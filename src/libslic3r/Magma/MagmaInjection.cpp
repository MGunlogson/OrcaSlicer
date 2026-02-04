#include "MagmaInjection.hpp"
#include "MagmaTubeMap.hpp"
#include "MagmaSpiralOffset.hpp"

#include "../GCode.hpp"
#include "../GCode/GCodeProcessor.hpp"
#include "../Print.hpp"

#include <algorithm>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace Slic3r {
namespace magma {

std::vector<InjectionPoint> collect_injection_points(
    const MagmaTubeMap& tube_map,
    int layer_id)
{
    std::vector<InjectionPoint> points;

    const auto& pairs = tube_map.u_tube_pairs();
    for (int i = 0; i < static_cast<int>(pairs.size()); ++i) {
        const auto& pair = pairs[i];
        if (pair.pair_end_layer != layer_id)
            continue;
        if (pair.volume_mm3 <= 0)
            continue;

        InjectionPoint pt;
        // Use spiral-offset lattice for this layer so injection aligns with drawn pattern
        TriangleLattice lattice = lattice_for_layer(
            tube_map.cell_spacing(), tube_map.spiral_params(), layer_id);
        pt.position   = lattice.cell_center(pair.cell_a);
        pt.volume_mm3 = pair.volume_mm3;
        pt.pair_index = i;
        points.push_back(pt);
    }

    // Serpentine sort for minimal travel: sort by Y rows, alternate X direction.
    // Points on a triangle lattice naturally form rows; serpentine avoids the
    // long return sweep that a pure left-to-right raster would require.
    std::sort(points.begin(), points.end(), [](const InjectionPoint& a, const InjectionPoint& b) {
        if (std::abs(a.position.y() - b.position.y()) > 0.5)
            return a.position.y() < b.position.y();
        return a.position.x() < b.position.x();
    });

    // Assign row indices based on Y proximity, then reverse even rows
    if (points.size() > 1) {
        int row = 0;
        for (size_t i = 1; i < points.size(); ++i) {
            if (std::abs(points[i].position.y() - points[i - 1].position.y()) > 0.5)
                ++row;
        }
        // Re-walk and reverse alternating rows in-place
        size_t row_start = 0;
        row = 0;
        for (size_t i = 1; i <= points.size(); ++i) {
            bool end_of_row = (i == points.size()) ||
                              (std::abs(points[i].position.y() - points[i - 1].position.y()) > 0.5);
            if (end_of_row) {
                if (row % 2 == 1)
                    std::reverse(points.begin() + row_start, points.begin() + i);
                row_start = i;
                ++row;
            }
        }
    }

    return points;
}

std::string generate_injection_gcode(
    GCode& gcodegen,
    const MagmaTubeMap& tube_map,
    const std::vector<InjectionPoint>& points,
    double layer_z)
{
    if (points.empty())
        return {};

    std::string gcode;
    const auto& config = gcodegen.config();

    // --- Config values ---
    int injection_temp   = config.magma_injection_temp.value;
    double injection_speed_vol = config.magma_injection_speed.value; // mm^3/s, 0 = auto
    bool iron_tube_ends  = config.magma_iron_tube_ends.value;
    bool park_enabled    = config.magma_injection_park.value;
    int dwell_ms         = config.magma_injection_dwell.value;
    double fill_factor   = config.magma_tube_fill_factor.value;

    // Get current extruder info
    unsigned int extruder_id = gcodegen.writer().filament()->id();
    double filament_diameter = config.filament_diameter.get_at(extruder_id);
    double filament_area = (M_PI / 4.0) * filament_diameter * filament_diameter;

    // Injection volumetric speed.
    // Default (0) uses the infill volumetric rate — this represents the
    // hotend's proven melt capacity during normal printing.
    // Tube depth is calculated separately using a coupled pressure+thermal
    // model (see MagmaTubeMap::build). The extruder self-limits against
    // tube back-pressure: if commanded speed exceeds what pressure allows,
    // the extruder pushes slower. Adjust magma_tube_fill_factor if tubes
    // are underfilled due to this effect.
    double vol_speed;
    if (injection_speed_vol > 0) {
        vol_speed = injection_speed_vol;
    } else {
        double max_vol = config.filament_max_volumetric_speed.get_at(extruder_id);
        if (max_vol > 0) {
            vol_speed = max_vol;
        } else {
            // Derive from infill speed — proven melt rate for this hotend
            float nozzle_d = config.nozzle_diameter.get_at(extruder_id);
            float line_width = config.sparse_infill_line_width.get_abs_value(nozzle_d);
            if (line_width <= 0) line_width = nozzle_d;
            float lh = config.layer_height.value;
            double mm3_per_mm = lh * (line_width - lh * (1.0 - M_PI / 4.0));
            vol_speed = config.sparse_infill_speed.value * mm3_per_mm;
        }
    }

    // Cap at hotend melt rate
    double max_vol = config.filament_max_volumetric_speed.get_at(extruder_id);
    if (max_vol > 0)
        vol_speed = std::min(vol_speed, max_vol);

    // Clamp injection temp to safety max
    if (injection_temp > 0) {
        int max_temp = config.nozzle_temperature_range_high.get_at(extruder_id);
        if (max_temp > 0)
            injection_temp = std::min(injection_temp, max_temp);
    }

    // Get current printing temperature for restore
    int print_temp = config.nozzle_temperature.get_at(extruder_id);

    // Convert volumetric speed to filament feedrate
    double feedrate_mms = vol_speed / filament_area;   // mm/s of filament
    double feedrate_mmmin = feedrate_mms * 60.0;       // for G-code F parameter

    // --- Heat up ---
    bool temp_changed = (injection_temp > 0 && injection_temp != print_temp);
    if (temp_changed) {
        gcode += "; Magma injection: heating\n";
        gcode += gcodegen.retract(false, false);

        if (park_enabled) {
            // Park at front-left of bed (simple approach)
            // A more sophisticated version could use OozePrevention logic
            gcode += "; parking for temperature change\n";
        }

        gcode += gcodegen.writer().set_temperature(injection_temp, true); // M109, wait
    }

    // --- Injection ---
    // Emit role tag for GCode processor
    char buf[256];
    sprintf(buf, ";%s%s\n",
            GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Role).c_str(),
            ExtrusionEntity::role_to_string(erMagmaInjection).c_str());
    gcode += buf;

    // Force display dimensions for stationary injection extrusion.
    // GCodeProcessor uses these instead of computing from volume/distance
    // (which would be division by zero for E-only moves).
    float display_width = tube_map.interior_width();
    float display_height = static_cast<float>(config.layer_height.value);
    sprintf(buf, ";%s%g\n",
            GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Height).c_str(),
            display_height);
    gcode += buf;
    sprintf(buf, ";%s%g\n",
            GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Width).c_str(),
            display_width);
    gcode += buf;

    for (const auto& pt : points) {
        double volume = pt.volume_mm3 * fill_factor;
        if (volume <= 0)
            continue;

        // Travel to injection point
        Point scaled_pos(scale_(pt.position.x()), scale_(pt.position.y()));
        gcode += gcodegen.travel_to(scaled_pos, erMagmaInjection, "move to injection point");

        // Unretract before stationary extrude
        gcode += gcodegen.unretract();

        // Calculate filament length from volume
        double filament_length = volume / filament_area;

        // Set injection feedrate
        gcode += gcodegen.writer().set_speed(feedrate_mmmin);

        // Stationary extrude at current position (E-only, no XY movement)
        // Uses extrude_to_xy with current pos to maintain proper E tracking
        Vec3d cur_pos = gcodegen.writer().get_position();
        sprintf(buf, "Magma injection %.2f mm3", volume);
        gcode += gcodegen.writer().extrude_to_xy(
            Vec2d(cur_pos.x(), cur_pos.y()), filament_length, std::string(buf));

        // Dwell for air displacement
        if (dwell_ms > 0) {
            sprintf(buf, "G4 P%d ; injection dwell\n", dwell_ms);
            gcode += buf;
        }
    }

    // Reset forced dimensions so subsequent extrusion calculates normally
    sprintf(buf, ";%s0\n",
            GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Height).c_str());
    gcode += buf;
    sprintf(buf, ";%s0\n",
            GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Width).c_str());
    gcode += buf;

    // --- Ironing of tube ends ---
    // TODO: Implement tube end ironing when magma_iron_tube_ends is enabled.
    // Requires proper coordinate transformation through GCode's origin offset system.
    // Will emit rectilinear ironing passes using erIroning role and user's ironing settings.
    (void)iron_tube_ends;

    // --- Cool down ---
    if (temp_changed) {
        gcode += "; Magma injection: cooling\n";
        gcode += gcodegen.retract(false, false);

        if (park_enabled) {
            gcode += "; parking for temperature restore\n";
        }

        gcode += gcodegen.writer().set_temperature(print_temp, true); // M109, wait
    }

    return gcode;
}

} // namespace magma
} // namespace Slic3r
