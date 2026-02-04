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
    char buf[256];

    // --- Config values ---
    int injection_temp       = config.magma_injection_temp.value;
    double injection_speed_vol = config.magma_injection_speed.value;
    bool iron_tube_ends      = config.magma_iron_tube_ends.value;
    bool park_enabled        = config.magma_injection_park.value;
    int dwell_ms             = config.magma_injection_dwell.value;
    double fill_factor       = config.magma_tube_fill_factor.value;
    bool z_slam              = config.magma_injection_z_slam.value;
    int inj_filament         = config.magma_injection_filament.value;

    // --- Filament switch (before any injection) ---
    unsigned int original_filament = gcodegen.writer().filament()->id();
    bool need_filament_switch = (inj_filament > 0 &&
                                 (unsigned int)(inj_filament - 1) != original_filament);

    if (need_filament_switch) {
        gcode += "; Magma injection: switching to injection filament\n";
        gcode += gcodegen.set_extruder(inj_filament - 1, layer_z);
    }

    // Get extruder info (may have changed after filament switch)
    unsigned int extruder_id = gcodegen.writer().filament()->id();
    double filament_diameter = config.filament_diameter.get_at(extruder_id);
    double filament_area = (M_PI / 4.0) * filament_diameter * filament_diameter;

    // Injection volumetric speed.
    // Default (0) uses the infill volumetric rate — this represents the
    // hotend's proven melt capacity during normal printing.
    double vol_speed;
    if (injection_speed_vol > 0) {
        vol_speed = injection_speed_vol;
    } else {
        double max_vol = config.filament_max_volumetric_speed.get_at(extruder_id);
        if (max_vol > 0) {
            vol_speed = max_vol;
        } else {
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
    double feedrate_mms = vol_speed / filament_area;
    double feedrate_mmmin = feedrate_mms * 60.0;

    // --- Heat up ---
    // Skip if filament switch already handled temperature via set_extruder()
    bool temp_changed = (injection_temp > 0 && injection_temp != print_temp
                         && !need_filament_switch);
    if (temp_changed) {
        gcode += "; Magma injection: heating\n";
        gcode += gcodegen.retract(false, false);

        if (park_enabled) {
            double park_z = layer_z + 2.0;
            gcode += gcodegen.writer().travel_to_z(park_z, "park z-hop for temp change");
        }

        gcode += gcodegen.writer().set_temperature(injection_temp, true);
    }

    // --- Injection loop ---
    // Emit role tag for GCode processor
    sprintf(buf, ";%s%s\n",
            GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Role).c_str(),
            ExtrusionEntity::role_to_string(erMagmaInjection).c_str());
    gcode += buf;

    // Force display dimensions for stationary injection extrusion
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

    constexpr double slam_depth = 0.1;  // mm, hardcoded for safety

    for (const auto& pt : points) {
        double volume = pt.volume_mm3 * fill_factor;
        if (volume <= 0)
            continue;

        // Travel to injection point
        Point scaled_pos(scale_(pt.position.x()), scale_(pt.position.y()));
        gcode += gcodegen.travel_to(scaled_pos, erMagmaInjection, "move to injection point");

        // Unretract before injection
        gcode += gcodegen.unretract();

        // Z-slam: lower nozzle into surface to seal against hole
        if (z_slam) {
            sprintf(buf, "G1 Z%.3f F600 ; z-slam seal\n", layer_z - slam_depth);
            gcode += buf;
        }

        // Stationary extrude (E-only, no XY movement)
        double filament_length = volume / filament_area;
        gcode += gcodegen.writer().set_speed(feedrate_mmmin);

        Vec3d cur_pos = gcodegen.writer().get_position();
        sprintf(buf, "Magma injection %.2f mm3", volume);
        gcode += gcodegen.writer().extrude_to_xy(
            Vec2d(cur_pos.x(), cur_pos.y()), filament_length, std::string(buf));

        // Dwell for air displacement
        if (dwell_ms > 0) {
            sprintf(buf, "G4 P%d ; injection dwell\n", dwell_ms);
            gcode += buf;
        }

        // Z-slam release: return to normal layer height
        if (z_slam) {
            sprintf(buf, "G1 Z%.3f F600 ; z-slam release\n", layer_z);
            gcode += buf;
        }

        // Retract after injection — relieves ooze pressure, triggers wipe
        gcode += gcodegen.retract(false, false);
    }

    // Reset forced dimensions
    sprintf(buf, ";%s0\n",
            GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Height).c_str());
    gcode += buf;
    sprintf(buf, ";%s0\n",
            GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Width).c_str());
    gcode += buf;

    // --- Ironing of tube ends ---
    if (iron_tube_ends) {
        // Emit ironing role tag
        sprintf(buf, ";%s%s\n",
                GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Role).c_str(),
                ExtrusionEntity::role_to_string(erIroning).c_str());
        gcode += buf;

        // Ironing parameters: use user's settings if ironing is configured,
        // otherwise use sensible defaults
        double ir_flow_pct = 0.10;
        double ir_spacing  = 0.1;
        double ir_speed    = 15.0 * 60.0;  // 15 mm/s → mm/min

        if (config.ironing_type.value != IroningType::NoIroning) {
            ir_flow_pct = config.ironing_flow.value / 100.0;
            ir_spacing  = config.ironing_spacing.value;
            ir_speed    = config.ironing_speed.value;
        }

        float nozzle_d = config.nozzle_diameter.get_at(extruder_id);
        double ir_height = config.layer_height.value * ir_flow_pct;
        double ir_flow_mm3_per_mm = nozzle_d * ir_height;
        double ir_e_per_mm = ir_flow_mm3_per_mm / filament_area;

        // Iron each injection hole with serpentine parallel lines
        for (const auto& pt : points) {
            double radius = tube_map.interior_width() / 2.0 - nozzle_d / 2.0;
            if (radius <= 0.05)
                continue;

            bool left_to_right = true;
            for (double dy = -radius; dy <= radius + 0.001; dy += ir_spacing) {
                double r2 = radius * radius - dy * dy;
                if (r2 <= 0)
                    continue;
                double half_chord = std::sqrt(r2);

                double x0 = pt.position.x() + (left_to_right ? -half_chord : half_chord);
                double x1 = pt.position.x() + (left_to_right ? half_chord : -half_chord);
                double y  = pt.position.y() + dy;

                // Travel to line start
                Point start_s(scale_(x0), scale_(y));
                gcode += gcodegen.travel_to(start_s, erIroning, "iron start");
                gcode += gcodegen.unretract();

                // Extrude ironing line
                double line_len = 2.0 * half_chord;
                double e_val = line_len * ir_e_per_mm;
                gcode += gcodegen.writer().set_speed(ir_speed);
                gcode += gcodegen.writer().extrude_to_xy(
                    Vec2d(x1, y), e_val, "tube iron");

                left_to_right = !left_to_right;
            }
            gcode += gcodegen.retract(false, false);
        }
    }

    // --- Cool down ---
    if (temp_changed) {
        gcode += "; Magma injection: cooling\n";
        gcode += gcodegen.retract(false, false);

        if (park_enabled) {
            double park_z = layer_z + 2.0;
            gcode += gcodegen.writer().travel_to_z(park_z, "park z-hop for temp restore");
        }

        gcode += gcodegen.writer().set_temperature(print_temp, true);
    }

    // --- Switch back to original filament ---
    if (need_filament_switch) {
        gcode += "; Magma injection: restoring print filament\n";
        gcode += gcodegen.set_extruder(original_filament, layer_z);
    }

    return gcode;
}

} // namespace magma
} // namespace Slic3r
