#ifndef slic3r_Magma_MagmaResolved_hpp_
#define slic3r_Magma_MagmaResolved_hpp_

#include "../PrintConfig.hpp"
#include "MagmaGeometry.hpp"

namespace Slic3r {
namespace magma {

// Every Magma value derived from config, resolved in one place so MagmaTubeMap::build,
// Print::validate and PresetHints all describe the same tube. Add a field here rather than
// re-deriving a quantity at a call site.
struct MagmaResolved
{
    // --- what the lattice is ---
    InfillPattern        pattern   = ipMagmaRectilinear;
    const MagmaGeometry *geometry  = nullptr;

    // --- inputs, after resolution ---
    int    sparse_extruder     = 0;    // 0-based index of the extruder printing the lattice
    double nozzle_diameter     = 0.0;  // that extruder's bore
    double line_width          = 0.0;  // deposited bead width ("auto" already resolved)
    // The injection nozzle's cone. There is no lattice-extruder tip flat: only the injection
    // nozzle is pressed into the print.
    double cone_half_angle_deg = 0.0;
    double seal_press          = 0.0;  // depth pressed past first contact, before injection

    // The injection extruder may differ from the lattice extruder. Sealing uses these fields;
    // printing the lattice uses sparse_extruder / nozzle_diameter above.
    int    injection_extruder       = 0;    // 0-based; falls back to sparse_extruder
    double injection_nozzle_diameter = 0.0;
    double injection_nozzle_flat     = 0.0;  // measured tip flat (or the bore-scaled stand-in)
    // True when no outer diameter was configured and injection_nozzle_flat is
    // resolve_nozzle_flat's bore-scaled guess. Print::validate refuses to slice in that state;
    // the GUI readout uses it to mark seal depths as estimates.
    bool   injection_flat_is_estimate = false;

    // --- derived geometry ---
    double interior_width   = 0.0;  // open face-to-face width of one cell
    double cell_spacing     = 0.0;  // centre-to-centre line spacing
    double opening_diameter = 0.0;  // circle the nozzle must cover to seal
    double bore_diameter    = 0.0;  // largest circle that fits inside the tube

    // --- tube height bounds ---
    // From the configured layer height, as MagmaTubeMap::build() computes them.
    double min_tube_height = 0.0;
    double max_tube_height = 0.0;

    // --- derived depths (nominal) ---
    // From the nominal opening. The G-code path recomputes these per tube from its clipped
    // opening, so use these for warnings and readouts only, never to emit a move.
    //
    // The nozzle drops to seal_depth before extruding, then sinks plunge_depth further during
    // the injection, reaching total_depth() as the last filament goes in.
    double seal_depth   = 0.0;  // derived from the tube width, not configured directly
    double plunge_depth = 0.0;  // additive, during the injection

    // Contact at the cell's corner, (press + plunge) * tan(theta): the corner is the last part
    // of the opening the cone covers, so this is the weakest point of the seal.
    double grip         = 0.0;
    // Cone diameter at seal depth (before plunge). Compared against the cell pitch only as a
    // geometry sanity bound (MAGMA_PITCH_ABSURD_RATIO), not as a quality predictor.
    double cone_at_seal = 0.0;

    // Plunge as configured, before clamp_plunge_depth(); lets the readout and
    // Print::validate report when it was reduced.
    double plunge_requested = 0.0;
    bool   plunge_clamped() const { return plunge_requested > plunge_depth + 1e-9; }

    // Ratio the grid-disruption warning fires on.
    double pitch_ratio() const {
        const double pitch = cell_spacing;
        return pitch > 0.0 ? cone_at_seal / pitch : 0.0;
    }

    // Depth below the print surface the nozzle reaches at the end of injection.
    double total_depth() const { return seal_depth + plunge_depth; }

};

// Returns false, leaving `out` untouched, when the effective pattern is not a Magma pattern.
// `region` supplies the pattern, lattice filament and line width; `object` the tube, nozzle
// tip and seal settings; `print` the nozzle diameters.
bool resolve_magma(const PrintRegionConfig &region,
                   const PrintObjectConfig &object,
                   const PrintConfig       &print,
                   MagmaResolved           &out);

} // namespace magma
} // namespace Slic3r

#endif // slic3r_Magma_MagmaResolved_hpp_
