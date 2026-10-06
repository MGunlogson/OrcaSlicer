#ifndef slic3r_PresetHints_hpp_
#define slic3r_PresetHints_hpp_

#include <string>

#include "libslic3r/PresetBundle.hpp"

namespace Slic3r {

// GUI utility functions to produce hint messages from the current profile.
class PresetHints
{
public:
    // Produce a textual description of the cooling logic of a currently active filament.
    static std::string cooling_description(const Preset &preset);
    
    // Produce a textual description of the maximum flow achived for the current configuration
    // (the current printer, filament and print settigns).
    // This description will be useful for getting a gut feeling for the maximum volumetric
    // print speed achievable with the extruder.
    static std::string maximum_volumetric_flow_description(const PresetBundle &preset_bundle);

    // Produce a textual description of a recommended thin wall thickness
    // from the provided number of perimeters and the external / internal perimeter width.
    static std::string recommended_thin_wall_thickness(const PresetBundle &preset_bundle);

    // Produce a textual explanation of the combined effects of the top/bottom_shell_layers
    // versus top/bottom_min_shell_thickness. Which of the two values wins depends
    // on the active layer height.
    static std::string top_bottom_shell_thickness_explanation(const PresetBundle &preset_bundle);

    // `text` is the terse "Label:  value" block shown inline on the settings page;
    // `tooltip` repeats the rows with an explanation of each.
    struct MagmaReadout {
        std::string text;
        std::string tooltip;
    };

    // Magma: resolved tube geometry (interior width, cell spacing, bore, seal opening, open area).
    static MagmaReadout magma_geometry_readout(const PresetBundle &preset_bundle);

    // Magma: resolved injection depths (seal depth, corner grip, plunge, total immersion), nozzle vs
    // cell pitch, estimated dose and injection time.
    static MagmaReadout magma_injection_readout(const PresetBundle &preset_bundle);
};

} // namespace Slic3r

#endif /* slic3r_PresetHints_hpp_ */
