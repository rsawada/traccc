/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// VecMem include(s).
#include <vecmem/memory/device_atomic_ref.hpp>

namespace traccc::device {

template <typename detector_t>
TRACCC_HOST_DEVICE inline void find_strip_pairs(
    const global_index_t globalIndex, typename detector_t::view det_view,
    const edm::measurement_collection<default_algebra>::const_view&
        measurements_view,
    const strip_measurement_surface_info_collection_types::const_view&
        surface_infos_view,
    const barrel_strip_pair_config& barrel_config,
    const endcap_strip_pair_config& endcap_config, unsigned int& pair_position,
    strip_pair_collection_types::view pairs_view) {

    const edm::measurement_collection<default_algebra>::const_device
        measurements(measurements_view);
    const strip_measurement_surface_info_collection_types::const_device
        surface_infos(surface_infos_view);
    (void)det_view;
    if (globalIndex >= measurements.size()) {
        return;
    }

    vecmem::device_vector<strip_pair> pairs(pairs_view);
    const edm::measurement inner_measurement = measurements.at(globalIndex);

    for (unsigned int other_index = 0u; other_index < measurements.size();
         ++other_index) {
        if (other_index == globalIndex) {
            continue;
        }

        const edm::measurement outer_measurement = measurements.at(other_index);
        strip_measurement_surface_info inner_info{};
        strip_measurement_surface_info outer_info{};
        if ((globalIndex < surface_infos.size()) &&
            (other_index < surface_infos.size())) {
            inner_info = surface_infos.at(globalIndex);
            outer_info = surface_infos.at(other_index);
        }

        const bool use_barrel_material =
            (inner_info.has_barrel_material != 0u) &&
            (outer_info.has_barrel_material != 0u);
        const bool use_endcap_material =
            (inner_info.has_endcap_material != 0u) &&
            (outer_info.has_endcap_material != 0u);
        if (use_barrel_material) {
            if (!details::is_compatible_barrel_strip_pair(
                    inner_measurement, outer_measurement, inner_info,
                    outer_info, barrel_config)) {
                continue;
            }
        } else if (use_endcap_material) {
            if (!details::is_compatible_endcap_strip_pair(
                    inner_measurement, outer_measurement, inner_info,
                    outer_info, endcap_config)) {
                continue;
            }
        } else {
            continue;
        }

        const point3 inner_strip_center =
            use_barrel_material ? inner_info.barrel_strip_center
                                : inner_info.endcap_strip_center;
        const point3 outer_strip_center =
            use_barrel_material ? outer_info.barrel_strip_center
                                : outer_info.endcap_strip_center;
        const scalar inner_strip_r =
            std::sqrt(inner_strip_center[0] * inner_strip_center[0] +
                      inner_strip_center[1] * inner_strip_center[1]);
        const scalar outer_strip_r =
            std::sqrt(outer_strip_center[0] * outer_strip_center[0] +
                      outer_strip_center[1] * outer_strip_center[1]);
        const scalar delta_x = inner_strip_center[0] - outer_strip_center[0];
        const scalar delta_y = inner_strip_center[1] - outer_strip_center[1];
        const scalar strip_center_delta_xy =
            std::sqrt(delta_x * delta_x + delta_y * delta_y);
        const scalar delta_z = inner_strip_center[2] - outer_strip_center[2];
        const scalar inner_strip_half_length = inner_info.strip_half_length;
        const scalar outer_strip_half_length = outer_info.strip_half_length;
        const vector3 inner_normal =
            use_barrel_material ? inner_info.barrel_strip_normal
                                : inner_info.endcap_strip_normal;
        const vector3 outer_normal =
            use_barrel_material ? outer_info.barrel_strip_normal
                                : outer_info.endcap_strip_normal;

        const scalar strip_length_gap_tolerance =
            details::g80_strip_length_gap_tolerance(inner_info, outer_info);
        vecmem::device_atomic_ref<unsigned int> next_position(pair_position);
        const unsigned int position = next_position.fetch_add(1u);
        if (position < pairs.size()) {
            pairs.at(position) = {static_cast<unsigned int>(globalIndex),
                                  other_index,
                                  inner_measurement.surface_link().value(),
                                  outer_measurement.surface_link().value(),
                                  outer_strip_r - inner_strip_r,
                                  strip_center_delta_xy,
                                  delta_z,
                                  inner_normal[0] * outer_normal[0] +
                                      inner_normal[1] * outer_normal[1] +
                                      inner_normal[2] * outer_normal[2],
                                  use_endcap_material ? 1u : 0u,
                                  inner_info.mid_r,
                                  outer_info.mid_r,
                                  inner_strip_half_length,
                                  outer_strip_half_length,
                                  strip_length_gap_tolerance};
        }
    }
}

}  // namespace traccc::device
