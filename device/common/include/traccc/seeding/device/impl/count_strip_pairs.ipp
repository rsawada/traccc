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
TRACCC_HOST_DEVICE inline void count_strip_pairs(
    const global_index_t globalIndex, typename detector_t::view det_view,
    const edm::measurement_collection<default_algebra>::const_view&
        measurements_view,
    const strip_measurement_surface_info_collection_types::const_view&
        surface_infos_view,
    const barrel_strip_pair_config& barrel_config,
    const endcap_strip_pair_config& endcap_config, unsigned int& n_pairs,
    unsigned int& n_barrel_pairs, unsigned int& n_endcap_pairs,
    unsigned int& n_endcap_boundary_pairs) {

    const edm::measurement_collection<default_algebra>::const_device
        measurements(measurements_view);
    const strip_measurement_surface_info_collection_types::const_device
        surface_infos(surface_infos_view);
    (void)det_view;
    if (globalIndex >= measurements.size()) {
        return;
    }

    const edm::measurement inner_measurement = measurements.at(globalIndex);
    unsigned int n_compatible_barrel_pairs = 0u;
    unsigned int n_compatible_endcap_pairs = 0u;
    unsigned int n_compatible_endcap_boundary_pairs = 0u;

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

        if ((inner_info.has_barrel_material != 0u) &&
            (outer_info.has_barrel_material != 0u)) {
            if (details::is_compatible_barrel_strip_pair(
                    inner_measurement, outer_measurement, inner_info,
                    outer_info, barrel_config)) {
                ++n_compatible_barrel_pairs;
            }
        } else if ((inner_info.has_endcap_material != 0u) &&
                   (outer_info.has_endcap_material != 0u)) {
            if (details::is_compatible_endcap_strip_pair(
                    inner_measurement, outer_measurement, inner_info,
                    outer_info, endcap_config)) {
                ++n_compatible_endcap_pairs;
                ++n_compatible_endcap_boundary_pairs;
            }
        }
    }

    const unsigned int n_compatible_pairs =
        n_compatible_barrel_pairs + n_compatible_endcap_pairs;
    if (n_compatible_pairs > 0u) {
        vecmem::device_atomic_ref<unsigned int> total_pairs(n_pairs);
        total_pairs.fetch_add(n_compatible_pairs);
    }
    if (n_compatible_barrel_pairs > 0u) {
        vecmem::device_atomic_ref<unsigned int> barrel_pairs(n_barrel_pairs);
        barrel_pairs.fetch_add(n_compatible_barrel_pairs);
    }
    if (n_compatible_endcap_pairs > 0u) {
        vecmem::device_atomic_ref<unsigned int> endcap_pairs(n_endcap_pairs);
        endcap_pairs.fetch_add(n_compatible_endcap_pairs);
    }
    if (n_compatible_endcap_boundary_pairs > 0u) {
        vecmem::device_atomic_ref<unsigned int> endcap_boundary_pairs(
            n_endcap_boundary_pairs);
        endcap_boundary_pairs.fetch_add(n_compatible_endcap_boundary_pairs);
    }
}

}  // namespace traccc::device
