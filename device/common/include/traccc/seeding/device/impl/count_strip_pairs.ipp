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
    if (globalIndex >= measurements.size()) {
        return;
    }

    typename detector_t::device det(det_view);
    const edm::measurement inner_measurement = measurements.at(globalIndex);
    unsigned int n_compatible_barrel_pairs = 0u;
    unsigned int n_compatible_endcap_pairs = 0u;
    unsigned int n_compatible_endcap_boundary_pairs = 0u;

    for (unsigned int other_index = 0u; other_index < measurements.size();
         ++other_index) {
        const edm::measurement outer_measurement = measurements.at(other_index);
        const detray::tracking_surface inner_surface{
            det, inner_measurement.surface_link()};
        const detray::tracking_surface outer_surface{
            det, outer_measurement.surface_link()};

        if ((static_cast<int>(inner_surface.shape_id()) == 0) &&
            (static_cast<int>(outer_surface.shape_id()) == 0)) {
            if (details::is_compatible_barrel_strip_pair(
                    det, inner_measurement, outer_measurement,
                    barrel_config)) {
                ++n_compatible_barrel_pairs;
            }
        } else {
            if (details::is_compatible_endcap_strip_pair(
                    det, inner_measurement, outer_measurement,
                    endcap_config)) {
                strip_measurement_surface_info inner_info{
                    inner_measurement.surface_link().value(), 0u, 0.f, 0.f,
                    0.f, 0.f};
                strip_measurement_surface_info outer_info{
                    outer_measurement.surface_link().value(), 0u, 0.f, 0.f,
                    0.f, 0.f};
                if ((globalIndex < surface_infos.size()) &&
                    (other_index < surface_infos.size())) {
                    inner_info = surface_infos.at(globalIndex);
                    outer_info = surface_infos.at(other_index);
                }

                if ((inner_info.is_endcap == 0u) ||
                    (outer_info.is_endcap == 0u)) {
                    continue;
                }

                const point2 inner_local_center{
                    inner_info.mid_r, inner_measurement.local_position()[1]};
                const point2 outer_local_center{
                    outer_info.mid_r, outer_measurement.local_position()[1]};
                const point3 inner_strip_center =
                    inner_surface.local_to_global({}, inner_local_center, {});
                const point3 outer_strip_center =
                    outer_surface.local_to_global({}, outer_local_center, {});
                const scalar delta_x =
                    inner_strip_center[0] - outer_strip_center[0];
                const scalar delta_y =
                    inner_strip_center[1] - outer_strip_center[1];
                const scalar strip_center_delta_xy =
                    std::sqrt(delta_x * delta_x + delta_y * delta_y);
                if (strip_center_delta_xy >=
                    endcap_config.max_strip_center_delta_xy) {
                    continue;
                }

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
