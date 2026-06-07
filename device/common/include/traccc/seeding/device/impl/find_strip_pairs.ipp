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
    if (globalIndex >= measurements.size()) {
        return;
    }

    typename detector_t::device det(det_view);
    vecmem::device_vector<strip_pair> pairs(pairs_view);
    const edm::measurement inner_measurement = measurements.at(globalIndex);

    for (unsigned int other_index = 0u; other_index < measurements.size();
         ++other_index) {
        const edm::measurement outer_measurement = measurements.at(other_index);
        if (!details::is_compatible_strip_pair(
                det, inner_measurement, outer_measurement, barrel_config,
                endcap_config)) {
            continue;
        }

        const detray::tracking_surface inner_surface{
            det, inner_measurement.surface_link()};
        const detray::tracking_surface outer_surface{
            det, outer_measurement.surface_link()};
        const point3 inner_surface_center = inner_surface.center({});
        const point3 outer_surface_center = outer_surface.center({});
        const scalar inner_surface_r =
            std::sqrt(inner_surface_center[0] * inner_surface_center[0] +
                      inner_surface_center[1] * inner_surface_center[1]);
        const scalar outer_surface_r =
            std::sqrt(outer_surface_center[0] * outer_surface_center[0] +
                      outer_surface_center[1] * outer_surface_center[1]);
        strip_measurement_surface_info inner_info{
            inner_measurement.surface_link().value(), 0u, 0.f, 0.f, 0.f};
        strip_measurement_surface_info outer_info{
            outer_measurement.surface_link().value(), 0u, 0.f, 0.f, 0.f};
        if ((globalIndex < surface_infos.size()) &&
            (other_index < surface_infos.size())) {
            inner_info = surface_infos.at(globalIndex);
            outer_info = surface_infos.at(other_index);
        }

        const bool use_endcap_mid_r =
            (inner_info.is_endcap != 0u) && (outer_info.is_endcap != 0u);
        const point2 inner_local_center{
            use_endcap_mid_r ? inner_info.mid_r
                             : inner_measurement.local_position()[0],
            use_endcap_mid_r ? inner_measurement.local_position()[1] : 0.f};
        const point2 outer_local_center{
            use_endcap_mid_r ? outer_info.mid_r
                             : outer_measurement.local_position()[0],
            use_endcap_mid_r ? outer_measurement.local_position()[1] : 0.f};
        const point3 inner_strip_center =
            inner_surface.local_to_global({}, inner_local_center, {});
        const point3 outer_strip_center =
            outer_surface.local_to_global({}, outer_local_center, {});
        const scalar delta_x =
            inner_strip_center[0] - outer_strip_center[0];
        const scalar delta_y =
            inner_strip_center[1] - outer_strip_center[1];
        const scalar delta_z =
            inner_strip_center[2] - outer_strip_center[2];
        const vector3 inner_normal =
            inner_surface.normal({}, inner_measurement.local_position());
        const vector3 outer_normal =
            outer_surface.normal({}, outer_measurement.local_position());

        vecmem::device_atomic_ref<unsigned int> next_position(pair_position);
        const unsigned int position = next_position.fetch_add(1u);
        if (position < pairs.size()) {
            pairs.at(position) = {static_cast<unsigned int>(globalIndex),
                                  other_index,
                                  inner_measurement.surface_link().value(),
                                  outer_measurement.surface_link().value(),
                                  outer_surface_r - inner_surface_r,
                                  std::sqrt(delta_x * delta_x +
                                            delta_y * delta_y),
                                  delta_z,
                                  inner_normal[0] * outer_normal[0] +
                                      inner_normal[1] * outer_normal[1] +
                                      inner_normal[2] * outer_normal[2],
                                  use_endcap_mid_r ? 1u : 0u,
                                  inner_info.mid_r,
                                  outer_info.mid_r};
        }
    }
}

}  // namespace traccc::device
