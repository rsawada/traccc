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
    const strip_pair_config& config, unsigned int& pair_position,
    strip_pair_collection_types::view pairs_view) {

    const edm::measurement_collection<default_algebra>::const_device
        measurements(measurements_view);
    if (globalIndex >= measurements.size()) {
        return;
    }

    typename detector_t::device det(det_view);
    vecmem::device_vector<strip_pair> pairs(pairs_view);
    const edm::measurement inner_measurement = measurements.at(globalIndex);

    for (unsigned int other_index = 0u; other_index < measurements.size();
         ++other_index) {
        const edm::measurement outer_measurement = measurements.at(other_index);
        if (!details::is_compatible_barrel_strip_pair(
                det, inner_measurement, outer_measurement, config)) {
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
        const point2 inner_local_center{
            inner_measurement.local_position()[0], 0.f};
        const point2 outer_local_center{
            outer_measurement.local_position()[0], 0.f};
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
                                      inner_normal[2] * outer_normal[2]};
        }
    }
}

}  // namespace traccc::device
