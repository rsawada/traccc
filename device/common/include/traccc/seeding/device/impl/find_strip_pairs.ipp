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
            inner_measurement.surface_link().value(), 0u, 0.f, 0.f, 0.f, 0.f};
        strip_measurement_surface_info outer_info{
            outer_measurement.surface_link().value(), 0u, 0.f, 0.f, 0.f, 0.f};
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
        const scalar strip_center_delta_xy =
            std::sqrt(delta_x * delta_x + delta_y * delta_y);
        if (use_endcap_mid_r &&
            (strip_center_delta_xy >=
             endcap_config.max_strip_center_delta_xy)) {
            continue;
        }
        const scalar delta_z =
            inner_strip_center[2] - outer_strip_center[2];
        const scalar inner_strip_half_length = inner_info.strip_half_length;
        const scalar outer_strip_half_length = outer_info.strip_half_length;
        const vector3 inner_normal =
            inner_surface.normal({}, inner_measurement.local_position());
        const vector3 outer_normal =
            outer_surface.normal({}, outer_measurement.local_position());

        scalar strip_length_gap_tolerance = 0.f;
        const scalar strip_gap_parameter =
            use_endcap_mid_r ? endcap_config.strip_gap_parameter
                             : barrel_config.strip_gap_parameter;
        if (strip_gap_parameter != 0.f) {
            // Match the offline offset() calculation: T1(?,0) and T2(?,0)
            // are the local-x axes of the two detector-element transforms,
            // expressed in global coordinates. They are not the per-strip
            // local0 probe directions.
            const vector3 inner_axis = inner_surface.transform({}).x();
            const vector3 outer_axis = outer_surface.transform({}).x();
            const scalar inner_axis_norm = std::sqrt(
                inner_axis[0] * inner_axis[0] + inner_axis[1] * inner_axis[1] +
                inner_axis[2] * inner_axis[2]);
            const scalar outer_axis_norm = std::sqrt(
                outer_axis[0] * outer_axis[0] + outer_axis[1] * outer_axis[1] +
                outer_axis[2] * outer_axis[2]);
            if ((inner_axis_norm > 0.f) && (outer_axis_norm > 0.f)) {
                const scalar x12 =
                    (inner_axis[0] * outer_axis[0] +
                     inner_axis[1] * outer_axis[1] +
                     inner_axis[2] * outer_axis[2]) /
                    (inner_axis_norm * outer_axis_norm);
                const scalar s =
                    (inner_surface_center[0] - outer_surface_center[0]) *
                        inner_normal[0] +
                    (inner_surface_center[1] - outer_surface_center[1]) *
                        inner_normal[1] +
                    (inner_surface_center[2] - outer_surface_center[2]) *
                        inner_normal[2];
                const scalar dm =
                    strip_gap_parameter * inner_surface_r * std::abs(s * x12);

                if (use_endcap_mid_r) {
                    strip_length_gap_tolerance = dm / 0.04f;
                } else {
                    const scalar denom2 = (1.f - x12) * (1.f + x12);
                    if (denom2 > 1e-12f) {
                        strip_length_gap_tolerance = dm / std::sqrt(denom2);
                    }
                }

                if ((strip_length_gap_tolerance > 0.f) &&
                    (std::abs(inner_normal[2]) > 0.7f) &&
                    (std::abs(inner_surface_center[2]) > 1e-6f)) {
                    strip_length_gap_tolerance *=
                        inner_surface_r / std::abs(inner_surface_center[2]);
                }
            }
        }

        vecmem::device_atomic_ref<unsigned int> next_position(pair_position);
        const unsigned int position = next_position.fetch_add(1u);
        if (position < pairs.size()) {
            pairs.at(position) = {static_cast<unsigned int>(globalIndex),
                                  other_index,
                                  inner_measurement.surface_link().value(),
                                  outer_measurement.surface_link().value(),
                                  outer_surface_r - inner_surface_r,
                                  strip_center_delta_xy,
                                  delta_z,
                                  inner_normal[0] * outer_normal[0] +
                                      inner_normal[1] * outer_normal[1] +
                                      inner_normal[2] * outer_normal[2],
                                  use_endcap_mid_r ? 1u : 0u,
                                  inner_info.mid_r,
                                  outer_info.mid_r,
                                  inner_strip_half_length,
                                  outer_strip_half_length,
                                  strip_length_gap_tolerance};
        }
    }
}

}  // namespace traccc::device
