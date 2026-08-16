/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/seeding/detail/spacepoint_formation.hpp"

// System include(s).
#include <cassert>

namespace traccc::device {

template <typename detector_t>
TRACCC_HOST_DEVICE inline void form_spacepoints(
    const global_index_t globalIndex, typename detector_t::view det_view,
    const edm::measurement_collection<default_algebra>::const_view&
        measurements_view,
    edm::spacepoint_collection::view spacepoints_view) {

    form_pixel_spacepoints<detector_t>(globalIndex, det_view, measurements_view,
                                       spacepoints_view);
}

template <typename detector_t>
TRACCC_HOST_DEVICE inline void form_pixel_spacepoints(
    const global_index_t globalIndex, typename detector_t::view det_view,
    const edm::measurement_collection<default_algebra>::const_view&
        measurements_view,
    edm::spacepoint_collection::view spacepoints_view) {

    // Set up the input container(s).
    const edm::measurement_collection<default_algebra>::const_device
        measurements(measurements_view);

    // Check if anything needs to be done
    if (globalIndex >= measurements.size()) {
        return;
    }

    // Create the tracking geometry
    typename detector_t::device det(det_view);

    // Set up the output container(s).
    edm::spacepoint_collection::device spacepoints(spacepoints_view);

    const edm::measurement meas = measurements.at(globalIndex);

    // Fill the spacepoint using the common function.
    if (details::is_valid_pixel_measurement(meas)) {
        const edm::spacepoint_collection::device::size_type i =
            spacepoints.push_back_default();
        edm::spacepoint_collection::device::proxy_type sp = spacepoints.at(i);
        traccc::details::fill_pixel_spacepoint(sp, det, meas);
        sp.measurement_index_1() = globalIndex;
        sp.measurement_index_2() =
            edm::spacepoint_collection::device::INVALID_MEASUREMENT_INDEX;
    }
}

template <typename detector_t>
TRACCC_HOST_DEVICE inline void form_strip_spacepoints(
    const global_index_t globalIndex, typename detector_t::view det_view,
    const edm::measurement_collection<default_algebra>::const_view&
        measurements_view,
    edm::spacepoint_collection::view spacepoints_view) {

    // Set up the input container(s).
    const edm::measurement_collection<default_algebra>::const_device
        measurements(measurements_view);

    // Check if anything needs to be done
    if (globalIndex >= measurements.size()) {
        return;
    }

    // Create the tracking geometry
    typename detector_t::device det(det_view);

    // Set up the output container(s).
    edm::spacepoint_collection::device spacepoints(spacepoints_view);

    const edm::measurement meas = measurements.at(globalIndex);

    // Simplified strip spacepoint: one 1D strip measurement creates one point.
    if (details::is_valid_strip_measurement(meas)) {
        const edm::spacepoint_collection::device::size_type i =
            spacepoints.push_back_default();
        edm::spacepoint_collection::device::proxy_type sp = spacepoints.at(i);
        traccc::details::fill_strip_spacepoint(sp, det, meas);
        sp.measurement_index_1() = globalIndex;
        sp.measurement_index_2() =
            edm::spacepoint_collection::device::INVALID_MEASUREMENT_INDEX;
    }
}

template <typename detector_t>
TRACCC_HOST_DEVICE inline void form_barrel_strip_spacepoints(
    const global_index_t globalIndex, typename detector_t::view det_view,
    const edm::measurement_collection<default_algebra>::const_view&
        measurements_view,
    const strip_pair_collection_types::const_view& pairs_view,
    const strip_measurement_surface_info_collection_types::const_view&
        surface_infos_view,
    const point3& beam_spot,
    edm::spacepoint_collection::view spacepoints_view) {

    const edm::measurement_collection<default_algebra>::const_device
        measurements(measurements_view);
    const strip_pair_collection_types::const_device pairs(pairs_view);
    const strip_measurement_surface_info_collection_types::const_device
        surface_infos(surface_infos_view);
    if (globalIndex >= pairs.size()) {
        return;
    }

    (void)det_view;
    edm::spacepoint_collection::device spacepoints(spacepoints_view);
    const strip_pair pair = pairs.at(globalIndex);
    const edm::measurement first_meas =
        measurements.at(pair.measurement_index_1);
    const edm::measurement second_meas =
        measurements.at(pair.measurement_index_2);

    strip_measurement_surface_info first_info{};
    strip_measurement_surface_info second_info{};
    if (!details::find_strip_surface_info(surface_infos, pair.surface_link_1,
                                          first_info) ||
        !details::find_strip_surface_info(surface_infos, pair.surface_link_2,
                                          second_info)) {
        return;
    }
    const details::strip_material first_material =
        details::make_strip_material(first_meas, first_info, beam_spot);
    const details::strip_material second_material =
        details::make_strip_material(second_meas, second_info, beam_spot);
    if ((first_material.valid == 0u) || (second_material.valid == 0u)) {
        return;
    }
    const scalar strip_length_gap_tolerance =
        details::g80_strip_length_gap_tolerance(first_info, second_info);
    point3 g80_spacepoint{};
    const bool accepted = traccc::details::make_g80_strip_spacepoint(
        g80_spacepoint, first_material.center, first_material.direction,
        second_material.direction, first_material.trajectory,
        second_material.trajectory, first_material.normal,
        second_material.normal, first_material.half_length,
        second_material.half_length, strip_length_gap_tolerance);
    if (!accepted) {
        return;
    }

    const edm::spacepoint_collection::device::size_type i =
        spacepoints.push_back_default();
    edm::spacepoint_collection::device::proxy_type sp = spacepoints.at(i);
    sp.x() = g80_spacepoint[0];
    sp.y() = g80_spacepoint[1];
    sp.z() = g80_spacepoint[2];
    sp.radius_variance() = 0.f;
    sp.z_variance() = 0.f;
    sp.measurement_index_1() = pair.measurement_index_1;
    sp.measurement_index_2() = pair.measurement_index_2;
}

}  // namespace traccc::device
