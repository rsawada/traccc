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
    edm::spacepoint_collection::view spacepoints_view) {

    const edm::measurement_collection<default_algebra>::const_device
        measurements(measurements_view);
    const strip_pair_collection_types::const_device pairs(pairs_view);
    const strip_measurement_surface_info_collection_types::const_device
        surface_infos(surface_infos_view);
    if (globalIndex >= pairs.size()) {
        return;
    }

    typename detector_t::device det(det_view);
    edm::spacepoint_collection::device spacepoints(spacepoints_view);
    const strip_pair pair = pairs.at(globalIndex);
    const edm::measurement first_meas =
        measurements.at(pair.measurement_index_1);
    const edm::measurement second_meas =
        measurements.at(pair.measurement_index_2);

    const strip_measurement_surface_info first_material =
        surface_infos.at(pair.measurement_index_1);
    const strip_measurement_surface_info second_material =
        surface_infos.at(pair.measurement_index_2);
    const scalar strip_length_gap_tolerance =
        details::g80_strip_length_gap_tolerance(first_material,
                                                second_material);
    point3 g80_spacepoint{};
    bool use_g80_spacepoint = false;
    bool accepted = true;
    if (pair.is_endcap == 0u) {
        if ((first_material.has_barrel_material != 0u) &&
            (second_material.has_barrel_material != 0u)) {
            use_g80_spacepoint = true;
            accepted = traccc::details::make_g80_strip_spacepoint(
                g80_spacepoint, first_material.barrel_strip_center,
                first_material.barrel_strip_direction,
                second_material.barrel_strip_direction,
                first_material.barrel_trajectory_direction,
                second_material.barrel_trajectory_direction,
                first_material.barrel_strip_normal,
                second_material.barrel_strip_normal,
                first_material.strip_half_length,
                second_material.strip_half_length,
                strip_length_gap_tolerance);
        }
    } else {
        if ((first_material.has_endcap_material != 0u) &&
            (second_material.has_endcap_material != 0u)) {
            use_g80_spacepoint = true;
            accepted = traccc::details::make_g80_strip_spacepoint(
                g80_spacepoint, first_material.endcap_strip_center,
                first_material.endcap_strip_direction,
                second_material.endcap_strip_direction,
                first_material.endcap_trajectory_direction,
                second_material.endcap_trajectory_direction,
                first_material.endcap_strip_normal,
                second_material.endcap_strip_normal,
                first_material.strip_half_length,
                second_material.strip_half_length,
                strip_length_gap_tolerance);
        }
    }

    if (!accepted) {
        return;
    }

    const edm::spacepoint_collection::device::size_type i =
        spacepoints.push_back_default();
    edm::spacepoint_collection::device::proxy_type sp = spacepoints.at(i);
    if (use_g80_spacepoint) {
        sp.x() = g80_spacepoint[0];
        sp.y() = g80_spacepoint[1];
        sp.z() = g80_spacepoint[2];
        sp.radius_variance() = 0.f;
        sp.z_variance() = 0.f;
    } else if (pair.is_endcap == 0u) {
        traccc::details::fill_barrel_strip_spacepoint(
            sp, det, first_meas, second_meas, {}, pair.strip_half_length_1,
            pair.strip_half_length_2, strip_length_gap_tolerance);
    } else {
        traccc::details::fill_endcap_strip_spacepoint(
            sp, det, first_meas, second_meas, {}, pair.surface_mid_r_1,
            pair.surface_mid_r_2, pair.strip_half_length_1,
            pair.strip_half_length_2, strip_length_gap_tolerance);
    }
    sp.measurement_index_1() = pair.measurement_index_1;
    sp.measurement_index_2() = pair.measurement_index_2;
}

}  // namespace traccc::device
