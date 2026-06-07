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
    edm::spacepoint_collection::view spacepoints_view) {

    const edm::measurement_collection<default_algebra>::const_device
        measurements(measurements_view);
    const strip_pair_collection_types::const_device pairs(pairs_view);
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

    const edm::spacepoint_collection::device::size_type i =
        spacepoints.push_back_default();
    edm::spacepoint_collection::device::proxy_type sp = spacepoints.at(i);

    const detray::tracking_surface first_surface{det,
                                                  first_meas.surface_link()};
    const detray::tracking_surface second_surface{det,
                                                   second_meas.surface_link()};
    if ((static_cast<int>(first_surface.shape_id()) == 0) &&
        (static_cast<int>(second_surface.shape_id()) == 0)) {
        traccc::details::fill_barrel_strip_spacepoint(sp, det, first_meas,
                                                       second_meas);
    } else {
        traccc::details::fill_endcap_strip_spacepoint(sp, det, first_meas,
                                                       second_meas);
    }
    sp.measurement_index_1() = pair.measurement_index_1;
    sp.measurement_index_2() = pair.measurement_index_2;
}

}  // namespace traccc::device
