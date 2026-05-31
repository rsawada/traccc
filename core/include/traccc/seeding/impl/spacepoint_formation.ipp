/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"

// Detray include(s).
#include <detray/geometry/tracking_surface.hpp>

namespace traccc::details {

template <typename measurement_backend_t>
TRACCC_HOST_DEVICE inline bool is_valid_measurement(
    const edm::measurement<measurement_backend_t>& meas) {
    return is_valid_pixel_measurement(meas);
}

template <typename measurement_backend_t>
TRACCC_HOST_DEVICE inline bool is_valid_pixel_measurement(
    const edm::measurement<measurement_backend_t>& meas) {
    return (meas.dimensions() == 2u);
}

template <typename measurement_backend_t>
TRACCC_HOST_DEVICE inline bool is_valid_strip_measurement(
    const edm::measurement<measurement_backend_t>& meas) {
    return (meas.dimensions() == 1u);
}

template <typename spacepoint_backend_t, typename detector_t,
          typename measurement_backend_t>
TRACCC_HOST_DEVICE inline void fill_pixel_spacepoint(
    edm::spacepoint<spacepoint_backend_t>& sp, const detector_t& det,
    const edm::measurement<measurement_backend_t>& meas,
    const typename detector_t::geometry_context gctx) {

    // Get the global position of this silicon pixel measurement.
    const detray::tracking_surface sf{det, meas.surface_link()};
    const point3 global = sf.local_to_global(gctx, meas.local_position(), {});

    // Fill the spacepoint with the global position and the measurement.
    sp.x() = global[0];
    sp.y() = global[1];
    sp.z() = global[2];
    sp.radius_variance() = 0.f;
    sp.z_variance() = 0.f;
}

template <typename spacepoint_backend_t, typename detector_t,
          typename measurement_backend_t>
TRACCC_HOST_DEVICE inline void fill_barrel_strip_spacepoint(
    edm::spacepoint<spacepoint_backend_t>& sp, const detector_t& det,
    const edm::measurement<measurement_backend_t>& first_meas,
    const edm::measurement<measurement_backend_t>& second_meas,
    const typename detector_t::geometry_context gctx) {

    const detray::tracking_surface first_surface{det,
                                                  first_meas.surface_link()};
    const detray::tracking_surface second_surface{det,
                                                   second_meas.surface_link()};
    const point2 first_local_center{first_meas.local_position()[0], 0.f};
    const point2 second_local_center{second_meas.local_position()[0], 0.f};
    const point2 first_local_plus{first_meas.local_position()[0], 1.f};
    const point2 second_local_plus{second_meas.local_position()[0], 1.f};
    const point3 first_center =
        first_surface.local_to_global(gctx, first_local_center, {});
    const point3 second_center =
        second_surface.local_to_global(gctx, second_local_center, {});
    const point3 first_plus =
        first_surface.local_to_global(gctx, first_local_plus, {});
    const point3 second_plus =
        second_surface.local_to_global(gctx, second_local_plus, {});
    const vector3 first_direction = first_plus - first_center;
    const vector3 second_direction = second_plus - second_center;
    const vector3 center_delta = first_center - second_center;

    const scalar a = first_direction[0] * first_direction[0] +
                     first_direction[1] * first_direction[1] +
                     first_direction[2] * first_direction[2];
    const scalar b = first_direction[0] * second_direction[0] +
                     first_direction[1] * second_direction[1] +
                     first_direction[2] * second_direction[2];
    const scalar c = second_direction[0] * second_direction[0] +
                     second_direction[1] * second_direction[1] +
                     second_direction[2] * second_direction[2];
    const scalar d = first_direction[0] * center_delta[0] +
                     first_direction[1] * center_delta[1] +
                     first_direction[2] * center_delta[2];
    const scalar e = second_direction[0] * center_delta[0] +
                     second_direction[1] * center_delta[1] +
                     second_direction[2] * center_delta[2];
    const scalar denominator = a * c - b * b;

    // Parallel strip lines do not define a stable intersection point. Keep a
    // provisional midpoint so that one candidate pair still produces one SP.
    if (std::abs(denominator) < 1e-12f) {
        sp.x() = 0.5f * (first_center[0] + second_center[0]);
        sp.y() = 0.5f * (first_center[1] + second_center[1]);
        sp.z() = 0.5f * (first_center[2] + second_center[2]);
        sp.radius_variance() = 0.f;
        sp.z_variance() = 0.f;
        return;
    }

    const scalar first_parameter = (b * e - c * d) / denominator;
    const scalar second_parameter = (a * e - b * d) / denominator;
    const point3 first_closest = first_center + first_parameter * first_direction;
    const point3 second_closest =
        second_center + second_parameter * second_direction;

    // Use the midpoint of the closest points on the two infinite strip lines.
    sp.x() = 0.5f * (first_closest[0] + second_closest[0]);
    sp.y() = 0.5f * (first_closest[1] + second_closest[1]);
    sp.z() = 0.5f * (first_closest[2] + second_closest[2]);
    sp.radius_variance() = 0.f;
    sp.z_variance() = 0.f;
}

template <typename spacepoint_backend_t, typename detector_t,
          typename measurement_backend_t>
TRACCC_HOST_DEVICE inline void fill_strip_spacepoint(
    edm::spacepoint<spacepoint_backend_t>& sp, const detector_t& det,
    const edm::measurement<measurement_backend_t>& meas,
    const typename detector_t::geometry_context gctx) {

    // Get the global position of this silicon strip measurement.
    const detray::tracking_surface sf{det, meas.surface_link()};
    const point3 global = sf.local_to_global(gctx, meas.local_position(), {});

    // Fill the spacepoint with the global position and the measurement.
    sp.x() = global[0];
    sp.y() = global[1];
    sp.z() = global[2];
    sp.radius_variance() = 0.f;
    sp.z_variance() = 0.f;
}

}  // namespace traccc::details
