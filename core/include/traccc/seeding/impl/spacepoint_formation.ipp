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
    // Offline-inspired construction: build the plane containing the second
    // strip and the beam spot, then intersect the first strip with that plane.
    // The prototype uses the nominal beam spot at the global origin.
    const vector3 second_trajectory{second_center[0], second_center[1],
                                    second_center[2]};
    const vector3 second_plane_normal{
        second_direction[1] * second_trajectory[2] -
            second_direction[2] * second_trajectory[1],
        second_direction[2] * second_trajectory[0] -
            second_direction[0] * second_trajectory[2],
        second_direction[0] * second_trajectory[1] -
            second_direction[1] * second_trajectory[0]};
    const scalar denominator = first_direction[0] * second_plane_normal[0] +
                               first_direction[1] * second_plane_normal[1] +
                               first_direction[2] * second_plane_normal[2];

    // If the first strip is close to parallel to the plane, the intersection is
    // unstable. Keep a provisional midpoint so that one candidate pair still
    // produces one SP.
    if (std::abs(denominator) < 1e-12f) {
        sp.x() = 0.5f * (first_center[0] + second_center[0]);
        sp.y() = 0.5f * (first_center[1] + second_center[1]);
        sp.z() = 0.5f * (first_center[2] + second_center[2]);
        sp.radius_variance() = 0.f;
        sp.z_variance() = 0.f;
        return;
    }

    const scalar first_parameter =
        -(first_center[0] * second_plane_normal[0] +
          first_center[1] * second_plane_normal[1] +
          first_center[2] * second_plane_normal[2]) /
        denominator;
    const point3 spacepoint = first_center + first_parameter * first_direction;

    sp.x() = spacepoint[0];
    sp.y() = spacepoint[1];
    sp.z() = spacepoint[2];
    sp.radius_variance() = 0.f;
    sp.z_variance() = 0.f;
}

template <typename spacepoint_backend_t, typename detector_t,
          typename measurement_backend_t>
TRACCC_HOST_DEVICE inline void fill_endcap_strip_spacepoint(
    edm::spacepoint<spacepoint_backend_t>& sp, const detector_t& det,
    const edm::measurement<measurement_backend_t>& first_meas,
    const edm::measurement<measurement_backend_t>& second_meas,
    const typename detector_t::geometry_context gctx) {

    const detray::tracking_surface first_surface{det,
                                                  first_meas.surface_link()};
    const detray::tracking_surface second_surface{det,
                                                   second_meas.surface_link()};

    const auto first_min_r_values = first_surface.boundary(0u);
    const auto first_max_r_values = first_surface.boundary(1u);
    const auto second_min_r_values = second_surface.boundary(0u);
    const auto second_max_r_values = second_surface.boundary(1u);

    const scalar first_min_r =
        !first_min_r_values.empty() ? first_min_r_values.front() : 0.f;
    const scalar first_max_r =
        !first_max_r_values.empty() ? first_max_r_values.front() : first_min_r;
    const scalar second_min_r =
        !second_min_r_values.empty() ? second_min_r_values.front() : 0.f;
    const scalar second_max_r =
        !second_max_r_values.empty() ? second_max_r_values.front() : second_min_r;

    // Prototype endpoint model: annulus strip is local1 fixed, local0 varying.
    // Use the radial-band midpoint as the line anchor and +1 in local0 as the
    // direction. A later version can solve the exact min/max-r intersections.
    const scalar first_mid_r = 0.5f * (first_min_r + first_max_r);
    const scalar second_mid_r = 0.5f * (second_min_r + second_max_r);
    const point2 first_local_center{first_mid_r, first_meas.local_position()[1]};
    const point2 second_local_center{second_mid_r,
                                     second_meas.local_position()[1]};
    const point2 first_local_plus{first_mid_r + 1.f,
                                  first_meas.local_position()[1]};
    const point2 second_local_plus{second_mid_r + 1.f,
                                   second_meas.local_position()[1]};

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
    // Offline-inspired construction: build the plane containing the second
    // strip and the beam spot, then intersect the first strip with that plane.
    // The prototype uses the nominal beam spot at the global origin.
    const vector3 second_trajectory{second_center[0], second_center[1],
                                    second_center[2]};
    const vector3 second_plane_normal{
        second_direction[1] * second_trajectory[2] -
            second_direction[2] * second_trajectory[1],
        second_direction[2] * second_trajectory[0] -
            second_direction[0] * second_trajectory[2],
        second_direction[0] * second_trajectory[1] -
            second_direction[1] * second_trajectory[0]};
    const scalar denominator = first_direction[0] * second_plane_normal[0] +
                               first_direction[1] * second_plane_normal[1] +
                               first_direction[2] * second_plane_normal[2];

    // If the first strip is close to parallel to the plane, the intersection is
    // unstable. Keep a provisional midpoint so that one candidate pair still
    // produces one SP.
    if (std::abs(denominator) < 1e-12f) {
        sp.x() = 0.5f * (first_center[0] + second_center[0]);
        sp.y() = 0.5f * (first_center[1] + second_center[1]);
        sp.z() = 0.5f * (first_center[2] + second_center[2]);
        sp.radius_variance() = 0.f;
        sp.z_variance() = 0.f;
        return;
    }

    const scalar first_parameter =
        -(first_center[0] * second_plane_normal[0] +
          first_center[1] * second_plane_normal[1] +
          first_center[2] * second_plane_normal[2]) /
        denominator;
    const point3 spacepoint = first_center + first_parameter * first_direction;

    sp.x() = spacepoint[0];
    sp.y() = spacepoint[1];
    sp.z() = spacepoint[2];
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
