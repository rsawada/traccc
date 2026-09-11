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

template <typename spacepoint_backend_t>
TRACCC_HOST_DEVICE inline void clear_detailed_strip_info(
    edm::spacepoint<spacepoint_backend_t>& sp) {
    sp.top_strip_vector() = vector3{0.f, 0.f, 0.f};
    sp.bottom_strip_vector() = vector3{0.f, 0.f, 0.f};
    sp.strip_center_distance() = vector3{0.f, 0.f, 0.f};
    sp.top_strip_center() = point3{0.f, 0.f, 0.f};
}

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

TRACCC_HOST_DEVICE inline scalar dot3(const vector3& a, const vector3& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/// Intersect strip lines with beam-spot planes, with explicit endpoint
/// allowances.
TRACCC_HOST_DEVICE inline bool make_strip_spacepoint(
    point3& spacepoint, const point3& first_center,
    const vector3& first_direction, const vector3& second_direction,
    const vector3& first_trajectory, const vector3& second_trajectory,
    const vector3& first_normal, const vector3& second_normal,
    const scalar first_half_length, const scalar second_half_length,
    const scalar strip_length_gap_tolerance,
    const scalar strip_length_tolerance) {

    const scalar strip_length_limit = scalar{1} + strip_length_tolerance;
    const scalar first_denominator = dot3(first_direction, second_normal);
    const scalar second_denominator = dot3(second_direction, first_normal);
    if ((first_half_length <= 0.f) || (second_half_length <= 0.f) ||
        (std::abs(first_denominator) < 1e-12f) ||
        (std::abs(second_denominator) < 1e-12f)) {
        return false;
    }

    const scalar a = -dot3(first_trajectory, second_normal);
    const scalar c = -dot3(second_trajectory, first_normal);
    const scalar first_one_over_strip = 0.5f / first_half_length;
    const scalar second_one_over_strip = 0.5f / second_half_length;
    const scalar first_pre_cut =
        strip_length_limit + first_one_over_strip * strip_length_gap_tolerance;
    const scalar second_pre_cut =
        strip_length_limit + second_one_over_strip * strip_length_gap_tolerance;

    if ((std::abs(a) > std::abs(first_denominator) * first_pre_cut) ||
        (std::abs(c) > std::abs(second_denominator) * second_pre_cut)) {
        return false;
    }

    scalar m = a / first_denominator;
    scalar n = c / second_denominator;

    if (strip_length_gap_tolerance != 0.f) {
        const scalar cs = dot3(first_direction, second_direction) *
                          first_one_over_strip * first_one_over_strip;
        if (std::abs(cs) < 1e-12f) {
            return false;
        }
        if ((m > strip_length_limit) || (n > strip_length_limit)) {
            scalar dm = m - 1.f;
            const scalar dmn = (n - 1.f) * cs;
            if (dmn > dm) {
                dm = dmn;
            }
            m -= dm;
            n -= dm / cs;
        } else if ((m < -strip_length_limit) || (n < -strip_length_limit)) {
            scalar dm = -(1.f + m);
            const scalar dmn = -(1.f + n) * cs;
            if (dmn > dm) {
                dm = dmn;
            }
            m += dm;
            n += dm / cs;
        }

        if ((std::abs(m) > strip_length_limit) ||
            (std::abs(n) > strip_length_limit)) {
            return false;
        }
    }

    spacepoint = first_center + (0.5f * m) * first_direction;
    return true;
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
    clear_detailed_strip_info(sp);
}

}  // namespace traccc::details
