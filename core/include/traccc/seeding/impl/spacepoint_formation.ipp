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


TRACCC_HOST_DEVICE inline scalar dot3(const vector3& a, const vector3& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

TRACCC_HOST_DEVICE inline scalar norm3(const vector3& a) {
    return std::sqrt(dot3(a, a));
}

TRACCC_HOST_DEVICE inline vector3 cross3(const vector3& a, const vector3& b) {
    return {a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]};
}

TRACCC_HOST_DEVICE inline scalar max_scalar(const scalar a, const scalar b) {
    return a > b ? a : b;
}

TRACCC_HOST_DEVICE inline scalar clamp_scalar(const scalar value,
                                             const scalar lo,
                                             const scalar hi) {
    return value < lo ? lo : (value > hi ? hi : value);
}

/// Calculate a strip spacepoint using the same cuts and line/plane
/// parameterisation as the G80 offline strip spacepoint implementation.
TRACCC_HOST_DEVICE inline bool make_g80_strip_spacepoint(
    point3& spacepoint, const point3& first_center,
    const vector3& first_direction, const vector3& second_direction,
    const vector3& first_trajectory,
    const vector3& second_trajectory, const vector3& first_normal,
    const vector3& second_normal, const scalar first_half_length,
    const scalar second_half_length,
    const scalar strip_length_gap_tolerance) {

    // Match G80's `double limit = 1. + float_tolerance` evaluation.
    constexpr scalar strip_length_tolerance = static_cast<scalar>(0.01f);
    constexpr scalar strip_length_limit =
        scalar{1} + strip_length_tolerance;
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
    const scalar first_pre_cut = strip_length_limit +
                                 first_one_over_strip *
                                     strip_length_gap_tolerance;
    const scalar second_pre_cut = strip_length_limit +
                                  second_one_over_strip *
                                      strip_length_gap_tolerance;

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
        } else if ((m < -strip_length_limit) ||
                   (n < -strip_length_limit)) {
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

template <typename spacepoint_backend_t>
TRACCC_HOST_DEVICE inline void fill_corrected_strip_spacepoint_from_lines(
    edm::spacepoint<spacepoint_backend_t>& sp, const point3& first_center,
    const point3& second_center, const vector3& first_direction,
    const vector3& second_direction, const scalar first_half_length,
    const scalar second_half_length, const scalar strip_length_gap_tolerance) {

    const vector3 first_trajectory{first_center[0], first_center[1],
                                   first_center[2]};
    const vector3 second_trajectory{second_center[0], second_center[1],
                                    second_center[2]};
    const vector3 first_plane_normal = cross3(first_direction, first_trajectory);
    const vector3 second_plane_normal =
        cross3(second_direction, second_trajectory);
    const scalar first_denominator = dot3(first_direction, second_plane_normal);

    if (std::abs(first_denominator) < 1e-12f) {
        sp.x() = 0.5f * (first_center[0] + second_center[0]);
        sp.y() = 0.5f * (first_center[1] + second_center[1]);
        sp.z() = 0.5f * (first_center[2] + second_center[2]);
        sp.radius_variance() = 0.f;
        sp.z_variance() = 0.f;
        return;
    }

    const vector3 first_center_vector{first_center[0], first_center[1],
                                      first_center[2]};
    scalar first_parameter =
        -dot3(first_center_vector, second_plane_normal) / first_denominator;

    const scalar second_denominator = dot3(second_direction, first_plane_normal);
    scalar second_parameter = 0.f;
    const bool has_second_parameter = std::abs(second_denominator) >= 1e-12f;
    if (has_second_parameter) {
        const vector3 second_center_vector{second_center[0], second_center[1],
                                           second_center[2]};
        second_parameter =
            -dot3(second_center_vector, first_plane_normal) / second_denominator;
    }

    // Offline-inspired strip-length correction. G80 works with a dimensionless
    // m/n parameter where +/-1 corresponds to the strip endpoints. Here the
    // line directions are local-coordinate steps, so parameter/half_length is
    // the corresponding dimensionless quantity.
    // Match G80's `double limit = 1. + float_tolerance` evaluation.
    constexpr scalar strip_length_tolerance = static_cast<scalar>(0.01f);
    constexpr scalar strip_length_limit =
        scalar{1} + strip_length_tolerance;
    if ((first_half_length > 0.f) && (second_half_length > 0.f)) {
        scalar first_q = first_parameter / first_half_length;
        scalar second_q = second_parameter / second_half_length;
        const scalar first_l =
            strip_length_limit +
            strip_length_gap_tolerance / (2.f * first_half_length);
        const scalar second_l =
            strip_length_limit +
            strip_length_gap_tolerance / (2.f * second_half_length);
        const scalar first_norm = norm3(first_direction);
        const scalar second_norm = norm3(second_direction);
        scalar correction_coupling = 0.f;
        if ((first_norm > 0.f) && (second_norm > 0.f)) {
            const scalar cos_angle = dot3(first_direction, second_direction) /
                                     (first_norm * second_norm);
            correction_coupling =
                cos_angle * (second_half_length / first_half_length);
        }

        if ((strip_length_gap_tolerance != 0.f) && has_second_parameter &&
            (std::abs(correction_coupling) > 1e-6f)) {
            if ((first_q > strip_length_limit) ||
                (second_q > strip_length_limit)) {
                scalar delta = first_q - 1.f;
                const scalar delta_second = (second_q - 1.f) * correction_coupling;
                if (delta_second > delta) {
                    delta = delta_second;
                }
                first_q -= delta;
                second_q -= delta / correction_coupling;
            } else if ((first_q < -strip_length_limit) ||
                       (second_q < -strip_length_limit)) {
                scalar delta = -(1.f + first_q);
                const scalar delta_second = -(1.f + second_q) * correction_coupling;
                if (delta_second > delta) {
                    delta = delta_second;
                }
                first_q += delta;
                second_q += delta / correction_coupling;
            }
        }

        // The current device kernel has already pushed an output slot, so it
        // cannot reject like G80. Keep the result finite and in the G80 pre-cut
        // acceptance envelope as a temporary fallback.
        first_q = clamp_scalar(first_q, -first_l, first_l);
        second_q = clamp_scalar(second_q, -second_l, second_l);
        first_parameter = first_q * first_half_length;
    } else if (first_half_length > 0.f) {
        first_parameter = clamp_scalar(first_parameter,
                                       -strip_length_limit * first_half_length,
                                       strip_length_limit * first_half_length);
    }

    const point3 spacepoint = first_center + first_parameter * first_direction;
    sp.x() = spacepoint[0];
    sp.y() = spacepoint[1];
    sp.z() = spacepoint[2];
    sp.radius_variance() = 0.f;
    sp.z_variance() = 0.f;
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
    const typename detector_t::geometry_context gctx,
    const scalar first_half_length, const scalar second_half_length,
    const scalar strip_length_gap_tolerance) {

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

    scalar corrected_first_half_length = first_half_length;
    scalar corrected_second_half_length = second_half_length;
    if (corrected_first_half_length <= 0.f) {
        const auto first_half_y_values = first_surface.boundary(1u);
        corrected_first_half_length = !first_half_y_values.empty()
                                          ? first_half_y_values.front()
                                          : 0.f;
    }
    if (corrected_second_half_length <= 0.f) {
        const auto second_half_y_values = second_surface.boundary(1u);
        corrected_second_half_length = !second_half_y_values.empty()
                                           ? second_half_y_values.front()
                                           : 0.f;
    }

    fill_corrected_strip_spacepoint_from_lines(
        sp, first_center, second_center, first_direction, second_direction,
        corrected_first_half_length, corrected_second_half_length,
        strip_length_gap_tolerance);
}

template <typename spacepoint_backend_t, typename detector_t,
          typename measurement_backend_t>
TRACCC_HOST_DEVICE inline void fill_endcap_strip_spacepoint(
    edm::spacepoint<spacepoint_backend_t>& sp, const detector_t& det,
    const edm::measurement<measurement_backend_t>& first_meas,
    const edm::measurement<measurement_backend_t>& second_meas,
    const typename detector_t::geometry_context gctx,
    const scalar first_mid_r_from_pair, const scalar second_mid_r_from_pair,
    const scalar first_half_length, const scalar second_half_length,
    const scalar strip_length_gap_tolerance) {

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

    scalar first_mid_r = first_mid_r_from_pair;
    scalar second_mid_r = second_mid_r_from_pair;
    if (first_mid_r <= 0.f) {
        first_mid_r = 0.5f * (first_min_r + first_max_r);
    }
    if (second_mid_r <= 0.f) {
        second_mid_r = 0.5f * (second_min_r + second_max_r);
    }

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

    // Apply the same G80-like endpoint correction as for the barrel. The
    // endcap half lengths are prepared on the host in the focal local0 frame,
    // after converting the beam-frame R boundaries for this measurement's
    // local1 value.
    fill_corrected_strip_spacepoint_from_lines(
        sp, first_center, second_center, first_direction, second_direction,
        first_half_length, second_half_length, strip_length_gap_tolerance);
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
