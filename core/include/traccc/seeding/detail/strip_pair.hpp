/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/container.hpp"
#include "traccc/edm/measurement_collection.hpp"

// Detray include(s).
#include <detray/geometry/tracking_surface.hpp>

// System include(s).
#include <cmath>
#include <cstdint>

namespace traccc {

/// Pair of strip measurements used to form one strip spacepoint.
struct strip_pair {
    /// Index of the measurement on the inner surface.
    unsigned int measurement_index_1;
    /// Index of the measurement on the outer surface.
    unsigned int measurement_index_2;
    /// Surface link of the measurement on the inner surface.
    std::uint64_t surface_link_1;
    /// Surface link of the measurement on the outer surface.
    std::uint64_t surface_link_2;
    /// Difference between the outer and inner surface radii.
    scalar surface_delta_r;
    /// Distance between the selected strip centers in the xy plane.
    scalar strip_center_delta_xy;
    /// Signed difference between the selected strip center z coordinates.
    scalar strip_center_delta_z;
    /// Dot product of the surface normal vectors.
    scalar normal_dot;
    /// Whether this pair uses endcap strip surface information.
    unsigned int is_endcap;
    /// Host-provided midpoint radius of the first strip surface.
    scalar surface_mid_r_1;
    /// Host-provided midpoint radius of the second strip surface.
    scalar surface_mid_r_2;
    /// Half length of the first strip direction used for the fill.
    scalar strip_half_length_1;
    /// Half length of the second strip direction used for the fill.
    scalar strip_half_length_2;
    /// G80-like strip-length gap tolerance used for m/n correction.
    scalar strip_length_gap_tolerance;
};

/// Declare all strip pair collection types.
using strip_pair_collection_types = collection_types<strip_pair>;

/// Per-measurement strip surface information prepared on the host.
struct strip_measurement_surface_info {
    /// Surface link associated with the measurement at the same index.
    std::uint64_t surface_link;
    /// Whether this entry corresponds to an endcap strip surface.
    unsigned int is_endcap;
    /// Minimum radius from the endcap annulus boundary.
    scalar min_r;
    /// Maximum radius from the endcap annulus boundary.
    scalar max_r;
    /// Radius used as a pseudo strip midpoint for pair diagnostics.
    scalar mid_r;
    /// Strip half length prepared on the host for device-side filling.
    scalar strip_half_length;
    /// G80-compatible center of the measured barrel strip.
    point3 barrel_strip_center;
    /// G80-compatible, unnormalised barrel strip direction (start - end).
    vector3 barrel_strip_direction;
    /// Direction from the beam spot to the strip center, multiplied by two.
    vector3 barrel_trajectory_direction;
    /// G80-compatible barrel plane normal, built from the strip and beam spot.
    vector3 barrel_strip_normal;
    /// Whether the G80-compatible barrel material is available.
    unsigned int has_barrel_material;
};

/// Declare all strip measurement surface information collection types.
using strip_measurement_surface_info_collection_types =
    collection_types<strip_measurement_surface_info>;

/// Configuration for the initial barrel strip pair search.
struct barrel_strip_pair_config {
    scalar max_surface_delta_z = 10.f;
    scalar min_surface_delta_r = 3.f;
    scalar max_surface_delta_r = 20.f;
    scalar max_strip_center_delta_xy = 20.f;
    scalar max_strip_center_delta_z = 20.f;
    scalar min_normal_dot = 0.995f;
    scalar strip_gap_parameter = 0.0015f;
};

/// Configuration for the initial endcap strip pair search.
struct endcap_strip_pair_config {
    scalar min_surface_delta_z = 3.f;
    scalar max_surface_delta_z = 8.f;
    scalar max_strip_center_delta_xy = 20.f;
    scalar min_normal_dot = 0.995f;
    scalar strip_gap_parameter = 0.0015f;
};

namespace details {

/// Return whether two barrel strip measurements are opposite-pair candidates.
template <typename detector_t, typename measurement_backend_t>
TRACCC_HOST_DEVICE inline bool is_compatible_barrel_strip_pair(
    const detector_t& det,
    const edm::measurement<measurement_backend_t>& inner_measurement,
    const edm::measurement<measurement_backend_t>& outer_measurement,
    const barrel_strip_pair_config& config) {

    if ((inner_measurement.dimensions() != 1u) ||
        (outer_measurement.dimensions() != 1u) ||
        (inner_measurement.surface_link().value() ==
         outer_measurement.surface_link().value())) {
        return false;
    }

    const detray::tracking_surface inner_surface{
        det, inner_measurement.surface_link()};
    const detray::tracking_surface outer_surface{
        det, outer_measurement.surface_link()};

    // The initial implementation handles barrel rectangle surfaces only.
    if ((static_cast<int>(inner_surface.shape_id()) != 0) ||
        (static_cast<int>(outer_surface.shape_id()) != 0)) {
        return false;
    }

    const point3 inner_surface_center = inner_surface.center({});
    const point3 outer_surface_center = outer_surface.center({});
    const scalar surface_delta_z =
        std::abs(outer_surface_center[2] - inner_surface_center[2]);
    if (surface_delta_z >= config.max_surface_delta_z) {
        return false;
    }

    const scalar inner_surface_r =
        std::sqrt(inner_surface_center[0] * inner_surface_center[0] +
                  inner_surface_center[1] * inner_surface_center[1]);
    const scalar outer_surface_r =
        std::sqrt(outer_surface_center[0] * outer_surface_center[0] +
                  outer_surface_center[1] * outer_surface_center[1]);
    const scalar surface_delta_r = outer_surface_r - inner_surface_r;

    // Only search inner-to-outer pairs. This also avoids duplicate pairs.
    if ((surface_delta_r <= config.min_surface_delta_r) ||
        (surface_delta_r >= config.max_surface_delta_r)) {
        return false;
    }

    const vector3 inner_normal =
        inner_surface.normal({}, inner_measurement.local_position());
    const vector3 outer_normal =
        outer_surface.normal({}, outer_measurement.local_position());
    const scalar normal_dot = inner_normal[0] * outer_normal[0] +
                              inner_normal[1] * outer_normal[1] +
                              inner_normal[2] * outer_normal[2];
    if (normal_dot <= config.min_normal_dot) {
        return false;
    }

    const point2 inner_local_center{inner_measurement.local_position()[0],
                                    0.f};
    const point2 outer_local_center{outer_measurement.local_position()[0],
                                    0.f};
    const point3 inner_strip_center =
        inner_surface.local_to_global({}, inner_local_center, {});
    const point3 outer_strip_center =
        outer_surface.local_to_global({}, outer_local_center, {});

    const scalar delta_x = inner_strip_center[0] - outer_strip_center[0];
    const scalar delta_y = inner_strip_center[1] - outer_strip_center[1];
    const scalar delta_z = inner_strip_center[2] - outer_strip_center[2];
    const scalar strip_center_delta_xy =
        std::sqrt(delta_x * delta_x + delta_y * delta_y);

    if ((strip_center_delta_xy > config.max_strip_center_delta_xy) ||
        (std::abs(delta_z) > config.max_strip_center_delta_z)) {
        return false;
    }

    return true;
}


/// Return whether two endcap strip measurements are pair candidates.
template <typename detector_t, typename measurement_backend_t>
TRACCC_HOST_DEVICE inline bool is_compatible_endcap_strip_pair(
    const detector_t& det,
    const edm::measurement<measurement_backend_t>& first_measurement,
    const edm::measurement<measurement_backend_t>& second_measurement,
    const endcap_strip_pair_config& config) {

    if ((first_measurement.dimensions() != 1u) ||
        (second_measurement.dimensions() != 1u) ||
        (first_measurement.surface_link().value() ==
         second_measurement.surface_link().value())) {
        return false;
    }

    const detray::tracking_surface first_surface{det,
                                                 first_measurement.surface_link()};
    const detray::tracking_surface second_surface{
        det, second_measurement.surface_link()};

    // Endcap strip modules are represented by annulus-like surfaces.
    if ((static_cast<int>(first_surface.shape_id()) == 0) ||
        (static_cast<int>(second_surface.shape_id()) == 0)) {
        return false;
    }

    // Endcap strip measurements use the annulus phi-like coordinate.
    if ((first_measurement.subspace()[0] != 1u) ||
        (second_measurement.subspace()[0] != 1u)) {
        return false;
    }

    const point3 first_center = first_surface.center({});
    const point3 second_center = second_surface.center({});
    const scalar surface_delta_z = std::abs(second_center[2]) -
                                   std::abs(first_center[2]);
    if ((surface_delta_z <= config.min_surface_delta_z) ||
        (surface_delta_z >= config.max_surface_delta_z)) {
        return false;
    }

    const vector3 first_normal =
        first_surface.normal({}, first_measurement.local_position());
    const vector3 second_normal =
        second_surface.normal({}, second_measurement.local_position());
    const scalar normal_dot = first_normal[0] * second_normal[0] +
                              first_normal[1] * second_normal[1] +
                              first_normal[2] * second_normal[2];
    if (normal_dot <= config.min_normal_dot) {
        return false;
    }

    // Debug step: do not use boundary() or pseudo-midpoint xy on device.
    // The host-side cutflow still prints these quantities separately.
    return true;
}

/// Return whether two strip measurements are pair candidates.
template <typename detector_t, typename measurement_backend_t>
TRACCC_HOST_DEVICE inline bool is_compatible_strip_pair(
    const detector_t& det,
    const edm::measurement<measurement_backend_t>& first_measurement,
    const edm::measurement<measurement_backend_t>& second_measurement,
    const barrel_strip_pair_config& barrel_config,
    const endcap_strip_pair_config& endcap_config) {

    const detray::tracking_surface first_surface{det,
                                                 first_measurement.surface_link()};
    const detray::tracking_surface second_surface{
        det, second_measurement.surface_link()};

    if ((static_cast<int>(first_surface.shape_id()) == 0) &&
        (static_cast<int>(second_surface.shape_id()) == 0)) {
        return is_compatible_barrel_strip_pair(det, first_measurement,
                                               second_measurement,
                                               barrel_config);
    }

    return is_compatible_endcap_strip_pair(det, first_measurement,
                                           second_measurement, endcap_config);
}

}  // namespace details
}  // namespace traccc
