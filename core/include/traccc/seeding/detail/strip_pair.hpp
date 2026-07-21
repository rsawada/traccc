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
    /// G80-compatible center of the measured endcap strip.
    point3 endcap_strip_center;
    /// G80-compatible, unnormalised endcap strip direction (start - end).
    vector3 endcap_strip_direction;
    /// Direction from the beam spot to the endcap strip center, times two.
    vector3 endcap_trajectory_direction;
    /// G80-compatible endcap plane normal.
    vector3 endcap_strip_normal;
    /// Whether the G80-compatible endcap material is available.
    unsigned int has_endcap_material;
    /// G80/Athena surface-frame origin, T(:,3), in global coordinates.
    point3 surface_origin;
    /// G80/Athena surface-frame local x axis, T(:,0), in global coordinates.
    vector3 surface_local_x;
    /// G80/Athena surface-frame local y axis, T(:,1), in global coordinates.
    vector3 surface_local_y;
    /// G80/Athena surface-frame normal, T(:,2), in global coordinates.
    vector3 surface_normal;
    /// Radius used by the G80 strip-gap calculation.
    scalar surface_reference_r;
    /// Whether the Athena detector element uses an annulus design.
    unsigned int is_annulus;
    /// Whether the complete Athena surface frame is available.
    unsigned int has_surface_frame;
};

/// Declare all strip measurement surface information collection types.
using strip_measurement_surface_info_collection_types =
    collection_types<strip_measurement_surface_info>;

/// Configuration for the initial barrel strip pair search.
struct barrel_strip_pair_config {
    scalar max_strip_center_delta_z = 10.f;
    scalar min_strip_center_delta_r = 3.f;
    scalar max_strip_center_delta_r = 10.f;
    scalar max_strip_center_delta_xy = 10.f;
};

/// Configuration for the initial endcap strip pair search.
struct endcap_strip_pair_config {
    scalar min_strip_center_delta_abs_z = 3.f;
    scalar max_strip_center_delta_abs_z = 8.f;
    scalar max_strip_center_delta_xy = 5.f;
};

namespace details {

/// Reproduce the G80 StripSpacePointFormationTool::offset calculation.
TRACCC_HOST_DEVICE inline scalar g80_strip_length_gap_tolerance(
    const strip_measurement_surface_info& first_info,
    const strip_measurement_surface_info& second_info,
    const scalar strip_gap_parameter = 0.0015f) {

    if ((first_info.has_surface_frame == 0u) ||
        (second_info.has_surface_frame == 0u) ||
        (strip_gap_parameter == 0.f)) {
        return 0.f;
    }

    const scalar x12 =
        first_info.surface_local_x[0] * second_info.surface_local_x[0] +
        first_info.surface_local_x[1] * second_info.surface_local_x[1] +
        first_info.surface_local_x[2] * second_info.surface_local_x[2];
    const scalar delta_origin_x =
        first_info.surface_origin[0] - second_info.surface_origin[0];
    const scalar delta_origin_y =
        first_info.surface_origin[1] - second_info.surface_origin[1];
    const scalar delta_origin_z =
        first_info.surface_origin[2] - second_info.surface_origin[2];
    const scalar surface_separation =
        delta_origin_x * first_info.surface_normal[0] +
        delta_origin_y * first_info.surface_normal[1] +
        delta_origin_z * first_info.surface_normal[2];
    const scalar dm = strip_gap_parameter * first_info.surface_reference_r *
                      std::abs(surface_separation * x12);

    scalar tolerance = 0.f;
    if (first_info.is_annulus != 0u) {
        tolerance = dm / 0.04f;
    } else {
        const scalar denominator2 = (1.f - x12) * (1.f + x12);
        if (denominator2 > 0.f) {
            tolerance = dm / std::sqrt(denominator2);
        }
    }

    if ((std::abs(first_info.surface_normal[2]) > 0.7f) &&
        (std::abs(first_info.surface_origin[2]) > 0.f)) {
        tolerance *= first_info.surface_reference_r /
                     std::abs(first_info.surface_origin[2]);
    }
    return tolerance;
}

/// Return whether two barrel strip measurements are opposite-pair candidates.
template <typename measurement_backend_t>
TRACCC_HOST_DEVICE inline bool is_compatible_barrel_strip_pair(
    const edm::measurement<measurement_backend_t>& inner_measurement,
    const edm::measurement<measurement_backend_t>& outer_measurement,
    const strip_measurement_surface_info& inner_info,
    const strip_measurement_surface_info& outer_info,
    const barrel_strip_pair_config& config) {

    if ((inner_measurement.dimensions() != 1u) ||
        (outer_measurement.dimensions() != 1u) ||
        (inner_info.has_barrel_material == 0u) ||
        (outer_info.has_barrel_material == 0u) ||
        (inner_info.is_endcap != 0u) || (outer_info.is_endcap != 0u)) {
        return false;
    }

    const point3& inner_center = inner_info.barrel_strip_center;
    const point3& outer_center = outer_info.barrel_strip_center;
    const scalar delta_z = outer_center[2] - inner_center[2];
    if (std::abs(delta_z) >= config.max_strip_center_delta_z) {
        return false;
    }

    const scalar inner_r =
        std::sqrt(inner_center[0] * inner_center[0] +
                  inner_center[1] * inner_center[1]);
    const scalar outer_r =
        std::sqrt(outer_center[0] * outer_center[0] +
                  outer_center[1] * outer_center[1]);
    const scalar delta_r = outer_r - inner_r;

    // Only search inner-to-outer pairs. This also avoids duplicate pairs.
    if ((delta_r <= config.min_strip_center_delta_r) ||
        (delta_r >= config.max_strip_center_delta_r)) {
        return false;
    }

    const scalar delta_x = outer_center[0] - inner_center[0];
    const scalar delta_y = outer_center[1] - inner_center[1];
    const scalar delta_xy2 = delta_x * delta_x + delta_y * delta_y;
    const scalar max_delta_xy2 = config.max_strip_center_delta_xy *
                                 config.max_strip_center_delta_xy;
    if (delta_xy2 >= max_delta_xy2) {
        return false;
    }

    return true;
}


/// Return whether two endcap strip measurements are pair candidates.
template <typename measurement_backend_t>
TRACCC_HOST_DEVICE inline bool is_compatible_endcap_strip_pair(
    const edm::measurement<measurement_backend_t>& first_measurement,
    const edm::measurement<measurement_backend_t>& second_measurement,
    const strip_measurement_surface_info& first_info,
    const strip_measurement_surface_info& second_info,
    const endcap_strip_pair_config& config) {

    if ((first_measurement.dimensions() != 1u) ||
        (second_measurement.dimensions() != 1u) ||
        (first_info.is_endcap == 0u) || (second_info.is_endcap == 0u) ||
        (first_info.has_endcap_material == 0u) ||
        (second_info.has_endcap_material == 0u)) {
        return false;
    }

    const point3& first_center = first_info.endcap_strip_center;
    const point3& second_center = second_info.endcap_strip_center;

    // Require both strips to be on the same endcap side.
    if (first_center[2] * second_center[2] <= 0.f) {
        return false;
    }

    // Search from smaller to larger absolute z to avoid duplicate pairs.
    const scalar delta_abs_z =
        std::abs(second_center[2]) - std::abs(first_center[2]);
    if ((delta_abs_z <= config.min_strip_center_delta_abs_z) ||
        (delta_abs_z >= config.max_strip_center_delta_abs_z)) {
        return false;
    }

    const scalar delta_x = second_center[0] - first_center[0];
    const scalar delta_y = second_center[1] - first_center[1];
    const scalar delta_xy2 = delta_x * delta_x + delta_y * delta_y;
    const scalar max_delta_xy2 = config.max_strip_center_delta_xy *
                                 config.max_strip_center_delta_xy;
    if (delta_xy2 >= max_delta_xy2) {
        return false;
    }

    return true;
}

}  // namespace details
}  // namespace traccc
