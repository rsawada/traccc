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

/// Relationship between the reference and candidate strip surfaces.
enum class strip_pair_relation : unsigned int {
    none = 0u,
    opposite = 1u,
    eta_minus = 2u,
    eta_plus = 3u,
    phi_minus = 4u,
    phi_plus = 5u
};

inline constexpr std::uint64_t invalid_strip_surface_link = UINT64_MAX;

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
    /// Active local coordinate used by the offline strip pairing.
    scalar active_local;
    /// Strip index used for ITk endcap phi-overlap matching.
    scalar strip_index;
    /// Signed barrel/endcap identifier (-2, 0, or +2).
    int barrel_ec;
    /// Whether this surface is the non-stereo reference side.
    unsigned int is_reference_surface;
    /// Detray surface links of the offline module relationships.
    std::uint64_t opposite_surface_link;
    std::uint64_t eta_minus_surface_link;
    std::uint64_t eta_plus_surface_link;
    std::uint64_t phi_minus_surface_link;
    std::uint64_t phi_plus_surface_link;
    /// Final offline compatibility ranges for opposite and eta neighbours.
    scalar opposite_min;
    scalar opposite_max;
    scalar eta_minus_min;
    scalar eta_minus_max;
    scalar eta_plus_min;
    scalar eta_plus_max;
    /// Final offline edge ranges for the reference and phi neighbour.
    scalar phi_minus_reference_min;
    scalar phi_minus_reference_max;
    scalar phi_minus_candidate_min;
    scalar phi_minus_candidate_max;
    scalar phi_plus_reference_min;
    scalar phi_plus_reference_max;
    scalar phi_plus_candidate_min;
    scalar phi_plus_candidate_max;
    /// Range in the flattened list of measurements on related surfaces.
    unsigned int candidate_measurement_begin;
    unsigned int candidate_measurement_count;
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

/// Return the offline relationship if two strip measurements are compatible.
template <typename measurement_backend_t>
TRACCC_HOST_DEVICE inline strip_pair_relation match_offline_strip_pair(
    const edm::measurement<measurement_backend_t>& reference_measurement,
    const edm::measurement<measurement_backend_t>& candidate_measurement,
    const strip_measurement_surface_info& reference_info,
    const strip_measurement_surface_info& candidate_info) {

    if ((reference_measurement.dimensions() != 1u) ||
        (candidate_measurement.dimensions() != 1u) ||
        (reference_info.is_reference_surface == 0u) ||
        (reference_info.is_endcap != candidate_info.is_endcap) ||
        ((reference_info.has_barrel_material == 0u) &&
         (reference_info.has_endcap_material == 0u)) ||
        ((candidate_info.has_barrel_material == 0u) &&
         (candidate_info.has_endcap_material == 0u))) {
        return strip_pair_relation::none;
    }

    const std::uint64_t candidate_link =
        candidate_measurement.surface_link().value();
    scalar min_value = 0.f;
    scalar max_value = 0.f;
    strip_pair_relation relation = strip_pair_relation::none;

    if (candidate_link == reference_info.opposite_surface_link) {
        relation = strip_pair_relation::opposite;
        min_value = reference_info.opposite_min;
        max_value = reference_info.opposite_max;
    } else if (candidate_link == reference_info.eta_minus_surface_link) {
        relation = strip_pair_relation::eta_minus;
        min_value = reference_info.eta_minus_min;
        max_value = reference_info.eta_minus_max;
    } else if (candidate_link == reference_info.eta_plus_surface_link) {
        relation = strip_pair_relation::eta_plus;
        min_value = reference_info.eta_plus_min;
        max_value = reference_info.eta_plus_max;
    } else if (candidate_link == reference_info.phi_minus_surface_link) {
        const scalar reference_value = reference_info.is_endcap != 0u
                                           ? reference_info.strip_index
                                           : reference_info.active_local;
        const scalar candidate_value = candidate_info.is_endcap != 0u
                                           ? candidate_info.strip_index
                                           : candidate_info.active_local;
        if ((reference_value < reference_info.phi_minus_reference_min) ||
            (reference_value > reference_info.phi_minus_reference_max) ||
            (candidate_value < reference_info.phi_minus_candidate_min) ||
            (candidate_value > reference_info.phi_minus_candidate_max)) {
            return strip_pair_relation::none;
        }
        return strip_pair_relation::phi_minus;
    } else if (candidate_link == reference_info.phi_plus_surface_link) {
        const scalar reference_value = reference_info.is_endcap != 0u
                                           ? reference_info.strip_index
                                           : reference_info.active_local;
        const scalar candidate_value = candidate_info.is_endcap != 0u
                                           ? candidate_info.strip_index
                                           : candidate_info.active_local;
        if ((reference_value < reference_info.phi_plus_reference_min) ||
            (reference_value > reference_info.phi_plus_reference_max) ||
            (candidate_value < reference_info.phi_plus_candidate_min) ||
            (candidate_value > reference_info.phi_plus_candidate_max)) {
            return strip_pair_relation::none;
        }
        return strip_pair_relation::phi_plus;
    } else {
        return strip_pair_relation::none;
    }

    scalar difference = candidate_info.active_local - reference_info.active_local;
    if (candidate_info.barrel_ec < 0) {
        difference = -difference;
    }
    return ((difference >= min_value) && (difference <= max_value))
               ? relation
               : strip_pair_relation::none;
}

TRACCC_HOST_DEVICE inline bool is_overlap_relation(
    const strip_pair_relation relation) {
    return (relation == strip_pair_relation::eta_minus) ||
           (relation == strip_pair_relation::eta_plus) ||
           (relation == strip_pair_relation::phi_minus) ||
           (relation == strip_pair_relation::phi_plus);
}

}  // namespace details
}  // namespace traccc
