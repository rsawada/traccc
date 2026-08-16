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
#include <detray/geometry/barcode.hpp>

// System include(s).
#include <array>
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
inline constexpr unsigned int max_strip_rows = 4u;

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

/// Static strip surface information prepared once on the host.
struct strip_measurement_surface_info {
    /// Surface link identifying this record. Records are sorted by this value.
    std::uint64_t surface_link;
    /// Whether this entry corresponds to an endcap strip surface.
    unsigned int is_endcap;
    /// Minimum radius from the endcap annulus boundary.
    double min_r;
    /// Maximum radius from the endcap annulus boundary.
    double max_r;
    /// Radius used as a pseudo strip midpoint for pair diagnostics.
    double mid_r;
    /// Whether the Athena detector element uses an annulus design.
    unsigned int is_annulus;
    /// Athena surface frame used by SiDetectorElement::globalPosition().
    std::array<double, 3u> athena_origin;
    std::array<double, 3u> athena_eta_axis;
    std::array<double, 3u> athena_phi_axis;
    std::array<double, 3u> athena_normal;
    /// Exact inputs used by StripSpacePointFormationTool::offset().
    std::array<double, 3u> athena_transform_local_x;
    std::array<double, 3u> athena_transform_normal;
    std::array<double, 3u> athena_transform_translation;
    std::array<double, 3u> athena_element_center;
    /// StripBoxDesign parameters used by cellIdOfPosition()/endsOfStrip().
    unsigned int barrel_n_rows;
    unsigned int barrel_n_strips;
    double barrel_pitch;
    double barrel_row_length;
    /// StripStereoAnnulusDesign parameters, stored for every strip row.
    unsigned int endcap_n_rows;
    std::array<unsigned int, max_strip_rows> endcap_n_strips;
    std::array<unsigned int, max_strip_rows> endcap_first_strip;
    std::array<double, max_strip_rows> endcap_phi_pitch;
    std::array<double, max_strip_rows> endcap_strip_start_radius;
    std::array<double, max_strip_rows> endcap_strip_end_radius;
    double endcap_stereo_angle;
    double endcap_wafer_center_r;
    double endcap_center_r;
    /// Radius used by the offline phi-overlap cellIdOfPosition() lookup.
    double endcap_overlap_lookup_r;
    /// Cached constants used verbatim by StripStereoAnnulusDesign.
    double endcap_phi_width;
    double endcap_length_bf;
    double endcap_sin_stereo;
    double endcap_cos_stereo;
    /// Whether the endcap surface uses polar local coordinates.
    unsigned int endcap_uses_polar_coordinates;
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

/// Geometry reconstructed for one measured strip on the device.
struct strip_material {
    point3 center{};
    vector3 direction{};
    vector3 trajectory{};
    vector3 normal{};
    scalar half_length{0.f};
    unsigned int valid{0u};
};

template <typename surface_info_collection_t>
TRACCC_HOST_DEVICE inline bool find_strip_surface_info(
    const surface_info_collection_t& surface_infos,
    const std::uint64_t surface_link,
    strip_measurement_surface_info& result) {
    unsigned int first = 0u;
    unsigned int last = static_cast<unsigned int>(surface_infos.size());
    while (first < last) {
        const unsigned int middle = first + (last - first) / 2u;
        const auto candidate = surface_infos.at(middle);
        if (candidate.surface_link < surface_link) {
            first = middle + 1u;
        } else {
            last = middle;
        }
    }
    if ((first >= surface_infos.size()) ||
        (surface_infos.at(first).surface_link != surface_link)) {
        return false;
    }
    result = surface_infos.at(first);
    return true;
}

template <typename measurement_collection_t>
TRACCC_HOST_DEVICE inline unsigned int lower_measurement_bound(
    const measurement_collection_t& measurements,
    const std::uint64_t surface_link) {
    // Measurement sorting compares detray barcodes by surface index, not by
    // their full encoded value. Use the same key for this binary search.
    const auto surface_index =
        detray::geometry::barcode{surface_link}.index();
    unsigned int first = 0u;
    unsigned int last = static_cast<unsigned int>(measurements.size());
    while (first < last) {
        const unsigned int middle = first + (last - first) / 2u;
        if (measurements.at(middle).surface_link().index() < surface_index) {
            first = middle + 1u;
        } else {
            last = middle;
        }
    }
    return first;
}

template <typename measurement_collection_t>
TRACCC_HOST_DEVICE inline unsigned int upper_measurement_bound(
    const measurement_collection_t& measurements,
    const std::uint64_t surface_link) {
    const auto surface_index =
        detray::geometry::barcode{surface_link}.index();
    unsigned int first = 0u;
    unsigned int last = static_cast<unsigned int>(measurements.size());
    while (first < last) {
        const unsigned int middle = first + (last - first) / 2u;
        if (measurements.at(middle).surface_link().index() <= surface_index) {
            first = middle + 1u;
        } else {
            last = middle;
        }
    }
    return first;
}

struct endcap_cell_id {
    int strip{-1};
    int row{-1};
    int strip_1d{-1};
    unsigned int valid{0u};
};

struct exact_point2 {
    double x_eta{0.};
    double x_phi{0.};
};

struct exact_point3 {
    double x{0.};
    double y{0.};
    double z{0.};
};

TRACCC_HOST_DEVICE inline exact_point3 athena_global_position(
    const strip_measurement_surface_info& info, const double x_eta,
    const double x_phi) {
    // SolidStateDetectorElementBase::globalPosition(Vector2D):
    // origin + local[distEta] * etaAxis + local[distPhi] * phiAxis.
    return {info.athena_origin[0] + x_eta * info.athena_eta_axis[0] +
                x_phi * info.athena_phi_axis[0],
            info.athena_origin[1] + x_eta * info.athena_eta_axis[1] +
                x_phi * info.athena_phi_axis[1],
            info.athena_origin[2] + x_eta * info.athena_eta_axis[2] +
                x_phi * info.athena_phi_axis[2]};
}

struct local_strip_endpoints {
    exact_point2 first{};
    exact_point2 second{};
    unsigned int valid{0u};
};

/// Device counterpart of StripBoxDesign::cellIdOfPosition()/endsOfStrip().
TRACCC_HOST_DEVICE inline local_strip_endpoints barrel_ends_of_strip(
    const double x_eta, const double x_phi,
    const strip_measurement_surface_info& info) {
    local_strip_endpoints result{};
    if ((info.barrel_n_rows == 0u) || (info.barrel_n_strips == 0u) ||
        (info.barrel_pitch <= 0.f) || (info.barrel_row_length <= 0.f)) {
        return result;
    }

    const int strip = static_cast<int>(
        std::floor(x_phi / info.barrel_pitch) +
        static_cast<double>(info.barrel_n_strips) * 0.5);
    int row = 0;
    if (info.barrel_n_rows > 1u) {
        row = static_cast<int>(
            std::floor(x_eta / info.barrel_row_length) +
            static_cast<double>(info.barrel_n_rows) * 0.5);
    }
    if ((strip < 0) ||
        (strip >= static_cast<int>(info.barrel_n_strips)) || (row < 0) ||
        (row >= static_cast<int>(info.barrel_n_rows))) {
        return result;
    }

    const double eta_start =
        (static_cast<double>(row) -
         static_cast<double>(info.barrel_n_rows) * 0.5) *
        info.barrel_row_length;
    const double eta_end = eta_start + info.barrel_row_length;
    const double phi =
        (static_cast<double>(strip) -
         static_cast<double>(info.barrel_n_strips) * 0.5 + 0.5) *
        info.barrel_pitch;
    result.first = {eta_start, phi};
    result.second = {eta_end, phi};
    result.valid = 1u;
    return result;
}

/// Exact device counterpart of StripStereoAnnulusDesign::beamToStrip().
TRACCC_HOST_DEVICE inline exact_point2 endcap_beam_to_strip(
    const exact_point2& position, const strip_measurement_surface_info& info) {
    double x_beam = position.x_eta;
    double y_beam = position.x_phi;
    if (info.endcap_uses_polar_coordinates != 0u) {
        x_beam = position.x_eta * std::cos(position.x_phi);
        y_beam = position.x_eta * std::sin(position.x_phi);
    }
    const double cos_stereo = info.endcap_cos_stereo;
    const double sin_stereo = info.endcap_sin_stereo;
    const double x_strip =
        cos_stereo * (x_beam - info.endcap_wafer_center_r) +
        sin_stereo * y_beam + info.endcap_wafer_center_r;
    const double y_strip =
        -sin_stereo * (x_beam - info.endcap_wafer_center_r) +
        cos_stereo * y_beam;
    if (info.endcap_uses_polar_coordinates != 0u) {
        return {std::sqrt(x_strip * x_strip + y_strip * y_strip),
                std::atan2(y_strip, x_strip)};
    }
    return {x_strip, y_strip};
}

/// Exact device counterpart of StripStereoAnnulusDesign::cellIdOfPosition().
template <typename measurement_backend_t>
TRACCC_HOST_DEVICE inline scalar strip_active_local(
    const edm::measurement<measurement_backend_t>& measurement,
    const strip_measurement_surface_info& info) {
    return info.is_endcap != 0u ? measurement.local_position()[1]
                                : measurement.local_position()[0];
}

TRACCC_HOST_DEVICE inline endcap_cell_id endcap_cell_id_of_position(
    const exact_point2& position, const strip_measurement_surface_info& info) {
    endcap_cell_id result{};
    if (info.endcap_n_rows == 0u) {
        return result;
    }
    const double radius = info.endcap_uses_polar_coordinates != 0u
                              ? position.x_eta
                              : std::sqrt(position.x_eta * position.x_eta +
                                          position.x_phi * position.x_phi);
    if ((radius < info.endcap_strip_start_radius[0]) ||
        (radius >=
         info.endcap_strip_end_radius[info.endcap_n_rows - 1u])) {
        return result;
    }
    unsigned int row = 0u;
    for (unsigned int candidate = 1u; candidate < info.endcap_n_rows;
         ++candidate) {
        if (radius >= info.endcap_strip_start_radius[candidate]) {
            row = candidate;
        }
    }
    const exact_point2 strip_position = endcap_beam_to_strip(position, info);
    const double phi_strip = info.endcap_uses_polar_coordinates != 0u
                                 ? strip_position.x_phi
                                 : std::atan2(strip_position.x_phi,
                                              strip_position.x_eta);
    const int strip = static_cast<int>(
        std::floor(phi_strip / info.endcap_phi_pitch[row]) +
        static_cast<double>(info.endcap_n_strips[row]) * 0.5);
    if ((strip < 0) ||
        (strip >= static_cast<int>(info.endcap_n_strips[row]))) {
        return result;
    }
    result.strip = strip;
    result.row = static_cast<int>(row);
    result.strip_1d =
        static_cast<int>(info.endcap_first_strip[row]) + strip;
    result.valid = 1u;
    return result;
}

template <typename measurement_backend_t>
TRACCC_HOST_DEVICE inline scalar endcap_strip_index(
    const edm::measurement<measurement_backend_t>& measurement,
    const strip_measurement_surface_info& info) {
    if ((info.endcap_n_rows == 0u) ||
        (info.endcap_n_strips[0] == 0u) ||
        (info.endcap_phi_width == 0.)) {
        return -1.f;
    }

    // Match StripSpacePointFormationTool::getStripEnds(). The offline
    // overlap code compares this strip index with limits obtained separately
    // through correctPolarRange()/cellIdOfPosition().
    const double phi_pitch_phi =
        info.endcap_phi_width /
        static_cast<double>(info.endcap_n_strips[0]);
    const double local_phi =
        static_cast<double>(measurement.local_position()[1]);
    const double strip_index_value =
        -std::floor(local_phi / phi_pitch_phi) +
        static_cast<double>(info.endcap_n_strips[0]) * 0.5 - 0.5;
    if ((strip_index_value < 0.) ||
        (strip_index_value >=
         static_cast<double>(info.endcap_n_strips[0]))) {
        return -1.f;
    }
    // Athena stores this value in std::size_t before using it, which
    // truncates the half-integer produced for rows with an even strip count.
    const unsigned int strip_index =
        static_cast<unsigned int>(strip_index_value);
    return static_cast<scalar>(strip_index);
}

/// Device counterpart of StripStereoAnnulusDesign::stripPosAtR().
TRACCC_HOST_DEVICE inline exact_point2 endcap_strip_position_at_r(
    const int strip_index, const unsigned int row, const double radius,
    const strip_measurement_surface_info& info) {
    const double phi_strip =
        (static_cast<double>(strip_index) -
         static_cast<double>(info.endcap_n_strips[row]) * 0.5 + 0.5) *
        info.endcap_phi_pitch[row];
    const double b = -2. * info.endcap_length_bf *
                     std::sin(0.5 * info.endcap_stereo_angle + phi_strip);
    const double c = info.endcap_length_bf * info.endcap_length_bf -
                     radius * radius;
    const double strip_r = 0.5 * (-b + std::sqrt(b * b - 4. * c));
    const double strip_x = strip_r * std::cos(phi_strip);
    const double strip_y = strip_r * std::sin(phi_strip);
    const double cos_stereo = info.endcap_cos_stereo;
    const double sin_stereo = info.endcap_sin_stereo;
    const double beam_x =
        cos_stereo * (strip_x - info.endcap_wafer_center_r) -
        sin_stereo * strip_y + info.endcap_wafer_center_r;
    const double beam_y =
        sin_stereo * (strip_x - info.endcap_wafer_center_r) +
        cos_stereo * strip_y;
    if (info.endcap_uses_polar_coordinates != 0u) {
        return {std::sqrt(beam_x * beam_x + beam_y * beam_y),
                std::atan2(beam_y, beam_x)};
    }
    return {beam_x, beam_y};
}

/// Device counterpart of the row-0 strip-index calculation used by the
/// validated Athena CPU prototype.
TRACCC_HOST_DEVICE inline int endcap_row0_strip_index(
    const double active_local, const strip_measurement_surface_info& info) {
    if ((info.endcap_n_rows == 0u) || (info.endcap_n_strips[0] == 0u) ||
        (info.endcap_phi_width <= 0.)) {
        return -1;
    }
    const double phi_pitch =
        info.endcap_phi_width /
        static_cast<double>(info.endcap_n_strips[0]);
    const double value =
        -std::floor(active_local / phi_pitch) +
        static_cast<double>(info.endcap_n_strips[0]) * 0.5 - 0.5;
    return ((value >= 0.) &&
            (value < static_cast<double>(info.endcap_n_strips[0])))
               ? static_cast<int>(value)
               : -1;
}

template <typename measurement_backend_t>
TRACCC_HOST_DEVICE inline strip_material make_strip_material(
    const edm::measurement<measurement_backend_t>& measurement,
    const strip_measurement_surface_info& info,
    const point3& beam_spot) {
    strip_material result{};
    exact_point3 first{};
    exact_point3 second{};
    if (info.is_endcap != 0u) {
        // This is the exact member-function chain used by the validated CPU
        // prototype: row-0 strip selection, stripPosAtR(minR/maxR), then
        // SolidStateDetectorElementBase::globalPosition(Vector2D).
        const int strip_index = endcap_row0_strip_index(
            static_cast<double>(measurement.local_position()[1]), info);
        if (strip_index < 0) {
            return result;
        }
        const exact_point2 first_local = endcap_strip_position_at_r(
            strip_index, 0u, info.min_r, info);
        const exact_point2 second_local = endcap_strip_position_at_r(
            strip_index, 0u, info.max_r, info);
        first = athena_global_position(info, first_local.x_eta,
                                       first_local.x_phi);
        second = athena_global_position(info, second_local.x_eta,
                                        second_local.x_phi);
    } else {
        // This exactly mirrors endsOfStrip(SiLocalPosition(0, local0, 0))
        // followed by SiDetectorElement::globalPosition().
        const local_strip_endpoints local_ends = barrel_ends_of_strip(
            0., static_cast<double>(measurement.local_position()[0]), info);
        if (local_ends.valid == 0u) {
            return result;
        }
        first = athena_global_position(info, local_ends.first.x_eta,
                                       local_ends.first.x_phi);
        second = athena_global_position(info, local_ends.second.x_eta,
                                        local_ends.second.x_phi);
    }

    // Athena computes these quantities in double precision. Convert only the
    // completed material to traccc scalar, matching the old CPU path.
    const double center_x = 0.5 * (first.x + second.x);
    const double center_y = 0.5 * (first.y + second.y);
    const double center_z = 0.5 * (first.z + second.z);
    const double direction_x = first.x - second.x;
    const double direction_y = first.y - second.y;
    const double direction_z = first.z - second.z;
    const double trajectory_x =
        2. * (center_x - static_cast<double>(beam_spot[0]));
    const double trajectory_y =
        2. * (center_y - static_cast<double>(beam_spot[1]));
    const double trajectory_z =
        2. * (center_z - static_cast<double>(beam_spot[2]));
    const double normal_x =
        direction_y * trajectory_z - direction_z * trajectory_y;
    const double normal_y =
        direction_z * trajectory_x - direction_x * trajectory_z;
    const double normal_z =
        direction_x * trajectory_y - direction_y * trajectory_x;
    const double length = std::sqrt(direction_x * direction_x +
                                    direction_y * direction_y +
                                    direction_z * direction_z);

    result.center = {static_cast<scalar>(center_x),
                     static_cast<scalar>(center_y),
                     static_cast<scalar>(center_z)};
    result.direction = {static_cast<scalar>(direction_x),
                        static_cast<scalar>(direction_y),
                        static_cast<scalar>(direction_z)};
    result.trajectory = {static_cast<scalar>(trajectory_x),
                         static_cast<scalar>(trajectory_y),
                         static_cast<scalar>(trajectory_z)};
    result.normal = {static_cast<scalar>(normal_x),
                     static_cast<scalar>(normal_y),
                     static_cast<scalar>(normal_z)};
    result.half_length = static_cast<scalar>(0.5 * length);
    result.valid = length > 0. ? 1u : 0u;
    return result;
}

/// Reproduce the ACTS CPU StripSpacePointFormationTool::offset calculation
/// using the exact transform columns and element center cached on the host.
TRACCC_HOST_DEVICE inline scalar g80_strip_length_gap_tolerance(
    const strip_measurement_surface_info& first_info,
    const strip_measurement_surface_info& second_info,
    const scalar strip_gap_parameter = 0.0015f) {

    if (strip_gap_parameter == 0.f) {
        return 0.f;
    }

    const double x12 =
        first_info.athena_transform_local_x[0] *
            second_info.athena_transform_local_x[0] +
        first_info.athena_transform_local_x[1] *
            second_info.athena_transform_local_x[1] +
        first_info.athena_transform_local_x[2] *
            second_info.athena_transform_local_x[2];
    const double dx = first_info.athena_transform_translation[0] -
                      second_info.athena_transform_translation[0];
    const double dy = first_info.athena_transform_translation[1] -
                      second_info.athena_transform_translation[1];
    const double dz = first_info.athena_transform_translation[2] -
                      second_info.athena_transform_translation[2];
    const double surface_separation =
        dx * first_info.athena_transform_normal[0] +
        dy * first_info.athena_transform_normal[1] +
        dz * first_info.athena_transform_normal[2];
    const double radius_x = first_info.is_annulus != 0u
                                ? first_info.athena_element_center[0]
                                : first_info.athena_transform_translation[0];
    const double radius_y = first_info.is_annulus != 0u
                                ? first_info.athena_element_center[1]
                                : first_info.athena_transform_translation[1];
    const double surface_reference_r =
        std::sqrt(radius_x * radius_x + radius_y * radius_y);
    const double dm = static_cast<double>(strip_gap_parameter) *
                      surface_reference_r *
                      std::abs(surface_separation * x12);

    double tolerance = 0.;
    if (first_info.is_annulus != 0u) {
        tolerance = dm / 0.04;
    } else {
        const double denominator2 = (1. - x12) * (1. + x12);
        if (denominator2 > 0.) {
            tolerance = dm / std::sqrt(denominator2);
        }
    }

    if ((std::abs(first_info.athena_transform_normal[2]) > 0.7) &&
        (std::abs(first_info.athena_transform_translation[2]) > 0.)) {
        tolerance *= surface_reference_r /
                     std::abs(first_info.athena_transform_translation[2]);
    }
    return static_cast<scalar>(tolerance);
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
        (reference_info.is_endcap != candidate_info.is_endcap)) {
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
        const scalar reference_value =
            reference_info.is_endcap != 0u
                ? endcap_strip_index(reference_measurement, reference_info)
                : strip_active_local(reference_measurement, reference_info);
        const scalar candidate_value =
            candidate_info.is_endcap != 0u
                ? endcap_strip_index(candidate_measurement, candidate_info)
                : strip_active_local(candidate_measurement, candidate_info);
        if ((reference_value < reference_info.phi_minus_reference_min) ||
            (reference_value > reference_info.phi_minus_reference_max) ||
            (candidate_value < reference_info.phi_minus_candidate_min) ||
            (candidate_value > reference_info.phi_minus_candidate_max)) {
            return strip_pair_relation::none;
        }
        return strip_pair_relation::phi_minus;
    } else if (candidate_link == reference_info.phi_plus_surface_link) {
        const scalar reference_value =
            reference_info.is_endcap != 0u
                ? endcap_strip_index(reference_measurement, reference_info)
                : strip_active_local(reference_measurement, reference_info);
        const scalar candidate_value =
            candidate_info.is_endcap != 0u
                ? endcap_strip_index(candidate_measurement, candidate_info)
                : strip_active_local(candidate_measurement, candidate_info);
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

    scalar difference = strip_active_local(candidate_measurement, candidate_info) -
                        strip_active_local(reference_measurement, reference_info);
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
