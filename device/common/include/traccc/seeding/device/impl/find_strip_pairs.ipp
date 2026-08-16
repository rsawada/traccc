/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <vecmem/memory/device_atomic_ref.hpp>

namespace traccc::device {

template <typename detector_t>
TRACCC_HOST_DEVICE inline void find_strip_pairs(
    const global_index_t globalIndex, typename detector_t::view det_view,
    const edm::measurement_collection<default_algebra>::const_view&
        measurements_view,
    const strip_measurement_surface_info_collection_types::const_view&
        surface_infos_view,
    const point3& beam_spot, unsigned int& opposite_position,
    unsigned int& overlap_position,
    strip_pair_collection_types::view opposite_pairs_view,
    strip_pair_collection_types::view overlap_pairs_view) {

    const edm::measurement_collection<default_algebra>::const_device
        measurements(measurements_view);
    const strip_measurement_surface_info_collection_types::const_device
        surface_infos(surface_infos_view);
    if (globalIndex >= measurements.size()) {
        return;
    }

    vecmem::device_vector<strip_pair> opposite_pairs(opposite_pairs_view);
    vecmem::device_vector<strip_pair> overlap_pairs(overlap_pairs_view);
    const edm::measurement reference_measurement =
        measurements.at(globalIndex);
    strip_measurement_surface_info reference_info{};
    if ((reference_measurement.dimensions() != 1u) ||
        !details::find_strip_surface_info(
            surface_infos, reference_measurement.surface_link().value(),
            reference_info) ||
        (reference_info.is_reference_surface == 0u)) {
        return;
    }

    (void)det_view;
    const details::strip_material reference_material =
        details::make_strip_material(reference_measurement, reference_info,
                                     beam_spot);
    if (reference_material.valid == 0u) {
        return;
    }

    const std::uint64_t related_links[5]{
        reference_info.opposite_surface_link,
        reference_info.eta_minus_surface_link,
        reference_info.eta_plus_surface_link,
        reference_info.phi_minus_surface_link,
        reference_info.phi_plus_surface_link};
    for (unsigned int relation_index = 0u; relation_index < 5u;
         ++relation_index) {
        const std::uint64_t candidate_link = related_links[relation_index];
        if (candidate_link == invalid_strip_surface_link) {
            continue;
        }
        bool duplicate = false;
        for (unsigned int previous = 0u; previous < relation_index;
             ++previous) {
            duplicate |= related_links[previous] == candidate_link;
        }
        if (duplicate) {
            continue;
        }

        strip_measurement_surface_info candidate_info{};
        if (!details::find_strip_surface_info(surface_infos, candidate_link,
                                              candidate_info)) {
            continue;
        }
        const unsigned int candidate_begin =
            details::lower_measurement_bound(measurements, candidate_link);
        const unsigned int candidate_end =
            details::upper_measurement_bound(measurements, candidate_link);
        for (unsigned int candidate_index = candidate_begin;
             candidate_index < candidate_end; ++candidate_index) {
            if (candidate_index == globalIndex) {
                continue;
            }
            const edm::measurement candidate_measurement =
                measurements.at(candidate_index);
            const auto relation = details::match_offline_strip_pair(
                reference_measurement, candidate_measurement, reference_info,
                candidate_info);
            if (relation == strip_pair_relation::none) {
                continue;
            }

            const details::strip_material candidate_material =
                details::make_strip_material(candidate_measurement,
                                             candidate_info, beam_spot);
            if (candidate_material.valid == 0u) {
                continue;
            }
            const scalar reference_r = vector::perp(reference_material.center);
            const scalar candidate_r = vector::perp(candidate_material.center);
            const scalar dx =
                reference_material.center[0] - candidate_material.center[0];
            const scalar dy =
                reference_material.center[1] - candidate_material.center[1];
            const strip_pair pair{
                static_cast<unsigned int>(globalIndex),
                candidate_index,
                reference_measurement.surface_link().value(),
                candidate_measurement.surface_link().value(),
                candidate_r - reference_r,
                std::sqrt(dx * dx + dy * dy),
                reference_material.center[2] -
                    candidate_material.center[2],
                vector::dot(reference_material.normal,
                            candidate_material.normal),
                reference_info.is_endcap,
                static_cast<scalar>(reference_info.mid_r),
                static_cast<scalar>(candidate_info.mid_r),
                reference_material.half_length,
                candidate_material.half_length,
                details::g80_strip_length_gap_tolerance(
                    reference_info, candidate_info)};

            if (relation == strip_pair_relation::opposite) {
                const unsigned int position =
                    vecmem::device_atomic_ref<unsigned int>(
                        opposite_position)
                        .fetch_add(1u);
                if (position < opposite_pairs.size()) {
                    opposite_pairs.at(position) = pair;
                }
            } else {
                const unsigned int position =
                    vecmem::device_atomic_ref<unsigned int>(
                        overlap_position)
                        .fetch_add(1u);
                if (position < overlap_pairs.size()) {
                    overlap_pairs.at(position) = pair;
                }
            }
        }
    }
}

}  // namespace traccc::device
