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
TRACCC_HOST_DEVICE inline void count_strip_pairs(
    const global_index_t globalIndex, typename detector_t::view det_view,
    const edm::measurement_collection<default_algebra>::const_view&
        measurements_view,
    const strip_measurement_surface_info_collection_types::const_view&
        surface_infos_view,
    unsigned int& n_opposite_pairs, unsigned int& n_overlap_pairs) {

    const edm::measurement_collection<default_algebra>::const_device
        measurements(measurements_view);
    const strip_measurement_surface_info_collection_types::const_device
        surface_infos(surface_infos_view);
    (void)det_view;
    if (globalIndex >= measurements.size()) {
        return;
    }

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
    // Keep the counting and finding passes identical. The finding pass skips
    // measurements whose strip endpoints cannot be reconstructed; counting
    // them here would leave unwritten entries in the allocated pair buffers.
    const details::strip_material reference_material =
        details::make_strip_material(reference_measurement, reference_info,
                                     point3{});
    if (reference_material.valid == 0u) {
        return;
    }

    unsigned int opposite_count = 0u;
    unsigned int overlap_count = 0u;
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
            const auto relation = details::match_offline_strip_pair(
                reference_measurement, measurements.at(candidate_index),
                reference_info, candidate_info);
            if (relation == strip_pair_relation::none) {
                continue;
            }
            const details::strip_material candidate_material =
                details::make_strip_material(measurements.at(candidate_index),
                                             candidate_info, point3{});
            if (candidate_material.valid == 0u) {
                continue;
            }
            if (relation == strip_pair_relation::opposite) {
                ++opposite_count;
            } else if (details::is_overlap_relation(relation)) {
                ++overlap_count;
            }
        }
    }

    if (opposite_count > 0u) {
        vecmem::device_atomic_ref<unsigned int>(n_opposite_pairs)
            .fetch_add(opposite_count);
    }
    if (overlap_count > 0u) {
        vecmem::device_atomic_ref<unsigned int>(n_overlap_pairs)
            .fetch_add(overlap_count);
    }
}

}  // namespace traccc::device
