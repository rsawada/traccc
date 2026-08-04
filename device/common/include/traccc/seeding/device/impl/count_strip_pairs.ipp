/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/device_atomic_ref.hpp>

namespace traccc::device {

template <typename detector_t>
TRACCC_HOST_DEVICE inline void count_strip_pairs(
    const global_index_t globalIndex, typename detector_t::view det_view,
    const edm::measurement_collection<default_algebra>::const_view& measurements_view,
    const strip_measurement_surface_info_collection_types::const_view& surface_infos_view,
    const vecmem::data::vector_view<unsigned int>& candidate_indices_view,
    unsigned int& n_opposite_pairs, unsigned int& n_overlap_pairs) {

    const edm::measurement_collection<default_algebra>::const_device measurements(measurements_view);
    const strip_measurement_surface_info_collection_types::const_device surface_infos(surface_infos_view);
    const vecmem::device_vector<unsigned int> candidate_indices(
        candidate_indices_view);
    (void)det_view;
    if ((globalIndex >= measurements.size()) || (globalIndex >= surface_infos.size())) {
        return;
    }

    const edm::measurement reference_measurement = measurements.at(globalIndex);
    const strip_measurement_surface_info reference_info = surface_infos.at(globalIndex);
    unsigned int opposite_count = 0u;
    unsigned int overlap_count = 0u;

    const unsigned int candidate_end =
        reference_info.candidate_measurement_begin +
        reference_info.candidate_measurement_count;
    for (unsigned int position = reference_info.candidate_measurement_begin;
         position < candidate_end; ++position) {
        if (position >= candidate_indices.size()) {
            break;
        }
        const unsigned int candidate_index = candidate_indices.at(position);
        if ((candidate_index == globalIndex) ||
            (candidate_index >= measurements.size()) ||
            (candidate_index >= surface_infos.size())) {
            continue;
        }
        const auto relation = details::match_offline_strip_pair(
            reference_measurement, measurements.at(candidate_index), reference_info,
            surface_infos.at(candidate_index));
        if (relation == strip_pair_relation::opposite) {
            ++opposite_count;
        } else if (details::is_overlap_relation(relation)) {
            ++overlap_count;
        }
    }

    if (opposite_count > 0u) {
        vecmem::device_atomic_ref<unsigned int>(n_opposite_pairs).fetch_add(opposite_count);
    }
    if (overlap_count > 0u) {
        vecmem::device_atomic_ref<unsigned int>(n_overlap_pairs).fetch_add(overlap_count);
    }
}

}  // namespace traccc::device
