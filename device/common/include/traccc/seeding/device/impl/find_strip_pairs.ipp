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
TRACCC_HOST_DEVICE inline void find_strip_pairs(
    const global_index_t globalIndex, typename detector_t::view det_view,
    const edm::measurement_collection<default_algebra>::const_view& measurements_view,
    const strip_measurement_surface_info_collection_types::const_view& surface_infos_view,
    const vecmem::data::vector_view<unsigned int>& candidate_indices_view,
    unsigned int& opposite_position, unsigned int& overlap_position,
    strip_pair_collection_types::view opposite_pairs_view,
    strip_pair_collection_types::view overlap_pairs_view) {

    const edm::measurement_collection<default_algebra>::const_device measurements(measurements_view);
    const strip_measurement_surface_info_collection_types::const_device surface_infos(surface_infos_view);
    const vecmem::device_vector<unsigned int> candidate_indices(
        candidate_indices_view);
    (void)det_view;
    if ((globalIndex >= measurements.size()) || (globalIndex >= surface_infos.size())) {
        return;
    }

    vecmem::device_vector<strip_pair> opposite_pairs(opposite_pairs_view);
    vecmem::device_vector<strip_pair> overlap_pairs(overlap_pairs_view);
    const edm::measurement reference_measurement = measurements.at(globalIndex);
    const strip_measurement_surface_info reference_info = surface_infos.at(globalIndex);

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
        const edm::measurement candidate_measurement = measurements.at(candidate_index);
        const strip_measurement_surface_info candidate_info = surface_infos.at(candidate_index);
        const auto relation = details::match_offline_strip_pair(
            reference_measurement, candidate_measurement, reference_info, candidate_info);
        if (relation == strip_pair_relation::none) {
            continue;
        }

        const bool barrel = reference_info.has_barrel_material != 0u;
        const point3 reference_center = barrel ? reference_info.barrel_strip_center
                                               : reference_info.endcap_strip_center;
        const point3 candidate_center = barrel ? candidate_info.barrel_strip_center
                                               : candidate_info.endcap_strip_center;
        const vector3 reference_normal = barrel ? reference_info.barrel_strip_normal
                                                : reference_info.endcap_strip_normal;
        const vector3 candidate_normal = barrel ? candidate_info.barrel_strip_normal
                                                : candidate_info.endcap_strip_normal;
        const scalar reference_r = std::sqrt(reference_center[0] * reference_center[0] +
                                             reference_center[1] * reference_center[1]);
        const scalar candidate_r = std::sqrt(candidate_center[0] * candidate_center[0] +
                                             candidate_center[1] * candidate_center[1]);
        const scalar dx = reference_center[0] - candidate_center[0];
        const scalar dy = reference_center[1] - candidate_center[1];
        const strip_pair pair{
            static_cast<unsigned int>(globalIndex), candidate_index,
            reference_measurement.surface_link().value(), candidate_measurement.surface_link().value(),
            candidate_r - reference_r, std::sqrt(dx * dx + dy * dy),
            reference_center[2] - candidate_center[2],
            reference_normal[0] * candidate_normal[0] +
                reference_normal[1] * candidate_normal[1] +
                reference_normal[2] * candidate_normal[2],
            reference_info.is_endcap, reference_info.mid_r, candidate_info.mid_r,
            reference_info.strip_half_length, candidate_info.strip_half_length,
            details::g80_strip_length_gap_tolerance(reference_info, candidate_info)};

        if (relation == strip_pair_relation::opposite) {
            const unsigned int position =
                vecmem::device_atomic_ref<unsigned int>(opposite_position).fetch_add(1u);
            if (position < opposite_pairs.size()) { opposite_pairs.at(position) = pair; }
        } else {
            const unsigned int position =
                vecmem::device_atomic_ref<unsigned int>(overlap_position).fetch_add(1u);
            if (position < overlap_pairs.size()) { overlap_pairs.at(position) = pair; }
        }
    }
}

}  // namespace traccc::device
