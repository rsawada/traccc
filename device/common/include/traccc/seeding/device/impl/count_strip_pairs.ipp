/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// VecMem include(s).
#include <vecmem/memory/device_atomic_ref.hpp>

namespace traccc::device {

template <typename detector_t>
TRACCC_HOST_DEVICE inline void count_strip_pairs(
    const global_index_t globalIndex, typename detector_t::view det_view,
    const edm::measurement_collection<default_algebra>::const_view&
        measurements_view,
    const strip_pair_config& config, unsigned int& n_pairs) {

    const edm::measurement_collection<default_algebra>::const_device
        measurements(measurements_view);
    if (globalIndex >= measurements.size()) {
        return;
    }

    typename detector_t::device det(det_view);
    const edm::measurement inner_measurement = measurements.at(globalIndex);
    unsigned int n_compatible_pairs = 0u;

    for (unsigned int other_index = 0u; other_index < measurements.size();
         ++other_index) {
        const edm::measurement outer_measurement = measurements.at(other_index);
        if (details::is_compatible_barrel_strip_pair(
                det, inner_measurement, outer_measurement, config)) {
            ++n_compatible_pairs;
        }
    }

    if (n_compatible_pairs > 0u) {
        vecmem::device_atomic_ref<unsigned int> total_pairs(n_pairs);
        total_pairs.fetch_add(n_compatible_pairs);
    }
}

}  // namespace traccc::device
