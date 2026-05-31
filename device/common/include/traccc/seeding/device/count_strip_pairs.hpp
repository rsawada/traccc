/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/global_index.hpp"

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/seeding/detail/strip_pair.hpp"

namespace traccc::device {

/// Count barrel strip measurement pairs that can form strip spacepoints.
template <typename detector_t>
TRACCC_HOST_DEVICE inline void count_strip_pairs(
    global_index_t globalIndex, typename detector_t::view det_view,
    const edm::measurement_collection<default_algebra>::const_view&
        measurements_view,
    const strip_pair_config& config, unsigned int& n_pairs);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/seeding/device/impl/count_strip_pairs.ipp"
