/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/seeding/device/silicon_strip_spacepoint_formation_algorithm.hpp"

// VecMem include(s).
#include <vecmem/containers/vector.hpp>
#include <vecmem/containers/data/vector_buffer.hpp>

namespace traccc::device {

silicon_strip_spacepoint_formation_algorithm::
    silicon_strip_spacepoint_formation_algorithm(
        const traccc::memory_resource& mr, vecmem::copy& copy,
        std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)), algorithm_base(mr, copy) {}

auto silicon_strip_spacepoint_formation_algorithm::operator()(
    const detector_buffer& det,
    const edm::measurement_collection<default_algebra>::const_view&
        measurements) const -> output_type {

    // Get the number of measurements. In an asynchronous way if possible.
    edm::measurement_collection<default_algebra>::const_view::size_type
        n_measurements = 0u;
    if (mr().host) {
        vecmem::async_size size = copy().get_size(measurements, *(mr().host));
        // Here we could give control back to the caller, once our code allows
        // for it. (coroutines...)
        n_measurements = size.get();
    } else {
        n_measurements = copy().get_size(measurements);
    }

    // If there are no measurements, return right away.
    if (n_measurements == 0) {
        return {};
    }

    // Count compatible barrel strip measurement pairs on the device.
    const strip_pair_config pair_config{};
    vecmem::data::vector_buffer<unsigned int> pair_counter_buffer(1u,
                                                                  mr().main);
    copy().setup(pair_counter_buffer)->ignore();
    copy().memset(pair_counter_buffer, 0)->ignore();
    count_strip_pairs_kernel({n_measurements, det, measurements, pair_config,
                              pair_counter_buffer.ptr()[0]});

    // Copy the pair count back to the host before allocating the pair buffer.
    vecmem::vector<unsigned int> pair_counter_host(
        mr().host ? mr().host : &(mr().main));
    copy()(pair_counter_buffer, pair_counter_host)->wait();
    const unsigned int n_pairs = pair_counter_host.at(0);
    TRACCC_INFO("compatible barrel strip pairs: " << n_pairs);

    // Fill the pair buffer using the same search conditions as the count pass.
    if (n_pairs == 0u) {
        return {};
    }
    strip_pair_collection_types::buffer pairs_buffer(n_pairs, mr().main);
    copy().setup(pairs_buffer)->ignore();
    copy().memset(pair_counter_buffer, 0)->ignore();
    find_strip_pairs_kernel({n_measurements, det, measurements, pair_config,
                             pair_counter_buffer.ptr()[0], pairs_buffer});

    // Copy a small sample to the host so the initial thresholds can be
    // checked while pair-based strip spacepoint filling is validated.
    strip_pair_collection_types::host pairs_host(
        mr().host ? mr().host : &(mr().main));
    copy()(pairs_buffer, pairs_host)->wait();
    const std::size_t n_pairs_to_print =
        (pairs_host.size() < 100u ? pairs_host.size() : 100u);
    TRACCC_INFO(
        "strip_pair_candidate_csv,index,measurement_index_1,"
        "measurement_index_2,surface_link_1,surface_link_2,"
        "surface_delta_r,strip_center_delta_xy,strip_center_delta_z,"
        "normal_dot");
    for (std::size_t i = 0u; i < n_pairs_to_print; ++i) {
        const strip_pair& pair = pairs_host.at(i);
        TRACCC_INFO("strip_pair_candidate_csv,"
                   << i << "," << pair.measurement_index_1 << ","
                   << pair.measurement_index_2 << "," << pair.surface_link_1
                   << "," << pair.surface_link_2 << ","
                   << pair.surface_delta_r << ","
                   << pair.strip_center_delta_xy << ","
                   << pair.strip_center_delta_z << "," << pair.normal_dot);
    }

    edm::spacepoint_collection::buffer spacepoints(
        n_pairs, mr().main, vecmem::data::buffer_type::resizable);
    copy().setup(spacepoints)->ignore();
    form_spacepoints_kernel(
        {n_pairs, det, measurements, pairs_buffer, spacepoints});

    // Return the reconstructed spacepoints.
    return spacepoints;
}

}  // namespace traccc::device
