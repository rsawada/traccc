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
        measurements,
    const strip_measurement_surface_info_collection_types::const_view&
        surface_infos) const -> output_type {

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

    // Count compatible strip measurement pairs on the device.
    const barrel_strip_pair_config barrel_pair_config{};
    const endcap_strip_pair_config endcap_pair_config{};
    vecmem::data::vector_buffer<unsigned int> pair_counter_buffer(4u,
                                                                  mr().main);
    copy().setup(pair_counter_buffer)->ignore();
    copy().memset(pair_counter_buffer, 0)->ignore();
    count_strip_pairs_kernel({n_measurements, det, measurements, surface_infos,
                              barrel_pair_config, endcap_pair_config,
                              pair_counter_buffer.ptr()[0],
                              pair_counter_buffer.ptr()[1],
                              pair_counter_buffer.ptr()[2],
                              pair_counter_buffer.ptr()[3]});

    // Copy the pair count back to the host before allocating the pair buffer.
    vecmem::vector<unsigned int> pair_counter_host(
        mr().host ? mr().host : &(mr().main));
    copy()(pair_counter_buffer, pair_counter_host)->wait();
    const unsigned int n_pairs = pair_counter_host.at(0);

    // Fill the pair buffer using the same search conditions as the count pass.
    if (n_pairs == 0u) {
        return {};
    }
    strip_pair_collection_types::buffer pairs_buffer(n_pairs, mr().main);
    copy().setup(pairs_buffer)->ignore();
    copy().memset(pair_counter_buffer, 0)->ignore();
    find_strip_pairs_kernel({n_measurements, det, measurements, surface_infos,
                             barrel_pair_config, endcap_pair_config,
                             pair_counter_buffer.ptr()[0], pairs_buffer});

    edm::spacepoint_collection::buffer spacepoints(
        n_pairs, mr().main, vecmem::data::buffer_type::resizable);
    copy().setup(spacepoints)->ignore();
    form_spacepoints_kernel(
        {n_pairs, det, measurements, pairs_buffer, surface_infos, spacepoints});

    // Return the reconstructed spacepoints.
    return spacepoints;
}

}  // namespace traccc::device
