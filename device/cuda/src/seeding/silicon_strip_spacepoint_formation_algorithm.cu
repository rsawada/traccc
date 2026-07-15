/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../utils/cuda_error_handling.hpp"
#include "../utils/global_index.hpp"
#include "../utils/utils.hpp"
#include "traccc/cuda/seeding/silicon_strip_spacepoint_formation_algorithm.hpp"

// Project include(s).
#include "traccc/geometry/detector.hpp"
#include "traccc/seeding/device/count_strip_pairs.hpp"
#include "traccc/seeding/device/find_strip_pairs.hpp"
#include "traccc/seeding/device/form_spacepoints.hpp"

namespace traccc::cuda {
namespace kernels {

/// Kernel wrapping @c device::count_strip_pairs.
template <typename detector_t>
__global__ void __launch_bounds__(1024, 1) count_strip_pairs_kernel(
    typename detector_t::view detector,
    typename edm::measurement_collection<
        typename detector_t::device::algebra_type>::const_view measurements,
    strip_measurement_surface_info_collection_types::const_view surface_infos,
    barrel_strip_pair_config barrel_config,
    endcap_strip_pair_config endcap_config, unsigned int& n_pairs,
    unsigned int& n_barrel_pairs, unsigned int& n_endcap_pairs,
    unsigned int& n_endcap_boundary_pairs)
    requires(traccc::is_detector_traits<detector_t>)
{
    device::count_strip_pairs<detector_t>(details::global_index1(), detector,
                                          measurements, surface_infos,
                                          barrel_config, endcap_config, n_pairs,
                                          n_barrel_pairs, n_endcap_pairs,
                                          n_endcap_boundary_pairs);
}

/// Kernel wrapping @c device::find_strip_pairs.
template <typename detector_t>
__global__ void __launch_bounds__(1024, 1) find_strip_pairs_kernel(
    typename detector_t::view detector,
    typename edm::measurement_collection<
        typename detector_t::device::algebra_type>::const_view measurements,
    strip_measurement_surface_info_collection_types::const_view surface_infos,
    barrel_strip_pair_config barrel_config,
    endcap_strip_pair_config endcap_config, unsigned int& pair_position,
    strip_pair_collection_types::view pairs)
    requires(traccc::is_detector_traits<detector_t>)
{
    device::find_strip_pairs<detector_t>(details::global_index1(), detector,
                                         measurements, surface_infos,
                                         barrel_config, endcap_config,
                                         pair_position, pairs);
}

/// Kernel wrapping @c device::form_barrel_strip_spacepoints.
template <typename detector_t>
__global__ void __launch_bounds__(1024, 1) form_barrel_strip_spacepoints_kernel(
    typename detector_t::view detector,
    typename edm::measurement_collection<
        typename detector_t::device::algebra_type>::const_view measurements,
    strip_pair_collection_types::const_view pairs,
    strip_measurement_surface_info_collection_types::const_view surface_infos,
    edm::spacepoint_collection::view spacepoints)
    requires(traccc::is_detector_traits<detector_t>)
{
    device::form_barrel_strip_spacepoints<detector_t>(
        details::global_index1(), detector, measurements, pairs, surface_infos,
        spacepoints);
}

}  // namespace kernels

silicon_strip_spacepoint_formation_algorithm::
    silicon_strip_spacepoint_formation_algorithm(
        const traccc::memory_resource& mr, vecmem::copy& copy,
        cuda::stream& str, std::unique_ptr<const Logger> logger)
    : device::silicon_strip_spacepoint_formation_algorithm(mr, copy,
                                                           std::move(logger)),
      cuda::algorithm_base(str) {}

void silicon_strip_spacepoint_formation_algorithm::count_strip_pairs_kernel(
    const count_strip_pairs_kernel_payload& payload) const {

    const unsigned int n_threads = warp_size() * 8;
    const unsigned int n_blocks =
        (payload.n_measurements + n_threads - 1) / n_threads;
    detector_buffer_visitor<detector_type_list>(
        payload.detector, [&]<typename detector_traits_t>(
                              const typename detector_traits_t::view& det) {
            kernels::count_strip_pairs_kernel<detector_traits_t>
                <<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(
                    det, payload.measurements, payload.surface_infos,
                    payload.barrel_config, payload.endcap_config, payload.n_pairs,
                    payload.n_barrel_pairs, payload.n_endcap_pairs,
                    payload.n_endcap_boundary_pairs);
        });
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void silicon_strip_spacepoint_formation_algorithm::find_strip_pairs_kernel(
    const find_strip_pairs_kernel_payload& payload) const {

    const unsigned int n_threads = warp_size() * 8;
    const unsigned int n_blocks =
        (payload.n_measurements + n_threads - 1) / n_threads;
    detector_buffer_visitor<detector_type_list>(
        payload.detector, [&]<typename detector_traits_t>(
                              const typename detector_traits_t::view& det) {
            kernels::find_strip_pairs_kernel<detector_traits_t>
                <<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(
                    det, payload.measurements, payload.surface_infos,
                    payload.barrel_config, payload.endcap_config,
                    payload.pair_position, payload.pairs);
        });
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void silicon_strip_spacepoint_formation_algorithm::form_spacepoints_kernel(
    const form_spacepoints_kernel_payload& payload) const {

    const unsigned int n_threads = warp_size() * 8;
    const unsigned int n_blocks =
        (payload.n_pairs + n_threads - 1) / n_threads;
    detector_buffer_visitor<detector_type_list>(
        payload.detector, [&]<typename detector_traits_t>(
                              const typename detector_traits_t::view& det) {
            kernels::form_barrel_strip_spacepoints_kernel<detector_traits_t>
                <<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(
                    det, payload.measurements, payload.pairs,
                    payload.surface_infos,
                    payload.spacepoints);
        });
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

}  // namespace traccc::cuda
