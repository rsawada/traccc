/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/algorithm_base.hpp"

// Project include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/geometry/detector_buffer.hpp"
#include "traccc/seeding/detail/strip_pair.hpp"
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/memory_resource.hpp"
#include "traccc/utils/messaging.hpp"

namespace traccc::device {

/// Algorithm forming space points out of measurements
///
/// This algorithm performs the local-to-global transformation of the 2D (strip)
/// measurements made on every detector module, into 3D spacepoint coordinates.
///
class silicon_strip_spacepoint_formation_algorithm
    : public algorithm<edm::spacepoint_collection::buffer(
          const detector_buffer&,
          const edm::measurement_collection<default_algebra>::const_view&,
          const strip_measurement_surface_info_collection_types::const_view&)>,
      public messaging,
      public algorithm_base {

    public:
    /// Constructor for spacepoint_formation algorithm
    ///
    /// @param mr The memory resource(s) to use in the algorithm
    /// @param copy The copy object to use for copying data between device
    ///             and host memory blocks
    /// @param logger The logger instance to use
    ///
    silicon_strip_spacepoint_formation_algorithm(
        const traccc::memory_resource& mr, vecmem::copy& copy,
        std::unique_ptr<const Logger> logger = getDummyLogger().clone());

    /// Construct spacepoints from 2D silicon strip measurements
    ///
    /// @param det Detector object
    /// @param measurements A collection of measurements
    /// @param surface_infos Per-measurement strip surface information
    /// @return A spacepoint buffer, with one spacepoint for every
    ///         silicon strip measurement
    ///
    output_type operator()(
        const detector_buffer& det,
        const edm::measurement_collection<default_algebra>::const_view&
            measurements,
        const strip_measurement_surface_info_collection_types::const_view&
            surface_infos) const override;

    protected:
    /// @name Function(s) to be implemented by derived classes
    /// @{

    /// Payload for the @c count_strip_pairs_kernel function.
    struct count_strip_pairs_kernel_payload {
        /// The number of measurements in the event.
        edm::measurement_collection<default_algebra>::const_view::size_type
            n_measurements;
        /// The detector object.
        const detector_buffer& detector;
        /// The input measurements.
        const edm::measurement_collection<default_algebra>::const_view&
            measurements;
        /// Configuration for the initial barrel strip pair search.
        const barrel_strip_pair_config& barrel_config;
        /// Configuration for the initial endcap strip pair search.
        const endcap_strip_pair_config& endcap_config;
        /// Total number of compatible strip pairs.
        unsigned int& n_pairs;
        /// Total number of compatible barrel strip pairs.
        unsigned int& n_barrel_pairs;
        /// Total number of compatible endcap strip pairs.
        unsigned int& n_endcap_pairs;
        /// Total number of compatible endcap pairs with boundary values.
        unsigned int& n_endcap_boundary_pairs;
    };

    /// Launch the strip pair counting kernel.
    virtual void count_strip_pairs_kernel(
        const count_strip_pairs_kernel_payload& payload) const = 0;

    /// Payload for the @c find_strip_pairs_kernel function.
    struct find_strip_pairs_kernel_payload {
        /// The number of measurements in the event.
        edm::measurement_collection<default_algebra>::const_view::size_type
            n_measurements;
        /// The detector object.
        const detector_buffer& detector;
        /// The input measurements.
        const edm::measurement_collection<default_algebra>::const_view&
            measurements;
        /// Per-measurement strip surface information.
        const strip_measurement_surface_info_collection_types::const_view&
            surface_infos;
        /// Configuration for the initial barrel strip pair search.
        const barrel_strip_pair_config& barrel_config;
        /// Configuration for the initial endcap strip pair search.
        const endcap_strip_pair_config& endcap_config;
        /// Next position in the output pair buffer.
        unsigned int& pair_position;
        /// Output strip pairs.
        strip_pair_collection_types::view& pairs;
    };

    /// Launch the strip pair finding kernel.
    virtual void find_strip_pairs_kernel(
        const find_strip_pairs_kernel_payload& payload) const = 0;

    /// Payload for the @c form_spacepoints_kernel function
    struct form_spacepoints_kernel_payload {
        /// The number of compatible barrel strip pairs in the event.
        strip_pair_collection_types::const_view::size_type n_pairs;
        /// The detector object.
        const detector_buffer& detector;
        /// The input measurements.
        const edm::measurement_collection<default_algebra>::const_view&
            measurements;
        /// The compatible barrel strip pairs.
        const strip_pair_collection_types::const_view& pairs;
        /// The output spacepoints.
        edm::spacepoint_collection::view& spacepoints;
    };

    /// Launch the spacepoint formation kernel
    ///
    /// @param payload The payload for the kernel
    ///
    virtual void form_spacepoints_kernel(
        const form_spacepoints_kernel_payload& payload) const = 0;

    /// @}

};  // class silicon_strip_spacepoint_formation_algorithm

}  // namespace traccc::device
