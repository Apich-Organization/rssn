//! Bincode-based FFI API for numerical topology.

use serde::Deserialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::topology::PersistenceDiagram;
use crate::numerical::topology::{
    self,
};

#[derive(Deserialize)]

struct BettiInput {
    points: Vec<Vec<f64>>,
    epsilon: f64,
    max_dim: usize,
}

/// Computes Betti numbers at a fixed radius for topological data analysis using bincode serialization.
///
/// Betti numbers characterize topological features: β₀ (connected components),
/// β₁ (holes/loops), β₂ (voids), etc., in the Vietoris-Rips complex.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `BettiInput` with:
///   - `points`: Point cloud data as vectors of coordinates
///   - `epsilon`: Radius parameter for Vietoris-Rips complex
///   - `max_dim`: Maximum homology dimension to compute
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<Vec<usize>, String>` with either:
/// - `ok`: Vector of Betti numbers [β₀, β₁, β₂, ...] up to `max_dim`
/// - `err`: Error message if computation failed
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_topology_betti_numbers_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : BettiInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<Vec<usize>, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
            )
        },
    };

    let pt_slices: Vec<&[f64]> = input
        .points
        .iter()
        .map(std::vec::Vec::as_slice)
        .collect();

    let res = topology::betti_numbers_at_radius(
        &pt_slices,
        input.epsilon,
        input.max_dim,
    );

    let ffi_res = FfiResult {
        ok: Some(res),
        err: None::<String>,
    };

    to_bincode_buffer(&ffi_res)
}

#[derive(Deserialize)]

struct PersistenceInput {
    points: Vec<Vec<f64>>,
    max_epsilon: f64,
    steps: usize,
    max_dim: usize,
}

/// Computes persistent homology for topological data analysis using bincode serialization.
///
/// Tracks the birth and death of topological features (components, holes, voids)
/// across multiple scales, producing persistence diagrams.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `PersistenceInput` with:
///   - `points`: Point cloud data as vectors of coordinates
///   - `max_epsilon`: Maximum radius to analyze
///   - `steps`: Number of radius values to sample
///   - `max_dim`: Maximum homology dimension to compute
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<Vec<PersistenceDiagram>, String>` with either:
/// - `ok`: Persistence diagrams for each dimension
/// - `err`: Error message if computation failed
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_topology_persistence_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : PersistenceInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<Vec<PersistenceDiagram>, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
            )
        },
    };

    let res =
        topology::compute_persistence(
            &input.points,
            input.max_epsilon,
            input.steps,
            input.max_dim,
        );

    let ffi_res = FfiResult {
        ok: Some(res),
        err: None::<String>,
    };

    to_bincode_buffer(&ffi_res)
}
