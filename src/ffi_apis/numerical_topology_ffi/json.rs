//! JSON-based FFI API for numerical topology.

use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
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

/// Computes Betti numbers at a fixed radius for topological data analysis via JSON serialization.
///
/// Betti numbers characterize topological features: β₀ (connected components),
/// β₁ (holes/loops), β₂ (voids), etc., in the Vietoris-Rips complex.
///
/// # Arguments
///
/// * `input_json` - A JSON string pointer containing:
///   - `points`: Point cloud data as arrays of coordinates
///   - `epsilon`: Radius parameter for Vietoris-Rips complex
///   - `max_dim`: Maximum homology dimension to compute
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<usize>, String>` with
/// a vector of Betti numbers [β₀, β₁, β₂, ...] up to `max_dim`.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_topology_betti_numbers_json(
    input_json: *const c_char
) -> *mut c_char {

    let input : BettiInput = match from_json_string(input_json) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<usize>, String> {
                        ok : None,
                        err : Some("Invalid JSON input".to_string()),
                    },
                )
                .unwrap(),
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

    to_c_string(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
}

#[derive(Deserialize)]

struct PersistenceInput {
    points: Vec<Vec<f64>>,
    max_epsilon: f64,
    steps: usize,
    max_dim: usize,
}

/// Computes persistent homology for topological data analysis via JSON serialization.
///
/// Tracks the birth and death of topological features (components, holes, voids)
/// across multiple scales, producing persistence diagrams.
///
/// # Arguments
///
/// * `input_json` - A JSON string pointer containing:
///   - `points`: Point cloud data as arrays of coordinates
///   - `max_epsilon`: Maximum radius to analyze
///   - `steps`: Number of radius values to sample
///   - `max_dim`: Maximum homology dimension to compute
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<PersistenceDiagram>, String>` with
/// persistence diagrams for each dimension.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_topology_persistence_json(
    input_json: *const c_char
) -> *mut c_char {

    let input : PersistenceInput = match from_json_string(input_json) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<PersistenceDiagram>, String> {
                        ok : None,
                        err : Some("Invalid JSON input".to_string()),
                    },
                )
                .unwrap(),
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

    to_c_string(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
}
