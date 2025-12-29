//! Bincode-based FFI API for physics MTM functions.

use serde::Deserialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_mtm;

#[derive(Deserialize)]

struct Multigrid2DInput {
    n: usize,
    f: Vec<f64>,
    num_cycles: usize,
}

/// Solves the 2D Poisson equation using Multigrid Method via bincode serialization.
///
/// The Poisson equation ∇²u = f is solved using the multigrid method, which achieves
/// optimal O(N) complexity through hierarchical coarse-grid correction.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `Multigrid2DInput` with:
///   - `n`: Grid size (n×n interior points)
///   - `f`: Right-hand side source term (flattened 2D array)
///   - `num_cycles`: Number of V-cycles or W-cycles to perform
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<Vec<f64>, String>` with either:
/// - `ok`: Solution vector u (flattened 2D array)
/// - `err`: Error message if computation failed
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_physics_mtm_solve_poisson_2d_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : Multigrid2DInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(&FfiResult::<
                Vec<f64>,
                String,
            >::err(
                "Invalid Bincode".to_string(),
            ))
        },
    };

    match physics_mtm::solve_poisson_2d_multigrid(
        input.n,
        &input.f,
        input.num_cycles,
    ) {
        | Ok(res) => {
            to_bincode_buffer(&FfiResult::<
                Vec<f64>,
                String,
            >::ok(
                res
            ))
        },
        | Err(e) => {
            to_bincode_buffer(&FfiResult::<
                Vec<f64>,
                String,
            >::err(
                e
            ))
        },
    }
}
