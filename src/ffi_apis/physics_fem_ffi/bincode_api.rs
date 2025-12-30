//! Bincode-based FFI API for physics FEM functions.

use serde::Deserialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_fem;

#[derive(Deserialize)]

struct Poisson1DInput {
    n_elements: usize,
    domain_length: f64,
}

/// Solves the 1D Poisson equation using Finite Element Method (FEM) via bincode serialization.
///
/// The Poisson equation -d²u/dx² = f(x) is solved using linear finite elements
/// with homogeneous Dirichlet boundary conditions (u = 0 at boundaries).
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `Poisson1DInput` with:
///   - `n_elements`: Number of finite elements in the mesh
///   - `domain_length`: Total length of the 1D domain
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<Vec<f64>, String>` with either:
/// - `ok`: Solution vector u at nodal points
/// - `err`: Error message if computation failed
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_physics_fem_solve_poisson_1d_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : Poisson1DInput = match from_bincode_buffer(&buffer) {
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

    match physics_fem::solve_poisson_1d(
        input.n_elements,
        input.domain_length,
        |_| 2.0,
    ) {
        | Ok(res) => {
            to_bincode_buffer(
                &FfiResult::<
                    Vec<f64>,
                    String,
                >::ok(
                    res
                ),
            )
        },
        | Err(e) => {
            to_bincode_buffer(
                &FfiResult::<
                    Vec<f64>,
                    String,
                >::err(
                    e
                ),
            )
        },
    }
}
