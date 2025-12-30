//! JSON-based FFI API for physics sim Schrodinger quantum functions.

use std::os::raw::c_char;

use num_complex::Complex;
use serde::Deserialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_sim::schrodinger_quantum::SchrodingerParameters;
use crate::physics::physics_sim::schrodinger_quantum::{
    self,
};

#[derive(Deserialize)]

struct SchrodingerInput {
    params: SchrodingerParameters,
    initial_psi_re: Vec<f64>,
    initial_psi_im: Vec<f64>,
}

/// Solves the time-dependent Schrödinger equation for quantum wavefunction evolution via JSON serialization.
///
/// The Schrödinger equation iℏ∂ψ/∂t = Ĥψ where Ĥ = -ℏ²∇²/(2m) + V(r) governs quantum
/// mechanical evolution of the wavefunction ψ(r,t) under a potential V. This solver uses
/// the Crank-Nicolson or split-operator method for unitary time evolution.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `params`: Schrödinger parameters including:
///     - `n_points`: Number of spatial grid points
///     - `dx`: Spatial discretization step
///     - `dt`: Time step size
///     - `steps`: Number of time steps
///     - `potential`: External potential V(x) values
///   - `initial_psi_re`: Real part of initial wavefunction ψ(x,0)
///   - `initial_psi_im`: Imaginary part of initial wavefunction ψ(x,0)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<f64>, String>` with
/// the final probability density |ψ(x,t)|² as a vector.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw C string pointer.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_physics_sim_schrodinger_run_json(
    input: *const c_char
) -> *mut c_char {

    let input : SchrodingerInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<f64>,
                    String,
                >::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        },
    };

    let mut initial_psi: Vec<
        Complex<f64>,
    > = input
        .initial_psi_re
        .iter()
        .zip(
            input
                .initial_psi_im
                .iter(),
        )
        .map(|(&r, &i)| {

            Complex::new(r, i)
        })
        .collect();

    match schrodinger_quantum::run_schrodinger_simulation(
        &input.params,
        &mut initial_psi,
    ) {
        | Ok(snapshots) => {
            if let Some(final_state) = snapshots.last() {

                to_c_string(
                    serde_json::to_string(&FfiResult::<
                        Vec<f64>,
                        String,
                    >::ok(
                        final_state
                            .clone()
                            .into_raw_vec_and_offset().0,
                    ))
                    .unwrap(),
                )
            } else {

                to_c_string(
                    serde_json::to_string(&FfiResult::<
                        Vec<f64>,
                        String,
                    >::err(
                        "No snapshots produced".to_string(),
                    ))
                    .unwrap(),
                )
            }
        },
        | Err(e) => {
            to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<f64>,
                    String,
                >::err(
                    e
                ))
                .unwrap(),
            )
        },
    }
}
