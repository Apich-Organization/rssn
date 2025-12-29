//! Bincode-based FFI API for physics sim Schrodinger quantum functions.

use num_complex::Complex;
use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
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

/// Solves the time-dependent Schrödinger equation for quantum wavefunction evolution via bincode serialization.
///
/// The Schrödinger equation iℏ∂ψ/∂t = Ĥψ where Ĥ = -ℏ²∇²/(2m) + V(r) governs quantum
/// mechanical evolution of the wavefunction ψ(r,t) under a potential V. This solver uses
/// the Crank-Nicolson or split-operator method for unitary time evolution.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `SchrodingerInput` with:
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
/// A bincode-encoded buffer containing `FfiResult<Vec<f64>, String>` with either:
/// - `ok`: Final probability density |ψ(x,t)|² as a vector
/// - `err`: Error message if computation failed or produced no snapshots
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_physics_sim_schrodinger_run_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : SchrodingerInput = match from_bincode_buffer(&buffer) {
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

                to_bincode_buffer(&FfiResult::<
                    Vec<f64>,
                    String,
                >::ok(
                    final_state
                        .clone()
                        .into_raw_vec(),
                ))
            } else {

                to_bincode_buffer(&FfiResult::<
                    Vec<f64>,
                    String,
                >::err(
                    "No snapshots".to_string(),
                ))
            }
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
