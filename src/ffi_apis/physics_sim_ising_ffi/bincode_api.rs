//! Bincode-based FFI API for physics sim Ising statistical functions.

use serde::Serialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_sim::ising_statistical::IsingParameters;
use crate::physics::physics_sim::ising_statistical::{
    self,
};

#[derive(Serialize)]

struct IsingOutput {
    pub grid: Vec<i8>,
    pub magnetization: f64,
}

/// Runs a 2D Ising model Monte Carlo simulation using the Metropolis algorithm via bincode serialization.
///
/// The Ising model with Hamiltonian H = -J∑⟨i,j⟩sᵢsⱼ - h∑ᵢsᵢ describes phase transitions
/// in magnetic systems. The simulation uses Metropolis-Hastings sampling to evolve spin
/// configurations toward thermal equilibrium.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `IsingParameters` with:
///   - `width`, `height`: Grid dimensions
///   - `temperature`: Temperature T in units of `J/k_B`
///   - `mc_steps`: Number of Monte Carlo sweeps to perform
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<IsingOutput, String>` with either:
/// - `ok`: Object containing:
///   - `grid`: Final spin configuration (±1 values)
///   - `magnetization`: Average magnetization M = ⟨∑ᵢsᵢ⟩/N
/// - `err`: Error message if computation failed
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_physics_sim_ising_run_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let params : IsingParameters = match from_bincode_buffer(&buffer) {
        | Some(p) => p,
        | None => {
            return to_bincode_buffer(&FfiResult::<
                IsingOutput,
                String,
            >::err(
                "Invalid Bincode".to_string(),
            ))
        },
    };

    let (grid, magnetization) = ising_statistical::run_ising_simulation(&params);

    let out = IsingOutput {
        grid,
        magnetization,
    };

    to_bincode_buffer(&FfiResult::<
        IsingOutput,
        String,
    >::ok(out))
}
