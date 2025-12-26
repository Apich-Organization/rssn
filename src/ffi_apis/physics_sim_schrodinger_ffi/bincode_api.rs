//! Bincode-based FFI API for physics sim Schrodinger quantum functions.

use crate::ffi_apis::common::{from_bincode_buffer, to_bincode_buffer, BincodeBuffer};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_sim::schrodinger_quantum::{self, SchrodingerParameters};
use num_complex::Complex;
use serde::{Deserialize, Serialize};

#[derive(Deserialize)]
struct SchrodingerInput {
    params: SchrodingerParameters,
    initial_psi_re: Vec<f64>,
    initial_psi_im: Vec<f64>,
}

#[no_mangle]
pub unsafe extern "C" fn rssn_physics_sim_schrodinger_run_bincode(buffer: BincodeBuffer) -> BincodeBuffer {
    let input: SchrodingerInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => return to_bincode_buffer(&FfiResult::<Vec<f64>, String>::err("Invalid Bincode".to_string())),
    };

    let mut initial_psi: Vec<Complex<f64>> = input.initial_psi_re.iter().zip(input.initial_psi_im.iter())
        .map(|(&r, &i)| Complex::new(r, i))
        .collect();

    match schrodinger_quantum::run_schrodinger_simulation(&input.params, &mut initial_psi) {
        Ok(snapshots) => {
            if let Some(final_state) = snapshots.last() {
                to_bincode_buffer(&FfiResult::<Vec<f64>, String>::ok(final_state.clone().into_raw_vec()))
            } else {
                to_bincode_buffer(&FfiResult::<Vec<f64>, String>::err("No snapshots".to_string()))
            }
        }
        Err(e) => to_bincode_buffer(&FfiResult::<Vec<f64>, String>::err(e)),
    }
}
