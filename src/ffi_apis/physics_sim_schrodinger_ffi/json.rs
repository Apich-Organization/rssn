//! JSON-based FFI API for physics sim Schrodinger quantum functions.

use crate::ffi_apis::common::{from_json_string, to_c_string};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_sim::schrodinger_quantum::{self, SchrodingerParameters};
use num_complex::Complex;
use serde::{Deserialize, Serialize};
use std::os::raw::c_char;

#[derive(Deserialize)]

struct SchrodingerInput {
    params: SchrodingerParameters,
    initial_psi_re: Vec<f64>,
    initial_psi_im: Vec<f64>,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_sim_schrodinger_run_json(
    input: *const c_char
) -> *mut c_char {

    let input: SchrodingerInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<f64>,
                    String,
                >::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        }
    };

    let mut initial_psi: Vec<Complex<f64>> = input
        .initial_psi_re
        .iter()
        .zip(
            input
                .initial_psi_im
                .iter(),
        )
        .map(|(&r, &i)| Complex::new(r, i))
        .collect();

    match schrodinger_quantum::run_schrodinger_simulation(
        &input.params,
        &mut initial_psi,
    ) {
        Ok(snapshots) => {
            if let Some(final_state) = snapshots.last() {

                to_c_string(
                    serde_json::to_string(&FfiResult::<
                        Vec<f64>,
                        String,
                    >::ok(
                        final_state
                            .clone()
                            .into_raw_vec(),
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
        }
        Err(e) => to_c_string(
            serde_json::to_string(&FfiResult::<
                Vec<f64>,
                String,
            >::err(
                e
            ))
            .unwrap(),
        ),
    }
}
