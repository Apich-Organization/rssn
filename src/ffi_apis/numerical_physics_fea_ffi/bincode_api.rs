//! Bincode-based FFI API for numerical FEA functions.

use crate::ffi_apis::common::{
    from_bincode_buffer,
    to_bincode_buffer,
    BincodeBuffer,
};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::physics_fea;
use serde::{
    Deserialize,
    Serialize,
};

#[derive(Deserialize)]

struct LinearElement1DInput {
    length: f64,
    youngs_modulus: f64,
    area: f64,
}

#[derive(Deserialize)]

struct StressInput {
    sx: f64,
    sy: f64,
    txy: f64,
}

#[derive(Serialize)]

struct PrincipalStressOutput {
    sigma1: f64,
    sigma2: f64,
    angle: f64,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_fea_linear_element_1d_stiffness_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: LinearElement1DInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok: None,
                    err: Some("Invalid Bincode".to_string()),
                },
            )
        }
    };

    let stiffness = input.youngs_modulus * input.area / input.length;

    to_bincode_buffer(&FfiResult {
        ok: Some(stiffness),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_fea_von_mises_stress_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: StressInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok: None,
                    err: Some("Invalid Bincode".to_string()),
                },
            )
        }
    };

    let vm = physics_fea::TriangleElement2D::von_mises_stress(&[
        input.sx, input.sy, input.txy,
    ]);

    to_bincode_buffer(&FfiResult {
        ok: Some(vm),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_fea_principal_stresses_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: StressInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(
                &FfiResult::<PrincipalStressOutput, String> {
                    ok: None,
                    err: Some("Invalid Bincode".to_string()),
                },
            )
        }
    };

    let (sigma1, sigma2, angle) = physics_fea::principal_stresses(&[
        input.sx, input.sy, input.txy,
    ]);

    let output = PrincipalStressOutput {
        sigma1,
        sigma2,
        angle,
    };

    to_bincode_buffer(&FfiResult {
        ok: Some(output),
        err: None::<String>,
    })
}
