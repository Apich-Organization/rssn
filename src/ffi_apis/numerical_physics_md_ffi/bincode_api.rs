//! Bincode-based FFI API for numerical MD functions.

use crate::ffi_apis::common::{from_bincode_buffer, to_bincode_buffer, BincodeBuffer};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::physics_md;
use serde::{Deserialize, Serialize};

#[derive(Deserialize)]
struct LennardJonesInput {
    p1_position: Vec<f64>,
    p2_position: Vec<f64>,
    epsilon: f64,
    sigma: f64,
}

#[derive(Serialize)]
struct InteractionOutput {
    potential: f64,
    force: Vec<f64>,
}

#[derive(Deserialize)]
struct PbcInput {
    position: Vec<f64>,
    box_size: Vec<f64>,
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_md_lennard_jones_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {
    let input: LennardJonesInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<InteractionOutput, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };
    
    let p1 = physics_md::Particle::new(0, 1.0, input.p1_position, vec![0.0, 0.0, 0.0]);
    let p2 = physics_md::Particle::new(1, 1.0, input.p2_position, vec![0.0, 0.0, 0.0]);
    
    match physics_md::lennard_jones_interaction(&p1, &p2, input.epsilon, input.sigma) {
        Ok((potential, force)) => {
            let output = InteractionOutput { potential, force };
            to_bincode_buffer(&FfiResult {
                ok: Some(output),
                err: None::<String>,
            })
        }
        Err(e) => to_bincode_buffer(&FfiResult::<InteractionOutput, String> {
            ok: None,
            err: Some(e),
        }),
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_md_apply_pbc_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {
    let input: PbcInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<Vec<f64>, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };
    
    let wrapped = physics_md::apply_pbc(&input.position, &input.box_size);
    
    to_bincode_buffer(&FfiResult {
        ok: Some(wrapped),
        err: None::<String>,
    })
}
