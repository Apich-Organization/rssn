//! Bincode-based FFI API for numerical topology.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::topology::PersistenceDiagram;
use crate::numerical::topology::{
    self,
};

#[derive(Deserialize)]

struct BettiInput {
    points : Vec<Vec<f64>>,
    epsilon : f64,
    max_dim : usize,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_topology_betti_numbers_bincode(
    buffer : BincodeBuffer
) -> BincodeBuffer {

    let input : BettiInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<Vec<usize>, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
            )
        },
    };

    let pt_slices : Vec<&[f64]> = input
        .points
        .iter()
        .map(|v| v.as_slice())
        .collect();

    let res = topology::betti_numbers_at_radius(
        &pt_slices,
        input.epsilon,
        input.max_dim,
    );

    let ffi_res = FfiResult {
        ok : Some(res),
        err : None::<String>,
    };

    to_bincode_buffer(&ffi_res)
}

#[derive(Deserialize)]

struct PersistenceInput {
    points : Vec<Vec<f64>>,
    max_epsilon : f64,
    steps : usize,
    max_dim : usize,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_topology_persistence_bincode(
    buffer : BincodeBuffer
) -> BincodeBuffer {

    let input : PersistenceInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<Vec<PersistenceDiagram>, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
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
        ok : Some(res),
        err : None::<String>,
    };

    to_bincode_buffer(&ffi_res)
}
