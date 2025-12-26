//! JSON-based FFI API for numerical topology.

use std::os::raw::c_char;

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::topology::PersistenceDiagram;
use crate::numerical::topology::{
    self,
};

#[derive(Deserialize)]

struct BettiInput {
    points: Vec<Vec<f64>>,
    epsilon: f64,
    max_dim: usize,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_topology_betti_numbers_json(
    input_json: *const c_char
) -> *mut c_char {

    let input : BettiInput = match from_json_string(input_json) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<usize>, String> {
                        ok : None,
                        err : Some("Invalid JSON input".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let pt_slices: Vec<&[f64]> = input
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
        ok: Some(res),
        err: None::<String>,
    };

    to_c_string(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
}

#[derive(Deserialize)]

struct PersistenceInput {
    points: Vec<Vec<f64>>,
    max_epsilon: f64,
    steps: usize,
    max_dim: usize,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_topology_persistence_json(
    input_json: *const c_char
) -> *mut c_char {

    let input : PersistenceInput = match from_json_string(input_json) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<PersistenceDiagram>, String> {
                        ok : None,
                        err : Some("Invalid JSON input".to_string()),
                    },
                )
                .unwrap(),
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
        ok: Some(res),
        err: None::<String>,
    };

    to_c_string(
        serde_json::to_string(&ffi_res)
            .unwrap(),
    )
}
