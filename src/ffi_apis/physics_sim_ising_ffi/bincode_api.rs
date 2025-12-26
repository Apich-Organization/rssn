//! Bincode-based FFI API for physics sim Ising statistical functions.

use serde::Deserialize;
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

#[no_mangle]

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
