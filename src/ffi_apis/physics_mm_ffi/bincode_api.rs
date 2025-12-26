//! Bincode-based FFI API for physics MM functions.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_mm::SPHSystem;

#[derive(Deserialize)]

struct SphInput {
    system : SPHSystem,
    dt : f64,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_mm_sph_update_bincode(
    buffer : BincodeBuffer
) -> BincodeBuffer {

    let mut input: SphInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    SPHSystem,
                    String,
                >::err(
                    "Invalid Bincode".to_string(),
                ))
            },
        };

    input
        .system
        .update(input.dt);

    to_bincode_buffer(&FfiResult::<
        SPHSystem,
        String,
    >::ok(
        input.system
    ))
}
