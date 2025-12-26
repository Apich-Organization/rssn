//! Bincode-based FFI API for physics sim FDTD electrodynamics functions.

use crate::ffi_apis::common::{
    from_bincode_buffer,
    to_bincode_buffer,
    BincodeBuffer,
};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_sim::fdtd_electrodynamics::{
    self,
    FdtdParameters,
};
use ndarray::Array2;

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_sim_fdtd_run_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let params: FdtdParameters = match from_bincode_buffer(&buffer) {
        | Some(p) => p,
        | None => {
            return to_bincode_buffer(&FfiResult::<
                Vec<f64>,
                String,
            >::err(
                "Invalid Bincode".to_string(),
            ))
        },
    };

    let snapshots = fdtd_electrodynamics::run_fdtd_simulation(&params);

    if let Some(final_ez) = snapshots.last() {

        to_bincode_buffer(&FfiResult::<
            Array2<f64>,
            String,
        >::ok(
            final_ez.clone()
        ))
    } else {

        to_bincode_buffer(&FfiResult::<
            Array2<f64>,
            String,
        >::err(
            "No snapshots".to_string(),
        ))
    }
}
