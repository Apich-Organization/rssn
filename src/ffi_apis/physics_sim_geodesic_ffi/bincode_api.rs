//! Bincode-based FFI API for physics sim geodesic relativity functions.

use crate::ffi_apis::common::{
    from_bincode_buffer,
    to_bincode_buffer,
    BincodeBuffer,
};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_sim::geodesic_relativity::{
    self,
    GeodesicParameters,
};

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_sim_geodesic_run_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let params: GeodesicParameters = match from_bincode_buffer(&buffer) {
        Some(p) => p,
        None => {
            return to_bincode_buffer(&FfiResult::<
                Vec<(f64, f64)>,
                String,
            >::err(
                "Invalid Bincode".to_string(),
            ))
        }
    };

    let path = geodesic_relativity::run_geodesic_simulation(&params);

    to_bincode_buffer(&FfiResult::<
        Vec<(f64, f64)>,
        String,
    >::ok(path))
}
