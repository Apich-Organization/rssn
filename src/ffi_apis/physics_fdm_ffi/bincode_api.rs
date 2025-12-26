//! Bincode-based FFI API for physics FDM functions.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_fdm::Dimensions;
use crate::physics::physics_fdm::FdmGrid;
use crate::physics::physics_fdm::{
    self,
};

#[derive(Deserialize)]

struct WaveEquationInput {
    width: usize,
    height: usize,
    c: f64,
    dx: f64,
    dy: f64,
    dt: f64,
    steps: usize,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_fdm_wave_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : WaveEquationInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(&FfiResult::<
                FdmGrid<f64>,
                String,
            >::err(
                "Invalid Bincode".to_string(),
            ))
        },
    };

    let result = physics_fdm::solve_wave_equation_2d(
        input.width,
        input.height,
        input.c,
        input.dx,
        input.dy,
        input.dt,
        input.steps,
        |x, y| {

            let dx_cen = x as f64 - (input.width / 2) as f64;

            let dy_cen = y as f64 - (input.height / 2) as f64;

            let dist2 = dx_cen.powi(2) + dy_cen.powi(2);

            (-dist2 / 20.0).exp()
        },
    );

    to_bincode_buffer(&FfiResult::<
        FdmGrid<f64>,
        String,
    >::ok(result))
}
