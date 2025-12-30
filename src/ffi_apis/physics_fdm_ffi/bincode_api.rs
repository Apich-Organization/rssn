//! Bincode-based FFI API for physics FDM functions.

use serde::Deserialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_fdm::FdmGrid;
use crate::physics::physics_fdm::FdmSolverConfig2D;
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

/// Solves the 2D wave equation using Finite Difference Method (FDM) via bincode serialization.
///
/// The wave equation ∂²u/∂t² = c²∇²u is solved using explicit finite differences with
/// a Gaussian initial condition centered at the grid midpoint.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `WaveEquationInput` with:
///   - `width`: Grid width (number of cells in x-direction)
///   - `height`: Grid height (number of cells in y-direction)
///   - `c`: Wave speed
///   - `dx`: Spatial step size in x-direction
///   - `dy`: Spatial step size in y-direction
///   - `dt`: Time step size
///   - `steps`: Number of time steps to simulate
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<FdmGrid<f64>, String>` with either:
/// - `ok`: Final wave field grid after time evolution
/// - `err`: Error message if computation failed
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

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

    let config = FdmSolverConfig2D {
        width: input.width,
        height: input.height,
        dx: input.dx,
        dy: input.dy,
        dt: input.dt,
        steps: input.steps,
    };

    let result = physics_fdm::solve_wave_equation_2d(
        &config,
        input.c,
        |x, y| {

            let dx_cen = x as f64 - (input.width / 2) as f64;

            let dy_cen = y as f64 - (input.height / 2) as f64;

            let dist2 = dy_cen.mul_add(dy_cen, dx_cen.powi(2));

            (-dist2 / 20.0).exp()
        },
    );

    to_bincode_buffer(&FfiResult::<
        FdmGrid<f64>,
        String,
    >::ok(result))
}
