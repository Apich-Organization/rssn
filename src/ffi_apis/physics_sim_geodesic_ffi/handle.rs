//! Handle-based FFI API for physics sim geodesic relativity functions.

use crate::numerical::matrix::Matrix;
use crate::physics::physics_sim::geodesic_relativity::GeodesicParameters;
use crate::physics::physics_sim::geodesic_relativity::{
    self,
};

/// Runs a geodesic simulation and returns the resulting path as a Matrix handle (Nx2).
#[unsafe(no_mangle)]

pub extern "C" fn rssn_physics_sim_geodesic_run(
    black_hole_mass: f64,
    r0: f64,
    rdot0: f64,
    phi0: f64,
    phidot0: f64,
    proper_time_end: f64,
    initial_dt: f64,
) -> *mut Matrix<f64> {

    let params = GeodesicParameters {
        black_hole_mass,
        initial_state: [
            r0,
            rdot0,
            phi0,
            phidot0,
        ],
        proper_time_end,
        initial_dt,
    };

    let path = geodesic_relativity::run_geodesic_simulation(&params);

    let n = path.len();

    let mut data =
        Vec::with_capacity(n * 2);

    for (x, y) in path {

        data.push(x);

        data.push(y);
    }

    Box::into_raw(Box::new(
        Matrix::new(n, 2, data),
    ))
}
