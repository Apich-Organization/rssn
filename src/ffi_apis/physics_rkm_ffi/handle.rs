//! Handle-based FFI API for physics RKM (Runge-Kutta Methods) functions.

use crate::numerical::matrix::Matrix;
use crate::physics::physics_rkm;

/// Simulates the Lorenz attractor scenario and returns the results as a Matrix handle.
#[no_mangle]
pub extern "C" fn rssn_physics_rkm_simulate_lorenz() -> *mut Matrix<f64> {
    let results = physics_rkm::simulate_lorenz_attractor_scenario();
    let rows = results.len();
    if rows == 0 {
        return std::ptr::null_mut();
    }
    let cols = results[0].1.len() + 1; // +1 for time
    let mut flattened = Vec::with_capacity(rows * cols);
    for (t, y) in results {
        flattened.push(t);
        flattened.extend(y);
    }
    Box::into_raw(Box::new(Matrix::new(rows, cols, flattened)))
}

/// Simulates the damped oscillator scenario and returns the results as a Matrix handle.
#[no_mangle]
pub extern "C" fn rssn_physics_rkm_simulate_damped_oscillator() -> *mut Matrix<f64> {
    let results = physics_rkm::simulate_damped_oscillator_scenario();
    let rows = results.len();
    if rows == 0 {
        return std::ptr::null_mut();
    }
    let cols = results[0].1.len() + 1;
    let mut flattened = Vec::with_capacity(rows * cols);
    for (t, y) in results {
        flattened.push(t);
        flattened.extend(y);
    }
    Box::into_raw(Box::new(Matrix::new(rows, cols, flattened)))
}

/// Simulates the Van der Pol oscillator scenario and returns the results as a Matrix handle.
#[no_mangle]
pub extern "C" fn rssn_physics_rkm_simulate_vanderpol() -> *mut Matrix<f64> {
    let results = physics_rkm::simulate_vanderpol_scenario();
    let rows = results.len();
    if rows == 0 {
        return std::ptr::null_mut();
    }
    let cols = results[0].1.len() + 1;
    let mut flattened = Vec::with_capacity(rows * cols);
    for (t, y) in results {
        flattened.push(t);
        flattened.extend(y);
    }
    Box::into_raw(Box::new(Matrix::new(rows, cols, flattened)))
}

/// Simulates the Lotka-Volterra predator-prey scenario and returns the results as a Matrix handle.
#[no_mangle]
pub extern "C" fn rssn_physics_rkm_simulate_lotka_volterra() -> *mut Matrix<f64> {
    let results = physics_rkm::simulate_lotka_volterra_scenario();
    let rows = results.len();
    if rows == 0 {
        return std::ptr::null_mut();
    }
    let cols = results[0].1.len() + 1;
    let mut flattened = Vec::with_capacity(rows * cols);
    for (t, y) in results {
        flattened.push(t);
        flattened.extend(y);
    }
    Box::into_raw(Box::new(Matrix::new(rows, cols, flattened)))
}
