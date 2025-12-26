//! Handle-based FFI API for physics EM (Euler Methods) functions.

use crate::numerical::matrix::Matrix;
use crate::physics::physics_em;

/// Simulates the oscillator forward Euler scenario and returns the results as a Matrix handle.
#[no_mangle]

pub extern "C" fn rssn_physics_em_simulate_oscillator_forward() -> *mut Matrix<f64> {

    let results = physics_em::simulate_oscillator_forward_euler_scenario();

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

    Box::into_raw(Box::new(
        Matrix::new(
            rows, cols, flattened,
        ),
    ))
}

/// Simulates the gravity semi-implicit Euler scenario and returns the results as a Matrix handle.
#[no_mangle]

pub extern "C" fn rssn_physics_em_simulate_gravity_semi_implicit() -> *mut Matrix<f64> {

    match physics_em::simulate_gravity_semi_implicit_euler_scenario() {
        | Ok(results) => {

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

            Box::into_raw(Box::new(
                Matrix::new(
                    rows, cols, flattened,
                ),
            ))
        },
        | Err(_) => std::ptr::null_mut(),
    }
}

/// Simulates the stiff decay scenario using backward Euler and returns a Matrix handle.
#[no_mangle]

pub extern "C" fn rssn_physics_em_simulate_stiff_decay_backward() -> *mut Matrix<f64> {

    match physics_em::simulate_stiff_decay_scenario() {
        | Ok(results) => {

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

            Box::into_raw(Box::new(
                Matrix::new(
                    rows, cols, flattened,
                ),
            ))
        },
        | Err(_) => std::ptr::null_mut(),
    }
}
