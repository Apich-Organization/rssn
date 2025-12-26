//! Handle-based FFI API for numerical fractal geometry and chaos functions.

use crate::numerical::fractal_geometry_and_chaos;

/// Computes the escape time for a single point in the Mandelbrot set.
#[no_mangle]

pub extern "C" fn rssn_num_fractal_mandelbrot_escape_time(
    c_real: f64,
    c_imag: f64,
    max_iter: u32,
) -> u32 {

    fractal_geometry_and_chaos::mandelbrot_escape_time(c_real, c_imag, max_iter)
}

/// Computes the escape time for a single point in a Julia set.
#[no_mangle]

pub extern "C" fn rssn_num_fractal_julia_escape_time(
    z_real: f64,
    z_imag: f64,
    c_real: f64,
    c_imag: f64,
    max_iter: u32,
) -> u32 {

    fractal_geometry_and_chaos::julia_escape_time(z_real, z_imag, c_real, c_imag, max_iter)
}

/// Computes the Lyapunov exponent for the logistic map.
#[no_mangle]

pub extern "C" fn rssn_num_fractal_lyapunov_logistic(
    r: f64,
    x0: f64,
    transient: usize,
    num_iterations: usize,
) -> f64 {

    fractal_geometry_and_chaos::lyapunov_exponent_logistic(r, x0, transient, num_iterations)
}

/// Computes the Lyapunov exponent for the Lorenz system.
#[no_mangle]

pub extern "C" fn rssn_num_fractal_lyapunov_lorenz(
    x0: f64,
    y0: f64,
    z0: f64,
    dt: f64,
    num_steps: usize,
    sigma: f64,
    rho: f64,
    beta: f64,
) -> f64 {

    fractal_geometry_and_chaos::lyapunov_exponent_lorenz(
        (x0, y0, z0),
        dt,
        num_steps,
        sigma,
        rho,
        beta,
    )
}

/// Generates Lorenz attractor points.
///
/// # Safety
/// `out_ptr` must be a valid pointer to at least `num_steps * 3` f64 values.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_lorenz_attractor(
    x0: f64,
    y0: f64,
    z0: f64,
    dt: f64,
    num_steps: usize,
    out_ptr: *mut f64,
) -> i32 {

    if out_ptr.is_null() {

        return -1;
    }

    let points = fractal_geometry_and_chaos::generate_lorenz_attractor((x0, y0, z0), dt, num_steps);

    for (i, (x, y, z)) in points.iter().enumerate() {

        *out_ptr.add(i * 3) = *x;

        *out_ptr.add(i * 3 + 1) = *y;

        *out_ptr.add(i * 3 + 2) = *z;
    }

    0
}

/// Generates Henon map points.
///
/// # Safety
/// `out_ptr` must be a valid pointer to at least `num_steps * 2` f64 values.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_henon_map(
    x0: f64,
    y0: f64,
    num_steps: usize,
    a: f64,
    b: f64,
    out_ptr: *mut f64,
) -> i32 {

    if out_ptr.is_null() {

        return -1;
    }

    let points = fractal_geometry_and_chaos::generate_henon_map((x0, y0), num_steps, a, b);

    for (i, (x, y)) in points.iter().enumerate() {

        *out_ptr.add(i * 2) = *x;

        *out_ptr.add(i * 2 + 1) = *y;
    }

    0
}

/// Iterates the logistic map.
///
/// # Safety
/// `out_ptr` must be a valid pointer to at least `num_steps + 1` f64 values.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_logistic_map(
    x0: f64,
    r: f64,
    num_steps: usize,
    out_ptr: *mut f64,
) -> i32 {

    if out_ptr.is_null() {

        return -1;
    }

    let orbit = fractal_geometry_and_chaos::logistic_map_iterate(x0, r, num_steps);

    for (i, &x) in orbit.iter().enumerate() {

        *out_ptr.add(i) = x;
    }

    0
}

/// Computes box-counting dimension.
///
/// # Safety
/// `points_ptr` must be a valid pointer to `num_points * 2` f64 values.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_box_counting_dim(
    points_ptr: *const f64,
    num_points: usize,
    num_scales: usize,
) -> f64 {

    if points_ptr.is_null() || num_points == 0 {

        return 0.0;
    }

    let mut points = Vec::with_capacity(num_points);

    for i in 0..num_points {

        let x = *points_ptr.add(i * 2);

        let y = *points_ptr.add(i * 2 + 1);

        points.push((x, y));
    }

    fractal_geometry_and_chaos::box_counting_dimension(&points, num_scales)
}

/// Computes correlation dimension.
///
/// # Safety
/// `points_ptr` must be a valid pointer to `num_points * 2` f64 values.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_correlation_dim(
    points_ptr: *const f64,
    num_points: usize,
    num_radii: usize,
) -> f64 {

    if points_ptr.is_null() || num_points == 0 {

        return 0.0;
    }

    let mut points = Vec::with_capacity(num_points);

    for i in 0..num_points {

        let x = *points_ptr.add(i * 2);

        let y = *points_ptr.add(i * 2 + 1);

        points.push((x, y));
    }

    fractal_geometry_and_chaos::correlation_dimension(&points, num_radii)
}
