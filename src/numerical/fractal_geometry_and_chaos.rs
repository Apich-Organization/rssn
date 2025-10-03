//! # Numerical Fractal Geometry and Chaos
//!
//! This module provides numerical tools for exploring fractal geometry and chaotic systems.
//! It includes functions for generating data for the Mandelbrot set and simulating
//! the Lorenz attractor, which are classic examples in these fields.

use num_complex::Complex;

/// Generates the data for a Mandelbrot set image.
///
/// The Mandelbrot set is a fractal defined by the iteration `z = z*z + c`.
/// Points `c` for which the iteration remains bounded form the set.
/// This function computes the escape time for each point in a given region.
///
/// # Arguments
/// * `width`, `height` - The dimensions of the output image.
/// * `x_range`, `y_range` - The region in the complex plane to plot.
/// * `max_iter` - The maximum number of iterations per point.
///
/// # Returns
/// A 2D vector where each element is the number of iterations it took for that point to escape.
pub fn generate_mandelbrot_set(
    width: usize,
    height: usize,
    x_range: (f64, f64),
    y_range: (f64, f64),
    max_iter: u32,
) -> Vec<Vec<u32>> {
    let mut data = vec![vec![0; width]; height];
    for (r, row) in data.iter_mut().enumerate() {
        for (c, val) in row.iter_mut().enumerate() {
            let x0 = x_range.0 + (x_range.1 - x_range.0) * (c as f64 / width as f64);
            let y0 = y_range.0 + (y_range.1 - y_range.0) * (r as f64 / height as f64);
            let mut z = Complex::new(0.0, 0.0);
            let c_val = Complex::new(x0, y0);
            let mut iter = 0;
            while z.norm_sqr() <= 4.0 && iter < max_iter {
                z = z * z + c_val;
                iter += 1;
            }
            *val = iter;
        }
    }
    data
}

/// Generates the points for a Lorenz attractor simulation.
///
/// The Lorenz attractor is a set of chaotic solutions for a simplified model of atmospheric
/// convection. This function numerically integrates the Lorenz system of differential equations
/// to produce a sequence of points that trace out the attractor.
///
/// # Arguments
/// * `start_point` - The initial `(x, y, z)` coordinates.
/// * `dt` - The time step for the integration.
/// * `num_steps` - The number of integration steps to perform.
///
/// # Returns
/// A `Vec` of `(f64, f64, f64)` tuples representing the trajectory of the attractor.
pub fn generate_lorenz_attractor(
    start_point: (f64, f64, f64),
    dt: f64,
    num_steps: usize,
) -> Vec<(f64, f64, f64)> {
    let (sigma, rho, beta) = (10.0, 28.0, 8.0 / 3.0);
    let mut points = Vec::with_capacity(num_steps);
    let (mut x, mut y, mut z) = start_point;

    for _ in 0..num_steps {
        let dx = sigma * (y - x);
        let dy = x * (rho - z) - y;
        let dz = x * y - beta * z;

        x += dx * dt;
        y += dy * dt;
        z += dz * dt;

        points.push((x, y, z));
    }
    points
}
