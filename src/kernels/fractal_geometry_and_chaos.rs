//! # Numerical Fractal Geometry and Chaos Theory
//!
//! This module provides numerical tools for exploring fractal geometry and chaotic systems.
//! It includes functions for generating fractal images, simulating famous attractors,
//! computing fractal dimensions, and analyzing chaotic dynamics through Lyapunov exponents.
//!
//! ## Features
//!
//! ### Fractal Generation
//! - **Mandelbrot Set**: `generate_mandelbrot_set` - Classic escape-time fractal
//! - **Julia Set**: `generate_julia_set` - Related to Mandelbrot with fixed parameter c
//! - **Burning Ship**: `generate_burning_ship` - Variant with absolute value operations
//! - **Multibrot**: `generate_multibrot` - Generalized Mandelbrot with power d
//! - **Newton Fractal**: `generate_newton_fractal` - Fractal from Newton's method
//!
//! ### Strange Attractors
//! - **Lorenz Attractor**: `generate_lorenz_attractor` - Atmospheric convection model
//! - **Rossler Attractor**: `generate_rossler_attractor` - Simpler chaotic flow
//! - **Henon Map**: `generate_henon_map` - Discrete chaotic map
//! - **Tinkerbell Map**: `generate_tinkerbell_map` - 2D discrete dynamical system
//!
//! ### Chaotic Maps
//! - **Logistic Map**: `logistic_map_iterate` - Population dynamics model
//! - **Bifurcation Diagram**: `logistic_bifurcation` - Visualize route to chaos
//!
//! ### Dimension Estimation
//! - **Box-counting Dimension**: `box_counting_dimension` - Estimate fractal dimension
//! - **Correlation Dimension**: `correlation_dimension` - Statistical dimension estimate
//!
//! ### Chaos Analysis
//! - **Lyapunov Exponent**: `lyapunov_exponent_logistic` - Measure of chaos
//! - **Orbit Analysis**: `orbit_density` - Distribution of orbits in phase space
//!
//! ## Examples
//!
//! ### Generating a Mandelbrot Set
//! ```
//! use rssn::kernels::fractal_geometry_and_chaos::generate_mandelbrot_set;
//!
//! let data = generate_mandelbrot_set(100, 100, (-2.0, 1.0), (-1.5, 1.5), 100);
//!
//! assert_eq!(data.len(), 100);
//!
//! assert_eq!(data[0].len(), 100);
//! ```
//!
//! ### Simulating the Lorenz Attractor
//! ```
//! use rssn::kernels::fractal_geometry_and_chaos::generate_lorenz_attractor;
//!
//! let trajectory = generate_lorenz_attractor((1.0, 1.0, 1.0), 0.01, 1000);
//!
//! assert_eq!(trajectory.len(), 1000);
//! ```

use num_complex::Complex;
use serde::Deserialize;
use serde::Serialize;

// ============================================================================
// Common Types
// ============================================================================

/// A 2D point in the plane
pub type Point2D = (f64, f64);

/// A 3D point in space
pub type Point3D = (f64, f64, f64);

/// Result of a fractal computation with escape time data
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FractalData {
    /// Width of the image
    pub width: usize,
    /// Height of the image
    pub height: usize,
    /// Escape time data (row-major order)
    pub data: Vec<u32>,
    /// Maximum iterations used
    pub max_iter: u32,
}

impl FractalData {
    /// Creates a new `FractalData` with given dimensions
    #[must_use]
    pub fn new(
        width: usize,
        height: usize,
        max_iter: u32,
    ) -> Self {
        Self {
            width,
            height,
            data: vec![0; width * height],
            max_iter,
        }
    }

    /// Gets the escape time at a specific pixel
    #[must_use]
    pub fn get(
        &self,
        x: usize,
        y: usize,
    ) -> Option<u32> {
        if x < self.width && y < self.height {
            Some(self.data[y * self.width + x])
        } else {
            None
        }
    }

    /// Sets the escape time at a specific pixel
    pub fn set(
        &mut self,
        x: usize,
        y: usize,
        value: u32,
    ) {
        if x < self.width && y < self.height {
            self.data[y * self.width + x] = value;
        }
    }
}

// ============================================================================
// Mandelbrot Set
// ============================================================================

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
///
/// # Example
/// ```
/// use rssn::kernels::fractal_geometry_and_chaos::generate_mandelbrot_set;
///
/// let data = generate_mandelbrot_set(100, 100, (-2.5, 1.0), (-1.5, 1.5), 50);
///
/// assert_eq!(data.len(), 100);
/// ```
#[must_use]
pub fn generate_mandelbrot_set(
    width: usize,
    height: usize,
    x_range: (f64, f64),
    y_range: (f64, f64),
    max_iter: u32,
) -> Vec<Vec<u32>> {
    let mut data = vec![vec![0; width]; height];

    let x_scale = (x_range.1 - x_range.0) / width as f64;

    let y_scale = (y_range.1 - y_range.0) / height as f64;

    for (r, row) in data.iter_mut().enumerate() {
        for (c, val) in row.iter_mut().enumerate() {
            let x0 = (c as f64).mul_add(x_scale, x_range.0);

            let y0 = (r as f64).mul_add(y_scale, y_range.0);

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

/// Computes the escape time for a single point in the Mandelbrot set.
///
/// # Arguments
/// * `c_real` - Real part of the complex number c
/// * `c_imag` - Imaginary part of the complex number c
/// * `max_iter` - Maximum iterations
///
/// # Returns
/// The number of iterations before escape, or `max_iter` if the point is in the set.
#[must_use]
pub fn mandelbrot_escape_time(
    c_real: f64,
    c_imag: f64,
    max_iter: u32,
) -> u32 {
    let c = Complex::new(c_real, c_imag);

    let mut z = Complex::new(0.0, 0.0);

    let mut iter = 0;

    while z.norm_sqr() <= 4.0 && iter < max_iter {
        z = z * z + c;

        iter += 1;
    }

    iter
}

// ============================================================================
// Julia Set
// ============================================================================

/// Generates the data for a Julia set image.
///
/// The Julia set is defined by iterating `z = z*z + c` for a fixed c parameter.
/// Unlike the Mandelbrot set, c is constant and z varies.
///
/// # Arguments
/// * `width`, `height` - The dimensions of the output image.
/// * `x_range`, `y_range` - The region in the complex plane to plot.
/// * `c` - The constant c parameter `(real, imaginary)`.
/// * `max_iter` - The maximum number of iterations per point.
///
/// # Returns
/// A 2D vector where each element is the escape time.
///
/// # Example
/// ```
/// use rssn::kernels::fractal_geometry_and_chaos::generate_julia_set;
///
/// // Julia set for c = -0.4 + 0.6i
/// let data = generate_julia_set(100, 100, (-2.0, 2.0), (-2.0, 2.0), (-0.4, 0.6), 50);
///
/// assert_eq!(data.len(), 100);
/// ```
#[must_use]
pub fn generate_julia_set(
    width: usize,
    height: usize,
    x_range: (f64, f64),
    y_range: (f64, f64),
    c: (f64, f64),
    max_iter: u32,
) -> Vec<Vec<u32>> {
    let mut data = vec![vec![0; width]; height];

    let c_val = Complex::new(c.0, c.1);

    let x_scale = (x_range.1 - x_range.0) / width as f64;

    let y_scale = (y_range.1 - y_range.0) / height as f64;

    for (r, row) in data.iter_mut().enumerate() {
        for (col, val) in row.iter_mut().enumerate() {
            let x0 = (col as f64).mul_add(x_scale, x_range.0);

            let y0 = (r as f64).mul_add(y_scale, y_range.0);

            let mut z = Complex::new(x0, y0);

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

/// Computes the escape time for a single point in a Julia set.
#[must_use]
pub fn julia_escape_time(
    z_real: f64,
    z_imag: f64,
    c_real: f64,
    c_imag: f64,
    max_iter: u32,
) -> u32 {
    let c = Complex::new(c_real, c_imag);

    let mut z = Complex::new(z_real, z_imag);

    let mut iter = 0;

    while z.norm_sqr() <= 4.0 && iter < max_iter {
        z = z * z + c;

        iter += 1;
    }

    iter
}

// ============================================================================
// Burning Ship Fractal
// ============================================================================

/// Generates the Burning Ship fractal.
///
/// Similar to the Mandelbrot set but uses absolute values:
/// `z = (|Re(z)| + i|Im(z)|)^2 + c`
///
/// The resulting fractal resembles a burning ship.
///
/// # Arguments
/// * `width`, `height` - The dimensions of the output image.
/// * `x_range`, `y_range` - The region in the complex plane to plot.
/// * `max_iter` - The maximum number of iterations per point.
///
/// # Returns
/// A 2D vector where each element is the escape time.
#[must_use]
pub fn generate_burning_ship(
    width: usize,
    height: usize,
    x_range: (f64, f64),
    y_range: (f64, f64),
    max_iter: u32,
) -> Vec<Vec<u32>> {
    let mut data = vec![vec![0; width]; height];

    let x_scale = (x_range.1 - x_range.0) / width as f64;

    let y_scale = (y_range.1 - y_range.0) / height as f64;

    for (r, row) in data.iter_mut().enumerate() {
        for (col, val) in row.iter_mut().enumerate() {
            let cx = (col as f64).mul_add(x_scale, x_range.0);

            let cy = (r as f64).mul_add(y_scale, y_range.0);

            let mut zx = 0.0;

            let mut zy = 0.0;

            let mut iter = 0;

            while zx * zx + zy * zy <= 4.0 && iter < max_iter {
                let xtemp = zx * zx - zy * zy + cx;

                zy = (2.0 * zx.abs()).mul_add(zy.abs(), cy);

                zx = xtemp;

                iter += 1;
            }

            *val = iter;
        }
    }

    data
}

// ============================================================================
// Multibrot Set
// ============================================================================

/// Generates the Multibrot set with exponent d.
///
/// The Multibrot set is a generalization where `z = z^d + c`.
/// d=2 gives the standard Mandelbrot set.
///
/// # Arguments
/// * `width`, `height` - The dimensions of the output image.
/// * `x_range`, `y_range` - The region in the complex plane to plot.
/// * `d` - The exponent (power).
/// * `max_iter` - The maximum number of iterations per point.
///
/// # Returns
/// A 2D vector where each element is the escape time.
#[must_use]
pub fn generate_multibrot(
    width: usize,
    height: usize,
    x_range: (f64, f64),
    y_range: (f64, f64),
    d: f64,
    max_iter: u32,
) -> Vec<Vec<u32>> {
    let mut data = vec![vec![0; width]; height];

    let x_scale = (x_range.1 - x_range.0) / width as f64;

    let y_scale = (y_range.1 - y_range.0) / height as f64;

    for (r, row) in data.iter_mut().enumerate() {
        for (col, val) in row.iter_mut().enumerate() {
            let x0 = (col as f64).mul_add(x_scale, x_range.0);

            let y0 = (r as f64).mul_add(y_scale, y_range.0);

            let mut z = Complex::new(0.0, 0.0);

            let c_val = Complex::new(x0, y0);

            let mut iter = 0;

            // Escape radius depends on d
            let escape_radius = 2.0_f64.max(c_val.norm().powf(1.0 / (d - 1.0)));

            while z.norm() <= escape_radius && iter < max_iter {
                z = z.powf(d) + c_val;

                iter += 1;
            }

            *val = iter;
        }
    }

    data
}

// ============================================================================
// Newton Fractal
// ============================================================================

/// Generates a Newton fractal for `z^3 - 1 = 0`.
///
/// Newton's method applied to finding roots of `f(z) = z^3 - 1` produces
/// a fractal pattern. Each pixel is colored based on which root it converges to.
///
/// # Arguments
/// * `width`, `height` - The dimensions of the output image.
/// * `x_range`, `y_range` - The region in the complex plane to plot.
/// * `max_iter` - The maximum number of iterations per point.
/// * `tolerance` - The convergence tolerance.
///
/// # Returns
/// A 2D vector where each element indicates which root (0, 1, or 2) was reached,
/// or 3 if no convergence.
#[must_use]
pub fn generate_newton_fractal(
    width: usize,
    height: usize,
    x_range: (f64, f64),
    y_range: (f64, f64),
    max_iter: u32,
    tolerance: f64,
) -> Vec<Vec<u32>> {
    let mut data = vec![vec![0; width]; height];

    // Roots of z^3 - 1 = 0
    let roots = [
        Complex::new(1.0, 0.0),
        Complex::new(-0.5, 3.0_f64.sqrt() / 2.0),
        Complex::new(-0.5, -3.0_f64.sqrt() / 2.0),
    ];

    let x_scale = (x_range.1 - x_range.0) / width as f64;

    let y_scale = (y_range.1 - y_range.0) / height as f64;

    for (r, row) in data.iter_mut().enumerate() {
        for (col, val) in row.iter_mut().enumerate() {
            let x0 = (col as f64).mul_add(x_scale, x_range.0);

            let y0 = (r as f64).mul_add(y_scale, y_range.0);

            let mut z = Complex::new(x0, y0);

            let mut root_index = 3u32; // No root found
            for _ in 0..max_iter {
                // f(z) = z^3 - 1
                // f'(z) = 3z^2
                let z2 = z * z;

                let z3 = z2 * z;

                let f = z3 - Complex::new(1.0, 0.0);

                let df = Complex::new(3.0, 0.0) * z2;

                if df.norm() < 1e-10 {
                    break;
                }

                z -= f / df;

                // Check if close to any root
                for (i, &root) in roots.iter().enumerate() {
                    if (z - root).norm() < tolerance {
                        root_index = i as u32;

                        break;
                    }
                }

                if root_index < 3 {
                    break;
                }
            }

            *val = root_index;
        }
    }

    data
}

// ============================================================================
// Lorenz Attractor
// ============================================================================

/// Generates the points for a Lorenz attractor simulation.
///
/// The Lorenz attractor is a set of chaotic solutions for a simplified model of atmospheric
/// convection. This function numerically integrates the Lorenz system of differential equations
/// to produce a sequence of points that trace out the attractor.
///
/// The system is defined by:
/// - dx/dt = σ(y - x)
/// - dy/dt = x(ρ - z) - y
/// - dz/dt = xy - βz
///
/// Uses standard parameters: σ = 10, ρ = 28, β = 8/3
///
/// # Arguments
/// * `start_point` - The initial `(x, y, z)` coordinates.
/// * `dt` - The time step for the integration.
/// * `num_steps` - The number of integration steps to perform.
///
/// # Returns
/// A `Vec` of `(f64, f64, f64)` tuples representing the trajectory of the attractor.
///
/// # Example
/// ```
/// use rssn::kernels::fractal_geometry_and_chaos::generate_lorenz_attractor;
///
/// let points = generate_lorenz_attractor((0.1, 0.0, 0.0), 0.01, 1000);
///
/// assert_eq!(points.len(), 1000);
/// ```
#[must_use]
pub fn generate_lorenz_attractor(
    start_point: (f64, f64, f64),
    dt: f64,
    num_steps: usize,
) -> Vec<(f64, f64, f64)> {
    generate_lorenz_attractor_custom(start_point, dt, num_steps, 10.0, 28.0, 8.0 / 3.0)
}

/// Generates a Lorenz attractor with custom parameters.
///
/// # Arguments
/// * `start_point` - The initial `(x, y, z)` coordinates.
/// * `dt` - The time step for the integration.
/// * `num_steps` - The number of integration steps.
/// * `sigma`, `rho`, `beta` - The Lorenz system parameters.
///
/// # Returns
/// A `Vec` of `(f64, f64, f64)` tuples representing the trajectory.
#[must_use]
pub fn generate_lorenz_attractor_custom(
    start_point: (f64, f64, f64),
    dt: f64,
    num_steps: usize,
    sigma: f64,
    rho: f64,
    beta: f64,
) -> Vec<(f64, f64, f64)> {
    let mut points = Vec::with_capacity(num_steps);

    let (mut x, mut y, mut z) = start_point;

    for _ in 0..num_steps {
        let dx = sigma * (y - x);

        let dy = x.mul_add(rho - z, -y);

        let dz = x.mul_add(y, -(beta * z));

        x += dx * dt;

        y += dy * dt;

        z += dz * dt;

        points.push((x, y, z));
    }

    points
}

// ============================================================================
// Rossler Attractor
// ============================================================================

/// Generates the Rossler attractor.
///
/// The Rossler attractor is a system of three non-linear ordinary differential equations
/// that exhibits chaotic dynamics. It is simpler than the Lorenz system.
///
/// The system is defined by:
/// - dx/dt = -y - z
/// - dy/dt = x + ay
/// - dz/dt = b + z(x - c)
///
/// Common parameter values: a = 0.2, b = 0.2, c = 5.7
///
/// # Arguments
/// * `start_point` - The initial `(x, y, z)` coordinates.
/// * `dt` - The time step for the integration.
/// * `num_steps` - The number of integration steps.
/// * `a`, `b`, `c` - The Rossler system parameters.
///
/// # Returns
/// A `Vec` of `(f64, f64, f64)` tuples representing the trajectory.
#[must_use]
pub fn generate_rossler_attractor(
    start_point: (f64, f64, f64),
    dt: f64,
    num_steps: usize,
    a: f64,
    b: f64,
    c: f64,
) -> Vec<(f64, f64, f64)> {
    let mut points = Vec::with_capacity(num_steps);

    let (mut x, mut y, mut z) = start_point;

    for _ in 0..num_steps {
        let dx = -y - z;

        let dy = a.mul_add(y, x);

        let dz = z.mul_add(x - c, b);

        x += dx * dt;

        y += dy * dt;

        z += dz * dt;

        points.push((x, y, z));
    }

    points
}

// ============================================================================
// Henon Map
// ============================================================================

/// Generates the Henon map trajectory.
///
/// The Henon map is a discrete-time dynamical system that exhibits chaotic behavior.
///
/// The map is:
/// - x_{n+1} = 1 - a*`x_n^2` + `y_n`
/// - y_{n+1} = b*`x_n`
///
/// Classic parameter values: a = 1.4, b = 0.3
///
/// # Arguments
/// * `start_point` - The initial `(x, y)` coordinates.
/// * `num_steps` - The number of iterations.
/// * `a`, `b` - The Henon map parameters.
///
/// # Returns
/// A `Vec` of `(f64, f64)` tuples representing the trajectory.
#[must_use]
pub fn generate_henon_map(
    start_point: (f64, f64),
    num_steps: usize,
    a: f64,
    b: f64,
) -> Vec<(f64, f64)> {
    let mut points = Vec::with_capacity(num_steps);

    let (mut x, mut y) = start_point;

    for _ in 0..num_steps {
        let x_new = (a * x).mul_add(-x, 1.0) + y;

        let y_new = b * x;

        x = x_new;

        y = y_new;

        points.push((x, y));
    }

    points
}

// ============================================================================
// Tinkerbell Map
// ============================================================================

/// Generates the Tinkerbell map trajectory.
///
/// The Tinkerbell map is a two-dimensional iterated map:
/// - x_{n+1} = `x_n^2` - `y_n^2` + `ax_n` + `by_n`
/// - y_{n+1} = `2x_n`*`y_n` + `cx_n` + `dy_n`
///
/// Classic parameter values: a = 0.9, b = -0.6013, c = 2.0, d = 0.5
///
/// # Arguments
/// * `start_point` - The initial `(x, y)` coordinates.
/// * `num_steps` - The number of iterations.
/// * `a`, `b`, `c`, `d` - The map parameters.
///
/// # Returns
/// A `Vec` of `(f64, f64)` tuples representing the trajectory.
#[must_use]
pub fn generate_tinkerbell_map(
    start_point: (f64, f64),
    num_steps: usize,
    a: f64,
    b: f64,
    c: f64,
    d: f64,
) -> Vec<(f64, f64)> {
    let mut points = Vec::with_capacity(num_steps);

    let (mut x, mut y) = start_point;

    for _ in 0..num_steps {
        let x_new = b.mul_add(y, a.mul_add(x, x.mul_add(x, -(y * y))));

        let y_new = d.mul_add(y, (2.0 * x).mul_add(y, c * x));

        x = x_new;

        y = y_new;

        points.push((x, y));
    }

    points
}

// ============================================================================
// Logistic Map
// ============================================================================

/// Iterates the logistic map.
///
/// The logistic map is defined by: x_{n+1} = r * `x_n` * (1 - `x_n`)
///
/// It is a classic example of how complex behavior can arise from simple nonlinear
/// dynamical equations. For r between 0 and 4, the dynamics range from stable
/// fixed points to period doubling to chaos.
///
/// # Arguments
/// * `x0` - Initial value (should be between 0 and 1).
/// * `r` - The growth rate parameter.
/// * `num_steps` - The number of iterations.
///
/// # Returns
/// A `Vec` of x values representing the orbit.
#[must_use]
pub fn logistic_map_iterate(
    x0: f64,
    r: f64,
    num_steps: usize,
) -> Vec<f64> {
    let mut orbit = Vec::with_capacity(num_steps + 1);

    let mut x = x0;

    orbit.push(x);

    for _ in 0..num_steps {
        x = r * x * (1.0 - x);

        orbit.push(x);
    }

    orbit
}

/// Generates data for a bifurcation diagram of the logistic map.
///
/// For each value of r, the map is iterated and the final values after
/// discarding transients are recorded.
///
/// # Arguments
/// * `r_range` - The range of r values `(r_min, r_max)`.
/// * `num_r_values` - Number of r values to sample.
/// * `transient` - Number of initial iterations to discard.
/// * `num_points` - Number of points to record after transient.
/// * `x0` - Initial x value.
///
/// # Returns
/// A `Vec` of `(r, x)` pairs for plotting the bifurcation diagram.
#[must_use]
pub fn logistic_bifurcation(
    r_range: (f64, f64),
    num_r_values: usize,
    transient: usize,
    num_points: usize,
    x0: f64,
) -> Vec<(f64, f64)> {
    let mut data = Vec::with_capacity(num_r_values * num_points);

    let r_step = (r_range.1 - r_range.0) / (num_r_values as f64 - 1.0);

    for i in 0..num_r_values {
        let r = (i as f64).mul_add(r_step, r_range.0);

        let mut x = x0;

        // Discard transient
        for _ in 0..transient {
            x = r * x * (1.0 - x);
        }

        // Record points
        for _ in 0..num_points {
            x = r * x * (1.0 - x);

            data.push((r, x));
        }
    }

    data
}

// ============================================================================
// Lyapunov Exponent
// ============================================================================

/// Computes the Lyapunov exponent for the logistic map.
///
/// The Lyapunov exponent λ measures the rate of separation of infinitesimally
/// close trajectories. For the logistic map:
/// λ = lim (1/n) * sum(ln|r(1-2x_i)|)
///
/// Positive λ indicates chaos.
///
/// # Arguments
/// * `r` - The growth rate parameter.
/// * `x0` - Initial value.
/// * `transient` - Number of iterations to discard.
/// * `num_iterations` - Number of iterations for averaging.
///
/// # Returns
/// The estimated Lyapunov exponent.
#[must_use]
pub fn lyapunov_exponent_logistic(
    r: f64,
    x0: f64,
    transient: usize,
    num_iterations: usize,
) -> f64 {
    let mut x = x0;

    // Discard transient
    for _ in 0..transient {
        x = r * x * (1.0 - x);
    }

    // Compute average
    let mut sum = 0.0;

    for _ in 0..num_iterations {
        // Derivative of f(x) = rx(1-x) is f'(x) = r(1-2x)
        let deriv = r * 2.0f64.mul_add(-x, 1.0);

        if deriv.abs() > 0.0 {
            sum += deriv.abs().ln();
        }

        x = r * x * (1.0 - x);
    }

    sum / num_iterations as f64
}

/// Computes the largest Lyapunov exponent for a 3D flow (e.g., Lorenz system).
///
/// Uses the method of tracking the evolution of a small perturbation.
///
/// # Arguments
/// * `start_point` - Initial point `(x, y, z)`.
/// * `dt` - Time step.
/// * `num_steps` - Number of integration steps.
/// * `sigma`, `rho`, `beta` - Lorenz system parameters.
///
/// # Returns
/// The estimated largest Lyapunov exponent.
#[must_use]
pub fn lyapunov_exponent_lorenz(
    start_point: (f64, f64, f64),
    dt: f64,
    num_steps: usize,
    sigma: f64,
    rho: f64,
    beta: f64,
) -> f64 {
    let (mut x, mut y, mut z) = start_point;

    let mut sum = 0.0;

    // Initial perturbation vector (normalized)
    let eps = 1e-8;

    let mut dx = eps;

    let mut dy = 0.0;

    let mut dz = 0.0;

    for _ in 0..num_steps {
        // Evolve the main trajectory
        let dx_main = sigma * (y - x);

        let dy_main = x.mul_add(rho - z, -y);

        let dz_main = x.mul_add(y, -(beta * z));

        x += dx_main * dt;

        y += dy_main * dt;

        z += dz_main * dt;

        // Evolve the perturbation using linearized equations
        let ddx = sigma * (dy - dx);

        let ddy = dx * (rho - z) - x * dz - dy;

        let ddz = beta.mul_add(-dz, dx * y + x * dy);

        dx += ddx * dt;

        dy += ddy * dt;

        dz += ddz * dt;

        // Compute and accumulate the growth rate
        let norm = (dx * dx + dy * dy + dz * dz).sqrt();

        if norm > 0.0 {
            sum += norm.ln();

            // Renormalize
            dx /= norm;

            dy /= norm;

            dz /= norm;

            dx *= eps;

            dy *= eps;

            dz *= eps;
        }
    }

    sum / (num_steps as f64 * dt) - eps.ln() / dt
}

// ============================================================================
// Dimension Estimation
// ============================================================================

/// Estimates the box-counting (Minkowski) dimension from a set of 2D points.
///
/// The box-counting dimension is estimated by:
/// D = lim(ε→0) [log(N(ε)) / log(1/ε)]
///
/// where N(ε) is the number of boxes of side ε needed to cover the set.
///
/// # Arguments
/// * `points` - A slice of 2D points.
/// * `num_scales` - Number of different box sizes to use.
///
/// # Returns
/// An estimate of the box-counting dimension.
#[must_use]
pub fn box_counting_dimension(
    points: &[(f64, f64)],
    num_scales: usize,
) -> f64 {
    if points.is_empty() || num_scales < 2 {
        return 0.0;
    }

    // Find bounding box
    let mut x_min = f64::MAX;

    let mut x_max = f64::MIN;

    let mut y_min = f64::MAX;

    let mut y_max = f64::MIN;

    for &(x, y) in points {
        x_min = x_min.min(x);

        x_max = x_max.max(x);

        y_min = y_min.min(y);

        y_max = y_max.max(y);
    }

    let extent = (x_max - x_min).max(y_max - y_min);

    if extent <= 0.0 {
        return 0.0;
    }

    // Collect (log(1/eps), log(N)) pairs
    let mut log_inv_eps = Vec::with_capacity(num_scales);

    let mut log_n = Vec::with_capacity(num_scales);

    for scale_idx in 0..num_scales {
        let num_boxes = 1 << (scale_idx + 1); // 2, 4, 8, 16, ...
        let eps = extent / f64::from(num_boxes);

        // Count occupied boxes using a hash set
        let mut occupied = std::collections::HashSet::new();

        for &(x, y) in points {
            let bx = ((x - x_min) / eps).floor() as i64;

            let by = ((y - y_min) / eps).floor() as i64;

            occupied.insert((bx, by));
        }

        let n = occupied.len();

        if n > 0 && eps > 0.0 {
            log_inv_eps.push((1.0 / eps).ln());

            log_n.push((n as f64).ln());
        }
    }

    // Linear regression to find slope
    if log_inv_eps.len() < 2 {
        return 0.0;
    }

    linear_regression_slope(&log_inv_eps, &log_n).unwrap_or(0.0)
}

/// Estimates the correlation dimension from the log-log slope of the fraction
/// of pairs closer than each of `num_radii` logarithmically spaced radii.
///
/// Returns zero when fewer than two scales can determine a slope, and NaN for
/// non-finite coordinates or distances outside the range of `f64`.
#[must_use]
pub fn correlation_dimension(
    points: &[Point2D],
    num_radii: usize,
) -> f64 {
    if points.len() < 2 || num_radii < 2 {
        return 0.0;
    }
    if !valid_points(points) {
        return f64::NAN;
    }

    let ranges = map_pair_blocks(points.len(), |rows| {
        let mut min = f64::INFINITY;
        let mut max: f64 = 0.0;
        for distance in pair_distances(points, rows) {
            if distance > 0.0 {
                min = min.min(distance);
                max = max.max(distance);
            }
        }
        (min, max)
    });
    let (min, max) = ranges.iter().fold((f64::INFINITY, 0.0_f64), |a, b| {
        (a.0.min(b.0), a.1.max(b.1))
    });
    if !max.is_finite() {
        return f64::NAN;
    }
    if max <= min {
        return 0.0;
    }

    let mut radii = Vec::new();
    if radii.try_reserve_exact(num_radii).is_err() {
        return f64::NAN;
    }
    let (start, end) = (min.ln(), max.ln());
    radii.push(min);
    for i in 1..num_radii - 1 {
        let log_radius = start + (end - start) * (i as f64 / (num_radii - 1) as f64);
        radii.push(log_radius.exp().clamp(min, max));
    }
    radii.push(max);
    correlation_slope(points, &radii).unwrap_or(0.0)
}

/// Estimates correlation dimension from caller-supplied positive, finite radii.
///
/// Radii can be supplied in any order. At radius `r`, only pairs strictly
/// closer than `r` contribute. At least two distinct radii with close pairs
/// are needed for a fit.
///
/// Returns `None` for non-finite coordinates or when the radii cannot determine
/// a slope.
#[must_use]
pub fn correlation_dimension_with_radii(
    points: &[Point2D],
    radii: &[f64],
) -> Option<f64> {
    if points.len() < 2
        || radii.len() < 2
        || !valid_points(points)
        || radii.iter().any(|&r| !r.is_finite() || r <= 0.0)
    {
        return None;
    }
    let mut radii = radii.to_vec();
    radii.sort_unstable_by(f64::total_cmp);
    correlation_slope(points, &radii)
}

fn valid_points(points: &[Point2D]) -> bool {
    points.iter().all(|(x, y)| x.is_finite() && y.is_finite())
}

fn map_pair_blocks<T: Send>(
    num_points: usize,
    f: impl Fn(std::ops::Range<usize>) -> T + Sync,
) -> Vec<T> {
    // Bound the number of scratch histograms independently of the point count
    // and of the number of worker threads.
    let block = if num_points >= 256 {
        num_points.div_ceil(16)
    } else {
        num_points
    };
    crate::kernels::par::map_indexed(num_points.div_ceil(block), 2, |i| {
        let start = i * block;
        f(start..(start + block).min(num_points))
    })
}

fn pair_distances(
    points: &[Point2D],
    rows: std::ops::Range<usize>,
) -> impl Iterator<Item = f64> {
    rows.flat_map(|j| {
        let (x, y) = points[j];
        points[j + 1..]
            .iter()
            .map(move |&(u, v)| (x - u).hypot(y - v))
    })
}

fn correlation_slope(
    points: &[Point2D],
    radii: &[f64],
) -> Option<f64> {
    let counts = map_pair_blocks(points.len(), |rows| {
        let mut histogram = vec![0_u128; radii.len()];
        for distance in pair_distances(points, rows) {
            // A pair contributes to the suffix of radii strictly above it.
            let first = radii.partition_point(|&r| r <= distance);
            if let Some(bucket) = histogram.get_mut(first) {
                *bucket += 1;
            }
        }
        histogram
    });
    let mut counts = counts.into_iter();
    let mut totals = counts.next()?;
    for histogram in counts {
        for (total, count) in totals.iter_mut().zip(histogram) {
            *total += count;
        }
    }
    // The normalising number of pairs is constant across all radii. Its
    // logarithm shifts the intercept and leaves the fitted slope unchanged.
    let mut cumulative = 0;
    let mut log_radii = Vec::with_capacity(radii.len());
    let mut log_counts = Vec::with_capacity(radii.len());
    for (&radius, count) in radii.iter().zip(totals) {
        cumulative += count;
        if cumulative > 0 {
            log_radii.push(radius.ln());
            log_counts.push((cumulative as f64).ln());
        }
    }
    linear_regression_slope(&log_radii, &log_counts)
}

fn linear_regression_slope(
    x: &[f64],
    y: &[f64],
) -> Option<f64> {
    if x.len() < 2 {
        return None;
    }
    let n = x.len() as f64;
    let mean_x = x.iter().sum::<f64>() / n;
    let mean_y = y.iter().sum::<f64>() / n;
    let mut variance = 0.0;
    let mut covariance = 0.0;
    for (&a, &b) in x.iter().zip(y) {
        let dx = a - mean_x;
        variance += dx * dx;
        covariance += dx * (b - mean_y);
    }
    // The subtraction from the mean can introduce residual variance even
    // when every x is the same, so also test for an actual distinct log scale.
    (variance > 0.0
        && x.windows(2)
            .any(|pair| !crate::kernels::exact_eq(pair[0], pair[1])))
    .then(|| covariance / variance)
    .filter(|slope| slope.is_finite())
}

// ============================================================================
// Orbit Analysis
// ============================================================================

/// Computes the density distribution of orbit points in a 2D histogram.
///
/// # Arguments
/// * `points` - The orbit trajectory.
/// * `x_bins` - Number of bins in x direction.
/// * `y_bins` - Number of bins in y direction.
/// * `x_range` - Range for x axis.
/// * `y_range` - Range for y axis.
///
/// # Returns
/// A 2D vector of counts for each bin.
#[must_use]
pub fn orbit_density(
    points: &[(f64, f64)],
    x_bins: usize,
    y_bins: usize,
    x_range: (f64, f64),
    y_range: (f64, f64),
) -> Vec<Vec<usize>> {
    let mut density = vec![vec![0usize; x_bins]; y_bins];

    let x_scale = (x_range.1 - x_range.0) / x_bins as f64;

    let y_scale = (y_range.1 - y_range.0) / y_bins as f64;

    for &(x, y) in points {
        if x >= x_range.0 && x < x_range.1 && y >= y_range.0 && y < y_range.1 {
            let xi: usize = (((x - x_range.0) / x_scale) as i64).try_into().unwrap_or(0);

            let yi: usize = (((y - y_range.0) / y_scale) as i64).try_into().unwrap_or(0);

            if xi < x_bins && yi < y_bins {
                density[yi][xi] += 1;
            }
        }
    }

    density
}

/// Computes the entropy of an orbit from its density distribution.
///
/// Uses Shannon entropy: H = -`sum(p_i` * `log(p_i)`)
///
/// # Arguments
/// * `density` - A 2D density histogram from `orbit_density`.
///
/// # Returns
/// The Shannon entropy of the distribution.
#[must_use]
pub fn orbit_entropy(density: &[Vec<usize>]) -> f64 {
    let total: usize = density.iter().flat_map(|row| row.iter()).sum();

    if total == 0 {
        return 0.0;
    }

    let mut entropy = 0.0;

    let total_f = total as f64;

    for row in density {
        for &count in row {
            if count > 0 {
                let p = count as f64 / total_f;

                entropy -= p * p.ln();
            }
        }
    }

    entropy
}

// ============================================================================
// IFS (Iterated Function System) Numerical Implementation
// ============================================================================

/// Represents an affine transformation in 2D.
///
/// The transformation is: (x', y') = (ax + by + e, cx + dy + f)
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct AffineTransform2D {
    /// Coefficients of the transformation
    pub a: f64,
    /// Coefficient b.
    pub b: f64,
    /// Coefficient c.
    pub c: f64,
    /// Coefficient d.
    pub d: f64,
    /// Coefficient e (translation in x).
    pub e: f64,
    /// Coefficient f (translation in y).
    pub f: f64,
}

impl AffineTransform2D {
    /// Creates a new affine transformation.
    #[must_use]
    pub const fn new(
        a: f64,
        b: f64,
        c: f64,
        d: f64,
        e: f64,
        f: f64,
    ) -> Self {
        Self { a, b, c, d, e, f }
    }

    /// Applies the transformation to a point.
    #[must_use]
    pub fn apply(
        &self,
        point: (f64, f64),
    ) -> (f64, f64) {
        let (x, y) = point;

        (
            self.a.mul_add(x, self.b * y) + self.e,
            self.c.mul_add(x, self.d * y) + self.f,
        )
    }
}

const IFS_DEFAULT_SEED: u64 = 12_345;
const IFS_DRAW_COUNT: f64 = 2_147_483_648.0;
const IFS_MAX_DRAW: u32 = (1 << 31) - 1;

/// Invalid input or an unrepresentable orbit for an affine IFS.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IfsError {
    /// The transform and weight counts differ or are zero.
    Shape,
    /// The weights must be finite, nonnegative and include a positive value.
    InvalidWeights,
    /// A coefficient or starting coordinate is nonfinite.
    NonFiniteInput,
    /// An affine step produced a nonfinite coordinate.
    NonFiniteOrbit,
    /// The work count overflows or a result allocation could not be reserved.
    WorkLimit,
}

impl std::fmt::Display for IfsError {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        f.write_str(match self {
            | Self::Shape => "IFS transforms and weights must have the same nonzero length",
            | Self::InvalidWeights => {
                "IFS weights must be finite and nonnegative with a positive total"
            },
            | Self::NonFiniteInput => "IFS coefficients and starting coordinates must be finite",
            | Self::NonFiniteOrbit => "IFS step produced a nonfinite coordinate",
            | Self::WorkLimit => "IFS work count or allocation exceeds capacity",
        })
    }
}

impl std::error::Error for IfsError {}

/// A prepared affine iterated function system with fixed selection weights.
///
/// Preparation keeps the reachable integer cut points of the 31-bit drawing
/// sequence. Reuse the same plan for independent starting points and seeds.
#[derive(Debug, Clone)]
pub struct AffineIfs {
    transforms: Box<[AffineTransform2D]>,
    cut_points: Box<[u32]>,
}

impl AffineIfs {
    pub(crate) fn heap_bytes(&self) -> usize {
        std::mem::size_of_val(self.transforms.as_ref())
            .saturating_add(std::mem::size_of_val(self.cut_points.as_ref()))
    }

    /// Prepares an affine IFS from finite coefficients and relative weights.
    ///
    /// Zero-weight transformations are excluded. Bounds are inclusive at
    /// the upper end of each draw interval.
    ///
    /// # Errors
    /// Invalid shape, weights, coefficients or a failed allocation.
    pub fn new(
        transforms: &[AffineTransform2D],
        weights: &[f64],
    ) -> Result<Self, IfsError> {
        if transforms.is_empty() || transforms.len() != weights.len() {
            return Err(IfsError::Shape);
        }
        if transforms.iter().any(|map| {
            ![map.a, map.b, map.c, map.d, map.e, map.f]
                .into_iter()
                .all(f64::is_finite)
        }) {
            return Err(IfsError::NonFiniteInput);
        }
        let (mut total, mut largest, mut positive) = (0.0, 0.0_f64, 0);
        for &weight in weights {
            if !weight.is_finite() || weight < 0.0 {
                return Err(IfsError::InvalidWeights);
            }
            total += weight;
            largest = largest.max(weight);
            positive += usize::from(weight > 0.0);
        }
        if positive == 0 {
            return Err(IfsError::InvalidWeights);
        }
        // Keep the old summation and normalization when the total is
        // representable. Otherwise scale the entire distribution first.
        let scale = if total.is_finite() {
            1.0
        } else {
            let power_of_two = f64::from_bits(largest.to_bits() & 0x7ff0_0000_0000_0000);
            total = weights.iter().map(|&weight| weight / power_of_two).sum();
            power_of_two
        };
        let mut selected = Vec::new();
        let mut cut_points = Vec::new();
        selected
            .try_reserve_exact(positive)
            .map_err(|_| IfsError::WorkLimit)?;
        cut_points
            .try_reserve_exact(positive)
            .map_err(|_| IfsError::WorkLimit)?;
        let mut cumulative = 0.0;
        for (&map, &weight) in transforms.iter().zip(weights) {
            if weight == 0.0 {
                continue;
            }
            cumulative += weight / scale;
            let cutoff = ((cumulative / total) * IFS_DRAW_COUNT) as u32;
            let cutoff = cutoff.min(IFS_MAX_DRAW);
            if cut_points.last().is_none_or(|&previous| previous < cutoff) {
                cut_points.push(cutoff);
                selected.push(map);
            }
        }
        if let Some(last) = cut_points.last_mut() {
            *last = IFS_MAX_DRAW;
        }
        Ok(Self {
            transforms: selected.into_boxed_slice(),
            cut_points: cut_points.into_boxed_slice(),
        })
    }

    /// Generates `num_points` transformed points after `skip` burn-in steps.
    ///
    /// The seed selects the same draw stream on every call. The starting
    /// point itself is excluded. A request for zero points skips the burn-in.
    ///
    /// # Errors
    /// A nonfinite starting point or step, work count overflow, or allocation
    /// failure.
    pub fn generate(
        &self,
        start: Point2D,
        num_points: usize,
        skip: usize,
        seed: u64,
    ) -> Result<Vec<Point2D>, IfsError> {
        if !start.0.is_finite() || !start.1.is_finite() {
            return Err(IfsError::NonFiniteInput);
        }
        let mut points = Vec::new();
        if num_points == 0 {
            return Ok(points);
        }
        skip.checked_add(num_points).ok_or(IfsError::WorkLimit)?;
        points
            .try_reserve_exact(num_points)
            .map_err(|_| IfsError::WorkLimit)?;
        let thresholds = self.cut_points.len() - 1;
        if thresholds == 0 {
            self.generate_with(start, num_points, skip, seed, points, |_| 0)
        } else if thresholds <= 3 {
            let mut small = [IFS_MAX_DRAW; 3];
            small[..thresholds].copy_from_slice(&self.cut_points[..thresholds]);
            let [a, b, c] = small;
            self.generate_with(start, num_points, skip, seed, points, |draw| {
                usize::from(a < draw) + usize::from(b < draw) + usize::from(c < draw)
            })
        } else {
            self.generate_with(start, num_points, skip, seed, points, |draw| {
                self.cut_points.partition_point(|&cutoff| cutoff < draw)
            })
        }
    }

    fn generate_with(
        &self,
        mut point: Point2D,
        num_points: usize,
        skip: usize,
        mut seed: u64,
        mut points: Vec<Point2D>,
        select: impl Fn(u32) -> usize,
    ) -> Result<Vec<Point2D>, IfsError> {
        const CHECK_INTERVAL: usize = 64;
        let mut advance = |point| {
            seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            self.transforms[select((seed >> 33) as u32)].apply(point)
        };
        let validate = |(x, y): Point2D| {
            if x.is_finite() && y.is_finite() {
                Ok(())
            } else {
                Err(IfsError::NonFiniteOrbit)
            }
        };
        // Each affine coordinate uses both inputs. A nonfinite coordinate
        // persists in subsequent steps, so a block-end check detects any
        // failure while bounding how much extra work a failure causes.
        for offset in (0..skip).step_by(CHECK_INTERVAL) {
            let block = (skip - offset).min(CHECK_INTERVAL);
            for _ in 0..block {
                point = advance(point);
            }
            validate(point)?;
        }
        for offset in (0..num_points).step_by(CHECK_INTERVAL) {
            let block = (num_points - offset).min(CHECK_INTERVAL);
            points.extend((0..block).map(|_| {
                point = advance(point);
                point
            }));
            validate(point)?;
        }
        Ok(points)
    }
}

/// Generates an IFS fractal using the chaos game algorithm.
///
/// # Arguments
/// * `transforms` - The affine transformations.
/// * `probabilities` - Probabilities for selecting each transformation.
/// * `start_point` - Initial point.
/// * `num_points` - Number of points to generate.
/// * `skip` - Number of initial points to skip (for convergence).
///
/// # Returns
/// A vector of 2D points on the fractal. Invalid input or a nonfinite orbit
/// yields an empty vector. Use [`AffineIfs`] to reuse preparation or choose a seed.
#[must_use]
pub fn generate_ifs_fractal(
    transforms: &[AffineTransform2D],
    probabilities: &[f64],
    start_point: (f64, f64),
    num_points: usize,
    skip: usize,
) -> Vec<(f64, f64)> {
    generate_ifs_checked(transforms, probabilities, start_point, num_points, skip)
        .unwrap_or_default()
}

pub(crate) fn generate_ifs_checked(
    transforms: &[AffineTransform2D],
    weights: &[f64],
    start: Point2D,
    points: usize,
    skip: usize,
) -> Result<Vec<Point2D>, IfsError> {
    AffineIfs::new(transforms, weights)?.generate(start, points, skip, IFS_DEFAULT_SEED)
}

/// Creates the Sierpinski triangle IFS.
#[must_use]
pub fn sierpinski_triangle_ifs() -> (Vec<AffineTransform2D>, Vec<f64>) {
    let transforms = vec![
        AffineTransform2D::new(0.5, 0.0, 0.0, 0.5, 0.0, 0.0),
        AffineTransform2D::new(0.5, 0.0, 0.0, 0.5, 0.5, 0.0),
        AffineTransform2D::new(0.5, 0.0, 0.0, 0.5, 0.25, 0.5),
    ];

    let probabilities = vec![1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];

    (transforms, probabilities)
}

/// Creates the Barnsley fern IFS.
#[must_use]
pub fn barnsley_fern_ifs() -> (Vec<AffineTransform2D>, Vec<f64>) {
    let transforms = vec![
        AffineTransform2D::new(0.0, 0.0, 0.0, 0.16, 0.0, 0.0),
        AffineTransform2D::new(0.85, 0.04, -0.04, 0.85, 0.0, 1.6),
        AffineTransform2D::new(0.2, -0.26, 0.23, 0.22, 0.0, 1.6),
        AffineTransform2D::new(-0.15, 0.28, 0.26, 0.24, 0.0, 0.44),
    ];

    let probabilities = vec![0.01, 0.85, 0.07, 0.07];

    (transforms, probabilities)
}
