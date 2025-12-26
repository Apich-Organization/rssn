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
//! 
//! use rssn::numerical::fractal_geometry_and_chaos::generate_mandelbrot_set;
//!
//! let data = generate_mandelbrot_set(
//!     100,
//!     100,
//!     (-2.0, 1.0),
//!     (-1.5, 1.5),
//!     100,
//! );
//!
//! assert_eq!(data.len(), 100);
//!
//! assert_eq!(data[0].len(), 100);
//! ```
//!
//! ### Simulating the Lorenz Attractor
//! ```
//! 
//! use rssn::numerical::fractal_geometry_and_chaos::generate_lorenz_attractor;
//!
//! let trajectory = generate_lorenz_attractor(
//!     (1.0, 1.0, 1.0),
//!     0.01,
//!     1000,
//! );
//!
//! assert_eq!(
//!     trajectory.len(),
//!     1000
//! );
//! ```

use num_complex::Complex;
use serde::{Deserialize, Serialize};

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
/// 
/// use rssn::numerical::fractal_geometry_and_chaos::generate_mandelbrot_set;
///
/// let data = generate_mandelbrot_set(
///     100,
///     100,
///     (-2.5, 1.0),
///     (-1.5, 1.5),
///     50,
/// );
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

    for (r, row) in data
        .iter_mut()
        .enumerate()
    {

        for (c, val) in row
            .iter_mut()
            .enumerate()
        {

            let x0 = x_range.0 + c as f64 * x_scale;

            let y0 = y_range.0 + r as f64 * y_scale;

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
/// 
/// use rssn::numerical::fractal_geometry_and_chaos::generate_julia_set;
///
/// // Julia set for c = -0.4 + 0.6i
/// let data = generate_julia_set(
///     100,
///     100,
///     (-2.0, 2.0),
///     (-2.0, 2.0),
///     (-0.4, 0.6),
///     50,
/// );
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

    for (r, row) in data
        .iter_mut()
        .enumerate()
    {

        for (col, val) in row
            .iter_mut()
            .enumerate()
        {

            let x0 = x_range.0 + col as f64 * x_scale;

            let y0 = y_range.0 + r as f64 * y_scale;

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

    for (r, row) in data
        .iter_mut()
        .enumerate()
    {

        for (col, val) in row
            .iter_mut()
            .enumerate()
        {

            let cx = x_range.0 + col as f64 * x_scale;

            let cy = y_range.0 + r as f64 * y_scale;

            let mut zx = 0.0;

            let mut zy = 0.0;

            let mut iter = 0;

            while zx * zx + zy * zy <= 4.0 && iter < max_iter {

                let xtemp = zx * zx - zy * zy + cx;

                zy = (2.0 * zx.abs() * zy.abs()) + cy;

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

    for (r, row) in data
        .iter_mut()
        .enumerate()
    {

        for (col, val) in row
            .iter_mut()
            .enumerate()
        {

            let x0 = x_range.0 + col as f64 * x_scale;

            let y0 = y_range.0 + r as f64 * y_scale;

            let mut z = Complex::new(0.0, 0.0);

            let c_val = Complex::new(x0, y0);

            let mut iter = 0;

            // Escape radius depends on d
            let escape_radius = 2.0_f64.max(
                c_val
                    .norm()
                    .powf(1.0 / (d - 1.0)),
            );

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
        Complex::new(
            -0.5,
            3.0_f64.sqrt() / 2.0,
        ),
        Complex::new(
            -0.5,
            -3.0_f64.sqrt() / 2.0,
        ),
    ];

    let x_scale = (x_range.1 - x_range.0) / width as f64;

    let y_scale = (y_range.1 - y_range.0) / height as f64;

    for (r, row) in data
        .iter_mut()
        .enumerate()
    {

        for (col, val) in row
            .iter_mut()
            .enumerate()
        {

            let x0 = x_range.0 + col as f64 * x_scale;

            let y0 = y_range.0 + r as f64 * y_scale;

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

                z = z - f / df;

                // Check if close to any root
                for (i, &root) in roots
                    .iter()
                    .enumerate()
                {

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
/// 
/// use rssn::numerical::fractal_geometry_and_chaos::generate_lorenz_attractor;
///
/// let points = generate_lorenz_attractor(
///     (0.1, 0.0, 0.0),
///     0.01,
///     1000,
/// );
///
/// assert_eq!(points.len(), 1000);
/// ```
#[must_use]

pub fn generate_lorenz_attractor(
    start_point: (f64, f64, f64),
    dt: f64,
    num_steps: usize,
) -> Vec<(f64, f64, f64)> {

    generate_lorenz_attractor_custom(
        start_point,
        dt,
        num_steps,
        10.0,
        28.0,
        8.0 / 3.0,
    )
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

        let dy = x * (rho - z) - y;

        let dz = x * y - beta * z;

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

        let dy = x + a * y;

        let dz = b + z * (x - c);

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
/// - x_{n+1} = 1 - a*x_n^2 + y_n
/// - y_{n+1} = b*x_n
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

        let x_new = 1.0 - a * x * x + y;

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
/// - x_{n+1} = x_n^2 - y_n^2 + ax_n + by_n
/// - y_{n+1} = 2x_n*y_n + cx_n + dy_n
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

        let x_new = x * x - y * y + a * x + b * y;

        let y_new = 2.0 * x * y + c * x + d * y;

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
/// The logistic map is defined by: x_{n+1} = r * x_n * (1 - x_n)
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

        let r = r_range.0 + i as f64 * r_step;

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
        let deriv = r * (1.0 - 2.0 * x);

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

        let dy_main = x * (rho - z) - y;

        let dz_main = x * y - beta * z;

        x += dx_main * dt;

        y += dy_main * dt;

        z += dz_main * dt;

        // Evolve the perturbation using linearized equations
        let ddx = sigma * (dy - dx);

        let ddy = dx * (rho - z) - x * dz - dy;

        let ddz = dx * y + x * dy - beta * dz;

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
        let eps = extent / num_boxes as f64;

        // Count occupied boxes using a hash set
        use std::collections::HashSet;

        let mut occupied = HashSet::new();

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

    linear_regression_slope(&log_inv_eps, &log_n)
}

/// Estimates the correlation dimension from a set of points.
///
/// The correlation dimension is estimated from the correlation integral:
/// C(r) ~ r^D
///
/// # Arguments
/// * `points` - A slice of 2D points.
/// * `num_radii` - Number of radius values to sample.
///
/// # Returns
/// An estimate of the correlation dimension.
#[must_use]

pub fn correlation_dimension(
    points: &[(f64, f64)],
    num_radii: usize,
) -> f64 {

    if points.len() < 2 || num_radii < 2 {

        return 0.0;
    }

    // Compute all pairwise distances
    let n = points.len();

    let mut distances: Vec<f64> = Vec::with_capacity(n * (n - 1) / 2);

    for i in 0..n {

        for j in (i + 1)..n {

            let dx = points[i].0 - points[j].0;

            let dy = points[i].1 - points[j].1;

            distances.push((dx * dx + dy * dy).sqrt());
        }
    }

    if distances.is_empty() {

        return 0.0;
    }

    distances.sort_by(|a, b| {

        a.partial_cmp(b)
            .unwrap()
    });

    let r_min = distances[0].max(1e-10);

    let r_max = distances[distances.len() - 1];

    if r_max <= r_min {

        return 0.0;
    }

    // Sample radii logarithmically
    let mut log_r = Vec::with_capacity(num_radii);

    let mut log_c = Vec::with_capacity(num_radii);

    for i in 0..num_radii {

        let log_r_val =
            r_min.ln() + (r_max.ln() - r_min.ln()) * (i as f64) / (num_radii as f64 - 1.0);

        let r = log_r_val.exp();

        // Count pairs with distance < r
        let count = distances.partition_point(|&d| d < r);

        let c_r = (2.0 * count as f64) / ((n * (n - 1)) as f64);

        if c_r > 0.0 {

            log_r.push(r.ln());

            log_c.push(c_r.ln());
        }
    }

    if log_r.len() < 2 {

        return 0.0;
    }

    linear_regression_slope(&log_r, &log_c)
}

/// Helper function for linear regression slope estimation.

fn linear_regression_slope(
    x: &[f64],
    y: &[f64],
) -> f64 {

    let n = x.len() as f64;

    let sum_x: f64 = x.iter().sum();

    let sum_y: f64 = y.iter().sum();

    let sum_xy: f64 = x
        .iter()
        .zip(y.iter())
        .map(|(xi, yi)| xi * yi)
        .sum();

    let sum_xx: f64 = x
        .iter()
        .map(|xi| xi * xi)
        .sum();

    let denom = n * sum_xx - sum_x * sum_x;

    if denom.abs() < 1e-10 {

        return 0.0;
    }

    (n * sum_xy - sum_x * sum_y) / denom
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

            let xi = ((x - x_range.0) / x_scale) as usize;

            let yi = ((y - y_range.0) / y_scale) as usize;

            if xi < x_bins && yi < y_bins {

                density[yi][xi] += 1;
            }
        }
    }

    density
}

/// Computes the entropy of an orbit from its density distribution.
///
/// Uses Shannon entropy: H = -sum(p_i * log(p_i))
///
/// # Arguments
/// * `density` - A 2D density histogram from `orbit_density`.
///
/// # Returns
/// The Shannon entropy of the distribution.
#[must_use]

pub fn orbit_entropy(density: &[Vec<usize>]) -> f64 {

    let total: usize = density
        .iter()
        .flat_map(|row| row.iter())
        .sum();

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
    pub b: f64,
    pub c: f64,
    pub d: f64,
    pub e: f64,
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
            self.a * x + self.b * y + self.e,
            self.c * x + self.d * y + self.f,
        )
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
/// A vector of 2D points on the fractal.
#[must_use]

pub fn generate_ifs_fractal(
    transforms: &[AffineTransform2D],
    probabilities: &[f64],
    start_point: (f64, f64),
    num_points: usize,
    skip: usize,
) -> Vec<(f64, f64)> {

    if transforms.is_empty() || transforms.len() != probabilities.len() {

        return vec![];
    }

    // Build cumulative probabilities
    let mut cumulative = Vec::with_capacity(probabilities.len());

    let mut sum = 0.0;

    for &p in probabilities {

        sum += p;

        cumulative.push(sum);
    }

    // Normalize
    for c in &mut cumulative {

        *c /= sum;
    }

    let mut points = Vec::with_capacity(num_points);

    let mut current = start_point;

    // Simple LCG random number generator for reproducibility
    let mut rng_state = 12345u64;

    let lcg_next = |state: &mut u64| -> f64 {

        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1);

        (*state >> 33) as f64 / (u32::MAX as f64)
    };

    for i in 0..(num_points + skip) {

        let r = lcg_next(&mut rng_state);

        let idx = cumulative
            .iter()
            .position(|&c| r <= c)
            .unwrap_or(0);

        current = transforms[idx].apply(current);

        if i >= skip {

            points.push(current);
        }
    }

    points
}

/// Creates the Sierpinski triangle IFS.
#[must_use]

pub fn sierpinski_triangle_ifs() -> (
    Vec<AffineTransform2D>,
    Vec<f64>,
) {

    let transforms = vec![
        AffineTransform2D::new(
            0.5, 0.0, 0.0, 0.5, 0.0, 0.0,
        ),
        AffineTransform2D::new(
            0.5, 0.0, 0.0, 0.5, 0.5, 0.0,
        ),
        AffineTransform2D::new(
            0.5, 0.0, 0.0, 0.5, 0.25, 0.5,
        ),
    ];

    let probabilities = vec![
        1.0 / 3.0,
        1.0 / 3.0,
        1.0 / 3.0,
    ];

    (
        transforms,
        probabilities,
    )
}

/// Creates the Barnsley fern IFS.
#[must_use]

pub fn barnsley_fern_ifs() -> (
    Vec<AffineTransform2D>,
    Vec<f64>,
) {

    let transforms = vec![
        AffineTransform2D::new(
            0.0, 0.0, 0.0, 0.16, 0.0, 0.0,
        ),
        AffineTransform2D::new(
            0.85, 0.04, -0.04, 0.85, 0.0, 1.6,
        ),
        AffineTransform2D::new(
            0.2, -0.26, 0.23, 0.22, 0.0, 1.6,
        ),
        AffineTransform2D::new(
            -0.15, 0.28, 0.26, 0.24, 0.0, 0.44,
        ),
    ];

    let probabilities = vec![
        0.01, 0.85, 0.07, 0.07,
    ];

    (
        transforms,
        probabilities,
    )
}
