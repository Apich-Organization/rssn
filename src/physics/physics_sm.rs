// src/physics/physics_sm.rs
// Implementation of Fourier Spectral Methods for solving PDEs.

use crate::numerical::transforms::{fft, fft_slice, ifft, ifft_slice};
use num_complex::Complex;
use rayon::prelude::*;

// --- FFT Utilities (2D & 3D) ---

/// Transposes a 2D matrix represented as a flat Vec.
pub(crate) fn transpose<T: Clone + Default>(data: &[T], width: usize, height: usize) -> Vec<T> {
    let mut transposed = vec![T::default(); width * height];
    for i in 0..height {
        for j in 0..width {
            transposed[j * height + i] = data[i * width + j].clone();
        }
    }
    transposed
}

/// Performs a 2D FFT by applying 1D FFT along rows and then columns.
pub fn fft2d(data: &mut Vec<Complex<f64>>, width: usize, height: usize) {
    // FFT on rows
    data.par_chunks_mut(width).for_each(fft_slice);

    // Transpose and FFT on columns
    let mut transposed = transpose(data, width, height);
    transposed
        .par_chunks_mut(height)
        .for_each(fft_slice);

    // Transpose back
    *data = transpose(&transposed, height, width);
}

/// Performs a 2D IFFT.
pub fn ifft2d(data: &mut Vec<Complex<f64>>, width: usize, height: usize) {
    // IFFT on rows
    data.par_chunks_mut(width).for_each(ifft_slice);

    // Transpose and IFFT on columns
    let mut transposed = transpose(data, width, height);
    transposed
        .par_chunks_mut(height)
        .for_each(ifft_slice);

    // Transpose back
    *data = transpose(&transposed, height, width);
}

// --- Wavenumber Generation ---

/// Creates a 1D wavenumber grid for FFT.
pub(crate) fn create_k_grid(n: usize, dx: f64) -> Vec<f64> {
    let dk = 2.0 * std::f64::consts::PI / (n as f64 * dx);
    let mut k: Vec<f64> = (0..n / 2).map(|i| i as f64 * dk).collect();
    let mut k_neg: Vec<f64> = (0..n / 2)
        .map(|i| -(n as f64 / 2.0 - i as f64) * dk)
        .collect();
    k.append(&mut k_neg);
    k
}

// --- 1D Solver ---

/// Solves the 1D advection-diffusion equation (u_t + c*u_x = D*u_xx) using a Fourier spectral method.
pub fn solve_advection_diffusion_1d(
    initial_condition: &[f64],
    dx: f64,
    c: f64, // Advection speed
    d: f64, // Diffusion coefficient
    dt: f64,
    steps: usize,
) -> Vec<f64> {
    let n = initial_condition.len();
    let k = create_k_grid(n, dx);

    // Convert initial condition to complex and FFT
    let mut u_hat: Vec<Complex<f64>> = initial_condition
        .iter()
        .map(|&v| Complex::new(v, 0.0))
        .collect();
    fft(&mut u_hat);

    // Time-stepping in Fourier space
    for _ in 0..steps {
        u_hat.par_iter_mut().enumerate().for_each(|(i, val)| {
            let ki = k[i];
            let advection_term = Complex::new(0.0, -c * ki);
            let diffusion_term = Complex::new(-d * ki * ki, 0.0);
            let rhs = (advection_term + diffusion_term) * *val;
            *val += rhs * dt; // Forward Euler
        });
    }

    // Inverse FFT to get solution in physical space
    ifft(&mut u_hat);
    u_hat.into_iter().map(|v| v.re).collect()
}

/// Example scenario for the 1D spectral solver.
pub fn simulate_1d_advection_diffusion_scenario() -> Vec<f64> {
    const N: usize = 128;
    const L: f64 = 2.0 * std::f64::consts::PI;
    let dx = L / N as f64;
    let x: Vec<f64> = (0..N).map(|i| i as f64 * dx).collect();
    let initial_condition: Vec<f64> = x
        .iter()
        .map(|&v| (-(v - L / 2.0).powi(2) / 0.5).exp())
        .collect();

    solve_advection_diffusion_1d(&initial_condition, dx, 1.0, 0.01, 0.01, 200)
}

// --- 2D Solver ---

pub struct AdvectionDiffusionConfig {
	pub width: usize,
    pub height: usize,
    pub dx: f64,
    pub dy: f64,
    pub c: (f64, f64), // Advection speeds (cx, cy)
    pub d: f64,        // Diffusion coefficient
    pub dt: f64,
    pub steps: usize,
}

/// Solves the 2D advection-diffusion equation using a Fourier spectral method.
pub fn solve_advection_diffusion_2d(
    initial_condition: &[f64],
	config: &AdvectionDiffusionConfig,
) -> Vec<f64> {
    let kx = create_k_grid(config.width, config.dx);
    let ky = create_k_grid(config.height, config.dy);

    let mut u_hat: Vec<Complex<f64>> = initial_condition
        .iter()
        .map(|&v| Complex::new(v, 0.0))
        .collect();
    fft2d(&mut u_hat, config.width, config.height);

    for _ in 0..config.steps {
        u_hat.par_iter_mut().enumerate().for_each(|(idx, val)| {
            let i = idx % config.width;
            let j = idx / config.width;
            let kxi = kx[i];
            let kyj = ky[j];

            let advection_term = Complex::new(0.0, -config.c.0 * kxi - config.c.1 * kyj);
            let diffusion_term = Complex::new(-config.d * (kxi * kxi + kyj * kyj), 0.0);
            let rhs = (advection_term + diffusion_term) * *val;
            *val += rhs * config.dt; // Forward Euler
        });
    }

    ifft2d(&mut u_hat, config.width, config.height);
    u_hat.into_iter().map(|v| v.re).collect()
}

/// Example scenario for the 2D spectrclsal solver.
pub fn simulate_2d_advection_diffusion_scenario() -> Vec<f64> {
    const WIDTH: usize = 64;
    const HEIGHT: usize = 64;
    const L: f64 = 2.0 * std::f64::consts::PI;
    const NU: f64 = 0.01;      // Example: Diffusion coefficient (e.g., nu)
    const U: f64 = 0.5;        // Example: Advection velocity component (x-direction)
    const V: f64 = 0.5;        // Example: Advection velocity component (y-direction)
    const DT: f64 = 0.001;
	const STEPS: usize = 100; // Example: Number of time steps
    let dx = L / WIDTH as f64;
    let dy = L / HEIGHT as f64;

	let config = AdvectionDiffusionConfig {
		width: WIDTH,
		height: HEIGHT,
		dx,
		dy,
		c: (U, V),
		d: NU,
		dt: DT,
		steps: STEPS,
	};

    let mut initial_condition = vec![0.0; WIDTH * HEIGHT];
    for j in 0..HEIGHT {
        for i in 0..WIDTH {
            let x = i as f64 * dx;
            let y = j as f64 * dy;
            let val = (-((x - L / 2.0).powi(2) + (y - L / 2.0).powi(2)) / 0.5).exp();
            initial_condition[j * WIDTH + i] = val;
        }
    }

    solve_advection_diffusion_2d(
        &initial_condition,
        &config,
    )
}

// --- 3D Solver ---
// Note: 3D FFT is memory and compute intensive.

/// Performs a 3D FFT.
pub fn fft3d(data: &mut Vec<Complex<f64>>, width: usize, height: usize, depth: usize) {
    let plane_size = width * height;
    // FFT along x-axis
    data.par_chunks_mut(width).for_each(fft_slice);

    // FFT along y-axis
    let mut transposed_xy = vec![Complex::default(); data.len()];
    for k in 0..depth {
        let plane_slice = &data[k * plane_size..(k + 1) * plane_size];
        let mut transposed_plane = transpose(plane_slice, width, height);
        transposed_plane
            .par_chunks_mut(height)
            .for_each(fft_slice);
        let retransposed_plane = transpose(&transposed_plane, height, width);
        transposed_xy[k * plane_size..(k + 1) * plane_size].copy_from_slice(&retransposed_plane);
    }
    *data = transposed_xy;

    // FFT along z-axis
    let transposed_z: Vec<Complex<f64>> = (0..plane_size)
        // 1. Parallelize over the outer loop index `i`
        .into_par_iter()
        // 2. Map each column index `i` to the processed column data
        .flat_map(|i| {
            // Build the current column `z_col` by reading from `data`
            let mut z_col: Vec<_> = (0..depth).map(|k| data[k * plane_size + i]).collect();

            // Perform the FFT calculation
            fft(&mut z_col);

            // 3. Return the processed column data, flattened into the final Vec.
            //    The order here ensures correct reassembly for the transpose.
            z_col
        })
        // 4. Rayon collects all the processed columns, stitching them together
        //    into the final `transposed_z` Vec.
        .collect();

    // Note: You no longer need `let mut transposed_z = vec![Complex::default(); data.len()];`
    // as the result of the collect is the new transposed_z.
    *data = transposed_z;
}

/// Performs a 3D IFFT.
pub fn ifft3d(data: &mut Vec<Complex<f64>>, width: usize, height: usize, depth: usize) {
    // Inverse of the 3D FFT process
    let plane_size = width * height;
    // FFT along z-axis
    let transposed_z: Vec<Complex<f64>> = (0..plane_size)
        // 1. Parallelize over the outer loop index `i`
        .into_par_iter()
        // 2. Map each column index `i` to the processed column data
        .flat_map(|i| {
            // Build the current column `z_col` by reading from `data`
            let mut z_col: Vec<_> = (0..depth).map(|k| data[k * plane_size + i]).collect();

            // Perform the FFT calculation
            ifft(&mut z_col);

            // 3. Return the processed column data, flattened into the final Vec.
            //    The order here ensures correct reassembly for the transpose.
            z_col
        })
        // 4. Rayon collects all the processed columns, stitching them together
        //    into the final `transposed_z` Vec.
        .collect();
    *data = transposed_z;

    let mut transposed_xy = vec![Complex::default(); data.len()];
    for k in 0..depth {
        let plane_slice = &data[k * plane_size..(k + 1) * plane_size];
        let mut transposed_plane = transpose(plane_slice, width, height);
        transposed_plane
            .par_chunks_mut(height)
            .for_each(ifft_slice);
        let retransposed_plane = transpose(&transposed_plane, height, width);
        transposed_xy[k * plane_size..(k + 1) * plane_size].copy_from_slice(&retransposed_plane);
    }
    *data = transposed_xy;

    data.par_chunks_mut(width).for_each(ifft_slice);
}

pub struct AdvectionDiffusionConfig3d{
	pub width: usize,
    pub height: usize,
    pub depth: usize,
    pub dx: f64,
    pub dy: f64,
    pub dz: f64,
    pub c: (f64, f64, f64), // Advection speeds (cx, cy, cz)
    pub d: f64,             // Diffusion coefficient
    pub dt: f64,
    pub steps: usize,
}

/// Solves the 3D advection-diffusion equation.
pub fn solve_advection_diffusion_3d(
    initial_condition: &[f64],
	config: &AdvectionDiffusionConfig3d,
) -> Vec<f64> {
    let kx = create_k_grid(config.width, config.dx);
    let ky = create_k_grid(config.height, config.dy);
    let kz = create_k_grid(config.depth, config.dz);

    let mut u_hat: Vec<Complex<f64>> = initial_condition
        .iter()
        .map(|&v| Complex::new(v, 0.0))
        .collect();
    fft3d(&mut u_hat, config.width, config.height, config.depth);

    let plane_size = config.width * config.height;
    for _ in 0..config.steps {
        u_hat.par_iter_mut().enumerate().for_each(|(idx, val)| {
            let i = idx % config.width;
            let j = (idx / config.width) % config.height;
            let k = idx / plane_size;
            let kxi = kx[i];
            let kyj = ky[j];
            let kzk = kz[k];

            let advection = Complex::new(0.0, -config.c.0 * kxi - config.c.1 * kyj - config.c.2 * kzk);
            let diffusion = Complex::new(-config.d * (kxi.powi(2) + kyj.powi(2) + kzk.powi(2)), 0.0);
            let rhs = (advection + diffusion) * *val;
            *val += rhs * config.dt; // Forward Euler
        });
    }

    ifft3d(&mut u_hat, config.width, config.height, config.depth);
    u_hat.into_iter().map(|v| v.re).collect()
}

/// Example scenario for the 3D spectral solver.
pub fn simulate_3d_advection_diffusion_scenario() -> Vec<f64> {
    const N: usize = 16; // Keep N small for 3D
    const L: f64 = 2.0 * std::f64::consts::PI;
    let d = L / N as f64;
	
	let config = AdvectionDiffusionConfig3d {
		width: N,
		height: N,
		depth: N,
		dx: d,
		dy: d,
		dz: d,
		c: (1.0, 0.5, 0.2), // Advection speeds (cx, cy, cz)
		d: 0.01,             // Diffusion coefficient
		dt: 0.005,
		steps: 200,
	};
	
    let mut initial_condition = vec![0.0; N * N * N];
    for k in 0..N {
        for j in 0..N {
            for i in 0..N {
                let x = i as f64 * d;
                let y = j as f64 * d;
                let z = k as f64 * d;
                let val =
                    (-((x - L / 2.0).powi(2) + (y - L / 2.0).powi(2) + (z - L / 2.0).powi(2))
                        / 0.5)
                        .exp();
                initial_condition[(k * N + j) * N + i] = val;
            }
        }
    }

    solve_advection_diffusion_3d(
        &initial_condition,
        &config
    )
}
