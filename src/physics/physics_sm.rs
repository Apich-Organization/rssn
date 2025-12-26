use crate::numerical::transforms::{
    fft,
    fft_slice,
    ifft,
    ifft_slice,
};
use num_complex::Complex;
use rayon::prelude::*;
use serde::{
    Deserialize,
    Serialize,
};

/// Transposes a 2D matrix represented as a flat Vec.

pub(crate) fn transpose<T: Clone + Default + Send + Sync>(
    data: &[T],
    width: usize,
    height: usize,
) -> Vec<T> {

    let mut transposed = vec![T::default(); width * height];

    transposed
        .par_iter_mut()
        .enumerate()
        .for_each(|(idx, val)| {

            let i = idx / height;

            let j = idx % height;

            // transposed[i * height + j] comes from data[j * width + i]
            *val = data[j * width + i].clone();
        });

    transposed
}

/// Performs a 2D FFT by applying 1D FFT along rows and then columns.

pub fn fft2d(
    data: &mut Vec<Complex<f64>>,
    width: usize,
    height: usize,
) {

    data.par_chunks_mut(width)
        .for_each(fft_slice);

    let mut transposed = transpose(data, width, height);

    transposed
        .par_chunks_mut(height)
        .for_each(fft_slice);

    *data = transpose(
        &transposed,
        height,
        width,
    );
}

/// Performs a 2D IFFT.

pub fn ifft2d(
    data: &mut Vec<Complex<f64>>,
    width: usize,
    height: usize,
) {

    data.par_chunks_mut(width)
        .for_each(ifft_slice);

    let mut transposed = transpose(data, width, height);

    transposed
        .par_chunks_mut(height)
        .for_each(ifft_slice);

    *data = transpose(
        &transposed,
        height,
        width,
    );
}

/// Creates a 1D wavenumber grid for FFT.

pub(crate) fn create_k_grid(
    n: usize,
    dx: f64,
) -> Vec<f64> {

    let dk = 2.0 * std::f64::consts::PI / (n as f64 * dx);

    let mut k: Vec<f64> = (0..n / 2)
        .map(|i| i as f64 * dk)
        .collect();

    let mut k_neg: Vec<f64> = (0..n / 2)
        .map(|i| -(n as f64 / 2.0 - i as f64) * dk)
        .collect();

    k.append(&mut k_neg);

    k
}

/// Solves the 1D advection-diffusion equation (u_t + c*u_x = D*u_xx) using a Fourier spectral method.

pub fn solve_advection_diffusion_1d(
    initial_condition: &[f64],
    dx: f64,
    c: f64,
    d: f64,
    dt: f64,
    steps: usize,
) -> Vec<f64> {

    let n = initial_condition.len();

    let k = create_k_grid(n, dx);

    let mut u_hat: Vec<Complex<f64>> = initial_condition
        .iter()
        .map(|&v| Complex::new(v, 0.0))
        .collect();

    fft(&mut u_hat);

    for _ in 0..steps {

        u_hat
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, val)| {

                let ki = k[i];

                let advection_term = Complex::new(0.0, -c * ki);

                let diffusion_term = Complex::new(-d * ki * ki, 0.0);

                let rhs = (advection_term + diffusion_term) * *val;

                *val += rhs * dt;
            });
    }

    ifft(&mut u_hat);

    u_hat
        .into_iter()
        .map(|v| v.re)
        .collect()
}

/// Example scenario for the 1D spectral solver.

pub fn simulate_1d_advection_diffusion_scenario() -> Vec<f64> {

    const N: usize = 128;

    const L: f64 = 2.0 * std::f64::consts::PI;

    let dx = L / N as f64;

    let x: Vec<f64> = (0..N)
        .map(|i| i as f64 * dx)
        .collect();

    let initial_condition: Vec<f64> = x
        .iter()
        .map(|&v| (-(v - L / 2.0).powi(2) / 0.5).exp())
        .collect();

    solve_advection_diffusion_1d(
        &initial_condition,
        dx,
        1.0,
        0.01,
        0.01,
        200,
    )
}

#[derive(Debug, Clone, Serialize, Deserialize)]

pub struct AdvectionDiffusionConfig {
    pub width: usize,
    pub height: usize,
    pub dx: f64,
    pub dy: f64,
    pub c: (f64, f64),
    pub d: f64,
    pub dt: f64,
    pub steps: usize,
}

/// Solves the 2D advection-diffusion equation using a Fourier spectral method.

pub fn solve_advection_diffusion_2d(
    initial_condition: &[f64],
    config: &AdvectionDiffusionConfig,
) -> Vec<f64> {

    let kx = create_k_grid(
        config.width,
        config.dx,
    );

    let ky = create_k_grid(
        config.height,
        config.dy,
    );

    let mut u_hat: Vec<Complex<f64>> = initial_condition
        .iter()
        .map(|&v| Complex::new(v, 0.0))
        .collect();

    fft2d(
        &mut u_hat,
        config.width,
        config.height,
    );

    for _ in 0..config.steps {

        u_hat
            .par_iter_mut()
            .enumerate()
            .for_each(|(idx, val)| {

                let i = idx % config.width;

                let j = idx / config.width;

                let kxi = kx[i];

                let kyj = ky[j];

                let advection_term = Complex::new(
                    0.0,
                    -config.c.0 * kxi - config.c.1 * kyj,
                );

                let diffusion_term = Complex::new(
                    -config.d * (kxi * kxi + kyj * kyj),
                    0.0,
                );

                let rhs = (advection_term + diffusion_term) * *val;

                *val += rhs * config.dt;
            });
    }

    ifft2d(
        &mut u_hat,
        config.width,
        config.height,
    );

    u_hat
        .into_iter()
        .map(|v| v.re)
        .collect()
}

/// Example scenario for the 2D spectrclsal solver.

pub fn simulate_2d_advection_diffusion_scenario() -> Vec<f64> {

    const WIDTH: usize = 64;

    const HEIGHT: usize = 64;

    const L: f64 = 2.0 * std::f64::consts::PI;

    const NU: f64 = 0.01;

    const U: f64 = 0.5;

    const V: f64 = 0.5;

    const DT: f64 = 0.001;

    const STEPS: usize = 100;

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

/// Performs a 3D FFT.

pub fn fft3d(
    data: &mut Vec<Complex<f64>>,
    width: usize,
    height: usize,
    depth: usize,
) {

    let plane_size = width * height;

    data.par_chunks_mut(width)
        .for_each(fft_slice);

    let mut transposed_xy = vec![Complex::default(); data.len()];

    for k in 0..depth {

        let plane_slice = &data[k * plane_size..(k + 1) * plane_size];

        let mut transposed_plane = transpose(
            plane_slice,
            width,
            height,
        );

        transposed_plane
            .par_chunks_mut(height)
            .for_each(fft_slice);

        let retransposed_plane = transpose(
            &transposed_plane,
            height,
            width,
        );

        transposed_xy[k * plane_size..(k + 1) * plane_size].copy_from_slice(&retransposed_plane);
    }

    *data = transposed_xy;

    let transposed_z: Vec<Complex<f64>> = (0..plane_size)
        .into_par_iter()
        .flat_map(|i| {

            let mut z_col: Vec<_> = (0..depth)
                .map(|k| data[k * plane_size + i])
                .collect();

            fft(&mut z_col);

            z_col
        })
        .collect();

    *data = transposed_z;
}

/// Performs a 3D IFFT.

pub fn ifft3d(
    data: &mut Vec<Complex<f64>>,
    width: usize,
    height: usize,
    depth: usize,
) {

    let plane_size = width * height;

    let transposed_z: Vec<Complex<f64>> = (0..plane_size)
        .into_par_iter()
        .flat_map(|i| {

            let mut z_col: Vec<_> = (0..depth)
                .map(|k| data[k * plane_size + i])
                .collect();

            ifft(&mut z_col);

            z_col
        })
        .collect();

    *data = transposed_z;

    let mut transposed_xy = vec![Complex::default(); data.len()];

    for k in 0..depth {

        let plane_slice = &data[k * plane_size..(k + 1) * plane_size];

        let mut transposed_plane = transpose(
            plane_slice,
            width,
            height,
        );

        transposed_plane
            .par_chunks_mut(height)
            .for_each(ifft_slice);

        let retransposed_plane = transpose(
            &transposed_plane,
            height,
            width,
        );

        transposed_xy[k * plane_size..(k + 1) * plane_size].copy_from_slice(&retransposed_plane);
    }

    *data = transposed_xy;

    data.par_chunks_mut(width)
        .for_each(ifft_slice);
}

#[derive(Debug, Clone, Serialize, Deserialize)]

pub struct AdvectionDiffusionConfig3d {
    pub width: usize,
    pub height: usize,
    pub depth: usize,
    pub dx: f64,
    pub dy: f64,
    pub dz: f64,
    pub c: (f64, f64, f64),
    pub d: f64,
    pub dt: f64,
    pub steps: usize,
}

/// Solves the 3D advection-diffusion equation.

pub fn solve_advection_diffusion_3d(
    initial_condition: &[f64],
    config: &AdvectionDiffusionConfig3d,
) -> Vec<f64> {

    let kx = create_k_grid(
        config.width,
        config.dx,
    );

    let ky = create_k_grid(
        config.height,
        config.dy,
    );

    let kz = create_k_grid(
        config.depth,
        config.dz,
    );

    let mut u_hat: Vec<Complex<f64>> = initial_condition
        .iter()
        .map(|&v| Complex::new(v, 0.0))
        .collect();

    fft3d(
        &mut u_hat,
        config.width,
        config.height,
        config.depth,
    );

    let plane_size = config.width * config.height;

    for _ in 0..config.steps {

        u_hat
            .par_iter_mut()
            .enumerate()
            .for_each(|(idx, val)| {

                let i = idx % config.width;

                let j = (idx / config.width) % config.height;

                let k = idx / plane_size;

                let kxi = kx[i];

                let kyj = ky[j];

                let kzk = kz[k];

                let advection = Complex::new(
                    0.0,
                    -config.c.0 * kxi - config.c.1 * kyj - config.c.2 * kzk,
                );

                let diffusion = Complex::new(
                    -config.d * (kxi.powi(2) + kyj.powi(2) + kzk.powi(2)),
                    0.0,
                );

                let rhs = (advection + diffusion) * *val;

                *val += rhs * config.dt;
            });
    }

    ifft3d(
        &mut u_hat,
        config.width,
        config.height,
        config.depth,
    );

    u_hat
        .into_iter()
        .map(|v| v.re)
        .collect()
}

/// Example scenario for the 3D spectral solver.

pub fn simulate_3d_advection_diffusion_scenario() -> Vec<f64> {

    const N: usize = 16;

    const L: f64 = 2.0 * std::f64::consts::PI;

    let d = L / N as f64;

    let config = AdvectionDiffusionConfig3d {
        width: N,
        height: N,
        depth: N,
        dx: d,
        dy: d,
        dz: d,
        c: (1.0, 0.5, 0.2),
        d: 0.01,
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
        &config,
    )
}
