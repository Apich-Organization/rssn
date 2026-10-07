//! Pseudo-spectral (Fourier, Chebyshev) tests.

use num_complex::Complex;
use rssn::sim::boundary::{BcValue, BoundaryCondition as Bc, Face, FaceBcs};
use rssn::sim::spectral::{ChebGrid, SpectralGrid, etdrk4, fft_nd};
use std::f64::consts::PI;

#[test]
fn fft_nd_roundtrip_and_known_transform() {
    let n = [8, 4, 16];
    let total: usize = n.iter().product();
    let data: Vec<Complex<f64>> = (0..total).map(|i| Complex::new((i as f64 * 0.37).sin(), (i as f64 * 0.11).cos())).collect();
    let mut a = data.clone();
    fft_nd(&mut a, &n, false).unwrap();
    fft_nd(&mut a, &n, true).unwrap();
    for (x, y) in a.iter().zip(&data) {
        assert!((x - y).norm() < 1e-12);
    }
    // single plane wave goes to a single mode
    let mut w: Vec<Complex<f64>> = (0..total)
        .map(|p| {
            let (i, j, k) = (p % 8, (p / 8) % 4, p / 32);
            Complex::new(0.0, 2.0 * PI * (2.0 * i as f64 / 8.0 + 1.0 * j as f64 / 4.0 + 3.0 * k as f64 / 16.0)).exp()
        })
        .collect();
    fft_nd(&mut w, &n, false).unwrap();
    let peak = 2 + 8 * (1 + 4 * 3);
    for (p, v) in w.iter().enumerate() {
        if p == peak {
            assert!((v.re - total as f64).abs() < 1e-9);
        } else {
            assert!(v.norm() < 1e-9);
        }
    }
    assert!(fft_nd(&mut [Complex::new(0.0, 0.0); 6], &[6], false).is_err());
}

#[test]
fn fourier_derivative_and_helmholtz_converge_spectrally() {
    for dim in 1..=3 {
        let g = SpectralGrid::new(&vec![16; dim], &vec![2.0 * PI; dim]).unwrap();
        let u = |x: [f64; 3]| (0..dim).map(|a| (x[a] + 0.3 * a as f64).sin()).sum::<f64>();
        let uh = g.to_spec(&(0..g.size()).map(|p| u(g.coord(p))).collect::<Vec<_>>());
        let dx = g.to_phys(&g.deriv(&uh, 0));
        for (p, v) in dx.iter().enumerate() {
            assert!((v - (g.coord(p)[0]).cos()).abs() < 1e-12);
        }
    }
    let err = |n: usize| {
        let g = SpectralGrid::new(&[n, n], &[2.0 * PI; 2]).unwrap();
        let sigma = 1.5;
        let u = |x: [f64; 3]| (x[0].cos() + x[1].cos()).exp();
        let f = |x: [f64; 3]| u(x) * (-x[0].sin().powi(2) + x[0].cos() - x[1].sin().powi(2) + x[1].cos() + sigma);
        let rhs: Vec<f64> = (0..g.size()).map(|p| f(g.coord(p))).collect();
        let sol = g.helmholtz(sigma, &rhs).unwrap();
        sol.iter().enumerate().map(|(p, v)| (v - u(g.coord(p))).abs()).fold(0.0, f64::max)
    };
    let (e16, e32) = (err(16), err(32));
    assert!(e16 < 1e-6 && e32 < 1e-13, "{e16} {e32}");
    let g = SpectralGrid::new(&[8], &[1.0]).unwrap();
    assert!(g.helmholtz(0.0, &vec![1.0; 8]).is_err());
    assert!(g.helmholtz(0.0, &(0..8).map(|i| (2.0 * PI * i as f64 / 8.0).sin()).collect::<Vec<_>>()).is_ok());
}

/// ETDRK4 on `u_t = ν Δu - u (∇·1) u + f` with `u = 0.5 sin(Σx - t)` in `dim` dimensions.
fn etd_error(dim: usize, n: usize, steps: usize) -> f64 {
    let nu = 0.1;
    let g = SpectralGrid::new(&vec![n; dim], &vec![2.0 * PI; dim]).unwrap();
    let d = dim as f64;
    let phase = move |x: [f64; 3], t: f64| (0..dim).map(|a| x[a]).sum::<f64>() - t;
    let u_ex = move |x: [f64; 3], t: f64| 0.5 * phase(x, t).sin();
    let forcing = move |x: [f64; 3], t: f64| {
        let ph = phase(x, t);
        -0.5 * ph.cos() + 0.125 * d * (2.0 * ph).sin() + nu * d * 0.5 * ph.sin()
    };
    let lin: Vec<Complex<f64>> = (0..g.size()).map(|p| Complex::new(-nu * g.k2(p), 0.0)).collect();
    let mask = g.dealias_mask();
    let gg = g.clone();
    let nonlinear = move |t: f64, uh: &[Complex<f64>]| -> Vec<Complex<f64>> {
        let u = gg.to_phys(uh);
        let mut grad_sum = vec![Complex::new(0.0, 0.0); uh.len()];
        for a in 0..dim {
            let da = gg.deriv(uh, a);
            for (s, v) in grad_sum.iter_mut().zip(&da) {
                *s += v;
            }
        }
        let ux = gg.to_phys(&grad_sum);
        let f: Vec<f64> = (0..u.len()).map(|p| forcing(gg.coord(p), t)).collect();
        let nl: Vec<f64> = (0..u.len()).map(|p| -u[p] * ux[p] + f[p]).collect();
        let mut nh = gg.to_spec(&nl);
        for (v, &keep) in nh.iter_mut().zip(&mask) {
            if !keep {
                *v = Complex::new(0.0, 0.0);
            }
        }
        nh
    };
    let u0: Vec<f64> = (0..g.size()).map(|p| u_ex(g.coord(p), 0.0)).collect();
    let uh = etdrk4(g.to_spec(&u0), &lin, &nonlinear, 1.0, steps);
    let u = g.to_phys(&uh);
    u.iter().enumerate().map(|(p, v)| (v - u_ex(g.coord(p), 1.0)).abs()).fold(0.0, f64::max)
}

#[test]
fn etdrk4_is_fourth_order_in_1d_2d_3d() {
    for (dim, n) in [(1, 16), (2, 8), (3, 8)] {
        let e: Vec<f64> = [10, 20, 40].iter().map(|&s| etd_error(dim, n, s)).collect();
        let o = (e[0] / e[1]).log2().min((e[1] / e[2]).log2());
        assert!(o > 3.5 || e[2] < 1e-11, "dim {dim}: order {o}, {e:?}");
    }
}

fn u2(x: [f64; 3]) -> f64 {
    (x[0] + x[1]).exp()
}

#[test]
fn chebyshev_helmholtz_exponential_convergence_1d_2d_3d() {
    // 1D: -u'' + 2u = f, u = exp(x) on [0,1]; Neumann at 0, Robin at 1
    let u1 = |x: [f64; 3]| x[0].exp();
    let bcs1 = FaceBcs::dirichlet_zero()
        .with(Face::XMin, Bc::Neumann(BcValue::space(move |x| -u1(x))))
        .with(Face::XMax, Bc::Robin { alpha: 2.0, beta: 1.0, g: BcValue::space(move |x| 3.0 * u1(x)) });
    let run1 = |n: usize| {
        let g = ChebGrid::new(&[n], &[0.0], &[1.0]).unwrap();
        let u = g.solve_helmholtz(&bcs1, 2.0, |x| -u1(x) + 2.0 * u1(x)).unwrap();
        g.max_error(&u, u1)
    };
    assert!(run1(6) < 1e-4 && run1(12) < 1e-11, "{} {}", run1(6), run1(12));
    // 2D mixed
    let bcs2 = FaceBcs::dirichlet_zero()
        .with(Face::XMin, Bc::Dirichlet(BcValue::space(u2)))
        .with(Face::XMax, Bc::Dirichlet(BcValue::space(u2)))
        .with(Face::YMin, Bc::Neumann(BcValue::space(move |x| -u2(x))))
        .with(Face::YMax, Bc::Robin { alpha: 2.0, beta: 1.0, g: BcValue::space(move |x| 3.0 * u2(x)) });
    let run2 = |n: usize| {
        let g = ChebGrid::new(&[n, n], &[0.0; 2], &[1.0; 2]).unwrap();
        let u = g.solve_helmholtz(&bcs2, 0.0, |x| -2.0 * u2(x)).unwrap();
        g.max_error(&u, u2)
    };
    let (a, b) = (run2(6), run2(12));
    assert!(a < 1e-3 && b < 1e-9, "{a} {b}");
    // 3D
    let u3 = |x: [f64; 3]| (x[0] + x[1] + x[2]).exp();
    let bcs3 = FaceBcs::uniform(&Bc::Dirichlet(BcValue::space(u3)))
        .with(Face::ZMax, Bc::Neumann(BcValue::space(u3)));
    let g = ChebGrid::new(&[8, 8, 8], &[0.0; 3], &[1.0; 3]).unwrap();
    let u = g.solve_helmholtz(&bcs3, 0.0, |x| -3.0 * u3(x)).unwrap();
    assert!(g.max_error(&u, u3) < 1e-6);
    // unsupported
    let bad = FaceBcs::uniform(&Bc::Periodic);
    assert!(g.solve_helmholtz(&bad, 0.0, |_| 0.0).is_err());
}

#[test]
fn chebyshev_heat_equation_second_order_in_time_with_spectral_space() {
    let bcs = FaceBcs::dirichlet_zero();
    let exact = |x: [f64; 3], t: f64| (-PI * PI * t).exp() * (PI * x[0]).sin();
    let run = |steps: usize| {
        let g = ChebGrid::new(&[16], &[0.0], &[1.0]).unwrap();
        let u = g.heat(&bcs, 1.0, |x| exact(x, 0.0), None, 0.2, steps, 0.5).unwrap();
        g.max_error(&u, |x| exact(x, 0.2))
    };
    let (e1, e2) = (run(20), run(40));
    assert!((e1 / e2).log2() > 1.9 && e2 < 1e-4, "{e1} {e2}");
    // 2D Neumann cosine mode
    let bcs = FaceBcs::uniform(&Bc::neumann(0.0));
    let g = ChebGrid::new(&[14, 14], &[0.0; 2], &[1.0; 2]).unwrap();
    let ex = |x: [f64; 3], t: f64| (-2.0 * PI * PI * t).exp() * (PI * x[0]).cos() * (PI * x[1]).cos();
    let u = g.heat(&bcs, 1.0, |x| ex(x, 0.0), None, 0.1, 100, 0.5).unwrap();
    assert!(g.max_error(&u, |x| ex(x, 0.1)) < 1e-4);
    // time dependent Dirichlet data with a source
    let bcs = FaceBcs::uniform(&Bc::Dirichlet(BcValue::spacetime(|t, x| t * x[0])));
    let g = ChebGrid::new(&[8], &[0.0], &[1.0]).unwrap();
    let src = |_t: f64, x: [f64; 3]| x[0];
    let u = g.heat(&bcs, 1.0, |_| 0.0, Some(&src), 1.0, 50, 0.5).unwrap();
    assert!(g.max_error(&u, |x| x[0]) < 1e-10);
}
