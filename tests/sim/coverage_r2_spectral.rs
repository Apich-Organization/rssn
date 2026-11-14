//! Coverage round 2: spectral helpers (wavenumbers, complex transforms, differentiation
//! matrices, dealiasing, grid validation).

use std::f64::consts::PI;

use num_complex::Complex;
use rssn::sim::boundary::FaceBcs;
use rssn::sim::spectral::ChebGrid;
use rssn::sim::spectral::SpectralGrid;
use rssn::sim::spectral::cheb_diff_matrix;
use rssn::sim::spectral::fft_nd;

fn c(
    re: f64,
    im: f64,
) -> Complex<f64> {
    Complex::new(re, im)
}

#[test]
fn grids_reject_invalid_shapes() {
    assert!(SpectralGrid::new(&[], &[]).is_err());
    assert!(SpectralGrid::new(&[8, 8, 8, 8], &[1.0; 4]).is_err());
    assert!(SpectralGrid::new(&[8, 8], &[1.0]).is_err());
    assert!(
        SpectralGrid::new(&[6], &[1.0]).is_err(),
        "not a power of two"
    );
    assert!(
        SpectralGrid::new(&[2], &[1.0]).is_err(),
        "fewer than four points"
    );
    assert!(SpectralGrid::new(&[8, 3], &[1.0, 1.0]).is_err());
    assert!(SpectralGrid::new(&[4, 8, 16], &[1.0, 2.0, 3.0]).is_ok());
    assert!(ChebGrid::new(&[], &[], &[]).is_err());
    assert!(
        ChebGrid::new(&[1], &[0.0], &[1.0]).is_err(),
        "degree below 2"
    );
    assert!(ChebGrid::new(&[4], &[0.0, 0.0], &[1.0]).is_err());
    assert!(ChebGrid::new(&[4, 4, 4, 4], &[0.0; 4], &[1.0; 4]).is_err());
    let g = ChebGrid::new(&[4, 6], &[0.0, -1.0], &[2.0, 3.0]).unwrap();
    assert_eq!(g.size(), 5 * 7);
}

#[test]
fn fft_nd_validates_extents_and_is_unitary_up_to_scale() {
    assert!(
        fft_nd(&mut [c(0.0, 0.0); 8], &[4], false).is_err(),
        "length mismatch"
    );
    assert!(
        fft_nd(&mut [c(0.0, 0.0); 12], &[3, 4], true).is_err(),
        "3 is not a power of two"
    );
    // forward then inverse is the identity
    let n = [8usize, 4];
    let orig: Vec<Complex<f64>> = (0..32)
        .map(|p| c((p as f64).sin(), (p as f64 * 0.7).cos()))
        .collect();
    let mut w = orig.clone();
    fft_nd(&mut w, &n, false).unwrap();
    // Parseval: sum |F|^2 = N sum |f|^2
    let e0: f64 = orig.iter().map(|v| v.norm_sqr()).sum();
    let e1: f64 = w.iter().map(|v| v.norm_sqr()).sum();
    assert!((e1 - 32.0 * e0).abs() < 1e-9 * e1);
    fft_nd(&mut w, &n, true).unwrap();
    for (a, b) in w.iter().zip(&orig) {
        assert!((a - b).norm() < 1e-12);
    }
    // linearity
    let x: Vec<Complex<f64>> = (0..16).map(|p| c(p as f64, 1.0)).collect();
    let y: Vec<Complex<f64>> = (0..16).map(|p| c(1.0, (p * p) as f64 * 0.01)).collect();
    let (mut fx, mut fy, mut fz) = (
        x.clone(),
        y.clone(),
        x.iter()
            .zip(&y)
            .map(|(a, b)| a * 2.0 + b * c(0.0, 3.0))
            .collect::<Vec<_>>(),
    );
    fft_nd(&mut fx, &[16], false).unwrap();
    fft_nd(&mut fy, &[16], false).unwrap();
    fft_nd(&mut fz, &[16], false).unwrap();
    for k in 0..16 {
        assert!((fz[k] - (fx[k] * 2.0 + fy[k] * c(0.0, 3.0))).norm() < 1e-9);
    }
}

#[test]
fn wavenumbers_wrap_to_negative_frequencies_and_kvec_agrees() {
    let g = SpectralGrid::new(&[8, 4], &[2.0 * PI, 4.0 * PI]).unwrap();
    let w: Vec<f64> = (0..8).map(|j| g.wavenumber(0, j)).collect();
    assert_eq!(w, vec![0.0, 1.0, 2.0, 3.0, 4.0, -3.0, -2.0, -1.0]);
    // second axis has length 4 pi so the spacing of k is 1/2
    let w1: Vec<f64> = (0..4).map(|j| g.wavenumber(1, j)).collect();
    assert_eq!(w1, vec![0.0, 0.5, 1.0, -0.5]);
    for p in 0..g.size() {
        let k = g.kvec(p);
        let (i, j) = (p % 8, (p / 8) % 4);
        assert_eq!(k[0], g.wavenumber(0, i));
        assert_eq!(k[1], g.wavenumber(1, j));
        assert_eq!(k[2], 0.0, "unused axis has zero wavenumber");
        assert!((g.k2(p) - (k[0] * k[0] + k[1] * k[1])).abs() < 1e-15);
        let x = g.coord(p);
        assert!(
            (x[0] - i as f64 * g.h(0)).abs() < 1e-15 && (x[1] - j as f64 * g.h(1)).abs() < 1e-15
        );
    }
    assert_eq!(g.size(), 32);
}

#[test]
fn complex_transforms_round_trip_and_real_input_has_hermitian_spectrum() {
    let g = SpectralGrid::new(&[8, 8], &[1.0, 1.0]).unwrap();
    let u: Vec<Complex<f64>> = (0..64)
        .map(|p| c((p as f64 * 0.3).sin(), (p as f64 * 0.11).cos() - 0.5))
        .collect();
    let back = g.to_phys_c(&g.to_spec_c(&u));
    for (a, b) in back.iter().zip(&u) {
        assert!((a - b).norm() < 1e-12);
    }
    // real input: hat u(-k) = conj(hat u(k))
    let r: Vec<f64> = (0..64).map(|p| ((p * 7 % 13) as f64).sin()).collect();
    let rh = g.to_spec(&r);
    for j in 0..8 {
        for i in 0..8 {
            let (mi, mj) = ((8 - i) % 8, (8 - j) % 8);
            assert!((rh[i + 8 * j] - rh[mi + 8 * mj].conj()).norm() < 1e-10);
        }
    }
    // to_phys of a real signal's spectrum returns the signal
    let rr = g.to_phys(&rh);
    for (a, b) in rr.iter().zip(&r) {
        assert!((a - b).abs() < 1e-12);
    }
    // complex derivative: d/dx of exp(i m x) is i m exp(i m x)
    let m = 3.0;
    let e: Vec<Complex<f64>> = (0..64)
        .map(|p| Complex::new(0.0, m * 2.0 * PI * g.coord(p)[0]).exp())
        .collect();
    let de = g.to_phys_c(&g.deriv(&g.to_spec_c(&e), 0));
    for (p, v) in de.iter().enumerate() {
        let want = c(0.0, m * 2.0 * PI) * e[p];
        assert!((v - want).norm() < 1e-9, "p = {p}");
    }
}

#[test]
fn derivative_kills_the_nyquist_mode_and_the_mean() {
    let g = SpectralGrid::new(&[8], &[2.0 * PI]).unwrap();
    // the Nyquist mode (-1)^j has no well-defined derivative on the grid: it is zeroed
    let nyq: Vec<f64> = (0..8)
        .map(|j| if j % 2 == 0 { 1.0 } else { -1.0 })
        .collect();
    let d = g.to_phys(&g.deriv(&g.to_spec(&nyq), 0));
    assert!(d.iter().all(|v| v.abs() < 1e-12));
    let k = vec![5.0; 8];
    let d = g.to_phys(&g.deriv(&g.to_spec(&k), 0));
    assert!(d.iter().all(|v| v.abs() < 1e-12));
    // the second derivative is minus k^2 on a resolved mode
    let u: Vec<f64> = (0..8).map(|j| (2.0 * g.coord(j)[0]).sin()).collect();
    let uh = g.to_spec(&u);
    let d2 = g.to_phys(&g.deriv(&g.deriv(&uh, 0), 0));
    for j in 0..8 {
        assert!((d2[j] + 4.0 * u[j]).abs() < 1e-11);
    }
}

#[test]
fn dealias_mask_keeps_the_lowest_two_thirds_of_each_axis() {
    let g = SpectralGrid::new(&[16], &[1.0]).unwrap();
    let mask = g.dealias_mask();
    assert_eq!(mask.len(), 16);
    // |m| < 16/3 -> m in 0..=5, plus mirrored negative modes
    let kept: Vec<usize> = (0..16).filter(|&j| mask[j]).collect();
    assert_eq!(kept, vec![0, 1, 2, 3, 4, 5, 11, 12, 13, 14, 15]);
    // mask is symmetric under k -> -k
    for j in 1..16 {
        assert_eq!(mask[j], mask[16 - j]);
    }
    let g2 = SpectralGrid::new(&[8, 8], &[1.0, 1.0]).unwrap();
    let m2 = g2.dealias_mask();
    assert_eq!(
        m2.iter().filter(|&&b| b).count(),
        5 * 5,
        "tensor product of 5 retained modes"
    );
}

#[test]
fn helmholtz_solves_manufactured_solutions_and_rejects_singular_forcing() {
    let g = SpectralGrid::new(&[16, 16], &[2.0 * PI, 2.0 * PI]).unwrap();
    let u = |x: [f64; 3]| (2.0 * x[0]).sin() * (x[1]).cos() + 0.5 * (3.0 * x[1]).sin();
    let sigma = 1.5;
    // (k^2 + sigma) uhat = fhat: -lap u + sigma u = f
    let f: Vec<f64> = (0..g.size())
        .map(|p| {
            let x = g.coord(p);
            let lapu = -4.0 * (2.0 * x[0]).sin() * x[1].cos()
                - (2.0 * x[0]).sin() * x[1].cos()
                - 4.5 * (3.0 * x[1]).sin();
            -lapu + sigma * u(x)
        })
        .collect();
    let sol = g.helmholtz(sigma, &f).unwrap();
    for (p, v) in sol.iter().enumerate() {
        assert!((v - u(g.coord(p))).abs() < 1e-10);
    }
    // sigma = 0: the constant mode is singular. Zero-mean forcing is fine and returns the
    // zero-mean solution; a forcing with mean is an error.
    let zero_mean: Vec<f64> = (0..g.size()).map(|p| (g.coord(p)[0]).sin()).collect();
    let s0 = g.helmholtz(0.0, &zero_mean).unwrap();
    for (p, v) in s0.iter().enumerate() {
        assert!((v - g.coord(p)[0].sin()).abs() < 1e-10);
    }
    assert!(g.helmholtz(0.0, &vec![1.0; g.size()]).is_err());
    // negative sigma = -k^2 of a resolved mode is singular for that mode
    let mode: Vec<f64> = (0..g.size()).map(|p| g.coord(p)[0].cos()).collect();
    assert!(g.helmholtz(-1.0, &mode).is_err());
}

#[test]
fn chebyshev_matrix_differentiates_polynomials_exactly() {
    for &n in &[4usize, 8, 12] {
        for &len in &[2.0, 5.0] {
            let d = cheb_diff_matrix(n, len);
            assert_eq!(d.len(), n + 1);
            // nodes are y_i = (len/2) * (-cos(pi i / n))
            let y: Vec<f64> = (0..=n)
                .map(|i| -0.5 * len * (PI * i as f64 / n as f64).cos())
                .collect();
            for row in &d {
                assert_eq!(row.len(), n + 1);
                let s: f64 = row.iter().sum();
                assert!(
                    s.abs() < 1e-9 * n as f64 * n as f64,
                    "rows sum to zero: {s}"
                );
            }
            for k in 0..=n.min(7) {
                let f: Vec<f64> = y.iter().map(|v| v.powi(k as i32)).collect();
                for i in 0..=n {
                    let got: f64 = (0..=n).map(|j| d[i][j] * f[j]).sum();
                    let want = if k == 0 {
                        0.0
                    } else {
                        k as f64 * y[i].powi(k as i32 - 1)
                    };
                    let scale = 1.0 + want.abs();
                    assert!(
                        (got - want).abs() < 1e-8 * scale * (n * n) as f64,
                        "n={n} len={len} k={k} i={i}: {got} vs {want}"
                    );
                }
            }
            // antisymmetry of the matrix under node reflection: D[i][j] = -D[n-i][n-j]
            for i in 0..=n {
                for j in 0..=n {
                    assert!((d[i][j] + d[n - i][n - j]).abs() < 1e-8 * (n * n) as f64);
                }
            }
        }
    }
    // the exponential is differentiated with spectral accuracy
    let n = 14;
    let d = cheb_diff_matrix(n, 2.0);
    let y: Vec<f64> = (0..=n).map(|i| -(PI * i as f64 / n as f64).cos()).collect();
    for i in 0..=n {
        let got: f64 = (0..=n).map(|j| d[i][j] * y[j].exp()).sum();
        assert!((got - y[i].exp()).abs() < 1e-10);
    }
}

#[test]
fn chebyshev_grid_nodes_sample_and_error_norm() {
    let g = ChebGrid::new(&[6, 4], &[-1.0, 2.0], &[3.0, 6.0]).unwrap();
    // endpoints are nodes, nodes ascend and mirror about the centre
    assert!((g.node(0, 0) + 1.0).abs() < 1e-14 && (g.node(0, 6) - 3.0).abs() < 1e-14);
    assert!((g.node(1, 0) - 2.0).abs() < 1e-14 && (g.node(1, 4) - 6.0).abs() < 1e-14);
    for i in 0..6 {
        assert!(g.node(0, i) < g.node(0, i + 1));
        assert!((g.node(0, i) + g.node(0, 6 - i) - 2.0).abs() < 1e-13);
    }
    let f = |x: [f64; 3]| x[0] * x[1] + 1.0;
    let u = g.sample(f);
    assert_eq!(u.len(), g.size());
    assert_eq!(g.max_error(&u, f), 0.0);
    let off = |x: [f64; 3]| f(x) + 0.25;
    assert!((g.max_error(&u, off) - 0.25).abs() < 1e-14);
    // coordinates fill the tensor grid in x-fastest order
    assert_eq!(g.coord(0)[0], g.node(0, 0));
    assert_eq!(g.coord(1)[0], g.node(0, 1));
    assert_eq!(g.coord(7)[1], g.node(1, 1));
}

#[test]
fn chebyshev_heat_flow_obeys_the_maximum_principle_and_decays() {
    let g = ChebGrid::new(&[16], &[0.0], &[1.0]).unwrap();
    let bcs = FaceBcs::dirichlet_zero();
    let init = |x: [f64; 3]| (PI * x[0]).sin() + 0.3 * (3.0 * PI * x[0]).sin();
    let max0 = g.sample(init).iter().fold(0.0f64, |a, b| a.max(b.abs()));
    let mut prev = max0;
    for &t in &[0.01, 0.05, 0.2] {
        let u = g.heat(&bcs, 1.0, init, None, t, 40, 1.0).unwrap();
        let m = u.iter().fold(0.0f64, |a, b| a.max(b.abs()));
        assert!(m <= prev + 1e-12, "monotone decay: {m} > {prev}");
        prev = m;
        // boundaries hold
        assert!(u[0].abs() < 1e-12 && u[16].abs() < 1e-12);
    }
    // the slowest mode dominates: exp(-pi^2 t) amplitude
    let u = g.heat(&bcs, 1.0, init, None, 0.2, 200, 0.5).unwrap();
    let mid = u[8];
    let want =
        (-PI * PI * 0.2).exp() * 1.0 + 0.3 * (-9.0 * PI * PI * 0.2).exp() * (3.0 * PI * 0.5).sin();
    assert!((mid - want).abs() < 1e-3, "{mid} vs {want}");
    // unsupported boundary type is an error, not a panic
    assert!(
        g.heat(
            &FaceBcs::uniform(&rssn::sim::boundary::BoundaryCondition::Periodic),
            1.0,
            init,
            None,
            0.1,
            4,
            1.0
        )
        .is_err()
    );
}
