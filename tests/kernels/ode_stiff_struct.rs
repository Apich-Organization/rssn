//! Tests for banded/sparse Jacobians, colouring finite differences and
//! index-2/3 DAEs in the Radau and BDF integrators.

use rssn::kernels::dense::Mat;
use rssn::kernels::jacobian::{
    BandLu, JacStructure, band_pattern, color_columns, fd_jacobian_colored,
};
use rssn::kernels::krylov::Csr;
use rssn::kernels::ode_adaptive::{OdeOptions, dopri5};
use rssn::kernels::ode_stiff::{
    BdfOptions, DaeIndex, RadauOptions, bdf_with, radau_with,
};

fn heat(n: usize) -> impl Fn(f64, &[f64], &mut [f64]) {
    let h = 1.0 / (n as f64 + 1.0);
    move |_t, u, du| {
        for i in 0..n {
            let l = if i > 0 { u[i - 1] } else { 0.0 };
            let r = if i + 1 < n { u[i + 1] } else { 0.0 };
            du[i] = (l - 2.0 * u[i] + r) / (h * h);
        }
    }
}

fn heat_exact(n: usize, t: f64) -> Vec<f64> {
    let h = 1.0 / (n as f64 + 1.0);
    let lam = (2.0 - 2.0 * (std::f64::consts::PI * h).cos()) / (h * h);
    (1..=n).map(|i| (-lam * t).exp() * (std::f64::consts::PI * i as f64 * h).sin()).collect()
}

fn heat_u0(n: usize) -> Vec<f64> {
    heat_exact(n, 0.0)
}

#[test]
fn band_lu_matches_dense_solve() {
    let n = 40;
    let (kl, ku) = (3, 2);
    let mut trip = Vec::new();
    let mut m = Mat::zeros(n, n);
    let mut seed = 12345u64;
    let mut rnd = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((seed >> 33) as f64) / f64::from(1u32 << 31) - 0.5
    };
    for i in 0..n {
        for j in i.saturating_sub(kl)..=(i + ku).min(n - 1) {
            let v = rnd();
            trip.push((i, j, v));
            m.set(i, j, v);
        }
    }
    let b: Vec<f64> = (0..n).map(|i| (i as f64).sin()).collect();
    let lu = BandLu::factor(n, kl, ku, &trip).unwrap();
    let x = lu.solve(&b);
    let want = rssn::kernels::dense::solve(&m, &b).unwrap();
    for i in 0..n {
        assert!((x[i] - want[i]).abs() < 1e-9 * (1.0 + want[i].abs()), "{i}");
    }
    // outside the band: rejected
    assert!(BandLu::factor(4, 1, 1, &[(0, 3, 1.0)]).is_none());
}

#[test]
fn colouring_uses_band_width_colours() {
    let pat = band_pattern(50, 1, 1);
    let (colors, k) = color_columns(&pat);
    assert_eq!(k, 3);
    // columns of one colour share no row
    for i in 0..50usize {
        let cols: Vec<usize> = (i.saturating_sub(1)..=(i + 1).min(49)).collect();
        for a in &cols {
            for b in &cols {
                assert!(a == b || colors[*a] != colors[*b]);
            }
        }
    }
    let f = heat(50);
    let y = heat_u0(50);
    let mut f0 = vec![0.0; 50];
    f(0.0, &y, &mut f0);
    let mut nfev = 0;
    let jac = fd_jacobian_colored(&f, 0.0, &y, &f0, &pat, &(colors, k), &mut nfev);
    assert_eq!(nfev, 3);
    let h2 = (1.0 / 51.0_f64).powi(2);
    for i in 0..50 {
        for p in jac.indptr[i]..jac.indptr[i + 1] {
            let j = jac.indices[p];
            let want = if i == j { -2.0 / h2 } else { 1.0 / h2 };
            assert!((jac.values[p] - want).abs() < 1e-3 * want.abs(), "({i},{j})");
        }
    }
}

#[test]
fn colouring_of_irregular_pattern() {
    // arrow matrix: dense first row and column force many colours only for the
    // first column; the rest share colours
    let n = 20;
    let mut trip = Vec::new();
    for i in 0..n {
        trip.push((i, i, 1.0));
        trip.push((i, 0, 1.0));
        trip.push((0, i, 1.0));
    }
    let pat = Csr::from_triplets(n, n, &trip);
    let (colors, k) = color_columns(&pat);
    // every column meets row 0, so all columns need distinct colours
    assert_eq!(k, n);
    assert_eq!(colors.len(), n);
}

#[test]
fn radau_banded_and_sparse_match_dense() {
    let n = 300;
    let f = heat(n);
    let u0 = heat_u0(n);
    let t1 = 0.002;
    let opts = OdeOptions { rtol: 1e-9, atol: 1e-12, ..OdeOptions::default() };
    let three = RadauOptions { max_stages: 3, ..RadauOptions::default() };
    let dense = radau_with(&f, 0.0, t1, &u0, &opts, &three).unwrap();
    let banded_opts = RadauOptions { structure: JacStructure::Banded { kl: 1, ku: 1 }, ..three };
    let banded = radau_with(&f, 0.0, t1, &u0, &opts, &banded_opts).unwrap();
    let pat = band_pattern(n, 1, 1);
    let sparse_opts =
        RadauOptions { structure: JacStructure::Sparse(&pat), ..three };
    let sparse = radau_with(&f, 0.0, t1, &u0, &opts, &sparse_opts).unwrap();
    let want = heat_exact(n, t1);
    for sol in [&dense, &banded, &sparse] {
        let y = sol.y.last().unwrap();
        for i in 0..n {
            assert!((y[i] - want[i]).abs() < 1e-7, "{i}: {} vs {}", y[i], want[i]);
        }
    }
    // the colouring needs 3 evaluations per Jacobian instead of n
    assert!(banded.nfev * 2 < dense.nfev, "{} vs {}", banded.nfev, dense.nfev);
    assert!(sparse.nfev * 2 < dense.nfev);
}

#[test]
fn radau_user_sparse_jacobian_and_large_system() {
    let n = 2000;
    let f = heat(n);
    let h2 = (1.0 / (n as f64 + 1.0)).powi(2);
    let sj = move |_t: f64, _y: &[f64]| {
        let mut trip = Vec::new();
        for i in 0..n {
            trip.push((i, i, -2.0 / h2));
            if i > 0 {
                trip.push((i, i - 1, 1.0 / h2));
            }
            if i + 1 < n {
                trip.push((i, i + 1, 1.0 / h2));
            }
        }
        Csr::from_triplets(n, n, &trip)
    };
    let u0 = heat_u0(n);
    let opts = OdeOptions { rtol: 1e-8, atol: 1e-11, ..OdeOptions::default() };
    let ro = RadauOptions { sparse_jac: Some(&sj), max_stages: 3, ..RadauOptions::default() };
    let sol = radau_with(&f, 0.0, 1e-4, &u0, &opts, &ro).unwrap();
    let want = heat_exact(n, 1e-4);
    let y = sol.y.last().unwrap();
    let err = (0..n).map(|i| (y[i] - want[i]).abs()).fold(0.0, f64::max);
    assert!(err < 1e-6, "err {err}");
}

#[test]
fn bdf_banded_and_sparse() {
    let n = 80;
    let f = heat(n);
    let u0 = heat_u0(n);
    let t1 = 0.01;
    let opts = OdeOptions { rtol: 1e-8, atol: 1e-11, ..OdeOptions::default() };
    let want = heat_exact(n, t1);
    let dense = bdf_with(&f, 0.0, t1, &u0, &opts, &BdfOptions::default()).unwrap();
    let bo = BdfOptions { structure: JacStructure::Banded { kl: 1, ku: 1 }, ..BdfOptions::default() };
    let banded = bdf_with(&f, 0.0, t1, &u0, &opts, &bo).unwrap();
    let pat = band_pattern(n, 1, 1);
    let so = BdfOptions { structure: JacStructure::Sparse(&pat), ..BdfOptions::default() };
    let sparse = bdf_with(&f, 0.0, t1, &u0, &opts, &so).unwrap();
    for sol in [&dense, &banded, &sparse] {
        let y = sol.y.last().unwrap();
        for i in 0..n {
            assert!((y[i] - want[i]).abs() < 5e-6, "{i}: {} vs {}", y[i], want[i]);
        }
    }
    assert!(banded.nfev < dense.nfev);
}

const G: f64 = 9.81;

/// Pendulum with unit mass and length: position `(x, y)`, velocity `(u, v)`,
/// multiplier `lambda`; `index` selects the constraint level.
fn pendulum(index3: bool) -> impl Fn(f64, &[f64], &mut [f64]) {
    move |_t, w, dw| {
        let (x, y, u, v, l) = (w[0], w[1], w[2], w[3], w[4]);
        dw[0] = u;
        dw[1] = v;
        dw[2] = -2.0 * x * l;
        dw[3] = -G - 2.0 * y * l;
        dw[4] = if index3 { x * x + y * y - 1.0 } else { x * u + y * v };
    }
}

fn pendulum_reference(t: f64) -> [f64; 4] {
    // angle formulation theta'' = -g cos(theta), x = cos, y = sin
    let ang = |_t: f64, s: &[f64], ds: &mut [f64]| {
        ds[0] = s[1];
        ds[1] = -G * s[0].cos();
    };
    let opts = OdeOptions { rtol: 1e-13, atol: 1e-13, ..OdeOptions::default() };
    let sol = dopri5(ang, 0.0, t, &[0.0, 0.0], &opts, &[]).unwrap();
    let s = sol.y.last().unwrap();
    let (th, om) = (s[0], s[1]);
    [th.cos(), th.sin(), -om * th.sin(), om * th.cos()]
}

fn pendulum_mass() -> Mat {
    let mut m = Mat::identity(5);
    m.set(4, 4, 0.0);
    m
}

#[test]
fn radau_pendulum_index3() {
    let m = pendulum_mass();
    let f = pendulum(true);
    let opts = OdeOptions { rtol: 1e-8, atol: 1e-8, ..OdeOptions::default() };
    let y0 = [1.0, 0.0, 0.0, 0.0, 0.0];
    let ro = RadauOptions {
        mass: Some(&m),
        index: DaeIndex { n2: 2, n3: 1 },
        ..RadauOptions::default()
    };
    let sol = radau_with(&f, 0.0, 1.0, &y0, &opts, &ro).unwrap();
    let y = sol.y.last().unwrap();
    let r = pendulum_reference(1.0);
    for i in 0..4 {
        assert!((y[i] - r[i]).abs() < 2e-5, "component {i}: {} vs {}", y[i], r[i]);
    }
    assert!((y[0] * y[0] + y[1] * y[1] - 1.0).abs() < 1e-8);
    // the multiplier equals (u^2 + v^2 - g y) / 2 on the manifold
    let lam = (y[2] * y[2] + y[3] * y[3] - G * y[1]) / 2.0;
    assert!((y[4] - lam).abs() < 1e-3 * lam.abs().max(1.0), "{} vs {lam}", y[4]);
    // index-aware error control takes far fewer rejected steps than treating
    // all variables as index 1
    let plain = RadauOptions { mass: Some(&m), ..RadauOptions::default() };
    match radau_with(&f, 0.0, 1.0, &y0, &opts, &plain) {
        Ok(s) => assert!(
            s.rejected > sol.rejected || s.t.len() > sol.t.len(),
            "aware: {} steps/{} rej; plain: {} steps/{} rej",
            sol.t.len(),
            sol.rejected,
            s.t.len(),
            s.rejected
        ),
        Err(_) => {}
    }
}

#[test]
fn radau_pendulum_index2() {
    let m = pendulum_mass();
    let f = pendulum(false);
    let opts = OdeOptions { rtol: 1e-9, atol: 1e-9, ..OdeOptions::default() };
    let y0 = [1.0, 0.0, 0.0, 0.0, 0.0];
    let ro = RadauOptions {
        mass: Some(&m),
        index: DaeIndex { n2: 1, n3: 0 },
        ..RadauOptions::default()
    };
    let sol = radau_with(&f, 0.0, 0.5, &y0, &opts, &ro).unwrap();
    let y = sol.y.last().unwrap();
    let r = pendulum_reference(0.5);
    for i in 0..4 {
        assert!((y[i] - r[i]).abs() < 5e-5, "component {i}: {} vs {}", y[i], r[i]);
    }
}

#[test]
fn radau_hessenberg_index2_linear() {
    // y' = -y + z, 0 = y - cos t  =>  y = cos t, z = cos t - sin t
    let f = |t: f64, w: &[f64], dw: &mut [f64]| {
        dw[0] = -w[0] + w[1];
        dw[1] = w[0] - t.cos();
    };
    let m = Mat::from_rows(&[vec![1.0, 0.0], vec![0.0, 0.0]]).unwrap();
    let opts = OdeOptions { rtol: 1e-9, atol: 1e-10, ..OdeOptions::default() };
    let ro = RadauOptions {
        mass: Some(&m),
        index: DaeIndex { n2: 1, n3: 0 },
        ..RadauOptions::default()
    };
    let sol = radau_with(f, 0.0, 3.0, &[1.0, 1.0], &opts, &ro).unwrap();
    for (k, t) in sol.t.iter().enumerate() {
        let w = &sol.y[k];
        assert!((w[0] - t.cos()).abs() < 1e-7);
        assert!((w[1] - (t.cos() - t.sin())).abs() < 1e-6, "t={t}: {}", w[1]);
    }
}

#[test]
fn radau_sparse_mass_dae() {
    let f = |_t: f64, u: &[f64], du: &mut [f64]| {
        du[0] = -u[0] + u[1];
        du[1] = u[1] - u[0] * u[0];
    };
    let m = Csr::from_triplets(2, 2, &[(0, 0, 1.0)]);
    let opts = OdeOptions { rtol: 1e-9, atol: 1e-11, ..OdeOptions::default() };
    let pat = Csr::from_triplets(2, 2, &[(0, 0, 1.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 1.0)]);
    let ro = RadauOptions {
        sparse_mass: Some(&m),
        structure: JacStructure::Sparse(&pat),
        ..RadauOptions::default()
    };
    let sol = radau_with(f, 0.0, 4.0, &[0.5, 0.25], &opts, &ro).unwrap();
    let last = sol.y.last().unwrap();
    let exact = 1.0 / (1.0 + 4.0_f64.exp());
    assert!((last[0] - exact).abs() < 1e-7);
}

#[test]
fn invalid_structures_are_rejected() {
    let f = |_t: f64, y: &[f64], dy: &mut [f64]| dy[0] = -y[0];
    let opts = OdeOptions::default();
    let pat = band_pattern(3, 1, 1);
    let ro = RadauOptions { structure: JacStructure::Sparse(&pat), ..RadauOptions::default() };
    assert!(radau_with(f, 0.0, 1.0, &[1.0], &opts, &ro).is_err());
    let ro = RadauOptions { index: DaeIndex { n2: 2, n3: 0 }, ..RadauOptions::default() };
    assert!(radau_with(f, 0.0, 1.0, &[1.0], &opts, &ro).is_err());
}
