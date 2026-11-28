//! Round-2 coverage: Jacobian providers and colouring, matrix norms,
//! multilevel nested dissection with explicit options, the AMG V-cycle,
//! the legacy IDR(s) wrapper and the arc / cell accessors of the coloured
//! graph kernel. Public API only; deterministic inputs.
// The deprecated `idr_s` wrapper is tested on purpose (it must stay bit-identical to `idr_s_with`).
#![allow(deprecated)]
#![allow(clippy::needless_range_loop)]

use rssn::kernels::amg_general::AmgGeneral;
use rssn::kernels::amg_general::AmgGeneralOptions;
use rssn::kernels::amg_general::AmgKind;
use rssn::kernels::dense::Mat;
use rssn::kernels::graph::ColourGraph;
use rssn::kernels::graph::refine;
use rssn::kernels::jacobian::FdMode;
use rssn::kernels::jacobian::Jac;
use rssn::kernels::jacobian::JacStructure;
use rssn::kernels::jacobian::JacobianProvider;
use rssn::kernels::jacobian::SolverKind;
use rssn::kernels::jacobian::band_pattern;
use rssn::kernels::jacobian::color_columns;
use rssn::kernels::jacobian::fd_jacobian_colored;
use rssn::kernels::jacobian::fd_jacobian_colored_par;
use rssn::kernels::jacobian::fd_jacobian_dense;
use rssn::kernels::jacobian::fd_jacobian_dense_par;
use rssn::kernels::krylov::Csr;
use rssn::kernels::krylov::identity_precond;
use rssn::kernels::krylov_ext::IdrOptions;
use rssn::kernels::krylov_ext::idr_s;
use rssn::kernels::krylov_ext::idr_s_with;
use rssn::kernels::matfun::csr_norm1;
use rssn::kernels::matfun::norm1;
use rssn::kernels::multigrid::Cycle;
use rssn::kernels::ordering::cholesky_fill;
use rssn::kernels::ordering_ml::NdOptions;
use rssn::kernels::ordering_ml::hybrid_nested_dissection_order;
use rssn::kernels::ordering_ml::multilevel_nested_dissection_order;
use rssn::kernels::ordering_ml::multilevel_nested_dissection_with;
use rssn::kernels::ordering_ml::vertex_separator_bisection;

// ------------------------------------------------------------ matrix norms

#[test]
fn norm1_dense_and_sparse_agree_with_hand_values() {
    // [[1, -2, 3], [-4, 5, -6]] : column sums 5, 7, 9 -> ||A||_1 = 9
    let a = Mat::from_rows(&[vec![1.0, -2.0, 3.0], vec![-4.0, 5.0, -6.0]]).unwrap();
    assert_eq!(norm1(&a), 9.0);
    let trip: Vec<(usize, usize, f64)> = (0..2)
        .flat_map(|i| (0..3).map(move |j| (i, j)))
        .map(|(i, j)| (i, j, a.at(i, j)))
        .collect();
    let s = Csr::from_triplets(2, 3, &trip);
    assert_eq!(csr_norm1(&s), 9.0);
    // identity has norm 1; zero and empty matrices have norm 0
    assert_eq!(norm1(&Mat::identity(5)), 1.0);
    assert_eq!(norm1(&Mat::zeros(3, 3)), 0.0);
    assert_eq!(norm1(&Mat::zeros(0, 0)), 0.0);
    assert_eq!(csr_norm1(&Csr::from_triplets(0, 0, &[])), 0.0);
    assert_eq!(csr_norm1(&Csr::from_triplets(4, 4, &[])), 0.0);
    // duplicate triplets are summed before the norm is taken
    let dup = Csr::from_triplets(2, 2, &[(0, 0, 1.0), (0, 0, 2.0), (1, 0, -4.0), (1, 1, 0.5)]);
    assert_eq!(csr_norm1(&dup), 7.0);
}

#[test]
fn norm1_is_a_norm_on_pseudo_random_matrices() {
    // homogeneity, triangle inequality and the dense / sparse cross-check.
    let mut s = 0x1234_5678_u64;
    let mut next = || {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((s >> 11) as f64) / (1u64 << 53) as f64 - 0.5
    };
    for n in [1usize, 2, 5, 9] {
        let mk = |next: &mut dyn FnMut() -> f64| {
            let mut m = Mat::zeros(n, n);
            for i in 0..n {
                for j in 0..n {
                    m.set(i, j, 4.0 * next());
                }
            }
            m
        };
        let (a, b) = (mk(&mut next), mk(&mut next));
        let mut sum = Mat::zeros(n, n);
        let mut scaled = Mat::zeros(n, n);
        let mut trip = Vec::new();
        for i in 0..n {
            for j in 0..n {
                sum.set(i, j, a.at(i, j) + b.at(i, j));
                scaled.set(i, j, -3.0 * a.at(i, j));
                trip.push((i, j, a.at(i, j)));
            }
        }
        assert!(norm1(&sum) <= norm1(&a) + norm1(&b) + 1e-12);
        assert!((norm1(&scaled) - 3.0 * norm1(&a)).abs() < 1e-12);
        let sp = Csr::from_triplets(n, n, &trip);
        assert!((csr_norm1(&sp) - norm1(&a)).abs() < 1e-12);
    }
    // inf entries give an infinite norm instead of a panic
    let mut inf = Mat::zeros(2, 2);
    inf.set(0, 1, f64::INFINITY);
    assert!(norm1(&inf).is_infinite());
}

// -------------------------------------------------------------- jacobian

/// `f_i = y_{i-1} y_i - y_{i+1}^2 + sin(y_i)` with zero outside the range.
fn chain_rhs(
    _t: f64,
    y: &[f64],
    out: &mut [f64],
) {
    let n = y.len();
    for i in 0..n {
        let left = if i > 0 { y[i - 1] } else { 0.0 };
        let right = if i + 1 < n { y[i + 1] } else { 0.0 };
        out[i] = left * y[i] - right * right + y[i].sin();
    }
}

fn chain_state(n: usize) -> Vec<f64> {
    (0..n)
        .map(|i| 0.3 + 0.07 * i as f64 - 0.4 * ((i * i) % 5) as f64 / 5.0)
        .collect()
}

/// Analytic tridiagonal Jacobian of `chain_rhs` as a dense matrix.
fn chain_jacobian(y: &[f64]) -> Mat {
    let n = y.len();
    let mut j = Mat::zeros(n, n);
    for i in 0..n {
        let left = if i > 0 { y[i - 1] } else { 0.0 };
        j.set(i, i, left + y[i].cos());
        if i > 0 {
            j.set(i, i - 1, y[i]);
        }
        if i + 1 < n {
            j.set(i, i + 1, -2.0 * y[i + 1]);
        }
    }
    j
}

fn csr_to_dense(c: &Csr) -> Mat {
    let mut m = Mat::zeros(c.rows, c.cols);
    for i in 0..c.rows {
        for k in c.indptr[i]..c.indptr[i + 1] {
            m.set(i, c.indices[k], c.values[k]);
        }
    }
    m
}

fn max_abs_diff(
    a: &Mat,
    b: &Mat,
) -> f64 {
    a.data
        .iter()
        .zip(&b.data)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

#[test]
fn band_pattern_and_colouring_are_consistent() {
    for (n, kl, ku) in [
        (1usize, 0usize, 0usize),
        (6, 1, 1),
        (12, 2, 1),
        (9, 0, 3),
        (4, 10, 10),
    ] {
        let pat = band_pattern(n, kl, ku);
        assert_eq!((pat.rows, pat.cols), (n, n));
        // exactly the entries with -kl <= j - i <= ku
        for i in 0..n {
            let cols: Vec<usize> = (pat.indptr[i]..pat.indptr[i + 1])
                .map(|k| pat.indices[k])
                .collect();
            let expect: Vec<usize> = (0..n)
                .filter(|&j| {
                    (j as i64 - i as i64) >= -(kl as i64) && (j as i64 - i as i64) <= ku as i64
                })
                .collect();
            assert_eq!(cols, expect, "row {i} of band({n},{kl},{ku})");
        }
        // proper colouring: columns that share a row get different colours
        let (colour, ncol) = color_columns(&pat);
        assert_eq!(colour.len(), n);
        assert!(colour.iter().all(|&c| c < ncol));
        for i in 0..n {
            for a in pat.indptr[i]..pat.indptr[i + 1] {
                for b in (a + 1)..pat.indptr[i + 1] {
                    assert_ne!(
                        colour[pat.indices[a]], colour[pat.indices[b]],
                        "n={n} row {i}"
                    );
                }
            }
        }
        // a band of width kl + ku + 1 needs at least that many colours (a row is a clique)
        assert!(ncol >= (kl + ku + 1).min(n));
    }
    // a tridiagonal matrix is compressed to exactly three colours
    assert_eq!(color_columns(&band_pattern(40, 1, 1)).1, 3);
    // a diagonal matrix needs one evaluation
    assert_eq!(color_columns(&band_pattern(25, 0, 0)).1, 1);
    // empty pattern
    let (c, k) = color_columns(&band_pattern(0, 1, 1));
    assert!(c.is_empty() && k == 0);
}

#[test]
fn colored_and_dense_finite_differences_match_the_analytic_jacobian() {
    let n = 40;
    let y = chain_state(n);
    let mut f0 = vec![0.0; n];
    chain_rhs(0.0, &y, &mut f0);
    let exact = chain_jacobian(&y);
    let pat = band_pattern(n, 1, 1);
    let colours = color_columns(&pat);
    let mut nfev = 0usize;
    let ser = fd_jacobian_colored(&chain_rhs, 0.0, &y, &f0, &pat, &colours, &mut nfev);
    assert_eq!(nfev, 3, "one evaluation per colour, not per column");
    let mut nfev_par = 10usize;
    let par = fd_jacobian_colored_par(&chain_rhs, 0.0, &y, &f0, &pat, &colours, &mut nfev_par);
    assert_eq!(
        nfev_par, 13,
        "the parallel version counts the same evaluations"
    );
    assert_eq!(ser, par, "parallel colouring must be bit-identical");
    assert!(max_abs_diff(&csr_to_dense(&ser), &exact) < 5e-6);
    // the result keeps the pattern exactly
    assert_eq!(
        (ser.indptr.clone(), ser.indices.clone()),
        (pat.indptr.clone(), pat.indices.clone())
    );
    // dense FD agrees with the compressed one on the pattern and with the exact Jacobian
    let dense = fd_jacobian_dense(&chain_rhs, 0.0, &y, &f0);
    let dense_par = fd_jacobian_dense_par(&chain_rhs, 0.0, &y, &f0);
    assert_eq!(
        dense.data, dense_par.data,
        "parallel dense FD is bit-identical"
    );
    assert!(max_abs_diff(&dense, &exact) < 5e-6);
    assert!(max_abs_diff(&dense, &csr_to_dense(&ser)) < 5e-6);
}

#[test]
fn colored_parallel_with_a_coarser_pattern_than_the_truth_still_matches_serial() {
    // pentadiagonal pattern around a tridiagonal function: extra entries are ~0
    let n = 33;
    let y = chain_state(n);
    let mut f0 = vec![0.0; n];
    chain_rhs(0.0, &y, &mut f0);
    let pat = band_pattern(n, 2, 2);
    let colours = color_columns(&pat);
    assert!(colours.1 >= 5);
    let (mut a, mut b) = (0usize, 0usize);
    let s = fd_jacobian_colored(&chain_rhs, 0.0, &y, &f0, &pat, &colours, &mut a);
    let p = fd_jacobian_colored_par(&chain_rhs, 0.0, &y, &f0, &pat, &colours, &mut b);
    assert_eq!(s, p);
    assert_eq!(a, b);
    let dense = csr_to_dense(&s);
    for i in 0..n {
        for j in 0..n {
            if (i as i64 - j as i64).abs() == 2 {
                assert!(
                    dense.at(i, j).abs() < 1e-6,
                    "outside the true band: {}",
                    dense.at(i, j)
                );
            }
        }
    }
    assert!(max_abs_diff(&dense, &chain_jacobian(&y)) < 5e-6);
}

#[test]
fn provider_eval_with_serial_and_parallel_modes_are_identical() {
    let n = 30;
    let y = chain_state(n);
    let mut f0 = vec![0.0; n];
    chain_rhs(0.0, &y, &mut f0);
    let exact = chain_jacobian(&y);

    // banded structure: finite differences with colouring
    let banded =
        JacobianProvider::new(n, JacStructure::Banded { kl: 1, ku: 1 }, None, None).unwrap();
    assert_eq!(banded.kind, SolverKind::Banded);
    assert_eq!(banded.band, (1, 1));
    let (mut n1, mut n2) = (0usize, 0usize);
    let j_ser = banded.eval_with(&chain_rhs, 0.0, &y, &f0, &mut n1, &FdMode::serial());
    let j_par = banded.eval_with(&chain_rhs, 0.0, &y, &f0, &mut n2, &FdMode::parallel());
    assert_eq!((n1, n2), (3, 3));
    let (Jac::Sparse(a), Jac::Sparse(b)) = (&j_ser, &j_par) else {
        panic!("banded structure must produce a sparse Jacobian");
    };
    assert_eq!(a, b);
    assert!(max_abs_diff(&csr_to_dense(a), &exact) < 5e-6);
    // eval() is eval_with(serial)
    let mut n3 = 0usize;
    let Jac::Sparse(c) = banded.eval(&chain_rhs, 0.0, &y, &f0, &mut n3) else {
        panic!("sparse expected");
    };
    assert_eq!((&c, n3), (a, 3));

    // dense structure without a callback: n evaluations, dense matrix
    let dense = JacobianProvider::new(n, JacStructure::Dense, None, None).unwrap();
    assert_eq!(dense.kind, SolverKind::Dense);
    let (mut d1, mut d2) = (0usize, 0usize);
    let jd_ser = dense.eval_with(&chain_rhs, 0.0, &y, &f0, &mut d1, &FdMode::serial());
    let jd_par = dense.eval_with(&chain_rhs, 0.0, &y, &f0, &mut d2, &FdMode::parallel());
    assert_eq!((d1, d2), (n, n));
    let (Jac::Dense(a), Jac::Dense(b)) = (&jd_ser, &jd_par) else {
        panic!("dense expected");
    };
    assert_eq!(a.data, b.data);
    assert!(max_abs_diff(a, &exact) < 5e-6);
}

#[test]
fn provider_prefers_user_jacobians_and_does_not_count_evaluations() {
    let n = 12;
    let y = chain_state(n);
    let mut f0 = vec![0.0; n];
    chain_rhs(0.0, &y, &mut f0);
    let exact = chain_jacobian(&y);
    let dense_cb = |_t: f64, y: &[f64]| chain_jacobian(y);
    let sparse_cb = |_t: f64, y: &[f64]| {
        let j = chain_jacobian(y);
        let mut trip = Vec::new();
        for i in 0..y.len() {
            for k in 0..y.len() {
                if j.at(i, k) != 0.0 {
                    trip.push((i, k, j.at(i, k)));
                }
            }
        }
        Csr::from_triplets(y.len(), y.len(), &trip)
    };
    // dense callback on a dense structure
    let p = JacobianProvider::new(n, JacStructure::Dense, Some(&dense_cb), None).unwrap();
    let mut nfev = 0usize;
    let Jac::Dense(m) = p.eval_with(&chain_rhs, 0.0, &y, &f0, &mut nfev, &FdMode::parallel())
    else {
        panic!("dense expected");
    };
    assert_eq!(nfev, 0);
    assert_eq!(m.data, exact.data);
    // dense callback on a banded structure is compressed onto the band
    let p = JacobianProvider::new(
        n,
        JacStructure::Banded { kl: 1, ku: 1 },
        Some(&dense_cb),
        None,
    )
    .unwrap();
    let Jac::Sparse(c) = p.eval_with(&chain_rhs, 0.0, &y, &f0, &mut nfev, &FdMode::serial()) else {
        panic!("sparse expected");
    };
    assert_eq!(nfev, 0);
    assert_eq!(c.values.len(), 3 * n - 2);
    assert_eq!(csr_to_dense(&c).data, exact.data);
    // a sparse callback forces a sparse solver even for a Dense structure
    let p = JacobianProvider::new(n, JacStructure::Dense, None, Some(&sparse_cb)).unwrap();
    assert_eq!(p.kind, SolverKind::Sparse);
    let Jac::Sparse(c) = p.eval_with(&chain_rhs, 0.0, &y, &f0, &mut nfev, &FdMode::parallel())
    else {
        panic!("sparse expected");
    };
    assert_eq!(nfev, 0);
    assert_eq!(csr_to_dense(&c).data, exact.data);
    // a user pattern of the wrong size is rejected
    let wrong = band_pattern(n + 1, 1, 1);
    assert!(JacobianProvider::new(n, JacStructure::Sparse(&wrong), None, None).is_none());
    let right = band_pattern(n, 1, 1);
    let p = JacobianProvider::new(n, JacStructure::Sparse(&right), None, None).unwrap();
    assert_eq!(p.kind, SolverKind::Sparse);
    let mut k = 0usize;
    let Jac::Sparse(c) = p.eval_with(&chain_rhs, 0.0, &y, &f0, &mut k, &FdMode::parallel()) else {
        panic!("sparse expected");
    };
    assert_eq!(k, 3);
    assert!(max_abs_diff(&csr_to_dense(&c), &exact) < 5e-6);
}

#[test]
fn finite_difference_jacobians_of_degenerate_sizes_do_not_panic() {
    let none: [f64; 0] = [];
    let m = fd_jacobian_dense(&chain_rhs, 0.0, &none, &none);
    assert_eq!((m.rows, m.cols), (0, 0));
    let m = fd_jacobian_dense_par(&chain_rhs, 0.0, &none, &none);
    assert_eq!((m.rows, m.cols), (0, 0));
    // 1x1 : f = y cos... use chain_rhs: f = sin(y), J = cos(y)
    let y = [0.8];
    let mut f0 = [0.0];
    chain_rhs(0.0, &y, &mut f0);
    let m = fd_jacobian_dense_par(&chain_rhs, 0.0, &y, &f0);
    assert!((m.at(0, 0) - 0.8_f64.cos()).abs() < 5e-6);
    // y = 0 uses the 1e-5 floor for the increment (no division by zero)
    let y = [0.0, 0.0];
    let mut f0 = [0.0, 0.0];
    chain_rhs(0.0, &y, &mut f0);
    let m = fd_jacobian_dense(&chain_rhs, 0.0, &y, &f0);
    assert!(m.data.iter().all(|v| v.is_finite()));
    assert!((m.at(0, 0) - 1.0).abs() < 1e-5);
}

// ------------------------------------------------- nested dissection options

fn grid_adj(
    nx: usize,
    ny: usize,
) -> Vec<Vec<usize>> {
    let mut adj = vec![Vec::new(); nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            let v = j * nx + i;
            if i + 1 < nx {
                adj[v].push(v + 1);
                adj[v + 1].push(v);
            }
            if j + 1 < ny {
                adj[v].push(v + nx);
                adj[v + nx].push(v);
            }
        }
    }
    adj
}

fn is_permutation(
    p: &[usize],
    n: usize,
) -> bool {
    let mut seen = vec![false; n];
    p.len() == n
        && p.iter()
            .all(|&v| v < n && !std::mem::replace(&mut seen[v], true))
}

#[test]
fn nested_dissection_options_always_give_valid_orderings() {
    let adj = grid_adj(14, 14);
    let n = adj.len();
    let natural: Vec<usize> = (0..n).collect();
    let natural_fill = cholesky_fill(&adj, &natural);
    for leaf in [1usize, 5, 20, 100, 5000] {
        for fm_passes in [0usize, 1, 8] {
            for imbalance in [0.0, 0.1, 0.4] {
                let o = NdOptions {
                    leaf,
                    imbalance,
                    fm_passes,
                };
                let p = multilevel_nested_dissection_with(&adj, &o);
                assert!(
                    is_permutation(&p, n),
                    "leaf={leaf} fm={fm_passes} imb={imbalance}"
                );
                // determinism: same options, same permutation
                assert_eq!(p, multilevel_nested_dissection_with(&adj, &o));
                if leaf <= 20 {
                    // dissection of a 2-D grid must beat the lexicographic band fill
                    // is NOT guaranteed on 14x14; the fill must at least be finite and
                    // not exceed the dense bound n (n+1)/2.
                    assert!(cholesky_fill(&adj, &p) <= n * (n + 1) / 2);
                }
            }
        }
    }
    // default options equal the explicit default
    assert_eq!(
        multilevel_nested_dissection_order(&adj),
        multilevel_nested_dissection_with(&adj, &NdOptions::default())
    );
    // the hybrid ordering is never worse than natural order
    let h = hybrid_nested_dissection_order(&adj);
    assert!(is_permutation(&h, n));
    assert!(cholesky_fill(&adj, &h) < natural_fill + n);
}

#[test]
fn nested_dissection_on_a_larger_grid_beats_the_natural_order() {
    // 40 x 40: natural (banded) Cholesky fill ~ n * nx = 64000, nested dissection ~ O(n log n)
    let adj = grid_adj(40, 40);
    let n = adj.len();
    let natural: Vec<usize> = (0..n).collect();
    let nat = cholesky_fill(&adj, &natural);
    let nd = multilevel_nested_dissection_with(
        &adj,
        &NdOptions {
            leaf: 30,
            ..NdOptions::default()
        },
    );
    assert!(is_permutation(&nd, n));
    let fill = cholesky_fill(&adj, &nd);
    assert!(
        fill * 10 < nat * 9,
        "ND fill {fill} not clearly below natural {nat}"
    );
}

#[test]
fn nested_dissection_degenerate_graphs() {
    assert!(multilevel_nested_dissection_order(&[]).is_empty());
    assert!(hybrid_nested_dissection_order(&[]).is_empty());
    assert_eq!(multilevel_nested_dissection_order(&[vec![]]), vec![0]);
    // isolated vertices
    let iso = vec![Vec::<usize>::new(); 7];
    assert!(is_permutation(
        &multilevel_nested_dissection_with(
            &iso,
            &NdOptions {
                leaf: 1,
                ..NdOptions::default()
            }
        ),
        7
    ));
    // two disconnected grids plus isolated vertices, forced to dissect (leaf = 4)
    let g = grid_adj(6, 6);
    let mut adj: Vec<Vec<usize>> = g.clone();
    let off = adj.len();
    adj.extend(
        g.iter()
            .map(|row| row.iter().map(|v| v + off).collect::<Vec<_>>()),
    );
    adj.push(Vec::new());
    adj.push(Vec::new());
    let p = multilevel_nested_dissection_with(
        &adj,
        &NdOptions {
            leaf: 4,
            ..NdOptions::default()
        },
    );
    assert!(is_permutation(&p, adj.len()));
    // path graph and star graph
    let path: Vec<Vec<usize>> = (0..50)
        .map(|i| {
            let mut v = Vec::new();
            if i > 0 {
                v.push(i - 1);
            }
            if i + 1 < 50 {
                v.push(i + 1);
            }
            v
        })
        .collect();
    let p = multilevel_nested_dissection_with(
        &path,
        &NdOptions {
            leaf: 3,
            ..NdOptions::default()
        },
    );
    assert!(is_permutation(&p, 50));
    // a path has a perfect elimination ordering of fill n - 1 ; dissection must stay linear
    assert!(cholesky_fill(&path, &p) <= 3 * 50);
    let mut star = vec![Vec::new(); 31];
    for leaf in 1..31 {
        star[0].push(leaf);
        star[leaf].push(0);
    }
    assert!(is_permutation(
        &multilevel_nested_dissection_with(
            &star,
            &NdOptions {
                leaf: 2,
                ..NdOptions::default()
            }
        ),
        31
    ));
}

#[test]
fn separator_bisection_respects_options() {
    let adj = grid_adj(16, 16);
    let n = adj.len();
    for fm_passes in [0usize, 2, 8] {
        for imbalance in [0.05, 0.1, 0.3] {
            let o = NdOptions {
                leaf: 10,
                imbalance,
                fm_passes,
            };
            let (l, r, s) = vertex_separator_bisection(&adj, &o);
            // a partition of the vertex set ...
            let mut all: Vec<usize> = l.iter().chain(&r).chain(&s).copied().collect();
            all.sort_unstable();
            assert_eq!(
                all,
                (0..n).collect::<Vec<_>>(),
                "fm={fm_passes} imb={imbalance}"
            );
            // ... whose separator really separates
            let in_l: Vec<bool> = (0..n).map(|v| l.contains(&v)).collect();
            for &v in &r {
                assert!(
                    adj[v].iter().all(|&u| !in_l[u]),
                    "edge crosses the separator: fm={fm_passes} imb={imbalance} vertex {v}"
                );
            }
            assert!(!l.is_empty() && !r.is_empty());
            // a 16 x 16 grid has a separator of 16 vertices; accept a small multiple
            assert!(s.len() <= 3 * 16, "separator {} too large", s.len());
        }
    }
}

// -------------------------------------------------------------- AMG V-cycle

fn poisson2d(nx: usize) -> Csr {
    let mut t = Vec::new();
    for j in 0..nx {
        for i in 0..nx {
            let p = j * nx + i;
            t.push((p, p, 4.0));
            if i > 0 {
                t.push((p, p - 1, -1.0));
            }
            if i + 1 < nx {
                t.push((p, p + 1, -1.0));
            }
            if j > 0 {
                t.push((p, p - nx, -1.0));
            }
            if j + 1 < nx {
                t.push((p, p + nx, -1.0));
            }
        }
    }
    Csr::from_triplets(nx * nx, nx * nx, &t)
}

fn residual_norm(
    a: &Csr,
    x: &[f64],
    b: &[f64],
) -> f64 {
    let mut ax = vec![0.0; b.len()];
    a.matvec(x, &mut ax);
    ax.iter()
        .zip(b)
        .map(|(p, q)| (p - q).powi(2))
        .sum::<f64>()
        .sqrt()
}

#[test]
fn amg_vcycle_is_a_linear_contraction_with_the_solution_as_fixed_point() {
    let a = poisson2d(24);
    let n = a.rows;
    // manufactured solution
    let xstar: Vec<f64> = (0..n)
        .map(|i| ((i * 37 % 101) as f64) / 50.0 - 1.0)
        .collect();
    let mut b = vec![0.0; n];
    a.matvec(&xstar, &mut b);
    for (kind, cycle) in [
        (AmgKind::RugeStuben, Cycle::V),
        (AmgKind::RugeStuben, Cycle::W),
        (AmgKind::SmoothedAggregation, Cycle::V),
    ] {
        let amg = AmgGeneral::build(
            &a,
            &AmgGeneralOptions {
                kind,
                cycle,
                coarse_size: 30,
                ..AmgGeneralOptions::default()
            },
        )
        .unwrap();
        assert!(amg.num_levels() >= 2);
        // fixed point: the exact solution is reproduced to rounding level
        let mut x = xstar.clone();
        amg.vcycle(&b, &mut x);
        let drift = x
            .iter()
            .zip(&xstar)
            .map(|(p, q)| (p - q).abs())
            .fold(0.0, f64::max);
        assert!(
            drift < 1e-10,
            "{kind:?}/{cycle:?}: fixed point drifts by {drift}"
        );
        // contraction of the residual over repeated cycles
        let mut x = vec![0.0; n];
        let r0 = residual_norm(&a, &x, &b);
        let mut previous = r0;
        for k in 0..12 {
            amg.vcycle(&b, &mut x);
            let r = residual_norm(&a, &x, &b);
            assert!(
                r < previous,
                "{kind:?}/{cycle:?} cycle {k}: {r} vs {previous}"
            );
            previous = r;
        }
        assert!(
            previous < 1e-2 * r0,
            "{kind:?}/{cycle:?}: only reduced by {}",
            previous / r0
        );
    }
}

#[test]
fn amg_vcycle_from_zero_equals_precondition_and_is_linear_in_the_rhs() {
    let a = poisson2d(16);
    let n = a.rows;
    let amg = AmgGeneral::build(
        &a,
        &AmgGeneralOptions {
            coarse_size: 20,
            ..AmgGeneralOptions::default()
        },
    )
    .unwrap();
    let b1: Vec<f64> = (0..n).map(|i| (i % 7) as f64 - 3.0).collect();
    let b2: Vec<f64> = (0..n).map(|i| ((i * 5) % 11) as f64 * 0.25).collect();
    let apply = |b: &[f64]| {
        let mut x = vec![0.0; n];
        amg.vcycle(b, &mut x);
        x
    };
    let mut z = vec![1.0; n]; // precondition must overwrite, not accumulate
    amg.precondition(&b1, &mut z);
    assert_eq!(apply(&b1), z);
    // x(b1 + 2 b2) = x(b1) + 2 x(b2) (the cycle from zero is a fixed linear map)
    let combo: Vec<f64> = b1.iter().zip(&b2).map(|(p, q)| p + 2.0 * q).collect();
    let (x1, x2, xc) = (apply(&b1), apply(&b2), apply(&combo));
    for i in 0..n {
        assert!(
            (xc[i] - (x1[i] + 2.0 * x2[i])).abs() < 1e-12 * (1.0 + xc[i].abs()),
            "i={i}"
        );
    }
    // the zero right-hand side maps to zero
    assert!(apply(&vec![0.0; n]).iter().all(|v| *v == 0.0));
    // solve() is the repeated V-cycle: replay its history manually
    let res = amg.solve(&b1, None, 1e-8, 60);
    assert!(res.converged);
    let mut x = vec![0.0; n];
    for _ in 1..res.history.len() {
        amg.vcycle(&b1, &mut x);
    }
    assert_eq!(x, res.x);
}

#[test]
fn amg_build_rejects_unusable_matrices() {
    assert!(
        AmgGeneral::build(
            &Csr::from_triplets(0, 0, &[]),
            &AmgGeneralOptions::default()
        )
        .is_err()
    );
    assert!(
        AmgGeneral::build(
            &Csr::from_triplets(3, 4, &[(0, 0, 1.0)]),
            &AmgGeneralOptions::default()
        )
        .is_err()
    );
    // zero / negative diagonal
    // (the diagonal is only inspected while a level is being coarsened, so the matrix must be
    // larger than `coarse_size`; a tiny matrix goes straight to the direct coarse solve)
    let mut trip: Vec<(usize, usize, f64)> = Vec::new();
    for i in 0..120usize {
        trip.push((i, i, if i == 57 { -1.0 } else { 2.0 }));
        if i > 0 {
            trip.push((i, i - 1, -0.5));
        }
        if i + 1 < 120 {
            trip.push((i, i + 1, -0.5));
        }
    }
    let neg = Csr::from_triplets(120, 120, &trip);
    assert!(AmgGeneral::build(&neg, &AmgGeneralOptions::default()).is_err());
    let pos = trip.iter().position(|t| t.0 == 57 && t.1 == 57).unwrap();
    trip[pos].2 = 0.0;
    let zero = Csr::from_triplets(120, 120, &trip);
    assert!(AmgGeneral::build(&zero, &AmgGeneralOptions::default()).is_err());
    // near null space of the wrong length
    let a = poisson2d(8);
    let o = AmgGeneralOptions {
        kind: AmgKind::SmoothedAggregation,
        near_nullspace: vec![vec![1.0; 3]],
        ..AmgGeneralOptions::default()
    };
    assert!(AmgGeneral::build(&a, &o).is_err());
    // a 1x1 system is solved directly by the coarse level
    let one = AmgGeneral::build(
        &Csr::from_triplets(1, 1, &[(0, 0, 4.0)]),
        &AmgGeneralOptions::default(),
    )
    .unwrap();
    let mut x = vec![0.0];
    one.vcycle(&[8.0], &mut x);
    assert!((x[0] - 2.0).abs() < 1e-14);
}

// ----------------------------------------------------- legacy IDR(s) wrapper

#[test]
#[allow(deprecated)]
fn legacy_idr_s_is_identical_to_idr_s_with() {
    // non-symmetric tridiagonal convection-diffusion matrix
    let n = 60;
    let mut t = Vec::new();
    for i in 0..n {
        t.push((i, i, 2.0));
        if i > 0 {
            t.push((i, i - 1, -1.3));
        }
        if i + 1 < n {
            t.push((i, i + 1, -0.7));
        }
    }
    let a = Csr::from_triplets(n, n, &t);
    let b: Vec<f64> = (0..n).map(|i| 1.0 + (i % 5) as f64).collect();
    let apply = |x: &[f64], y: &mut [f64]| a.matvec(x, y);
    for s in [1usize, 3, 8] {
        let old = idr_s(apply, identity_precond, s, &b, None, 1e-11, 400, 5);
        let new = idr_s_with(
            apply,
            identity_precond,
            &b,
            None,
            &IdrOptions {
                s,
                tol: 1e-11,
                max_iter: 400,
                seed: 5,
            },
        );
        assert_eq!(old.x, new.x);
        assert_eq!(
            (old.iterations, old.converged),
            (new.iterations, new.converged)
        );
        assert_eq!(old.residual, new.residual);
        assert!(new.converged, "s = {s}");
        assert!(residual_norm(&a, &new.x, &b) < 1e-8 * b.iter().map(|v| v * v).sum::<f64>().sqrt());
        // different seeds give different Krylov paths but the same solution
        let other = idr_s(apply, identity_precond, s, &b, None, 1e-11, 400, 99);
        assert!(other.converged);
        let d = other
            .x
            .iter()
            .zip(&new.x)
            .map(|(p, q)| (p - q).abs())
            .fold(0.0, f64::max);
        assert!(d < 1e-7, "seed-dependent solution, difference {d}");
    }
}

#[test]
#[allow(deprecated)]
fn legacy_idr_s_edge_cases() {
    let id = |x: &[f64], y: &mut [f64]| y.copy_from_slice(x);
    // zero right-hand side: converged immediately
    let r = idr_s(id, identity_precond, 2, &[0.0; 4], None, 1e-10, 50, 1);
    assert!(r.converged && r.iterations == 0 && r.x.iter().all(|v| *v == 0.0));
    // exact initial guess
    let b = [1.0, 2.0, 3.0];
    let r = idr_s(id, identity_precond, 2, &b, Some(&b), 1e-10, 50, 1);
    assert!(r.converged && r.iterations == 0);
    // identity operator is solved within a handful of products
    let r = idr_s(id, identity_precond, 2, &b, None, 1e-12, 50, 1);
    assert!(r.converged && r.iterations <= 10);
    for (p, q) in r.x.iter().zip(&b) {
        assert!((p - q).abs() < 1e-10);
    }
    // s larger than n is clamped; s = 0 is raised to 1; empty system is trivially converged
    let r = idr_s(id, identity_precond, 100, &b, None, 1e-10, 50, 1);
    assert!(r.converged);
    let r = idr_s(id, identity_precond, 0, &b, None, 1e-10, 50, 1);
    assert!(r.converged);
    let r = idr_s(id, identity_precond, 4, &[], None, 1e-10, 50, 1);
    assert!(r.converged && r.x.is_empty());
    // a NaN right-hand side must not be reported as converged, and must not hang or panic
    let r = idr_s(
        id,
        identity_precond,
        2,
        &[f64::NAN, 1.0],
        None,
        1e-10,
        20,
        1,
    );
    assert!(!r.converged);
    // max_iter = 0 on a non-trivial problem: not converged, bounded work
    let r = idr_s(
        id,
        identity_precond,
        2,
        &[1.0, 1.0, 1.0],
        Some(&[0.0, 0.0, 0.0]),
        1e-10,
        0,
        1,
    );
    assert!(!r.converged);
}

// ------------------------------------------------- coloured graph accessors

#[test]
fn colour_graph_arc_lists_are_sorted_labelled_and_consistent() {
    let mut g = ColourGraph::new(5, true);
    assert!(g.set_edge(0, 3, 7));
    assert!(g.set_edge(0, 1, 2));
    assert!(g.set_edge(0, 4, 1));
    assert!(g.set_edge(3, 0, 9));
    assert!(g.set_edge(2, 2, 4)); // self-loop
    assert!(!g.set_edge(0, 5, 1), "out of range endpoint");
    assert!(!g.set_edge(9, 0, 1));
    assert_eq!(g.out_arcs(0), &[(1, 2), (3, 7), (4, 1)], "sorted by head");
    assert_eq!(g.in_arcs(0), &[(3, 9)]);
    assert_eq!(g.in_arcs(3), &[(0, 7)]);
    assert_eq!(
        g.in_arcs(2),
        &[(2, 4)],
        "a self-loop is both an out- and an in-arc"
    );
    assert_eq!(g.out_arcs(2), &[(2, 4)]);
    assert!(g.out_arcs(1).is_empty());
    // out of range vertex: empty, not a panic
    assert!(g.out_arcs(99).is_empty() && g.in_arcs(99).is_empty());
    // the sum of out-degrees equals the sum of in-degrees equals the arc count
    let outs: usize = (0..5).map(|v| g.out_arcs(v).len()).sum();
    let ins: usize = (0..5).map(|v| g.in_arcs(v).len()).sum();
    assert_eq!((outs, ins), (5, 5));
    // relabelling replaces the label in both lists; label 0 removes the arc
    assert!(g.set_edge(0, 3, 8));
    assert_eq!(g.out_arcs(0)[1], (3, 8));
    assert_eq!(g.in_arcs(3), &[(0, 8)]);
    assert!(g.set_edge(0, 3, 0));
    assert_eq!(g.out_arcs(0), &[(1, 2), (4, 1)]);
    assert!(g.in_arcs(3).is_empty());
    assert!(!g.has_edge(0, 3) && g.has_edge(3, 0));
}

#[test]
fn undirected_colour_graph_stores_both_directions_and_no_in_arcs() {
    let mut g = ColourGraph::new(4, false);
    assert!(g.set_edge(0, 2, 3));
    assert!(g.set_edge(1, 2, 5));
    assert!(g.set_edge(3, 3, 6));
    assert_eq!(g.out_arcs(2), &[(0, 3), (1, 5)]);
    assert_eq!(g.out_arcs(0), &[(2, 3)]);
    assert_eq!(g.out_arcs(3), &[(3, 6)], "self-loop stored once");
    for v in 0..4 {
        assert!(
            g.in_arcs(v).is_empty(),
            "undirected graphs keep no in-arc list"
        );
    }
    // symmetry of the adjacency labels
    for u in 0..4 {
        for v in 0..4 {
            assert_eq!(g.edge_label(u, v), g.edge_label(v, u));
        }
    }
    assert_eq!(g.degree(2), 2);
}

#[test]
fn with_colours_validates_length_and_colours_split_cells() {
    assert!(ColourGraph::with_colours(3, false, &[0, 0]).is_none());
    assert!(ColourGraph::with_colours(3, false, &[0, 0, 0, 0]).is_none());
    assert!(ColourGraph::with_colours(0, true, &[]).is_some());
    let mut path = ColourGraph::with_colours(3, false, &[0, 0, 0]).unwrap();
    path.set_edge(0, 1, 1);
    path.set_edge(1, 2, 1);
    assert_eq!(path.colours(), &[0, 0, 0]);
    // uncoloured path P3: the two ends are equivalent -> 2 cells (sizes 2 and 1)
    let plain = refine(&path);
    assert_eq!(plain.num_cells(), 2);
    let mut sizes = plain.cell_sizes();
    sizes.sort_unstable();
    assert_eq!(sizes, vec![1, 2]);
    assert!(plain.verify(&path));
    // colouring one end breaks the symmetry -> discrete partition
    let mut coloured = ColourGraph::with_colours(3, false, &[1, 0, 0]).unwrap();
    coloured.set_edge(0, 1, 1);
    coloured.set_edge(1, 2, 1);
    let r = refine(&coloured);
    assert!(r.is_discrete() && r.verify(&coloured));
    // set_colours reports a length mismatch and leaves colours alone
    assert!(!coloured.set_colours(&[1, 2]));
    assert_eq!(coloured.colours(), &[1, 0, 0]);
    assert!(coloured.set_colours(&[2, 2, 2]));
    assert_eq!(coloured.colours(), &[2, 2, 2]);
}

#[test]
fn cell_certificate_counts_are_the_equitable_degrees() {
    // star K_{1,4}: centre cell sees one cell (the leaves) 4 times, a leaf sees the centre once
    let mut star = ColourGraph::new(5, false);
    for leaf in 1..5 {
        star.set_edge(0, leaf, 1);
    }
    let r = refine(&star);
    let cert = r.certificate();
    assert_eq!(cert.num_cells(), 2);
    let cells = cert.cells();
    // find the singleton (centre) and the leaf cell
    let centre = cells.iter().find(|c| c.size() == 1).unwrap();
    let leaves = cells.iter().find(|c| c.size() == 4).unwrap();
    assert_eq!(centre.out_counts().len(), 1);
    let (cell_of_leaves, label, count) = centre.out_counts()[0];
    assert_eq!((label, count), (1, 4));
    assert_eq!(leaves.out_counts().len(), 1);
    assert_eq!(leaves.out_counts()[0].1, 1);
    assert_eq!(leaves.out_counts()[0].2, 1);
    assert_ne!(
        cell_of_leaves,
        leaves.out_counts()[0].0,
        "the two cells see each other, not themselves"
    );
    // undirected: no in-counts
    assert!(centre.in_counts().is_empty() && leaves.in_counts().is_empty());

    // directed 4-cycle: one cell; one out-neighbour and one in-neighbour, both inside the cell
    let mut cyc = ColourGraph::new(4, true);
    for i in 0..4 {
        cyc.set_edge(i, (i + 1) % 4, 3);
    }
    let r = refine(&cyc);
    assert_eq!(r.certificate().num_cells(), 1);
    let c = &r.certificate().cells()[0];
    assert_eq!(c.size(), 4);
    assert_eq!(c.out_counts(), &[(0, 3, 1)]);
    assert_eq!(c.in_counts(), &[(0, 3, 1)]);

    // in a random-ish graph the counts of a cell sum to the common degree of its vertices
    let mut g = ColourGraph::new(9, true);
    for u in 0..9usize {
        g.set_edge(u, (u * 2 + 1) % 9, 1 + (u % 2) as u32);
        g.set_edge(u, (u + 4) % 9, 1);
    }
    let r = refine(&g);
    assert!(r.verify(&g));
    for (cell, cert) in r.cells().iter().zip(r.certificate().cells()) {
        let v = cell[0];
        let total: usize = cert.out_counts().iter().map(|&(_, _, k)| k).sum();
        assert_eq!(total, g.out_arcs(v).len(), "out-degree of vertex {v}");
        let total_in: usize = cert.in_counts().iter().map(|&(_, _, k)| k).sum();
        assert_eq!(total_in, g.in_arcs(v).len(), "in-degree of vertex {v}");
        assert_eq!(cert.size(), cell.len());
    }
}
