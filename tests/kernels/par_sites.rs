//! Determinism and timing checks of the in-algorithm parallel sites
//! (Part 6 of `docs/ALGORITHM_COMPLETION_CHECKLIST.md`).
//!
//! Every site is run once inside a one-thread rayon pool (the serial
//! behaviour) and once in the global pool; the two results must be equal
//! (`Debug` text, so floats are compared bit-exactly). The timings are
//! printed, never asserted: run with `--nocapture` and grep `PAR`.

use std::time::Instant;

use rssn::kernels::homotopy::HomotopyOptions;
use rssn::kernels::homotopy::Precision;
use rssn::kernels::homotopy::homotopy_solve;
use rssn::kernels::homotopy_cert::certify_endpoints;
use rssn::kernels::homotopy_cert::certify_endpoints_rigorous;
use rssn::kernels::homotopy_cert::homotopy_solve_verified;
use rssn::kernels::rootfind::PolyTerm;

/// Runs `f` serially (one-thread pool) and in the global pool, asserts the
/// results agree exactly and prints one timing row.
pub fn compare<T: std::fmt::Debug + Send>(
    site: &str,
    size: &str,
    f: impl Fn() -> T + Send + Sync,
) {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let t0 = Instant::now();
    let a = pool.install(&f);
    let ts = t0.elapsed().as_secs_f64();
    let t1 = Instant::now();
    let b = f();
    let tp = t1.elapsed().as_secs_f64();
    assert_eq!(
        format!("{a:?}"),
        format!("{b:?}"),
        "{site}: parallel result differs from serial"
    );
    println!(
        "PAR | {site} | {size} | serial {:.1} ms | parallel {:.1} ms | x{:.2} | threads {}",
        ts * 1e3,
        tp * 1e3,
        ts / tp.max(1e-12),
        rayon::current_num_threads()
    );
}

fn t(
    coef: f64,
    exps: &[u32],
) -> PolyTerm {
    PolyTerm {
        coef,
        exps: exps.to_vec(),
    }
}

/// Cyclic n-roots system.
fn cyclic(n: usize) -> Vec<Vec<PolyTerm>> {
    let mut sys = Vec::new();
    for k in 1..n {
        let mut eq = Vec::new();
        for i in 0..n {
            let mut e = vec![0u32; n];
            for j in 0..k {
                e[(i + j) % n] = 1;
            }
            eq.push(t(1.0, &e));
        }
        sys.push(eq);
    }
    let mut e = vec![1u32; n];
    eq_last(&mut sys, &mut e);
    sys
}

fn eq_last(
    sys: &mut Vec<Vec<PolyTerm>>,
    e: &mut [u32],
) {
    sys.push(vec![t(1.0, e), t(-1.0, &vec![0u32; e.len()])]);
}

#[test]
fn par_homotopy_paths_and_certification() {
    let sys = cyclic(5); // Bezout number 120
    let opts = HomotopyOptions::default();
    compare("P1 homotopy_solve (double)", "cyclic-5, 120 paths", || {
        homotopy_solve(&sys, &opts).unwrap()
    });
    let sys4 = cyclic(3);
    let adaptive = HomotopyOptions {
        precision: Precision::Adaptive,
        ..HomotopyOptions::default()
    };
    compare("P1 homotopy_solve (adaptive)", "cyclic-3, 6 paths", || {
        homotopy_solve(&sys4, &adaptive).unwrap()
    });
    let ends = homotopy_solve(&sys, &opts).unwrap();
    compare("P2 certify_endpoints", "120 endpoints", || {
        certify_endpoints(&sys, &ends, 1.0)
    });
    compare("P2 certify_endpoints_rigorous", "120 endpoints", || {
        certify_endpoints_rigorous(&sys, &ends)
    });
    compare("P1+P2 homotopy_solve_verified", "cyclic-3, 6 paths", || {
        homotopy_solve_verified(&sys4, &opts).unwrap()
    });
}

// ---------------------------------------------------------------- P3: Groebner

mod groebner_sites {
    use super::compare;
    use num_rational::BigRational;
    use rssn::rules::poly::groebner::GPoly;
    use rssn::rules::poly::groebner::GroebnerLimits;
    use rssn::rules::poly::groebner::Order;
    use rssn::rules::poly::groebner::fglm;
    use rssn::rules::poly::groebner::groebner;

    fn q(n: i64) -> BigRational {
        BigRational::from_integer(n.into())
    }
    fn mono(
        n: usize,
        powers: &[(usize, u32)],
    ) -> Vec<u32> {
        let mut m = vec![0; n];
        for &(i, e) in powers {
            m[i] += e;
        }
        m
    }
    fn poly(
        terms: &[(i64, Vec<u32>)],
        order: Order,
    ) -> GPoly {
        GPoly::new(
            terms.iter().map(|(c, e)| (e.clone(), q(*c))).collect(),
            order,
        )
    }
    pub fn katsura(
        m: usize,
        order: Order,
    ) -> Vec<GPoly> {
        let n = m + 1;
        let mut out = Vec::new();
        let mut linear = vec![(1, mono(n, &[(0, 1)])), (-1, mono(n, &[]))];
        for l in 1..=m {
            linear.push((2, mono(n, &[(l, 1)])));
        }
        out.push(poly(&linear, order));
        for k in 0..m {
            let mut terms = vec![(-1, mono(n, &[(k, 1)]))];
            for l in -(m as i64)..=(m as i64) {
                let (a, b) = (
                    l.unsigned_abs() as usize,
                    (k as i64 - l).unsigned_abs() as usize,
                );
                if a <= m && b <= m {
                    terms.push((1, mono(n, &[(a, 1), (b, 1)])));
                }
            }
            out.push(poly(&terms, order));
        }
        out
    }
    pub fn cyclic(
        n: usize,
        order: Order,
    ) -> Vec<GPoly> {
        let mut out = Vec::new();
        for d in 1..n {
            let terms: Vec<(i64, Vec<u32>)> = (0..n)
                .map(|start| {
                    (
                        1,
                        mono(n, &(0..d).map(|j| ((start + j) % n, 1)).collect::<Vec<_>>()),
                    )
                })
                .collect();
            out.push(poly(&terms, order));
        }
        out.push(poly(&[(1, vec![1; n]), (-1, mono(n, &[]))], order));
        out
    }

    #[test]
    fn par_groebner_and_fglm() {
        let lim = GroebnerLimits {
            reductions: 200_000,
            basis: 2000,
        };
        for m in [5, 6] {
            let g = katsura(m, Order::GradedRevLex);
            compare("P3 groebner", &format!("katsura-{m} grevlex"), || {
                groebner(&g, Order::GradedRevLex, lim)
            });
        }
        let c = cyclic(5, Order::GradedRevLex);
        compare("P3 groebner", "cyclic-5 grevlex", || {
            groebner(&c, Order::GradedRevLex, lim)
        });
        let basis = groebner(&c, Order::GradedRevLex, lim).unwrap();
        // FGLM stays sequential (the multiplication columns are not its bottleneck);
        // checked here only for equality of the results under both pools
        compare(
            "P3 fglm (serial, no site)",
            "cyclic-5 -> lex (dim 70)",
            || fglm(&basis, Order::GradedRevLex, Order::Lex),
        );
    }
}

// ---------------------------------------------------------------- P4: factorisation

mod factor_sites {
    use super::compare;
    use num_bigint::BigInt;
    use num_rational::BigRational;
    use rssn::rules::poly::univariate::factor;

    /// Swinnerton-Dyer polynomial of the first `k` primes (degree `2^k`,
    /// irreducible over `Q`, reducible modulo every prime): the worst case
    /// of Zassenhaus recombination.
    pub fn swinnerton_dyer(primes: &[i64]) -> Vec<BigInt> {
        let mut p: Vec<BigInt> = vec![BigInt::from(0), BigInt::from(1)]; // x
        for &q in primes {
            // P(x + s) P(x - s) = A^2 - q B^2 with s^2 = q
            let n = p.len();
            let mut a = vec![BigInt::from(0); n];
            let mut b = vec![BigInt::from(0); n];
            let mut binom: Vec<Vec<BigInt>> = vec![vec![BigInt::from(1)]];
            for j in 1..n {
                let prev = &binom[j - 1];
                let mut row = vec![BigInt::from(1)];
                for i in 1..j {
                    row.push(&prev[i - 1] + &prev[i]);
                }
                row.push(BigInt::from(1));
                binom.push(row);
            }
            for (j, c) in p.iter().enumerate() {
                for m in 0..=j {
                    // c * C(j,m) x^(j-m) s^m
                    let term = c * &binom[j][m] * BigInt::from(q).pow((m / 2) as u32);
                    if m % 2 == 0 {
                        a[j - m] += term;
                    } else {
                        b[j - m] += term;
                    }
                }
            }
            let mul = |u: &[BigInt], v: &[BigInt]| {
                let mut out = vec![BigInt::from(0); u.len() + v.len() - 1];
                for (i, x) in u.iter().enumerate() {
                    for (j, y) in v.iter().enumerate() {
                        out[i + j] += x * y;
                    }
                }
                out
            };
            let a2 = mul(&a, &a);
            let b2 = mul(&b, &b);
            let mut next = a2;
            for (i, c) in b2.iter().enumerate() {
                next[i] -= BigInt::from(q) * c;
            }
            while next.last().is_some_and(|c| *c == BigInt::from(0)) {
                next.pop();
            }
            p = next;
        }
        p
    }

    #[test]
    fn par_zassenhaus_recombination() {
        for primes in [&[2, 3, 5, 7][..], &[2, 3, 5, 7, 11][..]] {
            let sd = swinnerton_dyer(primes);
            let q: Vec<BigRational> = sd
                .iter()
                .map(|c| BigRational::from_integer(c.clone()))
                .collect();
            compare(
                "P4 univariate::factor",
                &format!("Swinnerton-Dyer deg {}", q.len() - 1),
                || factor(&q),
            );
            assert_eq!(factor(&q).1.len(), 1, "irreducible");
        }
    }
}

// ---------------------------------------------------------------- P5: exact elimination

mod qlinalg_sites {
    use super::compare;
    use rssn::kernels::qlinalg::QMat;
    use rssn::kernels::qlinalg::nullspace;
    use rssn::kernels::qlinalg::q;
    use rssn::kernels::qlinalg::rref;

    fn lcg(state: &mut u64) -> i64 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*state >> 33) % 7) as i64 - 3
    }

    #[test]
    fn par_rational_elimination() {
        let mut st = 7_u64;
        // 80 x 80 matrix of small integers with rank deficiency (sums of rows)
        let mut m: QMat = (0..60)
            .map(|_| (0..80).map(|_| q(lcg(&mut st))).collect())
            .collect();
        for i in 0..20 {
            let row: Vec<_> = m[i].iter().zip(&m[i + 20]).map(|(a, b)| a + b).collect();
            m.push(row);
        }
        compare("P5 qlinalg::nullspace", "80 x 80 over Q", || {
            nullspace(&m, 80)
        });
        compare("P5 qlinalg::rref", "80 x 80 over Q", || {
            let mut a = m.clone();
            let p = rref(&mut a, 80);
            (p, a)
        });
    }
}

// ---------------------------------------------------------------- P8-P10: dense / sparse / FD / quadrature

mod numeric_sites {
    use super::compare;
    use rssn::kernels::dense::Mat;
    use rssn::kernels::dense::lu_factor;
    use rssn::kernels::jacobian::fd_jacobian_dense;
    use rssn::kernels::jacobian::fd_jacobian_dense_par;
    use rssn::kernels::krylov::Csr;
    use rssn::kernels::quadrature::GkRule;
    use rssn::kernels::quadrature::LowDiscrepancy;
    use rssn::kernels::quadrature::cubature_nested;
    use rssn::kernels::quadrature::cubature_nested_par;
    use rssn::kernels::quadrature::integrate_gk;
    use rssn::kernels::quadrature::integrate_gk_par;
    use rssn::kernels::quadrature::monte_carlo_par;
    use rssn::kernels::quadrature::quasi_monte_carlo;
    use rssn::kernels::quadrature::quasi_monte_carlo_par;

    fn lcg(state: &mut u64) -> f64 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (*state >> 11) as f64 / (1u64 << 53) as f64
    }

    fn random_mat(
        n: usize,
        st: &mut u64,
    ) -> Mat {
        let mut m = Mat::zeros(n, n);
        for v in &mut m.data {
            *v = lcg(st) - 0.5;
        }
        for i in 0..n {
            m.data[i * n + i] += n as f64 * 0.1;
        }
        m
    }

    #[test]
    fn par_dense_and_sparse() {
        let mut st = 99_u64;
        let a = random_mat(500, &mut st);
        let b = random_mat(500, &mut st);
        compare("P8 Mat::matmul", "500 x 500", || a.matmul(&b).unwrap());
        compare(
            "P8 lu_factor (sequential, no site: parallel update did not pay)",
            "700 x 700",
            || {
                let big = random_mat(700, &mut 5_u64.clone());
                let lu = lu_factor(&big).unwrap();
                lu.det()
            },
        );
        let x: Vec<f64> = (0..a.cols).map(|i| (i as f64).sin()).collect();
        compare("P8 Mat::matvec", "500 x 500 (x200)", || {
            (0..200).map(|_| a.matvec(&x)).next_back().unwrap()
        });
        // 2-D five-point Laplacian on a 400 x 400 grid
        let n = 400;
        let mut trip = Vec::new();
        for i in 0..n {
            for j in 0..n {
                let r = i * n + j;
                trip.push((r, r, 4.0));
                if i > 0 {
                    trip.push((r, r - n, -1.0));
                }
                if i + 1 < n {
                    trip.push((r, r + n, -1.0));
                }
                if j > 0 {
                    trip.push((r, r - 1, -1.0));
                }
                if j + 1 < n {
                    trip.push((r, r + 1, -1.0));
                }
            }
        }
        let csr = Csr::from_triplets(n * n, n * n, &trip);
        let v: Vec<f64> = (0..n * n).map(|i| ((i % 97) as f64) * 0.01).collect();
        compare(
            "P8 Csr::matvec",
            "160000 rows, 5-pt Laplacian (x50)",
            || {
                let mut y = vec![0.0; n * n];
                for _ in 0..50 {
                    csr.matvec(&v, &mut y);
                }
                y
            },
        );
    }

    #[test]
    fn par_fd_jacobian_and_quadrature() {
        // an expensive right-hand side: 200 unknowns, inner loop of transcendental calls
        let n = 200;
        let f = |_t: f64, y: &[f64], out: &mut [f64]| {
            for i in 0..y.len() {
                let mut s = 0.0;
                for k in 0..400 {
                    s += (y[(i + k) % y.len()] * (k as f64 + 1.0)).sin();
                }
                out[i] = s;
            }
        };
        let y: Vec<f64> = (0..n).map(|i| 0.01 * i as f64).collect();
        let mut f0 = vec![0.0; n];
        f(0.0, &y, &mut f0);
        compare(
            "P9 fd_jacobian_dense_par",
            "n = 200, 400 sin per component",
            || fd_jacobian_dense_par(&f, 0.0, &y, &f0),
        );
        let serial = fd_jacobian_dense(&f, 0.0, &y, &f0);
        assert_eq!(
            serial,
            fd_jacobian_dense_par(&f, 0.0, &y, &f0),
            "parallel Jacobian differs from the serial function"
        );
        // expensive nested integrand: a 3-D cubature whose innermost evaluation is costly
        let g = |p: &[f64]| {
            (p[0] * p[1]).sin() * (p[1] + p[2]).cos()
                + (0..200)
                    .map(|k| (p[2] * k as f64).sin() * 1e-6)
                    .sum::<f64>()
        };
        let lo = [0.0; 3];
        let hi = [1.0; 3];
        compare("P10 cubature_nested_par", "3-D, tol 1e-6", || {
            cubature_nested_par(g, &lo, &hi, 1e-6)
        });
        let s = cubature_nested(g, &lo, &hi, 1e-6);
        let p = cubature_nested_par(g, &lo, &hi, 1e-6);
        assert_eq!(
            (s.value, s.evaluations),
            (p.value, p.evaluations),
            "nested cubature differs from the serial one"
        );
        let h = |x: f64| {
            (0..2000)
                .map(|k| (x * (1.0 + k as f64 * 1e-3)).sin())
                .sum::<f64>()
                .abs()
                .sqrt()
        };
        compare("P10 integrate_gk_par", "[0,4], 2000 sin per node", || {
            integrate_gk_par(h, 0.0, 4.0, GkRule::G7K15, 1e-12, 1e-10, 400)
        });
        let s = integrate_gk(h, 0.0, 4.0, GkRule::G7K15, 1e-12, 1e-10, 400);
        let p = integrate_gk_par(h, 0.0, 4.0, GkRule::G7K15, 1e-12, 1e-10, 400);
        assert_eq!(s, p, "parallel GK differs from the serial one");
        let w = |x: &[f64]| (x[0] + x[1] * x[2]).sin() * (x[3] - x[0]).cos();
        compare(
            "P10 quasi_monte_carlo_par (Sobol)",
            "4-D, 4M points",
            || quasi_monte_carlo_par(w, &[0.0; 4], &[1.0; 4], 4_000_000, LowDiscrepancy::Sobol),
        );
        compare("P10 monte_carlo_par", "4-D, 4M points", || {
            monte_carlo_par(w, &[0.0; 4], &[1.0; 4], 4_000_000, 42)
        });
        // the parallel QMC agrees with the serial one up to summation order
        let a = quasi_monte_carlo(w, &[0.0; 4], &[1.0; 4], 1 << 16, LowDiscrepancy::Sobol);
        let b = quasi_monte_carlo_par(w, &[0.0; 4], &[1.0; 4], 1 << 16, LowDiscrepancy::Sobol);
        assert!((a.value - b.value).abs() < 1e-12, "{a:?} {b:?}");
        let a = quasi_monte_carlo(w, &[0.0; 4], &[1.0; 4], 100_000, LowDiscrepancy::Halton);
        let b = quasi_monte_carlo_par(w, &[0.0; 4], &[1.0; 4], 100_000, LowDiscrepancy::Halton);
        assert!((a.value - b.value).abs() < 1e-12, "{a:?} {b:?}");
    }
}

#[test]
fn par_sde_ensemble() {
    use rssn::kernels::sde::EnsembleScheme;
    use rssn::kernels::sde::sde_ensemble;
    let drift = |_t: f64, y: &[f64], o: &mut [f64]| {
        o[0] = 0.8 * y[0];
        o[1] = -y[1];
    };
    let diff = |_t: f64, y: &[f64], o: &mut [f64]| {
        o[0] = 0.6 * y[0];
        o[1] = 0.0;
        o[2] = 0.0;
        o[3] = 0.3 * y[1].abs().sqrt();
    };
    compare(
        "A4 sde_ensemble (Milstein)",
        "20000 paths x 200 steps, n=2",
        || {
            sde_ensemble(
                drift,
                diff,
                2,
                (0.0, 0.005, 200),
                &[1.0, 1.0],
                20_000,
                5,
                EnsembleScheme::EulerMaruyama,
            )
            .unwrap()
            .iter()
            .map(|v| v[0] + v[1])
            .sum::<f64>()
        },
    );
}

#[test]
fn par_stiff_solvers_with_parallel_jacobians() {
    use rssn::kernels::ode_adaptive::OdeOptions;
    use rssn::kernels::ode_stiff::BdfOptions;
    use rssn::kernels::ode_stiff::bdf_with;
    use rssn::kernels::ode_stiff::bdf_with_par;
    // Stiff reaction-diffusion chain; every component runs a short sine series.
    let make = |n: usize| {
        move |_t: f64, y: &[f64], out: &mut [f64]| {
            for i in 0..n {
                let left = if i == 0 { y[i] } else { y[i - 1] };
                let right = if i + 1 == n { y[i] } else { y[i + 1] };
                let mut react = 0.0;
                for k in 1..=40 {
                    react += (f64::from(k) * y[i]).sin() / f64::from(k * k);
                }
                out[i] = 400.0 * (left - 2.0 * y[i] + right) - 0.5 * react * y[i];
            }
        }
    };
    let opts = OdeOptions {
        rtol: 1e-6,
        atol: 1e-9,
        ..OdeOptions::default()
    };
    let summary =
        |s: rssn::kernels::ode_adaptive::OdeSolution| (s.t.len(), s.nfev, s.y.last().cloned());
    for n in [100usize, 200, 400] {
        let rhs = make(n);
        let y0: Vec<f64> = (0..n).map(|i| 1.0 + 0.1 * (i as f64 * 0.3).sin()).collect();
        let t1 = if n > 200 { 0.05 } else { 0.1 };
        let a = bdf_with(rhs, 0.0, t1, &y0, &opts, &BdfOptions::default()).unwrap();
        let b = bdf_with_par(rhs, 0.0, t1, &y0, &opts, &BdfOptions::default()).unwrap();
        assert_eq!(format!("{:?}", summary(a)), format!("{:?}", summary(b)));
        let size = format!("n = {n}, expensive f");
        compare("A5 bdf_with_par", &size, || {
            summary(bdf_with_par(rhs, 0.0, t1, &y0, &opts, &BdfOptions::default()).unwrap())
        });
    }
}

// ---------------------------------------------------------------- A2: several heavy kernels

mod job_level {
    use std::time::Instant;

    use num_rational::BigRational;
    use rayon::prelude::*;
    use rssn::rules::poly::groebner::GroebnerLimits;
    use rssn::rules::poly::groebner::Order;
    use rssn::rules::poly::groebner::groebner;
    use rssn::rules::poly::univariate::factor;

    /// What a two-phase kernel would do: the heavy parts of several nodes of one request run as
    /// parallel jobs (each of which is itself parallel inside) instead of one after the other.
    fn jobs_row<T: std::fmt::Debug + Send + Sync>(
        site: &str,
        size: &str,
        inputs: &[T],
        job: impl Fn(&T) -> String + Sync,
    ) {
        let t0 = Instant::now();
        let a: Vec<String> = inputs.iter().map(&job).collect();
        let ts = t0.elapsed().as_secs_f64();
        let t1 = Instant::now();
        let b: Vec<String> = inputs.par_iter().map(&job).collect();
        let tp = t1.elapsed().as_secs_f64();
        assert_eq!(a, b, "{site}");
        println!(
            "PAR | {site} | {size} | one after the other {:.1} ms | as parallel jobs {:.1} ms | x{:.2} | threads {}",
            ts * 1e3,
            tp * 1e3,
            ts / tp.max(1e-12),
            rayon::current_num_threads()
        );
    }

    #[test]
    fn par_two_phase_candidates() {
        let sets: Vec<Vec<BigRational>> = [
            &[2, 3, 5, 7, 11][..],
            &[2, 3, 5, 7, 13],
            &[2, 3, 5, 11, 13],
            &[2, 3, 7, 11, 13],
        ]
        .iter()
        .map(|p| {
            super::factor_sites::swinnerton_dyer(p)
                .iter()
                .map(|c| BigRational::from_integer(c.clone()))
                .collect()
        })
        .collect();
        jobs_row(
            "A2 factor jobs",
            "4 Swinnerton-Dyer deg 32 in one request",
            &sets,
            |q| format!("{:?}", factor(q)),
        );
        let lim = GroebnerLimits {
            reductions: 200_000,
            basis: 2000,
        };
        let systems: Vec<_> = (0..4)
            .map(|k| super::groebner_sites::katsura(4 + (k % 2), Order::GradedRevLex))
            .collect();
        jobs_row(
            "A2 groebner jobs",
            "4 Katsura-4/5 systems in one request",
            &systems,
            |g| format!("{:?}", groebner(g, Order::GradedRevLex, lim)),
        );
    }
}
