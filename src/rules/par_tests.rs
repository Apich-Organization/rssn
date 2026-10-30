//! Determinism and timing checks of the symbolic parallel sites (Part 6 of
//! `docs/ALGORITHM_COMPLETION_CHECKLIST.md`). Each site runs in a one-thread
//! pool and in the global pool; the results must be equal. The timings are
//! printed (`--nocapture`, grep `PAR`), never asserted.

use std::time::Instant;

use super::number_theory::number_theory;
use super::testing::simplify;

/// Runs `f` serially (one-thread pool) and in the global pool, asserts the
/// results agree and prints one timing row.
fn compare<T: std::fmt::Debug + Send>(
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

#[test]
fn par_ecm_factorint() {
    for src in [
        "factorint(2^128 + 1)",
        "factorint(2^101 - 1)",
        "factorint((2^107 - 1) * (2^61 - 1))",
    ] {
        compare("P7 factorint (ECM curves)", src, || {
            simplify(&[number_theory()], src)
        });
    }
}

#[test]
fn par_soundness_check_program() {
    use crate::graph::Engine;
    use crate::graph::Graph;
    use crate::graph::soundness::check_program;
    let mut graph = Graph::new();
    let engine = Engine::install(
        &mut graph,
        &[
            super::arith(),
            super::elementary(),
            super::poly(),
            super::complex(),
        ],
    )
    .unwrap();
    compare(
        "P15 check_program",
        "arith+elementary+poly+complex rewrites, 60 samples",
        || {
            check_program(&graph, engine.program(), 60)
                .unwrap()
                .into_iter()
                .map(|(name, v)| (name, v.agreed, v.inconclusive))
                .collect::<Vec<_>>()
        },
    );
}

#[test]
fn par_galois_resolvents() {
    use super::poly::galois_high::analyse;
    use num_bigint::BigInt;
    use num_rational::BigRational;
    let poly = |c: &[i64]| -> Vec<BigRational> {
        c.iter()
            .map(|&v| BigRational::from_integer(BigInt::from(v)))
            .collect()
    };
    // S7 (x^7 - x - 1): the whole coset list of 5040 permutations is searched
    // for every candidate group; PSL(3,2) (Trinks): a later group in the list
    for (name, coefficients) in [
        ("x^7 - x - 1 (S7)", [-1, -1, 0, 0, 0, 0, 0, 1]),
        ("x^7 - 7x + 3 (PSL(3,2))", [3, -7, 0, 0, 0, 0, 0, 1]),
    ] {
        let c = poly(&coefficients);
        compare("P6 galois_high::analyse", name, || {
            analyse(&c).map(|a| (a.group.label, a.group.order))
        });
    }
}

#[test]
fn par_zeilberger_systems() {
    let rules = super::standard();
    for f in [
        "binomial(n, k)^2*binomial(n + k, k)^2",
        "(-1)^k*binomial(n, k)^3",
    ] {
        let request = format!("zeilberger({f}, k, n)");
        compare("P5 zeilberger (not parallelised: no win)", f, || {
            super::testing::reduce_with(&rules, &request, &[]).0
        });
    }
}

#[test]
fn par_engine_batched_matching() {
    let rules = super::standard();
    for src in [
        "diff(sin(x)^2*exp(x^2)*ln(1 + x^2), x)",
        "(x + 1)^6 - (x - 1)^6 + (x^2 + 2*x + 1)*(x^2 - 2*x + 1)",
        "integral(x^3*exp(x) + sin(x)*cos(x) + 1/(x^2 + 1), x)",
        "diff(((x^2 + 1)^5*sin(x)^3*exp(x)*ln(x^2 + 2))^2, x)",
        "(x + y + z)^6 - (x - y + z)^6 + (x*y + y*z + z*x)^3 - (x + y)*(y + z)*(z + x)*(x*y*z)",
        "diff(diff(sin(x)*cos(x)*exp(x)*ln(x)*atan(x)*sinh(x), x), x)",
        "diff(sin(x)^3*cos(x)^2*exp(2*x)*ln(x + 1)*atan(x)^2, x)",
        "diff(exp(x)*sin(x)*cos(x)*ln(x)*sinh(x)*atan(x)*tan(x), x)",
    ] {
        compare("A1 Engine::apply_tier batched e-matching", src, || {
            super::testing::reduce_with(&rules, src, &[]).0
        });
    }
}

#[test]
fn par_session_compute_all() {
    use crate::api::Config;
    use crate::api::Session;
    let texts = [
        "diff(sin(x)^2*exp(x^2)*ln(1 + x^2), x)",
        "(x + 1)^6 - (x - 1)^6 + (x^2 + 2*x + 1)*(x^2 - 2*x + 1)",
        "integral(x^3*exp(x) + sin(x)*cos(x) + 1/(x^2 + 1), x)",
        "diff(x^5*atan(x)/(1 + x^2), x)",
        "integral(x^2*sin(x), x)",
        "diff(exp(sin(x))*cos(x^3), x)",
        "(x^2 - 1)^3 - (x - 1)^3*(x + 1)^3",
        "integral(exp(2*x)*cos(x), x)",
    ];
    // Returns the answers and the seconds spent in the requests (the
    // template installation is excluded).
    let job = |batch: bool| -> (Vec<String>, f64) {
        let s = Session::new();
        let config = Config::new();
        let _ = s.compute(s.parse("x + x").unwrap(), &config);
        let terms: Vec<_> = texts.iter().map(|t| s.parse(t).unwrap()).collect();
        let t0 = Instant::now();
        let answers: Vec<_> = if batch {
            s.compute_all(&terms, &config)
        } else {
            terms.iter().map(|&t| s.compute(t, &config)).collect()
        };
        let dt = t0.elapsed().as_secs_f64();
        let out = answers
            .into_iter()
            .map(|a| a.map_or_else(|e| format!("{e}"), |a| a.term.to_string()))
            .collect();
        (out, dt)
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let (looped, tl) = job(false);
    let (serial, ts) = pool.install(|| job(true));
    let (parallel, tp) = job(true);
    assert_eq!(serial, looped);
    assert_eq!(parallel, looped);
    println!(
        "PAR | A3 Session::compute_all | 8 requests, standard rules | loop {:.1} ms | serial {:.1} ms | parallel {:.1} ms | x{:.2} | threads {}",
        tl * 1e3,
        ts * 1e3,
        tp * 1e3,
        ts / tp.max(1e-12),
        rayon::current_num_threads()
    );
}

#[test]
fn par_backend_plot_and_verify_sampling() {
    use crate::backend::call_batch_par;
    use crate::backend::compile;
    use crate::graph::Env;
    use crate::graph::Graph;
    let mut g = Graph::new();
    let _engine = crate::graph::Engine::install(&mut g, &super::standard()).unwrap();
    let term = g
        .parse(
            "sin(x)*exp(-x^2/9) + cos(y)*(x^2 + y^2 + 1) + x*y/(1 + x^2) + sin(3*y)*cos(x)*cos(x)",
        )
        .unwrap();
    let (x, y) = (g.interner_mut().symbol("x"), g.interner_mut().symbol("y"));
    let compiled = compile(&g, term, &[x, y]).unwrap_or_else(|e| panic!("{e}"));
    let n = 1_000_000;
    let xs: Vec<f64> = (0..n).map(|i| -5.0 + 10.0 * i as f64 / n as f64).collect();
    let ys: Vec<f64> = (0..n).map(|i| 3.0 * (i as f64 * 1e-3).sin()).collect();
    compare("A6 backend::call_batch_par", "1M rows, 2 inputs", || {
        let mut out = vec![0.0; n];
        call_batch_par(&*compiled, &[&xs, &ys], &mut out);
        out
    });
    // plot sampling through Graph::eval
    let mut g2 = Graph::new();
    let _engine2 = crate::graph::Engine::install(&mut g2, &super::standard()).unwrap();
    let t2 = g2
        .parse(
            "sin(x)*exp(-x^2/9) + cos(3*x)*sqrt(x^2 + 1) + atan(x)*ln(x^2 + 2) + sin(x)^2*cos(x)^3",
        )
        .unwrap();
    let dir = std::env::temp_dir();
    compare(
        "A6 io::plot_term_2d",
        "20000 samples (evaluation + SVG)",
        || {
            let path = dir.join(format!("rssn_plot_{}.svg", std::process::id()));
            super::super::io::plot::plot_term_2d(
                &mut g2.clone(),
                t2,
                "x",
                (-8.0, 8.0),
                20_000,
                &path,
            )
            .unwrap();
            let text = std::fs::read_to_string(&path).unwrap();
            let _ = std::fs::remove_file(&path);
            text.len()
        },
    );
    // the sampling pattern of rules/verify.rs: 12 evaluations of the simplified term
    let sym = g2.interner_mut().symbol("x");
    let g2 = g2;
    let eval12 = |parallel: bool| -> Vec<Option<u64>> {
        let one = |s: u32| {
            let mut env = Env::numeric(0.0);
            env.bind(sym, 0.31 + 0.47 * f64::from(s));
            g2.eval(t2, &env).map(f64::to_bits)
        };
        if parallel {
            crate::kernels::par::map_indexed(12, 2, |s| one(s as u32))
        } else {
            (0..12_u32).map(one).collect()
        }
    };
    let reps = 2000;
    let t0 = Instant::now();
    for _ in 0..reps {
        std::hint::black_box(eval12(false));
    }
    let ts = t0.elapsed().as_secs_f64() / f64::from(reps);
    let t1 = Instant::now();
    for _ in 0..reps {
        std::hint::black_box(eval12(true));
    }
    let tp = t1.elapsed().as_secs_f64() / f64::from(reps);
    assert_eq!(eval12(false), eval12(true));
    println!(
        "PAR | A6 rules/verify.rs vanishes (12 samples) | one term | serial {:.1} us | parallel {:.1} us | x{:.2}",
        ts * 1e6,
        tp * 1e6,
        ts / tp.max(1e-12)
    );
}

/// Katsura-`m` system in the variables `prefix0..prefix{m}` as a `groebner(...)` request.
fn katsura_request(
    prefix: &str,
    m: usize,
) -> String {
    let u = |l: i64| format!("{prefix}{}", l.unsigned_abs());
    let mut eqs = vec![{
        let mut terms = vec![u(0)];
        terms.extend((1..=m as i64).map(|l| format!("2*{}", u(l))));
        format!("{} - 1", terms.join(" + "))
    }];
    for k in 0..m as i64 {
        let mut terms = Vec::new();
        for l in -(m as i64)..=(m as i64) {
            if (k - l).unsigned_abs() as usize <= m {
                terms.push(format!("{}*{}", u(l), u(k - l)));
            }
        }
        eqs.push(format!("{} - {}", terms.join(" + "), u(k)));
    }
    let vars: Vec<String> = (0..=m as i64).map(u).collect();
    format!(
        "groebner(list({}), list({}), grevlex)",
        eqs.join(", "),
        vars.join(", ")
    )
}

#[test]
fn par_two_phase_groebner_requests() {
    let rules = super::standard();
    for (label, m, copies) in [("4 Katsura-4", 4, 4), ("4 Katsura-5", 5, 4)] {
        let request = format!(
            "list({})",
            (0..copies)
                .map(|c| katsura_request(&format!("{}", (b'a' + c as u8) as char), m))
                .collect::<Vec<_>>()
                .join(", ")
        );
        compare("A2 two-phase groebner kernel", label, || {
            super::testing::reduce_with(&rules, &request, &[])
        });
    }
}
