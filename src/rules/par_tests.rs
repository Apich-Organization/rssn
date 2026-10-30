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
