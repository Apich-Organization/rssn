//! Sequence acceleration on classical slowly converging series.

use std::f64::consts::LN_2;
use std::f64::consts::PI;

use rssn::kernels::convergence::LevinVariant;
use rssn::kernels::convergence::aitken_acceleration;
use rssn::kernels::convergence::find_sequence_limit;
use rssn::kernels::convergence::levin_condition_number;
use rssn::kernels::convergence::levin_sum;
use rssn::kernels::convergence::levin_sum_lazy;
use rssn::kernels::convergence::levin_transform;
use rssn::kernels::convergence::levin_transform_terms;
use rssn::kernels::convergence::richardson_extrapolation;
use rssn::kernels::convergence::richardson_extrapolation_with;
use rssn::kernels::convergence::sum_series_numerical;
use rssn::kernels::convergence::wynn_epsilon;

fn partial_sums(
    term: impl Fn(f64) -> f64,
    first: f64,
    count: usize,
) -> Vec<f64> {
    let mut sum = 0.0;
    (0..count)
        .map(|i| {
            sum += term(first + i as f64);
            sum
        })
        .collect()
}

fn last(v: &[f64]) -> f64 {
    v.last().copied().unwrap_or(f64::NAN)
}

#[test]
fn aitken_on_leibniz_and_alternating_harmonic() {
    let leibniz = partial_sums(|k| (-1.0_f64).powf(k) / (2.0 * k + 1.0), 0.0, 20);
    let acc = aitken_acceleration(&leibniz);
    assert_eq!(acc.len(), 18);
    assert!((last(&acc) - PI / 4.0).abs() < (leibniz[19] - PI / 4.0).abs() / 100.0);
    let harmonic = partial_sums(|k| (-1.0_f64).powf(k + 1.0) / k, 1.0, 20);
    assert!((last(&aitken_acceleration(&harmonic)) - LN_2).abs() < 1e-3);
    assert!(aitken_acceleration(&[1.0, 2.0]).is_empty());
}

#[test]
fn aitken_is_exact_for_geometric_error() {
    let seq: Vec<f64> = (0..6).map(|n| 3.0 + 0.5_f64.powi(n)).collect();
    for v in aitken_acceleration(&seq) {
        assert!((v - 3.0).abs() < 1e-12, "{v}");
    }
}

#[test]
fn wynn_epsilon_on_alternating_series() {
    let leibniz = partial_sums(|k| (-1.0_f64).powf(k) / (2.0 * k + 1.0), 0.0, 13);
    let est = wynn_epsilon(&leibniz);
    assert!((last(&est) - PI / 4.0).abs() < 1e-8, "{est:?}");
    assert_eq!(est[0], leibniz[12]);
    let harmonic = partial_sums(|k| (-1.0_f64).powf(k + 1.0) / k, 1.0, 13);
    assert!((last(&wynn_epsilon(&harmonic)) - LN_2).abs() < 1e-8);
    assert!(wynn_epsilon(&[]).is_empty());
    // Exactly converged sequences are returned as they are.
    assert_eq!(last(&wynn_epsilon(&[1.0, 1.0, 1.0])), 1.0);
}

#[test]
fn richardson_on_basel_partial_sums_at_doubling_counts() {
    let counts = [8_u32, 16, 32, 64, 128, 256, 512];
    let sums: Vec<f64> = counts
        .iter()
        .map(|&n| (1..=n).map(|k| 1.0 / f64::from(k * k)).sum())
        .collect();
    let rich = richardson_extrapolation_with(&sums, 2.0);
    assert!(
        (last(&rich) - PI * PI / 6.0).abs() < 1e-8,
        "{}",
        last(&rich)
    );
    assert!((last(&rich) - PI * PI / 6.0).abs() < (sums[6] - PI * PI / 6.0).abs() * 1e-5);
}

#[test]
fn richardson_removes_h_squared_errors() {
    // Trapezoid-like approximations of 1 with error c2 h^2 + c4 h^4.
    let seq: Vec<f64> = (0..4)
        .map(|i| {
            let h = 0.5_f64.powi(i);
            1.0 + 3.0 * h * h - 2.0 * h.powi(4)
        })
        .collect();
    let rich = richardson_extrapolation(&seq);
    assert_eq!(rich[0], seq[0]);
    assert!((rich[2] - 1.0).abs() < 1e-12, "{rich:?}");
    assert!(richardson_extrapolation(&[]).is_empty());
}

#[test]
fn find_sequence_limit_finds_known_limits() {
    let l =
        find_sequence_limit(|n| 2.0 + 0.5_f64.powf(n), 20, 1e-8).unwrap_or_else(|e| panic!("{e}"));
    assert!((l - 2.0).abs() < 1e-8);
    // Partial sums of the Leibniz series.
    let leibniz = partial_sums(|k| (-1.0_f64).powf(k) / (2.0 * k + 1.0), 0.0, 40);
    let l =
        find_sequence_limit(|n| leibniz[n as usize], 40, 1e-6).unwrap_or_else(|e| panic!("{e}"));
    assert!((l - PI / 4.0).abs() < 1e-5, "{l}");
    assert!(find_sequence_limit(|_| f64::NAN, 10, 1e-8).is_err());
    assert!(find_sequence_limit(|n| n, 10, 1e-8).is_err());
}

#[test]
fn sum_series_numerical_stops_at_small_terms() {
    let s = sum_series_numerical(|n| 0.5_f64.powf(n), 0, 200, 1e-14);
    assert!((s - 2.0).abs() < 1e-12);
    let s = sum_series_numerical(|n| n, 1, 10, 0.0);
    assert_eq!(s, 55.0);
}

fn ps(
    term: impl Fn(usize) -> f64,
    count: usize,
) -> Vec<f64> {
    let mut s = 0.0;
    (0..count)
        .map(|i| {
            s += term(i);
            s
        })
        .collect()
}

#[test]
fn levin_alternating_ln2() {
    // sum (-1)^k / (k+1) = ln 2
    let t = |k: usize| if k.is_multiple_of(2) { 1.0 } else { -1.0 } / (k as f64 + 1.0);
    let sums = ps(t, 16);
    for variant in [LevinVariant::T, LevinVariant::U, LevinVariant::V] {
        let est = levin_transform(&sums, variant, 1.0);
        let best = est
            .iter()
            .map(|e| (e - LN_2).abs())
            .fold(f64::INFINITY, f64::min);
        assert!(best < 1e-10, "{variant:?}: {best}");
        // far better than the raw partial sum
        assert!((sums[15] - LN_2).abs() > 1e-2);
    }
    let (v, err) = levin_sum(t, 20, LevinVariant::T).unwrap();
    assert!((v - LN_2).abs() < 1e-12, "{v}");
    assert!(err < 1e-9);
}

#[test]
fn levin_logarithmic_zeta2() {
    // sum 1/(k+1)^2 = pi^2/6: logarithmic convergence, Aitken is poor
    let t = |k: usize| 1.0 / ((k + 1) as f64).powi(2);
    let target = PI * PI / 6.0;
    let (v, err) = levin_sum(t, 24, LevinVariant::U).unwrap();
    assert!(
        (v - target).abs() < 1e-10,
        "u: {} err {err}",
        (v - target).abs()
    );
    // zeta(3) (Apery) with the u transform
    let (v, _) = levin_sum(|k| 1.0 / ((k + 1) as f64).powi(3), 24, LevinVariant::U).unwrap();
    assert!((v - 1.202_056_903_159_594_2).abs() < 1e-10);
    // compare: Aitken on the same data is much worse
    let sums = ps(t, 24);
    let aitken = aitken_acceleration(&sums);
    let a_err = (aitken.last().unwrap() - target).abs();
    assert!(a_err > 1e-6, "aitken {a_err}");
}

#[test]
fn levin_edge_cases() {
    assert!(levin_transform(&[], LevinVariant::U, 1.0).is_empty());
    assert!(levin_sum(|_| 1.0, 2, LevinVariant::U).is_err());
    // geometric series is summed essentially exactly
    let (v, _) = levin_sum(|k| 0.5_f64.powi(k as i32), 12, LevinVariant::T).unwrap();
    assert!((v - 2.0).abs() < 1e-10);
}

#[test]
fn levin_u_on_terms_beats_partial_sums_on_hard_series() {
    let zeta2 = |k: usize| 1.0 / ((k + 1) as f64).powi(2);
    let ln2 = |k: usize| {
        let sign = if k.is_multiple_of(2) {
            1.0
        } else {
            -1.0
        };
        sign / ((k + 1) as f64)
    };
    let zeta3 = |k: usize| 1.0 / ((k + 1) as f64).powi(3);
    let slow_tail = |k: usize| {
        let sign = if k.is_multiple_of(2) {
            1.0
        } else {
            -1.0
        };
        sign / ((k + 1) as f64).sqrt()
    };

    let cases: [(fn(usize) -> f64, f64, &str); 4] = [
        (zeta2, PI * PI / 6.0, "zeta2"),
        (ln2, LN_2, "ln2"),
        (zeta3, 1.202_056_903_159_594_2, "zeta3"),
        (slow_tail, 0.604_898_643_421_630_4, "eta(1/2)"),
    ];
    for (term, exact, name) in cases {
        let sums = ps(term, 32);
        let raw = (sums[31] - exact).abs();
        let (value, _) = levin_sum_lazy((0..32).map(term), 32, LevinVariant::U, 1.0)
            .unwrap_or_else(|e| panic!("{name}: {e}"));
        let err = (value - exact).abs();
        println!(
            "{name:10} raw {raw:.3e}  levin {err:.3e}  gain {:.1e}",
            raw / err
        );
        // Every one of these converges too slowly to reach 1e-9 by summing
        // 32 terms, and Levin reaches it.
        assert!(raw > 1e-4, "{name}: {raw}");
        assert!(err < 1e-9, "{name}: {err}");
    }
    // The logarithmic case: nine orders of magnitude on 32 terms.
    let sums = ps(zeta2, 32);
    let raw = (sums[31] - PI * PI / 6.0).abs();
    let (value, _) = levin_sum_lazy((0..32).map(zeta2), 32, LevinVariant::U, 1.0).unwrap();
    assert!(raw / (value - PI * PI / 6.0).abs() > 1e8);
}

#[test]
fn levin_terms_and_partial_sums_agree_low_order_and_differ_high_order() {
    let zeta2 = |k: usize| 1.0 / ((k + 1) as f64).powi(2);
    let terms: Vec<f64> = (0..17).map(zeta2).collect();
    let sums = ps(zeta2, 17);
    let by_terms = levin_transform_terms(&terms, LevinVariant::U, 1.0);
    let by_sums = levin_transform(&sums, LevinVariant::U, 1.0);
    // Order 8 is well below the conditioning limit: both paths agree.
    assert!((by_terms[8] - by_sums[8]).abs() < 1e-11);
    // At order 16 the differenced terms have lost the digits the terms keep.
    let exact = PI * PI / 6.0;
    let terms_err = (by_terms[16] - exact).abs();
    let sums_err = (by_sums[16] - exact).abs();
    println!("order 16: terms {terms_err:.3e}  differenced sums {sums_err:.3e}");
    assert!(terms_err < sums_err / 10.0, "{terms_err} vs {sums_err}");
}

#[test]
fn levin_condition_number_says_where_the_order_stops_helping() {
    let zeta2: Vec<f64> = (0..20).map(|k| 1.0 / ((k + 1) as f64).powi(2)).collect();
    let k11 = levin_condition_number(&zeta2[..12], LevinVariant::U, 1.0);
    let k16 = levin_condition_number(&zeta2[..17], LevinVariant::U, 1.0);
    assert!((k11 / 4.6e5 - 1.0).abs() < 0.1, "{k11}");
    assert!((k16 / 2.3e8 - 1.0).abs() < 0.1, "{k16}");
    // Roughly a factor 3.5 per order, so the floor kappa * eps rises from
    // 1e-10 at order 11 to 5e-8 at order 16.
    assert!(k11 * f64::EPSILON > 5e-11 && k11 * f64::EPSILON < 2e-10);
    assert!(k16 * f64::EPSILON > 4e-8 && k16 * f64::EPSILON < 6e-8);

    // The measured error follows that floor: best at order 11, worse at 19.
    let est = levin_transform_terms(&zeta2, LevinVariant::U, 1.0);
    let exact = PI * PI / 6.0;
    let e11 = (est[11] - exact).abs();
    let e14 = (est[14] - exact).abs();
    let e19 = (est[19] - exact).abs();
    println!("order 11 {e11:.3e}  14 {e14:.3e}  19 {e19:.3e}");
    assert!(e11 < 1e-10);
    assert!(e14 > e11);
    assert!(e19 > 1e-7);

    // Alternating terms are a convex combination: no conditioning limit.
    let alt: Vec<f64> = (0..17_usize)
        .map(|k| {
            let sign = if k.is_multiple_of(2) {
                1.0
            } else {
                -1.0
            };
            sign / ((k + 1) as f64)
        })
        .collect();
    assert!((levin_condition_number(&alt, LevinVariant::U, 1.0) - 1.0).abs() < 1e-12);
}

#[test]
fn levin_sums_an_abel_divergent_series() {
    // sum_{n>=1} (-1)^n n = -1/4 in the Abel sense.
    let terms = (1..25_usize).map(|n| {
        let sign = if n % 2 == 1 { -1.0 } else { 1.0 };
        sign * (n as f64)
    });
    let (value, _) = levin_sum_lazy(terms, 24, LevinVariant::U, 1.0).unwrap();
    assert!((value + 0.25).abs() < 1e-12, "{value}");
}
