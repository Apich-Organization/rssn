//! Symbolic and native conditional stochastic integrals.

use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::value::TensorValue;
use rssn::kernels::random::Rng;
use rssn::kernels::sde::levy_areas;

fn session() -> Session {
    Session::with_rules(&[rssn::rules::sde()]).unwrap()
}

fn sample_in_graph(
    session: &Session,
    seed: u64,
    h: f64,
    dw: &[f64],
    terms: usize,
) -> Value {
    session
        .call_value(
            "sde_levy_areas",
            &[
                Value::int(seed),
                Value::Float(h),
                Value::Tensor(TensorValue::vector(
                    dw.iter().copied().map(Value::Float).collect(),
                )),
                Value::int(terms),
            ],
            &CallCtx::new(),
        )
        .unwrap()
}

#[test]
fn symbolic_conditional_areas_decline_nonfinite_or_excessive_requests() {
    let session = session();
    for increments in [[0.1, f64::NAN], [0.1, f64::INFINITY], [0.1, 1e200]] {
        assert!(matches!(
            sample_in_graph(&session, 7, 0.01, &increments, 8),
            Value::Term(_)
        ));
    }
    // An independent bound applies to each of the output size and work.
    assert!(matches!(
        sample_in_graph(&session, 7, 0.01, &vec![0.0; 2049], 1),
        Value::Term(_)
    ));
    assert!(matches!(
        sample_in_graph(&session, 7, 0.01, &vec![0.0; 368], 1_000),
        Value::Term(_)
    ));

    // A representable physical integral can have an overflowing unscaled square.
    let Value::Tensor(answer) = sample_in_graph(&session, 7, 1.5e308, &[1.5e154], 1) else {
        panic!("the scalar area is representable")
    };
    assert!((answer.data[0].as_f64().unwrap() / 1e307 - 3.75).abs() < 1e-13);
}

#[test]
fn symbolic_and_native_conditional_areas_share_the_same_seeded_sample() {
    let seed = u64::MAX;
    let h = 0.375;
    let increments = [0.125, -0.25, 0.5];
    let terms = 23;
    let Value::Tensor(answer) = sample_in_graph(&session(), seed, h, &increments, terms) else {
        panic!("expected a flat list of areas")
    };
    assert_eq!(answer.shape, [9]);
    let expected = levy_areas(&mut Rng::new(seed), increments.len(), h, &increments, terms);
    for (actual, expected) in answer.data.iter().zip(expected) {
        assert_eq!(actual.as_f64().unwrap().to_bits(), expected.to_bits());
    }
}

#[test]
fn native_conditional_areas_preserve_the_rng_and_results_across_worker_counts() {
    let widths = [(1, 7), (48, 257)];
    let compute = || {
        let mut results = Vec::new();
        for (m, terms) in widths {
            let h = 0.125_f64;
            let dw: Vec<f64> = (0..m).map(|j| (j as f64 * 0.75).sin() * h.sqrt()).collect();
            let mut rng = Rng::new(3321);
            let mut expected_rng = rng.clone();
            for _ in 0..2 * m * terms + m + m * (m - 1) / 2 {
                let _ = expected_rng.normal();
            }
            let area = levy_areas(&mut rng, m, h, &dw, terms);
            assert_eq!(area.len(), m * m);
            for j in 0..m {
                assert_eq!(area[j * m + j], 0.5 * (dw[j] * dw[j] - h));
                for l in j + 1..m {
                    assert!((area[j * m + l] + area[l * m + j] - dw[j] * dw[l]).abs() < 1e-15);
                }
            }
            assert_eq!(rng.normal().to_bits(), expected_rng.normal().to_bits());
            assert_eq!(rng.next_u64(), expected_rng.next_u64());
            results.extend(area.into_iter().map(f64::to_bits));
        }
        results
    };
    let one = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let four = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap();
    assert_eq!(one.install(compute), four.install(compute));
}
