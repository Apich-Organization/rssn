//! Seeded resampling through the symbolic operator and value bridge.

use rssn::graph::DimensionMode;
use rssn::kernel::CallCtx;
use rssn::kernel::DiagKind;
use rssn::kernel::TensorValue;
use rssn::kernel::Value;
use rssn::kernel::signature::DimRuleInfo;
use rssn::kernel::signature::TypeExpr;
use rssn::kernels::robust::bootstrap_basic_ci;
use rssn::kernels::robust::bootstrap_ci;
use rssn::kernels::robust::iqr;
use rssn::kernels::robust::mad;
use rssn::kernels::robust::median;
use rssn::prelude::ComputeError;
use rssn::prelude::Config;
use rssn::prelude::Session;

fn bootstrap_args(
    samples: &[f64],
    statistic: &str,
    resamples: u32,
) -> [Value; 5] {
    [
        Value::Tensor(TensorValue::vector(
            samples.iter().map(|&x| Value::Float(x)).collect(),
        )),
        Value::Symbol(statistic.into()),
        Value::Float(0.1),
        Value::Int(resamples.into()),
        Value::Int(17.into()),
    ]
}

fn bounds(result: Value) -> [u64; 2] {
    let Value::Tensor(tensor) = result else {
        panic!("expected a tensor, got {result:?}")
    };
    assert_eq!(tensor.shape, [2]);
    let [lower, upper] = tensor.data.as_slice() else {
        panic!("expected two bounds")
    };
    let bits = |value: &Value| {
        match value {
            | Value::Float(x) => x.to_bits(),
            | Value::Real(x) => x.to_f64().to_bits(),
            | Value::Int(x) => x.to_string().parse::<f64>().unwrap().to_bits(),
            | other => panic!("expected a numerical bound, got {other:?}"),
        }
    };
    [bits(lower), bits(upper)]
}

#[test]
fn order_statistics_on_mutable_rows_match_the_original_seeded_statistic() {
    let data: Vec<_> = (0..299)
        .map(|i| f64::from((i * 7919) % 563) - 271.25)
        .collect();
    let statistics = [
        ("median", median as fn(&[f64]) -> f64),
        ("mad", |x| mad(x, false)),
        ("iqr", iqr),
    ];
    for workers in [1, 4] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .unwrap();
        pool.install(|| {
            let session = Session::new();
            let context = CallCtx::new();
            for (label, statistic) in statistics {
                for (op, (lo, hi)) in [
                    (
                        "bootstrap_ci",
                        bootstrap_ci(&data, statistic, 1293, 0.1, 17),
                    ),
                    (
                        "bootstrap_basic_ci",
                        bootstrap_basic_ci(&data, statistic, 1293, 0.1, 17),
                    ),
                ] {
                    let actual = session
                        .apply_diag(op, &bootstrap_args(&data, label, 1293), &context)
                        .unwrap();
                    assert_eq!(
                        bounds(actual),
                        [lo.to_bits(), hi.to_bits()],
                        "{label} {op}, {workers} workers"
                    );
                }
            }
        });
    }
}

#[test]
fn basic_intervals_advertise_symbol_and_dimensionless_arguments() {
    let session = Session::new();
    let manifest = session.manifest();
    let op = manifest.op("bootstrap_basic_ci").unwrap();
    assert_eq!(op.dim_rule, Some(DimRuleInfo::Dimensionless));
    assert!(op.doc.contains("reverse percentile"));
    assert!(
        op.overloads
            .iter()
            .any(|sig| sig.params.get(1) == Some(&TypeExpr::Symbol))
    );

    let config = Config::new().dimensions(DimensionMode::Strict);
    let good = session
        .parse("bootstrap_basic_ci(list(3.0, 3.0), median, 0.1, 10, 7)")
        .unwrap();
    let answer = session.compute(good, &config).unwrap();
    assert_eq!(answer.dimension.unwrap().to_string(), "1");
    assert_eq!(
        bounds(answer.term.to_value().unwrap()),
        [3.0_f64.to_bits(); 2]
    );
    let measured = session
        .parse("bootstrap_basic_ci(list(quantity(3, m)), median, 0.1, 10, 7)")
        .unwrap();
    assert!(matches!(
        session.compute(measured, &config),
        Err(ComputeError::Dimension(_))
    ));
}

#[test]
fn rejected_bootstrap_arguments_identify_their_slot() {
    let session = Session::new();
    for (index, argument) in [
        (1, Value::Symbol("unknown".into())),
        (2, Value::Float(1.0)),
        (3, Value::Int(0.into())),
        (4, Value::Int((-1).into())),
    ] {
        let mut arguments = bootstrap_args(&[1.0, 3.0, 9.0], "median", 20);
        arguments[index] = argument;
        let error = session
            .apply_diag("bootstrap_basic_ci", &arguments, &CallCtx::new())
            .unwrap_err();
        assert_eq!(error.kind, DiagKind::Domain);
        assert_eq!(error.arg, Some(index));
    }
}
