//! Regressions for forward-backward filtering and the measured value contract.

use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::overload::matches_type;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::TensorValue;
use rssn::kernels::signal;

fn tensor(samples: &[f64]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![samples.len()],
        data: samples.iter().copied().map(Value::Float).collect(),
        elem: TypeExpr::Float,
    })
}

#[test]
fn forward_backward_matches_the_difference_equation_with_odd_padding() {
    let b = [1.0, 0.5, 0.125];
    let a = [2.0, -0.5, 0.125];
    let samples = [1.0, -2.0, 0.5, 3.0, -0.25];
    // Independent rational direct-form-I recurrence, with a constant
    // prehistory at the outer edge of each reflected pass.
    let expected = [
        1.003382862847502,
        0.010159065350308083,
        0.48553389334119856,
        0.9020575687754899,
        -0.25463422294706106,
    ];
    let filtered = signal::checked_filtfilt(&b, &a, &samples).unwrap();
    assert_eq!(filtered.len(), samples.len());
    for (&observed, reference) in filtered.iter().zip(expected) {
        assert!(
            (observed - reference).abs() < 4e-15,
            "{observed} vs {reference}"
        );
    }
    assert_eq!(signal::filtfilt(&b, &a, &samples), filtered);
}

#[test]
fn singleton_and_constant_signals_share_the_steady_state() {
    for samples in [vec![3.0], vec![3.0; 20]] {
        assert_eq!(
            signal::checked_filtfilt(&[0.25], &[0.5, -0.25], &samples).unwrap(),
            samples
        );
    }
    // A finite reflected endpoint should survive even when 2 * endpoint
    // would overflow if evaluated separately before subtracting the sample.
    let samples = [f64::MAX; 3];
    assert_eq!(
        signal::checked_filtfilt(&[0.5, 0.0], &[1.0], &samples).unwrap(),
        [0.25 * f64::MAX; 3]
    );
    assert!(
        signal::checked_filtfilt(&[1.0], &[1.0], &[])
            .unwrap()
            .is_empty()
    );
}

#[test]
fn checked_filter_reports_invalid_parameters_and_unrepresentable_results() {
    use signal::FilterError;
    use signal::SignalError;
    for (b, a) in [
        (&[][..], &[1.0][..]),
        (&[1.0], &[0.0]),
        (&[1.0], &[1.0, -1.0]),
    ] {
        assert_eq!(
            signal::checked_filtfilt(b, a, &[2.0]),
            Err(FilterError::InvalidParameters)
        );
    }
    assert_eq!(
        signal::checked_filtfilt(&[1.0], &[1.0], &[f64::NAN]),
        Err(FilterError::Signal(SignalError::NonFiniteInput))
    );
    assert_eq!(
        signal::checked_filtfilt(&[2.0], &[1.0], &[f64::MAX]),
        Err(FilterError::Signal(SignalError::NonFiniteResult))
    );
}

#[test]
fn numerical_filter_keeps_sample_units_and_requires_unitless_coefficients() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let dim = DimExpr::base("Length");
    let measured_tensor = |samples: &[f64]| {
        Value::Tensor(TensorValue {
            shape: vec![samples.len()],
            data: samples.iter().copied().map(Value::Float).collect(),
            elem: TypeExpr::Float.with_dim(dim.clone()),
        })
    };
    let measured = measured_tensor(&[3.0, 3.0, 3.0]);
    let b = tensor(&[0.25]);
    let a = tensor(&[0.5, -0.25]);
    let observed = session
        .call_value(
            "num_signal_filtfilt",
            &[b, a.clone(), measured.clone()],
            &ctx,
        )
        .unwrap();
    let Value::Tensor(observed) = observed else {
        panic!("expected a measured tensor, got {observed:?}");
    };
    assert_eq!(observed.shape, [3]);
    assert!(matches_type(&TypeExpr::Float.with_dim(dim.clone()), &observed.elem).accepted());
    for cell in observed.data {
        let Value::Dyn(cell) = cell else {
            panic!("expected a measured sample, got {cell:?}");
        };
        assert_eq!(cell.dim, dim);
        assert_eq!(*cell.value, Value::Float(3.0));
    }

    assert!(
        session
            .call_value(
                "num_signal_filtfilt",
                &[measured_tensor(&[0.25]), a, measured],
                &ctx
            )
            .is_err()
    );
    assert!(matches!(
        session.call_value(
            "num_signal_filtfilt",
            &[tensor(&[1.0]), tensor(&[1.0, -1.0]), tensor(&[2.0])],
            &ctx
        ),
        Err(OpError::Domain(_))
    ));
}
