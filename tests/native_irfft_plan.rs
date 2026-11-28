//! Prepared inverse real FFTs through the typed, measured native ABI.

use std::f64::consts::TAU;

use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::manifest::SignatureSource;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::NatExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::ComplexValue;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::QuantityValue;
use rssn::kernel::value::ScaleValue;
use rssn::kernel::value::TensorValue;
use rssn::kernel::value::UnitValue;

const NEW: &str = "num_native_irfft_plan_new";
const APPLY: &str = "num_native_irfft_plan_apply";
const OUTPUT_LEN: &str = "num_native_irfft_plan_output_len";
const SPECTRUM_LEN: &str = "num_native_irfft_plan_spectrum_len";
const NATIVE_KIND: &str = "rssn/native/fft/irfft-plan/v1";

fn invoke(
    s: &Session,
    name: &str,
    args: &[Value],
) -> Result<Value, OpError> {
    s.call_value(name, args, &CallCtx::new())
}

fn call(
    s: &Session,
    name: &str,
    args: &[Value],
) -> Value {
    invoke(s, name, args).unwrap_or_else(|error| panic!("{name}: {error:?}"))
}

fn plan(
    s: &Session,
    n: usize,
) -> Value {
    call(s, NEW, &[Value::int(n)])
}

fn complex(
    re: f64,
    im: f64,
) -> Value {
    Value::Complex(ComplexValue::Float { re, im })
}

fn bins(
    values: Vec<Value>,
    dim: DimExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        data: values,
        elem: TypeExpr::Complex.with_dim(dim),
    })
}

fn check_samples(
    value: &Value,
    expected: &[f64],
    dim: DimExpr,
) {
    let Value::Tensor(tensor) = value else {
        panic!("expected a measured sample tensor")
    };
    assert_eq!(tensor.shape, [expected.len()]);
    assert_eq!(tensor.elem, TypeExpr::Float.with_dim(dim));
    assert_eq!(tensor.data.len(), expected.len());
    for (actual, &expected) in tensor.data.iter().zip(expected) {
        let Value::Float(actual) = actual else {
            panic!("expected a real sample")
        };
        assert!((actual - expected).abs() < 2e-14, "{actual} vs {expected}");
    }
}

fn retained_bytes(s: &Session) -> f64 {
    let Value::Record(usage) = call(s, "num_native_usage", &[]) else {
        panic!("expected a native usage record");
    };
    usage.get("bytes").and_then(Value::as_f64).unwrap()
}

#[test]
fn inverse_uses_the_declared_real_length_and_preserves_each_input_dimension() {
    let s = Session::new();
    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let even = plan(&s, 4);
    let odd = plan(&s, 5);
    let even_bins = [complex(1.0, 17.0), complex(0.0, -1.0), complex(-1.0, -22.0)];
    let even_expected = [0.0, 1.0, 0.0, 0.0];
    check_samples(
        &call(
            &s,
            APPLY,
            &[even.clone(), bins(even_bins.to_vec(), length.clone())],
        ),
        &even_expected,
        length,
    );
    check_samples(
        &call(&s, APPLY, &[even, bins(even_bins.to_vec(), time.clone())]),
        &even_expected,
        time.clone(),
    );

    let odd_values = vec![complex(2.5, 8.0), complex(0.0, 0.0), complex(0.0, 2.5)];
    let odd_expected: Vec<f64> = (0..5)
        .map(|j| 0.5 - (TAU * 2.0 * f64::from(j) / 5.0).sin())
        .collect();
    check_samples(
        &call(&s, APPLY, &[odd, bins(odd_values, time.clone())]),
        &odd_expected,
        time,
    );
}

#[test]
fn complex_quantities_decode_to_coherent_units_and_tiny_samples_survive() {
    let s = Session::new();
    let length = DimExpr::base("Length");
    let centimeter = UnitValue {
        name: "cm".into(),
        display: "cm".into(),
        dimension: length.clone(),
        scale: ScaleValue {
            numer: 1.into(),
            denom: 100.into(),
            pi: 0,
        },
        offset: None,
        affine: false,
        aliases: vec![],
        prefixable: false,
        category: "length".into(),
        since: "beta".into(),
    };
    let measured = (0..3)
        .map(|_| {
            Value::Quantity(QuantityValue {
                value: Box::new(complex(200.0, 0.0)),
                unit: centimeter.clone(),
            })
        })
        .collect();
    check_samples(
        &call(&s, APPLY, &[plan(&s, 4), bins(measured, length.clone())]),
        &[2.0, 0.0, 0.0, 0.0],
        length.clone(),
    );

    let tiny = f64::from_bits(1);
    let expected = Value::Tensor(TensorValue {
        shape: vec![2],
        data: vec![Value::Float(tiny), Value::Float(0.0)],
        elem: TypeExpr::Float.with_dim(length.clone()),
    });
    assert_eq!(
        call(
            &s,
            APPLY,
            &[plan(&s, 2), bins(vec![complex(tiny, 0.0); 2], length)]
        ),
        expected,
    );
}

#[test]
fn rejects_bad_shape_mixed_units_nonfinite_and_foreign_or_released_plans() {
    let s = Session::new();
    let foreign = Session::new();
    let length = DimExpr::base("Length");
    let two = plan(&s, 2);
    let good = bins(vec![complex(1.0, 0.0); 2], length.clone());
    let rank_two = Value::Tensor(TensorValue {
        shape: vec![1, 2],
        data: vec![complex(1.0, 0.0); 2],
        elem: TypeExpr::Complex.with_dim(length.clone()),
    });
    let mixed = bins(
        vec![
            complex(1.0, 0.0),
            Value::Dyn(DynValue {
                value: Box::new(complex(1.0, 0.0)),
                dim: DimExpr::base("Time"),
            }),
        ],
        length.clone(),
    );
    for invalid in [
        rank_two,
        mixed,
        bins(vec![complex(1.0, 0.0)], length.clone()),
        bins(vec![complex(1.0, f64::NAN); 2], length),
    ] {
        assert!(invoke(&s, APPLY, &[two.clone(), invalid]).is_err());
    }
    assert!(invoke(&foreign, APPLY, &[two.clone(), good.clone()]).is_err());
    assert!(invoke(&s, APPLY, &[two.clone(), good.clone()]).is_ok());
    call(&s, "num_native_release", std::slice::from_ref(&two));
    assert!(invoke(&s, APPLY, &[two, good]).is_err());
}

#[test]
fn holds_no_sample_objects_and_checks_length_before_planning() {
    let s = Session::new();
    let dimension = DimExpr::base("Temperature");
    let empty = plan(&s, 0);
    let empty_bins = bins(vec![], dimension.clone());
    let before = retained_bytes(&s);
    check_samples(
        &call(&s, APPLY, &[empty.clone(), empty_bins]),
        &[],
        dimension.clone(),
    );
    assert_eq!(retained_bytes(&s), before);
    assert!(
        invoke(
            &s,
            APPLY,
            &[empty, bins(vec![complex(1.0, 0.0)], dimension)]
        )
        .is_err()
    );
    for invalid in [Value::int(-1), Value::int(262_144), Value::int(usize::MAX)] {
        assert!(invoke(&s, NEW, &[invalid]).is_err());
    }
    assert_eq!(retained_bytes(&s), before);

    let limit = plan(&s, 262_143);
    call(&s, "num_native_release", &[limit]);
    let large = plan(&s, 131_072);
    let retained = retained_bytes(&s) - before;
    assert!((1_048_576.0..=1_049_000.0).contains(&retained));
    call(&s, "num_native_release", &[large]);
    assert_eq!(retained_bytes(&s), before);
}

#[test]
fn reports_real_and_half_spectrum_lengths_for_live_plans() {
    let s = Session::new();
    let other = Session::new();
    for (n, spectrum_len) in [(0, 0), (1, 1), (4, 3), (5, 3)] {
        let plan = plan(&s, n);
        let args = std::slice::from_ref(&plan);
        assert_eq!(call(&s, OUTPUT_LEN, args), Value::int(n));
        assert_eq!(call(&s, SPECTRUM_LEN, args), Value::int(spectrum_len));
        assert!(invoke(&other, OUTPUT_LEN, args).is_err());
        call(&s, "num_native_release", args);
        assert!(invoke(&s, SPECTRUM_LEN, args).is_err());
    }
}

#[test]
fn manifest_maps_the_checked_plan_methods_and_binds_one_shared_dimension() {
    let s = Session::new();
    let manifest = s.manifest();
    let new = manifest.op(NEW).unwrap();
    let apply = manifest.op(APPLY).unwrap();
    for op in [new, apply] {
        assert_eq!(op.signature_source, SignatureSource::Declared);
        assert_eq!(op.overloads.len(), 1);
    }
    let opaque = TypeExpr::Opaque(NATIVE_KIND.into());
    for (name, path) in [
        (
            OUTPUT_LEN,
            "crate::kernels::transforms::IrfftPlan::output_len",
        ),
        (
            SPECTRUM_LEN,
            "crate::kernels::transforms::IrfftPlan::spectrum_len",
        ),
    ] {
        let op = manifest.op(name).unwrap();
        assert_eq!(op.signature_source, SignatureSource::Declared);
        assert_eq!(op.rust_path.as_deref(), Some(path));
        assert_eq!(op.overloads.len(), 1);
        assert_eq!(
            op.overloads[0].params.as_slice(),
            std::slice::from_ref(&opaque)
        );
        assert_eq!(op.overloads[0].ret, TypeExpr::Int);
    }
    assert_eq!(
        new.rust_path.as_deref(),
        Some("crate::kernels::transforms::IrfftPlan::new")
    );
    assert_eq!(
        apply.rust_path.as_deref(),
        Some("crate::kernels::transforms::IrfftPlan::apply_into")
    );
    assert_eq!(new.overloads[0].params, [TypeExpr::Int]);
    assert_eq!(new.overloads[0].ret, opaque);
    let dimension = DimExpr::Var("D".into());
    let vector = |elem: TypeExpr| {
        TypeExpr::Tensor {
            rank: NatExpr::Lit(1),
            shape: None,
            elem: Box::new(elem.with_dim(dimension.clone())),
        }
    };
    let signature = &apply.overloads[0];
    assert_eq!(signature.params, [opaque, vector(TypeExpr::Complex)]);
    assert_eq!(signature.ret, vector(TypeExpr::Float));
    assert_eq!(signature.type_params.len(), 1);
    assert_eq!(signature.type_params[0].name, "D");
    assert!(apply.doc.contains("caller-owned"));
}
