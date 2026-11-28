//! Branch-cut evaluation through the checked value ABI and symbolic graph.

use std::collections::HashMap;
use std::f64::consts::PI;
use std::f64::consts::TAU;

use num_complex::Complex64;
use rssn::api::Session;
use rssn::graph::Engine;
use rssn::graph::Graph;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::manifest::SignatureSource;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::DimRuleInfo;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::ComplexValue;
use rssn::kernel::value::DynValue;
use rssn::rules::complex::complex;

const ARG: &str = "num_complex_coords_arg_branch";
const LOG: &str = "num_complex_coords_log_cut";
const SQRT: &str = "num_complex_coords_sqrt_cut";
const POW: &str = "num_complex_coords_pow_cut";

const fn value(z: Complex64) -> Value {
    Value::Complex(ComplexValue::Float { re: z.re, im: z.im })
}

fn as_complex(value: &Value) -> Complex64 {
    let Value::Complex(number) = value else {
        panic!("expected a complex result, got {value:?}");
    };
    let (real, imaginary) = number.to_f64();
    Complex64::new(real, imaginary)
}

#[test]
fn checked_values_select_the_upper_endpoint_of_the_cut() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let theta = Value::Float(0.0);
    for imaginary in [0.0, -0.0] {
        let z = value(Complex64::new(4.0, imaginary));
        let argument = session
            .call_value(ARG, &[z.clone(), theta.clone()], &ctx)
            .unwrap();
        assert_eq!(argument.as_f64(), Some(TAU));
        let log = session
            .call_value(LOG, &[z.clone(), theta.clone()], &ctx)
            .unwrap();
        assert_eq!(as_complex(&log), Complex64::new(4.0_f64.ln(), TAU));
        for (name, args) in [
            (SQRT, vec![z.clone(), theta.clone()]),
            (POW, vec![z, value(Complex64::new(0.5, 0.0)), theta.clone()]),
        ] {
            let result = session.call_value(name, &args, &ctx).unwrap();
            assert!((as_complex(&result) + 2.0).norm() < 1e-14, "{name}");
        }
    }
    for (imaginary, upper) in [(1e-12, false), (-1e-12, true)] {
        let z = value(Complex64::new(1.0, imaginary));
        let phase = session.call_value(ARG, &[z, theta.clone()], &ctx).unwrap();
        assert_eq!(phase.as_f64().unwrap() > PI, upper);
    }
}

#[test]
fn checked_values_validate_zero_angles_and_dimensions() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let z = value(Complex64::new(0.0, 0.0));
    let theta = Value::Float(0.0);
    let phase = session
        .call_value(ARG, &[z.clone(), theta.clone()], &ctx)
        .unwrap();
    assert_eq!(phase.as_f64(), Some(TAU));
    let root = session
        .call_value(SQRT, &[z.clone(), theta.clone()], &ctx)
        .unwrap();
    assert_eq!(as_complex(&root), Complex64::new(0.0, 0.0));
    for (name, args) in [
        (LOG, vec![z.clone(), theta.clone()]),
        (POW, vec![z, value(Complex64::new(0.5, 0.0)), theta.clone()]),
    ] {
        assert!(matches!(
            session.call_value(name, &args, &ctx),
            Err(OpError::Domain(_))
        ));
    }

    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let measured = Value::Dyn(DynValue {
        value: Box::new(value(Complex64::new(4.0, 0.0))),
        dim: length.clone().pow(2_i64),
    });
    let root = session
        .call_value(SQRT, &[measured.clone(), theta.clone()], &ctx)
        .unwrap();
    let Value::Dyn(root) = root else {
        panic!("a dimensioned radicand must produce a dimensioned result");
    };
    assert_eq!(root.dim, length);
    assert!((as_complex(&root.value) + 2.0).norm() < 1e-14);

    let measured_angle = Value::Dyn(DynValue {
        value: Box::new(theta),
        dim: time,
    });
    assert!(
        session
            .call_value(SQRT, &[measured.clone(), measured_angle], &ctx)
            .is_err()
    );
    for op in [ARG, LOG] {
        assert!(
            session
                .call_value(op, &[measured.clone(), Value::Float(0.0)], &ctx)
                .is_err()
        );
    }
    assert!(
        session
            .call_value(
                SQRT,
                &[value(Complex64::new(4.0, 0.0)), Value::Float(f64::INFINITY)],
                &ctx,
            )
            .is_err()
    );
}

#[test]
fn symbolic_cut_operators_keep_their_keys_and_use_the_same_branch() {
    let session = Session::new();
    let manifest = session.manifest();
    let ctx = CallCtx::new();
    let z = session.value_of(session.parse("z").unwrap()).unwrap();
    for (name, native, base, arity) in [
        ("arg_branch", ARG, TypeExpr::Float, 2),
        ("log_cut", LOG, TypeExpr::Complex, 2),
        ("sqrt_cut", SQRT, TypeExpr::Complex, 2),
        ("pow_cut", POW, TypeExpr::Complex, 3),
    ] {
        let ordinary = manifest.op(name).unwrap();
        assert_eq!(ordinary.signature_source, SignatureSource::Overlay);
        assert_eq!(ordinary.keys.len(), 1);
        assert!(matches!(&ordinary.overloads[0].ret, TypeExpr::Dyn(_)));
        let mut args = vec![z.clone(), Value::Float(0.0)];
        if arity == 3 {
            args.insert(1, Value::Float(0.5));
        }
        assert!(matches!(
            session.call_key(&ordinary.keys[0], &args, &ctx),
            Ok(Value::Term(_))
        ));

        let checked = manifest.op(native).unwrap();
        assert_eq!(
            checked.rust_path.as_deref(),
            Some(format!("crate::kernels::complex_coords::{name}").as_str())
        );
        assert_eq!(checked.signature_source, SignatureSource::Declared);
        assert_eq!(checked.overloads[0].arity(), arity);
        if native == SQRT {
            assert_eq!(checked.dim_rule, Some(DimRuleInfo::Root(2)));
            assert!(
                matches!(&checked.overloads[0].ret, TypeExpr::Dim(inner, _) if **inner == base)
            );
        } else {
            assert_eq!(checked.dim_rule, Some(DimRuleInfo::Dimensionless));
            assert_eq!(checked.overloads[0].ret, base);
        }
    }

    let mut graph = Graph::new();
    let _engine = Engine::install(&mut graph, &[complex()]).unwrap();
    let z = graph.sym("z");
    let symbol = graph.symbol_of(z).unwrap();
    let bindings = HashMap::from([(symbol, Complex64::new(4.0, -0.0))]);
    for (expression, expected) in [
        ("arg_branch(z,0)", Complex64::new(TAU, 0.0)),
        ("log_cut(z,0)", Complex64::new(4.0_f64.ln(), TAU)),
        ("sqrt_cut(z,0)", Complex64::new(-2.0, 0.0)),
        ("pow_cut(z,1/2,0)", Complex64::new(-2.0, 0.0)),
    ] {
        let term = graph.parse(expression).unwrap();
        let result = graph.eval_complex(term, &bindings).unwrap();
        assert!((result - expected).norm() < 1e-14, "{expression}: {result}");
    }
}
