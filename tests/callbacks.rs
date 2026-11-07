//! End-to-end higher-order dispatch, signature evidence and serial quadrature.
use rssn::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::HostCallback;
use rssn::kernel::Kernel;
use rssn::kernel::NativeCallback;
use rssn::kernel::OpError;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::Signature;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::QuantityValue;
use rssn::kernel::value::RecordValue;
use rssn::kernel::value::ScaleValue;
use rssn::kernel::value::TensorValue;
use rssn::kernel::value::UnitValue;
use rssn::kernel::value::Value;
use std::sync::Arc;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;

fn number(value: &Value) -> f64 {
    match value {
        | Value::Dyn(v) => number(&v.value),
        | Value::Quantity(q) => number(&q.value),
        | _ => value.as_f64().unwrap(),
    }
}
fn field<'a>(
    value: &'a Value,
    name: &str,
) -> &'a Value {
    let Value::Record(record) = value else {
        panic!("record")
    };
    &record.fields.iter().find(|(n, _)| n == name).unwrap().1
}
fn args() -> Vec<Value> {
    vec![
        Value::Float(0.0),
        Value::Float(1.0),
        Value::Float(1e-10),
        Value::Float(1e-10),
        Value::int(64),
    ]
}
fn sig() -> Signature {
    Signature::new(vec![TypeExpr::Float], TypeExpr::Float)
}
fn dimensioned(
    value: f64,
    dim: DimExpr,
) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(Value::Float(value)),
        dim,
    })
}

#[test]
fn trailing_callback_handles_reach_the_same_name_and_pinned_operator() {
    let session = Session::new();
    let mut context = CallCtx::new();
    let unused_calls = Arc::new(AtomicUsize::new(0));
    let count = Arc::clone(&unused_calls);
    context
        .add_callback(Arc::new(NativeCallback::new(
            move |_: &[Value]| {
                count.fetch_add(1, Ordering::SeqCst);
                Ok(Value::String("wrong callback".into()))
            },
            sig(),
        )))
        .unwrap();
    let callback = context
        .add_callback(Arc::new(NativeCallback::new(
            |x: &[Value]| Ok(Value::Float(number(&x[0]).powi(2))),
            sig(),
        )))
        .unwrap();
    let ordinary = args();
    let mut combined = ordinary.clone();
    combined.push(Value::Function(callback));
    let key = session
        .resolve_overload_with_context("num_integrate", &combined, &context)
        .unwrap();
    let results = [
        session
            .call_ho("num_integrate", &ordinary, &[callback], &context)
            .unwrap(),
        session.call_key(&key, &combined, &context).unwrap(),
        Kernel::call_ho(&session, &key, &ordinary, &[callback], &context).unwrap(),
    ];
    for result in results {
        assert!((number(field(&result, "value")) - 1.0 / 3.0).abs() < 1e-12);
        assert_eq!(field(&result, "converged"), &Value::Bool(true));
    }
    assert_eq!(unused_calls.load(Ordering::SeqCst), 0);
    assert!(
        session
            .call_ho("num_integrate", &ordinary, &[], &context)
            .is_err()
    );
}

struct ChangingSignature(Arc<AtomicUsize>);
impl HostCallback for ChangingSignature {
    fn signature(&self) -> Signature {
        if self.0.fetch_add(1, Ordering::SeqCst) == 0 {
            sig()
        } else {
            Signature::new(vec![TypeExpr::String], TypeExpr::String)
        }
    }

    fn call(
        &self,
        args: &[Value],
    ) -> Result<Value, OpError> {
        Ok(Value::Float(number(&args[0]) + 1.0))
    }
}
#[test]
fn signatures_are_snapshotted_once_and_unscoped_handles_have_no_fabricated_type() {
    let queries = Arc::new(AtomicUsize::new(0));
    let mut context = CallCtx::new();
    let handle = context
        .add_callback(Arc::new(ChangingSignature(Arc::clone(&queries))))
        .unwrap();
    assert!(!matches!(
        Value::Function(handle).type_expr(),
        TypeExpr::Function(..)
    ));
    assert!(matches!(
        context.value_type(&Value::Function(handle)).unwrap(),
        TypeExpr::Function(..)
    ));
    for _ in 0..3 {
        assert_eq!(
            context
                .invoke_callback(handle, &[Value::Float(2.0)])
                .unwrap(),
            Value::Float(3.0)
        );
    }
    assert_eq!(queries.load(Ordering::SeqCst), 1);
    assert!(CallCtx::new().value_type(&Value::Function(handle)).is_err());
    let generic = Signature::new(vec![TypeExpr::Var("A".into())], TypeExpr::Var("A".into()));
    assert!(
        context
            .add_callback(Arc::new(NativeCallback::new(
                |_: &[Value]| Ok(Value::Float(1.0)),
                generic
            )))
            .is_err()
    );
}

#[test]
fn callback_function_matching_checks_arity_and_input_contravariance() {
    use rssn::kernel::overload::matches_type;
    let f = |a, b| TypeExpr::Function(vec![a], Box::new(b));
    assert!(
        matches_type(
            &f(TypeExpr::Float, TypeExpr::Real),
            &f(TypeExpr::Real, TypeExpr::Float)
        )
        .accepted()
    );
    assert!(
        !matches_type(
            &f(TypeExpr::Real, TypeExpr::Float),
            &f(TypeExpr::Float, TypeExpr::Real)
        )
        .accepted()
    );
    assert!(
        !matches_type(
            &TypeExpr::Function(vec![], Box::new(TypeExpr::Real)),
            &f(TypeExpr::Real, TypeExpr::Real)
        )
        .accepted()
    );
}

#[test]
fn wrong_callback_input_is_rejected_before_host_execution() {
    let calls = Arc::new(AtomicUsize::new(0));
    let counted = Arc::clone(&calls);
    let mut context = CallCtx::new();
    let signature = Signature::new(vec![TypeExpr::String], TypeExpr::Float);
    let handle = context
        .add_callback(Arc::new(NativeCallback::new(
            move |_: &[Value]| {
                counted.fetch_add(1, Ordering::SeqCst);
                Ok(Value::Float(1.0))
            },
            signature,
        )))
        .unwrap();
    assert!(
        Session::new()
            .call_ho("num_integrate", &args(), &[handle], &context)
            .is_err()
    );
    assert_eq!(calls.load(Ordering::SeqCst), 0);
    assert!(matches!(
        context.invoke_callback(handle, &[Value::Float(1.0)]),
        Err(OpError::CallbackType {
            argument: Some(0),
            ..
        })
    ));
}

#[test]
fn wrong_results_and_original_host_diagnostics_are_not_replaced_by_nan() {
    let mut context = CallCtx::new();
    let bad = context
        .add_callback(Arc::new(NativeCallback::new(
            |_: &[Value]| Ok(Value::String("bad".into())),
            sig(),
        )))
        .unwrap();
    assert!(matches!(
        Session::new().call_ho("num_integrate", &args(), &[bad], &context),
        Err(OpError::CallbackType { argument: None, .. })
    ));
    let calls = Arc::new(AtomicUsize::new(0));
    let counted = Arc::clone(&calls);
    let error = context
        .add_callback(Arc::new(NativeCallback::new(
            move |_: &[Value]| {
                counted.fetch_add(1, Ordering::SeqCst);
                Err(OpError::HostDiagnostic {
                    code: "E0805",
                    message: "callback domain".into(),
                })
            },
            sig(),
        )))
        .unwrap();
    assert_eq!(
        Session::new()
            .call_ho("num_integrate", &args(), &[error], &context)
            .unwrap_err(),
        OpError::HostDiagnostic {
            code: "E0805",
            message: "callback domain".into()
        }
    );
    assert_eq!(calls.load(Ordering::SeqCst), 1);
}

#[test]
fn callback_tensor_results_validate_data_lengths_elements_and_duplicate_fields() {
    let mut context = CallCtx::new();
    let tensor_type = TypeExpr::Tensor {
        rank: 1_u64.into(),
        shape: None,
        elem: Box::new(TypeExpr::Float),
    };
    for value in [
        Value::Tensor(TensorValue {
            shape: vec![3],
            data: vec![Value::Float(1.0)],
            elem: TypeExpr::Float,
        }),
        Value::Tensor(TensorValue {
            shape: vec![1],
            data: vec![Value::String("bad".into())],
            elem: TypeExpr::Float,
        }),
        Value::Tensor(TensorValue {
            shape: vec![usize::MAX, 2],
            data: vec![],
            elem: TypeExpr::Float,
        }),
    ] {
        let handle = context
            .add_callback(Arc::new(NativeCallback::new(
                move |_: &[Value]| Ok(value.clone()),
                Signature::new(vec![], tensor_type.clone()),
            )))
            .unwrap();
        assert!(matches!(
            context.invoke_callback(handle, &[]),
            Err(OpError::CallbackType { argument: None, .. })
        ));
    }
    let value = Value::Record(RecordValue {
        fields: vec![
            ("a".into(), Value::Float(1.0)),
            ("a".into(), Value::Float(2.0)),
        ],
    });
    let handle = context
        .add_callback(Arc::new(NativeCallback::new(
            move |_: &[Value]| Ok(value.clone()),
            Signature::new(
                vec![],
                TypeExpr::Record(vec![("a".into(), TypeExpr::Float)]),
            ),
        )))
        .unwrap();
    assert!(context.invoke_callback(handle, &[]).is_err());
}

#[test]
fn dimensional_bounds_use_coherent_units_and_integral_carries_domain_times_codomain() {
    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let output = DimExpr::Mul(vec![length.clone(), time.clone()]);
    let mut context = CallCtx::new();
    let return_dim = time.clone();
    let callback = context
        .add_callback(Arc::new(NativeCallback::new(
            move |values: &[Value]| Ok(dimensioned(number(&values[0]), return_dim.clone())),
            Signature::new(
                vec![TypeExpr::Float.with_dim(length.clone())],
                TypeExpr::Float.with_dim(time),
            ),
        )))
        .unwrap();
    let kilometers = UnitValue {
        name: "km".into(),
        display: "km".into(),
        dimension: length.clone(),
        scale: ScaleValue {
            numer: 1000.into(),
            denom: 1.into(),
            pi: 0,
        },
        offset: None,
        affine: false,
        aliases: vec![],
        prefixable: false,
        category: "test".into(),
        since: "test".into(),
    };
    let ordinary = vec![
        Value::Quantity(QuantityValue {
            value: Box::new(Value::Float(1.0)),
            unit: kilometers,
        }),
        dimensioned(2000.0, length),
        dimensioned(1e-6, output.clone()),
        Value::Float(1e-12),
        Value::int(64),
    ];
    let result = Session::new()
        .call_ho("num_integrate", &ordinary, &[callback], &context)
        .unwrap();
    assert!((number(field(&result, "value")) - 1.5e6).abs() < 1e-7);
    let Value::Dyn(value) = field(&result, "value") else {
        panic!("dimensioned result")
    };
    assert_eq!(value.dim, output);
    let mut wrong = ordinary.clone();
    wrong[2] = Value::Float(1.0);
    assert!(matches!(
        Session::new().call_ho("num_integrate", &wrong, &[callback], &context),
        Err(OpError::Type { arg: 2, .. })
    ));
    wrong[2] = Value::Float(0.0);
    assert!(
        Session::new()
            .call_ho("num_integrate", &wrong, &[callback], &context)
            .is_ok()
    );
}

#[test]
fn zero_and_reversed_intervals_and_budget_status_are_explicit() {
    let calls = Arc::new(AtomicUsize::new(0));
    let counted = Arc::clone(&calls);
    let mut context = CallCtx::new();
    let callback = context
        .add_callback(Arc::new(NativeCallback::new(
            move |values: &[Value]| {
                counted.fetch_add(1, Ordering::SeqCst);
                Ok(Value::Float(number(&values[0]).powi(2)))
            },
            sig(),
        )))
        .unwrap();
    let mut ordinary = args();
    ordinary[1] = Value::Float(0.0);
    let result = Session::new()
        .call_ho("num_integrate", &ordinary, &[callback], &context)
        .unwrap();
    assert_eq!(number(field(&result, "value")), 0.0);
    assert_eq!(calls.load(Ordering::SeqCst), 0);
    ordinary[0] = Value::Float(1.0);
    ordinary[1] = Value::Float(0.0);
    let result = Session::new()
        .call_ho("num_integrate", &ordinary, &[callback], &context)
        .unwrap();
    assert!((number(field(&result, "value")) + 1.0 / 3.0).abs() < 1e-12);
    let oscillatory = context
        .add_callback(Arc::new(NativeCallback::new(
            |v: &[Value]| Ok(Value::Float((1000.0 * number(&v[0])).sin())),
            sig(),
        )))
        .unwrap();
    let mut ordinary = args();
    ordinary[2] = Value::Float(1e-14);
    ordinary[3] = Value::Float(1e-14);
    ordinary[4] = Value::int(1);
    let result = Session::new()
        .call_ho("num_integrate", &ordinary, &[oscillatory], &context)
        .unwrap();
    assert_eq!(field(&result, "converged"), &Value::Bool(false));
}

thread_local! {static NESTED_SESSION:Session=Session::new();}
#[test]
fn callback_can_reenter_the_same_session_without_a_graph_borrow() {
    NESTED_SESSION.with(|session| {
        let mut context = CallCtx::new();
        let callback = context
            .add_callback(Arc::new(NativeCallback::new(
                |_: &[Value]| {
                    NESTED_SESSION.with(|session| {
                        let result = session.call_value(
                            "num_factorial",
                            &[Value::int(3)],
                            &CallCtx::new(),
                        )?;
                        Ok(Value::Float(number(&result)))
                    })
                },
                sig(),
            )))
            .unwrap();
        let result = session
            .call_ho("num_integrate", &args(), &[callback], &context)
            .unwrap();
        assert!((number(field(&result, "value")) - 6.0).abs() < 1e-12);
    });
}

#[test]
fn manifest_provenance_effects_api_version_and_binary_layout_are_explicit() {
    use rssn::kernel::KernelManifest;
    use rssn::kernel::OpInfo;
    use rssn::kernel::SignatureSource;
    use rssn::kernel::signature::Determinism;
    use rssn::kernel::signature::EvalKind;
    use rssn::kernel::signature::Purity;
    use rssn::kernel::signature::ThreadSafety;
    let manifest = Session::new().manifest();
    assert_eq!(manifest.api_version, 2);
    let op = manifest
        .ops
        .iter()
        .find(|op| op.name == "num_integrate")
        .unwrap();
    assert_eq!(
        op.rust_path.as_deref(),
        Some("crate::kernels::quadrature::integrate_gk")
    );
    assert_eq!(op.eval, EvalKind::HostCallback);
    assert_eq!(op.purity, Purity::Effectful);
    assert_eq!(op.determinism, Determinism::EnvDependent);
    assert_eq!(op.thread_safety, ThreadSafety::SessionBound);
    assert_eq!(op.signature_source, SignatureSource::Declared);
    assert!(
        op.higher_order
            && op.backend_support.interpreter
            && !op.backend_support.ffi
            && !op.backend_support.jit
    );
    let bytes = manifest.to_bytes().unwrap();
    assert_eq!(KernelManifest::from_bytes(&bytes).unwrap(), *manifest);
    let mut legacy = serde_json::to_value(op).unwrap();
    legacy.as_object_mut().unwrap().remove("rust_path");
    let legacy: OpInfo = serde_json::from_value(legacy).unwrap();
    assert_eq!(legacy.rust_path, None);
}

#[test]
fn api_v2_manifest_crosses_the_c_string_boundary() {
    let pointer = rssn::ffi::rssn_manifest_json();
    assert!(!pointer.is_null());
    // SAFETY: this pointer is the live owned CString returned by RSSN.
    let text = unsafe { std::ffi::CStr::from_ptr(pointer) }
        .to_str()
        .unwrap()
        .to_owned();
    // SAFETY: the returned CString is released exactly once.
    unsafe { rssn::ffi::rssn_string_free(pointer) };
    let manifest: serde_json::Value = serde_json::from_str(&text).unwrap();
    assert_eq!(manifest["api_version"].as_u64(), Some(2));
    assert!(manifest["ops"].as_array().unwrap().iter().any(|op| {
        op["name"] == "num_integrate"
            && op["rust_path"] == "crate::kernels::quadrature::integrate_gk"
    }));
}

#[test]
fn equivalent_concrete_dimensions_work_for_callback_name_and_key_dispatch() {
    use rssn::kernel::signature::Exp;
    use rssn::kernel::signature::Rational;
    let session = Session::new();
    let length = DimExpr::Base("Length".into());
    let square = DimExpr::Pow(Box::new(length.clone()), Exp::Rat(Rational::new(2, 1)));
    let product = DimExpr::Mul(vec![length.clone(), length.clone()]);
    let time = DimExpr::Base("Time".into());
    let mut context = CallCtx::new();
    let callback = context
        .add_callback(Arc::new(NativeCallback::new(
            {
                let time = time.clone();
                move |_: &[Value]| Ok(dimensioned(1.0, time.clone()))
            },
            Signature::new(
                vec![TypeExpr::Dim(Box::new(TypeExpr::Float), square.clone())],
                TypeExpr::Dim(Box::new(TypeExpr::Float), time.clone()),
            ),
        )))
        .unwrap();
    let mut ordinary = args();
    ordinary[0] = dimensioned(0.0, product.clone());
    ordinary[1] = dimensioned(1.0, product.clone());
    ordinary[2] = Value::Float(0.0);
    let mut combined = ordinary.clone();
    combined.push(Value::Function(callback));
    let key = session
        .resolve_overload_with_context("num_integrate", &combined, &context)
        .unwrap();
    for result in [
        session.call_ho("num_integrate", &ordinary, &[callback], &context),
        session.call_key(&key, &combined, &context),
    ] {
        let result = result.unwrap();
        assert!((number(field(&result, "value")) - 1.0).abs() < 1e-12);
    }
    ordinary[1] = dimensioned(1.0, time.clone());
    assert!(
        session
            .call_ho("num_integrate", &ordinary, &[callback], &context)
            .is_err()
    );
    combined[1] = dimensioned(1.0, time);
    assert!(session.call_key(&key, &combined, &context).is_err());
}

#[test]
fn callback_tensor_metadata_uses_algebraic_dimension_equality() {
    use rssn::kernel::signature::Exp;
    use rssn::kernel::signature::NatExpr;
    use rssn::kernel::signature::Rational;
    let length = DimExpr::Base("Length".into());
    let square = DimExpr::Pow(Box::new(length.clone()), Exp::Rat(Rational::new(2, 1)));
    let product = DimExpr::Mul(vec![length.clone(), length]);
    let tensor_ty = |dim| {
        TypeExpr::Tensor {
            rank: NatExpr::Lit(1),
            shape: None,
            elem: Box::new(TypeExpr::Dim(Box::new(TypeExpr::Float), dim)),
        }
    };
    let mut context = CallCtx::new();
    let callback = context
        .add_callback(Arc::new(NativeCallback::new(
            |values: &[Value]| Ok(values[0].clone()),
            Signature::new(vec![tensor_ty(square.clone())], tensor_ty(square)),
        )))
        .unwrap();
    let tensor = Value::Tensor(TensorValue {
        shape: vec![1],
        data: vec![Value::Float(2.0)],
        elem: TypeExpr::Dim(Box::new(TypeExpr::Float), product),
    });
    assert_eq!(
        context
            .invoke_callback(callback, std::slice::from_ref(&tensor))
            .unwrap(),
        tensor
    );
}

#[test]
fn callback_schema_checks_directional_shapes_and_nested_records() {
    use rssn::kernel::overload::matches_type;
    use rssn::kernel::signature::NatExpr;
    let tensor = |shape: Option<Vec<NatExpr>>| {
        TypeExpr::Tensor {
            rank: NatExpr::Lit(1),
            shape,
            elem: Box::new(TypeExpr::Float),
        }
    };
    let fixed = tensor(Some(vec![NatExpr::Lit(3)]));
    let other = tensor(Some(vec![NatExpr::Lit(4)]));
    let open = tensor(None);
    let variable = tensor(Some(vec![NatExpr::Var("N".into())]));
    let fun = |input, output| TypeExpr::Function(vec![input], Box::new(output));
    let accepts =
        |expected: &TypeExpr, actual: &TypeExpr| matches_type(expected, actual).accepted();
    assert!(!accepts(
        &fun(open.clone(), fixed.clone()),
        &fun(fixed.clone(), fixed.clone())
    ));
    assert!(accepts(
        &fun(fixed.clone(), fixed.clone()),
        &fun(open.clone(), fixed.clone())
    ));
    assert!(!accepts(
        &fun(fixed.clone(), fixed.clone()),
        &fun(open.clone(), open.clone())
    ));
    assert!(accepts(
        &fun(variable.clone(), variable.clone()),
        &fun(fixed.clone(), fixed.clone())
    ));
    assert!(!accepts(
        &fun(fixed.clone(), fixed.clone()),
        &fun(other, fixed.clone())
    ));
    let record = |v| TypeExpr::Record(vec![("state".into(), v)]);
    assert!(accepts(
        &fun(record(fixed.clone()), record(fixed.clone())),
        &fun(record(open), record(fixed))
    ));
    let duplicate = TypeExpr::Record(vec![
        ("x".into(), TypeExpr::Float),
        ("x".into(), TypeExpr::Float),
    ]);
    assert!(!accepts(&duplicate, &duplicate));
    assert!(!accepts(
        &fun(TypeExpr::Float, TypeExpr::Float),
        &fun(
            TypeExpr::Dyn(rssn::kernel::signature::RepExpr {
                repr: "unknown".into(),
                dim: None
            }),
            TypeExpr::Float
        )
    ));
}

#[test]
fn callback_registration_rejects_contradictory_tensor_schema() {
    use rssn::kernel::signature::NatExpr;
    let invalid = TypeExpr::Tensor {
        rank: NatExpr::Lit(1),
        shape: Some(vec![NatExpr::Lit(2), NatExpr::Lit(3)]),
        elem: Box::new(TypeExpr::Float),
    };
    let mut context = CallCtx::new();
    let result = context.add_callback(Arc::new(NativeCallback::new(
        |_: &[Value]| Ok(Value::Float(0.0)),
        Signature::new(vec![invalid], TypeExpr::Float),
    )));
    assert!(result.is_err());
}
