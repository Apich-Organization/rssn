//! Checked and reusable machine-word modular reducers over the public value ABI.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::manifest::SignatureSource;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::DimRuleInfo;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::VariantValue;

fn call(
    session: &Session,
    name: &str,
    args: &[Value],
) -> Value {
    session
        .call_value(&format!("num_modular_{name}"), args, &CallCtx::new())
        .unwrap_or_else(|error| panic!("{name}: {error:?}"))
}

fn optional(value: Option<u64>) -> Value {
    let (tag, payload) = match value {
        | Some(word) => ("Some", vec![Value::int(word)]),
        | None => ("None", vec![]),
    };
    Value::Variant(VariantValue {
        tag: tag.into(),
        payload,
    })
}

#[test]
fn reducers_keep_full_word_products_at_large_odd_and_even_moduli() {
    let session = Session::new();
    for (kind, moduli) in [
        ("montgomery", [15_u64, 65_537, u64::MAX]),
        ("barrett", [2, 12, u64::MAX - 1]),
    ] {
        for modulus in moduli {
            let plan = call(&session, &format!("{kind}_new"), &[Value::int(modulus)]);
            let modulus_wide = u128::from(modulus);
            let reduce = |x: u128| Value::int(x % modulus_wide);
            for left in [0, 1, 14, u64::MAX / 2, u64::MAX - 2, u64::MAX] {
                for right in [0, 1, 13, u64::MAX - 1] {
                    let a = u128::from(left);
                    let b = u128::from(right);
                    for (operation, expected) in [
                        ("add", reduce(a + b)),
                        (
                            "sub",
                            reduce(modulus_wide + a % modulus_wide - b % modulus_wide),
                        ),
                        ("mul", reduce(a * b)),
                    ] {
                        assert_eq!(
                            call(
                                &session,
                                &format!("{kind}_{operation}"),
                                &[plan.clone(), Value::int(left), Value::int(right)]
                            ),
                            expected,
                            "{kind} {operation} mod {modulus}"
                        );
                    }
                }
            }
            let x = Value::int(u64::MAX - 2);
            let one = call(
                &session,
                &format!("{kind}_pow"),
                &[plan.clone(), x.clone(), Value::int(0)],
            );
            let squared = call(
                &session,
                &format!("{kind}_pow"),
                &[plan.clone(), x, Value::int(2)],
            );
            assert_eq!(one, Value::int(1));
            assert_eq!(squared, reduce(u128::from(u64::MAX - 2).pow(2)));
        }
    }
}

#[test]
fn montgomery_helpers_keep_encoded_and_ordinary_representations_separate() {
    let session = Session::new();
    let modulus = 17_u64;
    let plan = call(&session, "montgomery_new", &[Value::int(modulus)]);
    let with = |name: &str, xs: &[Value]| {
        let args = [vec![plan.clone()], xs.to_vec()].concat();
        call(&session, &format!("montgomery_{name}"), &args)
    };
    let a = with("to_mont", &[Value::int(15)]);
    let b = with("to_mont_u128", &[Value::int(u128::MAX - 1)]);
    let decode = |value| with("from_mont", &[value]);
    assert_eq!(decode(a.clone()), Value::int(15));
    assert_eq!(decode(b.clone()), Value::int(16));
    assert_eq!(
        decode(with("add_mont", &[a.clone(), b.clone()])),
        Value::int(14)
    );
    assert_eq!(
        decode(with("sub_mont", &[a.clone(), b.clone()])),
        Value::int(16)
    );
    assert_eq!(
        decode(with("neg_mont", std::slice::from_ref(&a))),
        Value::int(2)
    );
    assert_eq!(
        decode(with("mul_mont", &[a.clone(), b.clone()])),
        Value::int(2)
    );
    assert_eq!(
        decode(with("fma_mont", &[a.clone(), b, a.clone()])),
        Value::int(0)
    );
    assert_eq!(
        decode(with("pow_mont", &[a.clone(), Value::int(3)])),
        Value::int(9)
    );
    assert_eq!(decode(with("half", &[a])), Value::int(16));
    assert_eq!(with("one", &[]), Value::int(1));
    assert_eq!(with("one_mont", &[]), with("r_mod", &[]));
    assert_eq!(with("r2_mod", &[]), Value::int(1));
    let Value::Int(multiplier) = with("n_prime", &[]) else {
        panic!("integer multiplier")
    };
    let multiplier = u64::try_from(multiplier).unwrap();
    assert_eq!(modulus.wrapping_mul(multiplier), u64::MAX);
    assert_eq!(
        with(
            "redc",
            &[Value::int(u128::from(modulus) * (1_u128 << 64) - 1)]
        ),
        Value::int(16)
    );
}

#[test]
fn composite_moduli_and_wide_reduction_have_total_checked_results() {
    let session = Session::new();
    for (kind, modulus) in [("montgomery", 15_u64), ("barrett", 12)] {
        let plan = call(&session, &format!("{kind}_new"), &[Value::int(modulus)]);
        let with = |name: &str, xs: &[Value]| {
            let args = [vec![plan.clone()], xs.to_vec()].concat();
            call(&session, &format!("{kind}_{name}"), &args)
        };
        assert_eq!(with("modulus", &[]), Value::int(modulus));
        assert_eq!(with("inv", &[Value::int(3)]), optional(None));
        assert_eq!(with("div", &[Value::int(7), Value::int(3)]), optional(None));
        assert_eq!(
            with("inv", &[Value::int(modulus - 1)]),
            optional(Some(modulus - 1))
        );
        assert_eq!(
            with("div", &[Value::int(1), Value::int(modulus - 1)]),
            optional(Some(modulus - 1))
        );
        assert_eq!(
            with("neg", &[Value::int(modulus + 1)]),
            Value::int(modulus - 1)
        );
    }
    let modulus = u64::MAX - 1;
    let plan = call(&session, "barrett_new", &[Value::int(modulus)]);
    let with = |name, value| call(&session, name, &[plan.clone(), Value::int(value)]);
    assert_eq!(
        with("barrett_reduce_u128", u128::MAX),
        Value::int(u128::MAX % u128::from(modulus))
    );
    assert_eq!(
        with("barrett_reduce_u64", u128::from(u64::MAX)),
        Value::int(1)
    );
    assert_eq!(call(&session, "barrett_mu", &[plan]), Value::int(1));
    assert_eq!(
        call(
            &session,
            "barrett_modulus",
            &[call(&session, "barrett_new_prime", &[Value::int(2)])]
        ),
        Value::int(2)
    );
    assert_eq!(
        call(
            &session,
            "montgomery_modulus",
            &[call(&session, "montgomery_new_prime", &[Value::int(97)])]
        ),
        Value::int(97)
    );
}

#[test]
fn invalid_inputs_and_handles_are_rejected_before_the_kernel_runs() {
    let session = Session::new();
    let foreign = Session::new();
    let ctx = CallCtx::new();
    let plan = call(&session, "montgomery_new", &[Value::int(15)]);
    let barrett = call(&session, "barrett_new", &[Value::int(12)]);
    let invalid = [
        ("num_modular_montgomery_new", vec![Value::int(2)]),
        ("num_modular_barrett_new", vec![Value::int(0)]),
        ("num_modular_montgomery_new_prime", vec![Value::int(15)]),
        ("num_modular_barrett_new_prime", vec![Value::int(12)]),
        (
            "num_modular_montgomery_add_mont",
            vec![plan.clone(), Value::int(15), Value::int(1)],
        ),
        (
            "num_modular_montgomery_mul_mont",
            vec![plan.clone(), Value::int(1), Value::int(15)],
        ),
        (
            "num_modular_montgomery_fma_mont",
            vec![plan.clone(), Value::int(1), Value::int(2), Value::int(15)],
        ),
        (
            "num_modular_montgomery_redc",
            vec![plan.clone(), Value::int(15_u128 << 64)],
        ),
        (
            "num_modular_barrett_reduce_u64",
            vec![barrett.clone(), Value::int(-1)],
        ),
        (
            "num_modular_barrett_reduce_u64",
            vec![barrett, Value::int(1_u128 << 64)],
        ),
        (
            "num_modular_montgomery_half",
            vec![plan.clone(), Value::int(15)],
        ),
    ];
    for (name, args) in invalid {
        assert!(
            matches!(
                session.call_value(name, &args, &ctx),
                Err(OpError::Domain(_))
            ),
            "{name}"
        );
    }
    let measured = Value::Dyn(DynValue {
        value: Box::new(Value::int(15)),
        dim: DimExpr::base("Length"),
    });
    assert!(
        session
            .call_value("num_modular_montgomery_new", &[measured], &ctx)
            .is_err()
    );
    assert!(
        session
            .call_value(
                "num_modular_barrett_modulus",
                std::slice::from_ref(&plan),
                &ctx
            )
            .is_err()
    );
    assert!(
        foreign
            .call_value(
                "num_modular_montgomery_modulus",
                std::slice::from_ref(&plan),
                &ctx
            )
            .is_err()
    );
    assert_eq!(
        session
            .call_value("num_native_release", std::slice::from_ref(&plan), &ctx)
            .unwrap(),
        Value::Bool(true)
    );
    assert!(
        session
            .call_value("num_modular_montgomery_modulus", &[plan], &ctx)
            .is_err()
    );
}

#[test]
fn every_reducer_method_has_one_typed_dimensionless_source_descriptor() {
    let manifest = Session::new().manifest();
    let methods: Vec<_> = manifest
        .ops
        .iter()
        .filter(|op| {
            op.name.starts_with("num_modular_montgomery_")
                || op.name.starts_with("num_modular_barrett_")
        })
        .collect();
    assert_eq!(methods.len(), 39);
    for op in methods {
        let rest = op.name.strip_prefix("num_modular_").unwrap();
        let (kind, method) = rest.split_once('_').unwrap();
        let source_kind = if kind == "barrett" {
            "Barrett"
        } else {
            "Montgomery"
        };
        assert_eq!(
            op.rust_path.as_deref(),
            Some(format!("crate::kernels::modular::{source_kind}::{method}").as_str())
        );
        assert_eq!(op.dim_rule, Some(DimRuleInfo::Dimensionless));
        assert_eq!(op.signature_source, SignatureSource::Declared);
        assert_eq!(op.overloads.len(), 1);
        let signature = &op.overloads[0];
        let constructor = method.starts_with("new");
        let handle = TypeExpr::Opaque(format!("rssn/native/modular/{kind}/v1"));
        if constructor {
            assert_eq!(signature.params, [TypeExpr::Int]);
            assert_eq!(signature.ret, handle);
        } else {
            assert_eq!(signature.params.first(), Some(&handle));
            assert!(
                signature
                    .params
                    .iter()
                    .skip(1)
                    .all(|ty| *ty == TypeExpr::Int)
            );
            if matches!(method, "inv" | "div") {
                assert_eq!(
                    signature.ret,
                    TypeExpr::Variant(vec![
                        ("None".into(), vec![]),
                        ("Some".into(), vec![TypeExpr::Int])
                    ])
                );
            } else {
                assert_eq!(signature.ret, TypeExpr::Int);
            }
        }
    }
}
