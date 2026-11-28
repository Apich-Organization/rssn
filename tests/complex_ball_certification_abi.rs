//! Closed-form checks and whole-disk callback contracts for complex certification.

use num_complex::Complex64 as C;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::NativeCallback;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::Signature;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::ComplexValue;
use rssn::kernel::value::RecordValue;
use rssn::kernel::value::TensorValue;
use rssn::kernel::value::VariantValue;
use rssn::kernels::complex_ball::ComplexBall as B;
use std::sync::Arc;

fn tensor(
    shape: Vec<usize>,
    values: impl IntoIterator<Item = f64>,
) -> Value {
    Value::Tensor(TensorValue {
        shape,
        elem: TypeExpr::Float,
        data: values.into_iter().map(Value::Float).collect(),
    })
}

fn raw_disk(ball: B) -> [f64; 3] {
    [ball.mid.re, ball.mid.im, ball.rad]
}

fn disk(ball: B) -> Value {
    tensor(vec![3], raw_disk(ball))
}

fn decode_disk(value: &Value) -> B {
    let Value::Tensor(tensor) = value else {
        panic!("expected disk tensor")
    };
    let [re, im, rad] = tensor.data.as_slice() else {
        panic!("expected three disk coordinates")
    };
    B::new(
        C::new(re.as_f64().unwrap(), im.as_f64().unwrap()),
        rad.as_f64().unwrap(),
    )
    .unwrap()
}

fn disks(value: &Value) -> Vec<B> {
    let Value::Tensor(tensor_value) = value else {
        panic!("expected disk box tensor")
    };
    let (cells, _) = tensor_value.data.as_chunks::<3>();
    cells
        .iter()
        .map(|cell| decode_disk(&tensor(vec![3], cell.iter().map(|v| v.as_f64().unwrap()))))
        .collect()
}

fn field<'a>(
    value: &'a Value,
    name: &str,
) -> &'a Value {
    let Value::Record(record) = value else {
        panic!("expected jet")
    };
    record.get(name).unwrap()
}

fn some(value: Value) -> Value {
    Value::Variant(VariantValue {
        tag: "Some".into(),
        payload: vec![value],
    })
}

fn unwrap_some(value: &Value) -> &Value {
    let Value::Variant(option) = value else {
        panic!("expected certificate option")
    };
    assert_eq!(option.tag, "Some");
    let [certificate] = option.payload.as_slice() else {
        panic!("expected a certificate")
    };
    certificate
}

fn callback_type(
    session: &Session,
    op: &str,
    slot: usize,
) -> Signature {
    let manifest = session.manifest();
    let parameter = &manifest
        .ops
        .iter()
        .find(|s| s.name == op)
        .unwrap()
        .overloads[0]
        .params[slot];
    let TypeExpr::Function(args, result) = parameter else {
        panic!("expected callback signature")
    };
    Signature::new(args.clone(), result.as_ref().clone())
}

#[test]
fn disk_special_functions_cover_known_values_and_branch_boundaries() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let origin = disk(B::new(C::new(0.0, 0.0), 0.01).unwrap());
    let erf = session
        .call_value("num_complex_disk_erf", &[origin], &ctx)
        .unwrap();
    assert!(decode_disk(unwrap_some(&erf)).contains(C::new(0.0, 0.011_284_167_808_628_218)));
    let point = disk(B::point(1.0, 0.0));
    let j0 = session
        .call_value(
            "num_complex_disk_bessel_j",
            &[point, Value::Float(0.0)],
            &ctx,
        )
        .unwrap();
    assert!(decode_disk(unwrap_some(&j0)).contains(C::new(0.765_197_686_557_966_6, 0.0)));
    let near_three = disk(B::new(C::new(3.0, 0.0), 1e-4).unwrap());
    let gamma = session
        .call_value("num_complex_disk_gamma", &[near_three], &ctx)
        .unwrap();
    assert!(decode_disk(unwrap_some(&gamma)).contains(C::new(2.0, 0.0)));
    let log_gamma = session
        .call_value(
            "num_complex_disk_log_gamma",
            &[disk(B::point(3.0, 0.0))],
            &ctx,
        )
        .unwrap();
    assert!(decode_disk(unwrap_some(&log_gamma)).contains(C::new(std::f64::consts::LN_2, 0.0)));
    let cut = disk(B::new(C::new(-2.0, 0.0), 1e-4).unwrap());
    for (name, args) in [
        ("log_gamma", vec![cut.clone()]),
        ("bessel_j", vec![cut, Value::Float(0.5)]),
    ] {
        let result = session
            .call_value(&format!("num_complex_disk_{name}"), &args, &ctx)
            .unwrap();
        assert!(matches!(result,Value::Variant(option) if option.tag=="None"));
    }
    let unbounded = disk(B::new(C::new(0.0, 0.0), f64::INFINITY).unwrap());
    assert_eq!(
        session
            .call_value("num_complex_disk_is_finite", &[unbounded], &ctx)
            .unwrap(),
        Value::Bool(false)
    );
}

#[test]
fn scalar_krawczyk_accepts_a_jet_enclosing_a_complete_disk() {
    let session = Session::new();
    let mut ctx = CallCtx::new();
    let name = "num_complex_disk_certify_root";
    let signature = callback_type(&session, name, 2);
    let callback = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |args: &[Value]| {
                let x = decode_disk(field(&args[0], "value"));
                let dx = decode_disk(field(&args[0], "derivative"));
                Ok(some(Value::Record(RecordValue {
                    fields: vec![
                        ("value".into(), disk(x.sqr() - B::point(2.0, 0.0))),
                        ("derivative".into(), disk(x.scale_real(2.0) * dx)),
                    ],
                })))
            },
            signature,
        )))
        .unwrap();
    let start = B::new(C::new(1.41, 0.0), 0.02).unwrap();
    let result = session
        .call_ho(name, &[disk(start), Value::int(2)], &[callback], &ctx)
        .unwrap();
    let root = decode_disk(unwrap_some(&result));
    assert!(root.contains(C::new(2.0_f64.sqrt(), 0.0)) && root.inside_strictly(&start));
    assert!(matches!(
        session.call_ho(name, &[disk(start), Value::int(1)], &[callback], &ctx),
        Err(OpError::BudgetExceeded)
    ));
}

#[test]
fn multidimensional_krawczyk_uses_a_whole_box_jacobian() {
    let session = Session::new();
    let mut ctx = CallCtx::new();
    let name = "num_complex_disk_certify_system";
    let value_callback = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |args: &[Value]| {
                let values = disks(&args[0]);
                let outputs = values
                    .iter()
                    .zip([2.0, 3.0])
                    .map(|(x, target)| x.sqr() - B::point(target, 0.0));
                Ok(some(tensor(vec![2, 3], outputs.flat_map(raw_disk))))
            },
            callback_type(&session, name, 2),
        )))
        .unwrap();
    let jac = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |args: &[Value]| {
                let values = disks(&args[0]);
                Ok(some(tensor(
                    vec![2, 2, 3],
                    [
                        values[0].scale_real(2.0),
                        B::zero(),
                        B::zero(),
                        values[1].scale_real(2.0),
                    ]
                    .into_iter()
                    .flat_map(raw_disk),
                )))
            },
            callback_type(&session, name, 3),
        )))
        .unwrap();
    let start = [
        B::new(C::new(1.41, 0.0), 0.02).unwrap(),
        B::new(C::new(1.73, 0.0), 0.02).unwrap(),
    ];
    let result = session
        .call_ho(
            name,
            &[
                tensor(vec![2, 3], start.into_iter().flat_map(raw_disk)),
                Value::int(3),
            ],
            &[value_callback, jac],
            &ctx,
        )
        .unwrap();
    let root = disks(unwrap_some(&result));
    for ((out, start), expected) in root.iter().zip(start).zip([2.0_f64.sqrt(), 3.0_f64.sqrt()]) {
        assert!(out.contains(C::new(expected, 0.0)) && out.inside_strictly(&start));
    }
}

#[test]
fn rouche_proof_uses_covering_disks_and_fails_closed_on_coarse_samples() {
    let session = Session::new();
    let mut ctx = CallCtx::new();
    let name = "num_complex_disk_certify_rouche";
    let f = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |args: &[Value]| Ok(some(disk(decode_disk(&args[0]).scale_real(3.0)))),
            callback_type(&session, name, 4),
        )))
        .unwrap();
    let g = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |args: &[Value]| {
                Ok(some(disk(
                    decode_disk(&args[0]).powi(5).unwrap() + B::one(),
                )))
            },
            callback_type(&session, name, 5),
        )))
        .unwrap();
    let verify = |segments, budget| {
        session.call_ho(
            name,
            &[
                Value::Complex(ComplexValue::Float { re: 0.0, im: 0.0 }),
                Value::Float(1.0),
                Value::int(segments),
                Value::int(budget),
            ],
            &[f, g],
            &ctx,
        )
    };
    assert_eq!(verify(64, 128).unwrap(), Value::Bool(true));
    assert_eq!(verify(4, 8).unwrap(), Value::Bool(false));
    assert!(matches!(verify(64, 127), Err(OpError::BudgetExceeded)));
}

#[test]
fn closed_disks_cover_complex_high_precision_references() {
    // 80-digit references at interior sample points, rounded to binary64.
    for (center, radius, references) in [
        (
            C::new(0.3, 0.4),
            0.01,
            [
                C::new(0.3899603698167964, 0.43302525766198674),
                C::new(0.9102980398049176, -1.341103147751163),
                C::new(0.4829593532874888, -0.9744625315018476),
                C::new(0.3107061171161256, 0.2391152508150076),
            ],
        ),
        (
            C::new(1.25, -0.75),
            0.001,
            [
                C::new(1.0826540568489909, -0.09491226458197148),
                C::new(0.664623668659425, 0.06009002425536292),
                C::new(-0.4044637505428701, 0.09016695191166137),
                C::new(0.7466061726062971, -0.09276772743606214),
            ],
        ),
        (
            C::new(-1.6, 0.45),
            0.005,
            [
                C::new(-1.0025765007841507, 0.02862004582424221),
                C::new(1.0638221876091796, 0.07352166614760036),
                C::new(0.06425073033109763, -6.214184164694161),
                C::new(-0.48624395495726586, 0.45885581235775047),
            ],
        ),
        (
            C::new(-2.2, -0.6),
            0.001,
            [
                C::new(-1.0025170969572996, -0.0007045846009643734),
                C::new(-0.019040696846734557, 0.4233446838467155),
                C::new(-0.858558140904525, 7.898928160769499),
                C::new(-0.5380071035116917, -0.28433174603545897),
            ],
        ),
        (
            C::new(15.2, 2.9),
            0.01,
            [
                C::new(1.0, 9.19178947666843e-86),
                C::new(3509055013.15981, 113893202505.11443),
                C::new(25.459001430512064, 7.823181329532258),
                C::new(1.3159226362442258, -1.2953282450429418),
            ],
        ),
    ] {
        let ball = B::new(center, radius).unwrap();
        let results = [
            ball.erf(),
            ball.gamma(),
            ball.ln_gamma(),
            ball.bessel_j(0.75),
        ];
        for (result, expected) in results.into_iter().zip(references) {
            let result = result.unwrap_or_else(|| panic!("missing enclosure at {center}"));
            assert!(
                result.contains(expected),
                "{expected} lies outside {result:?} at {center}"
            );
        }
    }
}
