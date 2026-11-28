//! Portable NTT plans preserve their exact DFT convention and reuse cached tables.

use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::manifest::SignatureSource;
use rssn::kernel::signature::DimRuleInfo;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::RecordValue;
use rssn::kernel::value::TensorValue;
use rssn::kernels::finite_field;

const GOLDILOCKS: u64 = 18_446_744_069_414_584_321;

fn vector(values: &[u64]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        data: values.iter().copied().map(Value::int).collect(),
        elem: TypeExpr::Int,
    })
}

#[test]
fn forward_matches_the_exact_small_dft() {
    const MODULUS: u64 = 17;
    const ROOT: u64 = 9;
    let input = [1, 2, 3, 4, u64::MAX, 12, 6, 9];
    let mut lookup = [1; 8];
    for index in 1..lookup.len() {
        lookup[index] = lookup[index - 1] * ROOT % MODULUS;
    }
    let mut output = [0; 8];
    for (frequency, coefficient) in output.iter_mut().enumerate() {
        *coefficient = input.iter().enumerate().fold(0, |acc, (index, value)| {
            (acc + (value % MODULUS) * lookup[(index * frequency) % input.len()]) % MODULUS
        });
    }
    let session = Session::new();
    let ctx = CallCtx::new();
    let plan = session
        .call_value(
            "num_finite_field_ntt_plan",
            &[Value::int(MODULUS), Value::int(8)],
            &ctx,
        )
        .unwrap();
    let spectrum = session
        .call_value(
            "num_finite_field_ntt_transform",
            &[vector(&input), plan.clone(), Value::Bool(false)],
            &ctx,
        )
        .unwrap();
    assert_eq!(spectrum, vector(&output));
    assert_eq!(
        session
            .call_value(
                "num_finite_field_ntt_transform",
                &[spectrum, plan, Value::Bool(true)],
                &ctx
            )
            .unwrap(),
        vector(&input.map(|value| value % MODULUS))
    );
}

#[test]
fn inverse_recomposes_full_width_coefficients_across_serial_and_parallel() {
    for length in [1, 2, 8, 32_768, 65_536] {
        let plan = finite_field::ntt_plan(GOLDILOCKS, length).unwrap();
        let again = finite_field::ntt_plan(GOLDILOCKS, length).unwrap();
        assert!(std::sync::Arc::ptr_eq(&plan, &again));
        let input: Vec<_> = (0..length)
            .map(|index| u64::MAX - (index as u64 * 891_683_231))
            .map(|value| plan.mont.to_mont(value))
            .collect();
        let run = |threads| {
            let mut data = input.clone();
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            let freq = pool.install(|| {
                finite_field::ntt_in_place_mont(&mut data, &plan, false);
                let freq = data.clone();
                finite_field::ntt_in_place_mont(&mut data, &plan, true);
                freq
            });
            assert_eq!(input, data);
            freq
        };
        assert_eq!(run(1), run(4));
    }
}

#[test]
fn rejects_invalid_or_forged_plans_before_transformation() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let plan_op = "num_finite_field_ntt_plan";
    let transform_op = "num_finite_field_ntt_transform";
    let manifest = session.manifest();
    for name in [plan_op, transform_op] {
        let op = manifest.op(name).unwrap();
        assert_eq!(op.signature_source, SignatureSource::Declared);
        assert_eq!(op.dim_rule, Some(DimRuleInfo::Dimensionless));
    }
    for (p, length) in [
        (2, 1),
        (341, 1),
        (341, 4),
        (17, 0),
        (17, 6),
        (17, 32),
        (GOLDILOCKS, 262_144),
    ] {
        assert!(
            session
                .call_value(plan_op, &[Value::int(p), Value::int(length)], &ctx)
                .is_err()
        );
    }
    let plan = session
        .call_value(plan_op, &[Value::int(17), Value::int(4)], &ctx)
        .unwrap();
    assert!(
        session
            .call_value(
                transform_op,
                &[vector(&[1, 2]), plan, Value::Bool(false)],
                &ctx
            )
            .is_err()
    );
    let forged = Value::Record(RecordValue {
        fields: vec![
            ("length".into(), Value::int(4)),
            ("modulus".into(), Value::int(341)),
        ],
    });
    assert!(
        session
            .call_value(
                transform_op,
                &[vector(&[1, 2, 3, 4]), forged, Value::Bool(true)],
                &ctx
            )
            .is_err()
    );
}
