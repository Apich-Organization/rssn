//! Native homology objects execute real kernel algorithms with checked lifetimes.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::RecordValue;
use rssn::kernel::value::TensorValue;
fn call(
    s: &Session,
    name: &str,
    args: &[Value],
) -> Result<Value, OpError> {
    s.call_value(name, args, &CallCtx::new())
}
fn ints(data: &[i64]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![data.len()],
        data: data.iter().copied().map(Value::int).collect(),
        elem: TypeExpr::Int,
    })
}
fn measured(data: &[f64]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![data.len()],
        data: data.iter().copied().map(Value::Float).collect(),
        elem: TypeExpr::Float.with_dim(DimExpr::Base("Length".into())),
    })
}
fn record(fields: Vec<(&str, Value)>) -> Value {
    Value::Record(RecordValue {
        fields: fields
            .into_iter()
            .map(|(key, value)| (key.into(), value))
            .collect(),
    })
}
fn field<'a>(
    v: &'a Value,
    key: &str,
) -> &'a Value {
    let Value::Record(r) = v else {
        panic!("record")
    };
    r.get(key).unwrap()
}
#[test]
fn retained_projective_plane_chains_have_exact_torsion_and_checked_session_identity() {
    let s = Session::new();
    let foreign = Session::new();
    let manifest = s.manifest();
    for (name, purity) in [
        (
            "num_homology_complex_chains",
            rssn::kernel::signature::Purity::Effectful,
        ),
        (
            "num_homology_chain_homology",
            rssn::kernel::signature::Purity::ReadOnly,
        ),
    ] {
        let info = manifest.op(name).unwrap();
        assert_eq!(info.purity, purity);
        assert_eq!(
            info.thread_safety,
            rssn::kernel::signature::ThreadSafety::SessionBound
        );
        assert_eq!(
            info.determinism,
            rssn::kernel::signature::Determinism::EnvDependent
        );
        assert!(info.backend_support.interpreter);
        assert!(
            !info.backend_support.jit && !info.backend_support.ffi && !info.backend_support.gpu
        );
        assert!(
            info.overloads[0]
                .effects
                .effects
                .contains(&rssn::kernel::signature::Effect::Stateful)
        );
    }
    let complex = call(&s, "num_homology_projective_plane", &[Value::int(0)]).unwrap();
    let chains = call(&s, "num_homology_complex_chains", &[complex.clone()]).unwrap();
    let group = call(
        &s,
        "num_homology_chain_homology",
        &[chains.clone(), Value::int(1)],
    )
    .unwrap();
    assert_eq!(field(&group, "rank"), &Value::int(0));
    assert_eq!(field(&group, "torsion"), &ints(&[2]));
    let cohomology = call(
        &s,
        "num_homology_chain_cohomology",
        &[chains.clone(), Value::int(2)],
    )
    .unwrap();
    assert_eq!(field(&cohomology, "torsion"), &ints(&[2]));
    assert_eq!(
        call(
            &s,
            "num_homology_chain_homology_mod_p",
            &[chains.clone(), Value::int(1), Value::int(2)]
        )
        .unwrap(),
        Value::int(1)
    );
    assert!(
        call(
            &s,
            "num_homology_chain_homology_mod_p",
            &[chains.clone(), Value::int(1), Value::int(4)]
        )
        .is_err()
    );
    assert!(
        call(
            &foreign,
            "num_homology_chain_homology",
            &[chains.clone(), Value::int(1)]
        )
        .is_err()
    );
    assert!(call(&s, "num_homology_chain_homology", &[complex, Value::int(1)]).is_err());
    call(&s, "num_native_release", &[chains.clone()]).unwrap();
    assert!(call(&s, "num_homology_chain_homology", &[chains, Value::int(1)]).is_err());
}
#[test]
fn filtration_objects_preserve_radius_units_and_explicit_essential_intervals() {
    let s = Session::new();
    let simplices = record(vec![
        ("offsets", ints(&[0, 1, 2, 4])),
        ("indices", ints(&[0, 1, 0, 1])),
    ]);
    let input = record(vec![
        ("radii", measured(&[0.0, 0.0, 2.0])),
        ("simplices", simplices),
    ]);
    let model = call(&s, "num_homology_filtration", &[input]).unwrap();
    let bars = call(
        &s,
        "num_homology_persistence_bars",
        &[model.clone(), Value::int(0)],
    )
    .unwrap();
    let Value::Variant(some) = bars else {
        panic!("optional bars")
    };
    assert_eq!(some.tag, "Some");
    let bars = &some.payload[0];
    let Value::Tensor(births) = field(bars, "births") else {
        panic!("births")
    };
    assert!(matches!(&births.elem,TypeExpr::Dim(_,DimExpr::Base(name))if name=="Length"));
    let Value::Tensor(essential) = field(bars, "essential") else {
        panic!("essential")
    };
    assert_eq!(essential.data, [Value::Bool(false), Value::Bool(true)]);
    let mut forged = model.clone();
    let Value::Record(r) = &mut forged else {
        panic!("model")
    };
    r.fields
        .iter_mut()
        .find(|(key, _)| key == "value_unit")
        .unwrap()
        .1 = Value::Dyn(DynValue {
        value: Box::new(Value::Float(1.0)),
        dim: DimExpr::Base("Time".into()),
    });
    assert!(
        call(
            &s,
            "num_homology_persistence_bars",
            &[forged, Value::int(0)]
        )
        .is_err()
    );
    let diagram = record(vec![
        ("births", field(bars, "births").clone()),
        ("deaths", field(bars, "deaths").clone()),
        ("essential", field(bars, "essential").clone()),
    ]);
    let distance = call(
        &s,
        "num_homology_bottleneck_distance",
        &[diagram.clone(), diagram],
    )
    .unwrap();
    let Value::Dyn(distance) = distance else {
        panic!("measured distance")
    };
    assert_eq!(*distance.value, Value::Float(0.0));
    assert_eq!(distance.dim, DimExpr::Base("Length".into()));
}
#[test]
fn exact_matrix_and_complex_constructions_reject_invalid_budgets_before_math() {
    let s = Session::new();
    let matrix = Value::Tensor(TensorValue {
        shape: vec![2, 2],
        data: vec![Value::int(2), Value::int(0), Value::int(0), Value::int(4)],
        elem: TypeExpr::Int,
    });
    assert_eq!(
        call(&s, "num_homology_smith_invariants", &[matrix.clone()]).unwrap(),
        ints(&[2, 4])
    );
    assert!(call(&s, "num_homology_rank_mod_p", &[matrix, Value::int(0)]).is_err());
    assert!(call(&s, "num_homology_full_simplex", &[Value::int(usize::MAX)]).is_err());
    let sphere = call(&s, "num_homology_sphere", &[Value::int(1)]).unwrap();
    assert_eq!(
        call(&s, "num_homology_betti_list", &[sphere]).unwrap(),
        ints(&[1, 1])
    );
}
