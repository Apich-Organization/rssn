//! Source-inventory families execute through registered codecs, not metadata.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::RecordValue;
use rssn::kernel::value::TensorValue;
fn call(
    s: &Session,
    name: &str,
    args: &[Value],
) -> Result<Value, OpError> {
    s.call_value(name, args, &CallCtx::new())
}
fn ints(xs: &[i64]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![xs.len()],
        data: xs.iter().copied().map(Value::int).collect(),
        elem: TypeExpr::Int,
    })
}
fn measured(
    shape: Vec<usize>,
    xs: &[f64],
    unit: &str,
) -> Value {
    Value::Tensor(TensorValue {
        shape,
        data: xs.iter().copied().map(Value::Float).collect(),
        elem: TypeExpr::Float.with_dim(DimExpr::Base(unit.into())),
    })
}
fn field<'a>(
    value: &'a Value,
    name: &str,
) -> &'a Value {
    let Value::Record(r) = value else {
        panic!("record")
    };
    r.get(name).unwrap()
}
fn csr(
    offsets: &[i64],
    indices: &[i64],
) -> Value {
    Value::Record(RecordValue {
        fields: vec![
            ("offsets".into(), ints(offsets)),
            ("indices".into(), ints(indices)),
        ],
    })
}
#[test]
fn coding_and_finite_field_codecs_execute_and_reject_out_of_range_values() {
    let s = Session::new();
    assert_eq!(
        call(
            &s,
            "num_error_correction_crc32_compute_numerical",
            &[ints(&[49, 50, 51, 52, 53, 54, 55, 56, 57])]
        )
        .unwrap(),
        Value::int(3_421_780_262u64)
    );
    let inv = call(&s, "num_finite_field_gf256_inv", &[Value::int(83)]).unwrap();
    assert_eq!(
        call(&s, "num_finite_field_gf256_mul", &[Value::int(83), inv]).unwrap(),
        Value::int(1)
    );
    assert!(
        call(
            &s,
            "num_finite_field_gf256_mul",
            &[Value::int(256), Value::int(1)]
        )
        .is_err()
    );
    let encoded = call(
        &s,
        "num_error_correction_reed_solomon_encode",
        &[ints(&[1, 2, 3]), Value::int(4)],
    )
    .unwrap();
    let mut damaged = encoded.clone();
    let Value::Tensor(t) = &mut damaged else {
        panic!("block")
    };
    t.data[1] = Value::int(77);
    assert_eq!(
        call(
            &s,
            "num_error_correction_reed_solomon_decode",
            &[damaged, Value::int(4)]
        )
        .unwrap(),
        encoded
    );
    for depth in [1_048_577usize, usize::MAX] {
        assert!(
            call(
                &s,
                "num_error_correction_interleave",
                &[ints(&[1, 2]), Value::int(depth)]
            )
            .is_err()
        );
        assert!(
            call(
                &s,
                "num_error_correction_deinterleave",
                &[ints(&[1, 2]), Value::int(depth)]
            )
            .is_err()
        );
    }
}
#[test]
fn topology_csr_and_persistence_preserve_topology_radius_dimensions_and_essential_classes() {
    let s = Session::new();
    let circle = csr(&[0, 1, 2, 3, 5, 7, 9], &[0, 1, 2, 0, 1, 1, 2, 0, 2]);
    assert_eq!(
        call(
            &s,
            "num_topology_betti_number",
            &[circle.clone(), Value::int(1)]
        )
        .unwrap(),
        Value::int(1)
    );
    assert_eq!(
        call(&s, "num_topology_euler_characteristic", &[circle]).unwrap(),
        Value::int(0)
    );
    let points = measured(vec![2, 1], &[0.0, 2.0], "Length");
    let radius = Value::Dyn(rssn::kernel::value::DynValue {
        value: Box::new(Value::Float(3.0)),
        dim: DimExpr::Base("Length".into()),
    });
    let result = call(
        &s,
        "num_topology_persistent_homology_rips",
        &[points.clone(), radius, Value::int(1), Value::int(100)],
    )
    .unwrap();
    let Value::Tensor(births) = field(&result, "births") else {
        panic!("births")
    };
    assert!(matches!(&births.elem,TypeExpr::Dim(_,DimExpr::Base(name)) if name=="Length"));
    let Value::Tensor(essential) = field(&result, "essential") else {
        panic!("essential")
    };
    assert_eq!(
        essential
            .data
            .iter()
            .filter(|v| matches!(v, Value::Bool(true)))
            .count(),
        1
    );
    let Value::Tensor(deaths) = field(&result, "deaths") else {
        panic!("deaths")
    };
    assert!(deaths.data.iter().all(|v| {
        match v {
            | Value::Dyn(d) => d.value.as_f64().is_some_and(f64::is_finite),
            | _ => v.as_f64().is_some_and(f64::is_finite),
        }
    }));
    let bad_radius = Value::Dyn(rssn::kernel::value::DynValue {
        value: Box::new(Value::Float(3.0)),
        dim: DimExpr::Base("Time".into()),
    });
    assert!(
        call(
            &s,
            "num_topology_vietoris_rips_complex",
            &[points, bad_radius, Value::int(1)]
        )
        .is_err()
    );
    assert!(call(&s, "num_topology_close_complex", &[csr(&[0, 2, 1], &[0])]).is_err());
}
#[test]
fn normalized_fractals_and_exact_statistics_reach_real_registered_kernels() {
    let s = Session::new();
    let orbit = call(
        &s,
        "num_fractal_geometry_and_chaos_logistic_map_iterate",
        &[Value::Float(0.5), Value::Float(4.0), Value::int(2)],
    )
    .unwrap();
    let Value::Tensor(t) = orbit else {
        panic!("orbit")
    };
    assert_eq!(
        t.data,
        [Value::Float(0.5), Value::Float(1.0), Value::Float(0.0)]
    );
    let fern = call(&s, "num_fractal_geometry_and_chaos_barnsley_fern_ifs", &[]).unwrap();
    let origin = Value::Record(RecordValue {
        fields: vec![
            ("item0".into(), Value::Float(0.0)),
            ("item1".into(), Value::Float(0.0)),
        ],
    });
    let points = call(
        &s,
        "num_fractal_geometry_and_chaos_generate_ifs_fractal",
        &[
            field(&fern, "item0").clone(),
            field(&fern, "item1").clone(),
            origin,
            Value::int(10),
            Value::int(5),
        ],
    )
    .unwrap();
    let Value::Tensor(points) = points else {
        panic!("IFS")
    };
    assert_eq!(points.shape, [10, 2]);
    let rational = num_rational::BigRational::new(1.into(), 2.into());
    let expected = num_rational::BigRational::new((-1).into(), 8.into());
    assert_eq!(
        call(
            &s,
            "num_special_exact_legendre_p_exact",
            &[Value::int(2), Value::rational(rational)]
        )
        .unwrap(),
        Value::rational(expected)
    );
    assert!(
        call(
            &s,
            "num_fractal_geometry_and_chaos_logistic_map_iterate",
            &[Value::Float(0.5), Value::Float(4.0), Value::int(usize::MAX)]
        )
        .is_err()
    );
}

#[test]
fn permutation_algorithms_validate_bijections_and_preserve_ragged_cycles() {
    let s = Session::new();
    let permutation = ints(&[1, 0, 2]);
    assert_eq!(
        call(&s, "num_permgroup_is_permutation", &[permutation.clone()]).unwrap(),
        Value::Bool(true)
    );
    let cycles = call(&s, "num_permgroup_cycles", &[permutation.clone()]).unwrap();
    assert_eq!(field(&cycles, "offsets"), &ints(&[0, 2, 3]));
    assert_eq!(field(&cycles, "indices"), &ints(&[0, 1, 2]));
    assert_eq!(
        call(&s, "num_permgroup_element_order", &[permutation.clone()]).unwrap(),
        Value::int(2)
    );
    assert_eq!(
        call(&s, "num_permgroup_power", &[permutation, Value::int(2)]).unwrap(),
        ints(&[0, 1, 2])
    );
    assert!(call(&s, "num_permgroup_cycles", &[ints(&[1, 1, 2])]).is_err());
    assert!(
        call(
            &s,
            "num_permgroup_conjugate",
            &[ints(&[0, 1]), ints(&[0, 1, 2])]
        )
        .is_err()
    );
    let generators = Value::Tensor(TensorValue {
        shape: vec![1, 3],
        data: vec![Value::int(1), Value::int(2), Value::int(0)],
        elem: TypeExpr::Int,
    });
    let orbits = call(
        &s,
        "num_permgroup_orbits",
        &[Value::int(3), generators.clone()],
    )
    .unwrap();
    assert_eq!(field(&orbits, "offsets"), &ints(&[0, 3]));
    assert_eq!(
        call(
            &s,
            "num_permgroup_is_primitive",
            &[Value::int(3), generators.clone()]
        )
        .unwrap(),
        Value::Bool(true)
    );
    assert!(
        call(
            &s,
            "num_permgroup_minimal_blocks",
            &[Value::int(3), generators, Value::int(0), Value::int(3)]
        )
        .is_err()
    );
    let empty_oversized = Value::Tensor(TensorValue {
        shape: vec![usize::MAX, 0],
        data: vec![],
        elem: TypeExpr::Int,
    });
    assert!(
        call(
            &s,
            "num_permgroup_orbits",
            &[Value::int(0), empty_oversized]
        )
        .is_err()
    );
}
