//! Small physical contracts for boundary elements and mixed-unit continua.
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::Exp;
use rssn::kernel::signature::Rational;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
fn call(
    s: &rssn::Session,
    name: &str,
    args: &[Value],
) -> Value {
    s.call_value(name, args, &CallCtx::new())
        .unwrap_or_else(|e| panic!("{name}: {e:?}"))
}
fn f(v: f64) -> Value {
    Value::Float(v)
}
fn i(v: usize) -> Value {
    Value::int(v)
}
fn t(
    shape: Vec<usize>,
    xs: &[f64],
) -> Value {
    Value::Tensor(TensorValue {
        shape,
        elem: TypeExpr::Float,
        data: xs.iter().copied().map(f).collect(),
    })
}
fn q(
    v: f64,
    d: DimExpr,
) -> Value {
    Value::Dyn(DynValue {
        dim: d,
        value: Box::new(f(v)),
    })
}
fn base(n: &str) -> DimExpr {
    DimExpr::base(n)
}
fn power(
    n: &str,
    p: i64,
) -> DimExpr {
    DimExpr::Pow(Box::new(base(n)), Exp::Rat(Rational::new(p, 1)))
}
fn scalar(v: &Value) -> f64 {
    match v {
        | Value::Dyn(v) => scalar(&v.value),
        | _ => v.as_f64().unwrap(),
    }
}
fn field<'a>(
    v: &'a Value,
    n: &str,
) -> &'a Value {
    let Value::Record(r) = v else {
        panic!("record")
    };
    r.get(n).unwrap()
}
fn refs(s: &rssn::Session) -> Value {
    call(
        s,
        "num_grid_mechanical_references",
        &[
            q(2., base("Length")),
            q(4., base("Time")),
            q(3., base("Mass")),
            q(1., base("Temperature")),
        ],
    )
}
fn frame(
    s: &rssn::Session,
    u: Value,
) -> Value {
    call(
        s,
        "num_grid_frame_new",
        &[q(2., base("Length")), q(4., base("Time")), u],
    )
}
#[test]
fn native_bem_geometry_retains_length_area_volume_and_bounds() {
    let s = rssn::Session::new();
    let frame = frame(&s, q(3., base("Temperature")));
    let mesh = call(
        &s,
        "num_grid_bem_mesh_cube",
        &[frame.clone(), f(1.), t(vec![3], &[0.; 3]), i(1)],
    );
    assert!((scalar(&call(&s, "num_grid_bem_mesh_volume", &[mesh.clone()])) - 8.).abs() < 1e-12);
    assert!((scalar(&call(&s, "num_grid_bem_mesh_total_area", &[mesh])) - 24.).abs() < 1e-12);
    let v = call(
        &s,
        "num_grid_bem_vector2_new",
        &[frame.clone(), t(vec![2], &[3., 4.])],
    );
    assert_eq!(scalar(&call(&s, "num_grid_bem_vector2_norm", &[v])), 10.);
    assert!(
        s.call_value(
            "num_grid_bem_mesh_icosphere",
            &[frame, f(1.), t(vec![3], &[0.; 3]), i(usize::MAX)],
            &CallCtx::new()
        )
        .is_err()
    );
}
#[test]
fn native_hyperelastic_identity_stress_and_lame_references() {
    let s = rssn::Session::new();
    let refs = refs(&s);
    let material = call(
        &s,
        "num_grid_hyper_material",
        &[
            refs,
            Value::String("NeoHookean".into()),
            f(10.),
            f(0.25),
            f(1.),
        ],
    );
    let lame = call(&s, "num_grid_hyper_lame", &[material.clone()]);
    let Value::Tensor(values) = field(&lame, "normalized") else {
        panic!("vector")
    };
    assert_eq!(
        values.data.iter().map(scalar).collect::<Vec<_>>(),
        vec![4., 4.]
    );
    assert!((scalar(field(&lame, "reference")) - 3. / 32.).abs() < 1e-15);
    let result = call(
        &s,
        "num_grid_hyper_piola",
        &[
            material,
            t(vec![3, 3], &[1., 0., 0., 0., 1., 0., 0., 0., 1.]),
        ],
    );
    let Value::Tensor(entries) = field(field(&result, "entries"), "normalized") else {
        panic!("entries")
    };
    assert!(entries.data.iter().all(|v| scalar(v) == 0.));
}
#[test]
fn native_gas_flux_separates_mass_momentum_and_energy_dimensions() {
    let s = rssn::Session::new();
    let result = call(
        &s,
        "num_grid_gas_hllc",
        &[
            refs(&s),
            t(vec![3], &[1., 2., 1.]),
            t(vec![3], &[1., 2., 1.]),
            f(1.4),
        ],
    );
    assert!((scalar(field(&result, "mass_flux")) - 0.375).abs() < 1e-12);
    assert!((scalar(field(&result, "momentum_flux")) - 0.46875).abs() < 1e-12);
    assert!((scalar(field(&result, "energy_flux")) - 0.515625).abs() < 1e-12);
    let Value::Dyn(mass) = field(&result, "mass_flux") else {
        panic!("measured mass flux")
    };
    let expected = DimExpr::Mul(vec![base("Mass"), power("Length", -2), power("Time", -1)]);
    assert!(
        rssn::kernel::overload::matches_type(
            &TypeExpr::Float.with_dim(expected),
            &TypeExpr::Float.with_dim(mass.dim.clone())
        )
        .accepted()
    );
}
#[test]
fn native_acoustics_reduced_energy_and_invalid_work_budget() {
    let s = rssn::Session::new();
    let pressure = DimExpr::Mul(vec![base("Mass"), base("Length"), power("Time", -2)]);
    let frame = frame(&s, q(0.375, pressure));
    let neumann = call(&s, "num_grid_neumann", &[frame.clone(), f(0.)]);
    let faces = call(&s, "num_grid_faces_uniform", &[neumann]);
    let config = call(
        &s,
        "num_grid_acoustics_config",
        &[
            refs(&s),
            faces,
            Value::Tensor(TensorValue {
                shape: vec![1],
                elem: TypeExpr::Int,
                data: vec![i(4)],
            }),
            t(vec![1], &[1.]),
            f(1.),
            f(1.),
            f(0.5),
            i(0),
            f(3.),
        ],
    );
    let solver = call(&s, "num_grid_acoustics_new", &[config]);
    call(&s, "num_grid_acoustics_rest", &[solver.clone()]);
    call(&s, "num_grid_acoustics_run", &[solver.clone(), i(3)]);
    assert_eq!(
        scalar(&call(&s, "num_grid_acoustics_energy", &[solver.clone()])),
        0.
    );
    assert!(
        s.call_value(
            "num_grid_acoustics_run",
            &[solver.clone(), i(usize::MAX)],
            &CallCtx::new()
        )
        .is_err()
    );
    assert_eq!(
        scalar(&call(&s, "num_grid_acoustics_energy", &[solver])),
        0.
    );
    let bad_neumann = call(&s, "num_grid_neumann", &[frame, f(1.)]);
    let faces = call(&s, "num_grid_faces_uniform", &[bad_neumann]);
    assert!(
        s.call_value(
            "num_grid_acoustics_config",
            &[
                refs(&s),
                faces,
                Value::Tensor(TensorValue {
                    shape: vec![1],
                    elem: TypeExpr::Int,
                    data: vec![i(4)]
                }),
                t(vec![1], &[1.]),
                f(1.),
                f(1.),
                f(0.5),
                i(0),
                f(3.)
            ],
            &CallCtx::new()
        )
        .is_err()
    );
}
#[test]
fn native_fem_scenario_returns_its_measured_field_reference() {
    let s = rssn::Session::new();
    let frame = frame(&s, q(3., base("Temperature")));
    let result = call(&s, "num_grid_fem_legacy_scenario1", &[frame]);
    let Value::Tensor(data) = field(&result, "data") else {
        panic!("field tensor")
    };
    assert_eq!(data.data.len(), 51);
    assert!((scalar(&data.data[25]) - 0.25).abs() < 1e-9);
    assert_eq!(scalar(field(&result, "value_reference")), 3.);
}
#[test]
fn native_quadrilateral_checks_jacobian_before_assembly() {
    let s = rssn::Session::new();
    let refs = refs(&s);
    let result = call(
        &s,
        "num_grid_quad_element",
        &[
            refs.clone(),
            t(vec![4, 2], &[0., 0., 1., 0., 1., 1., 0., 1.]),
            f(10.),
            f(0.25),
        ],
    );
    let Value::Tensor(values) = field(&result, "entries") else {
        panic!("matrix")
    };
    assert_eq!(values.data.len(), 64);
    assert!(
        s.call_value(
            "num_grid_quad_element",
            &[
                refs,
                t(vec![4, 2], &[0., 0., 0., 1., 1., 1., 1., 0.]),
                f(10.),
                f(0.25)
            ],
            &CallCtx::new()
        )
        .is_err()
    );
}

#[test]
fn native_bem_dense_solution_keeps_potential_and_flux_reference_separate() {
    let s = rssn::Session::new();
    let frame = frame(&s, q(3., base("Temperature")));
    let mesh = call(
        &s,
        "num_grid_bem_mesh_cube",
        &[frame, f(1.), t(vec![3], &[0.; 3]), i(1)],
    );
    let kinds = Value::Tensor(TensorValue {
        shape: vec![12],
        elem: TypeExpr::String,
        data: vec![Value::String("Potential".into()); 12],
    });
    let solution = call(
        &s,
        "num_grid_bem_solve3",
        &[mesh, kinds, t(vec![12], &[1.; 12])],
    );
    let fields = call(&s, "num_grid_bem_solution_fields", &[solution.clone()]);
    assert_eq!(scalar(field(&fields, "potential_reference")), 3.);
    assert_eq!(scalar(field(&fields, "flux_reference")), 1.5);
    let value = call(&s, "num_grid_bem_eval3", &[solution, t(vec![3], &[0.; 3])]);
    assert!(
        (scalar(&value) - 3.).abs() < 0.03,
        "constant harmonic potential: {value:?}"
    );
}
