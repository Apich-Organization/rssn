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
    assert!(
        (scalar(&call(
            &s,
            "num_grid_bem_mesh_volume",
            std::slice::from_ref(&mesh)
        )) - 8.)
            .abs()
            < 1e-12
    );
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
    let lame = call(&s, "num_grid_hyper_lame", std::slice::from_ref(&material));
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
fn prepared_hyperelastic_point_preserves_pressure_shape_and_tangent_action() {
    let session = rssn::Session::new();
    let material = call(
        &session,
        "num_grid_hyper_material",
        &[
            refs(&session),
            Value::String("NeoHookean".into()),
            f(10.0),
            f(0.25),
            f(1.0),
        ],
    );
    let identity = t(vec![3, 3], &[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]);
    let inverted = t(vec![3, 3], &[-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]);
    assert!(
        session
            .call_value(
                "num_grid_hyper_prepare_point",
                &[material.clone(), inverted],
                &CallCtx::new()
            )
            .is_err()
    );
    let point = call(
        &session,
        "num_grid_hyper_prepare_point",
        &[material.clone(), identity],
    );
    call(&session, "num_native_release", &[material]);
    let query = |name| call(&session, name, std::slice::from_ref(&point));
    let energy = query("num_grid_hyper_point_energy_density");
    assert_eq!(scalar(&energy), 0.0);
    let Value::Tensor(stress) = query("num_grid_hyper_point_piola") else {
        panic!("stress tensor")
    };
    assert_eq!(stress.shape, [3, 3]);
    assert!(stress.data.iter().all(|value| scalar(value) == 0.0));
    let Value::Tensor(tangent) = query("num_grid_hyper_point_tangent") else {
        panic!("tangent tensor")
    };
    assert_eq!(tangent.shape, [3, 3, 3, 3]);
    let pressure = DimExpr::Mul(vec![base("Mass"), power("Length", -1), power("Time", -2)]);
    assert!(
        rssn::kernel::overload::matches_type(&TypeExpr::Float.with_dim(pressure), &tangent.elem)
            .accepted()
    );
    let direction = t(vec![3, 3], &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    let Value::Tensor(action) = call(
        &session,
        "num_grid_hyper_point_tangent_action",
        &[point, direction],
    ) else {
        panic!("applied tangent tensor");
    };
    assert_eq!(action.shape, [3, 3]);
    let action = action.data.iter().map(scalar).collect::<Vec<_>>();
    assert_eq!(action, [1.125, 0.0, 0.0, 0.0, 0.375, 0.0, 0.0, 0.0, 0.375]);
    assert_eq!(
        action,
        tangent
            .data
            .iter()
            .step_by(9)
            .map(scalar)
            .collect::<Vec<_>>()
    );
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
    call(&s, "num_grid_acoustics_rest", std::slice::from_ref(&solver));
    call(&s, "num_grid_acoustics_run", &[solver.clone(), i(3)]);
    assert_eq!(
        scalar(&call(
            &s,
            "num_grid_acoustics_energy",
            std::slice::from_ref(&solver)
        )),
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
fn native_quadrilateral_solves_an_affine_patch_with_a_measured_reference() {
    let session = rssn::Session::new();
    let op = "num_grid_quad_solve";
    let fixed = |nodes: Vec<Value>| {
        Value::Tensor(TensorValue {
            shape: vec![nodes.len()],
            data: nodes,
            elem: TypeExpr::Int,
        })
    };
    let mut args = [
        refs(&session),
        t(vec![4, 2], &[0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0]),
        t(vec![1, 4], &[0.0, 1.0, 2.0, 3.0]),
        f(2.0),
        f(0.0),
        fixed(vec![i(0), i(3)]),
        t(vec![2, 3], &[1.0, 0.5, 0.0, 2.0, 0.5, 0.0]),
    ];
    let result = call(&session, op, &args);
    let Value::Dyn(reference) = field(&result, "reference") else {
        panic!("expected a measured Length reference");
    };
    assert_eq!(reference.dim, base("Length"));
    assert_eq!(scalar(&reference.value), 2.0);
    let Value::Tensor(displacement) = field(&result, "normalized_displacement") else {
        panic!("expected a displacement vector");
    };
    assert_eq!(displacement.shape, [8]);
    // Linear u = x/2, v = 0 gives unit stress and consistent right-edge nodal forces.
    for (node, expected) in displacement
        .data
        .iter()
        .zip([0., 0., 0.5, 0., 0.5, 0., 0., 0.])
    {
        assert!((scalar(node) - expected).abs() < 1e-12);
    }

    args[5] = fixed([i(3), i(0)].into_iter().cycle().take(4096).collect());
    assert_eq!(call(&session, op, &args), result);
    args[5] = fixed(vec![i(0), i(4)]);
    assert!(session.call_value(op, &args, &CallCtx::new()).is_err());
    args[5] = fixed(vec![i(0), i(3)]);
    args[2] = t(vec![1, 4], &[0.0, 3.0, 2.0, 1.0]);
    assert!(session.call_value(op, &args, &CallCtx::new()).is_err());
}

#[test]
fn native_quadrilateral_cantilever_writes_both_coordinate_tables() {
    let session = rssn::Session::new();
    let path = std::env::temp_dir().join(format!("rssn_native_q4_{}", std::process::id()));
    std::fs::create_dir(&path).unwrap();
    let original = path.join("beam_original.csv");
    let deformed = path.join("beam_deformed.csv");
    let result = session.call_value(
        "num_grid_quad_cantilever",
        &[Value::String(path.to_str().unwrap().into())],
        &CallCtx::new(),
    );
    let before = std::fs::read_to_string(&original);
    let after = std::fs::read_to_string(&deformed);
    if before.is_ok() {
        std::fs::remove_file(original).unwrap();
    }
    if after.is_ok() {
        std::fs::remove_file(deformed).unwrap();
    }
    std::fs::remove_dir(path).unwrap();
    assert_eq!(result.unwrap(), Value::Bool(true));
    let before = before.unwrap();
    let after = after.unwrap();
    assert_eq!(before.lines().next(), Some("x,y"));
    assert_eq!(after.lines().next(), Some("x,y"));
    assert_eq!(before.lines().count(), 106);
    assert_eq!(after.lines().count(), 106);
    for (index, (before, after)) in before
        .lines()
        .skip(1)
        .zip(after.lines().skip(1))
        .enumerate()
    {
        if index % 21 == 0 {
            assert_eq!(before, after);
        }
        if index == 62 {
            let y = |row: &str| row.rsplit_once(',').unwrap().1.parse::<f64>().unwrap();
            assert!(y(after) < y(before));
        }
    }
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
    let fields = call(
        &s,
        "num_grid_bem_solution_fields",
        std::slice::from_ref(&solution),
    );
    assert_eq!(scalar(field(&fields, "potential_reference")), 3.);
    assert_eq!(scalar(field(&fields, "flux_reference")), 1.5);
    let value = call(&s, "num_grid_bem_eval3", &[solution, t(vec![3], &[0.; 3])]);
    assert!(
        (scalar(&value) - 3.).abs() < 0.03,
        "constant harmonic potential: {value:?}"
    );
}

#[test]
fn native_bem_exterior_retains_outward_flux_and_the_solved_domain() {
    let session = rssn::Session::new();
    let temperature = base("Temperature");
    let mesh = call(
        &session,
        "num_grid_bem_mesh_icosphere",
        &[
            frame(&session, q(3.0, temperature.clone())),
            f(1.0),
            t(vec![3], &[0.0; 3]),
            i(2),
        ],
    );
    let panel_values = |kind: &str| {
        Value::Tensor(TensorValue {
            shape: vec![320],
            elem: TypeExpr::String,
            data: vec![Value::String(kind.into()); 320],
        })
    };
    let args = [
        mesh,
        panel_values("Potential"),
        t(vec![320], &[1.0; 320]),
        i(8),
        f(1.5),
        f(1e-8),
        i(32),
        i(128),
    ];
    let exterior = call(&session, "num_grid_bem_exterior", &args);
    let fields = call(
        &session,
        "num_grid_bem_solution_fields",
        std::slice::from_ref(&exterior),
    );
    assert_eq!(scalar(field(&fields, "potential_reference")), 3.0);
    assert_eq!(scalar(field(&fields, "flux_reference")), 1.5);
    let Value::Tensor(flux) = field(&fields, "flux") else {
        panic!("expected panel flux")
    };
    assert_eq!(flux.shape, [320]);
    assert!(flux.data.iter().all(|q| (scalar(q) + 1.0).abs() < 0.15));
    for radius in [2.0, 4.0] {
        let value = call(
            &session,
            "num_grid_bem_eval3_exterior",
            &[exterior.clone(), t(vec![3], &[radius, 0.0, 0.0])],
        );
        let Value::Dyn(value) = value else {
            panic!("expected measured potential")
        };
        assert_eq!(value.dim, temperature);
        assert!(
            (scalar(&value.value) - 3.0 / radius).abs() < 0.035,
            "radius {radius}: {value:?}"
        );
    }
    let wrong_domain = session.call_value(
        "num_grid_bem_eval3",
        &[exterior, t(vec![3], &[0.0; 3])],
        &CallCtx::new(),
    );
    assert!(wrong_domain.is_err());
}

#[test]
fn native_bem_aca_matches_dense_and_reports_real_compression() {
    let session = rssn::Session::new();
    // Two well-separated surfaces give the hierarchy an unambiguous far block.
    let mut source = rssn::sim::physics_bem::SurfaceMesh3D::icosphere(1.0, [0.0; 3], 1);
    let other = rssn::sim::physics_bem::SurfaceMesh3D::icosphere(1.0, [16.0, 0.0, 0.0], 1);
    let offset = source.vertices.len();
    source.vertices.extend(other.vertices);
    source.triangles.extend(
        other
            .triangles
            .into_iter()
            .map(|face| face.map(|node| node + offset)),
    );
    let num_panels = source.triangles.len();
    let nodes: Vec<_> = source.vertices.iter().flatten().copied().collect();
    let faces: Vec<_> = source
        .triangles
        .iter()
        .flatten()
        .map(|&node| node as f64)
        .collect();
    let mesh = call(
        &session,
        "num_grid_bem_mesh_new",
        &[
            frame(&session, q(3.0, base("Temperature"))),
            t(vec![source.vertices.len(), 3], &nodes),
            t(vec![num_panels, 3], &faces),
        ],
    );
    let kinds = Value::Tensor(TensorValue {
        shape: vec![num_panels],
        elem: TypeExpr::String,
        data: vec![Value::String("Potential".into()); num_panels],
    });
    let values: Vec<_> = (0..num_panels)
        .map(|panel| source.centroid(panel).unwrap()[0])
        .collect();
    let mut args = [
        mesh,
        kinds,
        t(vec![num_panels], &values),
        i(8),
        f(2.0),
        f(1e-8),
        i(32),
        i(128),
    ];
    let dense = call(&session, "num_grid_bem_solve3", &args[..3]);
    let aca = call(&session, "num_grid_bem_aca", &args);
    let eval_point = t(vec![3], &[0.2, 0.15, -0.1]);
    let dense_value = call(
        &session,
        "num_grid_bem_eval3",
        &[dense.clone(), eval_point.clone()],
    );
    let aca_value = call(&session, "num_grid_bem_eval3", &[aca.clone(), eval_point]);
    assert!((scalar(&dense_value) - 0.6).abs() < 0.04);
    assert!((scalar(&dense_value) - scalar(&aca_value)).abs() < 1e-6);
    let ratio = call(
        &session,
        "num_grid_bem_compression",
        std::slice::from_ref(&aca),
    );
    assert!(
        (0.0..1.0).contains(&scalar(&ratio)),
        "compression {ratio:?}"
    );
    assert!(
        session
            .call_value("num_grid_bem_compression", &[dense], &CallCtx::new())
            .is_err()
    );
    args[7] = i(4097);
    let error = session
        .call_value("num_grid_bem_aca", &args, &CallCtx::new())
        .unwrap_err();
    assert!(
        error.to_string().contains("ACA options/work cap"),
        "{error}"
    );
    assert!(
        session
            .call_value(
                "num_grid_bem_eval3_exterior",
                &[aca, t(vec![3], &[2.0, 0.0, 0.0])],
                &CallCtx::new()
            )
            .is_err()
    );
}

#[test]
fn prepared_native_bem_reuses_panel_types_and_returns_ordinary_measured_solutions() {
    let s = rssn::Session::new();
    let temperature = base("Temperature");
    let frame = frame(&s, q(3.0, temperature.clone()));
    let mesh = call(
        &s,
        "num_grid_bem_mesh_cube",
        &[frame, f(1.0), t(vec![3], &[0.0; 3]), i(1)],
    );
    let kinds = |kind: &str| {
        Value::Tensor(TensorValue {
            shape: vec![12],
            elem: TypeExpr::String,
            data: vec![Value::String(kind.into()); 12],
        })
    };
    let plan = call(
        &s,
        "num_grid_bem_prepare3",
        &[mesh.clone(), kinds("Potential")],
    );
    assert_eq!(
        call(
            &s,
            "num_grid_bem_plan_panel_count",
            std::slice::from_ref(&plan)
        ),
        i(12)
    );
    for input in [1.0, -2.0] {
        let solution = call(
            &s,
            "num_grid_bem_solve3_prepared",
            &[plan.clone(), t(vec![12], &[input; 12])],
        );
        let value = call(&s, "num_grid_bem_eval3", &[solution, t(vec![3], &[0.0; 3])]);
        assert_eq!(value.dimension(), Some(&temperature));
        assert!((scalar(&value) - 3.0 * input).abs() < 0.03);
    }
    let flux = call(&s, "num_grid_bem_prepare3", &[mesh, kinds("Flux")]);
    let solve = "num_grid_bem_solve3_prepared";
    assert!(
        s.call_value(
            solve,
            &[flux.clone(), t(vec![12], &[1.0; 12])],
            &CallCtx::new()
        )
        .is_err()
    );
    assert!(
        s.call_value(solve, &[plan, t(vec![11], &[1.0; 11])], &CallCtx::new())
            .is_err()
    );
    let zero = call(&s, solve, &[flux, t(vec![12], &[0.0; 12])]);
    let value = call(&s, "num_grid_bem_eval3", &[zero, t(vec![3], &[0.0; 3])]);
    assert_eq!(scalar(&value), 0.0);
}
