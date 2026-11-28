//! Concrete source-inventory objects and callbacks execute through session APIs.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::NativeCallback;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::Exp;
use rssn::kernel::signature::NatExpr;
use rssn::kernel::signature::Rational;
use rssn::kernel::signature::Signature;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::ComplexValue;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::RecordValue;
use rssn::kernel::value::TensorValue;
use std::sync::Arc;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;
fn call(
    s: &Session,
    n: &str,
    args: &[Value],
) -> Value {
    s.call_value(n, args, &CallCtx::new())
        .unwrap_or_else(|e| panic!("{n}: {e:?}"))
}
fn dim(n: &str) -> DimExpr {
    DimExpr::Base(n.into())
}
fn scalar(
    x: f64,
    d: DimExpr,
) -> Value {
    if d == DimExpr::Dimless {
        Value::Float(x)
    } else {
        Value::Dyn(DynValue {
            dim: d,
            value: Box::new(Value::Float(x)),
        })
    }
}
fn number(x: &Value) -> f64 {
    match x {
        | Value::Dyn(x) => number(&x.value),
        | _ => x.as_f64().unwrap(),
    }
}

#[test]
fn exact_matrix_lie_basis_reuses_a_single_coordinate_factorization() {
    let session = Session::new();
    let matrix = |coefficients: [i64; 4]| {
        call(
            &session,
            "num_native_qmatrix_new",
            &[
                tensor(
                    vec![2, 2],
                    coefficients.into_iter().map(Value::int).collect(),
                    TypeExpr::Int,
                ),
                Value::int(1),
            ],
        )
    };
    let h = matrix([1, 0, 0, -1]);
    let e = matrix([0, 1, 0, 0]);
    let f = matrix([0, 0, 1, 0]);
    let kind = h.type_expr();
    let basis = |left, right| tensor(vec![2], vec![left, right], kind.clone());
    let result = call(
        &session,
        "num_native_lie_from_matrices",
        &[basis(h, e.clone())],
    );
    let Value::Variant(lie) = result else {
        panic!("expected an optional Lie algebra")
    };
    assert_eq!(lie.tag, "Some");
    let coordinate = |entries: [i64; 2]| {
        tensor(
            vec![2],
            entries.into_iter().map(Value::int).collect(),
            TypeExpr::Int,
        )
    };
    for (left, right, expected) in [([1, 0], [0, 1], 2.0), ([0, 1], [1, 0], -2.0)] {
        let bracket = call(
            &session,
            "num_native_lie_bracket",
            &[lie.payload[0].clone(), coordinate(left), coordinate(right)],
        );
        let Value::Tensor(bracket) = bracket else {
            panic!("expected exact bracket coordinates")
        };
        assert_eq!(
            bracket.data.iter().map(number).collect::<Vec<_>>(),
            [0.0, expected]
        );
    }
    let Value::Variant(not_closed) = call(&session, "num_native_lie_from_matrices", &[basis(e, f)])
    else {
        panic!("expected an optional Lie algebra for a non-closed basis")
    };
    assert_eq!(not_closed.tag, "None");
}

#[test]
fn combined_ricci_restores_tensor_and_scalar_dimensions_and_propagates_budget_failure() {
    let session = Session::new();
    let mut ctx = CallCtx::new();
    let length = dim("Length");
    let metric_dimension = dim("Temperature");
    let metric_element = TypeExpr::Float.with_dim(metric_dimension.clone());
    let metric = ctx
        .add_callback(Arc::new(NativeCallback::new(
            {
                let metric_element = metric_element.clone();
                move |args: &[Value]| {
                    let Value::Tensor(point) = &args[0] else {
                        panic!("metric point");
                    };
                    let sine = (number(&point.data[0]) / 2.0).sin();
                    Ok(tensor(
                        vec![2, 2],
                        [3.0, 0.0, 0.0, 3.0 * sine * sine]
                            .map(Value::Float)
                            .to_vec(),
                        metric_element.clone(),
                    ))
                }
            },
            Signature::new(
                vec![TypeExpr::Tensor {
                    rank: NatExpr::Lit(1),
                    shape: None,
                    elem: Box::new(TypeExpr::Float.with_dim(length.clone())),
                }],
                TypeExpr::Tensor {
                    rank: NatExpr::Lit(2),
                    shape: None,
                    elem: Box::new(metric_element),
                },
            ),
        )))
        .unwrap();
    let run = |budget| {
        session.call_value(
            "num_native_geometry_ricci_curvature",
            &[
                Value::Function(metric),
                floats(&[std::f64::consts::PI, 0.4], length.clone()),
                scalar(2.0, length.clone()),
                scalar(3.0, metric_dimension.clone()),
                Value::int(budget),
            ],
            &ctx,
        )
    };
    let Value::Record(out) = run(128).unwrap() else {
        panic!("Ricci record");
    };
    assert!(run(1).unwrap_err().to_string().contains("budget"));
    let get = |name| &out.fields.iter().find(|(key, _)| key == name).unwrap().1;
    let Value::Tensor(tensor) = get("tensor") else {
        panic!("covariant tensor");
    };
    let inverse_square = length.pow(Exp::Rat(Rational::int(-2)));
    assert_eq!(tensor.shape, [2, 2]);
    assert_eq!(
        tensor.elem,
        TypeExpr::Float.with_dim(inverse_square.clone())
    );
    for (entry, expected) in tensor.data.iter().zip([0.25, 0.0, 0.0, 0.25]) {
        assert!((number(entry) - expected).abs() < 1e-7);
    }
    let Value::Dyn(value) = get("scalar") else {
        panic!("measured scalar curvature");
    };
    assert_eq!(
        value.dim,
        DimExpr::Mul(vec![
            inverse_square,
            metric_dimension.pow(Exp::Rat(Rational::int(-1)))
        ])
    );
    assert!((number(&value.value) - 1.0 / 6.0).abs() < 1e-7);
}
fn tensor(
    shape: Vec<usize>,
    data: Vec<Value>,
    elem: TypeExpr,
) -> Value {
    Value::Tensor(TensorValue { shape, data, elem })
}
fn ints(data: &[i64]) -> Value {
    tensor(
        vec![data.len()],
        data.iter().copied().map(Value::int).collect(),
        TypeExpr::Int,
    )
}

#[test]
fn fresnel_interval_accepts_improper_limits_and_resolves_adjacent_bounds() {
    let session = Session::new();
    let run = |a, b| {
        let result = call(
            &session,
            "num_fresnel_sc_interval",
            &[Value::Float(a), Value::Float(b)],
        );
        let Value::Record(result) = result else {
            panic!("Fresnel result");
        };
        let get = |name| number(&result.fields.iter().find(|(key, _)| key == name).unwrap().1);
        (get("item0"), get("item1"))
    };
    let near = run(1.0, 1.0_f64.next_up());
    assert_eq!(near.0, f64::EPSILON);
    assert!(near.1 < 0.0);
    assert_eq!(run(1.0_f64.next_up(), 1.0), (-near.0, -near.1));
    assert_eq!(run(f64::NEG_INFINITY, f64::INFINITY), (1.0, 1.0));
    assert_eq!(run(f64::INFINITY, f64::INFINITY), (0.0, 0.0));
    assert!(
        session
            .call_value(
                "num_fresnel_sc_interval",
                &[Value::Float(f64::NAN), Value::Float(1.0)],
                &CallCtx::new(),
            )
            .is_err()
    );
}

#[test]
fn elliptic_pi_plan_reuses_periods_and_keeps_batch_order() {
    let session = Session::new();
    let plan = call(
        &session,
        "num_special_elliptic_pi_plan_new",
        &[Value::Float(0.0), Value::Float(0.0)],
    );
    let complete = call(
        &session,
        "num_special_elliptic_pi_plan_complete",
        std::slice::from_ref(&plan),
    );
    assert!((number(&complete) - std::f64::consts::FRAC_PI_2).abs() < 1e-14);
    let amplitudes = [std::f64::consts::PI, 0.0, -std::f64::consts::PI, 0.25];
    let batch = call(
        &session,
        "num_special_elliptic_pi_plan_eval_many",
        &[plan.clone(), floats(&amplitudes, DimExpr::Dimless)],
    );
    let Value::Tensor(batch) = batch else {
        panic!("elliptic integral tensor");
    };
    assert_eq!(batch.shape, [4]);
    for (value, amplitude) in batch.data.iter().zip(amplitudes) {
        assert!((number(value) - amplitude).abs() < 1e-14);
        assert_eq!(
            value,
            &call(
                &session,
                "num_special_elliptic_pi_plan_eval",
                &[plan.clone(), Value::Float(amplitude)],
            )
        );
    }
    assert_eq!(
        call(
            &session,
            "num_special_elliptic_pi_plan_eval_many",
            &[plan, floats(&[], DimExpr::Dimless)],
        ),
        Value::Tensor(TensorValue {
            shape: vec![0],
            elem: TypeExpr::Float,
            data: vec![]
        })
    );
    assert!(
        session
            .call_value(
                "num_special_elliptic_pi_plan_new",
                &[Value::Float(1.0), Value::Float(0.0)],
                &CallCtx::new(),
            )
            .is_err()
    );
}

#[test]
fn gf256_batch_scale_and_mul_add_preserve_the_field_and_byte_contract() {
    let session = Session::new();
    let src = ints(&[0, 1, 128, 255]);
    let dst = ints(&[1, 3, 0, 255]);
    let scale = "num_finite_field_gf256_scale";
    let add = "num_finite_field_gf256_mul_add";
    assert_eq!(
        call(&session, scale, &[Value::int(2), src.clone()]),
        ints(&[0, 2, 29, 227]),
    );
    assert_eq!(
        call(&session, add, &[Value::int(2), src.clone(), dst.clone()]),
        ints(&[1, 1, 29, 28]),
    );
    assert_eq!(
        call(&session, add, &[Value::int(0), src.clone(), dst.clone()]),
        dst,
    );
    let ctx = CallCtx::new();
    assert!(
        session
            .call_value(add, &[Value::int(1), src.clone(), ints(&[0])], &ctx)
            .is_err()
    );
    assert!(
        session
            .call_value(scale, &[Value::int(256), src], &ctx)
            .is_err()
    );
    assert!(
        session
            .call_value(scale, &[Value::int(1), ints(&[0, -1])], &ctx)
            .is_err()
    );
    assert_eq!(
        call(&session, scale, &[Value::int(255), ints(&[])]),
        ints(&[])
    );
}
fn floats(
    data: &[f64],
    d: DimExpr,
) -> Value {
    tensor(
        vec![data.len()],
        data.iter().copied().map(Value::Float).collect(),
        TypeExpr::Float.with_dim(d),
    )
}
fn some(x: &Value) -> &Value {
    let Value::Variant(v) = x else {
        panic!("option")
    };
    assert_eq!(v.tag, "Some");
    &v.payload[0]
}

#[test]
fn native_weighted_graph_operators_keep_dimensions_and_manifest_shapes() {
    let s = Session::new();
    let length = dim("Length");
    let endpoint_matrix = tensor(
        vec![4, 2],
        [0, 1, 0, 2, 1, 3, 2, 3]
            .into_iter()
            .map(Value::int)
            .collect(),
        TypeExpr::Int,
    );
    let graph = call(
        &s,
        "num_native_graph_from_edges",
        &[
            Value::int(5),
            endpoint_matrix.clone(),
            floats(&[2.0, 5.0, 4.0, 2.0], length.clone()),
        ],
    );
    let handle = std::slice::from_ref(&graph);
    assert_eq!(
        call(&s, "num_native_graph_num_edges", handle),
        Value::int(4)
    );
    assert_eq!(
        call(&s, "num_native_graph_weakly_connected_components", handle),
        ints(&[0, 0, 0, 0, 1])
    );

    let ordered = |name: &str| {
        let Value::Tensor(t) = call(&s, name, handle) else {
            panic!("{name} should return a vertex vector");
        };
        t.data
            .iter()
            .map(|v| v.as_f64().unwrap() as usize)
            .collect::<Vec<_>>()
    };
    let topological = ordered("num_native_graph_topological_sort");
    let scc = ordered("num_native_graph_strongly_connected_components");
    for (u, v) in [(0, 1), (0, 2), (1, 3), (2, 3)] {
        let position = |vertex| topological.iter().position(|&w| w == vertex).unwrap();
        assert!(position(u) < position(v));
        assert!(scc[u] < scc[v]);
    }

    let manifest = s.manifest();
    for name in [
        "num_native_graph_from_edges",
        "num_native_graph_num_edges",
        "num_native_graph_topological_sort",
        "num_native_graph_strongly_connected_components",
        "num_native_graph_weakly_connected_components",
        "num_native_graph_condensation",
        "num_native_graph_undirected_cuts",
        "num_native_graph_multi_source_bfs",
        "num_native_graph_bellman_ford",
        "num_native_graph_dag_shortest_paths",
        "num_native_graph_maximum_flow",
    ] {
        let op = manifest.op(name).expect(name);
        assert_eq!(op.overloads.len(), 1);
        assert!(op.rust_path.is_some());
        assert!(!op.doc.is_empty());
        assert_eq!(
            op.signature_source,
            rssn::kernel::manifest::SignatureSource::Declared
        );
    }

    let Value::Record(condensation) = call(&s, "num_native_graph_condensation", handle) else {
        panic!("expected a condensation record");
    };
    assert_eq!(condensation.fields[0].0, "component_of");
    assert_eq!(condensation.fields[1].0, "arcs");
    let Value::Tensor(arcs) = &condensation.fields[1].1 else {
        panic!("expected a condensation arc matrix");
    };
    assert_eq!(arcs.shape, [4, 2]);
    assert!(
        arcs.data
            .as_chunks::<2>()
            .0
            .iter()
            .all(|edge| number(&edge[0]) < number(&edge[1]))
    );

    let endpoints = tensor(
        vec![5, 2],
        [0, 1, 0, 1, 1, 2, 1, 3, 3, 3]
            .into_iter()
            .map(Value::int)
            .collect(),
        TypeExpr::Int,
    );
    let Value::Record(cuts) = call(
        &s,
        "num_native_graph_undirected_cuts",
        &[Value::int(4), endpoints.clone()],
    ) else {
        panic!("expected a cut record");
    };
    assert_eq!(cuts.fields[0].0, "bridges");
    assert_eq!(
        cuts.fields[1],
        ("articulation_points".to_owned(), ints(&[1]))
    );
    let Value::Tensor(bridges) = &cuts.fields[0].1 else {
        panic!("expected a bridge matrix");
    };
    assert_eq!(bridges.shape, [2, 2]);
    assert_eq!(bridges.data, [1, 2, 1, 3].map(Value::int));
    assert!(
        s.call_value(
            "num_native_graph_undirected_cuts",
            &[Value::int(3), endpoints],
            &CallCtx::new()
        )
        .is_err()
    );

    let Value::Tensor(unweighted) = call(
        &s,
        "num_native_graph_multi_source_bfs",
        &[graph.clone(), ints(&[2, 1, 2])],
    ) else {
        panic!("expected an option tensor of distances");
    };
    assert_eq!(unweighted.shape, [5]);
    assert!(matches!(&unweighted.data[0], Value::Variant(v) if v.tag == "None"));
    assert_eq!(some(&unweighted.data[1]), &Value::int(0));
    assert_eq!(some(&unweighted.data[2]), &Value::int(0));
    assert_eq!(some(&unweighted.data[3]), &Value::int(1));
    assert!(matches!(&unweighted.data[4], Value::Variant(v) if v.tag == "None"));
    assert!(
        s.call_value(
            "num_native_graph_multi_source_bfs",
            &[graph.clone(), ints(&[5])],
            &CallCtx::new()
        )
        .is_err()
    );

    for name in [
        "num_native_graph_bellman_ford",
        "num_native_graph_dag_shortest_paths",
    ] {
        let schema = &manifest.op(name).unwrap().overloads[0].ret;
        assert!(matches!(schema, TypeExpr::Record(fields) if fields.len() == 2));
        let Value::Record(result) = call(&s, name, &[graph.clone(), Value::int(0)]) else {
            panic!("{name} should return a path record");
        };
        let Value::Tensor(distances) = &result.fields[0].1 else {
            panic!("expected distance vector");
        };
        assert_eq!(distances.elem, TypeExpr::Float.with_dim(length.clone()));
        assert_eq!(
            distances.data.iter().map(number).collect::<Vec<_>>(),
            [0.0, 2.0, 5.0, 6.0, f64::INFINITY]
        );
        let Value::Tensor(previous) = &result.fields[1].1 else {
            panic!("expected predecessor vector");
        };
        assert_eq!(some(&previous.data[3]), &Value::int(1));
        assert!(matches!(&previous.data[4], Value::Variant(v) if v.tag == "None"));
    }

    let Value::Record(flow) = call(
        &s,
        "num_native_graph_maximum_flow",
        &[graph.clone(), Value::int(0), Value::int(3)],
    ) else {
        panic!("expected a flow record");
    };
    assert_eq!(flow.fields[0].1, scalar(4.0, length.clone()));
    let Value::Tensor(source_side) = &flow.fields[1].1 else {
        panic!("expected a cut mask");
    };
    assert_eq!(
        source_side.data,
        [true, false, true, false, false].map(Value::Bool)
    );
    assert!(
        s.call_value(
            "num_native_graph_from_edges",
            &[
                Value::int(5),
                endpoint_matrix,
                floats(&[2.0], length.clone())
            ],
            &CallCtx::new()
        )
        .is_err()
    );
    assert!(
        s.call_value(
            "num_native_graph_add_edge",
            &[
                graph.clone(),
                Value::int(0),
                Value::int(4),
                Value::Float(1.0)
            ],
            &CallCtx::new()
        )
        .is_err()
    );
    assert_eq!(
        call(&s, "num_native_graph_num_edges", handle),
        Value::int(4)
    );
}

#[test]
fn native_group_and_exact_lie_objects_are_actual_checked_algorithms() {
    let s = Session::new();
    let generators = tensor(
        vec![2, 3],
        vec![1, 0, 2, 1, 2, 0].into_iter().map(Value::int).collect(),
        TypeExpr::Int,
    );
    let g = call(&s, "num_native_group_new", &[Value::int(3), generators]);
    assert_eq!(
        call(&s, "num_native_group_order", std::slice::from_ref(&g)),
        Value::int(6)
    );
    let center = call(&s, "num_native_group_center", std::slice::from_ref(&g));
    assert_eq!(call(&s, "num_native_group_order", &[center]), Value::int(1));
    assert!(
        s.call_value(
            "num_native_group_sylow_subgroup",
            &[g, Value::int(1), Value::int(1), Value::int(10)],
            &CallCtx::new()
        )
        .is_err()
    );
    let c = call(
        &s,
        "num_native_lie_parse_type",
        &[Value::String("A2".into())],
    );
    let r = call(&s, "num_native_root_system_new", &[some(&c).clone()]);
    let r = some(&r).clone();
    assert_eq!(
        call(
            &s,
            "num_native_root_system_dimension",
            std::slice::from_ref(&r)
        ),
        Value::int(8)
    );
    assert_eq!(
        some(&call(
            &s,
            "num_native_root_system_weyl_dim",
            &[r.clone(), ints(&[1, 0])]
        )),
        &Value::int(3)
    );
    assert_eq!(
        call(&s, "num_native_root_system_weyl_order", &[r]),
        Value::int(6)
    );
}
#[test]
fn native_geometry_frames_and_ball_witnesses_retain_physical_dimensions() {
    let s = Session::new();
    let m = |x| scalar(x, dim("Length"));
    let t = |x| scalar(x, dim("Time"));
    let p = call(&s, "num_native_graphics_point3_new", &[m(1.), m(2.), m(3.)]);
    let transform = call(
        &s,
        "num_native_graphics_translation_matrix",
        &[m(4.), m(0.), m(0.), m(2.)],
    );
    let p = call(&s, "num_native_graphics_transform_point", &[transform, p]);
    let v = call(&s, "num_native_graphics_point3_components", &[p]);
    let Value::Tensor(v) = v else {
        panic!("tensor")
    };
    assert_eq!(v.elem, TypeExpr::Float.with_dim(dim("Length")));
    assert_eq!(
        v.data.iter().map(number).collect::<Vec<_>>(),
        vec![5., 2., 3.]
    );
    let b = call(&s, "num_native_cball_square", &[m(2.), m(0.), m(0.25)]);
    assert_eq!(
        call(
            &s,
            "num_native_cball_contains",
            &[b.clone(), scalar(2., dim("Length"))]
        ),
        Value::Bool(true)
    );
    assert!(
        s.call_value(
            "num_native_cball_contains",
            &[b.clone(), t(2.)],
            &CallCtx::new()
        )
        .is_err()
    );
    let squared = call(&s, "num_native_cball_powi", &[b, Value::int(2)]);
    let width = call(&s, "num_native_cball_half_width", &[squared]);
    let Value::Dyn(witness) = &width else {
        panic!("measured ball width must retain its dimension witness")
    };
    let expected = TypeExpr::Float.with_dim(dim("Length").pow(rssn::kernel::signature::Exp::Rat(
        rssn::kernel::signature::Rational::int(2),
    )));
    assert!(
        rssn::kernel::overload::matches_type(
            &expected,
            &TypeExpr::Float.with_dim(witness.dim.clone())
        )
        .accepted()
    );
    assert!((number(&width) - 1.125).abs() < 1e-12);
}
#[test]
fn native_polynomial_certification_and_repeated_roots_are_scale_consistent() {
    let s = Session::new();
    let system = call(
        &s,
        "num_native_polynomial_system_new",
        &[
            floats(&[1., -2.], DimExpr::Dimless),
            tensor(
                vec![2, 1],
                vec![Value::int(2), Value::int(0)],
                TypeExpr::Int,
            ),
            ints(&[0, 2]),
            scalar(2., dim("Length")),
            scalar(3., dim("Time")),
        ],
    );
    let p = floats(&[2. * 2f64.sqrt()], dim("Length"));
    let proof = call(&s, "num_native_homotopy_certify_zero", &[system, p]);
    assert!(matches!(proof,Value::Variant(ref v)if v.tag=="Some"));
    let poly = call(
        &s,
        "num_native_polynomial_normalized_new",
        &[floats(&[1., 0., -4., 0., 4.], DimExpr::Dimless)],
    );
    let roots = call(
        &s,
        "num_native_real_roots_find_roots",
        &[
            poly,
            scalar(1e-10, dim("Length")),
            scalar(1., dim("Length")),
        ],
    );
    let Value::Tensor(roots) = roots else {
        panic!("roots")
    };
    assert_eq!(roots.data.len(), 2);
    assert!((number(&roots.data[0]) + 2f64.sqrt()).abs() < 1e-9);
    assert!((number(&roots.data[1]) - 2f64.sqrt()).abs() < 1e-9);
}
#[test]
fn native_parallel_callbacks_run_on_owner_thread_preserve_units_and_failures() {
    let s = Session::new();
    let mut ctx = CallCtx::new();
    let owner = std::thread::current().id();
    let calls = Arc::new(AtomicUsize::new(0));
    let seen = calls.clone();
    let h = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |args: &[Value]| {
                assert_eq!(std::thread::current().id(), owner);
                seen.fetch_add(1, Ordering::Relaxed);
                Ok(scalar(args[0].as_f64().unwrap() + 1., dim("Length")))
            },
            Signature::new(vec![TypeExpr::Int], TypeExpr::Float.with_dim(dim("Length"))),
        )))
        .unwrap();
    let out = s
        .call_value(
            "num_native_par_map_indexed",
            &[
                Value::int(100),
                Value::int(2),
                Value::Function(h),
                Value::int(100),
            ],
            &ctx,
        )
        .unwrap();
    let Value::Tensor(out) = out else {
        panic!("tensor")
    };
    assert_eq!(out.data.len(), 100);
    assert_eq!(number(&out.data[99]), 100.);
    assert_eq!(calls.load(Ordering::Relaxed), 100);
    let err = s
        .call_value(
            "num_native_par_sum_fixed",
            &[
                Value::int(100),
                Value::int(4),
                Value::int(2),
                Value::Function(h),
                Value::int(2),
            ],
            &ctx,
        )
        .unwrap_err();
    assert!(err.to_string().contains("budget"));
    let mut failed = CallCtx::new();
    let h = failed
        .add_callback(Arc::new(NativeCallback::new(
            |_: &[Value]| Err(OpError::Host("parallel-sentinel".into())),
            Signature::new(vec![TypeExpr::Int], TypeExpr::Float),
        )))
        .unwrap();
    let err = s
        .call_value(
            "num_native_par_map_indexed",
            &[
                Value::int(100),
                Value::int(1),
                Value::Function(h),
                Value::int(100),
            ],
            &failed,
        )
        .unwrap_err();
    assert!(err.to_string().contains("parallel-sentinel"));
}
#[test]
fn native_metric_and_bootstrap_callbacks_keep_independent_input_output_dimensions() {
    let s = Session::new();
    let mut ctx = CallCtx::new();
    let vector_ty = |d| {
        TypeExpr::Tensor {
            rank: NatExpr::Lit(1),
            shape: None,
            elem: Box::new(TypeExpr::Float.with_dim(d)),
        }
    };
    let statistic = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |args: &[Value]| {
                let Value::Tensor(t) = &args[0] else {
                    panic!("sample")
                };
                Ok(scalar(
                    t.data.iter().map(number).sum::<f64>() / t.data.len() as f64,
                    dim("Time"),
                ))
            },
            Signature::new(
                vec![vector_ty(dim("Length"))],
                TypeExpr::Float.with_dim(dim("Time")),
            ),
        )))
        .unwrap();
    let out = s
        .call_value(
            "num_native_robust_bootstrap_ci",
            &[
                floats(&[2., 2., 2.], dim("Length")),
                Value::Function(statistic),
                Value::int(8),
                Value::Float(0.1),
                Value::int(1),
                Value::int(8),
            ],
            &ctx,
        )
        .unwrap();
    let Value::Tensor(out) = out else {
        panic!("CI")
    };
    assert_eq!(out.elem, TypeExpr::Float.with_dim(dim("Time")));
    assert_eq!(number(&out.data[0]), 2.);
    let metric = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |_: &[Value]| {
                Ok(tensor(
                    vec![2, 2],
                    vec![
                        Value::Float(1.),
                        Value::Float(0.),
                        Value::Float(0.),
                        Value::Float(1.),
                    ],
                    TypeExpr::Float,
                ))
            },
            Signature::new(
                vec![vector_ty(dim("Length"))],
                TypeExpr::Tensor {
                    rank: NatExpr::Lit(2),
                    shape: None,
                    elem: Box::new(TypeExpr::Float),
                },
            ),
        )))
        .unwrap();
    let out = s
        .call_value(
            "num_native_geometry_ricci_scalar",
            &[
                Value::Function(metric),
                floats(&[1., 2.], dim("Length")),
                scalar(1., dim("Length")),
                Value::Float(1.),
                Value::int(1000),
            ],
            &ctx,
        )
        .unwrap();
    assert!(number(&out).abs() < 1e-9);
}
#[test]
fn native_rng_and_online_covariance_are_mutable_session_objects() {
    let s = Session::new();
    let a = call(&s, "num_native_rng_new", &[Value::int(42)]);
    let b = call(&s, "num_native_rng_new", &[Value::int(42)]);
    assert_eq!(
        call(&s, "num_native_rng_next_u64", std::slice::from_ref(&a)),
        call(&s, "num_native_rng_next_u64", &[b])
    );
    assert!(
        Session::new()
            .call_value("num_native_rng_uniform", &[a], &CallCtx::new())
            .is_err()
    );
    let c = call(
        &s,
        "num_native_online_covariance_new",
        &[Value::int(2), scalar(1., dim("Length"))],
    );
    for p in [[1., 3.], [3., 7.]] {
        call(
            &s,
            "num_native_online_covariance_push",
            &[c.clone(), floats(&p, dim("Length"))],
        );
    }
    assert_eq!(
        call(
            &s,
            "num_native_online_covariance_count",
            std::slice::from_ref(&c)
        ),
        Value::int(2)
    );
    let cov = call(&s, "num_native_online_covariance_sample_covariance", &[c]);
    let Value::Tensor(cov) = some(&cov) else {
        panic!("covariance")
    };
    assert_eq!(
        cov.data.iter().map(number).collect::<Vec<_>>(),
        vec![2., 4., 4., 8.]
    );
}

#[test]
fn online_covariance_batch_preserves_units_and_returns_dimensionless_correlation() {
    let s = Session::new();
    let length = dim("Length");
    let covariance = call(
        &s,
        "num_native_online_covariance_new",
        &[Value::int(2), scalar(1., length.clone())],
    );
    let data = [1., 3., 3., 7., 5., 11.].map(Value::Float).to_vec();
    for (shape, dimension) in [(vec![3, 2], dim("Time")), (vec![2, 3], length.clone())] {
        let invalid = tensor(shape, data.clone(), TypeExpr::Float.with_dim(dimension));
        assert!(
            s.call_value(
                "num_native_online_covariance_push_batch",
                &[covariance.clone(), invalid],
                &CallCtx::new(),
            )
            .is_err()
        );
    }
    let valid = tensor(vec![3, 2], data, TypeExpr::Float.with_dim(length));
    assert_eq!(
        call(
            &s,
            "num_native_online_covariance_push_batch",
            &[covariance.clone(), valid]
        ),
        Value::Bool(true)
    );
    let handle = std::slice::from_ref(&covariance);
    assert_eq!(
        call(&s, "num_native_online_covariance_count", handle),
        Value::int(3)
    );
    let correlation = call(&s, "num_native_online_covariance_correlation", handle);
    let Value::Tensor(result) = some(&correlation) else {
        panic!("correlation matrix");
    };
    assert_eq!(result.shape, vec![2, 2]);
    assert_eq!(result.elem, TypeExpr::Float);
    assert!(
        result
            .data
            .iter()
            .all(|value| (number(value) - 1.).abs() < 1e-15)
    );
}

#[test]
fn online_covariance_frequency_batch_retains_units_and_is_atomic() {
    let s = Session::new();
    let length = dim("Length");
    let squared_length = length.clone().pow(Exp::Rat(Rational::int(2)));
    let covariance = call(
        &s,
        "num_native_online_covariance_new",
        &[Value::int(2), scalar(1., length.clone())],
    );
    let observations = tensor(
        vec![2, 2],
        [1., 3., 5., 11.].map(Value::Float).to_vec(),
        TypeExpr::Float.with_dim(length),
    );
    let op = "num_native_online_covariance_push_frequency_batch";
    assert_eq!(
        call(
            &s,
            op,
            &[covariance.clone(), observations.clone(), ints(&[3, 1])]
        ),
        Value::Bool(true)
    );
    for invalid in [
        ints(&[1]),
        tensor(vec![2], vec![Value::int(-1), Value::int(1)], TypeExpr::Int),
        tensor(
            vec![2],
            vec![Value::int(1), Value::int(1)],
            TypeExpr::Int.with_dim(dim("Time")),
        ),
        tensor(
            vec![2],
            vec![Value::int(1_u64 << 53), Value::int(0)],
            TypeExpr::Int,
        ),
    ] {
        assert!(
            s.call_value(
                op,
                &[covariance.clone(), observations.clone(), invalid],
                &CallCtx::new()
            )
            .is_err()
        );
    }
    let handle = std::slice::from_ref(&covariance);
    assert_eq!(
        call(&s, "num_native_online_covariance_count", handle),
        Value::int(4)
    );
    let result = call(&s, "num_native_online_covariance_sample_covariance", handle);
    let Value::Tensor(result) = some(&result) else {
        panic!("sample covariance matrix");
    };
    assert_eq!(result.shape, [2, 2]);
    assert!(
        rssn::kernel::overload::matches_type(
            &TypeExpr::Float.with_dim(squared_length),
            &result.elem,
        )
        .accepted()
    );
    assert_eq!(
        result.data.iter().map(number).collect::<Vec<_>>(),
        [4., 8., 8., 16.]
    );
}

#[test]
fn ordering_metrics_accept_ragged_offsets_and_checked_permutations() {
    let session = Session::new();
    let op = "num_native_ordering_metrics";
    let offsets = ints(&[0, 2, 3, 3, 3]);
    let neighbors = ints(&[3, 3, 3]);
    for (order, bandwidth, profile) in [([0, 1, 2, 3], 3, 3), ([0, 3, 1, 2], 1, 2)] {
        let Value::Record(result) = call(
            &session,
            op,
            &[offsets.clone(), neighbors.clone(), ints(&order)],
        ) else {
            panic!("expected an ordering metrics record");
        };
        assert_eq!(
            result.fields,
            [
                ("half_bandwidth".to_string(), Value::int(bandwidth)),
                ("profile".to_string(), Value::int(profile)),
            ]
        );
    }
    let Value::Record(empty) = call(&session, op, &[ints(&[0]), ints(&[]), ints(&[])]) else {
        panic!("expected an empty pattern's ordering metrics record");
    };
    assert!(
        empty
            .fields
            .iter()
            .all(|(_, value)| *value == Value::int(0))
    );
    for (offsets, neighbors, order) in [
        (
            ints(&[0, 3, 2, 3, 3]),
            neighbors.clone(),
            ints(&[0, 1, 2, 3]),
        ),
        (offsets.clone(), ints(&[3, 4, 3]), ints(&[0, 1, 2, 3])),
        (offsets, neighbors, ints(&[0, 1, 2, 2])),
    ] {
        assert!(
            session
                .call_value(op, &[offsets, neighbors, order], &CallCtx::new())
                .is_err()
        );
    }
}

#[test]
fn general_multivector_inverse_covers_nonhomogeneous_and_tiny_elements() {
    use rssn::kernels::geometric_algebra::Multivector3D;
    let a = Multivector3D {
        s: 2.,
        v1: 1.,
        ..Default::default()
    };
    let inverse = a.inv().unwrap();
    assert!((inverse.s - 2. / 3.).abs() < 1e-15);
    assert!((inverse.v1 + 1. / 3.).abs() < 1e-15);
    assert!(
        Multivector3D {
            s: 1e-300,
            ..Default::default()
        }
        .inv()
        .is_some()
    );
    assert!(
        Multivector3D {
            s: 1.,
            v1: 1.,
            ..Default::default()
        }
        .inv()
        .is_none()
    );
}

#[test]
fn native_multivector_central_product_combines_both_physical_dimensions() {
    let session = Session::new();
    let length = dim("Length");
    let time = dim("Time");
    let coefficients =
        [2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0].map(|value| scalar(value, length.clone()));
    let multivector = call(&session, "num_native_multivector_new", &coefficients);
    let name = "num_native_multivector_mul_central";
    let result = call(
        &session,
        name,
        &[
            multivector.clone(),
            scalar(4.0, time.clone()),
            scalar(5.0, time.clone()),
        ],
    );
    let squared = call(&session, "num_native_multivector_norm_sq", &[result]);
    let dimension = DimExpr::Mul(vec![length.clone().pow(2), time.clone().pow(2)]);
    assert!(
        rssn::kernel::overload::matches_type(
            &TypeExpr::Float.with_dim(dimension),
            &squared.type_expr()
        )
        .accepted()
    );
    assert_eq!(number(&squared), 533.0);
    assert!(
        session
            .call_value(
                name,
                &[multivector, scalar(4.0, time), scalar(5.0, length)],
                &CallCtx::new()
            )
            .is_err()
    );
}

#[test]
fn native_double_double_norms_preserve_the_low_component_and_dimension() {
    let s = Session::new();
    let length = dim("Length");
    let measured = scalar(1.0, length.clone());
    let x = call(
        &s,
        "num_native_dd_from_f64",
        std::slice::from_ref(&measured),
    );
    let hypot = call(&s, "num_native_dd_hypot", &[x.clone(), x.clone()]);
    let complex = call(&s, "num_native_cdd_new", &[measured.clone(), measured]);
    let modulus = call(&s, "num_native_cdd_norm_dd", &[complex]);

    let parts = call(&s, "num_native_dd_parts", &[hypot]);
    assert_eq!(parts, call(&s, "num_native_dd_parts", &[modulus]));
    let Value::Record(parts) = parts else {
        panic!("expected double-double components")
    };
    let high = parts.get("item0").unwrap();
    let low = parts.get("item1").unwrap();
    let measured_float = TypeExpr::Float.with_dim(length);
    for value in [high, low] {
        assert!(
            rssn::kernel::overload::matches_type(&measured_float, &value.type_expr()).accepted()
        );
    }
    assert_eq!(number(high), std::f64::consts::SQRT_2);
    assert!((number(low) + 9.667_293_313_452_913e-17).abs() <= 2e-31);

    let wrong = call(&s, "num_native_dd_from_f64", &[scalar(1.0, dim("Time"))]);
    let foreign = Session::new();
    for other in [
        wrong,
        call(&foreign, "num_native_dd_from_f64", &[Value::Float(1.0)]),
    ] {
        assert!(
            s.call_value("num_native_dd_hypot", &[x.clone(), other], &CallCtx::new())
                .is_err()
        );
    }
}

#[test]
fn native_complex_double_double_reciprocal_handles_extreme_scales_and_inverse_units() {
    let s = Session::new();
    let length = dim("Length");
    for (a, b, expected) in [(3e200, 4e200, 2e-201), (3e-200, 4e-200, 2e199)] {
        let z = call(
            &s,
            "num_native_cdd_new",
            &[scalar(a, length.clone()), scalar(b, length.clone())],
        );
        let inverse = call(&s, "num_native_cdd_recip", &[z]);
        let modulus = call(&s, "num_native_cdd_norm", &[inverse]);
        let ty = TypeExpr::Float.with_dim(DimExpr::Inverse(Box::new(length.clone())));
        assert!(rssn::kernel::overload::matches_type(&ty, &modulus.type_expr()).accepted());
        assert!((number(&modulus) / expected - 1.0).abs() < 1e-14);
    }
    let zero = call(
        &s,
        "num_native_cdd_new",
        &[scalar(0.0, length.clone()), scalar(0.0, length)],
    );
    assert!(
        s.call_value("num_native_cdd_recip", &[zero], &CallCtx::new())
            .is_err()
    );
}

#[test]
fn native_polygon_winding_has_a_measured_complex_input_and_integer_result() {
    let s = Session::new();
    let name = "num_complex_winding_number";
    let length = dim("Length");
    let elem = TypeExpr::Complex.with_dim(length);
    let complex = |re, im| Value::Complex(ComplexValue::Float { re, im });
    let r = 1e-200;
    let vertices = vec![
        complex(r, 0.0),
        complex(0.0, r),
        complex(-r, 0.0),
        complex(0.0, -r),
    ];
    assert_eq!(
        call(&s, name, &[tensor(vec![4], vertices.clone(), elem.clone())]),
        Value::int(1)
    );
    let mut mixed = vertices;
    mixed[1] = Value::Dyn(DynValue {
        dim: dim("Time"),
        value: Box::new(complex(0.0, r)),
    });
    for invalid in [
        tensor(vec![4], mixed, elem.clone()),
        tensor(
            vec![2],
            vec![complex(r, 0.0), complex(-r, 0.0)],
            elem.clone(),
        ),
        tensor(vec![0], vec![], elem),
    ] {
        assert!(s.call_value(name, &[invalid], &CallCtx::new()).is_err());
    }
}

#[test]
fn native_quadratic_sieve_returns_a_sorted_exact_factor_pair() {
    let s = Session::new();
    let limits = call(&s, "num_native_qs_limits_for_digits", &[Value::int(2)]);
    let name = "num_native_qs_factor_pair";
    let Value::Variant(answer) = call(&s, name, &[Value::int(77), limits.clone()]) else {
        panic!("expected an optional factor pair")
    };
    assert_eq!(answer.tag, "Some");
    let Value::Record(pair) = &answer.payload[0] else {
        panic!("expected two ordered factors")
    };
    assert_eq!(pair.get("item0"), Some(&Value::int(7)));
    assert_eq!(pair.get("item1"), Some(&Value::int(11)));

    let Value::Variant(missing) = call(&s, name, &[Value::int(79), limits]) else {
        panic!("expected a missing factor pair for a prime")
    };
    assert_eq!(missing.tag, "None");
}

#[test]
fn prepared_reed_solomon_encoders_keep_exact_bytes_and_check_block_extents() {
    let s = Session::new();
    let parity = Value::int(4);
    let encoder = call(
        &s,
        "num_native_rs_encoder_new",
        std::slice::from_ref(&parity),
    );
    assert_eq!(
        call(
            &s,
            "num_native_rs_encoder_parity_symbols",
            std::slice::from_ref(&encoder)
        ),
        parity
    );
    let message = ints(&[0, 127, 255, 37, 11]);
    let expected = call(
        &s,
        "num_error_correction_reed_solomon_encode",
        &[message.clone(), parity],
    );
    let encode = "num_native_rs_encoder_encode";
    for _ in 0..3 {
        assert_eq!(
            call(&s, encode, &[encoder.clone(), message.clone()]),
            expected
        );
    }
    for invalid in [ints(&[256]), ints(&[0; 252])] {
        assert!(
            s.call_value(encode, &[encoder.clone(), invalid], &CallCtx::new())
                .is_err()
        );
    }
    assert!(
        s.call_value(
            "num_native_rs_encoder_new",
            &[Value::int(256)],
            &CallCtx::new()
        )
        .is_err()
    );
    assert!(
        Session::new()
            .call_value(encode, &[encoder, message], &CallCtx::new())
            .is_err()
    );

    let max = call(&s, "num_native_rs_encoder_new", &[Value::int(255)]);
    assert_eq!(call(&s, encode, &[max, ints(&[])]), ints(&[0; 255]));
}

#[test]
fn symbolic_reed_solomon_error_count_checks_block_extents() {
    let s = Session::new();
    let config = rssn::api::Config::new();
    for (length, parity) in [(256, 1), (4, 5)] {
        let symbols = vec!["0"; length].join(", ");
        let input = format!("rs_error_count(list({symbols}), {parity})");
        let term = s.parse(&input).unwrap();
        assert_eq!(s.compute(term, &config).unwrap().term.to_string(), input);
    }
    let valid = s.parse("rs_error_count(list(0, 0, 0, 0), 1)").unwrap();
    assert_eq!(s.compute(valid, &config).unwrap().term.to_string(), "0");
}

#[test]
fn native_ordered_selections_preserve_exact_word_boundaries() {
    use num_bigint::BigInt;
    let s = Session::new();
    let name = "num_native_number_exact_falling_factorial";
    let expected = BigInt::from(i64::MAX) * (i64::MAX - 1) * (i64::MAX - 2);
    assert_eq!(
        call(&s, name, &[Value::int(i64::MAX), Value::int(3)]),
        Value::Int(expected)
    );
    assert_eq!(
        call(&s, name, &[Value::int(2), Value::int(3)]),
        Value::int(0)
    );
    assert!(
        s.call_value(name, &[Value::int(5000), Value::int(4097)], &CallCtx::new())
            .is_err()
    );
}

#[test]
fn affine_ifs_reuses_maps_and_checks_inputs_at_the_native_boundary() {
    let s = Session::new();
    let maps = tensor(
        vec![2, 6],
        [0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 7.0, -7.0]
            .map(Value::Float)
            .to_vec(),
        TypeExpr::Float,
    );
    let weights = floats(&[f64::MAX, f64::MAX], DimExpr::Dimless);
    let origin = Value::Record(RecordValue {
        fields: vec![
            ("item0".into(), Value::Float(0.0)),
            ("item1".into(), Value::Float(0.0)),
        ],
    });
    let point_tensor = |columns: &[f64]| {
        tensor(
            vec![columns.len() / 2, 2],
            columns.iter().copied().map(Value::Float).collect(),
            TypeExpr::Float,
        )
    };
    let new = "num_native_affine_ifs_new";
    let generate = "num_native_affine_ifs_generate";
    let plan = call(&s, new, &[maps.clone(), weights.clone()]);
    let args = |points: i64, skip: i64, seed: Value| {
        [
            plan.clone(),
            origin.clone(),
            Value::int(points),
            Value::int(skip),
            seed,
        ]
    };
    let historical = point_tensor(&[7.0, -7.0, 1.0, -1.0, 1.0, -1.0]);
    for _ in 0..2 {
        assert_eq!(
            call(&s, generate, &args(3, 1, Value::int(12345))),
            historical
        );
    }
    assert_eq!(
        call(
            &s,
            "num_fractal_geometry_and_chaos_generate_ifs_fractal",
            &[
                maps.clone(),
                weights,
                origin.clone(),
                Value::int(3),
                Value::int(1)
            ],
        ),
        historical
    );
    assert_eq!(
        call(&s, generate, &args(2, 0, Value::int(u64::MAX))),
        point_tensor(&[7.0, -7.0, 7.0, -7.0])
    );
    assert_eq!(
        call(&s, generate, &args(0, 49_999_999, Value::int(0))),
        point_tensor(&[])
    );
    for invalid in [
        args(524_289, 0, Value::int(0)),
        args(1, 50_000_000, Value::int(0)),
        args(1, 0, Value::int(-1)),
    ] {
        assert!(s.call_value(generate, &invalid, &CallCtx::new()).is_err());
    }
    assert!(
        Session::new()
            .call_value(generate, &args(3, 0, Value::int(0)), &CallCtx::new())
            .is_err()
    );
    for invalid in [
        floats(&[0.0, 0.0], DimExpr::Dimless),
        floats(&[1.0], DimExpr::Dimless),
        floats(&[1.0, 1.0], dim("Length")),
    ] {
        assert!(
            s.call_value(new, &[maps.clone(), invalid], &CallCtx::new())
                .is_err()
        );
    }
}
