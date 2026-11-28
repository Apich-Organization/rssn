//! Acoustic boundary and source programs keep their budgets after detachment.
use rssn::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;

fn call(
    s: &Session,
    name: &str,
    args: &[Value],
) -> Result<Value, OpError> {
    s.call_value(name, args, &CallCtx::new())
}

fn invoke(
    s: &Session,
    name: &str,
    args: &[Value],
) -> Value {
    call(s, name, args).unwrap_or_else(|error| panic!("{name}: {error}"))
}

fn release(
    s: &Session,
    object: &Value,
) -> Result<Value, OpError> {
    call(s, "num_backend_release", std::slice::from_ref(object))
}

fn vector(
    values: &[Value],
    elem: TypeExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        data: values.to_vec(),
        elem,
    })
}

fn measured(
    n: f64,
    dim: DimExpr,
) -> Value {
    Value::Dyn(DynValue {
        dim,
        value: Box::new(Value::Float(n)),
    })
}

fn length() -> DimExpr {
    DimExpr::base("Length")
}
fn time() -> DimExpr {
    DimExpr::base("Time")
}

fn pressure() -> DimExpr {
    DimExpr::Mul(vec![DimExpr::base("Mass"), length(), time().pow(-2)])
}

fn reference(
    s: &Session,
    dim: DimExpr,
) -> Value {
    invoke(s, "num_backend_reference", &[measured(1., dim)])
}

fn compile(
    s: &Session,
    expression: &str,
    output: DimExpr,
) -> Value {
    let t = reference(s, time());
    let output = reference(s, output);
    invoke(
        s,
        "num_backend_compile_current",
        &[
            Value::String(expression.into()),
            vector(&[Value::String("t".into())], TypeExpr::String),
            vector(
                &[t],
                TypeExpr::Opaque("rssn/native/backend-scalar/v1".into()),
            ),
            output,
            Value::int(128),
        ],
    )
}

struct Fixture {
    session: Session,
    refs: Value,
    frame: Value,
    uniform_faces: Value,
}

impl Fixture {
    fn new() -> Self {
        let session = Session::new();
        let refs = invoke(
            &session,
            "num_grid_mechanical_references",
            &[
                measured(1., length()),
                measured(1., time()),
                measured(1., DimExpr::base("Mass")),
                measured(1., DimExpr::base("Temperature")),
            ],
        );
        let frame = invoke(
            &session,
            "num_grid_frame_new",
            &[
                measured(1., length()),
                measured(1., time()),
                measured(1., pressure()),
            ],
        );
        let uniform_faces = invoke(
            &session,
            "num_grid_faces_dirichlet_zero",
            std::slice::from_ref(&frame),
        );
        Self {
            session,
            refs,
            frame,
            uniform_faces,
        }
    }

    fn config(
        &self,
        faces: &Value,
    ) -> Result<Value, OpError> {
        call(
            &self.session,
            "num_grid_acoustics_config",
            &[
                self.refs.clone(),
                faces.clone(),
                vector(&[Value::int(4)], TypeExpr::Int),
                vector(&[Value::Float(1.)], TypeExpr::Float),
                Value::Float(1.),
                Value::Float(1.),
                Value::Float(0.5),
                Value::int(0),
                Value::Float(2.),
            ],
        )
    }

    fn source(
        &self,
        solver: &Value,
        program: &Value,
        budget: i64,
    ) -> Result<Value, OpError> {
        call(
            &self.session,
            "num_grid_acoustics_source",
            &[
                solver.clone(),
                Value::int(0),
                program.clone(),
                Value::int(budget),
            ],
        )
    }

    fn faces_with_datum(
        &self,
        base: &Value,
        program: &Value,
        face: i64,
    ) -> (Value, Value, Value) {
        let s = &self.session;
        let datum = invoke(
            s,
            "num_grid_datum_time",
            &[
                self.frame.clone(),
                program.clone(),
                Value::int(16),
                measured(1., pressure()),
            ],
        );
        let boundary = invoke(
            s,
            "num_grid_boundary_from_datum",
            &[
                datum.clone(),
                Value::String("Dirichlet".into()),
                Value::Float(0.),
                Value::Float(0.),
            ],
        );
        let faces = invoke(
            s,
            "num_grid_faces_with",
            &[base.clone(), boundary.clone(), Value::int(face)],
        );
        (faces, boundary, datum)
    }
}

fn field<'a>(
    record: &'a Value,
    name: &str,
) -> &'a Value {
    let Value::Record(record) = record else {
        panic!("expected record, found {record:?}")
    };
    record.get(name).unwrap()
}

fn first_pressure(
    s: &Session,
    solver: &Value,
) -> f64 {
    let state = invoke(s, "num_grid_acoustics_fields", std::slice::from_ref(solver));
    let data = field(field(&state, "pressure"), "normalized");
    let Value::Tensor(tensor) = data else {
        panic!("expected pressure vector")
    };
    tensor.data[0].as_f64().unwrap()
}

#[test]
fn detached_boundary_uses_the_replacement_and_allows_both_programs_to_be_released() {
    let model = Fixture::new();
    let s = &model.session;
    let program_a = compile(s, "1", pressure());
    let program_b = compile(s, "2", pressure());
    let (old_faces, old_boundary, old_datum) =
        model.faces_with_datum(&model.uniform_faces, &program_a, 0);
    let (new_faces, new_boundary, new_datum) = model.faces_with_datum(&old_faces, &program_b, 0);
    // An unused Y face may hold the old datum. A 1D solver keeps only its
    // two active X faces regardless of the history of the user-facing Faces.
    let inactive = invoke(
        s,
        "num_grid_faces_with",
        &[new_faces.clone(), old_boundary.clone(), Value::int(2)],
    );
    let config = model.config(&inactive).unwrap();
    let solver = invoke(s, "num_grid_acoustics_new", std::slice::from_ref(&config));

    for value in [
        config,
        inactive,
        new_faces,
        old_faces,
        new_boundary,
        old_boundary,
        new_datum,
        old_datum,
        program_a,
        program_b,
    ] {
        assert_eq!(release(s, &value).unwrap(), Value::Bool(true));
    }
    invoke(s, "num_grid_acoustics_step", std::slice::from_ref(&solver));
    assert_eq!(first_pressure(s, &solver), 1.);
}

#[test]
fn detached_sources_use_independent_budgets_and_remain_additive_in_one_cell() {
    let model = Fixture::new();
    let s = &model.session;
    let config = model.config(&model.uniform_faces).unwrap();
    let solver = invoke(s, "num_grid_acoustics_new", std::slice::from_ref(&config));
    let rate = DimExpr::Mul(vec![pressure(), time().pow(-1)]);
    let program = compile(s, "1", rate);
    let invalid = model.source(&solver, &program, 0).unwrap_err();
    assert!(
        invalid.to_string().contains("evaluation limit"),
        "{invalid}"
    );
    assert_eq!(
        model.source(&solver, &program, 1).unwrap(),
        Value::Bool(true)
    );
    assert_eq!(
        model.source(&solver, &program, 1).unwrap(),
        Value::Bool(true)
    );
    assert_eq!(release(s, &program).unwrap(), Value::Bool(true));
    assert!(
        model.source(&solver, &program, 1).is_err(),
        "released aliases cannot publish a new source"
    );
    invoke(s, "num_grid_acoustics_step", std::slice::from_ref(&solver));
    assert_eq!(first_pressure(s, &solver), 1.);
}

#[test]
fn a_scoped_program_can_be_dropped_while_its_earlier_solver_retains_a_source() {
    let model = Fixture::new();
    let s = &model.session;
    let config = model.config(&model.uniform_faces).unwrap();
    let solver = invoke(s, "num_grid_acoustics_new", std::slice::from_ref(&config));
    s.with_native_object_scope(|| {
        let rate = DimExpr::Mul(vec![pressure(), time().pow(-1)]);
        let program = compile(s, "1", rate);
        assert_eq!(
            model.source(&solver, &program, 1).unwrap(),
            Value::Bool(true)
        );
    })
    .unwrap();
    invoke(s, "num_grid_acoustics_step", std::slice::from_ref(&solver));
    assert_eq!(first_pressure(s, &solver), 0.5);
}

#[test]
fn exhaustion_during_step_or_run_invalidates_the_old_solver_without_exposing_partial_state() {
    let model = Fixture::new();
    let s = &model.session;
    let config = model.config(&model.uniform_faces).unwrap();
    for in_one_call in [false, true] {
        let solver = invoke(s, "num_grid_acoustics_new", std::slice::from_ref(&config));
        let rate = DimExpr::Mul(vec![pressure(), time().pow(-1)]);
        let program = compile(s, "1", rate);
        assert_eq!(
            model.source(&solver, &program, 1).unwrap(),
            Value::Bool(true)
        );
        assert_eq!(release(s, &program).unwrap(), Value::Bool(true));
        let source_error = if in_one_call {
            call(
                s,
                "num_grid_acoustics_run",
                &[solver.clone(), Value::int(2)],
            )
        } else {
            assert_eq!(
                invoke(s, "num_grid_acoustics_step", std::slice::from_ref(&solver)),
                Value::Bool(true)
            );
            assert_eq!(first_pressure(s, &solver), 0.5);
            call(s, "num_grid_acoustics_step", std::slice::from_ref(&solver))
        }
        .unwrap_err();
        assert!(
            matches!(source_error, OpError::Domain(ref message) if message == "compiled function evaluation budget exhausted"),
            "{source_error}"
        );
        let state_error = call(
            s,
            "num_grid_acoustics_fields",
            std::slice::from_ref(&solver),
        )
        .unwrap_err();
        assert!(
            matches!(state_error, OpError::Domain(ref message) if message == "acoustic solver invalidated by a failed mutation"),
            "{state_error}"
        );
        assert_eq!(release(s, &solver).unwrap(), Value::Bool(true));
    }

    let replacement = invoke(s, "num_grid_acoustics_new", std::slice::from_ref(&config));
    invoke(
        s,
        "num_grid_acoustics_step",
        std::slice::from_ref(&replacement),
    );
    assert_eq!(first_pressure(s, &replacement), 0.);
}

#[test]
fn detached_point_callbacks_preserve_physics_across_cpu_backends() {
    for selector in ["cpu", "cpu-parallel", "tiered", "cranelift"] {
        let backend = match rssn::backend::named(selector) {
            | Ok(backend) => backend,
            | Err(_) if matches!(selector, "tiered" | "cranelift") => continue,
            | Err(error) => panic!("{selector}: {error}"),
        };
        let model = Fixture::new();
        let s = &model.session;
        let (boundary_program, source_program) = rssn::backend::with_backend(backend, || {
            let rate = DimExpr::Mul(vec![pressure(), time().pow(-1)]);
            (compile(s, "2", pressure()), compile(s, "1", rate))
        });
        let (faces, boundary, datum) =
            model.faces_with_datum(&model.uniform_faces, &boundary_program, 0);
        let config = model.config(&faces).unwrap();
        let solver = invoke(s, "num_grid_acoustics_new", std::slice::from_ref(&config));
        assert_eq!(
            model.source(&solver, &source_program, 1).unwrap(),
            Value::Bool(true)
        );
        for value in [
            config,
            faces,
            boundary,
            datum,
            boundary_program,
            source_program,
        ] {
            assert_eq!(release(s, &value).unwrap(), Value::Bool(true));
        }
        invoke(s, "num_grid_acoustics_step", std::slice::from_ref(&solver));
        assert_eq!(first_pressure(s, &solver), 1.5, "{selector}");
    }
}

#[test]
fn the_native_object_bank_still_checks_sixty_four_mib_and_recovers_after_release() {
    let model = Fixture::new();
    let s = &model.session;
    let config = model.config(&model.uniform_faces).unwrap();
    let mut solvers = (0..7)
        .map(|_| invoke(s, "num_grid_acoustics_new", std::slice::from_ref(&config)))
        .collect::<Vec<_>>();
    let error = call(s, "num_grid_acoustics_new", std::slice::from_ref(&config)).unwrap_err();
    assert!(
        error
            .to_string()
            .contains("native object memory budget exhausted"),
        "{error}"
    );
    assert_eq!(
        release(s, &solvers.pop().unwrap()).unwrap(),
        Value::Bool(true)
    );
    invoke(s, "num_grid_acoustics_new", std::slice::from_ref(&config));
}

mod byte_caps {
    use super::*;
    use rssn::backend::Backend;
    use rssn::backend::BackendError;
    use rssn::backend::Compiled;
    use rssn::backend::NativeProgramBudget;
    use rssn::backend::with_backend;
    use rssn::graph::Graph;
    use rssn::graph::id::NodeId;
    use rssn::graph::id::SymbolId;
    use std::sync::Arc;

    /// A conservative test backend bound, through the same native compile path.
    struct Constant(usize);
    impl Backend for Constant {
        fn uses_gpu(&self) -> bool {
            false
        }

        fn native_program_budget(&self) -> Option<NativeProgramBudget> {
            Some(NativeProgramBudget {
                instructions: 4096,
                extra_bytes: self.0,
            })
        }

        fn retained_bytes(&self) -> Option<usize> {
            Some(0)
        }

        fn name(&self) -> &'static str {
            "interpreter"
        }

        fn compile(
            &self,
            _: &Graph,
            _: NodeId,
            _: &[SymbolId],
        ) -> Result<Box<dyn Compiled>, BackendError> {
            Ok(Box::new(Self(self.0)))
        }
    }
    impl Compiled for Constant {
        fn retained_bytes(&self) -> Option<usize> {
            Some(self.0)
        }

        fn arity(&self) -> usize {
            1
        }

        fn call(
            &self,
            _: &[f64],
        ) -> f64 {
            1.
        }
    }

    #[test]
    fn over_limit_callbacks_do_not_commit_reservations_or_poison_retry() {
        let model = Fixture::new();
        let s = &model.session;
        let hard_boundary = with_backend(Arc::new(Constant(512 * 1024)), || {
            compile(s, "1", pressure())
        });
        let (faces, boundary, datum) =
            model.faces_with_datum(&model.uniform_faces, &hard_boundary, 0);
        let error = model.config(&faces).unwrap_err();
        assert!(
            error
                .to_string()
                .contains("acoustic boundary compiled storage cap exceeded"),
            "{error}"
        );
        for handle in [faces, boundary, datum, hard_boundary] {
            assert_eq!(release(s, &handle).unwrap(), Value::Bool(true));
        }

        let config = model.config(&model.uniform_faces).unwrap();
        let solver = invoke(s, "num_grid_acoustics_new", std::slice::from_ref(&config));
        let rate = DimExpr::Mul(vec![pressure(), time().pow(-1)]);
        let hard_source = with_backend(Arc::new(Constant(2 * 1024 * 1024)), || {
            compile(s, "1", rate.clone())
        });
        let error = model.source(&solver, &hard_source, 1).unwrap_err();
        assert!(
            matches!(error, OpError::Domain(ref message) if message == "acoustic source storage cap exceeded"),
            "{error}"
        );
        assert_eq!(release(s, &hard_source).unwrap(), Value::Bool(true));
        let small = compile(s, "1", rate);
        assert_eq!(model.source(&solver, &small, 1).unwrap(), Value::Bool(true));
        invoke(s, "num_grid_acoustics_step", std::slice::from_ref(&solver));
        assert_eq!(first_pressure(s, &solver), 0.5);
    }

    #[test]
    fn repeated_source_keeps_one_code_charge_and_rejects_a_distinct_program_past_cap() {
        let model = Fixture::new();
        let s = &model.session;
        let config = model.config(&model.uniform_faces).unwrap();
        let solver = invoke(s, "num_grid_acoustics_new", std::slice::from_ref(&config));
        let rate = DimExpr::Mul(vec![pressure(), time().pow(-1)]);
        let shared = with_backend(Arc::new(Constant(1_500_000)), || {
            compile(s, "1", rate.clone())
        });
        let other = with_backend(Arc::new(Constant(800_000)), || compile(s, "1", rate));
        assert_eq!(
            model.source(&solver, &shared, 1).unwrap(),
            Value::Bool(true)
        );
        assert_eq!(
            model.source(&solver, &shared, 1).unwrap(),
            Value::Bool(true)
        );
        let error = model.source(&solver, &other, 1).unwrap_err();
        assert!(
            error
                .to_string()
                .contains("acoustic source storage cap exceeded"),
            "{error}"
        );
        assert_eq!(release(s, &other).unwrap(), Value::Bool(true));
        assert_eq!(
            model.source(&solver, &shared, 1).unwrap(),
            Value::Bool(true)
        );
        assert_eq!(release(s, &shared).unwrap(), Value::Bool(true));
        invoke(s, "num_grid_acoustics_step", std::slice::from_ref(&solver));
        assert_eq!(first_pressure(s, &solver), 1.5);
    }
}
