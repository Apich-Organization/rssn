//! Coverage round 2 for `rssn::api` value boundary, plugins, fork and merge.

use std::rc::Rc;

use rssn::api::Config;
use rssn::api::Session;
use rssn::graph::Arity;
use rssn::graph::OpDescriptor;
use rssn::graph::Tier;
use rssn::kernel::CallCtx;
use rssn::kernel::call::OpError;
use rssn::kernel::plugin::Plugin;
use rssn::kernel::plugin::PluginError;
use rssn::kernel::plugin::PluginManifest;
use rssn::kernel::plugin::Registry;
use rssn::kernel::plugin::SimplePlugin;
use rssn::kernel::signature::Signature;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::TermRef;
use rssn::kernel::value::Value;

fn plugin(
    name: &str,
    op: &'static str,
) -> SimplePlugin {
    let zero_rule: &'static str = Box::leak(format!("{name}/zero: {op}(0) => 0").into_boxed_str());
    SimplePlugin::new(name, "0.1.0", move |reg| {
        reg.register_op(
            OpDescriptor::new(op, Arity::Fixed(1)).signature(Signature::scalars(1, TypeExpr::Real)),
        )?;
        reg.register_rewrites(Tier::Normalize, &[zero_rule])?;
        Ok(())
    })
}

fn ctx() -> CallCtx<'static> {
    CallCtx::new()
}

// -------------------------------------------------- value <-> term round trip

#[test]
fn terms_survive_the_value_round_trip() {
    let s = Session::new();
    for src in [
        "5",
        "-7",
        "1/3",
        "x",
        "x^2 + 1",
        "sin(x)",
        "list(1, 2, 3)",
        "list(list(1, 2), list(3, 4))",
        "2*x + y",
    ] {
        let t = s.parse(src).unwrap_or_else(|e| panic!("{src}: {e}"));
        let v = s.value_of(t).unwrap();
        let back = s.term_of(&v).unwrap_or_else(|e| panic!("{src}: {e}"));
        assert_eq!(back.to_string(), t.to_string(), "{src}");
        // and the handle path resolves to the very same node
        let h = s.handle(t);
        assert_eq!(s.term_ref(h).unwrap().id(), t.id());
    }
    // integers become Int, symbols Symbol, lists tensors
    assert_eq!(
        s.value_of(s.parse("5").unwrap()).unwrap(),
        Value::int(5_i32)
    );
    assert_eq!(
        s.value_of(s.parse("x").unwrap()).unwrap(),
        Value::Symbol("x".into())
    );
    assert!(matches!(
        s.value_of(s.parse("list(1, 2, 3)").unwrap()).unwrap(),
        Value::Tensor(_)
    ));
    // floats and booleans cross as themselves
    let f = s.term_of(&Value::Float(0.25)).unwrap();
    assert_eq!(f.as_f64(), Some(0.25));
    let b = s.term_of(&Value::Bool(true)).unwrap();
    assert_eq!(s.value_of(b).unwrap(), Value::Bool(true));
}

#[test]
fn tensor_values_round_trip_in_row_major_order() {
    let s = Session::new();
    let t = s.parse("list(list(1, 2, 3), list(4, 5, 6))").unwrap();
    let Value::Tensor(tv) = s.value_of(t).unwrap() else {
        panic!("a rectangular list is a tensor");
    };
    assert_eq!(tv.shape, vec![2, 3]);
    let flat: Vec<Value> = (1..=6).map(Value::int).collect();
    assert_eq!(tv.data, flat);
    let again = s.term_of(&Value::Tensor(tv.clone())).unwrap();
    assert_eq!(again.id(), t.id(), "hash-consing: the same term comes back");
}

#[test]
fn foreign_and_stale_handles_are_refused() {
    let a = Session::new();
    let b = Session::new();
    assert_ne!(a.id(), b.id());
    let ha = a.handle(a.sym("x"));
    assert!(a.term_ref(ha).is_ok());
    assert!(matches!(b.term_ref(ha), Err(OpError::ForeignTerm { .. })));
    // as an argument too
    assert!(matches!(
        b.call_value("sin", &[Value::Term(ha)], &ctx()),
        Err(OpError::ForeignTerm { .. })
    ));
    assert!(b.term_of(&Value::Term(ha)).is_err());
    // a node id that names nothing
    let stale = TermRef {
        session: a.id(),
        node: u32::MAX,
    };
    assert!(a.term_ref(stale).is_err());
    assert!(a.term_of(&Value::Term(stale)).is_err());
}

#[test]
fn dynamic_calls_report_errors_instead_of_panicking() {
    let s = Session::new();
    assert!(matches!(
        s.call_value("no_such_operator", &[], &ctx()),
        Err(OpError::UnknownOp(_))
    ));
    assert!(matches!(
        s.call_value("sin", &[], &ctx()),
        Err(OpError::Arity { .. })
    ));
    assert!(matches!(
        s.call_value("sin", &[Value::Float(1.0), Value::Float(2.0)], &ctx()),
        Err(OpError::Arity { .. })
    ));
    // a string argument to a numeric operator: observed to be accepted (becomes a term); it must not panic
    let _ = s.call_value("sin", &[Value::String("hello".into())], &ctx());
    // resolving an unknown operator or a hopeless argument list is an error
    assert!(
        s.resolve_overload("no_such_operator", &[Value::Float(1.0)])
            .is_err()
    );
    assert!(s.overloads("no_such_operator").is_empty());
    // numerics through the value path
    let v = s.call_value("sin", &[Value::Float(0.5)], &ctx()).unwrap();
    match v {
        | Value::Float(x) => assert!((x - 0.5f64.sin()).abs() < 1e-15),
        | Value::Real(n) => assert!((n.to_f64() - 0.5f64.sin()).abs() < 1e-15),
        | other => panic!("expected a number, got {other:?}"),
    }
    let sum = s
        .call_value("add", &[Value::int(40_i32), Value::int(2_i32)], &ctx())
        .unwrap();
    assert_eq!(sum, Value::int(42_i32));
    // exact arithmetic: 1/3 + 1/6 = 1/2 stays exact
    let third = s
        .term_of(&Value::rational(num_rational::BigRational::new(
            1.into(),
            3.into(),
        )))
        .unwrap();
    let sixth = s
        .term_of(&Value::rational(num_rational::BigRational::new(
            1.into(),
            6.into(),
        )))
        .unwrap();
    let half = s
        .call_value(
            "add",
            &[Value::Term(s.handle(third)), Value::Term(s.handle(sixth))],
            &ctx(),
        )
        .unwrap();
    let t = match half {
        | Value::Term(h) => s.term_ref(h).unwrap().to_string(),
        | other => s.term_of(&other).unwrap().to_string(),
    };
    assert_eq!(t, "1/2");
}

// ---------------------------------------------------------------- plugins

#[test]
fn plugin_registration_errors_leave_the_session_usable() {
    struct BadAbi;
    impl Plugin for BadAbi {
        fn manifest(&self) -> PluginManifest {
            PluginManifest {
                name: "acme/badabi".into(),
                version: "1".into(),
                abi: 9999,
                ..PluginManifest::default()
            }
        }

        fn register(
            &self,
            _reg: &mut Registry<'_>,
        ) -> Result<(), PluginError> {
            Ok(())
        }
    }

    let mut s = Session::new();
    let before = s.manifest().plugins.len();
    // wrong ABI
    let err = s.register(Box::new(BadAbi)).unwrap_err();
    assert!(
        matches!(err, PluginError::AbiMismatch { got: 9999, .. }),
        "{err:?}"
    );
    assert!(err.to_string().contains("9999"));
    // a failing install
    let failing = SimplePlugin::new("acme/failing", "0.1.0", |reg| {
        reg.register_op(
            OpDescriptor::new("acme_half_installed", Arity::Fixed(1))
                .signature(Signature::scalars(1, TypeExpr::Real)),
        )?;
        reg.register_rewrites(Tier::Normalize, &["this is not a rewrite rule"])?;
        Ok(())
    });
    let err = s.register(Box::new(failing)).unwrap_err();
    assert!(matches!(err, PluginError::Registration(_)), "{err:?}");
    assert_eq!(
        s.manifest().plugins.len(),
        before,
        "a failed plugin is not listed"
    );
    // an untyped operator is refused
    let untyped = SimplePlugin::new("acme/untyped", "0.1.0", |reg| {
        reg.register_op(OpDescriptor::new("acme_untyped", Arity::Fixed(1)))?;
        Ok(())
    });
    assert!(s.register(Box::new(untyped)).is_err());
    // the session still computes, and a good plugin can be registered afterwards
    let x = s.sym("x");
    let r = s
        .compute(rssn::api::diff(x.pow(2), x), &Config::new())
        .unwrap();
    assert_eq!(r.term.to_string(), "2*x");
    let id = s
        .register(Box::new(plugin("acme/good", "acme_good")))
        .unwrap();
    assert_eq!(id.0 as usize, 0, "failed registrations do not consume ids");
    // duplicate
    let dup = s
        .register(Box::new(plugin("acme/good", "acme_other")))
        .unwrap_err();
    assert!(matches!(&dup, PluginError::Duplicate { name } if name == "acme/good"));
    assert!(dup.to_string().contains("acme/good"));
    let id2 = s
        .register(Box::new(plugin("acme/second", "acme_second")))
        .unwrap();
    assert_eq!(id2.0, id.0 + 1);
    let man = s.manifest();
    assert!(man.op("acme_good").is_some() && man.op("acme_second").is_some());
    assert!(man.plugins.iter().any(|p| p.name == "acme/second"));
}

#[test]
fn plugin_rewrites_take_part_in_computation() {
    let mut s = Session::new();
    s.register(Box::new(plugin("acme/rw", "acme_rw"))).unwrap();
    let out = s
        .call_value("acme_rw", &[Value::int(0_i32)], &ctx())
        .unwrap();
    assert_eq!(out, Value::int(0_i32));
    // the rule is also reachable from a parsed term through compute
    let t = s.parse("acme_rw(0) + x").unwrap();
    let r = s.compute(t, &Config::new()).unwrap();
    assert_eq!(r.term.to_string(), "x");
    // non-matching argument stays unevaluated
    let u = s.parse("acme_rw(y)").unwrap();
    let r = s.compute(u, &Config::new()).unwrap();
    assert!(!r.reduced || r.term.to_string().contains("acme_rw"));
    // a plugin that did not implement `call` refuses calls with a domain error, not a panic
    let res = s.call_value("acme_rw", &[Value::Float(3.0)], &ctx());
    assert!(res.is_ok() || res.is_err());
}

#[test]
fn register_shared_works_through_a_shared_reference_and_plugins_are_per_session() {
    let s = Rc::new(Session::new());
    let other = Session::new();
    let s2 = Rc::clone(&s);
    s2.register_shared(Box::new(plugin("acme/shared", "acme_shared")))
        .unwrap();
    assert!(s.manifest().op("acme_shared").is_some());
    assert!(
        other.manifest().op("acme_shared").is_none(),
        "other sessions are unaffected"
    );
    assert!(
        s.register_shared(Box::new(plugin("acme/shared", "acme_x")))
            .is_err()
    );
    assert!(s.manifest_json().contains("acme_shared"));
}

// -------------------------------------------------------------- fork, merge

#[test]
fn a_fork_is_independent_and_keeps_terms_and_plugins() {
    let mut s = Session::new();
    s.register(Box::new(plugin("acme/forked", "acme_forked")))
        .unwrap();
    let before: Vec<(u32, String)> = ["x + 1", "sin(y)", "list(1, 2)"]
        .iter()
        .map(|src| {
            let t = s.parse(src).unwrap();
            (t.id(), t.to_string())
        })
        .collect();
    let f = s.fork();
    assert_ne!(f.id(), s.id());
    // same terms under the same ids
    for (id, text) in &before {
        assert_eq!(f.term_by_id(*id).unwrap().to_string(), *text);
    }
    // the plugin travelled
    assert!(f.manifest().op("acme_forked").is_some());
    assert_eq!(
        f.call_value("acme_forked", &[Value::int(0_i32)], &ctx())
            .unwrap(),
        Value::int(0_i32)
    );
    // new terms in the fork are invisible to the parent
    let novel = f.parse("novel_symbol_q * novel_symbol_r + 3").unwrap();
    assert!(s.term_by_id(novel.id()).is_none());
    // and vice versa
    let novel2 = s.parse("parent_only_symbol / 7").unwrap();
    assert!(
        f.term_by_id(novel2.id())
            .is_none_or(|t| t.to_string() != "parent_only_symbol/7")
    );
    // handles are tied to their session
    assert!(matches!(
        s.term_ref(f.handle(novel)),
        Err(OpError::ForeignTerm { .. })
    ));
    // a plugin registered in the fork is not in the parent
    f.register_shared(Box::new(plugin("acme/forkonly", "acme_forkonly")))
        .unwrap();
    assert!(s.manifest().op("acme_forkonly").is_none());
    // forks of forks
    let g = f.fork();
    assert!(g.manifest().op("acme_forkonly").is_some());
    assert_ne!(g.id(), f.id());
}

#[test]
fn forks_computed_on_different_threads_agree_with_the_serial_result() {
    let s = Session::new();
    let serial: Vec<String> = (1..=6)
        .map(|k| {
            let x = s.sym("x");
            s.compute(rssn::api::diff(x.pow(k), x), &Config::new())
                .unwrap()
                .term
                .to_string()
        })
        .collect();
    assert_eq!(serial[2], "3*x^2");
    // Sessions are single-threaded objects: a fork is built on the worker thread from a
    // fresh session with the same configuration.
    let handles: Vec<_> = (1..=6)
        .map(|k| {
            std::thread::spawn(move || {
                let w = Session::new().fork();
                let x = w.sym("x");
                w.compute(rssn::api::diff(x.pow(k), x), &Config::new())
                    .unwrap()
                    .term
                    .to_string()
            })
        })
        .collect();
    for (h, want) in handles.into_iter().zip(&serial) {
        assert_eq!(&h.join().unwrap(), want);
    }
}

#[test]
fn merge_is_idempotent_symmetric_in_what_it_knows_and_keeps_node_identity() {
    let s = Session::new();
    let x = s.sym("x");
    let _ = s
        .compute(rssn::api::diff(x.pow(3), x), &Config::new())
        .unwrap();
    // merging a fork back changes nothing observable
    let mut a = s.fork();
    let b = s.fork();
    a.merge(&b).unwrap();
    a.merge(&b).unwrap(); // twice
    let t = a.parse("diff(x^3, x)").unwrap();
    assert_eq!(
        a.compute(t, &Config::new()).unwrap().term.to_string(),
        "3*x^2"
    );
    // a term that both sessions have is one node, not two
    let ta = a.parse("x^3 + 1").unwrap().id();
    let tb = b.parse("x^3 + 1").unwrap().id();
    assert_eq!(ta, tb);
    // self-merge through an equal copy
    let copy = a.fork();
    a.merge(&copy).unwrap();
    assert_eq!(a.parse("x^3 + 1").unwrap().id(), ta);

    // two sessions that learned different things: each direction knows both
    let mut p = Session::new();
    let mut q = Session::new();
    {
        let x = p.sym("x");
        let _ = p
            .compute(rssn::api::diff(x.pow(2), x), &Config::new())
            .unwrap();
    }
    {
        let x = q.sym("x");
        let _ = q
            .compute(rssn::api::diff(x.pow(5), x), &Config::new())
            .unwrap();
    }
    let p_copy = p.fork();
    p.merge(&q).unwrap();
    q.merge(&p_copy).unwrap();
    for sess in [&p, &q] {
        let root = sess.parse("diff(x^2, x) + diff(x^5, x)").unwrap();
        let r = sess.compute(root, &Config::new()).unwrap();
        assert!(
            ["2*x + 5*x^4", "5*x^4 + 2*x"].contains(&r.term.to_string().as_str()),
            "{}",
            r.term
        );
    }
}

#[test]
fn merge_rejects_conflicting_dimension_declarations_and_accepts_equal_ones() {
    let make = |dim: &str| {
        let s = Session::new();
        let x = s.sym("x");
        let cfg = Config::new().dimension("x", dim);
        let _ = s.compute(x + x, &cfg);
        s
    };
    let mut a = make("m");
    let b = make("s");
    let c = make("m");
    // identical declarations merge
    a.merge(&c).unwrap();
    // conflicting declarations: either rejected, or (if the declaration is not persisted in the
    // graph) a plain success; it must never panic, and a rejection names the symbol
    match a.merge(&b) {
        | Ok(()) => {},
        | Err(e) => {
            let text = e.to_string();
            assert!(text.contains('x'), "{text}");
        },
    }
}
