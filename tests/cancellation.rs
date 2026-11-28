//! Host cancellation interrupts work without reporting partial numeric success.

use std::sync::Arc;
use std::sync::atomic::AtomicBool;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;

use rssn::api::ComputeError;
use rssn::api::Config;
use rssn::api::Session;
use rssn::cancellation;
use rssn::graph::Arity;
use rssn::graph::Budget;
use rssn::graph::Cx;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Graph;
use rssn::graph::Kernel;
use rssn::graph::NodeId;
use rssn::graph::OpDescriptor;
use rssn::graph::OpFlags;
use rssn::graph::OpId;
use rssn::graph::Outcome;
use rssn::graph::Prepared;
use rssn::graph::PreparedWork;
use rssn::graph::RuleSet;
use rssn::graph::Saturate;
use rssn::graph::Stop;
use rssn::graph::Tier;
use rssn::kernels::ode::OdeSolverMethod;
use rssn::kernels::ode::{
    self,
};
use rssn::kernels::ode_adaptive::OdeError;
use rssn::kernels::ode_adaptive::OdeOptions;
use rssn::kernels::ode_adaptive::{
    self,
};

#[test]
fn fixed_ode_stops_when_the_derivative_requests_cancellation() {
    let token = Arc::new(AtomicBool::new(false));
    let calls = AtomicUsize::new(0);
    {
        let _scope = cancellation::install(Arc::clone(&token));
        let result = ode::solve_fixed(
            |_, _, rate| {
                if calls.fetch_add(1, Ordering::Relaxed) == 2 {
                    token.store(true, Ordering::Relaxed);
                }
                rate[0] = 1.0;
            },
            &[0.0],
            (0.0, 10.0),
            10_000,
            OdeSolverMethod::Euler,
        );
        assert_eq!(result.unwrap_err(), "computation cancelled");
        assert_eq!(calls.load(Ordering::Relaxed), 3);
    }
    let next = ode::solve_fixed(
        |_, _, rate| rate[0] = 1.0,
        &[0.0],
        (0.0, 1.0),
        2,
        OdeSolverMethod::Euler,
    )
    .unwrap();
    assert_eq!(next.last(), &[1.0]);
}

#[test]
fn adaptive_ode_returns_typed_cancellation_after_work_begins() {
    let token = Arc::new(AtomicBool::new(false));
    let calls = AtomicUsize::new(0);
    let _scope = cancellation::install(Arc::clone(&token));
    let result = ode_adaptive::dopri5(
        |_, state, rate| {
            if calls.fetch_add(1, Ordering::Relaxed) == 2 {
                token.store(true, Ordering::Relaxed);
            }
            rate[0] = state[0];
        },
        0.0,
        100.0,
        &[1.0],
        &OdeOptions::default(),
        &[],
    );
    assert_eq!(result.unwrap_err(), OdeError::Cancelled);
    assert!(calls.load(Ordering::Relaxed) >= 3);
    assert!(calls.load(Ordering::Relaxed) < 10);
}

struct CancelFirst {
    op: OpId,
    token: Arc<AtomicBool>,
    calls: Arc<AtomicUsize>,
}
impl Kernel for CancelFirst {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        _: NodeId,
    ) -> Outcome {
        if self.calls.fetch_add(1, Ordering::Relaxed) == 0 {
            self.token.store(true, Ordering::Relaxed);
            Outcome::Pass
        } else {
            Outcome::Equal(cx.graph.int(7))
        }
    }
}

#[test]
fn cancelled_graph_answers_are_errors_and_do_not_poison_the_cache() {
    let token = Arc::new(AtomicBool::new(false));
    let calls = Arc::new(AtomicUsize::new(0));
    let rules = {
        let token = Arc::clone(&token);
        let calls = Arc::clone(&calls);
        RuleSet::new("cancel-first", move |installer| {
            let op = installer
                .op(OpDescriptor::new("interruptible", Arity::Fixed(1)).flags(OpFlags::HEAVY))?;
            installer.kernel(
                "cancel-first",
                Tier::Reduce,
                CancelFirst {
                    op,
                    token: Arc::clone(&token),
                    calls: Arc::clone(&calls),
                },
            );
            Ok(())
        })
    };
    let session = Session::with_rules(std::slice::from_ref(&rules)).unwrap();
    let request = session.parse("interruptible(x)").unwrap();
    let config = Config::empty().with(rules);
    let _scope = cancellation::install(Arc::clone(&token));
    assert!(matches!(
        session.compute(request, &config),
        Err(ComputeError::Cancelled)
    ));
    token.store(false, Ordering::Relaxed);
    let answer = session.compute(request, &config).unwrap();
    assert_eq!(answer.term.to_string(), "7");
    assert_eq!(calls.load(Ordering::Relaxed), 2);
}

struct ParallelCancel {
    op: OpId,
    token: Arc<AtomicBool>,
    inherited: Arc<AtomicBool>,
    commits: Arc<AtomicUsize>,
}
struct CancelWork {
    token: Arc<AtomicBool>,
    inherited: Arc<AtomicBool>,
}
impl PreparedWork for CancelWork {
    fn run(self: Box<Self>) -> Box<dyn std::any::Any + Send> {
        if cancellation::current().is_some_and(|current| Arc::ptr_eq(&current, &self.token)) {
            self.inherited.store(true, Ordering::Relaxed);
        }
        self.token.store(true, Ordering::Relaxed);
        Box::new(())
    }
}
impl Kernel for ParallelCancel {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        _: &mut Cx<'_>,
        _: NodeId,
    ) -> Outcome {
        Outcome::Pass
    }

    fn two_phase(&self) -> bool {
        true
    }

    fn prepare(
        &self,
        _: &mut Cx<'_>,
        _: NodeId,
    ) -> Prepared {
        Prepared::Work(Box::new(CancelWork {
            token: Arc::clone(&self.token),
            inherited: Arc::clone(&self.inherited),
        }))
    }

    fn commit(
        &self,
        cx: &mut Cx<'_>,
        _: NodeId,
        _: Box<dyn std::any::Any + Send>,
    ) -> Outcome {
        self.commits.fetch_add(1, Ordering::Relaxed);
        Outcome::Equal(cx.graph.int(1))
    }
}

#[test]
fn rayon_work_inherits_the_scope_and_cancelled_work_is_not_committed() {
    let token = Arc::new(AtomicBool::new(false));
    let inherited = Arc::new(AtomicBool::new(false));
    let commits = Arc::new(AtomicUsize::new(0));
    let rules = {
        let (token, inherited, commits) = (
            Arc::clone(&token),
            Arc::clone(&inherited),
            Arc::clone(&commits),
        );
        RuleSet::new("parallel-cancel", move |installer| {
            let op = installer
                .op(OpDescriptor::new("parallel_work", Arity::Fixed(1)).flags(OpFlags::HEAVY))?;
            installer.kernel(
                "parallel-cancel",
                Tier::Reduce,
                ParallelCancel {
                    op,
                    token: Arc::clone(&token),
                    inherited: Arc::clone(&inherited),
                    commits: Arc::clone(&commits),
                },
            );
            Ok(())
        })
    };
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[rules]).unwrap();
    let root = graph.parse("parallel_work(x) + parallel_work(y)").unwrap();
    let _scope = cancellation::install(Arc::clone(&token));
    let report = engine.run(
        &mut graph,
        &[root],
        &Env::symbolic(),
        &Saturate,
        &Budget::default(),
    );
    assert_eq!(report.stop, Stop::Cancelled);
    assert!(!report.solved);
    assert!(inherited.load(Ordering::Relaxed));
    assert_eq!(commits.load(Ordering::Relaxed), 0);
}

#[test]
fn nested_and_worker_scopes_restore_their_previous_tokens() {
    let outer = Arc::new(AtomicBool::new(false));
    {
        let _outer = cancellation::install(Arc::clone(&outer));
        {
            let _inner = cancellation::install(Arc::new(AtomicBool::new(true)));
            assert!(cancellation::checkpoint().is_err());
        }
        assert!(cancellation::checkpoint().is_ok());
        let inherited = cancellation::current();
        let token = Arc::clone(&outer);
        std::thread::spawn(move || {
            assert!(!cancellation::is_cancelled());
            {
                let _scope = inherited.map(cancellation::install);
                token.store(true, Ordering::Relaxed);
                assert!(cancellation::checkpoint().is_err());
            }
            assert!(!cancellation::is_cancelled());
        })
        .join()
        .unwrap();
        assert!(cancellation::is_cancelled());
        {
            let _unrelated = cancellation::install_optional(None);
            assert!(!cancellation::is_cancelled());
        }
        assert!(cancellation::is_cancelled());
    }
    assert!(!cancellation::is_cancelled());
}
