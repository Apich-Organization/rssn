//! The SDE noise rule set: the closure-free part of the stochastic solvers,
//! exposed as manifest operators.
//!
//! The integrators in [`crate::kernels::sde`] take closures for the drift and
//! the diffusion, so they cannot be manifest operators: an operator is called
//! with values, and a host callback is not a value. What *can* cross the
//! boundary is the noise itself — a Brownian sample is a pure function of a
//! seed and the shape of the grid — so this rule set exposes exactly those
//! primitives:
//!
//! | operator | value |
//! |---|---|
//! | `sde_wiener(seed, m, steps, h)` | the `steps x m` matrix of Brownian increments `dW[k][j] ~ N(0, h)` |
//! | `sde_levy_areas(seed, h, dw, terms)` | the `m*m` matrix (row-major) of iterated Ito integrals `I_jl` over one step, conditioned on the increments `dw` |
//! | `sde_scheme_info(name)` | `list(strong_order, ito, commutative_only, adaptive, levy_areas)` for one scheme name |
//!
//! The first two are pure functions of their arguments: the same seed and shape
//! give the same numbers bit for bit, so a host can reproduce a path after the
//! fact, and two hosts can compare paths without exchanging the increments. A
//! host that has its own increments can ask for the Lévy areas alone with
//! `sde_levy_areas`; the areas are *conditional* on those increments, exactly
//! as in [`levy_areas`](crate::kernels::sde::levy_areas), so the increments are
//! an argument rather than a second draw.
//!
//! `sde_scheme_info` answers the question "what does this solver promise?"
//! without running it: the strong order and four design flags
//! (`1`/`0`) for the eight schemes of [`crate::kernels::sde`] — the Ito versus
//! Stratonovich interpretation, whether the scheme needs commutative noise,
//! whether it adapts the step size, and whether it needs the Lévy areas. It is
//! metadata, not a solver, and an unknown name stays a request.
//!
//! # Dimensions
//!
//! `sde_wiener`'s dimension is `sqrt(dim(h))` — the root of argument 3, which
//! is why the rule set needs [`DimRule::RootOf`] rather than
//! [`DimRule::Root`] (that one is hard-wired to argument 0). A Brownian
//! increment over a time step has the dimension of the square root of time.
//! `sde_levy_areas` has the dimension of `h`: an iterated integral
//! `∫∫ dW dW` is `T` when the Wiener process is `sqrt(T)`.
//!
//! # Conventions
//!
//! The kernel encodes the Ito/Stratonovich interpretation in *function names*
//! (`milstein` is Ito, `milstein_stratonovich` is Stratonovich). The noise
//! crossing this rule set is convention-free — a Brownian increment is the
//! same object either way; the Stratonovich schemes use the same increments
//! and drop the `-h` term in the iterated integral instead. See the
//! [`crate::kernels::sde`] module header for the scheme-by-scheme table.
//!
//! # Limitations
//!
//! Only the noise is exposed. The solver loops themselves stay Rust-only: they
//! take host closures for `a(t, x)` and `g(t, x)`, and wrapping them as
//! operators would need a callback argument, which `Value` does not carry for
//! the ensembles (they share `&Fn + Sync` and never mutate the closure). A
//! host that wants the ensemble statistics computes them around the noise this
//! rule set produces.
//!
//! `m` is not a parameter of `sde_levy_areas`: the length of `dw` *is* the
//! number of Wiener processes, and a second copy could only disagree with it.
//! `terms` bounds the Fourier modes of the bridge (accuracy improves like
//! `terms^{-3/2}`, per `levy_areas`). A request whose arguments are symbolic,
//! non-integral where an integer is required, out of the documented bounds
//! (`1 <= m <= 4096`, `1 <= steps <= 1_000_000`, `steps*m <= 4_194_304`,
//! `1 <= terms <= 1000`), or whose `h` is not a positive finite number stays a
//! request ([`Outcome::Pass`]): the operator never clamps a parameter to make
//! it fit, because a silently resized grid is a different sample.
//!
//! `sde_scheme_info` takes its scheme name as a `Value::String`; rssn's surface
//! syntax has no string literal, so it is a host-facing operator (the name
//! crosses the dynamic ABI as a value) rather than one spelled in parsed text,
//! and a name it does not know stays a request.

use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::Tier;
use crate::graph::dim::DimRule;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::kernels::random::Rng;
use crate::kernels::sde::brownian_increments;
use crate::kernels::sde::levy_areas;

/// Largest number of Wiener processes one request may draw.
const MAX_WIENER: usize = 4096;
/// Largest number of time steps one request may draw.
const MAX_STEPS: usize = 1_000_000;
/// Largest `steps * m` one request may draw (bounds the returned tensor).
const MAX_SAMPLES: usize = 4_194_304;
/// Largest number of bridge modes one area request may use.
const MAX_TERMS: usize = 1000;

/// The SDE-noise rule set.
#[must_use]
pub fn sde() -> RuleSet {
    RuleSet::new("sde", install).about(env!("CARGO_PKG_VERSION"), "3.0.0", &["sde", "stochastic"])
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let wiener = i.op(OpDescriptor::new("sde_wiener", Arity::Fixed(4))
        .flags(OpFlags::HEAVY)
        .cost(100)
        .category("sde")
        .since("3.0.0")
        .doc(
            "`sde_wiener(seed, m, steps, h)`: the `steps x m` matrix of seeded Brownian \
                 increments `dW[k][j] ~ N(0, h)`, reproducible bit for bit from `seed`.",
        ))?;
    let areas = i.op(OpDescriptor::new("sde_levy_areas", Arity::Fixed(4))
        .flags(OpFlags::HEAVY)
        .cost(100)
        .category("sde")
        .since("3.0.0")
        .doc(
            "`sde_levy_areas(seed, h, dw, terms)`: the `m*m` row-major iterated Ito \
                 integrals over one step, conditioned on the `m` increments `dw`, with a \
                 seeded Fourier bridge truncated at `terms` modes.",
        ))?;
    // A Brownian increment over `h` is `sqrt(h)`: the rooted operand is the
    // fourth argument, not the first, so the explicit-operand rule is needed.
    i.dim(wiener, DimRule::RootOf { degree: 2, index: 3 });
    // `I_jl = int int dW_j dW_l` has the dimension of the step length.
    i.dim(areas, DimRule::Same(1));
    i.kernel("sde/wiener", Tier::Reduce, Wiener { op: wiener });
    i.kernel("sde/levy-areas", Tier::Reduce, LevyAreas { op: areas });
    let info = i.op(OpDescriptor::new("sde_scheme_info", Arity::Fixed(1))
        .flags(OpFlags::HEAVY)
        .cost(10)
        .category("sde")
        .since("3.0.0")
        .doc(
            "`sde_scheme_info(name)`: the metadata of one stochastic scheme as \
                 `list(strong_order, ito, commutative_only, adaptive, levy_areas)` — the \
                 strong order and four flags (1 or 0). Unknown names stay a request.",
        ))?;
    // A metadata vector is a label result, not a physical quantity: the string
    // argument is non-metric, and so is the answer, so a dimensional check never
    // tries to unify either with a scalar.
    i.dim(info, DimRule::NonMetric);
    i.kernel("sde/scheme-info", Tier::Reduce, SchemeInfo { op: info });
    Ok(())
}

/// The integer an exact number holds, if it is one.
///
/// Exact rationals that are integral are accepted as well as `Int` (a host can
/// hand the ABI `Number::Rat(3/1)`); a float is not an integer here, however
/// round it looks, because the request must be reproducible from exactly the
/// arguments it was given.
fn integer(
    graph: &Graph,
    node: NodeId,
) -> Option<u64> {
    use num_traits::ToPrimitive;

    let rational = graph.number_of(node)?.to_rational()?;
    if !rational.is_integer() {
        return None;
    }
    rational.to_integer().to_u64()
}

/// The real number an argument holds, exact or floating.
fn real(
    graph: &Graph,
    node: NodeId,
) -> Option<f64> {
    graph.number_of(node).map(Number::to_f64)
}

/// `sde_wiener`: a seeded Brownian sample.
struct Wiener {
    op: OpId,
}

impl Kernel for Wiener {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        let &[seed, m, steps, h] = args.as_slice() else {
            return Outcome::Pass;
        };
        let (Some(seed), Some(m), Some(steps), Some(h)) = (
            integer(cx.graph, seed),
            integer(cx.graph, m),
            integer(cx.graph, steps),
            real(cx.graph, h),
        ) else {
            return Outcome::Pass;
        };
        let (Ok(m), Ok(steps)) = (usize::try_from(m), usize::try_from(steps)) else {
            return Outcome::Pass;
        };
        if !h.is_finite() || h <= 0.0 {
            return Outcome::Pass;
        }
        if !(1..=MAX_WIENER).contains(&m) || !(1..=MAX_STEPS).contains(&steps) {
            return Outcome::Pass;
        }
        if m.checked_mul(steps)
            .is_none_or(|samples| samples > MAX_SAMPLES)
        {
            return Outcome::Pass;
        }
        let mut rng = Rng::new(seed);
        let increments = brownian_increments(&mut rng, m, steps, h);
        let rows: Vec<NodeId> = increments
            .iter()
            .map(|row| {
                let items: Vec<NodeId> = row.iter().map(|&value| cx.graph.float(value)).collect();
                cx.graph.node(core::LIST, &items)
            })
            .collect();
        Outcome::Equal(cx.graph.node(core::LIST, &rows))
    }
}

/// `sde_levy_areas`: the iterated Ito integrals over one step, given its
/// increments.
struct LevyAreas {
    op: OpId,
}

impl Kernel for LevyAreas {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        let &[seed, h, increments, terms] = args.as_slice() else {
            return Outcome::Pass;
        };
        let (Some(seed), Some(h), Some(terms)) = (
            integer(cx.graph, seed),
            real(cx.graph, h),
            integer(cx.graph, terms),
        ) else {
            return Outcome::Pass;
        };
        let Some(dw) = list_of_numbers(cx.graph, increments) else {
            return Outcome::Pass;
        };
        let Ok(terms) = usize::try_from(terms) else {
            return Outcome::Pass;
        };
        let m = dw.len();
        if !h.is_finite() || h <= 0.0 {
            return Outcome::Pass;
        }
        if !(1..=MAX_WIENER).contains(&m) || !(1..=MAX_TERMS).contains(&terms) {
            return Outcome::Pass;
        }
        let mut rng = Rng::new(seed);
        let areas = levy_areas(&mut rng, m, h, &dw, terms);
        let items: Vec<NodeId> = areas.iter().map(|&value| cx.graph.float(value)).collect();
        Outcome::Equal(cx.graph.node(core::LIST, &items))
    }
}

/// The scheme-name string an argument holds.
fn scheme_name(
    graph: &Graph,
    node: NodeId,
) -> Option<&str> {
    match graph.payload(node) {
        | Some(crate::graph::Payload::Str(text)) => Some(text),
        | _ => None,
    }
}

/// The metadata of one scheme, in the order `sde_scheme_info` returns it:
/// `[strong_order, ito, commutative_only, adaptive, levy_areas]`.
///
/// The strong orders are the ones the module header and
/// `tests/kernels/sde.rs` establish (Euler-Maruyama 1/2, everything else
/// strong order 1); the flags are design facts: `ito` is 0 for the two
/// Stratonovich schemes, `commutative_only` is 1 exactly where the header's
/// limitation names the scheme (the Heun scheme is the exception that handles
/// general noise), `adaptive` marks the step-controlled solvers and
/// `levy_areas` marks the non-commutative ones that need the iterated
/// integrals.
fn scheme_fields(name: &str) -> Option<[f64; 5]> {
    let fields = match name {
        | "euler_maruyama" => [0.5, 1.0, 0.0, 0.0, 0.0],
        | "milstein" | "srk_platen" => [1.0, 1.0, 1.0, 0.0, 0.0],
        | "milstein_stratonovich" => [1.0, 0.0, 1.0, 0.0, 0.0],
        | "milstein_noncommutative" => [1.0, 1.0, 0.0, 0.0, 1.0],
        | "heun_stratonovich" => [1.0, 0.0, 0.0, 0.0, 0.0],
        | "sde_adaptive" => [1.0, 1.0, 1.0, 1.0, 0.0],
        | "sde_adaptive_noncommutative" => [1.0, 1.0, 0.0, 1.0, 1.0],
        | _ => return None,
    };
    Some(fields)
}

/// `sde_scheme_info`: the metadata of one stochastic scheme.
struct SchemeInfo {
    op: OpId,
}

impl Kernel for SchemeInfo {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        let &[scheme] = args.as_slice() else {
            return Outcome::Pass;
        };
        let Some(fields) = scheme_name(cx.graph, scheme).and_then(scheme_fields) else {
            // A symbolic argument or an unknown scheme name: the caller asked a
            // question this operator cannot answer, so it stays a request
            // rather than getting a guessed default.
            return Outcome::Pass;
        };
        let items: Vec<NodeId> = fields.iter().map(|&value| cx.graph.float(value)).collect();
        Outcome::Equal(cx.graph.node(core::LIST, &items))
    }
}

/// The numbers a `list` term holds, in order.
///
/// `None` for anything that is not a list of numeric leaves: a symbolic
/// increment would make the areas incomplete, and completing it with zeros is
/// the kind of plausible wrong answer this rule set refuses.
#[allow(clippy::single_call_fn)] // one kernel; separate so the contract is documented
fn list_of_numbers(
    graph: &Graph,
    node: NodeId,
) -> Option<Vec<f64>> {
    if graph.op(node) != core::LIST {
        return None;
    }
    graph
        .children(node)
        .iter()
        .map(|&item| graph.number_of(item).map(Number::to_f64))
        .collect()
}

#[cfg(test)]
mod tests {
    use crate::api::Session;
    use crate::graph::Engine;
    use crate::graph::Graph;
    use crate::graph::RuleSet;
    use crate::kernel::CallCtx;
    use crate::kernel::KernelManifest;
    use crate::kernel::SignatureSource;
    use crate::kernel::infer::signatures;
    use crate::kernel::signature::DimExpr;
    use crate::kernel::signature::DimRuleInfo;
    use crate::kernel::signature::Exp;
    use crate::kernel::signature::NatExpr;
    use crate::kernel::signature::Rational;
    use crate::kernel::signature::TypeExpr;
    use crate::kernel::value::TensorValue;
    use crate::kernel::value::Value;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    use super::sde;

    fn installed() -> (Graph, KernelManifest) {
        let sets: Vec<RuleSet> = vec![sde()];
        let mut graph = Graph::new();
        Engine::install(&mut graph, &sets).unwrap_or_else(|e| panic!("{e}"));
        let manifest = KernelManifest::build(&graph, &sets, &[]);
        (graph, manifest)
    }

    fn real(value: f64) -> Value {
        Value::Float(value)
    }

    fn floats(value: &Value) -> Vec<f64> {
        let Value::Tensor(tensor) = value else {
            panic!("expected a tensor, got {value:?}");
        };
        tensor
            .data
            .iter()
            .map(|entry| {
                let Value::Real(number) = entry else {
                    panic!("expected a real, got {entry:?}");
                };
                number.to_f64()
            })
            .collect()
    }

    fn call(
        name: &str,
        args: &[Value],
    ) -> Value {
        let session = Session::new();
        session
            .call_value(name, args, &CallCtx::new())
            .unwrap_or_else(|e| panic!("{name}: {e}"))
    }

    fn bits(values: &[f64]) -> Vec<u64> {
        values.iter().map(|value| value.to_bits()).collect()
    }

    #[test]
    fn the_noise_operators_are_manifest_operators() {
        let (graph, manifest) = installed();
        assert!(
            manifest.features.sde,
            "the build reports the SDE capability"
        );
        let set = manifest
            .rule_sets
            .iter()
            .find(|set| set.name == "sde")
            .unwrap_or_else(|| panic!("the sde rule set is not described"));
        assert!(set.capabilities.iter().any(|c| c == "sde"));
        for (name, arity) in [
            ("sde_wiener", 4_u8),
            ("sde_levy_areas", 4),
            ("sde_scheme_info", 1),
        ] {
            let op = manifest
                .op(name)
                .unwrap_or_else(|| panic!("{name} is not in the manifest"));
            assert_eq!(op.arity.min, arity, "{name}");
            assert_eq!(op.category, "sde", "{name}");
            assert_eq!(op.since, "3.0.0", "{name}");
            assert!(!op.doc.is_empty(), "{name} has no documentation");
            let (overloads, source) = signatures(
                graph.ops().get(
                    graph
                        .ops()
                        .lookup(name)
                        .unwrap_or_else(|| panic!("{name} is not registered")),
                ),
            );
            assert_eq!(source, SignatureSource::Overlay, "{name}");
            assert_eq!(overloads.len(), 1, "{name}");
        }
        assert_eq!(
            manifest.op("sde_wiener").and_then(|op| op.dim_rule.clone()),
            Some(DimRuleInfo::RootOf { degree: 2, index: 3 })
        );
        assert_eq!(
            manifest
                .op("sde_levy_areas")
                .and_then(|op| op.dim_rule.clone()),
            Some(DimRuleInfo::Same(1))
        );
        assert_eq!(
            manifest
                .op("sde_scheme_info")
                .and_then(|op| op.dim_rule.clone()),
            Some(DimRuleInfo::NonMetric),
            "a metadata vector is a label result"
        );
        assert!(manifest.verify(), "the hash must match the contents");
    }

    #[test]
    fn the_overlay_signature_states_the_root_dimension() {
        let (graph, _) = installed();
        let op = graph
            .ops()
            .lookup("sde_wiener")
            .unwrap_or_else(|| panic!("sde_wiener is registered"));
        let (overloads, source) = signatures(graph.ops().get(op));
        assert_eq!(source, SignatureSource::Overlay);
        let signature = overloads.first().unwrap_or_else(|| panic!("one overload"));
        let half = TypeExpr::Real
            .with_dim(DimExpr::Var(String::from("T")).pow(Exp::Rat(Rational::new(1, 2))));
        let dimensioned_parameter = TypeExpr::Real.with_dim(DimExpr::Var(String::from("T")));
        assert_eq!(signature.params.len(), 4);
        assert_eq!(signature.params.get(3), Some(&dimensioned_parameter));
        let TypeExpr::Tensor { rank, elem, .. } = &signature.ret else {
            panic!("the result is not a tensor: {:?}", signature.ret);
        };
        assert_eq!(*rank, NatExpr::Lit(2));
        assert_eq!(**elem, half);

        let areas = graph
            .ops()
            .lookup("sde_levy_areas")
            .unwrap_or_else(|| panic!("sde_levy_areas is registered"));
        let (overloads, _) = signatures(graph.ops().get(areas));
        let signature = overloads.first().unwrap_or_else(|| panic!("one overload"));
        let TypeExpr::Tensor { rank, elem, .. } = &signature.ret else {
            panic!("the result is not a tensor: {:?}", signature.ret);
        };
        assert_eq!(*rank, NatExpr::Lit(1));
        assert_eq!(**elem, dimensioned_parameter);
    }

    #[test]
    fn wiener_increments_are_reproducible_and_shaped() {
        let first = call(
            "sde_wiener",
            &[
                Value::int(7_i32),
                Value::int(2_i32),
                Value::int(3_i32),
                real(0.01),
            ],
        );
        let second = call(
            "sde_wiener",
            &[
                Value::int(7_i32),
                Value::int(2_i32),
                Value::int(3_i32),
                real(0.01),
            ],
        );
        let Value::Tensor(tensor) = &first else {
            panic!("expected a tensor, got {first:?}");
        };
        assert_eq!(tensor.shape, vec![3, 2], "steps x m");
        assert_eq!(
            bits(&floats(&first)),
            bits(&floats(&second)),
            "the same seed must reproduce the sample bit for bit"
        );
        let other = call(
            "sde_wiener",
            &[
                Value::int(8_i32),
                Value::int(2_i32),
                Value::int(3_i32),
                real(0.01),
            ],
        );
        assert_ne!(
            bits(&floats(&first)),
            bits(&floats(&other)),
            "a different seed must give a different sample"
        );
    }

    #[test]
    fn levy_areas_are_conditional_on_the_increments() {
        let increments = [0.1_f64, -0.2];
        let h = 0.01_f64;
        let dw = Value::Tensor(TensorValue::vector(vec![real(0.1), real(-0.2)]));
        let args = [Value::int(3_i32), real(h), dw.clone(), Value::int(8_i32)];
        let first = call("sde_levy_areas", &args);
        let second = call("sde_levy_areas", &args);
        let areas = floats(&first);
        assert_eq!(areas.len(), 4, "m*m");
        assert_eq!(bits(&areas), bits(&floats(&second)));
        // I_jj = (dW_j^2 - h) / 2 and I_jl + I_lj = dW_j dW_l.
        let diagonal = |j: usize| 0.5 * (increments[j] * increments[j] - h);
        assert!((areas[0] - diagonal(0)).abs() < 1e-15, "{}", areas[0]);
        assert!((areas[3] - diagonal(1)).abs() < 1e-15, "{}", areas[3]);
        let cross = increments[0] * increments[1];
        assert!((areas[1] + areas[2] - cross).abs() < 1e-15);
        // The increments are not redrawn: the areas depend on the seed only
        // through the bridge, so the diagonal is exact regardless of it.
        let other = call(
            "sde_levy_areas",
            &[Value::int(4_i32), real(h), dw, Value::int(8_i32)],
        );
        let other = floats(&other);
        assert!((other[0] - diagonal(0)).abs() < 1e-15);
        assert_ne!(bits(&areas[1..2]), bits(&other[1..2]));
    }

    #[test]
    fn scheme_metadata_is_the_documented_table() {
        let table = [
            ("euler_maruyama", [0.5, 1.0, 0.0, 0.0, 0.0]),
            ("milstein", [1.0, 1.0, 1.0, 0.0, 0.0]),
            ("milstein_stratonovich", [1.0, 0.0, 1.0, 0.0, 0.0]),
            ("milstein_noncommutative", [1.0, 1.0, 0.0, 0.0, 1.0]),
            ("srk_platen", [1.0, 1.0, 1.0, 0.0, 0.0]),
            ("heun_stratonovich", [1.0, 0.0, 0.0, 0.0, 0.0]),
            ("sde_adaptive", [1.0, 1.0, 1.0, 1.0, 0.0]),
            ("sde_adaptive_noncommutative", [1.0, 1.0, 0.0, 1.0, 1.0]),
        ];
        for (name, want) in table {
            let answer = call("sde_scheme_info", &[Value::String(name.to_owned())]);
            assert_eq!(floats(&answer), want, "{name}");
        }
        // An unknown name leaves the request unreduced; a non-string argument is
        // a type error. (The scheme name crosses as `Value::String`; rssn's
        // surface syntax has no string literal, so this operator is called
        // through the value ABI, not from parsed text.)
        let session = Session::new();
        let unknown = session
            .call_value(
                "sde_scheme_info",
                &[Value::String(String::from("no_such_scheme"))],
                &CallCtx::new(),
            )
            .unwrap_or_else(|e| panic!("{e}"));
        assert!(
            matches!(unknown, Value::Term(_)),
            "an unknown name must stay a request, got {unknown:?}"
        );
        for argument in [Value::Symbol(String::from("n")), Value::int(1_i32)] {
            let error = session
                .call_value("sde_scheme_info", &[argument], &CallCtx::new())
                .err()
                .unwrap_or_else(|| panic!("a non-string scheme name must be a type error"));
            assert!(
                matches!(error, crate::kernel::call::OpError::Type { .. }),
                "{error:?}"
            );
        }
    }

    #[test]
    fn symbolic_or_out_of_bounds_requests_stay_requests() {
        for (request, why) in [
            ("sde_wiener(n, 2, 3, 1/100)", "symbolic seed"),
            ("sde_wiener(1, 0, 3, 1/100)", "no Wiener processes"),
            ("sde_wiener(1, 2, 0, 1/100)", "no steps"),
            ("sde_wiener(1, 2, 3, 0)", "zero step length"),
            ("sde_levy_areas(1, 1/100, list(), 8)", "empty increments"),
            (
                "sde_levy_areas(1, 1/100, list(a, b), 8)",
                "symbolic increments",
            ),
            ("sde_levy_areas(1, 1/100, list(1, 2), 0)", "no bridge modes"),
        ] {
            let (text, reduced) = reduce_with(&[sde()], request, &[]);
            assert!(
                !reduced,
                "{request} ({why}) should stay a request, got {text}"
            );
        }
        // A plain numeric request does reduce.
        assert!(simplify(&[sde()], "sde_wiener(1, 2, 3, 1/100)").starts_with("list("));
    }
}
