//! Partial differential equations.
//!
//! `pdsolve(equation, u(x, t))` asks for the solution of a PDE for the
//! undetermined function `u`; derivatives are written as for ordinary
//! equations, `diff(u(x, t), t)`, `diff(diff(u(x, t), x), x)`. A third
//! argument is a list of conditions:
//!
//! * `u(x, 0) = f`, `u(0, t) = 0`, `u(L, t) = 0` — values on a coordinate
//!   line (initial or Dirichlet boundary conditions);
//! * `at(diff(u(x, t), t), t, 0) = g`, `at(diff(u(x, t), x), x, 0) = 0` —
//!   values of a derivative there (initial velocity, Neumann conditions).
//!
//! The answer is `u(x, t) = ...`: a closed form when one exists (checked by
//! substitution where it contains no arbitrary functions), an integral
//! representation (Green's functions, heat kernels, Kirchhoff's formula)
//! that the integration kernels evaluate further when they can, or a
//! Fourier series `sum(..., n, 1, oo)` whose coefficients are computed in
//! closed form. When the initial data is a finite combination of
//! eigenfunctions the series is replaced by the finite sum.
//!
//! The time variable of an evolution equation is `t` if it is an argument
//! of `u`, otherwise the last argument.
//!
//! | operator | value |
//! |---|---|
//! | `pdsolve(eq, u(...)[, conditions])`, `solve_pde(...)` | the solution, by whichever method applies |
//! | `pde_classify(eq, u(...))` | `list(type, order, dimension, linear, homogeneous, character, list(methods...))` |
//! | `pde_order(eq, u(...))` | the order |
//! | `solve_pde_by_characteristics(eq, u(x, y))` | first-order linear and quasi-linear equations |
//! | `solve_pde_by_separation_of_variables(eq, u(x, t), conditions)` | heat and wave equations on an interval, Laplace's equation on a rectangle or box |
//! | `solve_pde_by_greens_function(eq, u(...))` | Poisson and Helmholtz equations in free space |
//! | `solve_with_fourier_transform(eq, u(x, t), conditions)` | heat and Schrödinger equations on the line |
//! | `solve_wave_equation_1d_dalembert`, `solve_heat_equation_1d`, `solve_heat_equation_3d`, `solve_wave_equation_3d`, `solve_laplace_equation_2d`, `solve_laplace_equation_3d`, `solve_poisson_equation_2d`, `solve_poisson_equation_3d`, `solve_helmholtz_equation`, `solve_schrodinger_equation`, `solve_klein_gordon_equation`, `solve_burgers_equation`, `solve_second_order_pde` | one method each, same arguments as `pdsolve` |
//! | `at(f, x, a)` | `f` with `x = a`, once `f` is concrete in `x` |

use std::collections::HashMap;

use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::Facts;
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
use crate::graph::SymbolId;
use crate::graph::Tier;
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::neg;
use crate::rules::complex::build::pow;
use crate::rules::complex::build::powi;
use crate::rules::complex::build::sub;
use crate::rules::poly::best;
use crate::rules::poly::repr::from_term;
use crate::rules::poly::repr::to_term;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;
use crate::rules::poly::repr::Poly;

use super::calculus::antiderivative;
use super::calculus::derivative;
use super::complex::complex;
use super::logic::logic;
use super::ode::ode;
use super::solve::as_expression;
use super::solve::solve_for;
use super::special::special;

/// The partial differential equation rule set.
#[must_use]
pub fn pde() -> RuleSet {
    RuleSet::new("pde", install).needs(ode()).needs(special()).needs(complex()).needs(logic())
}

/// Which solver a request asks for.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Method {
    Any,
    Characteristics,
    Separation,
    Green,
    Fourier,
    Dalembert,
    Heat1,
    Heat3,
    Wave3,
    Laplace2,
    Laplace3,
    Poisson2,
    Poisson3,
    Helmholtz,
    Schrodinger,
    KleinGordon,
    Burgers,
    SecondOrder,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Request {
    Solve(Method),
    Classify,
    Order,
    At,
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let heavy = |name: &str, arity: Arity| OpDescriptor::new(name, arity).flags(OpFlags::HEAVY).cost(100);
    let at = i.op(heavy("at", Arity::Fixed(3)).flags(OpFlags::OPAQUE_ON_APPLY))?;
    i.kernel("pde/at", Tier::Reduce, Pde { op: at, request: Request::At });
    let table = [
        ("pdsolve", Method::Any),
        ("solve_pde", Method::Any),
        ("solve_pde_by_characteristics", Method::Characteristics),
        ("solve_pde_by_separation_of_variables", Method::Separation),
        ("solve_pde_by_greens_function", Method::Green),
        ("solve_with_fourier_transform", Method::Fourier),
        ("solve_wave_equation_1d_dalembert", Method::Dalembert),
        ("solve_heat_equation_1d", Method::Heat1),
        ("solve_heat_equation_3d", Method::Heat3),
        ("solve_wave_equation_3d", Method::Wave3),
        ("solve_laplace_equation_2d", Method::Laplace2),
        ("solve_laplace_equation_3d", Method::Laplace3),
        ("solve_poisson_equation_2d", Method::Poisson2),
        ("solve_poisson_equation_3d", Method::Poisson3),
        ("solve_helmholtz_equation", Method::Helmholtz),
        ("solve_schrodinger_equation", Method::Schrodinger),
        ("solve_klein_gordon_equation", Method::KleinGordon),
        ("solve_burgers_equation", Method::Burgers),
        ("solve_second_order_pde", Method::SecondOrder),
    ];
    for (name, method) in table {
        let op = i.op(heavy(name, Arity::Variadic))?;
        i.kernel(&format!("pde/{name}"), Tier::Reduce, Pde { op, request: Request::Solve(method) });
    }
    let op = i.op(heavy("pde_classify", Arity::Fixed(2)))?;
    i.kernel("pde/pde_classify", Tier::Reduce, Pde { op, request: Request::Classify });
    let op = i.op(heavy("pde_order", Arity::Fixed(2)))?;
    i.kernel("pde/pde_order", Tier::Reduce, Pde { op, request: Request::Order });
    Ok(())
}

struct Pde {
    op: OpId,
    request: Request,
}

impl Kernel for Pde {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        match self.request {
            | Request::At => at(cx, &args).map_or(Outcome::Pass, Outcome::Equal),
            | Request::Order => {
                let Some(problem) = args.get(..2).and_then(|a| Problem::parse(cx, a[0], a[1])) else {
                    return Outcome::Pass;
                };
                Outcome::Equal(cx.graph.int(i64::from(problem.order())))
            },
            | Request::Classify => {
                let Some(problem) = args.get(..2).and_then(|a| Problem::parse(cx, a[0], a[1])) else {
                    return Outcome::Pass;
                };
                classify_term(cx, &problem).map_or(Outcome::Pass, Outcome::Equal)
            },
            | Request::Solve(method) => {
                let (equation, unknown, conditions) = match *args.as_slice() {
                    | [e, u] => (e, u, None),
                    | [e, u, c] => (e, u, Some(c)),
                    | _ => return Outcome::Pass,
                };
                let Some(problem) = Problem::parse(cx, equation, unknown) else {
                    return Outcome::Pass;
                };
                let conditions = match conditions {
                    | Some(c) => match Conditions::parse(cx, &problem, c) {
                        | Some(c) => c,
                        | None => return Outcome::Pass,
                    },
                    | None => Conditions::default(),
                };
                solve(cx, &problem, &conditions, method).map_or(Outcome::Pass, Outcome::Pinned)
            },
        }
    }
}

/// `at(f, x, a)`: substitution once `f` is free of requests involving `x`.
fn at(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[f, x, a] = args else {
        return None;
    };
    let symbol = cx.graph.symbol_of(x)?;
    let term = best(cx.graph, f)?;
    let mut stack = vec![term];
    while let Some(n) = stack.pop() {
        let heavy = cx.graph.ops().get(cx.graph.op(n)).flags.has(OpFlags::HEAVY);
        if heavy && cx.graph.depends_on(cx.graph.find(n), symbol) {
            return None;
        }
        stack.extend_from_slice(cx.graph.children(n));
    }
    let substituted = cx.graph.substitute(term, x, a);
    Some(cx.simplify(substituted))
}

// ----------------------------------------------------------------------
// The problem: a linear decomposition of the equation
// ----------------------------------------------------------------------

/// A derivative `∂^α u`, as the count of differentiations per variable.
type Index = Vec<u32>;

struct Problem {
    /// `u(x, t, ...)`.
    unknown: NodeId,
    /// The function symbol `u`.
    function: NodeId,
    /// Independent variables, in the order of the arguments of `u`.
    vars: Vec<NodeId>,
    /// `lhs - rhs`, best form.
    residual: NodeId,
    /// Coefficient of each derivative occurring linearly.
    linear: Vec<(Index, NodeId)>,
    /// Terms free of `u`: the equation reads `Σ c_α ∂^α u + source = 0`.
    source: NodeId,
    /// Whether some term is not linear in `u` and its derivatives.
    nonlinear: bool,
    /// The derivative jets of `u` occurring, with their nodes.
    jets: Vec<(Index, NodeId)>,
}

impl Problem {
    fn parse(
        cx: &mut Cx<'_>,
        equation: NodeId,
        unknown: NodeId,
    ) -> Option<Self> {
        let unknown = best(cx.graph, unknown)?;
        if cx.graph.op(unknown) != core::APPLY {
            return None;
        }
        let children = cx.graph.children(unknown).to_vec();
        let (&function, vars) = children.split_first()?;
        if vars.is_empty() || vars.iter().any(|&v| cx.graph.symbol_of(v).is_none()) {
            return None;
        }
        let vars = vars.to_vec();
        let equation = best(cx.graph, equation)?;
        let residual = as_expression(cx.graph, equation);
        let residual = cx.simplify(residual);
        let diff = cx.graph.ops().lookup("diff")?;
        let mut problem = Self {
            unknown,
            function,
            vars,
            residual,
            linear: Vec::new(),
            source: NodeId::NONE,
            nonlinear: false,
            jets: Vec::new(),
        };
        problem.decompose(cx, diff)?;
        Some(problem)
    }

    /// If `node` is a derivative of `u`, its multi-index.
    fn index_of(
        &self,
        graph: &Graph,
        diff: OpId,
        node: NodeId,
    ) -> Option<Index> {
        if node == self.unknown {
            return Some(vec![0; self.vars.len()]);
        }
        if graph.op(node) != diff {
            return None;
        }
        let &[inner, x] = graph.children(node) else {
            return None;
        };
        let mut index = self.index_of(graph, diff, inner)?;
        let k = self.vars.iter().position(|&v| graph.same(v, x))?;
        index[k] = index[k].checked_add(1)?;
        Some(index)
    }

    fn decompose(
        &mut self,
        cx: &mut Cx<'_>,
        diff: OpId,
    ) -> Option<()> {
        let graph = &mut *cx.graph;
        let mut gens = Gens::default();
        let poly = from_term(graph, &mut gens, self.residual, Limits { terms: 512, exponent: 8 })?;
        // Generators that are jets of u.
        let mut jet_of: HashMap<u32, Index> = HashMap::new();
        for g in 0..u32::try_from(gens.len()).ok()? {
            let Some(node) = gens.node(g) else {
                continue;
            };
            if let Some(index) = self.index_of(graph, diff, node) {
                jet_of.insert(g, index.clone());
                if !self.jets.iter().any(|(i, _)| *i == index) {
                    self.jets.push((index, node));
                }
            } else if graph.depends_on(graph.find(node), graph.symbol_of(self.function)?) {
                // u hidden inside a function: sin(u), exp(u_x), ...
                self.nonlinear = true;
            }
        }
        let mut linear: Vec<(Index, Poly)> = Vec::new();
        let mut source = Poly::zero();
        for (mono, coeff) in poly.terms() {
            let mut jet = None;
            let mut rest = Vec::new();
            let mut degree = 0_u32;
            for &(g, e) in mono {
                if let Some(index) = jet_of.get(&g) {
                    degree = degree.saturating_add(e);
                    jet = Some(index.clone());
                } else {
                    rest.push((g, e));
                }
            }
            let term = Poly::monomial(rest, coeff.clone());
            match (jet, degree) {
                | (None, _) => source = source.add(&term),
                | (Some(index), 1) => match linear.iter_mut().find(|(i, _)| *i == index) {
                    | Some((_, c)) => *c = c.add(&term),
                    | None => linear.push((index, term)),
                },
                | _ => self.nonlinear = true,
            }
        }
        for (index, coeff) in linear {
            let c = to_term(cx.graph, &gens, &coeff);
            let c = cx.simplify(c);
            self.linear.push((index, c));
        }
        self.source = to_term(cx.graph, &gens, &source);
        self.source = cx.simplify(self.source);
        Some(())
    }

    fn order(&self) -> u32 {
        self.jets.iter().map(|(i, _)| i.iter().sum::<u32>()).max().unwrap_or(0)
    }

    const fn dimension(&self) -> usize {
        self.vars.len()
    }

    /// Coefficient of `∂^α u` (zero when absent).
    fn coefficient(
        &self,
        graph: &mut Graph,
        index: &[u32],
    ) -> NodeId {
        self.linear.iter().find(|(i, _)| i == index).map_or_else(|| graph.int(0), |&(_, c)| c)
    }

    fn unit(
        &self,
        k: usize,
        order: u32,
    ) -> Index {
        let mut index = vec![0; self.vars.len()];
        index[k] = order;
        index
    }

    /// Index of the time variable: `t` if present, else the last one.
    fn time(
        &self,
        graph: &Graph,
    ) -> usize {
        self.vars
            .iter()
            .position(|&v| graph.symbol_of(v).is_some_and(|s| graph.interner().symbol_name(s) == "t"))
            .unwrap_or_else(|| self.vars.len().saturating_sub(1))
    }

    fn homogeneous(
        &self,
        graph: &Graph,
    ) -> bool {
        graph.number_of(self.source).is_some_and(Number::is_zero)
    }

    /// Whether only the given derivatives occur (others have no term).
    fn only(
        &self,
        graph: &Graph,
        allowed: &[Index],
    ) -> bool {
        self.linear.iter().all(|(i, c)| allowed.contains(i) || graph.number_of(*c).is_some_and(Number::is_zero))
    }

    /// Whether `node` is free of every independent variable.
    fn constant(
        &self,
        graph: &Graph,
        node: NodeId,
    ) -> bool {
        self.vars.iter().all(|&v| graph.symbol_of(v).is_none_or(|s| !graph.depends_on(graph.find(node), s)))
    }
}

/// The time derivative order and spatial Laplacian structure of an
/// evolution equation `a ∂_t^k u = b Δu + c u + source`.
struct Evolution {
    time: usize,
    /// Order of the time derivative (1 or 2).
    order: u32,
    /// `b / a`: diffusivity or squared wave speed (per spatial variable,
    /// all equal).
    speed: NodeId,
    /// `c / a`: coefficient of `u` itself.
    potential: NodeId,
}

fn evolution(
    cx: &mut Cx<'_>,
    p: &Problem,
) -> Option<Evolution> {
    let time = p.time(cx.graph);
    let (mut order, mut a) = (0, NodeId::NONE);
    for k in [1, 2] {
        let index = p.unit(time, k);
        let c = p.coefficient(cx.graph, &index);
        if !cx.graph.number_of(c).is_some_and(Number::is_zero) {
            if order != 0 {
                return None;
            }
            order = k;
            a = c;
        }
    }
    if order == 0 || p.dimension() < 2 {
        return None;
    }
    let mut allowed = vec![p.unit(time, order), vec![0; p.dimension()]];
    let mut speed = None;
    for k in (0..p.dimension()).filter(|&k| k != time) {
        let index = p.unit(k, 2);
        let b = p.coefficient(cx.graph, &index);
        // a u_t + b' u_xx = 0  =>  u_t = (-b'/a) u_xx
        let a_inv = powi(cx.graph, a, -1);
        let ratio = mul(cx.graph, &[b, a_inv]);
        let ratio = neg(cx.graph, ratio);
        let ratio = cx.simplify(ratio);
        if cx.graph.number_of(ratio).is_some_and(Number::is_zero) {
            return None;
        }
        match speed {
            | None => speed = Some(ratio),
            | Some(s) => {
                let difference = sub(cx.graph, s, ratio);
                if !cx.is_zero(difference) {
                    return None;
                }
            },
        }
        allowed.push(index);
    }
    if !p.only(cx.graph, &allowed) || p.nonlinear {
        return None;
    }
    let speed = speed?;
    if order == 2 && cx.graph.facts(speed).has(Facts::NEGATIVE) {
        // u_tt = -u_xx is Laplace's equation, not a wave.
        return None;
    }
    let c = p.coefficient(cx.graph, &vec![0; p.dimension()]);
    let a_inv = powi(cx.graph, a, -1);
    let potential = mul(cx.graph, &[c, a_inv]);
    let potential = neg(cx.graph, potential);
    let potential = cx.simplify(potential);
    Some(Evolution { time, order, speed, potential })
}

// ----------------------------------------------------------------------
// Conditions
// ----------------------------------------------------------------------

/// `∂^α u = value` on the hyperplane `vars[on] = point`.
#[derive(Clone, Debug)]
struct Condition {
    on: usize,
    point: NodeId,
    derivative: Index,
    value: NodeId,
}

#[derive(Clone, Debug, Default)]
struct Conditions(Vec<Condition>);

impl Conditions {
    fn parse(
        cx: &mut Cx<'_>,
        p: &Problem,
        list: NodeId,
    ) -> Option<Self> {
        let list = best(cx.graph, list)?;
        if cx.graph.op(list) != core::LIST {
            return None;
        }
        let diff = cx.graph.ops().lookup("diff")?;
        let at = cx.graph.ops().lookup("at")?;
        let mut out = Vec::new();
        for condition in cx.graph.children(list).to_vec() {
            let &[target, value] = cx.graph.children(condition) else {
                return None;
            };
            if cx.graph.op(condition) != core::EQ {
                return None;
            }
            if cx.graph.op(target) == core::APPLY {
                // u(..., a, ...) = value: exactly one argument differs.
                let args = cx.graph.children(target).to_vec();
                if args.first() != Some(&p.function) || args.len() != p.vars.len() + 1 {
                    return None;
                }
                let differing: Vec<usize> = (0..p.vars.len()).filter(|&k| !cx.graph.same(args[k + 1], p.vars[k])).collect();
                let &[on] = differing.as_slice() else {
                    return None;
                };
                out.push(Condition { on, point: args[on + 1], derivative: vec![0; p.vars.len()], value });
            } else if cx.graph.op(target) == at {
                let &[inner, x, point] = cx.graph.children(target) else {
                    return None;
                };
                let derivative = p.index_of(cx.graph, diff, inner)?;
                let on = p.vars.iter().position(|&v| cx.graph.same(v, x))?;
                out.push(Condition { on, point, derivative, value });
            } else {
                return None;
            }
        }
        Some(Self(out))
    }

    fn find(
        &self,
        graph: &Graph,
        on: usize,
        derivative: &[u32],
        point: Option<NodeId>,
    ) -> Option<&Condition> {
        self.0.iter().find(|c| c.on == on && c.derivative == derivative && point.is_none_or(|p| graph.same(p, c.point)))
    }

    const fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

// ----------------------------------------------------------------------
// Classification
// ----------------------------------------------------------------------

fn symbol_term(
    graph: &mut Graph,
    name: &str,
) -> NodeId {
    graph.sym(name)
}

fn boolean(
    graph: &mut Graph,
    value: bool,
) -> NodeId {
    let op = graph.ops().lookup(if value { "true" } else { "false" });
    match op {
        | Some(op) => graph.node(op, &[]),
        | None => graph.int(i64::from(value)),
    }
}

/// The type of the equation, by the usual names.
fn kind(
    cx: &mut Cx<'_>,
    p: &Problem,
) -> &'static str {
    let order = p.order();
    let has_i = |graph: &Graph, n: NodeId| {
        let i = graph.ops().lookup("I");
        let mut stack = vec![n];
        while let Some(m) = stack.pop() {
            if Some(graph.op(m)) == i {
                return true;
            }
            stack.extend_from_slice(graph.children(m));
        }
        false
    };
    if order == 1 {
        return if p.nonlinear { "burgers" } else { "transport" };
    }
    if order != 2 {
        return "unknown";
    }
    if p.nonlinear {
        return "nonlinear";
    }
    let zero_index = vec![0; p.dimension()];
    let has_u = !{
        let c = p.coefficient(cx.graph, &zero_index);
        cx.graph.number_of(c).is_some_and(Number::is_zero)
    };
    if let Some(e) = evolution(cx, p) {
        let time = p.time(cx.graph);
        let a = p.coefficient(cx.graph, &p.unit(time, e.order));
        return match (e.order, has_i(cx.graph, a) || has_i(cx.graph, e.speed)) {
            | (1, true) => "schrodinger",
            | (1, false) => "heat",
            | (_, _) if has_u => "klein_gordon",
            | _ => "wave",
        };
    }
    // All pure second derivatives with one sign, no mixed or first ones.
    let mut signs = Vec::new();
    let mut allowed = vec![zero_index];
    for k in 0..p.dimension() {
        let index = p.unit(k, 2);
        let c = p.coefficient(cx.graph, &index);
        let facts = cx.graph.facts(c);
        signs.push(if facts.has(Facts::POSITIVE) { 1 } else if facts.has(Facts::NEGATIVE) { -1 } else { 0 });
        allowed.push(index);
    }
    if p.only(cx.graph, &allowed) && signs.iter().all(|&s| s != 0 && s == signs[0]) {
        return match (has_u, p.homogeneous(cx.graph)) {
            | (true, _) => "helmholtz",
            | (false, true) => "laplace",
            | (false, false) => "poisson",
        };
    }
    "unknown"
}

/// Hyperbolic, parabolic or elliptic by the discriminant of the principal
/// part (two variables), or by the signs of the pure second derivatives.
fn character(
    cx: &mut Cx<'_>,
    p: &Problem,
) -> &'static str {
    match p.order() {
        | 1 => return "first_order",
        | 2 => {},
        | _ => return "higher_order",
    }
    let n = p.dimension();
    let mut mixed = false;
    let mut signs = Vec::with_capacity(n);
    for k in 0..n {
        let c = p.coefficient(cx.graph, &p.unit(k, 2));
        let facts = cx.graph.facts(c);
        signs.push(if facts.has(Facts::POSITIVE) {
            1
        } else if facts.has(Facts::NEGATIVE) {
            -1
        } else if cx.graph.number_of(c).is_some_and(Number::is_zero) {
            0
        } else {
            2
        });
        for j in k + 1..n {
            let mut index = vec![0; n];
            index[k] = 1;
            index[j] = 1;
            let c = p.coefficient(cx.graph, &index);
            mixed |= !cx.graph.number_of(c).is_some_and(Number::is_zero);
        }
    }
    if n == 2 && mixed {
        // B² - 4AC with A, B, C the coefficients of u_xx, u_xy, u_yy.
        let a = p.coefficient(cx.graph, &[2, 0]);
        let b = p.coefficient(cx.graph, &[1, 1]);
        let c = p.coefficient(cx.graph, &[0, 2]);
        let b2 = powi(cx.graph, b, 2);
        let four = cx.graph.int(-4);
        let ac = mul(cx.graph, &[four, a, c]);
        let d = add(cx.graph, &[b2, ac]);
        let d = cx.simplify(d);
        let facts = cx.graph.facts(d);
        return if cx.graph.number_of(d).is_some_and(Number::is_zero) {
            "parabolic"
        } else if facts.has(Facts::POSITIVE) {
            "hyperbolic"
        } else if facts.has(Facts::NEGATIVE) {
            "elliptic"
        } else {
            "mixed"
        };
    }
    if mixed || signs.contains(&2) {
        return "mixed";
    }
    let positive = signs.iter().filter(|&&s| s == 1).count();
    let negative = signs.iter().filter(|&&s| s == -1).count();
    let zero = signs.iter().filter(|&&s| s == 0).count();
    match (positive, negative, zero) {
        | (_, 0, 0) | (0, _, 0) => "elliptic",
        | (1, _, 0) | (_, 1, 0) => "hyperbolic",
        | (_, _, 1) => "parabolic",
        | _ => "mixed",
    }
}

fn classify_term(
    cx: &mut Cx<'_>,
    p: &Problem,
) -> Option<NodeId> {
    let kind = kind(cx, p);
    let character = character(cx, p);
    let methods: &[&str] = match kind {
        | "transport" => &["characteristics"],
        | "burgers" => &["characteristics", "hopf_cole"],
        | "heat" => &["separation_of_variables", "fourier_transform", "heat_kernel"],
        | "wave" => &["dalembert", "separation_of_variables", "kirchhoff"],
        | "laplace" => &["separation_of_variables", "greens_function"],
        | "poisson" => &["greens_function"],
        | "helmholtz" => &["greens_function", "separation_of_variables"],
        | "schrodinger" => &["fourier_transform", "separation_of_variables"],
        | "klein_gordon" => &["plane_waves"],
        | _ => &[],
    };
    let graph = &mut *cx.graph;
    let kind_node = symbol_term(graph, kind);
    let order = graph.int(i64::from(p.order()));
    let dimension = graph.int(i64::try_from(p.dimension()).ok()?);
    let linear = boolean(graph, !p.nonlinear);
    let homogeneous = boolean(graph, p.homogeneous(graph));
    let character = symbol_term(graph, character);
    let methods: Vec<NodeId> = methods.iter().map(|m| symbol_term(graph, m)).collect();
    let methods = graph.node(core::LIST, &methods);
    Some(graph.node(core::LIST, &[kind_node, order, dimension, linear, homogeneous, character, methods]))
}

// ----------------------------------------------------------------------
// Solving
// ----------------------------------------------------------------------

fn solve(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
    method: Method,
) -> Option<NodeId> {
    use Method::{Any, Characteristics, Burgers, Separation, SecondOrder, Dalembert, Wave3, Heat1, Fourier, Heat3, Schrodinger, KleinGordon, Laplace2, Laplace3, Green, Poisson2, Poisson3, Helmholtz};
    let kind = kind(cx, p);
    let try_method = |m: Method| method == Any || method == m;
    let mut solution = if p.order() == 1 && (try_method(Characteristics) || try_method(Burgers)) {
        characteristics(cx, p, conditions)
    } else {
        None
    };
    if solution.is_none() && kind == "wave" {
        if !conditions.is_empty() && (try_method(Separation) || try_method(SecondOrder)) {
            solution = separation(cx, p, conditions);
        }
        if solution.is_none() && p.dimension() == 2 && (try_method(Dalembert) || try_method(SecondOrder)) {
            solution = dalembert(cx, p, conditions);
        }
        if solution.is_none() && p.dimension() == 4 && try_method(Wave3) {
            solution = kirchhoff(cx, p, conditions);
        }
        if solution.is_none() && p.dimension() == 4 && method == Any {
            solution = kirchhoff(cx, p, conditions);
        }
    }
    if solution.is_none() && kind == "heat" {
        if (try_method(Separation) || try_method(Heat1)) && p.dimension() == 2 {
            solution = separation(cx, p, conditions);
        }
        if solution.is_none() && (try_method(Fourier) || try_method(Heat1) || try_method(Heat3)) {
            solution = heat_kernel(cx, p, conditions);
        }
    }
    if solution.is_none() && kind == "schrodinger" && (try_method(Fourier) || try_method(Schrodinger)) {
        solution = schrodinger(cx, p, conditions);
    }
    if solution.is_none() && kind == "klein_gordon" && try_method(KleinGordon) {
        solution = klein_gordon(cx, p);
    }
    if solution.is_none() && kind == "laplace" {
        if try_method(Separation) || try_method(Laplace2) || try_method(Laplace3) {
            solution = laplace_box(cx, p, conditions);
        }
    }
    if solution.is_none() && matches!(kind, "poisson" | "laplace" | "helmholtz")
        && (try_method(Green) || try_method(Poisson2) || try_method(Poisson3) || try_method(Helmholtz))
    {
        solution = green(cx, p);
    }
    let solution = solution?;
    Some(cx.graph.node(core::EQ, &[p.unknown, solution]))
}

/// A fresh symbol named `stem` if that name is free in `term`, else a
/// numbered variant.
fn dummy(
    cx: &mut Cx<'_>,
    p: &Problem,
    stem: &str,
) -> (NodeId, SymbolId) {
    let graph = &mut *cx.graph;
    let candidate = graph.interner_mut().symbol(stem);
    let used = graph.depends_on(graph.find(p.residual), candidate)
        || p.vars.iter().any(|&v| graph.symbol_of(v) == Some(candidate));
    let symbol = if used { graph.interner_mut().fresh_symbol(stem) } else { candidate };
    (graph.symbol_node(symbol), symbol)
}

fn apply_function(
    graph: &mut Graph,
    name: &str,
    argument: NodeId,
) -> NodeId {
    let f = graph.sym(name);
    graph.node(core::APPLY, &[f, argument])
}

/// First-order equations: `a u_x + b u_y + c u = f` with constant `a`, `b`,
/// `c`; `a(x, y) u_x + b(x, y) u_y = 0` through the characteristic ODE;
/// and quasi-linear `u_t + a(u) u_x = 0` (Burgers) by its implicit
/// solution.
fn characteristics(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    if p.dimension() != 2 {
        return None;
    }
    if p.nonlinear {
        return quasilinear(cx, p, conditions);
    }
    let (x, y) = (p.vars[0], p.vars[1]);
    let a = p.coefficient(cx.graph, &[1, 0]);
    let b = p.coefficient(cx.graph, &[0, 1]);
    let c = p.coefficient(cx.graph, &[0, 0]);
    if !p.only(cx.graph, &[vec![1, 0], vec![0, 1], vec![0, 0]]) {
        return None;
    }
    let f = neg(cx.graph, p.source);
    let f = cx.simplify(f);
    let zero = |graph: &Graph, n: NodeId| graph.number_of(n).is_some_and(Number::is_zero);
    if p.constant(cx.graph, a) && p.constant(cx.graph, b) && p.constant(cx.graph, c) {
        // Swap the roles of x and y if u_x is absent.
        let (a, b, x, y) = if zero(cx.graph, a) { (b, a, y, x) } else { (a, b, x, y) };
        if zero(cx.graph, a) {
            return None;
        }
        let (xi, xi_symbol) = dummy(cx, p, "xi");
        let _ = xi_symbol;
        // Along y = (b x - ξ)/a, du/dx + (c/a) u = f/a.
        let inverse_a = powi(cx.graph, a, -1);
        let bx = mul(cx.graph, &[b, x]);
        let numerator = sub(cx.graph, bx, xi);
        let y_on = mul(cx.graph, &[numerator, inverse_a]);
        let f_on = cx.graph.substitute(f, y, y_on);
        let rate = mul(cx.graph, &[c, inverse_a]);
        let rate = cx.simplify(rate);
        let exp = cx.graph.ops().lookup("exp")?;
        let rate_x = mul(cx.graph, &[rate, x]);
        let growth = cx.graph.node(exp, &[rate_x]);
        let decay = {
            let minus = neg(cx.graph, rate_x);
            cx.graph.node(exp, &[minus])
        };
        let integrand = mul(cx.graph, &[f_on, inverse_a, growth]);
        let integrand = cx.simplify(integrand);
        let particular = if zero(cx.graph, integrand) { cx.graph.int(0) } else { antiderivative(cx, integrand, x)? };
        // ξ = b x - a y
        let ay = mul(cx.graph, &[a, y]);
        let invariant = sub(cx.graph, bx, ay);
        let invariant = cx.simplify(invariant);
        let particular = cx.graph.substitute(particular, xi, invariant);
        let arbitrary = apply_function(cx.graph, "F", invariant);
        let inside = add(cx.graph, &[particular, arbitrary]);
        let solution = mul(cx.graph, &[decay, inside]);
        let solution = cx.simplify(solution);
        if !verified(cx, p, solution) {
            return None;
        }
        return apply_initial_data(cx, p, conditions, solution, invariant);
    }
    // Variable coefficients, homogeneous transport: dy/dx = b/a.
    if !zero(cx.graph, c) || !zero(cx.graph, f) || zero(cx.graph, a) {
        return None;
    }
    let a_inv = powi(cx.graph, a, -1);
    let slope = mul(cx.graph, &[b, a_inv]);
    let slope = cx.simplify(slope);
    let y_symbol = cx.graph.symbol_of(y)?;
    let big_y = {
        let name = format!("{}_char", cx.graph.interner().symbol_name(y_symbol));
        let function = cx.graph.sym(&name);
        cx.graph.node(core::APPLY, &[function, x])
    };
    let slope_in_y = cx.graph.substitute(slope, y, big_y);
    let diff = cx.graph.ops().lookup("diff")?;
    let dsolve = cx.graph.ops().lookup("dsolve")?;
    let dy = cx.graph.node(diff, &[big_y, x]);
    let ode = cx.graph.node(core::EQ, &[dy, slope_in_y]);
    let request = cx.graph.node(dsolve, &[ode, big_y]);
    // The ODE is solved by its own kernel in a nested run.
    let solved = cx.simplify(request);
    if cx.graph.op(solved) != core::EQ {
        return None;
    }
    let relation = as_expression(cx.graph, solved);
    let relation = cx.graph.replace_subterm(relation, big_y, y);
    let constant = cx.graph.sym("C1");
    let invariant = *solve_for(cx.graph, relation, constant, 0)?.first()?;
    let invariant = cx.simplify(invariant);
    if !p.constant(cx.graph, invariant) {
        let solution = apply_function(cx.graph, "F", invariant);
        return apply_initial_data(cx, p, conditions, solution, invariant);
    }
    None
}

/// With a condition `u(x, y0) = g(x)` on a line, replace the arbitrary
/// function `F(ξ)` by the value that matches it.
fn apply_initial_data(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
    solution: NodeId,
    invariant: NodeId,
) -> Option<NodeId> {
    let Some(condition) = conditions.0.first() else {
        return Some(solution);
    };
    if conditions.0.len() != 1 || condition.derivative.iter().any(|&d| d != 0) {
        return None;
    }
    // On the line vars[on] = point, the solution with F(ξ) unknown must
    // equal the data. Write F(s) = data at the point of the line where
    // ξ = s, found by solving ξ(.., point, ..) = s for the free variable.
    let free = 1 - condition.on;
    let (s, _) = dummy(cx, p, "s");
    let on_line = cx.graph.substitute(invariant, p.vars[condition.on], condition.point);
    let equation = sub(cx.graph, on_line, s);
    let where_ = *solve_for(cx.graph, equation, p.vars[free], 0)?.first()?;
    let f_symbol = cx.graph.sym("F");
    // u on the line, with F(s) as an unknown: solve for F(s).
    let marker = cx.graph.sym("F_value");
    let arbitrary = cx.graph.node(core::APPLY, &[f_symbol, invariant]);
    let with_marker = cx.graph.replace_subterm(solution, arbitrary, marker);
    let line = cx.graph.substitute(with_marker, p.vars[condition.on], condition.point);
    let line = cx.graph.substitute(line, p.vars[free], where_);
    let value = cx.graph.substitute(condition.value, p.vars[free], where_);
    let matched = sub(cx.graph, line, value);
    let f_of_s = *solve_for(cx.graph, matched, marker, 0)?.first()?;
    let f_of_invariant = cx.graph.substitute(f_of_s, s, invariant);
    let result = cx.graph.substitute(with_marker, marker, f_of_invariant);
    let result = cx.simplify(result);
    verified(cx, p, result).then_some(result)
}

/// `u_t + a(u) u_x = 0` (with any constant factor on `u_t`): the implicit
/// solution `u = f(x - a(u) t)`.
fn quasilinear(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    let time = p.time(cx.graph);
    let space = 1 - time;
    let (u_t, u_x) = (p.unit(time, 1), p.unit(space, 1));
    let t_node = p.jets.iter().find(|(i, _)| *i == u_t)?.1;
    let x_node = p.jets.iter().find(|(i, _)| *i == u_x)?.1;
    // residual = A u_t + B(u) u_x with A free of u and both derivatives
    // occurring to the first power.
    let d_t = derivative_by_node(cx, p.residual, t_node)?;
    let d_x = derivative_by_node(cx, p.residual, x_node)?;
    let u_symbol = cx.graph.symbol_of(p.function)?;
    let depends_u = |graph: &Graph, n: NodeId| graph.depends_on(graph.find(n), u_symbol);
    if depends_u(cx.graph, d_t) || !p.homogeneous(cx.graph) {
        return None;
    }
    let d_t_inv = powi(cx.graph, d_t, -1);
    let speed = mul(cx.graph, &[d_x, d_t_inv]);
    let speed = cx.simplify(speed);
    // The speed may depend on u only (not on its derivatives).
    let (x, t) = (p.vars[space], p.vars[time]);
    let rebuilt = {
        let a = mul(cx.graph, &[d_t, t_node]);
        let b = mul(cx.graph, &[d_x, x_node]);
        add(cx.graph, &[a, b])
    };
    let check = sub(cx.graph, rebuilt, p.residual);
    if !cx.is_zero(check) {
        return None;
    }
    if p.jets.iter().any(|(i, n)| *i != u_t && *i != u_x && *n != p.unknown && depends_u(cx.graph, speed)) {
        return None;
    }
    let st = mul(cx.graph, &[speed, t]);
    let characteristic = sub(cx.graph, x, st);
    let data = match conditions.0.as_slice() {
        | [] => None,
        | [c] if c.on == time && c.derivative.iter().all(|&d| d == 0) => Some(c),
        | _ => return None,
    };
    let rhs = match data {
        | Some(c) => {
            if !cx.graph.number_of(c.point).is_some_and(Number::is_zero) {
                return None;
            }
            cx.graph.substitute(c.value, x, characteristic)
        },
        | None => apply_function(cx.graph, "F", characteristic),
    };
    Some(cx.simplify(rhs))
}

/// `∂ term / ∂ node` for a compound `node`, by freezing it as a symbol.
fn derivative_by_node(
    cx: &mut Cx<'_>,
    term: NodeId,
    node: NodeId,
) -> Option<NodeId> {
    let fresh = cx.graph.interner_mut().fresh_symbol("jet");
    let s = cx.graph.symbol_node(fresh);
    let frozen = cx.graph.replace_subterm(term, node, s);
    let d = derivative(cx.graph, frozen, s)?;
    let d = cx.simplify(d);
    Some(cx.graph.substitute(d, s, node))
}

/// The time and space variables of a two-variable evolution equation.
fn time_space(
    p: &Problem,
    e: &Evolution,
) -> Option<(NodeId, NodeId, usize, usize)> {
    if p.dimension() != 2 {
        return None;
    }
    let space = 1 - e.time;
    Some((p.vars[e.time], p.vars[space], e.time, space))
}

/// `u_tt = c² u_xx` on the line: `F(x - ct) + G(x + ct)`, or d'Alembert's
/// formula with `u(x, 0) = f`, `u_t(x, 0) = g`.
fn dalembert(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    let e = evolution(cx, p)?;
    if e.order != 2 || !p.homogeneous(cx.graph) || !cx.graph.number_of(e.potential).is_some_and(Number::is_zero) {
        return None;
    }
    let (t, x, time, _) = time_space(p, &e)?;
    let half = cx.graph.num(Number::fraction(1, 2)?);
    let c = pow(cx.graph, e.speed, half);
    let c = cx.simplify(c);
    let ct = mul(cx.graph, &[c, t]);
    let left = sub(cx.graph, x, ct);
    let right = add(cx.graph, &[x, ct]);
    if conditions.is_empty() {
        let f = apply_function(cx.graph, "F", left);
        let g = apply_function(cx.graph, "G", right);
        return Some(add(cx.graph, &[f, g]));
    }
    let zero_index = vec![0; p.dimension()];
    let position = conditions.find(cx.graph, time, &zero_index, None)?;
    let velocity = conditions.find(cx.graph, time, &p.unit(time, 1), None);
    if !cx.graph.number_of(position.point).is_some_and(Number::is_zero) {
        return None;
    }
    let f = position.value;
    let at_left = cx.graph.substitute(f, x, left);
    let at_right = cx.graph.substitute(f, x, right);
    let sum = add(cx.graph, &[at_left, at_right]);
    let mut terms = vec![mul(cx.graph, &[half, sum])];
    if let Some(v) = velocity {
        let (s, _) = dummy(cx, p, "s");
        let g = cx.graph.substitute(v.value, x, s);
        let defint = cx.graph.ops().lookup("defint")?;
        let integral = cx.graph.node(defint, &[g, s, left, right]);
        let two = cx.graph.int(2);
        let two_c = mul(cx.graph, &[two, c]);
        let scale = powi(cx.graph, two_c, -1);
        terms.push(mul(cx.graph, &[scale, integral]));
    }
    let solution = add(cx.graph, &terms);
    Some(cx.simplify(solution))
}

/// The heat equation `u_t = k Δu` in all of space with `u(·, 0) = f`: the
/// convolution with the heat kernel. One, two or three space dimensions.
fn heat_kernel(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    let e = evolution(cx, p)?;
    if e.order != 1 || !p.homogeneous(cx.graph) || !cx.graph.number_of(e.potential).is_some_and(Number::is_zero) {
        return None;
    }
    let zero_index = vec![0; p.dimension()];
    let initial = conditions.find(cx.graph, e.time, &zero_index, None)?;
    if !cx.graph.number_of(initial.point).is_some_and(Number::is_zero) || conditions.0.len() != 1 {
        return None;
    }
    convolve_gaussian(cx, p, e.time, e.speed, initial.value, false)
}

/// `∫ f(s) K(x - s, t) ds` over all space for the Gaussian kernel
/// `K = (4π k t)^(-d/2) exp(-|x - s|² / (4 k t))`. With `quantum`, the
/// free Schrödinger propagator, whose `k` is imaginary.
fn convolve_gaussian(
    cx: &mut Cx<'_>,
    p: &Problem,
    time: usize,
    k: NodeId,
    f: NodeId,
    quantum: bool,
) -> Option<NodeId> {
    let t = p.vars[time];
    let space: Vec<usize> = (0..p.dimension()).filter(|&j| j != time).collect();
    let d = i64::try_from(space.len()).ok()?;
    let graph_pi = {
        let pi = cx.graph.ops().lookup("pi")?;
        cx.graph.node(pi, &[])
    };
    let exp = cx.graph.ops().lookup("exp")?;
    let defint = cx.graph.ops().lookup("defint")?;
    let infinity = cx.graph.ops().lookup("oo")?;
    let four = cx.graph.int(4);
    let four_kt = mul(cx.graph, &[four, k, t]);
    let mut integrand = f;
    let mut squares = Vec::new();
    let mut dummies = Vec::new();
    for (n, &j) in space.iter().enumerate() {
        let stem = ["s", "r", "q"][n.min(2)];
        let (s, _) = dummy(cx, p, stem);
        integrand = cx.graph.substitute(integrand, p.vars[j], s);
        let difference = sub(cx.graph, p.vars[j], s);
        squares.push(powi(cx.graph, difference, 2));
        dummies.push(s);
    }
    let r2 = add(cx.graph, &squares);
    let exponent = {
        let four_kt_inv = powi(cx.graph, four_kt, -1);
        let ratio = mul(cx.graph, &[r2, four_kt_inv]);
        neg(cx.graph, ratio)
    };
    let gaussian = cx.graph.node(exp, &[exponent]);
    let normalisation = {
        let base = mul(cx.graph, &[graph_pi, four_kt]);
        let e = cx.graph.num(Number::fraction(-d, 2)?);
        pow(cx.graph, base, e)
    };
    let _ = quantum;
    let mut body = mul(cx.graph, &[normalisation, gaussian, integrand]);
    let oo = cx.graph.node(infinity, &[]);
    let minus_oo = neg(cx.graph, oo);
    for &s in dummies.iter().rev() {
        body = cx.graph.node(defint, &[body, s, minus_oo, oo]);
    }
    Some(body)
}

/// `i ħ ψ_t = -ħ²/(2m) ψ_xx` (any constant in place of `ħ²/2m`) with
/// `ψ(x, 0) = f`: the free propagator; without data, the separated
/// stationary states `exp(i (k x - ω t))`.
fn schrodinger(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    let e = evolution(cx, p)?;
    if e.order != 1 || !p.homogeneous(cx.graph) || !cx.graph.number_of(e.potential).is_some_and(Number::is_zero) {
        return None;
    }
    let zero_index = vec![0; p.dimension()];
    match conditions.find(cx.graph, e.time, &zero_index, None) {
        | Some(initial) if conditions.0.len() == 1 => {
            if !cx.graph.number_of(initial.point).is_some_and(Number::is_zero) {
                return None;
            }
            convolve_gaussian(cx, p, e.time, e.speed, initial.value, true)
        },
        | Some(_) => None,
        | None => {
            // Plane wave: ψ_t = κ ψ_xx with ψ = exp(i k x + κ (i k)² t)
            if p.dimension() != 2 {
                return None;
            }
            let (t, x, _, _) = time_space(p, &e)?;
            let i = cx.graph.ops().lookup("I")?;
            let i = cx.graph.node(i, &[]);
            let (k, _) = dummy(cx, p, "k");
            let ik = mul(cx.graph, &[i, k]);
            let ikx = mul(cx.graph, &[ik, x]);
            let k2 = powi(cx.graph, k, 2);
            let minus_one = cx.graph.int(-1);
            let rate = mul(cx.graph, &[minus_one, e.speed, k2, t]);
            let exponent = add(cx.graph, &[ikx, rate]);
            let exp = cx.graph.ops().lookup("exp")?;
            let wave = cx.graph.node(exp, &[exponent]);
            let amplitude = apply_function(cx.graph, "A", k);
            let body = mul(cx.graph, &[amplitude, wave]);
            let defint = cx.graph.ops().lookup("defint")?;
            let infinity = cx.graph.ops().lookup("oo")?;
            let oo = cx.graph.node(infinity, &[]);
            let minus_oo = neg(cx.graph, oo);
            Some(cx.graph.node(defint, &[body, k, minus_oo, oo]))
        },
    }
}

/// `φ_tt = c² Δφ - μ² φ` (Klein–Gordon): the superposition of plane waves
/// with `ω² = c² k² + μ²`.
fn klein_gordon(
    cx: &mut Cx<'_>,
    p: &Problem,
) -> Option<NodeId> {
    let e = evolution(cx, p)?;
    if e.order != 2 || !p.homogeneous(cx.graph) || p.dimension() != 2 {
        return None;
    }
    let (t, x, _, _) = time_space(p, &e)?;
    let (k, _) = dummy(cx, p, "k");
    let k2 = powi(cx.graph, k, 2);
    let c2k2 = mul(cx.graph, &[e.speed, k2]);
    let mass = neg(cx.graph, e.potential);
    let omega2 = add(cx.graph, &[c2k2, mass]);
    let half = cx.graph.num(Number::fraction(1, 2)?);
    let omega = pow(cx.graph, omega2, half);
    let omega = cx.simplify(omega);
    let i = cx.graph.ops().lookup("I")?;
    let i = cx.graph.node(i, &[]);
    let exp = cx.graph.ops().lookup("exp")?;
    let kx = mul(cx.graph, &[k, x]);
    let wt = mul(cx.graph, &[omega, t]);
    let wave = |cx: &mut Cx<'_>, sign: i64, amplitude: &str| {
        let s = cx.graph.int(sign);
        let phase_t = mul(cx.graph, &[s, wt]);
        let phase = add(cx.graph, &[kx, phase_t]);
        let exponent = mul(cx.graph, &[i, phase]);
        let e = cx.graph.node(exp, &[exponent]);
        let a = apply_function(cx.graph, amplitude, k);
        mul(cx.graph, &[a, e])
    };
    let outgoing = wave(cx, -1, "A");
    let incoming = wave(cx, 1, "B");
    let body = add(cx.graph, &[outgoing, incoming]);
    let defint = cx.graph.ops().lookup("defint")?;
    let infinity = cx.graph.ops().lookup("oo")?;
    let oo = cx.graph.node(infinity, &[]);
    let minus_oo = neg(cx.graph, oo);
    Some(cx.graph.node(defint, &[body, k, minus_oo, oo]))
}

/// Kirchhoff's formula for `u_tt = c² Δu` in three dimensions with
/// `u(x, 0) = f`, `u_t(x, 0) = g`:
/// `u = ∂_t (t M[f](ct)) + t M[g](ct)` with the spherical mean
/// `M[h](r) = (1/4π) ∫∫ h(x + r ω) sin θ dφ dθ`.
fn kirchhoff(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    let e = evolution(cx, p)?;
    if e.order != 2 || p.dimension() != 4 || !p.homogeneous(cx.graph) {
        return None;
    }
    let zero_index = vec![0; p.dimension()];
    let f = conditions.find(cx.graph, e.time, &zero_index, None).map(|c| c.value);
    let g = conditions.find(cx.graph, e.time, &p.unit(e.time, 1), None).map(|c| c.value);
    if f.is_none() && g.is_none() {
        return None;
    }
    let t = p.vars[e.time];
    let space: Vec<NodeId> = (0..4).filter(|&j| j != e.time).map(|j| p.vars[j]).collect();
    let half = cx.graph.num(Number::fraction(1, 2)?);
    let c = pow(cx.graph, e.speed, half);
    let r = mul(cx.graph, &[c, t]);
    let (theta, _) = dummy(cx, p, "theta");
    let (phi, _) = dummy(cx, p, "phi");
    let (sin, cos) = (cx.graph.ops().lookup("sin")?, cx.graph.ops().lookup("cos")?);
    let pi = cx.graph.ops().lookup("pi")?;
    let pi = cx.graph.node(pi, &[]);
    let defint = cx.graph.ops().lookup("defint")?;
    let st = cx.graph.node(sin, &[theta]);
    let ct = cx.graph.node(cos, &[theta]);
    let sp = cx.graph.node(sin, &[phi]);
    let cp = cx.graph.node(cos, &[phi]);
    let directions = [mul(cx.graph, &[st, cp]), mul(cx.graph, &[st, sp]), ct];
    let mean = |cx: &mut Cx<'_>, h: NodeId| -> Option<NodeId> {
        let mut shifted = h;
        // Substitute all three coordinates at once through fresh symbols.
        let fresh: Vec<NodeId> = (0..3)
            .map(|_| {
                let s = cx.graph.interner_mut().fresh_symbol("w");
                cx.graph.symbol_node(s)
            })
            .collect();
        for (k, &x) in space.iter().enumerate() {
            shifted = cx.graph.substitute(shifted, x, fresh[k]);
        }
        for (k, &x) in space.iter().enumerate() {
            let step = mul(cx.graph, &[r, directions[k]]);
            let moved = add(cx.graph, &[x, step]);
            shifted = cx.graph.substitute(shifted, fresh[k], moved);
        }
        // Expanded, the integrand is a sum of trigonometric monomials.
        let mut gens = Gens::default();
        if let Some(poly) = from_term(cx.graph, &mut gens, shifted, Limits { terms: 4096, exponent: 16 }) {
            shifted = to_term(cx.graph, &gens, &poly);
        }
        let body = mul(cx.graph, &[shifted, st]);
        let two = cx.graph.int(2);
        let two_pi = mul(cx.graph, &[two, pi]);
        let zero = cx.graph.int(0);
        let inner = cx.graph.node(defint, &[body, phi, zero, two_pi]);
        let outer = cx.graph.node(defint, &[inner, theta, zero, pi]);
        let four = cx.graph.int(4);
        let four_pi = mul(cx.graph, &[four, pi]);
        let scale = powi(cx.graph, four_pi, -1);
        let value = mul(cx.graph, &[scale, outer]);
        Some(cx.simplify(value))
    };
    let mut terms = Vec::new();
    if let Some(f) = f {
        let m = mean(cx, f)?;
        let tm = mul(cx.graph, &[t, m]);
        let tm = cx.simplify(tm);
        terms.push(derivative(cx.graph, tm, t)?);
    }
    if let Some(g) = g {
        let m = mean(cx, g)?;
        terms.push(mul(cx.graph, &[t, m]));
    }
    let sum = add(cx.graph, &terms);
    Some(cx.simplify(sum))
}

// ----------------------------------------------------------------------
// Separation of variables
// ----------------------------------------------------------------------

/// Boundary type at one end of an interval.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Boundary {
    Dirichlet,
    Neumann,
}

/// The eigenfunctions of `-X'' = λ X` on `[0, L]` for the two boundary
/// types, as `(sqrt(λ_n), X_n)` in the mode number `n`, together with the
/// first mode number (0 for the constant Neumann mode).
fn eigenfunctions(
    cx: &mut Cx<'_>,
    left: Boundary,
    right: Boundary,
    length: NodeId,
    x: NodeId,
    n: NodeId,
) -> Option<(NodeId, NodeId, i64)> {
    let pi = cx.graph.ops().lookup("pi")?;
    let pi = cx.graph.node(pi, &[]);
    let (sin, cos) = (cx.graph.ops().lookup("sin")?, cx.graph.ops().lookup("cos")?);
    let inverse_l = powi(cx.graph, length, -1);
    let (index, function, first) = match (left, right) {
        | (Boundary::Dirichlet, Boundary::Dirichlet) => (n, sin, 1),
        | (Boundary::Neumann, Boundary::Neumann) => (n, cos, 0),
        | (Boundary::Dirichlet, Boundary::Neumann) => {
            let half = cx.graph.num(Number::fraction(-1, 2)?);
            (add(cx.graph, &[n, half]), sin, 1)
        },
        | (Boundary::Neumann, Boundary::Dirichlet) => {
            let half = cx.graph.num(Number::fraction(-1, 2)?);
            (add(cx.graph, &[n, half]), cos, 1)
        },
    };
    let k = mul(cx.graph, &[index, pi, inverse_l]);
    let kx = mul(cx.graph, &[k, x]);
    let mode = cx.graph.node(function, &[kx]);
    Some((k, mode, first))
}

/// Projection coefficient `∫ h X_n dx / ∫ X_n² dx` over `[0, L]`.
fn coefficient(
    cx: &mut Cx<'_>,
    h: NodeId,
    mode: NodeId,
    x: NodeId,
    length: NodeId,
) -> Option<NodeId> {
    let defint = cx.graph.ops().lookup("defint")?;
    let zero = cx.graph.int(0);
    let product = mul(cx.graph, &[h, mode]);
    let numerator = cx.graph.node(defint, &[product, x, zero, length]);
    let square = powi(cx.graph, mode, 2);
    let denominator = cx.graph.node(defint, &[square, x, zero, length]);
    let denominator_inv = powi(cx.graph, denominator, -1);
    let ratio = mul(cx.graph, &[numerator, denominator_inv]);
    Some(cx.simplify(ratio))
}

/// Heat and wave equations on `[0, L]` with Dirichlet or Neumann ends,
/// and Laplace's equation (see [`laplace_box`]).
fn separation(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    let e = evolution(cx, p)?;
    if !p.homogeneous(cx.graph) || !cx.graph.number_of(e.potential).is_some_and(Number::is_zero) {
        return None;
    }
    let (t, x, time, space) = time_space(p, &e)?;
    let zero_index = vec![0; p.dimension()];
    // Boundary conditions at x = 0 and x = L, all homogeneous.
    let mut ends = Vec::new();
    for c in &conditions.0 {
        if c.on != space {
            continue;
        }
        if !cx.graph.number_of(c.value).is_some_and(Number::is_zero) {
            return None;
        }
        let kind = if c.derivative == zero_index {
            Boundary::Dirichlet
        } else if c.derivative == p.unit(space, 1) {
            Boundary::Neumann
        } else {
            return None;
        };
        ends.push((c.point, kind));
    }
    let [(a, ka), (b, kb)] = ends.as_slice() else {
        return None;
    };
    let (left, right, length) = if cx.graph.number_of(*a).is_some_and(Number::is_zero) {
        (*ka, *kb, *b)
    } else if cx.graph.number_of(*b).is_some_and(Number::is_zero) {
        (*kb, *ka, *a)
    } else {
        return None;
    };
    let initial = conditions.find(cx.graph, time, &zero_index, None)?;
    if !cx.graph.number_of(initial.point).is_some_and(Number::is_zero) {
        return None;
    }
    let velocity = conditions.find(cx.graph, time, &p.unit(time, 1), None).map(|c| c.value);
    if e.order == 1 && velocity.is_some() {
        return None;
    }
    // Mode number: an integer symbol.
    let (n, n_symbol) = dummy(cx, p, "n");
    cx.graph.assume(n_symbol, Facts::INTEGER | Facts::NONNEGATIVE);
    let (k, mode, first) = eigenfunctions(cx, left, right, length, x, n)?;
    let exp = cx.graph.ops().lookup("exp")?;
    let (sin, cos) = (cx.graph.ops().lookup("sin")?, cx.graph.ops().lookup("cos")?);
    // Time factor for mode n with coefficients A (from f) and B (from g).
    let k2 = powi(cx.graph, k, 2);
    let time_factor = |cx: &mut Cx<'_>, a: NodeId, b: Option<NodeId>| -> Option<NodeId> {
        if e.order == 1 {
            let minus_one = cx.graph.int(-1);
            let rate = mul(cx.graph, &[minus_one, e.speed, k2, t]);
            let decay = cx.graph.node(exp, &[rate]);
            Some(mul(cx.graph, &[a, decay]))
        } else {
            let half = cx.graph.num(Number::fraction(1, 2)?);
            let c = pow(cx.graph, e.speed, half);
            let omega = mul(cx.graph, &[c, k]);
            let wt = mul(cx.graph, &[omega, t]);
            let cosine = cx.graph.node(cos, &[wt]);
            let mut terms = vec![mul(cx.graph, &[a, cosine])];
            if let Some(b) = b {
                let sine = cx.graph.node(sin, &[wt]);
                let scale = powi(cx.graph, omega, -1);
                terms.push(mul(cx.graph, &[b, scale, sine]));
            }
            Some(add(cx.graph, &terms))
        }
    };
    // Finite data: a combination of modes with literal mode numbers.
    if let Some(finite) = finite_modes(cx, p, initial.value, velocity, mode, n, first, &time_factor) {
        if verified(cx, p, finite) && satisfies(cx, p, conditions, finite) {
            return Some(finite);
        }
    }
    let term_for = |cx: &mut Cx<'_>, mode: NodeId| -> Option<NodeId> {
        let a = coefficient(cx, initial.value, mode, x, length)?;
        let b = match velocity {
            | Some(g) => Some(coefficient(cx, g, mode, x, length)?),
            | None => None,
        };
        let factor = time_factor(cx, a, b)?;
        let term = mul(cx.graph, &[factor, mode]);
        Some(cx.simplify(term))
    };
    let general = term_for(cx, mode)?;
    let exact = |cx: &mut Cx<'_>, m: i64| -> Option<NodeId> {
        let number = cx.graph.int(m);
        let mode_m = cx.graph.substitute(mode, n, number);
        let mode_m = cx.simplify(mode_m);
        let term = term_for(cx, mode_m)?;
        let term = cx.graph.substitute(term, n, number);
        Some(cx.simplify(term))
    };
    assemble(cx, n, first, general, &exact)
}

/// `Σ_{n ≥ first} term(n)` from the general term with a symbolic mode
/// number. Integrating with a symbolic `n` silently assumes `n` differs
/// from the mode numbers where a denominator such as `n - 1` vanishes;
/// the first few terms are therefore also computed with `n` substituted
/// *before* integrating, and wherever the two disagree the exact terms
/// are written out and the series starts after them.
fn assemble(
    cx: &mut Cx<'_>,
    n: NodeId,
    first: i64,
    general: NodeId,
    exact: &dyn Fn(&mut Cx<'_>, i64) -> Option<NodeId>,
) -> Option<NodeId> {
    const CHECK: i64 = 8;
    let mut explicit = Vec::new();
    let mut last_exceptional = None;
    for m in first..first + CHECK {
        let exact_m = exact(cx, m)?;
        let number = cx.graph.int(m);
        let general_m = cx.graph.substitute(general, n, number);
        let general_m = cx.simplify(general_m);
        if !agree(cx, general_m, exact_m) {
            last_exceptional = Some(m);
        }
        explicit.push(exact_m);
    }
    let sum = cx.graph.ops().lookup("sum")?;
    let infinity = cx.graph.ops().lookup("oo")?;
    let oo = cx.graph.node(infinity, &[]);
    let general_vanishes = cx.graph.number_of(general).is_some_and(Number::is_zero);
    let (mut terms, start) = match last_exceptional {
        | None => (Vec::new(), first),
        | Some(m) => (explicit[..=usize::try_from(m - first).ok()?].to_vec(), m + 1),
    };
    if !general_vanishes {
        let start = cx.graph.int(start);
        terms.push(cx.graph.node(sum, &[general, n, start, oo]));
    }
    let total = add(cx.graph, &terms);
    Some(cx.simplify(total))
}

/// Whether two terms agree at generic values of their symbols.
fn agree(
    cx: &mut Cx<'_>,
    a: NodeId,
    b: NodeId,
) -> bool {
    let difference = sub(cx.graph, a, b);
    let difference = cx.simplify(difference);
    if cx.graph.number_of(difference).is_some_and(Number::is_zero) {
        return true;
    }
    let symbols = cx.graph.free_symbols(cx.graph.find(difference)).to_vec();
    let mut env = Env::numeric(0.0);
    for &s in &symbols {
        env.bind(s, 0.37 + 0.19 * f64::from(s.raw() % 11));
    }
    let scale = cx.graph.eval(b, &env).map_or(1.0, |v| v.abs().max(1.0));
    cx.graph.eval(difference, &env).is_some_and(|v| v.is_finite() && v.abs() <= 1e-9 * scale)
}

/// Builds the time-dependent factor of a mode from its coefficients.
type TimeFactor<'a> = dyn Fn(&mut Cx<'_>, NodeId, Option<NodeId>) -> Option<NodeId> + 'a;

/// When the data are sums of eigenmodes with literal mode numbers, the
/// solution is the corresponding finite sum.
#[allow(clippy::too_many_arguments)]
fn finite_modes(
    cx: &mut Cx<'_>,
    p: &Problem,
    f: NodeId,
    g: Option<NodeId>,
    mode: NodeId,
    n: NodeId,
    first: i64,
    time_factor: &TimeFactor<'_>,
) -> Option<NodeId> {
    let _ = p;
    // Candidate modes: n = first .. first + 16; project the data on each by
    // matching terms structurally.
    let mut terms = Vec::new();
    let data = [Some(f), g];
    let mut remaining: Vec<Option<NodeId>> = data.iter().map(|d| d.map(|d| best(cx.graph, d).unwrap_or(d))).collect();
    for m in first..first + 17 {
        let number = cx.graph.int(m);
        let mode_m = cx.graph.substitute(mode, n, number);
        let mode_m = cx.simplify(mode_m);
        if cx.graph.number_of(mode_m).is_some_and(Number::is_zero) {
            continue;
        }
        let mut coefficients = [None, None];
        for (slot, value) in remaining.iter_mut().enumerate() {
            let Some(v) = *value else {
                continue;
            };
            // coefficient c with v = c * mode_m + rest, c free of the
            // variables: d(v)/d(mode) after freezing the mode.
            let c = derivative_by_node(cx, v, mode_m)?;
            if !p.constant(cx.graph, c) {
                return None;
            }
            if !cx.graph.number_of(c).is_some_and(Number::is_zero) {
                let cm = mul(cx.graph, &[c, mode_m]);
                let rest = sub(cx.graph, v, cm);
                *value = Some(cx.simplify(rest));
                coefficients[slot] = Some(c);
            }
        }
        if coefficients.iter().all(Option::is_none) {
            continue;
        }
        let zero = cx.graph.int(0);
        let a = coefficients[0].unwrap_or(zero);
        let factor = time_factor(cx, a, coefficients[1])?;
        let factor = cx.graph.substitute(factor, n, number);
        terms.push(mul(cx.graph, &[factor, mode_m]));
    }
    // Everything must have been accounted for.
    if remaining.iter().flatten().any(|&v| !cx.graph.number_of(v).is_some_and(Number::is_zero)) || terms.is_empty() {
        return None;
    }
    let sum = add(cx.graph, &terms);
    Some(cx.simplify(sum))
}

/// Laplace's equation on a rectangle `[0, a] × [0, b]` (or a box) with
/// homogeneous Dirichlet data on every face but the one at the far end of
/// the last variable, where `u = h`.
fn laplace_box(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    let d = p.dimension();
    if !(2..=3).contains(&d) || !p.homogeneous(cx.graph) {
        return None;
    }
    let zero_index = vec![0; d];
    if conditions.0.iter().any(|c| c.derivative != zero_index) {
        return None;
    }
    // For every variable, the two ends; the last variable carries the data.
    let mut lengths = Vec::with_capacity(d);
    let mut data = None;
    for k in 0..d {
        let on: Vec<&Condition> = conditions.0.iter().filter(|c| c.on == k).collect();
        let [first, second] = on.as_slice() else {
            return None;
        };
        let zero_end = |graph: &Graph, c: &Condition| graph.number_of(c.point).is_some_and(Number::is_zero);
        let (start, end) = if zero_end(cx.graph, first) { (*first, *second) } else { (*second, *first) };
        if !zero_end(cx.graph, start) || !cx.graph.number_of(start.value).is_some_and(Number::is_zero) {
            return None;
        }
        lengths.push(end.point);
        if k + 1 == d {
            data = Some(end.value);
        } else if !cx.graph.number_of(end.value).is_some_and(Number::is_zero) {
            return None;
        }
    }
    let h = data?;
    let (sin, sinh) = (cx.graph.ops().lookup("sin")?, cx.graph.ops().lookup("sinh")?);
    let pi = cx.graph.ops().lookup("pi")?;
    let pi = cx.graph.node(pi, &[]);
    let defint = cx.graph.ops().lookup("defint")?;
    let sum = cx.graph.ops().lookup("sum")?;
    let infinity = cx.graph.ops().lookup("oo")?;
    let oo = cx.graph.node(infinity, &[]);
    let one = cx.graph.int(1);
    let zero = cx.graph.int(0);
    let last = p.vars[d - 1];
    let mut indices = Vec::new();
    let mut modes = Vec::new();
    let mut k2 = Vec::new();
    for (j, stem) in ["n", "m"].iter().enumerate().take(d - 1) {
        let (n, symbol) = dummy(cx, p, stem);
        cx.graph.assume(symbol, Facts::INTEGER | Facts::POSITIVE);
        let lengths_inv = powi(cx.graph, lengths[j], -1);
        let k = mul(cx.graph, &[n, pi, lengths_inv]);
        let kx = mul(cx.graph, &[k, p.vars[j]]);
        modes.push(cx.graph.node(sin, &[kx]));
        k2.push(powi(cx.graph, k, 2));
        indices.push(n);
    }
    // The term for given mode numbers: (2/a)(2/b) ∫∫ h modes · sinh(κ y)/sinh(κ c) · modes.
    let term_for = |cx: &mut Cx<'_>, values: &[NodeId]| -> Option<NodeId> {
        let modes: Vec<NodeId> = modes.iter().map(|&m| substitute_all(cx.graph, m, &indices, values)).collect();
        let k2: Vec<NodeId> = k2.iter().map(|&k| substitute_all(cx.graph, k, &indices, values)).collect();
        let total = add(cx.graph, &k2);
        let half = cx.graph.num(Number::fraction(1, 2)?);
        let kappa = pow(cx.graph, total, half);
        let mut projection = mul(cx.graph, &[h]);
        for &m in &modes {
            projection = mul(cx.graph, &[projection, m]);
        }
        let mut factors = Vec::new();
        #[allow(clippy::needless_range_loop)] // index is used for more than one array / arithmetic; iterator form would not be clearer
        for j in 0..d - 1 {
            projection = cx.graph.node(defint, &[projection, p.vars[j], zero, lengths[j]]);
            let two = cx.graph.int(2);
            let lengths_inv = powi(cx.graph, lengths[j], -1);
            factors.push(mul(cx.graph, &[two, lengths_inv]));
        }
        let kc = mul(cx.graph, &[kappa, lengths[d - 1]]);
        let ky = mul(cx.graph, &[kappa, last]);
        let numerator = cx.graph.node(sinh, &[ky]);
        let denominator = cx.graph.node(sinh, &[kc]);
        factors.push(projection);
        factors.push(numerator);
        factors.push(powi(cx.graph, denominator, -1));
        factors.extend(modes);
        let term = mul(cx.graph, &factors);
        Some(cx.simplify(term))
    };
    let general = term_for(cx, &indices)?;
    if d == 2 {
        let exact = |cx: &mut Cx<'_>, m: i64| -> Option<NodeId> {
            let number = cx.graph.int(m);
            term_for(cx, &[number])
        };
        return assemble(cx, indices[0], 1, general, &exact);
    }
    // A double series: only when the generic coefficient does not vanish
    // (a vanishing one means the data are special modes, which this
    // general formula cannot represent).
    if cx.graph.number_of(general).is_some_and(Number::is_zero) {
        return None;
    }
    let mut series = general;
    for &n in indices.iter().rev() {
        series = cx.graph.node(sum, &[series, n, one, oo]);
    }
    Some(series)
}

/// `term` with each of `from` replaced by the matching `to`.
fn substitute_all(
    graph: &mut Graph,
    term: NodeId,
    from: &[NodeId],
    to: &[NodeId],
) -> NodeId {
    from.iter().zip(to).fold(term, |acc, (&f, &t)| graph.substitute(acc, f, t))
}

// ----------------------------------------------------------------------
// Green's functions
// ----------------------------------------------------------------------

/// `Δu + k² u = f` (Poisson for `k = 0`) in all of space, `d = 1, 2, 3`:
/// the convolution of the source with the free-space Green's function.
fn green(
    cx: &mut Cx<'_>,
    p: &Problem,
) -> Option<NodeId> {
    let d = p.dimension();
    // a (u_xx + u_yy + ...) + c u + source = 0, all a equal.
    let a = p.coefficient(cx.graph, &p.unit(0, 2));
    let mut allowed = vec![vec![0; d]];
    for k in 0..d {
        let index = p.unit(k, 2);
        let ak = p.coefficient(cx.graph, &index);
        let difference = sub(cx.graph, ak, a);
        if !cx.is_zero(difference) {
            return None;
        }
        allowed.push(index);
    }
    if !p.only(cx.graph, &allowed) || p.nonlinear || !p.constant(cx.graph, a) {
        return None;
    }
    let inverse_a = powi(cx.graph, a, -1);
    let c = p.coefficient(cx.graph, &vec![0; d]);
    let k2 = mul(cx.graph, &[c, inverse_a]);
    let k2 = cx.simplify(k2);
    if !p.constant(cx.graph, k2) {
        return None;
    }
    // Δu + k² u = f with f = -source/a.
    let minus_one = cx.graph.int(-1);
    let f = mul(cx.graph, &[minus_one, p.source, inverse_a]);
    let f = cx.simplify(f);
    let helmholtz = !cx.graph.number_of(k2).is_some_and(Number::is_zero);
    if !helmholtz && cx.graph.number_of(f).is_some_and(Number::is_zero) {
        return None;
    }
    let pi = cx.graph.ops().lookup("pi")?;
    let pi = cx.graph.node(pi, &[]);
    let (ln, exp) = (cx.graph.ops().lookup("ln")?, cx.graph.ops().lookup("exp")?);
    let i = cx.graph.ops().lookup("I")?;
    let i = cx.graph.node(i, &[]);
    let mut source = f;
    let mut squares = Vec::new();
    let mut dummies = Vec::new();
    for (j, stem) in ["s", "r", "q"].iter().enumerate().take(d) {
        let (s, _) = dummy(cx, p, stem);
        source = cx.graph.substitute(source, p.vars[j], s);
        let difference = sub(cx.graph, p.vars[j], s);
        squares.push(powi(cx.graph, difference, 2));
        dummies.push(s);
    }
    let r2 = add(cx.graph, &squares);
    let half = cx.graph.num(Number::fraction(1, 2)?);
    let r = pow(cx.graph, r2, half);
    let kernel = match (d, helmholtz) {
        // Δ G = δ
        | (1, false) => mul(cx.graph, &[half, r]),
        | (2, false) => {
            let log = cx.graph.node(ln, &[r]);
            let two = cx.graph.int(2);
            let two_pi = mul(cx.graph, &[two, pi]);
            let two_pi_inv = powi(cx.graph, two_pi, -1);
            mul(cx.graph, &[log, two_pi_inv])
        },
        | (3, false) => {
            let four = cx.graph.int(4);
            let four_pi_r = mul(cx.graph, &[four, pi, r]);
            let minus_one = cx.graph.int(-1);
            let four_pi_r_inv = powi(cx.graph, four_pi_r, -1);
            mul(cx.graph, &[minus_one, four_pi_r_inv])
        },
        // (Δ + k²) G = δ, outgoing
        | (1, true) => {
            let k = pow(cx.graph, k2, half);
            let ikr = mul(cx.graph, &[i, k, r]);
            let wave = cx.graph.node(exp, &[ikr]);
            let two = cx.graph.int(2);
            let two_ik = mul(cx.graph, &[two, i, k]);
            let two_ik_inv = powi(cx.graph, two_ik, -1);
            mul(cx.graph, &[wave, two_ik_inv])
        },
        | (3, true) => {
            let k = pow(cx.graph, k2, half);
            let ikr = mul(cx.graph, &[i, k, r]);
            let wave = cx.graph.node(exp, &[ikr]);
            let four = cx.graph.int(4);
            let four_pi_r = mul(cx.graph, &[four, pi, r]);
            let minus_one = cx.graph.int(-1);
            let four_pi_r_inv = powi(cx.graph, four_pi_r, -1);
            mul(cx.graph, &[minus_one, wave, four_pi_r_inv])
        },
        | _ => return None,
    };
    let defint = cx.graph.ops().lookup("defint")?;
    let infinity = cx.graph.ops().lookup("oo")?;
    let oo = cx.graph.node(infinity, &[]);
    let minus_oo = neg(cx.graph, oo);
    let mut body = mul(cx.graph, &[kernel, source]);
    for &s in dummies.iter().rev() {
        body = cx.graph.node(defint, &[body, s, minus_oo, oo]);
    }
    if helmholtz && cx.graph.number_of(f).is_some_and(Number::is_zero) {
        // Homogeneous Helmholtz: the separated plane waves.
        let mut phase = Vec::new();
        let mut directions = Vec::new();
        for (j, stem) in ["kx", "ky", "kz"].iter().enumerate().take(d) {
            let (kj, _) = dummy(cx, p, stem);
            phase.push(mul(cx.graph, &[kj, p.vars[j]]));
            directions.push(kj);
        }
        let _ = directions;
        let total = add(cx.graph, &phase);
        let exponent = mul(cx.graph, &[i, total]);
        let wave = cx.graph.node(exp, &[exponent]);
        return Some(wave);
    }
    Some(body)
}

// ----------------------------------------------------------------------
// Checks
// ----------------------------------------------------------------------

/// Whether `solution` satisfies the equation: the residual with `u`
/// replaced must simplify to zero, or vanish numerically at sample
/// points. Solutions with arbitrary functions are accepted unchecked.
fn verified(
    cx: &mut Cx<'_>,
    p: &Problem,
    solution: NodeId,
) -> bool {
    // Arbitrary functions (F, G, A, ...) cannot be evaluated.
    let mut stack = vec![solution];
    while let Some(n) = stack.pop() {
        if cx.graph.op(n) == core::APPLY {
            return true;
        }
        stack.extend_from_slice(cx.graph.children(n));
    }
    let substituted = cx.graph.replace_subterm(p.residual, p.unknown, solution);
    let residual = cx.simplify(substituted);
    if cx.graph.number_of(residual).is_some_and(Number::is_zero) {
        return true;
    }
    let symbols = cx.graph.free_symbols(cx.graph.find(residual)).to_vec();
    let mut evaluated = 0;
    for sample in 0..4_u32 {
        let mut env = Env::numeric(0.0);
        for &s in &symbols {
            env.bind(s, 0.3 + 0.23 * f64::from(sample) + 0.11 * f64::from(s.raw() % 7));
        }
        match cx.graph.eval(residual, &env) {
            | Some(v) if v.is_finite() => {
                if v.abs() > 1e-7 {
                    return false;
                }
                evaluated += 1;
            },
            | _ => {},
        }
    }
    evaluated > 0
}

/// Whether `solution` meets every condition.
fn satisfies(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
    solution: NodeId,
) -> bool {
    for c in &conditions.0 {
        let mut value = solution;
        for (k, &count) in c.derivative.iter().enumerate() {
            for _ in 0..count {
                match derivative(cx.graph, value, p.vars[k]) {
                    | Some(d) => value = d,
                    | None => return false,
                }
            }
        }
        let on = cx.graph.substitute(value, p.vars[c.on], c.point);
        let difference = sub(cx.graph, on, c.value);
        if !cx.is_zero(difference) {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rules::testing::numeric;
    use crate::rules::testing::reduce_with;

    const ASSUME: [(&str, Facts); 4] =
        [("c", Facts::POSITIVE), ("k", Facts::POSITIVE), ("L", Facts::POSITIVE), ("a", Facts::POSITIVE)];

    fn run(src: &str) -> String {
        let (text, reduced) = reduce_with(&[pde()], src, &ASSUME);
        assert!(reduced, "`{src}` was not fully reduced: {text}");
        text
    }

    fn any(src: &str) -> String {
        reduce_with(&[pde()], src, &ASSUME).0
    }

    #[test]
    fn classification() {
        assert_eq!(
            run("pde_classify(diff(u(x, t), t) - k*diff(diff(u(x, t), x), x), u(x, t))"),
            "list(heat, 2, 2, true, true, parabolic, list(separation_of_variables, fourier_transform, heat_kernel))"
        );
        assert!(run("pde_classify(diff(diff(u(x, t), t), t) - c^2*diff(diff(u(x, t), x), x), u(x, t))").starts_with("list(wave, 2, 2, true, true, hyperbolic"));
        assert!(run("pde_classify(diff(diff(u(x, y), x), x) + diff(diff(u(x, y), y), y), u(x, y))").starts_with("list(laplace, 2, 2, true, true, elliptic"));
        assert!(run("pde_classify(diff(diff(u(x, y), x), x) + diff(diff(u(x, y), y), y) - x*y, u(x, y))").starts_with("list(poisson, 2, 2, true, false, elliptic"));
        assert!(run("pde_classify(diff(u(x, t), t) + u(x, t)*diff(u(x, t), x), u(x, t))").starts_with("list(burgers, 1, 2, false"));
        assert!(run("pde_classify(diff(u(x, y), x) + 2*diff(u(x, y), y), u(x, y))").starts_with("list(transport, 1, 2, true, true, first_order"));
        assert!(run("pde_classify(I*diff(u(x, t), t) + diff(diff(u(x, t), x), x), u(x, t))").starts_with("list(schrodinger"));
        // Mixed derivatives: u_xx + 3 u_xy + u_yy has discriminant 5.
        assert!(run("pde_classify(diff(diff(u(x, y), x), x) + 3*diff(diff(u(x, y), x), y) + diff(diff(u(x, y), y), y), u(x, y))").contains("hyperbolic"));
        assert_eq!(run("pde_order(diff(diff(diff(u(x, t), x), x), x) + diff(u(x, t), t), u(x, t))"), "3");
    }

    #[test]
    fn first_order_equations() {
        // u_x + 2 u_y = 0: u = F(2x - y).
        assert_eq!(any("pdsolve(diff(u(x, y), x) + 2*diff(u(x, y), y) = 0, u(x, y))"), "u(x, y) = F(2*x - y)");
        // With a source and a decay term, checked by substitution.
        let solution = any("pdsolve(diff(u(x, y), x) + diff(u(x, y), y) + u(x, y) = x, u(x, y))");
        assert!(solution.contains("F(x - y)"), "{solution}");
        // Initial data along y = 0.
        assert_eq!(
            run("solve_pde_by_characteristics(diff(u(x, y), y) + 3*diff(u(x, y), x) = 0, u(x, y), list(u(x, 0) = sin(x)))"),
            "u(x, y) = sin(x - 3*y)"
        );
        // Variable coefficients: x u_x + y u_y = 0 is constant along rays.
        let rays = any("pdsolve(x*diff(u(x, y), x) + y*diff(u(x, y), y) = 0, u(x, y))");
        assert!(rays.contains("F(") && rays.contains("y/x"), "{rays}");
        // Burgers' equation: the implicit solution.
        assert_eq!(
            run("solve_burgers_equation(diff(u(x, t), t) + u(x, t)*diff(u(x, t), x) = 0, u(x, t), list(u(x, 0) = x))"),
            "u(x, t) = x - t*u(x, t)"
        );
    }

    #[test]
    fn wave_equation_on_the_line() {
        assert_eq!(any("pdsolve(diff(diff(u(x, t), t), t) = c^2*diff(diff(u(x, t), x), x), u(x, t))"), any("u(x, t) = F(x - c*t) + G(x + c*t)"));
        assert_eq!(
            run("solve_wave_equation_1d_dalembert(diff(diff(u(x, t), t), t) = 4*diff(diff(u(x, t), x), x), u(x, t), list(u(x, 0) = x^2, at(diff(u(x, t), t), t, 0) = 0))"),
            "u(x, t) = 4*t^2 + x^2"
        );
        assert_eq!(
            run("pdsolve(diff(diff(u(x, t), t), t) = diff(diff(u(x, t), x), x), u(x, t), list(u(x, 0) = 0, at(diff(u(x, t), t), t, 0) = cos(x)))"),
            "u(x, t) = 1/2*(sin(t + x) - sin(x - t))"
        );
    }

    #[test]
    fn heat_equation() {
        // A single mode on [0, pi] decays with rate k.
        assert_eq!(
            run("pdsolve(diff(u(x, t), t) = k*diff(diff(u(x, t), x), x), u(x, t), list(u(0, t) = 0, u(pi, t) = 0, u(x, 0) = 3*sin(2*x)))"),
            "u(x, t) = 3*exp(-4*k*t)*sin(2*x)"
        );
        // Neumann ends: cos modes, the mean is conserved.
        assert_eq!(
            run("solve_heat_equation_1d(diff(u(x, t), t) = diff(diff(u(x, t), x), x), u(x, t), list(at(diff(u(x, t), x), x, 0) = 0, at(diff(u(x, t), x), x, pi) = 0, u(x, 0) = 1 + cos(x)))"),
            any("u(x, t) = 1 + exp(-t)*cos(x)")
        );
        // General data: a Fourier series with closed-form coefficients.
        let series = any("pdsolve(diff(u(x, t), t) = diff(diff(u(x, t), x), x), u(x, t), list(u(0, t) = 0, u(1, t) = 0, u(x, 0) = x*(1 - x)))");
        assert!(series.starts_with("u(x, t) = sum(") && series.contains("exp(") && !series.contains("defint"), "{series}");
        // On the line: a Gaussian stays Gaussian.
        let line = run("solve_heat_equation_1d(diff(u(x, t), t) = diff(diff(u(x, t), x), x), u(x, t), list(u(x, 0) = exp(-x^2)))");
        let (v, _) = numeric(&[pde()], &format!("at(at({}, x, 0.5), t, 0.25)", line.trim_start_matches("u(x, t) = ")), &[], 1e-10);
        let expected = (-(0.25_f64) / 2.0).exp() / 2.0_f64.sqrt();
        assert!((v - expected).abs() < 1e-9, "{line}: {v} vs {expected}");
    }

    #[test]
    fn laplace_and_poisson() {
        // Data sin(pi x) on the top of the unit square.
        assert_eq!(
            run("pdsolve(diff(diff(u(x, y), x), x) + diff(diff(u(x, y), y), y) = 0, u(x, y), list(u(0, y) = 0, u(1, y) = 0, u(x, 0) = 0, u(x, 1) = sin(pi*x)))")
                .replace(' ', ""),
            any("u(x, y) = sinh(pi*y)*sin(pi*x)/sinh(pi)").replace(' ', "")
        );
        // Poisson in 3D: the Newtonian potential.
        let potential = any("solve_poisson_equation_3d(diff(diff(u(x, y, z), x), x) + diff(diff(u(x, y, z), y), y) + diff(diff(u(x, y, z), z), z) = f(x, y, z), u(x, y, z))");
        assert!(potential.contains("defint(") && potential.contains("pi"), "{potential}");
        let two = any("solve_poisson_equation_2d(diff(diff(u(x, y), x), x) + diff(diff(u(x, y), y), y) = rho(x, y), u(x, y))");
        assert!(two.contains("ln("), "{two}");
    }

    #[test]
    fn evolution_in_three_dimensions() {
        // Kirchhoff: data f = x^2 + y^2 + z^2, g = 0 gives r² + 3c²t².
        assert_eq!(
            run("solve_wave_equation_3d(diff(diff(u(x, y, z, t), t), t) = c^2*(diff(diff(u(x, y, z, t), x), x) + diff(diff(u(x, y, z, t), y), y) + diff(diff(u(x, y, z, t), z), z)), u(x, y, z, t), list(u(x, y, z, 0) = x^2 + y^2 + z^2, at(diff(u(x, y, z, t), t), t, 0) = 0))"),
            run("u(x, y, z, t) = x^2 + y^2 + z^2 + 3*c^2*t^2")
        );
        let heat = any("solve_heat_equation_3d(diff(u(x, y, z, t), t) = diff(diff(u(x, y, z, t), x), x) + diff(diff(u(x, y, z, t), y), y) + diff(diff(u(x, y, z, t), z), z), u(x, y, z, t), list(u(x, y, z, 0) = f(x, y, z)))");
        assert!(heat.matches("defint(").count() == 3, "{heat}");
    }

    #[test]
    fn quantum_and_relativistic_waves() {
        let free = any("solve_schrodinger_equation(I*diff(u(x, t), t) = -diff(diff(u(x, t), x), x), u(x, t))");
        assert!(free.contains("A(k)") && free.contains("exp("), "{free}");
        let kg = any("solve_klein_gordon_equation(diff(diff(u(x, t), t), t) - diff(diff(u(x, t), x), x) + a^2*u(x, t) = 0, u(x, t))");
        assert!(kg.contains("(a^2 + k^2)^(1/2)") || kg.contains("(k^2 + a^2)^(1/2)"), "{kg}");
        let helmholtz = any("solve_helmholtz_equation(diff(diff(u(x, y, z), x), x) + diff(diff(u(x, y, z), y), y) + diff(diff(u(x, y, z), z), z) + 4*u(x, y, z) = f(x, y, z), u(x, y, z))");
        assert!(helmholtz.contains("exp(") && helmholtz.contains("defint("), "{helmholtz}");
    }
}
