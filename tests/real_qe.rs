//! Integration tests for real quantifier elimination (`real_qe`,
//! `cad_cells`) and for the routing of the logic rule set's `exists` /
//! `forall` to it.
//!
//! The tests do three independent things:
//!
//! * they evaluate every returned quantifier-free formula at many rational
//!   parameter values and compare it with an independent one-variable
//!   decision procedure built from the crate's existing exact
//!   `real_root_intervals` and `count_real_roots` operators;
//! * they check that `cad_cells` reports sign-invariant cells by evaluating
//!   the input polynomial at the reported sample and at perturbed samples of
//!   open cells;
//! * they check that non-polynomial requests stay unevaluated instead of
//!   falling back to Boolean semantics.

use std::time::Instant;

use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::NodeId;
use rssn::graph::Number;
use rssn::graph::OpId;
use rssn::graph::Payload;
use rssn::graph::Saturate;
use rssn::graph::op::core;
use rssn::rules::standard;

/// Truth of a concrete formula: comparisons are evaluated from their two
/// arithmetic sides, which is needed because `=` has no scalar evaluator of
/// its own.
fn formula_truth(
    graph: &Graph,
    node: NodeId,
    env: &Env,
    depth: usize,
) -> Option<bool> {
    if depth > 256 {
        return None;
    }
    if let Some(Payload::Bool(value)) = graph.payload(node) {
        return Some(*value);
    }
    let op = graph.op(node);
    let children = graph.children(node);
    let name = graph.ops().get(op).name.as_ref();
    if name == "true" {
        return Some(true);
    }
    if name == "false" {
        return Some(false);
    }
    if name == "and" {
        let mut result = true;
        for &child in children {
            result &= formula_truth(graph, child, env, depth.saturating_add(1))?;
        }
        return Some(result);
    }
    if name == "or" {
        let mut result = false;
        for &child in children {
            result |= formula_truth(graph, child, env, depth.saturating_add(1))?;
        }
        return Some(result);
    }
    if name == "not" {
        return Some(!formula_truth(
            graph,
            *children.first()?,
            env,
            depth.saturating_add(1),
        )?);
    }
    if name == "implies" {
        let &[a, b] = children else { return None };
        return Some(
            !formula_truth(graph, a, env, depth.saturating_add(1))?
                || formula_truth(graph, b, env, depth.saturating_add(1))?,
        );
    }
    if name == "iff" {
        let &[a, b] = children else { return None };
        return Some(
            formula_truth(graph, a, env, depth.saturating_add(1))?
                == formula_truth(graph, b, env, depth.saturating_add(1))?,
        );
    }
    if children.len() == 2 && (op == core::EQ || matches!(name, "lt" | "le" | "gt" | "ge" | "ne")) {
        let a = graph.eval(*children.first()?, env)?;
        let b = graph.eval(*children.get(1)?, env)?;
        let tolerance = 1e-9 * (1.0 + a.abs() + b.abs());
        return Some(match name {
            | "lt" => a < b - tolerance,
            | "le" => a <= b + tolerance,
            | "gt" => a > b + tolerance,
            | "ge" => a >= b - tolerance,
            | "ne" => (a - b).abs() > tolerance,
            | _ => (a - b).abs() <= tolerance,
        });
    }
    None
}

/// A graph with every standard rule set installed, reused across requests.
struct Harness {
    graph: Graph,
    engine: Engine,
}

impl Harness {
    fn new() -> Self {
        let mut graph = Graph::new();
        let engine = Engine::install(&mut graph, &standard()).expect("standard rule sets install");
        Self { graph, engine }
    }

    /// Reduces `src`, returning the displayed answer and whether it is free
    /// of heavy operators.
    fn reduce(
        &mut self,
        src: &str,
    ) -> (String, bool) {
        let root = self
            .graph
            .parse(src)
            .unwrap_or_else(|e| panic!("cannot parse `{src}`: {e}"));
        self.engine.run(
            &mut self.graph,
            &[root],
            &Env::symbolic(),
            &Saturate,
            &Budget::default(),
        );
        assert!(
            self.graph.conflicts().is_empty(),
            "unsound reduction of `{src}`"
        );
        match Extractor::new(&self.graph, &[root], &ClosedForm).build(&mut self.graph, root) {
            | Some(node) => (self.graph.display(node), true),
            | None => (src.to_owned(), false),
        }
    }

    fn answer(
        &mut self,
        src: &str,
    ) -> String {
        let (text, reduced) = self.reduce(src);
        assert!(reduced, "`{src}` was not fully reduced: {text}");
        text
    }

    /// Evaluates a concrete formula under float bindings.
    fn evaluate(
        &mut self,
        src: &str,
        bindings: &[(&str, f64)],
    ) -> Option<bool> {
        let root = self.graph.parse(src).ok()?;
        if let Some(Payload::Bool(value)) = self.graph.payload(root) {
            return Some(*value);
        }
        let mut env = Env::numeric(0.0);
        for (name, value) in bindings {
            env.bind(self.graph.interner_mut().symbol(name), *value);
        }
        formula_truth(&self.graph, root, &env, 0)
    }

    /// Runs one request and extracts its answer node.
    fn request(
        &mut self,
        node: NodeId,
    ) -> Option<NodeId> {
        self.engine.run(
            &mut self.graph,
            &[node],
            &Env::symbolic(),
            &Saturate,
            &Budget::default(),
        );
        Extractor::new(&self.graph, &[node], &ClosedForm).build(&mut self.graph, node)
    }
}

fn harness() -> Harness {
    Harness::new()
}

/// The fully reduced answer to `src`.
fn answer(src: &str) -> String {
    harness().answer(src)
}

// ---------------------------------------------------------------------------
// The independent one-variable decision procedure.

/// Operators used by the oracle.
struct Ops {
    and: OpId,
    or: OpId,
    not: OpId,
    implies: OpId,
    iff: OpId,
    lt: OpId,
    le: OpId,
    gt: OpId,
    ge: OpId,
    ne: OpId,
    truth: OpId,
    falsity: OpId,
    root_intervals: OpId,
    count_roots: OpId,
}

impl Ops {
    fn lookup(graph: &Graph) -> Option<Self> {
        let of = |name: &str| graph.ops().lookup(name);
        Some(Self {
            and: of("and")?,
            or: of("or")?,
            not: of("not")?,
            implies: of("implies")?,
            iff: of("iff")?,
            lt: of("lt")?,
            le: of("le")?,
            gt: of("gt")?,
            ge: of("ge")?,
            ne: of("ne")?,
            truth: of("true")?,
            falsity: of("false")?,
            root_intervals: of("real_root_intervals")?,
            count_roots: of("count_real_roots")?,
        })
    }

    fn is_comparison(
        &self,
        op: OpId,
    ) -> bool {
        op == self.lt
            || op == self.le
            || op == self.gt
            || op == self.ge
            || op == self.ne
            || op == core::EQ
    }
}

/// An exact rational number for a float that is a small dyadic rational.
fn exact(value: f64) -> Option<Number> {
    let scaled = (value * 1024.0).round();
    if (scaled / 1024.0 - value).abs() > 0.0 {
        return None;
    }
    Number::fraction(scaled as i64, 1024)
}

/// Collects the comparison nodes of a Boolean formula.
fn comparisons(
    graph: &Graph,
    ops: &Ops,
    node: NodeId,
    out: &mut Vec<NodeId>,
) {
    let op = graph.op(node);
    if ops.is_comparison(op) {
        out.push(node);
        return;
    }
    for &child in graph.children(node) {
        comparisons(graph, ops, child, out);
    }
}

/// The term `lhs - rhs` of a comparison.
fn difference(
    graph: &mut Graph,
    node: NodeId,
) -> Option<NodeId> {
    let &[lhs, rhs] = graph.children(node) else {
        return None;
    };
    let minus_one = graph.int(-1);
    let negative = graph.node(core::MUL, &[minus_one, rhs]);
    Some(graph.node(core::ADD, &[lhs, negative]))
}

/// An isolating interval, with exact endpoint nodes and float values.
#[derive(Copy, Clone)]
struct Interval {
    low: NodeId,
    high: NodeId,
    low_value: f64,
    high_value: f64,
}

/// One sample of the sign-invariant decomposition of the `y` line.
enum Sample {
    Rational(f64),
    Root(Interval),
}

/// The truth of a Boolean combination at a sign oracle.
fn truth(
    graph: &Graph,
    ops: &Ops,
    node: NodeId,
    sign: &dyn Fn(NodeId) -> i32,
) -> Option<bool> {
    let op = graph.op(node);
    if op == ops.truth || matches!(graph.payload(node), Some(Payload::Bool(true))) {
        return Some(true);
    }
    if op == ops.falsity || matches!(graph.payload(node), Some(Payload::Bool(false))) {
        return Some(false);
    }
    if ops.is_comparison(op) {
        let sign = sign(node);
        return Some(if op == core::EQ {
            sign == 0
        } else if op == ops.ne {
            sign != 0
        } else if op == ops.lt {
            sign < 0
        } else if op == ops.le {
            sign <= 0
        } else if op == ops.gt {
            sign > 0
        } else {
            sign >= 0
        });
    }
    let kids: Vec<NodeId> = graph.children(node).to_vec();
    if op == ops.not {
        return Some(!truth(graph, ops, *kids.first()?, sign)?);
    }
    if op == ops.and {
        let mut result = true;
        for child in kids {
            result &= truth(graph, ops, child, sign)?;
        }
        return Some(result);
    }
    if op == ops.or {
        let mut result = false;
        for child in kids {
            result |= truth(graph, ops, child, sign)?;
        }
        return Some(result);
    }
    if let &[a, b] = kids.as_slice() {
        if op == ops.implies {
            return Some(!truth(graph, ops, a, sign)? || truth(graph, ops, b, sign)?);
        }
        if op == ops.iff {
            return Some(truth(graph, ops, a, sign)? == truth(graph, ops, b, sign)?);
        }
    }
    None
}

/// Sign of the polynomial of `atom` at the sample.
fn atom_sign(
    harness: &mut Harness,
    ops: &Ops,
    atom: NodeId,
    y: NodeId,
    sample: &Sample,
) -> Option<i32> {
    let expression = difference(&mut harness.graph, atom)?;
    match sample {
        | Sample::Rational(value) => {
            let mut env = Env::numeric(0.0);
            env.bind(harness.graph.as_symbol(y)?, *value);
            let residual = harness.graph.eval(expression, &env)?;
            Some(if residual.abs() < 1e-9 {
                0
            } else if residual < 0.0 {
                -1
            } else {
                1
            })
        },
        | Sample::Root(interval) => {
            let request = harness.graph.node(
                ops.count_roots,
                &[expression, y, interval.low, interval.high],
            );
            let result = harness.request(request)?;
            let count = harness.graph.as_number(result).and_then(Number::to_i64)?;
            if count > 0 {
                return Some(0);
            }
            let middle = 0.5 * (interval.low_value + interval.high_value);
            let mut env = Env::numeric(0.0);
            env.bind(harness.graph.as_symbol(y)?, middle);
            let residual = harness.graph.eval(expression, &env)?;
            Some(if residual.abs() < 1e-9 {
                0
            } else if residual < 0.0 {
                -1
            } else {
                1
            })
        },
    }
}

impl Harness {
    /// Decides `Q y . body` at a rational assignment of the parameters using
    /// exact Sturm isolation and root counting for the boundaries.
    fn decide(
        &mut self,
        body: &str,
        y_name: &str,
        parameters: &[(&str, f64)],
        universal: bool,
    ) -> Option<bool> {
        let ops = Ops::lookup(&self.graph)?;
        let mut body = self.graph.parse(body).ok()?;
        for (name, value) in parameters {
            let symbol = self.graph.interner_mut().symbol(name);
            let leaf = self.graph.symbol_node(symbol);
            let replacement = self.graph.num(exact(*value)?);
            body = self.graph.substitute(body, leaf, replacement);
        }
        let y_symbol = self.graph.interner_mut().symbol(y_name);
        let y = self.graph.symbol_node(y_symbol);
        let mut atoms = Vec::new();
        comparisons(&self.graph, &ops, body, &mut atoms);
        if atoms.is_empty() {
            return truth(&self.graph, &ops, body, &|_| 0);
        }
        let mut product = self.graph.int(1);
        for &atom in &atoms {
            let factor = difference(&mut self.graph, atom)?;
            product = self.graph.node(core::MUL, &[product, factor]);
        }
        let request = self.graph.node(ops.root_intervals, &[product, y]);
        let result = self.request(request)?;
        if self.graph.op(result) != core::LIST {
            return None;
        }
        let mut intervals: Vec<Interval> = Vec::new();
        for &item in self.graph.children(result) {
            if self.graph.op(item) != core::LIST {
                continue;
            }
            let pair = self.graph.children(item).to_vec();
            let &[low, high] = pair.as_slice() else {
                continue;
            };
            intervals.push(Interval {
                low,
                high,
                low_value: self.graph.as_number(low).map(Number::to_f64)?,
                high_value: self.graph.as_number(high).map(Number::to_f64)?,
            });
        }
        let mut samples: Vec<Sample> = Vec::new();
        let initial = match intervals.first() {
            | Some(first) => first.low_value - 1.0,
            | None => 0.0,
        };
        samples.push(Sample::Rational(initial));
        for (index, interval) in intervals.iter().enumerate() {
            if (interval.high_value - interval.low_value).abs() < 1e-12 {
                samples.push(Sample::Rational(interval.low_value));
            } else {
                samples.push(Sample::Root(*interval));
            }
            let next = match intervals.get(index + 1) {
                | Some(next) => 0.5 * (interval.high_value + next.low_value),
                | None => interval.high_value + 1.0,
            };
            samples.push(Sample::Rational(next));
        }
        let mut result = universal;
        for sample in &samples {
            let mut signs = Vec::with_capacity(atoms.len());
            for &atom in &atoms {
                signs.push(atom_sign(self, &ops, atom, y, &sample)?);
            }
            let sign = |atom: NodeId| -> i32 {
                atoms
                    .iter()
                    .position(|&a| a == atom)
                    .and_then(|position| signs.get(position).copied())
                    .unwrap_or(0)
            };
            let value = truth(&self.graph, &ops, body, &sign)?;
            result = if universal {
                result && value
            } else {
                result || value
            };
        }
        Some(result)
    }
}

// ---------------------------------------------------------------------------
// One-variable elimination.

#[test]
fn one_variable_examples() {
    assert_eq!(answer("real_qe(exists(x, x^2 - 2 = 0), list(x))"), "true");
    assert_eq!(
        answer("real_qe(forall(x, gt(x^2 + 1, 0)), list(x))"),
        "true"
    );
    assert_eq!(answer("real_qe(exists(x, lt(x^2, 0)), list(x))"), "false");
    assert_eq!(answer("real_qe(forall(x, ge(x^2, 0)), list(x))"), "true");
    // The listed variable may also be free; it is then existential.
    assert_eq!(answer("real_qe(x^2 - 2 = 0, list(x))"), "true");
}

#[test]
fn parameter_quantifiers_with_dense_samples() {
    let mut harness = harness();
    let cases = [
        ("lt(x^2 + y^2, 1)", vec!["x"], false),
        ("ge(y^2, x)", vec!["x"], true),
        ("y^2 = x", vec!["x"], false),
        ("x*y = 1", vec!["x"], false),
        ("and(y^2 = x, gt(y, 0))", vec!["x"], false),
        ("a*y = 1", vec!["a"], false),
    ];
    let parameter_values = [
        -2.0, -1.5, -1.0, -0.75, -0.5, -0.25, -0.125, 0.0, 0.125, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0,
        3.0,
    ];
    for (body, parameters, universal) in cases {
        let quantifier = if universal {
            "forall"
        } else {
            "exists"
        };
        let source = format!("real_qe({quantifier}(y, {body}), list(y))");
        let text = harness.answer(&source);
        for &value in &parameter_values {
            let bindings: Vec<(&str, f64)> = parameters.iter().map(|&p| (p, value)).collect();
            let expected = harness
                .decide(body, "y", &bindings, universal)
                .unwrap_or_else(|| panic!("cannot decide `{body}` at {value}"));
            let actual = harness
                .evaluate(&text, &bindings)
                .unwrap_or_else(|| panic!("cannot evaluate `{text}`"));
            assert_eq!(
                actual, expected,
                "`{body}` ({quantifier}) at {value}: `{text}`"
            );
        }
    }
}

#[test]
fn quadratic_case_split_against_the_oracle() {
    let mut harness = harness();
    let exists = harness.answer("real_qe(exists(x, a*x^2 + b*x + c = 0), list(x))");
    let universal = harness.answer("real_qe(forall(x, gt(a*x^2 + b*x + c, 0)), list(x))");
    let values = [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0];
    for &a in &values {
        for &b in &values {
            for &c in &values {
                let bindings = [("a", a), ("b", b), ("c", c)];
                let expected_exists = harness
                    .decide("a*x^2 + b*x + c = 0", "x", &bindings, false)
                    .unwrap_or(false);
                assert_eq!(
                    harness.evaluate(&exists, &bindings).unwrap_or(false),
                    expected_exists,
                    "exists at a={a} b={b} c={c}: `{exists}`"
                );
                let expected_universal = harness
                    .decide("gt(a*x^2 + b*x + c, 0)", "x", &bindings, true)
                    .unwrap_or(false);
                assert_eq!(
                    harness.evaluate(&universal, &bindings).unwrap_or(false),
                    expected_universal,
                    "forall at a={a} b={b} c={c}: `{universal}`"
                );
            }
        }
    }
}

#[test]
fn sphere_and_cylinder() {
    let mut harness = harness();
    let sphere = harness.answer("real_qe(exists(z, lt(x^2 + y^2 + z^2, 1)), list(z))");
    for &x in &[-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5] {
        for &y in &[-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5] {
            let expected = harness
                .decide("lt(x^2 + y^2 + z^2, 1)", "z", &[("x", x), ("y", y)], false)
                .unwrap_or(false);
            assert_eq!(
                harness
                    .evaluate(&sphere, &[("x", x), ("y", y)])
                    .unwrap_or(false),
                expected,
                "sphere at x={x} y={y}: `{sphere}`"
            );
        }
    }
    let cylinder = harness.answer("real_qe(exists(z, lt(x^2 + y^2, 1)), list(z))");
    assert!(
        harness
            .evaluate(&cylinder, &[("x", 0.5), ("y", 0.5)])
            .unwrap_or(false)
    );
    assert!(
        !harness
            .evaluate(&cylinder, &[("x", 1.5), ("y", 0.0)])
            .unwrap_or(false)
    );
}

#[test]
fn nested_quantifiers() {
    let mut harness = harness();
    assert_eq!(
        harness.answer("real_qe(forall(a, exists(x, gt(x^2 + a, 0))), list(a))"),
        "true"
    );
    assert_eq!(
        harness.answer("real_qe(exists(a, forall(x, ge(x^2, a))), list(a))"),
        "true"
    );
    // forall a exists x: a*x^2 >= b is equivalent to b <= 0.
    let bounded = harness.answer("real_qe(forall(a, exists(x, ge(a*x^2, b))), list(a))");
    assert!(harness.evaluate(&bounded, &[("b", -2.0)]).unwrap_or(false));
    assert!(harness.evaluate(&bounded, &[("b", 0.0)]).unwrap_or(false));
    assert!(!harness.evaluate(&bounded, &[("b", 1.0)]).unwrap_or(false));
}

#[test]
fn logic_exists_and_forall_are_routed() {
    let mut harness = harness();
    // Real polynomial formulas are decided by real QE.
    let ball = harness.answer("exists(y, lt(x^2 + y^2, 1))");
    assert!(harness.evaluate(&ball, &[("x", 0.0)]).unwrap_or(false));
    assert!(harness.evaluate(&ball, &[("x", -0.5)]).unwrap_or(false));
    assert!(!harness.evaluate(&ball, &[("x", 1.5)]).unwrap_or(false));
    let square = harness.answer("forall(y, ge(y^2, x))");
    assert!(harness.evaluate(&square, &[("x", -1.0)]).unwrap_or(false));
    assert!(!harness.evaluate(&square, &[("x", 1.0)]).unwrap_or(false));
    // Genuine Boolean atoms keep the Boolean semantics.
    assert_eq!(harness.answer("exists(y, and(x, y))"), "x");
    assert_eq!(harness.answer("forall(y, or(x, y))"), "x");
    // A non-polynomial real formula must stay unevaluated.
    let (text, reduced) = harness.reduce("exists(y, sin(y) = 0)");
    assert!(!reduced, "`{text}` must stay unevaluated");
}

// ---------------------------------------------------------------------------
// CAD cells.

/// One reported cell.
struct CadCell {
    sample: Vec<f64>,
    open: Vec<bool>,
    signs: Vec<i32>,
}

fn sign_of(value: f64) -> i32 {
    if value.abs() < 1e-9 {
        0
    } else if value < 0.0 {
        -1
    } else {
        1
    }
}

impl Harness {
    /// Reads the `cad_cells` result.
    fn cad_cells(
        &mut self,
        src: &str,
    ) -> Vec<CadCell> {
        let root = self.graph.parse(src).expect("parse");
        let node = self.request(root).expect("reduced");
        assert_eq!(self.graph.op(node), core::LIST);
        let mut out = Vec::new();
        for &cell in self.graph.children(node).to_vec().iter() {
            assert_eq!(self.graph.op(cell), core::LIST);
            let parts = self.graph.children(cell).to_vec();
            assert_eq!(parts.len(), 3);
            let sample: Vec<f64> = self
                .graph
                .children(parts[0])
                .iter()
                .map(|&c| {
                    self.graph
                        .as_number(c)
                        .map(Number::to_f64)
                        .or_else(|| self.graph.eval(c, &Env::numeric(0.0)))
                        .expect("numeric sample")
                })
                .collect();
            let open: Vec<bool> = self
                .graph
                .children(parts[1])
                .iter()
                .map(|&c| {
                    match self.graph.payload(c) {
                        | Some(Payload::Bool(value)) => *value,
                        | _ => panic!("open flag is not a boolean"),
                    }
                })
                .collect();
            let signs: Vec<i32> = self
                .graph
                .children(parts[2])
                .iter()
                .map(|&c| {
                    let value = self
                        .graph
                        .as_number(c)
                        .and_then(Number::to_i64)
                        .expect("sign");
                    i32::try_from(value).expect("small sign")
                })
                .collect();
            out.push(CadCell { sample, open, signs });
        }
        out
    }
}

#[test]
fn cad_cells_are_sign_invariant() {
    let mut harness = harness();
    let cells = harness.cad_cells("cad_cells(list(x^2 + y^2 - 1), list(x, y))");
    assert!(cells.len() > 9, "{}", cells.len());
    assert!(cells.iter().any(|c| c.open.iter().all(|&o| o)));
    assert!(cells.iter().any(|c| c.open.iter().any(|&o| !o)));
    let mut sections = 0;
    for cell in &cells {
        assert_eq!(cell.sample.len(), 2);
        assert_eq!(cell.signs.len(), 1);
        let (x, y) = (cell.sample[0], cell.sample[1]);
        let value = x * x + y * y - 1.0;
        assert_eq!(cell.signs[0], sign_of(value), "sample sign at ({x}, {y})");
        if cell.open.iter().all(|&o| o) {
            for (dx, dy) in [
                (1e-4, 0.0),
                (-1e-4, 0.0),
                (0.0, 1e-4),
                (0.0, -1e-4),
                (7e-5, -3e-5),
            ] {
                let (px, py) = (x + dx, y + dy);
                let perturbed = px * px + py * py - 1.0;
                assert_eq!(
                    sign_of(perturbed),
                    cell.signs[0],
                    "perturbation ({dx}, {dy}) of ({x}, {y})"
                );
            }
        } else {
            sections += 1;
        }
    }
    assert!(sections > 0);
    assert!(
        cells.iter().any(|c| c.signs[0] == 0),
        "the decomposition must contain cells on the circle"
    );
}

#[test]
fn cad_cells_reports_openness() {
    let mut harness = harness();
    let cells = harness.cad_cells("cad_cells(list(x^2 - 1), list(x))");
    assert_eq!(cells.len(), 5);
    let patterns: Vec<Vec<bool>> = cells.iter().map(|c| c.open.clone()).collect();
    assert_eq!(
        patterns,
        vec![vec![true], vec![false], vec![true], vec![false], vec![true]]
    );
    let signs: Vec<i32> = cells.iter().map(|c| c.signs[0]).collect();
    assert_eq!(signs, vec![1, 0, -1, 0, 1]);
}

#[test]
fn cad_cells_degenerate_polynomials_and_variable_order() {
    let mut harness = harness();
    // The zero polynomial vanishes on every cell.
    let zero = harness.cad_cells("cad_cells(list(0), list(x))");
    assert_eq!(zero.len(), 1);
    assert_eq!(zero[0].signs, vec![0]);
    assert_eq!(zero[0].open, vec![true]);

    // A non-zero constant has one sign-invariant cell.
    let constant = harness.cad_cells("cad_cells(list(2), list(x))");
    assert_eq!(constant.len(), 1);
    assert_eq!(constant[0].signs, vec![1]);

    // A repeated input polynomial keeps one sign entry per input and the
    // square-free decomposition of the underlying factor.
    let repeated = harness.cad_cells("cad_cells(list(x^2 - 1, x^2 - 1), list(x))");
    assert_eq!(repeated.len(), 5);
    assert!(repeated.iter().all(|cell| cell.signs.len() == 2));
    let signs: Vec<i32> = repeated.iter().map(|cell| cell.signs[0]).collect();
    assert_eq!(signs, vec![1, 0, -1, 0, 1]);

    // The order of the listed variables selects the lifting order; the result
    // is a sign-invariant decomposition either way.
    for source in [
        "cad_cells(list(a*x - 1), list(a, x))",
        "cad_cells(list(a*x - 1), list(x, a))",
    ] {
        let cells = harness.cad_cells(source);
        assert!(!cells.is_empty(), "{source}");
        for cell in &cells {
            let value = cell.sample[0] * cell.sample[1] - 1.0;
            assert_eq!(
                cell.signs[0],
                sign_of(value),
                "{source} at {:?}",
                cell.sample
            );
        }
    }

    // No polynomials: one full-dimensional cell with no sign entries.
    let empty = harness.cad_cells("cad_cells(list(), list(x))");
    assert_eq!(empty.len(), 1);
    assert!(empty[0].signs.is_empty());
}

#[test]
fn malformed_and_mixed_requests_stay_unevaluated() {
    let mut harness = harness();
    for request in [
        "real_qe(gt(x, 0), 5)",
        "real_qe(lt(x, 0), list(1))",
        "cad_cells(list(x), list())",
        "cad_cells(x, list(x))",
        "cad_cells(list(x, y), list(x))",
        // A real comparison mixed with a Boolean atom is not a real formula;
        // it must not be decided with Boolean semantics.
        "exists(x, and(gt(x, 0), sin(x)))",
        "real_qe(exists(x, and(gt(x, 0), sin(x) = 0)), list(x))",
    ] {
        let (text, reduced) = harness.reduce(request);
        assert!(!reduced, "`{request}` must stay unevaluated: {text}");
    }
}

// ---------------------------------------------------------------------------
// Five variables.

const FIVE_VARIABLES: &str = "x1, x2, x3, x4, x5";

/// The five-dimensional unit ball, a linear chain and a two-parameter square
/// root are decided by the five-variable CAD. The returned quantifier-free
/// formulas are checked at dense rational parameter samples against the
/// geometric truth; for the square root that truth is `a + b >= 0`, since the
/// remaining coordinates can be zero.
#[test]
fn five_variable_quantifier_elimination() {
    let mut harness = harness();
    let variables = FIVE_VARIABLES;
    let ball = harness.answer(&format!(
        "real_qe(exists(list({variables}), lt(x1^2 + x2^2 + x3^2 + x4^2 + x5^2, 1)), list({variables}))"
    ));
    assert_eq!(ball, "true");
    let empty = harness.answer(&format!(
        "real_qe(exists(list({variables}), lt(x1^2 + x2^2 + x3^2 + x4^2 + x5^2, 0)), list({variables}))"
    ));
    assert_eq!(empty, "false");

    let chain = harness.answer(&format!(
        "real_qe(exists(list({variables}), and(gt(x2, x1), gt(x3, x2), gt(x4, x3), gt(x5, x4))), list({variables}))"
    ));
    assert_eq!(chain, "true");
    let cyclic = harness.answer(&format!(
        "real_qe(exists(list({variables}), and(gt(x2, x1), gt(x3, x2), gt(x4, x3), gt(x5, x4), gt(x1, x5))), list({variables}))"
    ));
    assert_eq!(cyclic, "false");

    let rooted = harness.answer(&format!(
        "real_qe(exists(list({variables}), and(x1 = a, x2 = b, x3^2 = a + b, x4 = 0, x5 = 0)), list({variables}))"
    ));

    // Rationals on a quarter grid: sums and products are exact in double
    // precision, so the comparison itself is exact at every sample.
    let samples: Vec<f64> = (-6..=6).map(|k| f64::from(k) / 4.0).collect();
    for &a in &samples {
        for &b in &samples {
            let bindings = [("a", a), ("b", b)];
            assert_eq!(
                harness.evaluate(&rooted, &bindings).unwrap_or(false),
                a + b >= 0.0,
                "rooted at a={a} b={b}: `{rooted}`"
            );
        }
    }
}

#[test]
fn cad_cells_five_variables_are_sign_invariant() {
    let mut harness = harness();
    let source =
        format!("cad_cells(list(x1^2 + x2^2 + x3^2 + x4^2 + x5^2 - 1), list({FIVE_VARIABLES}))");
    let cells = harness.cad_cells(&source);
    assert!(cells.len() > 30, "{}", cells.len());
    let mut all_open = 0_usize;
    let mut sections = 0_usize;
    for cell in &cells {
        assert_eq!(cell.sample.len(), 5, "{:?}", cell.sample);
        assert_eq!(cell.open.len(), 5);
        assert_eq!(cell.signs.len(), 1);
        let value = cell.sample.iter().map(|x| x * x).sum::<f64>() - 1.0;
        assert_eq!(
            cell.signs[0],
            sign_of(value),
            "sample sign at {:?}",
            cell.sample
        );
        if cell.open.iter().all(|&open| open) {
            all_open = all_open.saturating_add(1);
            // Only a full-dimensional cell is an open neighbourhood; a mixed
            // cell has section coordinates that constrain the others, so it
            // must not be perturbed as a product.
            for (index, _) in cell.open.iter().enumerate() {
                for delta in [1e-4, -1e-4] {
                    let mut point = cell.sample.clone();
                    if let Some(coordinate) = point.get_mut(index) {
                        *coordinate += delta;
                    }
                    let value = point.iter().map(|x| x * x).sum::<f64>() - 1.0;
                    assert_eq!(
                        sign_of(value),
                        cell.signs[0],
                        "perturbation {delta} of coordinate {index} in {:?}",
                        cell.sample
                    );
                }
            }
        }
        if cell.open.iter().any(|&open| !open) {
            sections = sections.saturating_add(1);
        }
    }
    assert!(all_open > 0, "the decomposition must have open cells");
    assert!(sections > 0, "the decomposition must have sections");
    assert!(
        cells.iter().any(|cell| cell.signs[0] == 0),
        "the sphere itself must occur as a cell"
    );
}

/// A cubic section embedded in five variables: the sample point of the
/// section coordinate is an exact `root_of`, and lifting the other four
/// variables keeps the sign conditions.
#[test]
fn cad_cells_five_variables_use_root_of() {
    let mut harness = harness();
    let source = format!("cad_cells(list(x1^3 - 2), list({FIVE_VARIABLES}))");
    let text = harness.answer(&source);
    assert!(text.contains("root_of("), "exact algebraic sample: {text}");
    let cells = harness.cad_cells(&source);
    assert_eq!(cells.len(), 3, "one sector, one section, one sector");
    let signs: Vec<i32> = cells.iter().map(|cell| cell.signs[0]).collect();
    assert_eq!(signs, vec![-1, 0, 1]);
    let mut root = None;
    for cell in &cells {
        assert_eq!(cell.sample.len(), 5);
        assert_eq!(cell.open.len(), 5);
        assert_eq!(cell.signs.len(), 1);
        let x = cell.sample[0];
        assert_eq!(
            cell.signs[0],
            sign_of(x * x * x - 2.0),
            "sample sign at {:?}",
            cell.sample
        );
        if !cell.open[0] {
            root = Some(x);
        } else {
            // The outer sectors may be perturbed in the cubic coordinate.
            for delta in [1e-4, -1e-4] {
                let moved = x + delta;
                assert_eq!(
                    sign_of(moved * moved * moved - 2.0),
                    cell.signs[0],
                    "perturbation {delta} at {x}"
                );
            }
        }
        // The four lifted variables are sectors: perturbing them, or moving
        // inside the cubic sector, leaves the sign unchanged because the
        // polynomial does not involve them.
        for index in 1..5 {
            assert!(cell.open[index], "{:?}", cell.open);
            for delta in [1e-4, -1e-4] {
                let mut point = cell.sample.clone();
                if let Some(coordinate) = point.get_mut(index) {
                    *coordinate += delta;
                }
                let moved = point[0];
                assert_eq!(
                    sign_of(moved * moved * moved - 2.0),
                    cell.signs[0],
                    "perturbation of lifted coordinate {index}"
                );
            }
        }
    }
    let root = root.expect("the section cell");
    assert!((root - 2.0_f64.cbrt()).abs() < 1e-6, "root sample {root}");
}

#[test]
fn performance_is_bounded() {
    let mut harness = harness();
    let requests = [
        "real_qe(exists(y, lt(x^2 + y^2, 1)), list(y))",
        "real_qe(forall(y, ge(y^2, x)), list(y))",
        "real_qe(exists(x, a*x^2 + b*x + c = 0), list(x))",
        "real_qe(forall(a, exists(x, gt(x^2 + a, 0))), list(a))",
        "cad_cells(list(x^2 + y^2 - 1), list(x, y))",
    ];
    for request in requests {
        let start = Instant::now();
        let (text, reduced) = harness.reduce(request);
        let elapsed = start.elapsed();
        assert!(reduced, "`{request}`: {text}");
        assert!(elapsed.as_secs_f64() < 20.0, "`{request}` took {elapsed:?}");
    }
}

#[test]
fn non_polynomial_requests_stay_unevaluated() {
    let mut harness = harness();
    for request in [
        "real_qe(exists(x, sin(x) = 0), list(x))",
        "exists(x, gt(exp(x), 0))",
        "real_qe(exists(x, gt(1/x, 0)), list(x))",
    ] {
        let (text, reduced) = harness.reduce(request);
        assert!(!reduced, "`{request}` must stay unevaluated: {text}");
    }
}
