//! # Tape IR
//!
//! The common input of every compiling backend: a typed SSA program over
//! `f64` registers with input, parameter and output tables. It is lowered
//! straight from the graph's DAG ([`lower`]) and optimised once
//! ([`Tape::optimise`]); the Cranelift JIT and the reference interpreters
//! below all consume it.
//!
//! Lowering undoes the canonical form of the engine, which only has
//! n-ary `add`/`mul` and `pow`:
//!
//! | graph form | instruction |
//! |---|---|
//! | `add(a, mul(-c, b))` | `Sub(a, c·b)` |
//! | `mul(-1, a)` | `Neg(a)` |
//! | `mul(a, pow(b, -k))` | `Div(a, b^k)` |
//! | `pow(x, n)`, small integer `n` | `PowI` (square-and-multiply) |
//! | `pow(x, k/2)`, `pow(x, k/3)` | `Sqrt` / `Cbrt` then `PowI` |
//! | n-ary `add`/`mul` | balanced pairwise trees (bounded rounding error) |
//!
//! Operators are matched first through the [`Lowering`](super::Lowering)
//! attribute, then by name for the elementary functions, and otherwise
//! fall back to calling their [`EvalFn`] through a shim, so every
//! operator with scalar semantics compiles.
//!
//! Besides plain evaluation the same program can be interpreted over
//! dual numbers ([`Tape::gradient`], forward-mode automatic
//! differentiation) and over outward-rounded intervals ([`Tape::enclose`]),
//! which turns a numeric witness into a rigorous enclosure.

use std::collections::HashMap;
use std::hash::Hash;
use std::hash::Hasher;

use num_traits::Signed;

use super::BackendError;
use super::Intrinsic;
use super::Lowering;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::SymbolId;
use crate::graph::op::EvalFn;
use crate::graph::op::core;

/// A register: the index of the instruction that produced it.
pub type Reg = u32;

/// Comparisons; they produce `1.0` or `0.0`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Cmp {
    /// `a < b`.
    Lt,
    /// `a <= b`.
    Le,
    /// `a > b`.
    Gt,
    /// `a >= b`.
    Ge,
    /// `a == b`.
    Eq,
    /// `a != b` (true for NaN operands).
    Ne,
}

/// One instruction of the tape.
#[derive(Clone, Copy, Debug)]
pub enum Inst {
    /// A constant.
    Const(f64),
    /// Input `i`.
    Input(u32),
    /// Parameter `i` (constant during a batch).
    Param(u32),
    /// `a + b`.
    Add(Reg, Reg),
    /// `a - b`.
    Sub(Reg, Reg),
    /// `a * b`.
    Mul(Reg, Reg),
    /// `a / b`.
    Div(Reg, Reg),
    /// `-a`.
    Neg(Reg),
    /// `a^n` for an integer `n`.
    PowI(Reg, i32),
    /// `a^b` (`powf`).
    PowF(Reg, Reg),
    /// An elementary function of one argument.
    Unary(Intrinsic, Reg),
    /// `atan2(a, b)`.
    Atan2(Reg, Reg),
    /// A comparison.
    Cmp(Cmp, Reg, Reg),
    /// `if c != 0 { a } else { b }`.
    Select(Reg, Reg, Reg),
    /// A call of an operator's scalar semantics: the function and the
    /// range `args[start..start + len]` of argument registers.
    Call(EvalFn, u32, u32),
    /// An external C function of one argument.
    Extern1(extern "C" fn(f64) -> f64, Reg),
    /// An external C function of two arguments.
    Extern2(extern "C" fn(f64, f64) -> f64, Reg, Reg),
}

/// A lowered program: instructions in dependency order, argument lists of
/// calls, and the output registers.
#[derive(Clone, Debug, Default)]
pub struct Tape {
    /// The instructions; register `r` is the value of `insts[r]`.
    pub insts: Vec<Inst>,
    /// Argument registers of [`Inst::Call`].
    pub args: Vec<Reg>,
    /// Number of inputs.
    pub inputs: usize,
    /// Number of parameters.
    pub params: usize,
    /// Output registers.
    pub outputs: Vec<Reg>,
}

/// Floating-point relaxations the optimiser may apply (all off by default,
/// which keeps IEEE semantics and NaN propagation).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct FastMath {
    /// Reassociate n-ary sums and products freely.
    pub reassoc: bool,
    /// Replace `a / b` by `a * (1/b)` and share reciprocals.
    pub reciprocal: bool,
}

fn key_of(inst: &Inst) -> (u8, u64, u64, u64) {
    let r = |x: Reg| u64::from(x);
    match *inst {
        | Inst::Const(c) => (0, c.to_bits(), 0, 0),
        | Inst::Input(i) => (1, u64::from(i), 0, 0),
        | Inst::Param(i) => (2, u64::from(i), 0, 0),
        | Inst::Add(a, b) => (3, r(a.min(b)), r(a.max(b)), 0),
        | Inst::Sub(a, b) => (4, r(a), r(b), 0),
        | Inst::Mul(a, b) => (5, r(a.min(b)), r(a.max(b)), 0),
        | Inst::Div(a, b) => (6, r(a), r(b), 0),
        | Inst::Neg(a) => (7, r(a), 0, 0),
        | Inst::PowI(a, n) => (8, r(a), u64::from(n.cast_unsigned()), 0),
        | Inst::PowF(a, b) => (9, r(a), r(b), 0),
        | Inst::Unary(f, a) => (10, f as u64, r(a), 0),
        | Inst::Atan2(a, b) => (11, r(a), r(b), 0),
        | Inst::Cmp(c, a, b) => (12, c as u64, r(a), r(b)),
        | Inst::Select(c, a, b) => (13, r(c), r(a) << 32 | r(b), 0),
        | Inst::Call(f, s, l) => (14, f as usize as u64, u64::from(s), u64::from(l)),
        | Inst::Extern1(f, a) => (15, f as usize as u64, r(a), 0),
        | Inst::Extern2(f, a, b) => (16, f as usize as u64, r(a), r(b)),
    }
}

/// Builds a tape with hash-consing of identical instructions.
struct Builder {
    tape: Tape,
    seen: HashMap<(u8, u64, u64, u64), Reg>,
}

impl Builder {
    fn push(
        &mut self,
        inst: Inst,
    ) -> Reg {
        // Calls with argument lists are keyed by their arguments too.
        let mut key = key_of(&inst);
        if let Inst::Call(_, start, len) = inst {
            let mut h = std::collections::hash_map::DefaultHasher::new();
            self.tape.args.get(start as usize..(start as usize).saturating_add(len as usize)).hash(&mut h);
            key.2 = h.finish();
        }
        if let Some(&r) = self.seen.get(&key) {
            return r;
        }
        let r = Reg::try_from(self.tape.insts.len()).unwrap_or(Reg::MAX);
        self.tape.insts.push(inst);
        self.seen.insert(key, r);
        r
    }

    fn constant(
        &mut self,
        c: f64,
    ) -> Reg {
        self.push(Inst::Const(c))
    }

    /// Balanced pairwise reduction.
    fn tree(
        &mut self,
        mut items: Vec<Reg>,
        combine: fn(Reg, Reg) -> Inst,
    ) -> Option<Reg> {
        while items.len() > 1 {
            let mut next = Vec::with_capacity(items.len().div_ceil(2));
            for pair in items.chunks(2) {
                match *pair {
                    | [a, b] => next.push(self.push(combine(a, b))),
                    | [a] => next.push(a),
                    | _ => {},
                }
            }
            items = next;
        }
        items.first().copied()
    }

    fn pow_int(
        &mut self,
        base: Reg,
        n: i64,
    ) -> Reg {
        match n {
            | 0 => self.constant(1.0),
            | 1 => base,
            | _ => self.push(Inst::PowI(base, i32::try_from(n).unwrap_or(i32::MAX))),
        }
    }
}

fn by_name(name: &str) -> Option<Intrinsic> {
    Some(match name {
        | "exp" => Intrinsic::Exp,
        | "ln" => Intrinsic::Ln,
        | "sin" => Intrinsic::Sin,
        | "cos" => Intrinsic::Cos,
        | "tan" => Intrinsic::Tan,
        | "asin" => Intrinsic::Asin,
        | "acos" => Intrinsic::Acos,
        | "atan" => Intrinsic::Atan,
        | "sinh" => Intrinsic::Sinh,
        | "cosh" => Intrinsic::Cosh,
        | "tanh" => Intrinsic::Tanh,
        | "asinh" => Intrinsic::Asinh,
        | "acosh" => Intrinsic::Acosh,
        | "atanh" => Intrinsic::Atanh,
        | "sqrt" => Intrinsic::Sqrt,
        | "abs" => Intrinsic::Abs,
        | "floor" => Intrinsic::Floor,
        | "ceil" => Intrinsic::Ceil,
        | "cbrt" => Intrinsic::Cbrt,
        | _ => return None,
    })
}

fn comparison(name: &str) -> Option<Cmp> {
    Some(match name {
        | "lt" => Cmp::Lt,
        | "le" => Cmp::Le,
        | "gt" => Cmp::Gt,
        | "ge" => Cmp::Ge,
        | "ne" => Cmp::Ne,
        | _ => return None,
    })
}

/// Lowers the concrete terms `roots` to one tape with an output per root,
/// as functions of `inputs` and `params` (in that order).
///
/// # Errors
/// An unbound symbol, a non-numeric literal or an operator without scalar
/// semantics.
pub fn lower(
    graph: &Graph,
    roots: &[NodeId],
    inputs: &[SymbolId],
    params: &[SymbolId],
) -> Result<Tape, BackendError> {
    let mut b = Builder { tape: Tape { inputs: inputs.len(), params: params.len(), ..Tape::default() }, seen: HashMap::new() };
    let mut slot: HashMap<NodeId, Reg> = HashMap::new();
    for &root in roots {
        let r = lower_node(graph, root, inputs, params, &mut b, &mut slot)?;
        b.tape.outputs.push(r);
    }
    Ok(b.tape)
}

/// A numeric literal of `node`, if it is one.
fn literal(
    graph: &Graph,
    node: NodeId,
) -> Option<f64> {
    graph.as_number(node).map(crate::graph::Number::to_f64)
}

fn lower_node(
    graph: &Graph,
    root: NodeId,
    inputs: &[SymbolId],
    params: &[SymbolId],
    b: &mut Builder,
    slot: &mut HashMap<NodeId, Reg>,
) -> Result<Reg, BackendError> {
    // Iterative post-order so that deep terms cannot overflow the stack.
    let mut stack = vec![(root, false)];
    while let Some((node, expanded)) = stack.pop() {
        if slot.contains_key(&node) {
            continue;
        }
        let children = graph.children(node);
        if !expanded && !children.is_empty() {
            stack.push((node, true));
            stack.extend(children.iter().filter(|c| !slot.contains_key(c)).map(|&c| (c, false)));
            continue;
        }
        let reg = lower_one(graph, node, inputs, params, b, slot)?;
        slot.insert(node, reg);
    }
    slot.get(&root).copied().ok_or(BackendError::NotNumeric)
}

fn lower_one(
    graph: &Graph,
    node: NodeId,
    inputs: &[SymbolId],
    params: &[SymbolId],
    b: &mut Builder,
    slot: &HashMap<NodeId, Reg>,
) -> Result<Reg, BackendError> {
    if let Some(c) = literal(graph, node) {
        return Ok(b.constant(c));
    }
    if let Some(symbol) = graph.as_symbol(node) {
        if let Some(i) = inputs.iter().position(|&s| s == symbol) {
            return Ok(b.push(Inst::Input(u32::try_from(i).unwrap_or(u32::MAX))));
        }
        if let Some(i) = params.iter().position(|&s| s == symbol) {
            return Ok(b.push(Inst::Param(u32::try_from(i).unwrap_or(u32::MAX))));
        }
        return Err(BackendError::UnboundSymbol(graph.interner().symbol_name(symbol).to_owned()));
    }
    if graph.payload(node).is_some() {
        return Err(BackendError::NotNumeric);
    }
    let op = graph.op(node);
    let children = graph.children(node);
    let reg = |c: &NodeId| slot.get(c).copied().ok_or(BackendError::NotNumeric);
    let args: Vec<Reg> = children.iter().map(reg).collect::<Result<_, _>>()?;
    if op == core::ADD {
        let (mut positive, mut negative) = (Vec::new(), Vec::new());
        for &child in children {
            match negated(graph, child, b, slot)? {
                | Some(r) => negative.push(r),
                | None => positive.push(reg(&child)?),
            }
        }
        let p = b.tree(positive, Inst::Add);
        let n = b.tree(negative, Inst::Add);
        return Ok(match (p, n) {
            | (Some(p), Some(n)) => b.push(Inst::Sub(p, n)),
            | (Some(p), None) => p,
            | (None, Some(n)) => b.push(Inst::Neg(n)),
            | (None, None) => b.constant(0.0),
        });
    }
    if op == core::MUL {
        let (mut numerator, mut denominator, mut sign) = (Vec::new(), Vec::new(), false);
        for &child in children {
            if let Some(c) = literal(graph, child) {
                if c < 0.0 {
                    sign = !sign;
                    if (c + 1.0).abs() > 0.0 {
                        numerator.push(b.constant(-c));
                    }
                } else {
                    numerator.push(b.constant(c));
                }
                continue;
            }
            if let (true, &[base, e]) = (graph.op(child) == core::POW, graph.children(child))
                && let Some(k) = literal(graph, e).filter(|k| *k < 0.0) {
                    let base = reg(&base)?;
                    denominator.push(power(graph, b, base, e, -k));
                    continue;
                }
            numerator.push(reg(&child)?);
        }
        let n = b.tree(numerator, Inst::Mul).unwrap_or_else(|| b.constant(1.0));
        let value = match b.tree(denominator, Inst::Mul) {
            | Some(d) => b.push(Inst::Div(n, d)),
            | None => n,
        };
        return Ok(if sign { b.push(Inst::Neg(value)) } else { value });
    }
    if op == core::POW {
        let (&[base_node, e], [base, exponent]) = (children, args.as_slice()) else {
            return Err(BackendError::NoSemantics("pow".to_owned()));
        };
        let _ = base_node;
        return match literal(graph, e) {
            | Some(k) if k < 0.0 => {
                let positive = power(graph, b, *base, e, -k);
                let one = b.constant(1.0);
                Ok(b.push(Inst::Div(one, positive)))
            },
            | Some(k) => Ok(power(graph, b, *base, e, k)),
            | None => Ok(b.push(Inst::PowF(*base, *exponent))),
        };
    }
    let desc = graph.ops().get(op);
    if let Some(lowering) = graph.ops().attr::<Lowering>(op) {
        match (lowering, args.as_slice()) {
            | (Lowering::Intrinsic(f), &[a]) => return Ok(b.push(Inst::Unary(*f, a))),
            | (Lowering::Extern1(f), &[a]) => return Ok(b.push(Inst::Extern1(*f, a))),
            | (Lowering::Extern2(f), &[a, c]) => return Ok(b.push(Inst::Extern2(*f, a, c))),
            | _ => {},
        }
    }
    let name: &str = &desc.name;
    match (by_name(name), comparison(name), name, args.as_slice()) {
        | (Some(f), _, _, &[a]) => return Ok(b.push(Inst::Unary(f, a))),
        | (_, Some(c), _, &[x, y]) => return Ok(b.push(Inst::Cmp(c, x, y))),
        | (_, _, "atan2", &[y, x]) => return Ok(b.push(Inst::Atan2(y, x))),
        | _ => {},
    }
    if op == core::EQ
        && let &[x, y] = args.as_slice() {
            return Ok(b.push(Inst::Cmp(Cmp::Eq, x, y)));
        }
    let f = desc.eval.ok_or_else(|| BackendError::NoSemantics(desc.name.to_string()))?;
    if args.is_empty() {
        // Constants such as `pi`.
        return Ok(b.constant(f(&[])));
    }
    let start = u32::try_from(b.tape.args.len()).unwrap_or(u32::MAX);
    b.tape.args.extend_from_slice(&args);
    Ok(b.push(Inst::Call(f, start, u32::try_from(args.len()).unwrap_or(u32::MAX))))
}

/// `c · rest` for a summand `mul(c, rest...)` with a negative literal `c`:
/// the register of `|c| · rest`, else `None`.
fn negated(
    graph: &Graph,
    child: NodeId,
    b: &mut Builder,
    slot: &HashMap<NodeId, Reg>,
) -> Result<Option<Reg>, BackendError> {
    if let Some(c) = literal(graph, child) {
        return Ok((c < 0.0).then(|| b.constant(-c)));
    }
    if graph.op(child) != core::MUL {
        return Ok(None);
    }
    let factors = graph.children(child);
    let Some(c) = factors.iter().find_map(|&f| literal(graph, f)).filter(|c| *c < 0.0) else {
        return Ok(None);
    };
    let mut rest = Vec::new();
    let mut skipped = false;
    let mut denominator = Vec::new();
    for &f in factors {
        if !skipped && literal(graph, f) == Some(c) {
            skipped = true;
            continue;
        }
        if let (true, &[base, e]) = (graph.op(f) == core::POW, graph.children(f))
            && let Some(k) = literal(graph, e).filter(|k| *k < 0.0) {
                let base = slot.get(&base).copied().ok_or(BackendError::NotNumeric)?;
                denominator.push(power(graph, b, base, e, -k));
                continue;
            }
        rest.push(slot.get(&f).copied().ok_or(BackendError::NotNumeric)?);
    }
    if (c + 1.0).abs() > 0.0 {
        rest.push(b.constant(-c));
    }
    let n = b.tree(rest, Inst::Mul).unwrap_or_else(|| b.constant(1.0));
    Ok(Some(match b.tree(denominator, Inst::Mul) {
        | Some(d) => b.push(Inst::Div(n, d)),
        | None => n,
    }))
}

/// `base^k` for a literal `k >= 0` (the node `e` gives its exact value).
fn power(
    graph: &Graph,
    b: &mut Builder,
    base: Reg,
    e: NodeId,
    k: f64,
) -> Reg {
    let exact = graph.as_number(e).and_then(crate::graph::Number::to_rational).map(|r| r.abs());
    if let Some(r) = exact {
        let (numer, denom) = (r.numer().clone(), r.denom().clone());
        if let (Ok(p), Ok(q)) = (i64::try_from(numer), i64::try_from(denom))
            && p <= 64 {
                match q {
                    | 1 => return b.pow_int(base, p),
                    | 2 => {
                        let s = b.push(Inst::Unary(Intrinsic::Sqrt, base));
                        return b.pow_int(s, p);
                    },
                    | 3 => {
                        let s = b.push(Inst::Unary(Intrinsic::Cbrt, base));
                        return b.pow_int(s, p);
                    },
                    | _ => {},
                }
            }
    }
    let k = b.constant(k);
    b.push(Inst::PowF(base, k))
}

impl Tape {
    /// Constant folding, then dead-code elimination (registers renumbered).
    #[must_use]
    pub fn optimise(&self) -> Self {
        // Fold instructions whose operands are all constants.
        let mut folded = self.clone();
        let mut constants: Vec<Option<f64>> = Vec::with_capacity(folded.insts.len());
        for i in 0..folded.insts.len() {
            let inst = folded.insts.get(i).copied();
            let value = inst.and_then(|inst| {
                let operands = operands(&inst, &folded.args);
                if matches!(inst, Inst::Input(_) | Inst::Param(_)) {
                    return None;
                }
                if let Inst::Const(c) = inst {
                    return Some(c);
                }
                let values: Option<Vec<f64>> = operands.iter().map(|&r| constants.get(r as usize).copied().flatten()).collect();
                let values = values?;
                Some(apply::<f64>(&inst, &values))
            });
            if let (Some(v), Some(slot)) = (value, folded.insts.get_mut(i)) {
                *slot = Inst::Const(v);
            }
            constants.push(value);
        }
        // Liveness from the outputs.
        let mut live = vec![false; folded.insts.len()];
        let mut work: Vec<Reg> = folded.outputs.clone();
        while let Some(r) = work.pop() {
            match live.get_mut(r as usize) {
                | Some(flag) if !*flag => {
                    *flag = true;
                    if let Some(inst) = folded.insts.get(r as usize) {
                        work.extend(operands(inst, &folded.args));
                    }
                },
                | _ => {},
            }
        }
        let mut map = vec![Reg::MAX; folded.insts.len()];
        let mut out = Self { inputs: folded.inputs, params: folded.params, ..Self::default() };
        for (i, inst) in folded.insts.iter().enumerate() {
            if !live.get(i).copied().unwrap_or(false) {
                continue;
            }
            let m = |r: Reg| map.get(r as usize).copied().unwrap_or(Reg::MAX);
            let renamed = match *inst {
                | Inst::Add(a, b) => Inst::Add(m(a), m(b)),
                | Inst::Sub(a, b) => Inst::Sub(m(a), m(b)),
                | Inst::Mul(a, b) => Inst::Mul(m(a), m(b)),
                | Inst::Div(a, b) => Inst::Div(m(a), m(b)),
                | Inst::Neg(a) => Inst::Neg(m(a)),
                | Inst::PowI(a, n) => Inst::PowI(m(a), n),
                | Inst::PowF(a, b) => Inst::PowF(m(a), m(b)),
                | Inst::Unary(f, a) => Inst::Unary(f, m(a)),
                | Inst::Atan2(a, b) => Inst::Atan2(m(a), m(b)),
                | Inst::Cmp(c, a, b) => Inst::Cmp(c, m(a), m(b)),
                | Inst::Select(c, a, b) => Inst::Select(m(c), m(a), m(b)),
                | Inst::Extern1(f, a) => Inst::Extern1(f, m(a)),
                | Inst::Extern2(f, a, b) => Inst::Extern2(f, m(a), m(b)),
                | Inst::Call(f, start, len) => {
                    let new_start = u32::try_from(out.args.len()).unwrap_or(u32::MAX);
                    let range = start as usize..(start as usize).saturating_add(len as usize);
                    for &r in folded.args.get(range).unwrap_or(&[]) {
                        out.args.push(m(r));
                    }
                    Inst::Call(f, new_start, len)
                },
                | other => other,
            };
            if let Some(slot) = map.get_mut(i) {
                *slot = Reg::try_from(out.insts.len()).unwrap_or(Reg::MAX);
            }
            out.insts.push(renamed);
        }
        out.outputs = folded.outputs.iter().map(|&r| map.get(r as usize).copied().unwrap_or(Reg::MAX)).collect();
        out
    }

    /// A structural hash, independent of node ids (for compilation caches).
    #[must_use]
    pub fn structural_hash(&self) -> u64 {
        let mut h = std::collections::hash_map::DefaultHasher::new();
        self.inputs.hash(&mut h);
        self.params.hash(&mut h);
        for inst in &self.insts {
            key_of(inst).hash(&mut h);
        }
        self.args.hash(&mut h);
        self.outputs.hash(&mut h);
        h.finish()
    }

    /// Evaluates every output over an arbitrary number type.
    pub fn run<T: Scalar>(
        &self,
        inputs: &[T],
        params: &[T],
        values: &mut Vec<T>,
        out: &mut [T],
    ) {
        values.clear();
        let mut scratch: Vec<T> = Vec::new();
        for inst in &self.insts {
            let value = match *inst {
                | Inst::Input(i) => inputs.get(i as usize).copied().unwrap_or_else(T::nan),
                | Inst::Param(i) => params.get(i as usize).copied().unwrap_or_else(T::nan),
                | _ => {
                    scratch.clear();
                    scratch.extend(operands(inst, &self.args).iter().map(|&r| values.get(r as usize).copied().unwrap_or_else(T::nan)));
                    apply(inst, &scratch)
                },
            };
            values.push(value);
        }
        for (slot, &r) in out.iter_mut().zip(&self.outputs) {
            *slot = values.get(r as usize).copied().unwrap_or_else(T::nan);
        }
    }

    /// Plain evaluation of every output.
    pub fn eval(
        &self,
        inputs: &[f64],
        params: &[f64],
        out: &mut [f64],
    ) {
        let mut values = Vec::with_capacity(self.insts.len());
        self.run(inputs, params, &mut values, out);
    }

    /// The gradient of output `k` with respect to every input at `inputs`
    /// (forward mode, one pass per input).
    #[must_use]
    pub fn gradient(
        &self,
        k: usize,
        inputs: &[f64],
        params: &[f64],
    ) -> Vec<f64> {
        let params: Vec<Dual> = params.iter().map(|&p| Dual(p, 0.0)).collect();
        let mut values = Vec::with_capacity(self.insts.len());
        let mut out = vec![Dual(0.0, 0.0); self.outputs.len()];
        (0..inputs.len())
            .map(|i| {
                let seeded: Vec<Dual> = inputs.iter().enumerate().map(|(j, &x)| Dual(x, if i == j { 1.0 } else { 0.0 })).collect();
                self.run(&seeded, &params, &mut values, &mut out);
                out.get(k).map_or(f64::NAN, |d| d.1)
            })
            .collect()
    }

    /// Rigorous enclosures of every output over the input boxes `inputs`
    /// (interval arithmetic with outward rounding; calls of operator
    /// semantics without interval rules give the whole line).
    #[must_use]
    pub fn enclose(
        &self,
        inputs: &[Interval],
        params: &[Interval],
    ) -> Vec<Interval> {
        let mut values = Vec::with_capacity(self.insts.len());
        let mut out = vec![Interval::ENTIRE; self.outputs.len()];
        self.run(inputs, params, &mut values, &mut out);
        out
    }
}

fn operands(
    inst: &Inst,
    args: &[Reg],
) -> Vec<Reg> {
    match *inst {
        | Inst::Const(_) | Inst::Input(_) | Inst::Param(_) => Vec::new(),
        | Inst::Neg(a) | Inst::PowI(a, _) | Inst::Unary(_, a) | Inst::Extern1(_, a) => vec![a],
        | Inst::Add(a, b)
        | Inst::Sub(a, b)
        | Inst::Mul(a, b)
        | Inst::Div(a, b)
        | Inst::PowF(a, b)
        | Inst::Atan2(a, b)
        | Inst::Cmp(_, a, b)
        | Inst::Extern2(_, a, b) => vec![a, b],
        | Inst::Select(c, a, b) => vec![c, a, b],
        | Inst::Call(_, start, len) => args.get(start as usize..(start as usize).saturating_add(len as usize)).unwrap_or(&[]).to_vec(),
    }
}

/// Number types a tape can be interpreted over.
pub trait Scalar: Copy + FromExtern {
    /// NaN (or the whole line).
    #[must_use]
    fn nan() -> Self;
    /// A constant.
    #[must_use]
    fn constant(c: f64) -> Self;
    /// `a + b`.
    #[must_use]
    fn add(self, b: Self) -> Self;
    /// `a - b`.
    #[must_use]
    fn sub(self, b: Self) -> Self;
    /// `a * b`.
    #[must_use]
    fn mul(self, b: Self) -> Self;
    /// `a / b`.
    #[must_use]
    fn div(self, b: Self) -> Self;
    /// `-a`.
    #[must_use]
    fn neg(self) -> Self;
    /// `a^n`.
    #[must_use]
    fn powi(self, n: i32) -> Self;
    /// `a^b`.
    #[must_use]
    fn powf(self, b: Self) -> Self;
    /// An elementary function.
    #[must_use]
    fn unary(self, f: Intrinsic) -> Self;
    /// `atan2(self, x)`.
    #[must_use]
    fn atan2(self, x: Self) -> Self;
    /// A comparison (`1` or `0`).
    #[must_use]
    fn compare(self, c: Cmp, b: Self) -> Self;
    /// `if self != 0 { a } else { b }`.
    #[must_use]
    fn select(self, a: Self, b: Self) -> Self;
    /// A call of scalar semantics.
    fn call(f: EvalFn, args: &[Self]) -> Self;
}

fn apply<T: Scalar>(
    inst: &Inst,
    v: &[T],
) -> T {
    let get = |i: usize| v.get(i).copied().unwrap_or_else(T::nan);
    match *inst {
        | Inst::Const(c) => T::constant(c),
        | Inst::Input(_) | Inst::Param(_) => T::nan(),
        | Inst::Add(..) => get(0).add(get(1)),
        | Inst::Sub(..) => get(0).sub(get(1)),
        | Inst::Mul(..) => get(0).mul(get(1)),
        | Inst::Div(..) => get(0).div(get(1)),
        | Inst::Neg(_) => get(0).neg(),
        | Inst::PowI(_, n) => get(0).powi(n),
        | Inst::PowF(..) => get(0).powf(get(1)),
        | Inst::Unary(f, _) => get(0).unary(f),
        | Inst::Atan2(..) => get(0).atan2(get(1)),
        | Inst::Cmp(c, ..) => get(0).compare(c, get(1)),
        | Inst::Select(..) => get(0).select(get(1), get(2)),
        | Inst::Call(f, ..) => T::call(f, v),
        | Inst::Extern1(f, _) => T::from_extern1(f, get(0)),
        | Inst::Extern2(f, ..) => T::from_extern2(f, get(0), get(1)),
    }
}

/// How a number type evaluates external functions.
pub trait FromExtern: Sized {
    /// `f(a)`.
    fn from_extern1(f: extern "C" fn(f64) -> f64, a: Self) -> Self;
    /// `f(a, b)`.
    fn from_extern2(f: extern "C" fn(f64, f64) -> f64, a: Self, b: Self) -> Self;
}

impl Intrinsic {
    /// The function on `f64`.
    #[must_use]
    pub fn apply(
        self,
        x: f64,
    ) -> f64 {
        match self {
            | Self::Exp => x.exp(),
            | Self::Ln => x.ln(),
            | Self::Sin => x.sin(),
            | Self::Cos => x.cos(),
            | Self::Tan => x.tan(),
            | Self::Asin => x.asin(),
            | Self::Acos => x.acos(),
            | Self::Atan => x.atan(),
            | Self::Sinh => x.sinh(),
            | Self::Cosh => x.cosh(),
            | Self::Tanh => x.tanh(),
            | Self::Asinh => x.asinh(),
            | Self::Acosh => x.acosh(),
            | Self::Atanh => x.atanh(),
            | Self::Sqrt => x.sqrt(),
            | Self::Abs => x.abs(),
            | Self::Floor => x.floor(),
            | Self::Ceil => x.ceil(),
            | Self::Cbrt => x.cbrt(),
        }
    }

    /// The derivative on `f64`.
    #[must_use]
    pub fn derivative(
        self,
        x: f64,
    ) -> f64 {
        match self {
            | Self::Exp => x.exp(),
            | Self::Ln => 1.0 / x,
            | Self::Sin => x.cos(),
            | Self::Cos => -x.sin(),
            | Self::Tan => 1.0 / (x.cos() * x.cos()),
            | Self::Asin => 1.0 / (1.0 - x * x).sqrt(),
            | Self::Acos => -1.0 / (1.0 - x * x).sqrt(),
            | Self::Atan => 1.0 / (1.0 + x * x),
            | Self::Sinh => x.cosh(),
            | Self::Cosh => x.sinh(),
            | Self::Tanh => 1.0 / (x.cosh() * x.cosh()),
            | Self::Asinh => 1.0 / (x * x + 1.0).sqrt(),
            | Self::Acosh => 1.0 / (x * x - 1.0).sqrt(),
            | Self::Atanh => 1.0 / (1.0 - x * x),
            | Self::Sqrt => 0.5 / x.sqrt(),
            | Self::Abs => x.signum(),
            | Self::Floor | Self::Ceil => 0.0,
            | Self::Cbrt => 1.0 / (3.0 * x.cbrt() * x.cbrt()),
        }
    }
}

const fn bit(b: bool) -> f64 {
    if b { 1.0 } else { 0.0 }
}

fn compare_f64(
    a: f64,
    c: Cmp,
    b: f64,
) -> bool {
    match c {
        | Cmp::Lt => a < b,
        | Cmp::Le => a <= b,
        | Cmp::Gt => a > b,
        | Cmp::Ge => a >= b,
        | Cmp::Eq => a.partial_cmp(&b) == Some(std::cmp::Ordering::Equal),
        | Cmp::Ne => a.partial_cmp(&b) != Some(std::cmp::Ordering::Equal),
    }
}

impl Scalar for f64 {
    fn nan() -> Self {
        Self::NAN
    }

    fn constant(c: f64) -> Self {
        c
    }

    fn add(
        self,
        b: Self,
    ) -> Self {
        self + b
    }

    fn sub(
        self,
        b: Self,
    ) -> Self {
        self - b
    }

    fn mul(
        self,
        b: Self,
    ) -> Self {
        self * b
    }

    fn div(
        self,
        b: Self,
    ) -> Self {
        self / b
    }

    fn neg(self) -> Self {
        -self
    }

    fn powi(
        self,
        n: i32,
    ) -> Self {
        int_pow(self, n)
    }

    fn powf(
        self,
        b: Self,
    ) -> Self {
        Self::powf(self, b)
    }

    fn unary(
        self,
        f: Intrinsic,
    ) -> Self {
        f.apply(self)
    }

    fn atan2(
        self,
        x: Self,
    ) -> Self {
        Self::atan2(self, x)
    }

    fn compare(
        self,
        c: Cmp,
        b: Self,
    ) -> Self {
        bit(compare_f64(self, c, b))
    }

    fn select(
        self,
        a: Self,
        b: Self,
    ) -> Self {
        if (self - 0.0).abs() > 0.0 || self.is_nan() { a } else { b }
    }

    fn call(
        f: EvalFn,
        args: &[Self],
    ) -> Self {
        f(args)
    }
}

impl FromExtern for f64 {
    fn from_extern1(
        f: extern "C" fn(f64) -> f64,
        a: Self,
    ) -> Self {
        f(a)
    }

    fn from_extern2(
        f: extern "C" fn(f64, f64) -> f64,
        a: Self,
        b: Self,
    ) -> Self {
        f(a, b)
    }
}

/// `x^n` by square-and-multiply (the same sequence the JIT emits).
#[must_use]
pub fn int_pow(
    x: f64,
    n: i32,
) -> f64 {
    let mut e = n.unsigned_abs();
    let mut base = x;
    let mut acc: Option<f64> = None;
    while e > 0 {
        if e & 1 == 1 {
            acc = Some(acc.map_or(base, |a| a * base));
        }
        e >>= 1;
        if e > 0 {
            base *= base;
        }
    }
    let value = acc.unwrap_or(1.0);
    if n < 0 { 1.0 / value } else { value }
}

/// A dual number `a + b ε`, `ε² = 0`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Dual(pub f64, pub f64);

impl Scalar for Dual {
    fn nan() -> Self {
        Self(f64::NAN, f64::NAN)
    }

    fn constant(c: f64) -> Self {
        Self(c, 0.0)
    }

    fn add(
        self,
        b: Self,
    ) -> Self {
        Self(self.0 + b.0, self.1 + b.1)
    }

    fn sub(
        self,
        b: Self,
    ) -> Self {
        Self(self.0 - b.0, self.1 - b.1)
    }

    fn mul(
        self,
        b: Self,
    ) -> Self {
        Self(self.0 * b.0, self.0 * b.1 + self.1 * b.0)
    }

    fn div(
        self,
        b: Self,
    ) -> Self {
        Self(self.0 / b.0, (self.1 * b.0 - self.0 * b.1) / (b.0 * b.0))
    }

    fn neg(self) -> Self {
        Self(-self.0, -self.1)
    }

    fn powi(
        self,
        n: i32,
    ) -> Self {
        let value = int_pow(self.0, n);
        let slope = if n == 0 { 0.0 } else { f64::from(n) * int_pow(self.0, n - 1) };
        Self(value, slope * self.1)
    }

    fn powf(
        self,
        b: Self,
    ) -> Self {
        let value = self.0.powf(b.0);
        let d_base = if self.1 == 0.0 { 0.0 } else { b.0 * self.0.powf(b.0 - 1.0) * self.1 };
        let d_exp = if b.1 == 0.0 { 0.0 } else { value * self.0.ln() * b.1 };
        Self(value, d_base + d_exp)
    }

    fn unary(
        self,
        f: Intrinsic,
    ) -> Self {
        Self(f.apply(self.0), f.derivative(self.0) * self.1)
    }

    fn atan2(
        self,
        x: Self,
    ) -> Self {
        let r2 = self.0 * self.0 + x.0 * x.0;
        Self(self.0.atan2(x.0), (x.0 * self.1 - self.0 * x.1) / r2)
    }

    fn compare(
        self,
        c: Cmp,
        b: Self,
    ) -> Self {
        Self(bit(compare_f64(self.0, c, b.0)), 0.0)
    }

    fn select(
        self,
        a: Self,
        b: Self,
    ) -> Self {
        if (self.0 - 0.0).abs() > 0.0 || self.0.is_nan() { a } else { b }
    }

    fn call(
        f: EvalFn,
        args: &[Self],
    ) -> Self {
        // Value exactly; directional derivative by central differences.
        let point: Vec<f64> = args.iter().map(|d| d.0).collect();
        let value = f(&point);
        let mut slope = 0.0;
        for (k, d) in args.iter().enumerate() {
            if d.1 == 0.0 {
                continue;
            }
            let h = 1e-6 * (1.0 + d.0.abs());
            let (mut up, mut down) = (point.clone(), point.clone());
            if let (Some(u), Some(w)) = (up.get_mut(k), down.get_mut(k)) {
                *u += h;
                *w -= h;
            }
            slope += (f(&up) - f(&down)) / (2.0 * h) * d.1;
        }
        Self(value, slope)
    }
}

impl FromExtern for Dual {
    fn from_extern1(
        f: extern "C" fn(f64) -> f64,
        a: Self,
    ) -> Self {
        let h = 1e-6 * (1.0 + a.0.abs());
        Self(f(a.0), (f(a.0 + h) - f(a.0 - h)) / (2.0 * h) * a.1)
    }

    fn from_extern2(
        f: extern "C" fn(f64, f64) -> f64,
        a: Self,
        b: Self,
    ) -> Self {
        let (ha, hb) = (1e-6 * (1.0 + a.0.abs()), 1e-6 * (1.0 + b.0.abs()));
        let da = (f(a.0 + ha, b.0) - f(a.0 - ha, b.0)) / (2.0 * ha);
        let db = (f(a.0, b.0 + hb) - f(a.0, b.0 - hb)) / (2.0 * hb);
        Self(f(a.0, b.0), da * a.1 + db * b.1)
    }
}

/// A closed interval `[lo, hi]` (NaN bounds mean "no information").
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Interval {
    /// Lower bound.
    pub lo: f64,
    /// Upper bound.
    pub hi: f64,
}

impl Interval {
    /// The whole line.
    pub const ENTIRE: Self = Self { lo: f64::NEG_INFINITY, hi: f64::INFINITY };

    /// The degenerate interval `[x, x]`.
    #[must_use]
    pub const fn point(x: f64) -> Self {
        Self { lo: x, hi: x }
    }

    /// `[lo, hi]` widened by one ulp on each side.
    const fn outward(
        lo: f64,
        hi: f64,
    ) -> Self {
        if lo.is_nan() || hi.is_nan() {
            return Self::ENTIRE;
        }
        Self { lo: lo.next_down(), hi: hi.next_up() }
    }

    /// Whether `x` lies in the interval.
    #[must_use]
    pub fn contains(
        &self,
        x: f64,
    ) -> bool {
        self.lo <= x && x <= self.hi
    }

    fn monotone(
        self,
        f: fn(f64) -> f64,
        increasing: bool,
    ) -> Self {
        let (a, b) = (f(self.lo), f(self.hi));
        if increasing { Self::outward(a, b) } else { Self::outward(b, a) }
    }

    fn contains_point_of(
        self,
        offset: f64,
        period: f64,
    ) -> bool {
        // Whether some x = offset + k period lies in [lo, hi].
        let k = ((self.lo - offset) / period).ceil();
        offset + k * period <= self.hi
    }
}

impl Scalar for Interval {
    fn nan() -> Self {
        Self::ENTIRE
    }

    fn constant(c: f64) -> Self {
        if c.is_finite() { Self::point(c) } else { Self { lo: c, hi: c } }
    }

    fn add(
        self,
        b: Self,
    ) -> Self {
        Self::outward(self.lo + b.lo, self.hi + b.hi)
    }

    fn sub(
        self,
        b: Self,
    ) -> Self {
        Self::outward(self.lo - b.hi, self.hi - b.lo)
    }

    fn mul(
        self,
        b: Self,
    ) -> Self {
        let p = [self.lo * b.lo, self.lo * b.hi, self.hi * b.lo, self.hi * b.hi];
        if p.iter().any(|v| v.is_nan()) {
            return Self::ENTIRE;
        }
        Self::outward(p.iter().copied().fold(f64::INFINITY, f64::min), p.iter().copied().fold(f64::NEG_INFINITY, f64::max))
    }

    fn div(
        self,
        b: Self,
    ) -> Self {
        if b.lo <= 0.0 && b.hi >= 0.0 {
            return Self::ENTIRE;
        }
        self.mul(Self { lo: 1.0 / b.hi, hi: 1.0 / b.lo }.widen())
    }

    fn neg(self) -> Self {
        Self { lo: -self.hi, hi: -self.lo }
    }

    fn powi(
        self,
        n: i32,
    ) -> Self {
        if n == 0 {
            return Self::point(1.0);
        }
        if n < 0 {
            return Self::point(1.0).div(self.powi(-n));
        }
        let (a, b) = (int_pow(self.lo, n), int_pow(self.hi, n));
        if n % 2 == 0 {
            if self.lo <= 0.0 && self.hi >= 0.0 {
                Self::outward(0.0, a.max(b)).clamp_low(0.0)
            } else {
                Self::outward(a.min(b), a.max(b))
            }
        } else {
            Self::outward(a, b)
        }
    }

    fn powf(
        self,
        b: Self,
    ) -> Self {
        // x^y = exp(y ln x) for x > 0.
        if self.lo <= 0.0 {
            return Self::ENTIRE;
        }
        b.mul(self.unary(Intrinsic::Ln)).unary(Intrinsic::Exp)
    }

    fn unary(
        self,
        f: Intrinsic,
    ) -> Self {
        use std::f64::consts::FRAC_PI_2;
        use std::f64::consts::PI;
        match f {
            | Intrinsic::Exp => self.monotone(f64::exp, true).clamp_low(0.0),
            | Intrinsic::Ln => {
                if self.hi <= 0.0 {
                    Self::ENTIRE
                } else {
                    Self { lo: if self.lo <= 0.0 { f64::NEG_INFINITY } else { self.lo.ln() }, hi: self.hi.ln() }.widen()
                }
            },
            | Intrinsic::Sqrt => {
                if self.hi < 0.0 {
                    Self::ENTIRE
                } else {
                    Self::outward(self.lo.max(0.0).sqrt(), self.hi.sqrt()).clamp_low(0.0)
                }
            },
            | Intrinsic::Cbrt => self.monotone(f64::cbrt, true),
            | Intrinsic::Atan => self.monotone(f64::atan, true),
            | Intrinsic::Sinh => self.monotone(f64::sinh, true),
            | Intrinsic::Tanh => self.monotone(f64::tanh, true),
            | Intrinsic::Asinh => self.monotone(f64::asinh, true),
            | Intrinsic::Floor => Self { lo: self.lo.floor(), hi: self.hi.floor() },
            | Intrinsic::Ceil => Self { lo: self.lo.ceil(), hi: self.hi.ceil() },
            | Intrinsic::Asin if self.lo >= -1.0 && self.hi <= 1.0 => self.monotone(f64::asin, true),
            | Intrinsic::Acos if self.lo >= -1.0 && self.hi <= 1.0 => self.monotone(f64::acos, false),
            | Intrinsic::Acosh if self.lo >= 1.0 => self.monotone(f64::acosh, true),
            | Intrinsic::Atanh if self.lo > -1.0 && self.hi < 1.0 => self.monotone(f64::atanh, true),
            | Intrinsic::Abs => {
                if self.lo >= 0.0 {
                    self
                } else if self.hi <= 0.0 {
                    self.neg()
                } else {
                    Self { lo: 0.0, hi: (-self.lo).max(self.hi) }
                }
            },
            | Intrinsic::Cosh => {
                let (a, b) = (self.lo.cosh(), self.hi.cosh());
                let lo = if self.lo <= 0.0 && self.hi >= 0.0 { 1.0 } else { a.min(b) };
                Self::outward(lo, a.max(b))
            },
            | Intrinsic::Sin | Intrinsic::Cos => {
                if !(self.hi - self.lo).is_finite() || self.hi - self.lo >= 2.0 * PI {
                    return Self { lo: -1.0, hi: 1.0 };
                }
                let g = if f == Intrinsic::Sin { f64::sin } else { f64::cos };
                let (a, b) = (g(self.lo), g(self.hi));
                let (max_at, min_at) = if f == Intrinsic::Sin { (FRAC_PI_2, -FRAC_PI_2) } else { (0.0, PI) };
                let hi = if self.contains_point_of(max_at, 2.0 * PI) { 1.0 } else { a.max(b) };
                let lo = if self.contains_point_of(min_at, 2.0 * PI) { -1.0 } else { a.min(b) };
                Self::outward(lo, hi).clamp(-1.0, 1.0)
            },
            | _ => Self::ENTIRE,
        }
    }

    fn atan2(
        self,
        x: Self,
    ) -> Self {
        let _ = x;
        Self { lo: -std::f64::consts::PI, hi: std::f64::consts::PI }
    }

    fn compare(
        self,
        c: Cmp,
        b: Self,
    ) -> Self {
        let certain = match c {
            | Cmp::Lt => (self.hi < b.lo, self.lo >= b.hi),
            | Cmp::Le => (self.hi <= b.lo, self.lo > b.hi),
            | Cmp::Gt => (self.lo > b.hi, self.hi <= b.lo),
            | Cmp::Ge => (self.lo >= b.hi, self.hi < b.lo),
            | Cmp::Eq => (false, self.hi < b.lo || self.lo > b.hi),
            | Cmp::Ne => (self.hi < b.lo || self.lo > b.hi, false),
        };
        match certain {
            | (true, _) => Self::point(1.0),
            | (_, true) => Self::point(0.0),
            | _ => Self { lo: 0.0, hi: 1.0 },
        }
    }

    fn select(
        self,
        a: Self,
        b: Self,
    ) -> Self {
        if !self.contains(0.0) {
            a
        } else if self.lo == 0.0 && self.hi == 0.0 {
            b
        } else {
            Self { lo: a.lo.min(b.lo), hi: a.hi.max(b.hi) }
        }
    }

    fn call(
        f: EvalFn,
        args: &[Self],
    ) -> Self {
        // Degenerate arguments: the value, widened; otherwise no information.
        if args.iter().all(Self::is_point) {
            let point: Vec<f64> = args.iter().map(|a| a.lo).collect();
            let v = f(&point);
            return Self::outward(v, v);
        }
        Self::ENTIRE
    }
}

impl FromExtern for Interval {
    fn from_extern1(
        f: extern "C" fn(f64) -> f64,
        a: Self,
    ) -> Self {
        if a.is_point() {
            let v = f(a.lo);
            return Self::outward(v, v);
        }
        Self::ENTIRE
    }

    fn from_extern2(
        f: extern "C" fn(f64, f64) -> f64,
        a: Self,
        b: Self,
    ) -> Self {
        if a.is_point() && b.is_point() {
            let v = f(a.lo, b.lo);
            return Self::outward(v, v);
        }
        Self::ENTIRE
    }
}

impl Interval {
    /// Whether the interval is a single point.
    #[must_use]
    pub fn is_point(&self) -> bool {
        self.lo.total_cmp(&self.hi).is_eq()
    }

    const fn widen(self) -> Self {
        Self::outward(self.lo, self.hi)
    }

    const fn clamp_low(
        self,
        lo: f64,
    ) -> Self {
        Self { lo: self.lo.max(lo), hi: self.hi }
    }

    const fn clamp(
        self,
        lo: f64,
        hi: f64,
    ) -> Self {
        Self { lo: self.lo.max(lo), hi: self.hi.min(hi) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Engine;
    use crate::rules;

    fn tape(src: &str) -> Tape {
        let mut g = Graph::new();
        assert!(Engine::install(&mut g, &rules::standard()).is_ok());
        let root = g.parse(src).unwrap_or_else(|e| panic!("{e}"));
        let inputs = [g.interner_mut().symbol("x"), g.interner_mut().symbol("y")];
        lower(&g, &[root], &inputs, &[]).unwrap_or_else(|e| panic!("{e}")).optimise()
    }

    #[test]
    fn canonical_forms_become_native_operations() {
        let t = tape("x - y");
        assert!(t.insts.iter().any(|i| matches!(i, Inst::Sub(..))) && !t.insts.iter().any(|i| matches!(i, Inst::Mul(..))));
        let t = tape("x / y");
        assert!(t.insts.iter().any(|i| matches!(i, Inst::Div(..))) && !t.insts.iter().any(|i| matches!(i, Inst::PowF(..))));
        let t = tape("x^(1/2) + y^(3/2) + x^5");
        assert_eq!(t.insts.iter().filter(|i| matches!(i, Inst::Unary(Intrinsic::Sqrt, _))).count(), 2);
        assert!(!t.insts.iter().any(|i| matches!(i, Inst::PowF(..))));
        // Constants fold (`x + 0.0` stays: it maps -0.0 to +0.0).
        let t = tape("2^10*x + sin(0)");
        assert!(!t.insts.iter().any(|i| matches!(i, Inst::Unary(..) | Inst::PowI(..))), "{:?}", t.insts);
    }

    #[test]
    fn gradients_by_forward_mode() {
        let t = tape("x^3*y + sin(x*y) + exp(y)/x");
        let (x, y) = (1.3, -0.4);
        let g = t.gradient(0, &[x, y], &[]);
        let dx = 3.0 * x * x * y + y * (x * y).cos() - y.exp() / (x * x);
        let dy = x.powi(3) + x * (x * y).cos() + y.exp() / x;
        assert!((g[0] - dx).abs() < 1e-12 && (g[1] - dy).abs() < 1e-12, "{g:?}");
    }

    #[test]
    fn special_functions_lower_to_direct_calls() {
        let mut g = Graph::new();
        assert!(crate::graph::Engine::install(&mut g, &crate::rules::standard()).is_ok());
        let root = g.parse("gamma(x) + erf(x) + besselj(1, x) + floor(x)").unwrap_or_else(|e| panic!("{e}"));
        let x = g.interner_mut().symbol("x");
        let tape = lower(&g, &[root], &[x], &[]).unwrap_or_else(|e| panic!("{e}")).optimise();
        assert!(!tape.insts.iter().any(|i| matches!(i, Inst::Call(..))), "{:?}", tape.insts);
        let mut out = [0.0];
        tape.eval(&[2.5], &[], &mut out);
        let want = 1.329_340_388_179_137 + 0.999_593_047_982_555 + 0.497_094_102_464_274_4 + 2.0;
        assert!((out[0] - want).abs() < 1e-12, "{}", out[0]);
    }

    #[test]
    fn interval_enclosures_are_rigorous() {
        let t = tape("x^2 - 2*x*y + sin(x) + exp(-y^2)");
        let boxed = [Interval { lo: 0.9, hi: 1.1 }, Interval { lo: -0.2, hi: 0.3 }];
        let enclosure = t.enclose(&boxed, &[]);
        let e = enclosure[0];
        let mut out = [0.0];
        for i in 0..=20 {
            for j in 0..=20 {
                let (x, y) = (0.9 + 0.01 * f64::from(i), -0.2 + 0.025 * f64::from(j));
                t.eval(&[x, y], &[], &mut out);
                assert!(e.contains(out[0]), "{} not in {e:?}", out[0]);
            }
        }
        // A sign decided rigorously: x^2 + 1 > 0 on any box.
        let positive = tape("x^2 + 1").enclose(&[Interval { lo: -5.0, hi: 3.0 }, Interval::point(0.0)], &[]);
        assert!(positive[0].lo > 0.0);
    }
}
