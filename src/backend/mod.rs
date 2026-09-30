//! # Execution backends
//!
//! Once a request has been reduced to a closed form, evaluating it many
//! times — over the rows of a data set, inside an ODE stepper, at every
//! cell of a simulation grid — is a different problem from finding it. A
//! [`Backend`] turns a closed-form term into a [`Compiled`] function of its
//! free symbols.
//!
//! The contract is deliberately small so that very different executors fit
//! behind it:
//!
//! * [`Interpreter`] (this module) flattens the term's DAG into a tape and
//!   walks it; always available, no code generation;
//! * a JIT backend lowers the same DAG to machine code. Each operator's
//!   [`EvalFn`](crate::graph::op::EvalFn) is a plain `fn` pointer it can
//!   emit a direct call to, and richer lowerings can be attached to
//!   operators as typed attributes
//!   ([`OpTable::set_attr`](crate::graph::OpTable::set_attr)) without the
//!   kernel knowing about them.
//!
//! Numeric reduction kernels and the simulation steppers only see
//! `dyn Compiled`, so switching executors never touches them.

use std::collections::HashMap;
use std::fmt;

use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::SymbolId;
use crate::graph::op::EvalFn;

/// Why a term could not be compiled.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BackendError {
    /// The term depends on a symbol that is not among the inputs.
    UnboundSymbol(String),
    /// The term uses an operator without scalar semantics (for example an
    /// unevaluated heavy operator).
    NoSemantics(String),
    /// The term contains a non-numeric literal.
    NotNumeric,
}

impl fmt::Display for BackendError {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        match self {
            | Self::UnboundSymbol(name) => write!(
                f,
                "symbol `{name}` is not an input of the compiled function"
            ),
            | Self::NoSemantics(op) => write!(f, "operator `{op}` has no numeric semantics"),
            | Self::NotNumeric => write!(f, "the term contains a non-numeric literal"),
        }
    }
}

impl std::error::Error for BackendError {}

/// A term compiled to a function `f64^n -> f64`.
pub trait Compiled: Send + Sync {
    /// Number of inputs.
    fn arity(&self) -> usize;

    /// Evaluates at one point. `args` must hold [`Compiled::arity`] values;
    /// missing ones are read as NaN.
    fn call(
        &self,
        args: &[f64],
    ) -> f64;

    /// Evaluates at many points: `columns[i][row]` is input `i` of `row`.
    /// Writes one result per element of `out`.
    fn call_batch(
        &self,
        columns: &[&[f64]],
        out: &mut [f64],
    ) {
        let mut args = vec![f64::NAN; self.arity()];
        for (row, slot) in out.iter_mut().enumerate() {
            for (arg, column) in args.iter_mut().zip(columns) {
                *arg = column.get(row).copied().unwrap_or(f64::NAN);
            }
            *slot = self.call(&args);
        }
    }
}

/// Turns closed-form terms into executable functions.
pub trait Backend {
    /// Name for diagnostics.
    fn name(&self) -> &'static str;

    /// Compiles the concrete term `root` as a function of `inputs`, in that
    /// order.
    ///
    /// # Errors
    /// See [`BackendError`].
    fn compile(
        &self,
        graph: &Graph,
        root: NodeId,
        inputs: &[SymbolId],
    ) -> Result<Box<dyn Compiled>, BackendError>;
}

#[derive(Clone, Copy)]
enum Instr {
    Const(f64),
    Input(u32),
    Call { f: EvalFn, start: u32, len: u32 },
}

/// The reference backend: a register tape in topological order.
///
/// Every distinct subterm becomes one instruction, so shared subterms are
/// computed once per call, exactly as a JIT would keep them in a register.
#[derive(Copy, Clone, Debug, Default)]
pub struct Interpreter;

struct Tape {
    arity: usize,
    instrs: Vec<Instr>,
    args: Vec<u32>,
}

impl Backend for Interpreter {
    fn name(&self) -> &'static str {
        "interpreter"
    }

    fn compile(
        &self,
        graph: &Graph,
        root: NodeId,
        inputs: &[SymbolId],
    ) -> Result<Box<dyn Compiled>, BackendError> {
        let mut tape = Tape {
            arity: inputs.len(),
            instrs: Vec::new(),
            args: Vec::new(),
        };
        let mut slot: HashMap<NodeId, u32> = HashMap::new();
        let mut stack = vec![(root, false)];
        while let Some((node, expanded)) = stack.pop() {
            if slot.contains_key(&node) {
                continue;
            }
            let children = graph.children(node);
            if !expanded && !children.is_empty() {
                stack.push((node, true));
                stack.extend(
                    children
                        .iter()
                        .filter(|c| !slot.contains_key(c))
                        .map(|&c| (c, false)),
                );
                continue;
            }
            let instr = if let Some(n) = graph.as_number(node) {
                Instr::Const(n.to_f64())
            } else if let Some(symbol) = graph.as_symbol(node) {
                let index = inputs.iter().position(|&s| s == symbol).ok_or_else(|| {
                    BackendError::UnboundSymbol(graph.interner().symbol_name(symbol).to_owned())
                })?;
                Instr::Input(u32::try_from(index).unwrap_or(u32::MAX))
            } else if graph.payload(node).is_some() {
                return Err(BackendError::NotNumeric);
            } else {
                let desc = graph.ops().get(graph.op(node));
                let f = desc
                    .eval
                    .ok_or_else(|| BackendError::NoSemantics(desc.name.to_string()))?;
                let start = u32::try_from(tape.args.len()).unwrap_or(u32::MAX);
                for child in children {
                    tape.args.push(slot.get(child).copied().unwrap_or(u32::MAX));
                }
                Instr::Call {
                    f,
                    start,
                    len: u32::try_from(children.len()).unwrap_or(u32::MAX),
                }
            };
            slot.insert(node, u32::try_from(tape.instrs.len()).unwrap_or(u32::MAX));
            tape.instrs.push(instr);
        }
        Ok(Box::new(tape))
    }
}

impl Tape {
    fn run(
        &self,
        args: &[f64],
        values: &mut Vec<f64>,
        scratch: &mut Vec<f64>,
    ) -> f64 {
        values.clear();
        for instr in &self.instrs {
            let value = match *instr {
                | Instr::Const(c) => c,
                | Instr::Input(i) => args.get(i as usize).copied().unwrap_or(f64::NAN),
                | Instr::Call { f, start, len } => {
                    scratch.clear();
                    let range = start as usize..(start as usize).saturating_add(len as usize);
                    for &arg in self.args.get(range).unwrap_or(&[]) {
                        scratch.push(values.get(arg as usize).copied().unwrap_or(f64::NAN));
                    }
                    f(scratch)
                },
            };
            values.push(value);
        }
        values.last().copied().unwrap_or(f64::NAN)
    }
}

impl Compiled for Tape {
    fn arity(&self) -> usize {
        self.arity
    }

    fn call(
        &self,
        args: &[f64],
    ) -> f64 {
        self.run(
            args,
            &mut Vec::with_capacity(self.instrs.len()),
            &mut Vec::new(),
        )
    }

    fn call_batch(
        &self,
        columns: &[&[f64]],
        out: &mut [f64],
    ) {
        let mut values = Vec::with_capacity(self.instrs.len());
        let mut scratch = Vec::new();
        let mut args = vec![f64::NAN; self.arity];
        for (row, slot) in out.iter_mut().enumerate() {
            for (arg, column) in args.iter_mut().zip(columns) {
                *arg = column.get(row).copied().unwrap_or(f64::NAN);
            }
            *slot = self.run(&args, &mut values, &mut scratch);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Engine;
    use crate::graph::Env;
    use crate::rules;

    fn setup(src: &str) -> (Graph, NodeId, Vec<SymbolId>) {
        let mut g = Graph::new();
        assert!(Engine::install(&mut g, &rules::standard()).is_ok());
        let root = g.parse(src).unwrap_or_else(|e| panic!("{e}"));
        let inputs = vec![g.interner_mut().symbol("x"), g.interner_mut().symbol("y")];
        (g, root, inputs)
    }

    #[test]
    fn agrees_with_the_reference_evaluator() {
        let (g, root, inputs) = setup("sin(x)^2 * exp(-y) + (x - y)/(1 + x^2) + pi");
        let f = Interpreter
            .compile(&g, root, &inputs)
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(f.arity(), 2);
        for (x, y) in [(0.0, 0.0), (1.5, -2.0), (-0.3, 4.0)] {
            let mut env = Env::numeric(0.0);
            env.bind(inputs[0], x);
            env.bind(inputs[1], y);
            let want = g.eval(root, &env).unwrap_or(f64::NAN);
            assert_eq!(f.call(&[x, y]).to_bits(), want.to_bits());
        }
    }

    #[test]
    fn batch_evaluation() {
        let (g, root, inputs) = setup("x^2 + y");
        let f = Interpreter
            .compile(&g, root, &inputs)
            .unwrap_or_else(|e| panic!("{e}"));
        let xs = [1.0, 2.0, 3.0];
        let ys = [10.0, 20.0, 30.0];
        let mut out = [0.0; 3];
        f.call_batch(&[&xs, &ys], &mut out);
        assert_eq!(out, [11.0, 24.0, 39.0]);
    }

    #[test]
    fn shared_subterms_become_one_instruction() {
        let (mut g, _, inputs) = setup("x");
        let mut node = g.sym("x");
        for _ in 0..30 {
            node = g.node(crate::graph::op::core::POW, &[node, node]);
        }
        // 2^30 paths, 31 distinct subterms: must compile and run instantly.
        let f = Interpreter
            .compile(&g, node, &inputs)
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(f.call(&[1.0, 0.0]), 1.0);
    }

    #[test]
    fn errors() {
        let (g, root, inputs) = setup("x + z");
        assert_eq!(
            Interpreter.compile(&g, root, &inputs).err(),
            Some(BackendError::UnboundSymbol("z".to_owned()))
        );
        let (g, root, inputs) = setup("diff(x^2, x)");
        assert_eq!(
            Interpreter.compile(&g, root, &inputs).err(),
            Some(BackendError::NoSemantics("diff".to_owned()))
        );
        assert_eq!(Interpreter.name(), "interpreter");
    }
}
