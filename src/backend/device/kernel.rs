//! The backend-neutral kernel IR.
//!
//! A [`KernelIr`] is a [`Tape`] with exactly the instructions a data-parallel
//! device can run: arithmetic, the elementary [`Intrinsic`]s, comparisons,
//! `select`, and `gather` (a read of another input column of the same row).
//! The tape instructions that call back into the host (`Call` of an
//! operator's `EvalFn`, `Extern1`/`Extern2` C functions) have no device
//! equivalent and are rejected by [`KernelIr::from_tape`] with
//! [`DeviceError::Unsupported`]; a caller that must run such a program keeps it
//! on the host (see [`super::DeviceBackend`]).
//!
//! Semantics are those of the tape over `f64`: a comparison yields `1.0` or
//! `0.0`, `select(c, a, b)` takes `a` unless `c == 0` (NaN counts as true),
//! an out-of-range `gather` reads NaN. A device that computes in `f32` (see
//! [`super::Dtype`]) approximates these; the emitters document where.
//!
//! Programs are *data-parallel*: one instance of the program per row. Row `r`
//! sees `Input(i)` = element `r` of input column `i`; parameters are the same
//! for every row.

use std::collections::HashMap;
use std::hash::Hash;
use std::hash::Hasher;

use super::DeviceError;
use crate::backend::Intrinsic;
use crate::backend::tape::Cmp;
use crate::backend::tape::Inst;
use crate::backend::tape::Reg;
use crate::backend::tape::Tape;
use crate::backend::tape::int_pow;

/// One instruction; register `r` is the value of `ops[r]`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum KernelOp {
    /// A constant.
    Const(f64),
    /// Element `row` of input column `i`.
    Input(u32),
    /// Parameter `i` (uniform across rows).
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
    /// `a^n`, `n` an integer (square-and-multiply).
    PowI(Reg, i32),
    /// `a^b`.
    PowF(Reg, Reg),
    /// An elementary function.
    Unary(Intrinsic, Reg),
    /// `atan2(a, b)`.
    Atan2(Reg, Reg),
    /// A comparison producing `1.0` or `0.0`.
    Cmp(Cmp, Reg, Reg),
    /// `if c != 0 { a } else { b }` (NaN selects `a`).
    Select(Reg, Reg, Reg),
    /// `input[base + trunc(index)]` of the same row; NaN when out of range.
    Gather(u32, Reg),
}

impl KernelOp {
    /// The operand registers (at most three) and how many are used.
    #[must_use]
    pub const fn operands(&self) -> ([Reg; 3], usize) {
        match *self {
            | Self::Const(_) | Self::Input(_) | Self::Param(_) => ([0; 3], 0),
            | Self::Neg(a) | Self::PowI(a, _) | Self::Unary(_, a) | Self::Gather(_, a) => {
                ([a, 0, 0], 1)
            },
            | Self::Add(a, b)
            | Self::Sub(a, b)
            | Self::Mul(a, b)
            | Self::Div(a, b)
            | Self::PowF(a, b)
            | Self::Atan2(a, b)
            | Self::Cmp(_, a, b) => ([a, b, 0], 2),
            | Self::Select(c, a, b) => ([c, a, b], 3),
        }
    }
}

/// A data-parallel program.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct KernelIr {
    /// Instructions in dependency order.
    pub ops: Vec<KernelOp>,
    /// Number of input columns.
    pub inputs: usize,
    /// Number of uniform parameters.
    pub params: usize,
    /// Output registers; output `k` is written to column `k`.
    pub outputs: Vec<Reg>,
}

impl KernelIr {
    /// Lowers a tape.
    ///
    /// # Errors
    /// [`DeviceError::Unsupported`] when the tape contains a host call
    /// (`Call`, `Extern1`, `Extern2`); [`DeviceError::InvalidKernel`] when it
    /// is malformed.
    pub fn from_tape(tape: &Tape) -> Result<Self, DeviceError> {
        let mut ops = Vec::with_capacity(tape.insts.len());
        for inst in &tape.insts {
            ops.push(match *inst {
                | Inst::Const(c) => KernelOp::Const(c),
                | Inst::Input(i) => KernelOp::Input(i),
                | Inst::Param(i) => KernelOp::Param(i),
                | Inst::Add(a, b) => KernelOp::Add(a, b),
                | Inst::Sub(a, b) => KernelOp::Sub(a, b),
                | Inst::Mul(a, b) => KernelOp::Mul(a, b),
                | Inst::Div(a, b) => KernelOp::Div(a, b),
                | Inst::Neg(a) => KernelOp::Neg(a),
                | Inst::PowI(a, n) => KernelOp::PowI(a, n),
                | Inst::PowF(a, b) => KernelOp::PowF(a, b),
                | Inst::Unary(f, a) => KernelOp::Unary(f, a),
                | Inst::Atan2(a, b) => KernelOp::Atan2(a, b),
                | Inst::Cmp(c, a, b) => KernelOp::Cmp(c, a, b),
                | Inst::Select(c, a, b) => KernelOp::Select(c, a, b),
                | Inst::Gather(base, index) => KernelOp::Gather(base, index),
                | Inst::Call(..) | Inst::Extern1(..) | Inst::Extern2(..) => {
                    return Err(DeviceError::Unsupported(
                        "host calls (operator EvalFn / extern C functions) have no device \
                         lowering"
                            .into(),
                    ));
                },
            });
        }
        let ir = Self {
            ops,
            inputs: tape.inputs,
            params: tape.params,
            outputs: tape.outputs.clone(),
        };
        ir.validate()?;
        Ok(ir)
    }

    /// The same program as a tape (used by the reference checks).
    #[must_use]
    pub fn to_tape(&self) -> Tape {
        let insts = self
            .ops
            .iter()
            .map(|op| {
                match *op {
                    | KernelOp::Const(c) => Inst::Const(c),
                    | KernelOp::Input(i) => Inst::Input(i),
                    | KernelOp::Param(i) => Inst::Param(i),
                    | KernelOp::Add(a, b) => Inst::Add(a, b),
                    | KernelOp::Sub(a, b) => Inst::Sub(a, b),
                    | KernelOp::Mul(a, b) => Inst::Mul(a, b),
                    | KernelOp::Div(a, b) => Inst::Div(a, b),
                    | KernelOp::Neg(a) => Inst::Neg(a),
                    | KernelOp::PowI(a, n) => Inst::PowI(a, n),
                    | KernelOp::PowF(a, b) => Inst::PowF(a, b),
                    | KernelOp::Unary(f, a) => Inst::Unary(f, a),
                    | KernelOp::Atan2(a, b) => Inst::Atan2(a, b),
                    | KernelOp::Cmp(c, a, b) => Inst::Cmp(c, a, b),
                    | KernelOp::Select(c, a, b) => Inst::Select(c, a, b),
                    | KernelOp::Gather(base, i) => Inst::Gather(base, i),
                }
            })
            .collect();
        Tape {
            insts,
            args: Vec::new(),
            inputs: self.inputs,
            params: self.params,
            outputs: self.outputs.clone(),
        }
    }

    /// Checks SSA order and index ranges. Every device calls this before it
    /// runs or emits a program, so a malformed program cannot index outside a
    /// register file or a buffer.
    ///
    /// # Errors
    /// [`DeviceError::InvalidKernel`].
    pub fn validate(&self) -> Result<(), DeviceError> {
        let bad = |why: String| Err(DeviceError::InvalidKernel(why));
        if self.outputs.is_empty() {
            return bad("a kernel needs at least one output".into());
        }
        for (at, op) in self.ops.iter().enumerate() {
            let (operands, n) = op.operands();
            if operands.iter().take(n).any(|&r| r as usize >= at) {
                return bad(format!("op {at} reads a register defined later"));
            }
            match *op {
                | KernelOp::Input(i) if i as usize >= self.inputs => {
                    return bad(format!("op {at}: input {i} of {}", self.inputs));
                },
                | KernelOp::Param(i) if i as usize >= self.params => {
                    return bad(format!("op {at}: parameter {i} of {}", self.params));
                },
                | KernelOp::Gather(base, _) if base as usize > self.inputs => {
                    return bad(format!("op {at}: gather base {base} of {}", self.inputs));
                },
                | _ => {},
            }
        }
        if let Some(&r) = self.outputs.iter().find(|&&r| r as usize >= self.ops.len()) {
            return bad(format!("output register {r} is out of range"));
        }
        Ok(())
    }

    /// A structural hash (constants by bit pattern): the key under which a
    /// device caches its compiled form.
    #[must_use]
    pub fn structural_hash(&self) -> u64 {
        let mut h = std::collections::hash_map::DefaultHasher::new();
        self.inputs.hash(&mut h);
        self.params.hash(&mut h);
        self.outputs.hash(&mut h);
        for op in &self.ops {
            match *op {
                | KernelOp::Const(c) => (0u8, c.to_bits()).hash(&mut h),
                | KernelOp::Input(i) => (1u8, i).hash(&mut h),
                | KernelOp::Param(i) => (2u8, i).hash(&mut h),
                | KernelOp::Add(a, b) => (3u8, a, b).hash(&mut h),
                | KernelOp::Sub(a, b) => (4u8, a, b).hash(&mut h),
                | KernelOp::Mul(a, b) => (5u8, a, b).hash(&mut h),
                | KernelOp::Div(a, b) => (6u8, a, b).hash(&mut h),
                | KernelOp::Neg(a) => (7u8, a).hash(&mut h),
                | KernelOp::PowI(a, n) => (8u8, a, n).hash(&mut h),
                | KernelOp::PowF(a, b) => (9u8, a, b).hash(&mut h),
                | KernelOp::Unary(f, a) => (10u8, f, a).hash(&mut h),
                | KernelOp::Atan2(a, b) => (11u8, a, b).hash(&mut h),
                | KernelOp::Cmp(c, a, b) => (12u8, c, a, b).hash(&mut h),
                | KernelOp::Select(c, a, b) => (13u8, c, a, b).hash(&mut h),
                | KernelOp::Gather(base, i) => (14u8, base, i).hash(&mut h),
            }
        }
        h.finish()
    }

    /// Reference evaluation of one row (the semantics every device must
    /// reproduce, up to its precision). Slow; for tests and fallbacks.
    #[must_use]
    pub fn eval_row(
        &self,
        inputs: &[f64],
        params: &[f64],
    ) -> Vec<f64> {
        let mut v: Vec<f64> = Vec::with_capacity(self.ops.len());
        for op in &self.ops {
            let g = |r: Reg| v.get(r as usize).copied().unwrap_or(f64::NAN);
            let x = match *op {
                | KernelOp::Const(c) => c,
                | KernelOp::Input(i) => inputs.get(i as usize).copied().unwrap_or(f64::NAN),
                | KernelOp::Param(i) => params.get(i as usize).copied().unwrap_or(f64::NAN),
                | KernelOp::Add(a, b) => g(a) + g(b),
                | KernelOp::Sub(a, b) => g(a) - g(b),
                | KernelOp::Mul(a, b) => g(a) * g(b),
                | KernelOp::Div(a, b) => g(a) / g(b),
                | KernelOp::Neg(a) => -g(a),
                | KernelOp::PowI(a, n) => int_pow(g(a), n),
                | KernelOp::PowF(a, b) => g(a).powf(g(b)),
                | KernelOp::Unary(f, a) => f.apply(g(a)),
                | KernelOp::Atan2(a, b) => g(a).atan2(g(b)),
                | KernelOp::Cmp(c, a, b) => cmp_value(c, g(a), g(b)),
                | KernelOp::Select(c, a, b) => {
                    if select_is_true(g(c)) {
                        g(a)
                    } else {
                        g(b)
                    }
                },
                | KernelOp::Gather(base, i) => {
                    gather_slot(base, g(i), self.inputs)
                        .and_then(|s| inputs.get(s).copied())
                        .unwrap_or(f64::NAN)
                },
            };
            v.push(x);
        }
        self.outputs
            .iter()
            .map(|&r| v.get(r as usize).copied().unwrap_or(f64::NAN))
            .collect()
    }
}

/// `select` takes its first branch unless the condition is exactly zero.
#[inline]
#[must_use]
pub fn select_is_true(c: f64) -> bool {
    c.is_nan() || c != 0.0
}

/// A comparison as `1.0` / `0.0` (`!=` is true for NaN operands).
#[inline]
#[must_use]
pub fn cmp_value(
    c: Cmp,
    a: f64,
    b: f64,
) -> f64 {
    let r = match c {
        | Cmp::Lt => a < b,
        | Cmp::Le => a <= b,
        | Cmp::Gt => a > b,
        | Cmp::Ge => a >= b,
        | Cmp::Eq => a == b,
        | Cmp::Ne => a != b,
    };
    f64::from(u8::from(r))
}

/// The input column a gather reads, or `None` when the index is not finite,
/// is negative, or points past the declared inputs.
#[inline]
#[must_use]
pub fn gather_slot(
    base: u32,
    offset: f64,
    inputs: usize,
) -> Option<usize> {
    if !offset.is_finite() {
        return None;
    }
    let offset = offset.trunc();
    if !(0.0..=1.0e18).contains(&offset) {
        return None;
    }
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let slot = (base as usize).checked_add(offset as usize)?;
    (slot < inputs).then_some(slot)
}

/// Builds a [`KernelIr`] with leaf deduplication.
#[derive(Debug, Default)]
pub struct KernelBuilder {
    ops: Vec<KernelOp>,
    leaves: HashMap<(u8, u64), Reg>,
    inputs: usize,
    params: usize,
}

impl KernelBuilder {
    /// A builder for a kernel with the given input and parameter counts.
    #[must_use]
    pub fn new(
        inputs: usize,
        params: usize,
    ) -> Self {
        Self {
            ops: Vec::new(),
            leaves: HashMap::new(),
            inputs,
            params,
        }
    }

    /// Appends an operation and returns its register (leaves are shared).
    pub fn push(
        &mut self,
        op: KernelOp,
    ) -> Reg {
        let key = match op {
            | KernelOp::Const(c) => Some((0u8, c.to_bits())),
            | KernelOp::Input(i) => Some((1, u64::from(i))),
            | KernelOp::Param(i) => Some((2, u64::from(i))),
            | _ => None,
        };
        if let Some(r) = key.and_then(|k| self.leaves.get(&k).copied()) {
            return r;
        }
        let r = Reg::try_from(self.ops.len()).unwrap_or(Reg::MAX);
        self.ops.push(op);
        if let Some(k) = key {
            self.leaves.insert(k, r);
        }
        r
    }

    /// Splices `other` into this program with its input `i` bound to
    /// `inputs[i]` and its parameter `j` to `params[j]`; returns its outputs.
    /// `other` must not use `gather` (it addresses the caller's columns).
    pub fn splice(
        &mut self,
        other: &KernelIr,
        inputs: &[Reg],
        params: &[Reg],
    ) -> Vec<Reg> {
        let mut map: Vec<Reg> = Vec::with_capacity(other.ops.len());
        for op in &other.ops {
            let m = |r: Reg| map.get(r as usize).copied().unwrap_or(0);
            let mapped = match *op {
                | KernelOp::Const(c) => KernelOp::Const(c),
                | KernelOp::Input(i) => {
                    map.push(inputs.get(i as usize).copied().unwrap_or(0));
                    continue;
                },
                | KernelOp::Param(j) => {
                    map.push(params.get(j as usize).copied().unwrap_or(0));
                    continue;
                },
                | KernelOp::Add(a, b) => KernelOp::Add(m(a), m(b)),
                | KernelOp::Sub(a, b) => KernelOp::Sub(m(a), m(b)),
                | KernelOp::Mul(a, b) => KernelOp::Mul(m(a), m(b)),
                | KernelOp::Div(a, b) => KernelOp::Div(m(a), m(b)),
                | KernelOp::Neg(a) => KernelOp::Neg(m(a)),
                | KernelOp::PowI(a, n) => KernelOp::PowI(m(a), n),
                | KernelOp::PowF(a, b) => KernelOp::PowF(m(a), m(b)),
                | KernelOp::Unary(f, a) => KernelOp::Unary(f, m(a)),
                | KernelOp::Atan2(a, b) => KernelOp::Atan2(m(a), m(b)),
                | KernelOp::Cmp(c, a, b) => KernelOp::Cmp(c, m(a), m(b)),
                | KernelOp::Select(c, a, b) => KernelOp::Select(m(c), m(a), m(b)),
                | KernelOp::Gather(..) => KernelOp::Const(f64::NAN),
            };
            let r = self.push(mapped);
            map.push(r);
        }
        other
            .outputs
            .iter()
            .map(|&r| map.get(r as usize).copied().unwrap_or(0))
            .collect()
    }

    /// Finishes with the given outputs.
    ///
    /// # Errors
    /// [`DeviceError::InvalidKernel`] when the result does not validate.
    pub fn finish(
        self,
        outputs: Vec<Reg>,
    ) -> Result<KernelIr, DeviceError> {
        let ir = KernelIr {
            ops: self.ops,
            inputs: self.inputs,
            params: self.params,
            outputs,
        };
        ir.validate()?;
        Ok(ir)
    }
}

/// One classical Runge-Kutta 4 step of `y' = f([t,] y)` as a kernel.
///
/// `rhs` has `n_state` outputs and either `n_state` inputs (autonomous) or
/// `n_state + 1` with the time as input 0. The step kernel has `n_state`
/// inputs (the state columns of every trajectory), the parameters of `rhs`
/// followed by two more, `t` at index `rhs.params` and the step `h` at
/// `rhs.params + 1`, and `n_state` outputs (the next state). Launching it
/// `steps` times with a changing `t` integrates every row, that is every
/// initial condition, with the same program: the ODE-ensemble pattern.
///
/// # Errors
/// [`DeviceError::InvalidKernel`] when the shape of `rhs` does not match, or
/// when `rhs` uses `gather` (which addresses columns and cannot be composed).
pub fn rk4_step(
    rhs: &KernelIr,
    n_state: usize,
) -> Result<KernelIr, DeviceError> {
    let timed = if rhs.inputs == n_state + 1 {
        true
    } else if rhs.inputs == n_state {
        false
    } else {
        return Err(DeviceError::InvalidKernel(format!(
            "rhs has {} inputs, expected {n_state} or {}",
            rhs.inputs,
            n_state + 1
        )));
    };
    if rhs.outputs.len() != n_state || n_state == 0 {
        return Err(DeviceError::InvalidKernel(
            "rhs must have one output per state variable".into(),
        ));
    }
    if rhs.ops.iter().any(|o| matches!(o, KernelOp::Gather(..))) {
        return Err(DeviceError::InvalidKernel(
            "an ODE right-hand side cannot use gather".into(),
        ));
    }
    let p = rhs.params;
    let idx = |i: usize| u32::try_from(i).unwrap_or(u32::MAX);
    let mut b = KernelBuilder::new(n_state, p + 2);
    let params: Vec<Reg> = (0..p).map(|j| b.push(KernelOp::Param(idx(j)))).collect();
    let t = b.push(KernelOp::Param(idx(p)));
    let h = b.push(KernelOp::Param(idx(p + 1)));
    let y: Vec<Reg> = (0..n_state)
        .map(|i| b.push(KernelOp::Input(idx(i))))
        .collect();
    let half = b.push(KernelOp::Const(0.5));
    let two = b.push(KernelOp::Const(2.0));
    let sixth = b.push(KernelOp::Const(1.0 / 6.0));
    let hh = b.push(KernelOp::Mul(h, half));
    let t_mid = b.push(KernelOp::Add(t, hh));
    let t_end = b.push(KernelOp::Add(t, h));

    let eval = |b: &mut KernelBuilder, time: Reg, state: &[Reg]| -> Vec<Reg> {
        let mut ins = Vec::with_capacity(rhs.inputs);
        if timed {
            ins.push(time);
        }
        ins.extend_from_slice(state);
        b.splice(rhs, &ins, &params)
    };
    let stage = |b: &mut KernelBuilder, k: &[Reg], factor: Reg| -> Vec<Reg> {
        y.iter()
            .zip(k)
            .map(|(&yi, &ki)| {
                let d = b.push(KernelOp::Mul(factor, ki));
                b.push(KernelOp::Add(yi, d))
            })
            .collect()
    };
    let k1 = eval(&mut b, t, &y);
    let y2 = stage(&mut b, &k1, hh);
    let k2 = eval(&mut b, t_mid, &y2);
    let y3 = stage(&mut b, &k2, hh);
    let k3 = eval(&mut b, t_mid, &y3);
    let y4 = stage(&mut b, &k3, h);
    let k4 = eval(&mut b, t_end, &y4);
    let mut out = Vec::with_capacity(n_state);
    for i in 0..n_state {
        let g = |k: &[Reg]| k.get(i).copied().unwrap_or(0);
        let k2x = b.push(KernelOp::Mul(two, g(&k2)));
        let k3x = b.push(KernelOp::Mul(two, g(&k3)));
        let s1 = b.push(KernelOp::Add(g(&k1), k2x));
        let s2 = b.push(KernelOp::Add(s1, k3x));
        let s3 = b.push(KernelOp::Add(s2, g(&k4)));
        let hs = b.push(KernelOp::Mul(h, sixth));
        let inc = b.push(KernelOp::Mul(hs, s3));
        out.push(b.push(KernelOp::Add(y.get(i).copied().unwrap_or(0), inc)));
    }
    b.finish(out)
}
