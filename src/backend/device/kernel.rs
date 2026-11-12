//! The backend-neutral kernel IR.
//!
//! A [`KernelIr`] is a typed SSA program that a data-parallel device runs
//! once per *row*. It has two layers.
//!
//! **The elementwise core** is a [`Tape`] restricted to what a device can
//! run: arithmetic, the elementary [`Intrinsic`]s, comparisons, `select`, and
//! `gather`. The tape instructions that call back into the host (`Call` of an
//! operator's `EvalFn`, `Extern1`/`Extern2` C functions) have no device
//! equivalent and are rejected by [`KernelIr::from_tape`] with
//! [`DeviceError::Unsupported`].
//!
//! **The structured extension** is what turns "evaluate a formula per row"
//! into a general array-programming primitive without adding a second IR:
//!
//! * *index values*: [`KernelOp::Row`], [`KernelOp::Rows`],
//!   [`KernelOp::Len`], integer constants and integer arithmetic. Registers
//!   are typed ([`Ty::Num`] or [`Ty::Idx`]) and never mixed implicitly;
//! * *addressed reads*: [`KernelOp::Load`] reads element `idx` of an input
//!   column. Columns have independent lengths, so a kernel can read a table,
//!   a neighbour, or a CSR structure;
//! * *one fold*: a [`Fold`] runs a body `count` times with the inner index
//!   [`KernelOp::Inner`] and combines values into accumulators with a
//!   [`ReduceOp`]; [`KernelOp::Acc`] reads an accumulator afterwards.
//!
//! With these, a reduction pass, a matrix product, a stencil, a sparse
//! matrix-vector product and an all-pairs interaction are all *kernels*, and
//! a device that implements the one launch entry point
//! ([`Device::launch_cols`](super::Device::launch_cols)) gets every one of
//! them (see [`super::array`]). Nothing in this file knows about a GPU.
//!
//! ## Phases
//!
//! A program with a fold has three phases, determined by data flow rather
//! than by syntax ([`KernelIr::analyse`]): *pre* (does not depend on
//! `Inner` or `Acc`; runs once per row), *body* (depends on `Inner`; runs
//! `count` times) and *post* (depends on `Acc`; runs once after the loop). A
//! register that would need both `Inner` and `Acc` is rejected.
//!
//! ## Semantics
//!
//! Numbers are `f64`: a comparison yields `1.0` or `0.0`, `select(c, a, b)`
//! takes `a` unless `c == 0` (NaN counts as true), an out-of-range `gather`
//! or `Load` reads NaN. Indices are signed integers; division and remainder
//! truncate towards zero and yield `0` for a zero divisor. A device that
//! computes in `f32` and `i32` (see [`super::Dtype`]) approximates these;
//! the emitters document where.

use std::collections::HashMap;
use std::hash::Hash;
use std::hash::Hasher;

use serde::Deserialize;
use serde::Serialize;

use super::DeviceError;
use super::ReduceOp;
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
    /// `if c != 0 { a } else { b }` (NaN selects `a`). `a` and `b` have one
    /// type, numbers or indices, and the result has it too.
    Select(Reg, Reg, Reg),
    /// `input[base + trunc(index)]` of the same row; NaN when out of range.
    Gather(u32, Reg),
    /// The row being computed (index).
    Row,
    /// The number of rows of the launch (index).
    Rows,
    /// The length of input column `i` (index).
    Len(u32),
    /// The inner index of the [`Fold`], `0..count` (index; body only).
    Inner,
    /// An index constant (must fit 32 bits).
    IConst(i64),
    /// Index `a + b`.
    IAdd(Reg, Reg),
    /// Index `a - b`.
    ISub(Reg, Reg),
    /// Index `a * b`.
    IMul(Reg, Reg),
    /// Index `a / b`, truncated; `0` when `b == 0`.
    IDiv(Reg, Reg),
    /// Index `a % b` (sign of `a`); `0` when `b == 0`.
    IMod(Reg, Reg),
    /// The smaller index.
    IMin(Reg, Reg),
    /// The larger index.
    IMax(Reg, Reg),
    /// A comparison of two indices producing the number `1.0` or `0.0`.
    ICmp(Cmp, Reg, Reg),
    /// An index as a number.
    ToNum(Reg),
    /// A number truncated to an index (NaN gives `0`; saturating).
    ToIdx(Reg),
    /// Element `idx` of input column `col`; NaN when out of range.
    Load(u32, Reg),
    /// Accumulator `k` of the [`Fold`], after the loop.
    Acc(u32),
}

/// The type of a register.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Ty {
    /// A floating-point number.
    Num,
    /// A signed integer index.
    Idx,
}

/// When a register is computed (see the module documentation).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Phase {
    /// Once per row, before the fold.
    Pre,
    /// Once per inner index.
    Body,
    /// Once per row, after the fold.
    Post,
}

/// The loop of a kernel: `count` iterations, each combining one value into
/// each accumulator.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Fold {
    /// Number of iterations (an index register of the pre phase; a negative
    /// count runs no iteration).
    pub count: Reg,
    /// Per accumulator: how values combine, and the register combined (a
    /// number of the pre phase or the body). Accumulator `k` starts at
    /// [`ReduceOp::identity`].
    pub accs: Vec<(ReduceOp, Reg)>,
}

/// Types and phases of a validated program.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Analysis {
    /// Type of each register.
    pub ty: Vec<Ty>,
    /// Phase of each register.
    pub phase: Vec<Phase>,
    /// Registers of the pre phase, in order.
    pub pre: Vec<usize>,
    /// Registers of the body, in order.
    pub body: Vec<usize>,
    /// Registers of the post phase, in order.
    pub post: Vec<usize>,
}

impl KernelOp {
    /// The operand registers (at most three) and how many are used.
    #[must_use]
    pub const fn operands(&self) -> ([Reg; 3], usize) {
        match *self {
            | Self::Const(_)
            | Self::Input(_)
            | Self::Param(_)
            | Self::Row
            | Self::Rows
            | Self::Len(_)
            | Self::Inner
            | Self::IConst(_)
            | Self::Acc(_) => ([0; 3], 0),
            | Self::Neg(a)
            | Self::PowI(a, _)
            | Self::Unary(_, a)
            | Self::Gather(_, a)
            | Self::ToNum(a)
            | Self::ToIdx(a)
            | Self::Load(_, a) => ([a, 0, 0], 1),
            | Self::Add(a, b)
            | Self::Sub(a, b)
            | Self::Mul(a, b)
            | Self::Div(a, b)
            | Self::PowF(a, b)
            | Self::Atan2(a, b)
            | Self::Cmp(_, a, b)
            | Self::IAdd(a, b)
            | Self::ISub(a, b)
            | Self::IMul(a, b)
            | Self::IDiv(a, b)
            | Self::IMod(a, b)
            | Self::IMin(a, b)
            | Self::IMax(a, b)
            | Self::ICmp(_, a, b) => ([a, b, 0], 2),
            | Self::Select(c, a, b) => ([c, a, b], 3),
        }
    }

    /// The type of the value this instruction defines (`select` takes the
    /// type of its branches, which [`KernelIr::analyse`] resolves; here it
    /// reports a number).
    #[must_use]
    pub const fn ty(&self) -> Ty {
        match self {
            | Self::Row
            | Self::Rows
            | Self::Len(_)
            | Self::Inner
            | Self::IConst(_)
            | Self::IAdd(..)
            | Self::ISub(..)
            | Self::IMul(..)
            | Self::IDiv(..)
            | Self::IMod(..)
            | Self::IMin(..)
            | Self::IMax(..)
            | Self::ToIdx(_) => Ty::Idx,
            | _ => Ty::Num,
        }
    }

    /// The type each operand must have.
    #[must_use]
    pub const fn operand_ty(&self) -> Ty {
        match self {
            | Self::IAdd(..)
            | Self::ISub(..)
            | Self::IMul(..)
            | Self::IDiv(..)
            | Self::IMod(..)
            | Self::IMin(..)
            | Self::IMax(..)
            | Self::ICmp(..)
            | Self::ToNum(_)
            | Self::Load(..) => Ty::Idx,
            | _ => Ty::Num,
        }
    }

    /// Whether this is an instruction of the elementwise core (a tape
    /// instruction).
    #[must_use]
    pub const fn is_elementwise(&self) -> bool {
        matches!(
            self,
            Self::Const(_)
                | Self::Input(_)
                | Self::Param(_)
                | Self::Add(..)
                | Self::Sub(..)
                | Self::Mul(..)
                | Self::Div(..)
                | Self::Neg(_)
                | Self::PowI(..)
                | Self::PowF(..)
                | Self::Unary(..)
                | Self::Atan2(..)
                | Self::Cmp(..)
                | Self::Select(..)
                | Self::Gather(..)
        )
    }

    /// The same instruction with every operand register renamed by `m`.
    #[must_use]
    pub fn map_operands(
        &self,
        m: impl Fn(Reg) -> Reg,
    ) -> Self {
        match *self {
            | Self::Add(a, b) => Self::Add(m(a), m(b)),
            | Self::Sub(a, b) => Self::Sub(m(a), m(b)),
            | Self::Mul(a, b) => Self::Mul(m(a), m(b)),
            | Self::Div(a, b) => Self::Div(m(a), m(b)),
            | Self::Neg(a) => Self::Neg(m(a)),
            | Self::PowI(a, n) => Self::PowI(m(a), n),
            | Self::PowF(a, b) => Self::PowF(m(a), m(b)),
            | Self::Unary(f, a) => Self::Unary(f, m(a)),
            | Self::Atan2(a, b) => Self::Atan2(m(a), m(b)),
            | Self::Cmp(c, a, b) => Self::Cmp(c, m(a), m(b)),
            | Self::Select(c, a, b) => Self::Select(m(c), m(a), m(b)),
            | Self::Gather(base, a) => Self::Gather(base, m(a)),
            | Self::IAdd(a, b) => Self::IAdd(m(a), m(b)),
            | Self::ISub(a, b) => Self::ISub(m(a), m(b)),
            | Self::IMul(a, b) => Self::IMul(m(a), m(b)),
            | Self::IDiv(a, b) => Self::IDiv(m(a), m(b)),
            | Self::IMod(a, b) => Self::IMod(m(a), m(b)),
            | Self::IMin(a, b) => Self::IMin(m(a), m(b)),
            | Self::IMax(a, b) => Self::IMax(m(a), m(b)),
            | Self::ICmp(c, a, b) => Self::ICmp(c, m(a), m(b)),
            | Self::ToNum(a) => Self::ToNum(m(a)),
            | Self::ToIdx(a) => Self::ToIdx(m(a)),
            | Self::Load(c, a) => Self::Load(c, m(a)),
            | leaf => leaf,
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
    /// The loop, when the program has one.
    pub fold: Option<Fold>,
}

/// The input columns of a launch, as the reference evaluator sees them: one
/// flat buffer and each column's offset and length.
#[derive(Clone, Copy, Debug)]
pub struct Columns<'a> {
    /// All columns, concatenated.
    pub data: &'a [f64],
    /// Offset of each column in `data`.
    pub offsets: &'a [usize],
    /// Length of each column.
    pub lens: &'a [usize],
}

impl Columns<'_> {
    /// Element `index` of column `col`, NaN when either is out of range.
    #[inline]
    #[must_use]
    pub fn get(
        &self,
        col: usize,
        index: i64,
    ) -> f64 {
        let (Some(&off), Some(&len)) = (self.offsets.get(col), self.lens.get(col)) else {
            return f64::NAN;
        };
        match usize::try_from(index) {
            | Ok(i) if i < len => {
                self.data
                    .get(off.saturating_add(i))
                    .copied()
                    .unwrap_or(f64::NAN)
            },
            | _ => f64::NAN,
        }
    }
}

/// Prefix sums of column lengths: the offsets of a column-major buffer.
#[must_use]
pub fn column_offsets(lens: &[usize]) -> Vec<usize> {
    let mut at = 0_usize;
    lens.iter()
        .map(|&len| {
            let here = at;
            at = at.saturating_add(len);
            here
        })
        .collect()
}

/// Reusable register storage for [`KernelIr::eval_at`].
#[derive(Clone, Debug, Default)]
pub struct Scratch {
    num: Vec<f64>,
    idx: Vec<i64>,
    acc: Vec<f64>,
}

/// The value of one register.
enum Val {
    Num(f64),
    Idx(i64),
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
            fold: None,
        };
        ir.validate()?;
        Ok(ir)
    }

    /// Whether the program is in the elementwise core: no fold and only tape
    /// instructions.
    #[must_use]
    pub fn is_elementwise(&self) -> bool {
        self.fold.is_none() && self.ops.iter().all(KernelOp::is_elementwise)
    }

    /// The same program as a tape (used by the reference checks), or `None`
    /// when it uses the structured extension, which a tape cannot express.
    #[must_use]
    pub fn to_tape(&self) -> Option<Tape> {
        if !self.is_elementwise() {
            return None;
        }
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
                    // Unreachable: `is_elementwise` was checked above.
                    | _ => Inst::Const(f64::NAN),
                }
            })
            .collect();
        Some(Tape {
            insts,
            args: Vec::new(),
            inputs: self.inputs,
            params: self.params,
            outputs: self.outputs.clone(),
        })
    }

    /// Checks SSA order, index ranges, register types and phases. Every
    /// device calls this before it runs or emits a program, so a malformed
    /// program cannot index outside a register file or a buffer.
    ///
    /// # Errors
    /// [`DeviceError::InvalidKernel`].
    pub fn validate(&self) -> Result<(), DeviceError> {
        self.analyse().map(|_| ())
    }

    /// Validates the program and returns the type and phase of every
    /// register.
    ///
    /// # Errors
    /// [`DeviceError::InvalidKernel`].
    #[allow(clippy::too_many_lines)]
    pub fn analyse(&self) -> Result<Analysis, DeviceError> {
        let bad = |why: String| Err(DeviceError::InvalidKernel(why));
        if self.outputs.is_empty() {
            return bad("a kernel needs at least one output".into());
        }
        let n = self.ops.len();
        let mut ty: Vec<Ty> = Vec::with_capacity(n);
        let mut phase: Vec<Phase> = Vec::with_capacity(n);
        let accs = self.fold.as_ref().map_or(0, |f| f.accs.len());
        for (at, op) in self.ops.iter().enumerate() {
            let (operands, count) = op.operands();
            let mut here = match *op {
                | KernelOp::Inner => Phase::Body,
                | KernelOp::Acc(_) => Phase::Post,
                | _ => Phase::Pre,
            };
            for (slot, &r) in operands.iter().take(count).enumerate() {
                let r = r as usize;
                if r >= at {
                    return bad(format!("op {at} reads a register defined later"));
                }
                // `select` takes a numeric condition and two values of one
                // type; every other instruction has uniformly typed operands.
                let want = match *op {
                    | KernelOp::Select(_, a, _) if slot > 0 => {
                        ty.get(a as usize).copied().unwrap_or(Ty::Num)
                    },
                    | _ => op.operand_ty(),
                };
                let got = ty.get(r).copied().unwrap_or(Ty::Num);
                if got != want {
                    return bad(format!(
                        "op {at}: operand {slot} is {got:?}, expected {want:?}"
                    ));
                }
                here = match (here, phase.get(r).copied().unwrap_or(Phase::Pre)) {
                    | (p, Phase::Pre) => p,
                    | (Phase::Pre, p) => p,
                    | (p, q) if p == q => p,
                    | _ => {
                        return bad(format!(
                            "op {at} depends on both the inner index and an accumulator"
                        ));
                    },
                };
            }
            match *op {
                | KernelOp::Input(i) if i as usize >= self.inputs => {
                    return bad(format!("op {at}: input {i} of {}", self.inputs));
                },
                | KernelOp::Len(i) | KernelOp::Load(i, _) if i as usize >= self.inputs => {
                    return bad(format!("op {at}: column {i} of {}", self.inputs));
                },
                | KernelOp::Param(i) if i as usize >= self.params => {
                    return bad(format!("op {at}: parameter {i} of {}", self.params));
                },
                | KernelOp::Gather(base, _) if base as usize > self.inputs => {
                    return bad(format!("op {at}: gather base {base} of {}", self.inputs));
                },
                | KernelOp::Inner if self.fold.is_none() => {
                    return bad(format!("op {at}: the inner index needs a fold"));
                },
                | KernelOp::Acc(k) if k as usize >= accs => {
                    return bad(format!("op {at}: accumulator {k} of {accs}"));
                },
                | KernelOp::IConst(c) if i32::try_from(c).is_err() => {
                    return bad(format!("op {at}: index constant {c} does not fit 32 bits"));
                },
                | _ => {},
            }
            ty.push(match *op {
                | KernelOp::Select(_, a, _) => ty.get(a as usize).copied().unwrap_or(Ty::Num),
                | _ => op.ty(),
            });
            phase.push(here);
        }
        for &r in &self.outputs {
            let r = r as usize;
            let (Some(&t), Some(&p)) = (ty.get(r), phase.get(r)) else {
                return bad(format!("output register {r} is out of range"));
            };
            if t != Ty::Num {
                return bad(format!(
                    "output register {r} is an index; convert it with ToNum"
                ));
            }
            if p == Phase::Body {
                return bad(format!(
                    "output register {r} is computed inside the fold; reduce it into an \
                     accumulator"
                ));
            }
        }
        if let Some(fold) = &self.fold {
            let count = fold.count as usize;
            if ty.get(count) != Some(&Ty::Idx) || phase.get(count) != Some(&Phase::Pre) {
                return bad("the fold count must be an index computed before the loop".into());
            }
            if fold.accs.is_empty() {
                return bad("a fold needs at least one accumulator".into());
            }
            for &(_, r) in &fold.accs {
                let r = r as usize;
                if ty.get(r) != Some(&Ty::Num) || matches!(phase.get(r), None | Some(Phase::Post)) {
                    return bad("a fold combines numbers computed before or inside the loop".into());
                }
            }
        }
        let of = |want: Phase| {
            phase
                .iter()
                .enumerate()
                .filter(|&(_, p)| *p == want)
                .map(|(k, _)| k)
                .collect::<Vec<usize>>()
        };
        Ok(Analysis {
            pre: of(Phase::Pre),
            body: of(Phase::Body),
            post: of(Phase::Post),
            ty,
            phase,
        })
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
                | KernelOp::Row => 15u8.hash(&mut h),
                | KernelOp::Rows => 16u8.hash(&mut h),
                | KernelOp::Len(i) => (17u8, i).hash(&mut h),
                | KernelOp::Inner => 18u8.hash(&mut h),
                | KernelOp::IConst(c) => (19u8, c).hash(&mut h),
                | KernelOp::IAdd(a, b) => (20u8, a, b).hash(&mut h),
                | KernelOp::ISub(a, b) => (21u8, a, b).hash(&mut h),
                | KernelOp::IMul(a, b) => (22u8, a, b).hash(&mut h),
                | KernelOp::IDiv(a, b) => (23u8, a, b).hash(&mut h),
                | KernelOp::IMod(a, b) => (24u8, a, b).hash(&mut h),
                | KernelOp::IMin(a, b) => (25u8, a, b).hash(&mut h),
                | KernelOp::IMax(a, b) => (26u8, a, b).hash(&mut h),
                | KernelOp::ICmp(c, a, b) => (27u8, c, a, b).hash(&mut h),
                | KernelOp::ToNum(a) => (28u8, a).hash(&mut h),
                | KernelOp::ToIdx(a) => (29u8, a).hash(&mut h),
                | KernelOp::Load(c, i) => (30u8, c, i).hash(&mut h),
                | KernelOp::Acc(k) => (31u8, k).hash(&mut h),
            }
        }
        self.fold.hash(&mut h);
        h.finish()
    }

    /// Reference evaluation of one row of an elementwise launch (the
    /// semantics every device must reproduce, up to its precision): `inputs`
    /// holds that row's element of every column. Slow; for tests and
    /// fallbacks. A malformed program evaluates to NaN.
    #[must_use]
    pub fn eval_row(
        &self,
        inputs: &[f64],
        params: &[f64],
    ) -> Vec<f64> {
        let lens = vec![1_usize; inputs.len()];
        let offsets = column_offsets(&lens);
        let cols = Columns {
            data: inputs,
            offsets: &offsets,
            lens: &lens,
        };
        let mut out = vec![f64::NAN; self.outputs.len()];
        if let Ok(analysis) = self.analyse() {
            self.eval_at(
                &analysis,
                0,
                1,
                cols,
                params,
                &mut Scratch::default(),
                &mut out,
            );
        }
        out
    }

    /// Reference evaluation of row `row` of a launch of `rows` rows over
    /// `cols`: the semantics of the whole IR, fold included. `analysis` must
    /// come from [`KernelIr::analyse`] of this program; `out` receives one
    /// value per output.
    #[allow(clippy::too_many_arguments)]
    pub fn eval_at(
        &self,
        analysis: &Analysis,
        row: usize,
        rows: usize,
        cols: Columns<'_>,
        params: &[f64],
        scratch: &mut Scratch,
        out: &mut [f64],
    ) {
        let n = self.ops.len();
        scratch.num.clear();
        scratch.num.resize(n, f64::NAN);
        scratch.idx.clear();
        scratch.idx.resize(n, 0);
        let row_i = i64::try_from(row).unwrap_or(i64::MAX);
        let rows_i = i64::try_from(rows).unwrap_or(i64::MAX);
        for &d in &analysis.pre {
            self.step(analysis, d, row_i, rows_i, 0, cols, params, scratch);
        }
        if let Some(fold) = &self.fold {
            scratch.acc.clear();
            scratch
                .acc
                .extend(fold.accs.iter().map(|(op, _)| op.identity()));
            let count = scratch.idx.get(fold.count as usize).copied().unwrap_or(0);
            for inner in 0..count.max(0) {
                for &d in &analysis.body {
                    self.step(analysis, d, row_i, rows_i, inner, cols, params, scratch);
                }
                for (k, &(op, r)) in fold.accs.iter().enumerate() {
                    let v = scratch.num.get(r as usize).copied().unwrap_or(f64::NAN);
                    if let Some(acc) = scratch.acc.get_mut(k) {
                        *acc = op.combine(*acc, v);
                    }
                }
            }
            for &d in &analysis.post {
                self.step(analysis, d, row_i, rows_i, 0, cols, params, scratch);
            }
        }
        for (slot, &r) in out.iter_mut().zip(&self.outputs) {
            *slot = scratch.num.get(r as usize).copied().unwrap_or(f64::NAN);
        }
    }

    /// Computes register `d`.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    #[inline]
    fn step(
        &self,
        analysis: &Analysis,
        d: usize,
        row: i64,
        rows: i64,
        inner: i64,
        cols: Columns<'_>,
        params: &[f64],
        scratch: &mut Scratch,
    ) {
        let Some(op) = self.ops.get(d) else {
            return;
        };
        let g = |r: Reg| scratch.num.get(r as usize).copied().unwrap_or(f64::NAN);
        let i = |r: Reg| scratch.idx.get(r as usize).copied().unwrap_or(0);
        let value = match *op {
            | KernelOp::Const(c) => Val::Num(c),
            | KernelOp::Input(col) => Val::Num(cols.get(col as usize, row)),
            | KernelOp::Param(k) => Val::Num(params.get(k as usize).copied().unwrap_or(f64::NAN)),
            | KernelOp::Add(a, b) => Val::Num(g(a) + g(b)),
            | KernelOp::Sub(a, b) => Val::Num(g(a) - g(b)),
            | KernelOp::Mul(a, b) => Val::Num(g(a) * g(b)),
            | KernelOp::Div(a, b) => Val::Num(g(a) / g(b)),
            | KernelOp::Neg(a) => Val::Num(-g(a)),
            | KernelOp::PowI(a, e) => Val::Num(int_pow(g(a), e)),
            | KernelOp::PowF(a, b) => Val::Num(g(a).powf(g(b))),
            | KernelOp::Unary(f, a) => Val::Num(f.apply(g(a))),
            | KernelOp::Atan2(a, b) => Val::Num(g(a).atan2(g(b))),
            | KernelOp::Cmp(c, a, b) => Val::Num(cmp_value(c, g(a), g(b))),
            | KernelOp::Select(c, a, b) => {
                let first = select_is_true(g(c));
                if analysis.ty.get(d) == Some(&Ty::Idx) {
                    Val::Idx(if first { i(a) } else { i(b) })
                } else {
                    Val::Num(if first { g(a) } else { g(b) })
                }
            },
            | KernelOp::Gather(base, off) => {
                Val::Num(
                    gather_slot(base, g(off), self.inputs)
                        .map_or(f64::NAN, |slot| cols.get(slot, row)),
                )
            },
            | KernelOp::Row => Val::Idx(row),
            | KernelOp::Rows => Val::Idx(rows),
            | KernelOp::Len(col) => {
                Val::Idx(
                    cols.lens
                        .get(col as usize)
                        .and_then(|&l| i64::try_from(l).ok())
                        .unwrap_or(0),
                )
            },
            | KernelOp::Inner => Val::Idx(inner),
            | KernelOp::IConst(c) => Val::Idx(c),
            | KernelOp::IAdd(a, b) => Val::Idx(i(a).wrapping_add(i(b))),
            | KernelOp::ISub(a, b) => Val::Idx(i(a).wrapping_sub(i(b))),
            | KernelOp::IMul(a, b) => Val::Idx(i(a).wrapping_mul(i(b))),
            | KernelOp::IDiv(a, b) => Val::Idx(i(a).checked_div(i(b)).unwrap_or(0)),
            | KernelOp::IMod(a, b) => Val::Idx(i(a).checked_rem(i(b)).unwrap_or(0)),
            | KernelOp::IMin(a, b) => Val::Idx(i(a).min(i(b))),
            | KernelOp::IMax(a, b) => Val::Idx(i(a).max(i(b))),
            | KernelOp::ICmp(c, a, b) => {
                let (x, y) = (i(a), i(b));
                let r = match c {
                    | Cmp::Lt => x < y,
                    | Cmp::Le => x <= y,
                    | Cmp::Gt => x > y,
                    | Cmp::Ge => x >= y,
                    | Cmp::Eq => x == y,
                    | Cmp::Ne => x != y,
                };
                Val::Num(f64::from(u8::from(r)))
            },
            #[allow(clippy::cast_precision_loss)]
            | KernelOp::ToNum(a) => Val::Num(i(a) as f64),
            | KernelOp::ToIdx(a) => Val::Idx(num_to_idx(g(a))),
            | KernelOp::Load(col, at) => Val::Num(cols.get(col as usize, i(at))),
            | KernelOp::Acc(k) => {
                Val::Num(scratch.acc.get(k as usize).copied().unwrap_or(f64::NAN))
            },
        };
        match value {
            | Val::Num(x) => {
                if let Some(slot) = scratch.num.get_mut(d) {
                    *slot = x;
                }
            },
            | Val::Idx(x) => {
                if let Some(slot) = scratch.idx.get_mut(d) {
                    *slot = x;
                }
            },
        }
    }
}

/// A number truncated to an index: NaN gives `0`, the range saturates at the
/// 32-bit limits every device shares.
#[inline]
#[must_use]
pub fn num_to_idx(x: f64) -> i64 {
    if x.is_nan() {
        return 0;
    }
    let t = x.trunc().clamp(f64::from(i32::MIN), f64::from(i32::MAX));
    #[allow(clippy::cast_possible_truncation)]
    {
        t as i64
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
            | KernelOp::Row => Some((3, 0)),
            | KernelOp::Rows => Some((4, 0)),
            | KernelOp::Len(i) => Some((5, u64::from(i))),
            | KernelOp::Inner => Some((6, 0)),
            #[allow(clippy::cast_sign_loss)]
            | KernelOp::IConst(c) => Some((7, c as u64)),
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

    /// A numeric constant.
    pub fn num(
        &mut self,
        c: f64,
    ) -> Reg {
        self.push(KernelOp::Const(c))
    }

    /// An index constant.
    pub fn idx(
        &mut self,
        c: i64,
    ) -> Reg {
        self.push(KernelOp::IConst(c))
    }

    /// Splices `other` into this program with its input `i` bound to
    /// `inputs[i]` and its parameter `j` to `params[j]`; returns its outputs.
    /// `other` must not use `gather` (it addresses the caller's columns) and
    /// its own fold, if any, is not carried over: splice elementwise programs.
    pub fn splice(
        &mut self,
        other: &KernelIr,
        inputs: &[Reg],
        params: &[Reg],
    ) -> Vec<Reg> {
        let mut map: Vec<Reg> = Vec::with_capacity(other.ops.len());
        for op in &other.ops {
            let mapped = match *op {
                | KernelOp::Input(i) => {
                    map.push(inputs.get(i as usize).copied().unwrap_or(0));
                    continue;
                },
                | KernelOp::Param(j) => {
                    map.push(params.get(j as usize).copied().unwrap_or(0));
                    continue;
                },
                | KernelOp::Gather(..) => KernelOp::Const(f64::NAN),
                | other_op => other_op.map_operands(|r| map.get(r as usize).copied().unwrap_or(0)),
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
        self.finish_fold(outputs, None)
    }

    /// Finishes a program with a loop.
    ///
    /// # Errors
    /// [`DeviceError::InvalidKernel`] when the result does not validate.
    pub fn finish_fold(
        self,
        outputs: Vec<Reg>,
        fold: Option<Fold>,
    ) -> Result<KernelIr, DeviceError> {
        let ir = KernelIr {
            ops: self.ops,
            inputs: self.inputs,
            params: self.params,
            outputs,
            fold,
        };
        ir.validate()?;
        Ok(ir)
    }

    /// Number of instructions so far.
    #[must_use]
    pub fn len(&self) -> usize {
        self.ops.len()
    }

    /// Whether no instruction was pushed yet.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.ops.is_empty()
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
/// when `rhs` is not elementwise or uses `gather` (which addresses columns
/// and cannot be composed).
pub fn rk4_step(
    rhs: &KernelIr,
    n_state: usize,
) -> Result<KernelIr, DeviceError> {
    let timed = rhs_shape(rhs, n_state)?;
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

/// Checks that `rhs` is an ODE right-hand side for `n_state` variables and
/// reports whether it takes the time as input 0.
///
/// # Errors
/// [`DeviceError::InvalidKernel`].
pub fn rhs_shape(
    rhs: &KernelIr,
    n_state: usize,
) -> Result<bool, DeviceError> {
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
    if !rhs.is_elementwise() || rhs.ops.iter().any(|o| matches!(o, KernelOp::Gather(..))) {
        return Err(DeviceError::InvalidKernel(
            "an ODE right-hand side must be elementwise and cannot use gather".into(),
        ));
    }
    Ok(timed)
}
