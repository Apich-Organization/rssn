//! # LIR: the loop IR shared by every compiling front end
//!
//! The [tape](super::tape) is a straight-line expression program and the
//! [device kernel IR](super::device::KernelIr) is a per-row program with one
//! fold. Neither can say "loop until", "call that function", or "branch
//! without evaluating both arms", which is what a language front end (Aevia's
//! `@jit`), a batch evaluator ("for every row") or a solver kernel needs. LIR
//! is the one IR that can, and the others lower into it:
//!
//! ```text
//!  tape ──► batch / scalar LIR ─┐
//!  KernelIr ──► row-loop LIR ───┼──► LIR ──► interpreter (reference, always available)
//!  Aevia Core IR ──► LIR ───────┘          └► Cranelift (feature `jit`, `lir_jit`)
//! ```
//!
//! ## The model: a typed register machine with basic blocks
//!
//! A [`Module`] holds [`Function`]s. A function has typed *registers*
//! ([`Reg`]), basic [`Block`]s of [`Inst`]ructions, and one [`Term`]inator per
//! block. Registers are **mutable and zero-initialised**: an instruction
//! assigns its destination, a jump assigns the target block's parameters (in
//! parallel), and a register read before any assignment reads zero. This is
//! exactly the "variables" model of Cranelift's SSA builder, so code
//! generation is a direct transcription and needs no SSA construction of its
//! own; the interpreter gives the same semantics by construction.
//!
//! Types ([`Ty`]) are `F64`, `I64`, `Bool` and `Slice`. A slice is a fat
//! pointer to `f64`s (address and length) passed in by the host; [`Inst::Load`]
//! and [`Inst::Store`] are bounds-checked (an out-of-range load reads NaN, an
//! out-of-range store does nothing), so no LIR program can touch memory it was
//! not given. Slices cannot be created, offset or stored in LIR — only passed
//! along — which keeps the host-facing API safe.
//!
//! ## Semantics worth stating
//!
//! * Float arithmetic is IEEE `f64` with no reassociation; integer powers use
//!   the same square-and-multiply sequence as the tape ([`super::tape::int_pow`]),
//!   so a tape lowered to LIR returns the same bits.
//! * Integer arithmetic is *checked*: overflow, and division or remainder by
//!   zero, **bail** (see below). `IDiv`/`IRem` truncate toward zero.
//! * `FToI` truncates toward zero; NaN or out-of-range bails.
//! * [`Term::Bail`] abandons the whole call: the host gets
//!   [`Outcome::Bailed`] with the code, and no result. Front ends use it for
//!   "this case is not mine" (an overflow an arbitrary-precision language must
//!   handle, a guard that failed) and re-run the call in their interpreter.
//!   Because LIR has no effects other than stores into host slices, re-running
//!   is safe whenever the program wrote nothing it read.
//! * [`Inst::Call`] calls another function of the module directly; recursion
//!   and mutual recursion are allowed. A bail inside a callee bails the caller.
//!   Calls are limited to a fixed depth ([`MAX_CALL_DEPTH`]); exceeding it
//!   bails with [`BAIL_DEPTH`] in both executors.
//!
//! ## Why a register machine and not structured regions
//!
//! Structured control flow (MLIR's `scf`) is easier to emit for GPUs, but a
//! front end's control flow — Aevia's tail calls and `match` fall-through —
//! is not always structured, and Cranelift wants a CFG anyway. Device
//! programs keep their own structured IR; LIR is the CPU-side common
//! denominator.

use std::fmt;

use super::Intrinsic;
use super::tape::Cmp;
use crate::graph::op::EvalFn;

/// A register of a function.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Reg(pub u32);

/// A basic block of a function.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct BlockId(pub u32);

/// A function of a module.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct FuncId(pub u32);

/// The type of a register.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Ty {
    /// A double.
    F64,
    /// A signed 64-bit integer.
    I64,
    /// A boolean.
    Bool,
    /// A host-provided slice of `f64` (address and length).
    Slice,
}

/// One instruction; it assigns its destination register.
#[derive(Clone, Debug)]
pub enum Inst {
    /// A float constant.
    F(f64),
    /// An integer constant.
    I(i64),
    /// A boolean constant.
    B(bool),
    /// A copy of a register (any type).
    Copy(Reg),
    /// `a + b`.
    FAdd(Reg, Reg),
    /// `a - b`.
    FSub(Reg, Reg),
    /// `a * b`.
    FMul(Reg, Reg),
    /// `a / b`.
    FDiv(Reg, Reg),
    /// `-a`.
    FNeg(Reg),
    /// `a^n`, `n` an integer constant (square-and-multiply).
    FPowI(Reg, i32),
    /// `a^b`.
    FPow(Reg, Reg),
    /// An elementary function.
    FUnary(Intrinsic, Reg),
    /// `atan2(a, b)`.
    FAtan2(Reg, Reg),
    /// The smaller of two floats (NaN if either is NaN).
    FMin(Reg, Reg),
    /// The larger of two floats (NaN if either is NaN).
    FMax(Reg, Reg),
    /// Checked `a + b`.
    IAdd(Reg, Reg),
    /// Checked `a - b`.
    ISub(Reg, Reg),
    /// Checked `a * b`.
    IMul(Reg, Reg),
    /// Checked `a / b`, truncated.
    IDiv(Reg, Reg),
    /// Checked `a % b` (sign of `a`).
    IRem(Reg, Reg),
    /// Checked `-a`.
    INeg(Reg),
    /// The smaller integer.
    IMin(Reg, Reg),
    /// The larger integer.
    IMax(Reg, Reg),
    /// A float comparison (false for NaN operands except `Ne`).
    FCmp(Cmp, Reg, Reg),
    /// An integer comparison.
    ICmp(Cmp, Reg, Reg),
    /// Boolean negation.
    Not(Reg),
    /// Boolean and (both operands already evaluated).
    And(Reg, Reg),
    /// Boolean or (both operands already evaluated).
    Or(Reg, Reg),
    /// An integer as a float (rounded to nearest).
    IToF(Reg),
    /// A float truncated to an integer; NaN or out of range bails.
    FToI(Reg),
    /// `if c { a } else { b }`, both of one type.
    Select(Reg, Reg, Reg),
    /// The length of a slice.
    Len(Reg),
    /// `slice[index]`, NaN when out of range.
    Load(Reg, Reg),
    /// `slice[index] = value` when in range; assigns its destination `0.0`.
    Store(Reg, Reg, Reg),
    /// A direct call of a function of the module returning one value.
    Call(FuncId, Vec<Reg>),
    /// A call of an operator's scalar semantics (host code).
    Eval(EvalFn, Vec<Reg>),
    /// An external C function of one argument.
    Extern1(extern "C" fn(f64) -> f64, Reg),
    /// An external C function of two arguments.
    Extern2(extern "C" fn(f64, f64) -> f64, Reg, Reg),
}

/// How a block ends.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Term {
    /// Assign the target's parameters (in parallel) and continue there.
    Jump(BlockId, Vec<Reg>),
    /// `if cond` jump to the first target, else to the second.
    Branch(Reg, (BlockId, Vec<Reg>), (BlockId, Vec<Reg>)),
    /// Return these registers (as many as the function's result types).
    Return(Vec<Reg>),
    /// Abandon the call with a code (see the module documentation).
    Bail(u32),
}

/// A basic block.
#[derive(Clone, Debug)]
pub struct Block {
    /// Parameters, assigned by the jumps into the block.
    pub params: Vec<Reg>,
    /// Instructions and their destinations.
    pub insts: Vec<(Reg, Inst)>,
    /// The terminator.
    pub term: Term,
}

/// A function: the entry block's parameters are its arguments.
#[derive(Clone, Debug)]
pub struct Function {
    /// A name, for diagnostics and symbols.
    pub name: String,
    /// Register types; register `r` has type `regs[r]`.
    pub regs: Vec<Ty>,
    /// Blocks; block 0 is the entry.
    pub blocks: Vec<Block>,
    /// Result types.
    pub results: Vec<Ty>,
}

/// Several functions that may call each other.
#[derive(Clone, Debug, Default)]
pub struct Module {
    /// The functions, indexed by [`FuncId`].
    pub functions: Vec<Function>,
}

/// Bail code: integer overflow, or division or remainder by zero.
pub const BAIL_OVERFLOW: u32 = 1;
/// Bail code: a float that does not convert to an integer.
pub const BAIL_CONVERT: u32 = 2;
/// Bail code: the call depth limit was exceeded.
pub const BAIL_DEPTH: u32 = 3;
/// Bail code: the interpreter's step budget ran out (interpreter only).
pub const BAIL_STEPS: u32 = 4;
/// Deepest call nesting either executor allows before bailing.
pub const MAX_CALL_DEPTH: u32 = 512;

/// Why a module is malformed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LirError(pub String);

impl fmt::Display for LirError {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        write!(f, "malformed LIR: {}", self.0)
    }
}

impl std::error::Error for LirError {}

impl From<LirError> for super::BackendError {
    fn from(e: LirError) -> Self {
        Self::Codegen(e.to_string())
    }
}

fn bad<T>(message: impl Into<String>) -> Result<T, LirError> {
    Err(LirError(message.into()))
}

// ---------------------------------------------------------------------------
// Building
// ---------------------------------------------------------------------------

/// Builds one function block by block.
///
/// ```
/// use rssn::backend::lir::FunctionBuilder;
/// use rssn::backend::lir::Inst;
/// use rssn::backend::lir::Module;
/// use rssn::backend::lir::Ty;
/// use rssn::backend::lir::Value;
/// use rssn::backend::tape::Cmp;
/// // sum(n) = 0 + 1 + ... + (n - 1), as a loop.
/// let mut f = FunctionBuilder::new("sum", &[Ty::I64], &[Ty::I64]);
/// let n = f.param(0);
/// let (head, body, exit) = (f.new_block(), f.new_block(), f.new_block());
/// let (i, acc) = (f.block_param(head, Ty::I64), f.block_param(head, Ty::I64));
/// let zero = f.ins(Inst::I(0));
/// f.jump(head, &[zero, zero]);
/// f.switch_to(head);
/// let more = f.ins(Inst::ICmp(Cmp::Lt, i, n));
/// f.branch(more, body, &[], exit, &[]);
/// f.switch_to(body);
/// let one = f.ins(Inst::I(1));
/// let (acc2, i2) = (f.ins(Inst::IAdd(acc, i)), f.ins(Inst::IAdd(i, one)));
/// f.jump(head, &[i2, acc2]);
/// f.switch_to(exit);
/// f.ret(&[acc]);
/// let module = Module {
///     functions: vec![f.finish().unwrap()],
/// };
/// let out = module.interpret(rssn::backend::lir::FuncId(0), &mut [Value::I(10)], 1 << 20);
/// assert_eq!(out.unwrap().values(), Some(&[Value::I(45)][..]));
/// ```
#[derive(Debug)]
pub struct FunctionBuilder {
    function: Function,
    current: usize,
    terminated: Vec<bool>,
}

impl FunctionBuilder {
    /// A function with the given parameter and result types; the builder is
    /// positioned in the entry block.
    #[must_use]
    pub fn new(
        name: &str,
        params: &[Ty],
        results: &[Ty],
    ) -> Self {
        let mut function = Function {
            name: name.to_owned(),
            regs: Vec::new(),
            blocks: vec![Block {
                params: Vec::new(),
                insts: Vec::new(),
                term: Term::Bail(0),
            }],
            results: results.to_vec(),
        };
        for &ty in params {
            let reg = Reg(u32::try_from(function.regs.len()).unwrap_or(u32::MAX));
            function.regs.push(ty);
            if let Some(entry) = function.blocks.first_mut() {
                entry.params.push(reg);
            }
        }
        Self {
            function,
            current: 0,
            terminated: vec![false],
        }
    }

    /// Parameter `i` of the function.
    #[must_use]
    pub fn param(
        &self,
        i: usize,
    ) -> Reg {
        self.function
            .blocks
            .first()
            .and_then(|b| b.params.get(i).copied())
            .unwrap_or(Reg(u32::MAX))
    }

    /// A fresh register of type `ty`.
    pub fn reg(
        &mut self,
        ty: Ty,
    ) -> Reg {
        let reg = Reg(u32::try_from(self.function.regs.len()).unwrap_or(u32::MAX));
        self.function.regs.push(ty);
        reg
    }

    /// The type of a register.
    #[must_use]
    pub fn ty(
        &self,
        reg: Reg,
    ) -> Option<Ty> {
        self.function.regs.get(reg.0 as usize).copied()
    }

    /// A new, empty block.
    pub fn new_block(&mut self) -> BlockId {
        let id = BlockId(u32::try_from(self.function.blocks.len()).unwrap_or(u32::MAX));
        self.function.blocks.push(Block {
            params: Vec::new(),
            insts: Vec::new(),
            term: Term::Bail(0),
        });
        self.terminated.push(false);
        id
    }

    /// Adds a parameter of type `ty` to `block`.
    pub fn block_param(
        &mut self,
        block: BlockId,
        ty: Ty,
    ) -> Reg {
        let reg = self.reg(ty);
        if let Some(b) = self.function.blocks.get_mut(block.0 as usize) {
            b.params.push(reg);
        }
        reg
    }

    /// Continues emitting into `block`.
    pub const fn switch_to(
        &mut self,
        block: BlockId,
    ) {
        self.current = block.0 as usize;
    }

    /// The block being emitted into.
    #[must_use]
    pub const fn current(&self) -> BlockId {
        BlockId(self.current as u32)
    }

    /// Whether the current block already has its terminator.
    #[must_use]
    pub fn is_terminated(&self) -> bool {
        self.terminated.get(self.current).copied().unwrap_or(true)
    }

    /// Emits `inst` into a fresh register of its natural type and returns it.
    /// For [`Inst::Call`] and [`Inst::Copy`], whose type is not determined by
    /// the instruction alone, use [`FunctionBuilder::ins_typed`].
    pub fn ins(
        &mut self,
        inst: Inst,
    ) -> Reg {
        let ty = self.natural_type(&inst).unwrap_or(Ty::F64);
        self.ins_typed(inst, ty)
    }

    /// Emits `inst` into a fresh register of type `ty`.
    pub fn ins_typed(
        &mut self,
        inst: Inst,
        ty: Ty,
    ) -> Reg {
        let reg = self.reg(ty);
        self.assign(reg, inst);
        reg
    }

    /// Emits `inst` assigning an existing register.
    pub fn assign(
        &mut self,
        reg: Reg,
        inst: Inst,
    ) {
        if let Some(b) = self.function.blocks.get_mut(self.current) {
            b.insts.push((reg, inst));
        }
    }

    fn natural_type(
        &self,
        inst: &Inst,
    ) -> Option<Ty> {
        Some(match inst {
            | Inst::F(_)
            | Inst::FAdd(..)
            | Inst::FSub(..)
            | Inst::FMul(..)
            | Inst::FDiv(..)
            | Inst::FNeg(_)
            | Inst::FPowI(..)
            | Inst::FPow(..)
            | Inst::FUnary(..)
            | Inst::FAtan2(..)
            | Inst::FMin(..)
            | Inst::FMax(..)
            | Inst::IToF(_)
            | Inst::Load(..)
            | Inst::Store(..)
            | Inst::Eval(..)
            | Inst::Extern1(..)
            | Inst::Extern2(..) => Ty::F64,
            | Inst::I(_)
            | Inst::IAdd(..)
            | Inst::ISub(..)
            | Inst::IMul(..)
            | Inst::IDiv(..)
            | Inst::IRem(..)
            | Inst::INeg(_)
            | Inst::IMin(..)
            | Inst::IMax(..)
            | Inst::FToI(_)
            | Inst::Len(_) => Ty::I64,
            | Inst::B(_)
            | Inst::FCmp(..)
            | Inst::ICmp(..)
            | Inst::Not(_)
            | Inst::And(..)
            | Inst::Or(..) => Ty::Bool,
            | Inst::Copy(r) | Inst::Select(_, r, _) => self.ty(*r)?,
            | Inst::Call(..) => return None,
        })
    }

    fn terminate(
        &mut self,
        term: Term,
    ) {
        if let Some(b) = self.function.blocks.get_mut(self.current) {
            b.term = term;
        }
        if let Some(t) = self.terminated.get_mut(self.current) {
            *t = true;
        }
    }

    /// Ends the current block with a jump.
    pub fn jump(
        &mut self,
        target: BlockId,
        args: &[Reg],
    ) {
        self.terminate(Term::Jump(target, args.to_vec()));
    }

    /// Ends the current block with a conditional branch.
    pub fn branch(
        &mut self,
        cond: Reg,
        then_block: BlockId,
        then_args: &[Reg],
        else_block: BlockId,
        else_args: &[Reg],
    ) {
        self.terminate(Term::Branch(
            cond,
            (then_block, then_args.to_vec()),
            (else_block, else_args.to_vec()),
        ));
    }

    /// Ends the current block with a return.
    pub fn ret(
        &mut self,
        values: &[Reg],
    ) {
        self.terminate(Term::Return(values.to_vec()));
    }

    /// Ends the current block with a bail.
    pub fn bail(
        &mut self,
        code: u32,
    ) {
        self.terminate(Term::Bail(code));
    }

    /// The finished function.
    ///
    /// # Errors
    /// A block without a terminator.
    pub fn finish(self) -> Result<Function, LirError> {
        if let Some(open) = self.terminated.iter().position(|t| !t) {
            return bad(format!(
                "{}: block {open} has no terminator",
                self.function.name
            ));
        }
        Ok(self.function)
    }
}

// ---------------------------------------------------------------------------
// Verification
// ---------------------------------------------------------------------------

impl Module {
    /// Checks types, arities and references.
    ///
    /// # Errors
    /// The first problem found.
    pub fn verify(&self) -> Result<(), LirError> {
        for (index, f) in self.functions.iter().enumerate() {
            self.verify_function(f)
                .map_err(|e| LirError(format!("function {index} `{}`: {}", f.name, e.0)))?;
        }
        Ok(())
    }

    fn verify_function(
        &self,
        f: &Function,
    ) -> Result<(), LirError> {
        let ty = |r: Reg| -> Result<Ty, LirError> {
            f.regs
                .get(r.0 as usize)
                .copied()
                .ok_or_else(|| LirError(format!("unknown register r{}", r.0)))
        };
        let want = |r: Reg, t: Ty| -> Result<(), LirError> {
            let have = ty(r)?;
            if have == t {
                Ok(())
            } else {
                bad(format!("r{} is {have:?}, expected {t:?}", r.0))
            }
        };
        let jump = |target: BlockId, args: &[Reg]| -> Result<(), LirError> {
            let block = f
                .blocks
                .get(target.0 as usize)
                .ok_or_else(|| LirError(format!("unknown block b{}", target.0)))?;
            if block.params.len() != args.len() {
                return bad(format!(
                    "jump to b{} with {} arguments for {} parameters",
                    target.0,
                    args.len(),
                    block.params.len()
                ));
            }
            for (&p, &a) in block.params.iter().zip(args) {
                want(a, ty(p)?)?;
            }
            Ok(())
        };
        if f.blocks.is_empty() {
            return bad("no blocks");
        }
        for (b, block) in f.blocks.iter().enumerate() {
            for &p in &block.params {
                ty(p)?;
            }
            for (dst, inst) in &block.insts {
                let d = ty(*dst)?;
                let expect = |t: Ty| -> Result<(), LirError> {
                    if d == t {
                        Ok(())
                    } else {
                        bad(format!(
                            "b{b}: destination r{} is {d:?}, expected {t:?}",
                            dst.0
                        ))
                    }
                };
                match inst {
                    | Inst::F(_) => expect(Ty::F64)?,
                    | Inst::I(_) => expect(Ty::I64)?,
                    | Inst::B(_) => expect(Ty::Bool)?,
                    | Inst::Copy(a) => expect(ty(*a)?)?,
                    | Inst::FAdd(a, c)
                    | Inst::FSub(a, c)
                    | Inst::FMul(a, c)
                    | Inst::FDiv(a, c)
                    | Inst::FPow(a, c)
                    | Inst::FAtan2(a, c)
                    | Inst::FMin(a, c)
                    | Inst::FMax(a, c)
                    | Inst::Extern2(_, a, c) => {
                        want(*a, Ty::F64)?;
                        want(*c, Ty::F64)?;
                        expect(Ty::F64)?;
                    },
                    | Inst::FNeg(a)
                    | Inst::FPowI(a, _)
                    | Inst::FUnary(_, a)
                    | Inst::Extern1(_, a) => {
                        want(*a, Ty::F64)?;
                        expect(Ty::F64)?;
                    },
                    | Inst::IAdd(a, c)
                    | Inst::ISub(a, c)
                    | Inst::IMul(a, c)
                    | Inst::IDiv(a, c)
                    | Inst::IRem(a, c)
                    | Inst::IMin(a, c)
                    | Inst::IMax(a, c) => {
                        want(*a, Ty::I64)?;
                        want(*c, Ty::I64)?;
                        expect(Ty::I64)?;
                    },
                    | Inst::INeg(a) => {
                        want(*a, Ty::I64)?;
                        expect(Ty::I64)?;
                    },
                    | Inst::FCmp(_, a, c) => {
                        want(*a, Ty::F64)?;
                        want(*c, Ty::F64)?;
                        expect(Ty::Bool)?;
                    },
                    | Inst::ICmp(_, a, c) => {
                        want(*a, Ty::I64)?;
                        want(*c, Ty::I64)?;
                        expect(Ty::Bool)?;
                    },
                    | Inst::Not(a) => {
                        want(*a, Ty::Bool)?;
                        expect(Ty::Bool)?;
                    },
                    | Inst::And(a, c) | Inst::Or(a, c) => {
                        want(*a, Ty::Bool)?;
                        want(*c, Ty::Bool)?;
                        expect(Ty::Bool)?;
                    },
                    | Inst::IToF(a) => {
                        want(*a, Ty::I64)?;
                        expect(Ty::F64)?;
                    },
                    | Inst::FToI(a) => {
                        want(*a, Ty::F64)?;
                        expect(Ty::I64)?;
                    },
                    | Inst::Select(c, x, y) => {
                        want(*c, Ty::Bool)?;
                        let t = ty(*x)?;
                        if t == Ty::Slice {
                            return bad("select of slices");
                        }
                        want(*y, t)?;
                        expect(t)?;
                    },
                    | Inst::Len(s) => {
                        want(*s, Ty::Slice)?;
                        expect(Ty::I64)?;
                    },
                    | Inst::Load(s, i) => {
                        want(*s, Ty::Slice)?;
                        want(*i, Ty::I64)?;
                        expect(Ty::F64)?;
                    },
                    | Inst::Store(s, i, v) => {
                        want(*s, Ty::Slice)?;
                        want(*i, Ty::I64)?;
                        want(*v, Ty::F64)?;
                        expect(Ty::F64)?;
                    },
                    | Inst::Call(callee, args) => {
                        let g = self
                            .functions
                            .get(callee.0 as usize)
                            .ok_or_else(|| LirError(format!("unknown function f{}", callee.0)))?;
                        let params = g.blocks.first().map_or(&[][..], |e| &e.params[..]);
                        if params.len() != args.len() {
                            return bad(format!("call of `{}` with wrong arity", g.name));
                        }
                        for (&p, &a) in params.iter().zip(args) {
                            want(a, g.regs.get(p.0 as usize).copied().unwrap_or(Ty::F64))?;
                        }
                        match g.results.as_slice() {
                            | [t] => expect(*t)?,
                            | _ => return bad(format!("call of `{}`: not one result", g.name)),
                        }
                    },
                    | Inst::Eval(_, args) => {
                        for &a in args {
                            want(a, Ty::F64)?;
                        }
                        expect(Ty::F64)?;
                    },
                }
            }
            match &block.term {
                | Term::Jump(t, args) => jump(*t, args)?,
                | Term::Branch(c, (t, ta), (e, ea)) => {
                    want(*c, Ty::Bool)?;
                    jump(*t, ta)?;
                    jump(*e, ea)?;
                },
                | Term::Return(values) => {
                    if values.len() != f.results.len() {
                        return bad(format!("b{b}: returns {} values", values.len()));
                    }
                    for (&v, &t) in values.iter().zip(&f.results) {
                        want(v, t)?;
                    }
                },
                | Term::Bail(_) => {},
            }
        }
        Ok(())
    }

    /// Which parameters of `func` may be stored into, by `func` itself or by
    /// anything it calls: a slice argument passed for such a parameter must
    /// be writable, every other one may be read-only.
    ///
    /// Slices cannot be created in LIR, only passed along, so the analysis
    /// follows each slice register back to the parameters it may alias
    /// (through copies, block parameters and call arguments) and iterates
    /// over the whole module to a fixpoint (calls may recurse).
    #[must_use]
    pub fn written_params(
        &self,
        func: FuncId,
    ) -> Vec<bool> {
        let all = self.written_all();
        all.get(func.0 as usize).cloned().unwrap_or_default()
    }

    fn written_all(&self) -> Vec<Vec<bool>> {
        // Per function: which entry parameters each register may alias.
        let provenance: Vec<Vec<Vec<usize>>> = self
            .functions
            .iter()
            .map(|f| {
                let mut from: Vec<Vec<usize>> = vec![Vec::new(); f.regs.len()];
                for (k, p) in f
                    .blocks
                    .first()
                    .map_or(&[][..], |e| &e.params[..])
                    .iter()
                    .enumerate()
                {
                    if let Some(slot) = from.get_mut(p.0 as usize) {
                        slot.push(k);
                    }
                }
                let mut changed = true;
                while changed {
                    changed = false;
                    let mut flow = |dst: Reg, src: Reg, from: &mut Vec<Vec<usize>>| {
                        let add = from.get(src.0 as usize).cloned().unwrap_or_default();
                        if let Some(slot) = from.get_mut(dst.0 as usize) {
                            for k in add {
                                if !slot.contains(&k) {
                                    slot.push(k);
                                    changed = true;
                                }
                            }
                        }
                    };
                    for block in &f.blocks {
                        for (dst, inst) in &block.insts {
                            if let Inst::Copy(src) = inst {
                                flow(*dst, *src, &mut from);
                            }
                        }
                        let jumps: Vec<(&BlockId, &Vec<Reg>)> = match &block.term {
                            | Term::Jump(t, args) => vec![(t, args)],
                            | Term::Branch(_, (t, ta), (e, ea)) => vec![(t, ta), (e, ea)],
                            | Term::Return(_) | Term::Bail(_) => Vec::new(),
                        };
                        for (target, args) in jumps {
                            if let Some(tb) = f.blocks.get(target.0 as usize) {
                                for (&p, &a) in tb.params.iter().zip(args) {
                                    flow(p, a, &mut from);
                                }
                            }
                        }
                    }
                }
                from
            })
            .collect();
        let mut written: Vec<Vec<bool>> = self
            .functions
            .iter()
            .map(|f| vec![false; f.blocks.first().map_or(0, |e| e.params.len())])
            .collect();
        let mut changed = true;
        while changed {
            changed = false;
            for (i, f) in self.functions.iter().enumerate() {
                let Some(from) = provenance.get(i) else {
                    continue;
                };
                let mut mark: Vec<usize> = Vec::new();
                for (_, inst) in f.blocks.iter().flat_map(|b| &b.insts) {
                    match inst {
                        | Inst::Store(slice, ..) => {
                            mark.extend(from.get(slice.0 as usize).into_iter().flatten());
                        },
                        | Inst::Call(g, args) => {
                            let callee = written.get(g.0 as usize).cloned().unwrap_or_default();
                            for (j, a) in args.iter().enumerate() {
                                if callee.get(j).copied().unwrap_or(false) {
                                    mark.extend(from.get(a.0 as usize).into_iter().flatten());
                                }
                            }
                        },
                        | _ => {},
                    }
                }
                if let Some(w) = written.get_mut(i) {
                    for k in mark {
                        if let Some(flag) = w.get_mut(k)
                            && !*flag
                        {
                            *flag = true;
                            changed = true;
                        }
                    }
                }
            }
        }
        written
    }
}

// ---------------------------------------------------------------------------
// Values and the reference interpreter
// ---------------------------------------------------------------------------

/// A host value passed to or returned from a LIR function.
#[derive(Debug, PartialEq)]
pub enum Value<'a> {
    /// A double.
    F(f64),
    /// An integer.
    I(i64),
    /// A boolean.
    B(bool),
    /// A slice the function may read.
    Read(&'a [f64]),
    /// A slice the function may read and write.
    Write(&'a mut [f64]),
}

impl Value<'_> {
    /// The LIR type of the value.
    #[must_use]
    pub const fn ty(&self) -> Ty {
        match self {
            | Self::F(_) => Ty::F64,
            | Self::I(_) => Ty::I64,
            | Self::B(_) => Ty::Bool,
            | Self::Read(_) | Self::Write(_) => Ty::Slice,
        }
    }
}

/// The result of a call.
#[derive(Debug, PartialEq)]
pub enum Outcome {
    /// The returned values.
    Returned(Vec<Value<'static>>),
    /// The call bailed with this code.
    Bailed(u32),
}

impl Outcome {
    /// The returned values, if the call did not bail.
    #[must_use]
    pub fn values(&self) -> Option<&[Value<'static>]> {
        match self {
            | Self::Returned(v) => Some(v),
            | Self::Bailed(_) => None,
        }
    }
}

/// An interpreter register: scalars by value, slices by index into the
/// call's slice table.
#[derive(Clone, Copy, Debug)]
enum Slot {
    F(f64),
    I(i64),
    B(bool),
    S(usize),
}

impl Slot {
    const fn zero(ty: Ty) -> Self {
        match ty {
            | Ty::F64 => Self::F(0.0),
            | Ty::I64 => Self::I(0),
            | Ty::Bool => Self::B(false),
            | Ty::Slice => Self::S(usize::MAX),
        }
    }

    const fn f(self) -> f64 {
        match self {
            | Self::F(x) => x,
            | _ => f64::NAN,
        }
    }

    const fn i(self) -> i64 {
        match self {
            | Self::I(x) => x,
            | _ => 0,
        }
    }

    const fn b(self) -> bool {
        match self {
            | Self::B(x) => x,
            | _ => false,
        }
    }
}

enum Slices<'s, 'a> {
    Host(&'s mut [Value<'a>]),
}

impl Slices<'_, '_> {
    fn len(
        &self,
        s: usize,
    ) -> i64 {
        let Self::Host(values) = self;
        match values.get(s) {
            | Some(Value::Read(x)) => i64::try_from(x.len()).unwrap_or(i64::MAX),
            | Some(Value::Write(x)) => i64::try_from(x.len()).unwrap_or(i64::MAX),
            | _ => 0,
        }
    }

    fn load(
        &self,
        s: usize,
        i: i64,
    ) -> f64 {
        let Self::Host(values) = self;
        let Ok(i) = usize::try_from(i) else {
            return f64::NAN;
        };
        match values.get(s) {
            | Some(Value::Read(x)) => x.get(i).copied().unwrap_or(f64::NAN),
            | Some(Value::Write(x)) => x.get(i).copied().unwrap_or(f64::NAN),
            | _ => f64::NAN,
        }
    }

    fn store(
        &mut self,
        s: usize,
        i: i64,
        v: f64,
    ) {
        let Self::Host(values) = self;
        let Ok(i) = usize::try_from(i) else {
            return;
        };
        if let Some(Value::Write(x)) = values.get_mut(s)
            && let Some(slot) = x.get_mut(i)
        {
            *slot = v;
        }
    }
}

/// The comparison of two floats.
#[must_use]
pub fn fcmp(
    c: Cmp,
    a: f64,
    b: f64,
) -> bool {
    match c {
        | Cmp::Lt => a < b,
        | Cmp::Le => a <= b,
        | Cmp::Gt => a > b,
        | Cmp::Ge => a >= b,
        | Cmp::Eq => a == b,
        | Cmp::Ne => a != b,
    }
}

/// The comparison of two integers.
#[must_use]
pub const fn icmp(
    c: Cmp,
    a: i64,
    b: i64,
) -> bool {
    match c {
        | Cmp::Lt => a < b,
        | Cmp::Le => a <= b,
        | Cmp::Gt => a > b,
        | Cmp::Ge => a >= b,
        | Cmp::Eq => a == b,
        | Cmp::Ne => a != b,
    }
}

/// `x` truncated toward zero to an integer, or `None` for NaN and values
/// outside the `i64` range (the bail condition of [`Inst::FToI`]).
#[must_use]
pub fn f64_to_i64(x: f64) -> Option<i64> {
    // 2^63 is exactly representable; every float in (-2^63 - 1, 2^63) truncates into range.
    const LIMIT: f64 = 9_223_372_036_854_775_808.0;
    if x.is_nan() || x >= LIMIT || x < -LIMIT {
        return None;
    }
    #[allow(clippy::cast_possible_truncation)]
    Some(x.trunc() as i64)
}

/// The float minimum with NaN propagation (both executors agree on it).
#[must_use]
pub fn fmin(
    a: f64,
    b: f64,
) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        a.min(b)
    }
}

/// The float maximum with NaN propagation.
#[must_use]
pub fn fmax(
    a: f64,
    b: f64,
) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        a.max(b)
    }
}

struct Interp<'m> {
    module: &'m Module,
    steps: u64,
}

enum Flow {
    Return(Vec<Slot>),
    Bail(u32),
}

impl Module {
    /// Runs `func` in the reference interpreter with at most `max_steps`
    /// instructions (a runaway loop bails with [`BAIL_STEPS`]).
    ///
    /// # Errors
    /// The module does not verify, or the arguments do not match the
    /// function's parameters, or a read-only slice is passed to a function
    /// that writes.
    pub fn interpret(
        &self,
        func: FuncId,
        args: &mut [Value<'_>],
        max_steps: u64,
    ) -> Result<Outcome, LirError> {
        self.verify()?;
        self.check_args(func, args)?;
        let mut slots: Vec<Slot> = Vec::with_capacity(args.len());
        for (k, a) in args.iter().enumerate() {
            slots.push(match a {
                | Value::F(x) => Slot::F(*x),
                | Value::I(x) => Slot::I(*x),
                | Value::B(x) => Slot::B(*x),
                // A slice register holds the index of its argument.
                | Value::Read(_) | Value::Write(_) => Slot::S(k),
            });
        }
        let mut interp = Interp {
            module: self,
            steps: max_steps,
        };
        let mut slices = Slices::Host(args);
        Ok(match interp.run(func, &slots, &mut slices, 0) {
            | Flow::Return(values) => {
                let results = &self
                    .functions
                    .get(func.0 as usize)
                    .map_or(&[][..], |f| &f.results[..]);
                Outcome::Returned(
                    values
                        .iter()
                        .zip(results.iter())
                        .map(|(v, t)| {
                            match t {
                                | Ty::F64 => Value::F(v.f()),
                                | Ty::I64 => Value::I(v.i()),
                                | Ty::Bool | Ty::Slice => Value::B(v.b()),
                            }
                        })
                        .collect(),
                )
            },
            | Flow::Bail(code) => Outcome::Bailed(code),
        })
    }

    /// Checks host arguments against a function's parameters.
    ///
    /// # Errors
    /// Unknown function, wrong arity or types, or a read-only slice for a
    /// function that writes, or a function returning a slice.
    pub fn check_args(
        &self,
        func: FuncId,
        args: &[Value<'_>],
    ) -> Result<&Function, LirError> {
        let f = self
            .functions
            .get(func.0 as usize)
            .ok_or_else(|| LirError(format!("unknown function f{}", func.0)))?;
        let params = f.blocks.first().map_or(&[][..], |e| &e.params[..]);
        if params.len() != args.len() {
            return bad(format!("`{}` takes {} arguments", f.name, params.len()));
        }
        let written = self.written_params(func);
        for (k, (&p, a)) in params.iter().zip(args).enumerate() {
            let t = f.regs.get(p.0 as usize).copied().unwrap_or(Ty::F64);
            if t != a.ty() {
                return bad(format!("argument {k} of `{}` must be {t:?}", f.name));
            }
            if written.get(k).copied().unwrap_or(false) && matches!(a, Value::Read(_)) {
                return bad(format!(
                    "argument {k} of `{}`: the function may write it, so it must be writable",
                    f.name
                ));
            }
        }
        if f.results.contains(&Ty::Slice) {
            return bad(format!("`{}` returns a slice", f.name));
        }
        Ok(f)
    }
}

impl Interp<'_> {
    #[allow(clippy::too_many_lines)]
    fn run(
        &mut self,
        func: FuncId,
        args: &[Slot],
        slices: &mut Slices<'_, '_>,
        depth: u32,
    ) -> Flow {
        if depth >= MAX_CALL_DEPTH {
            return Flow::Bail(BAIL_DEPTH);
        }
        let Some(f) = self.module.functions.get(func.0 as usize) else {
            return Flow::Bail(0);
        };
        let mut r: Vec<Slot> = f.regs.iter().map(|&t| Slot::zero(t)).collect();
        let mut block = 0_usize;
        let mut incoming: Vec<Slot> = args.to_vec();
        loop {
            let Some(b) = f.blocks.get(block) else {
                return Flow::Bail(0);
            };
            for (&p, &v) in b.params.iter().zip(&incoming) {
                if let Some(slot) = r.get_mut(p.0 as usize) {
                    *slot = v;
                }
            }
            for (dst, inst) in &b.insts {
                if self.steps == 0 {
                    return Flow::Bail(BAIL_STEPS);
                }
                self.steps -= 1;
                let g = |x: &Reg| r.get(x.0 as usize).copied().unwrap_or(Slot::F(f64::NAN));
                let value = match inst {
                    | Inst::F(x) => Slot::F(*x),
                    | Inst::I(x) => Slot::I(*x),
                    | Inst::B(x) => Slot::B(*x),
                    | Inst::Copy(a) => g(a),
                    | Inst::FAdd(a, c) => Slot::F(g(a).f() + g(c).f()),
                    | Inst::FSub(a, c) => Slot::F(g(a).f() - g(c).f()),
                    | Inst::FMul(a, c) => Slot::F(g(a).f() * g(c).f()),
                    | Inst::FDiv(a, c) => Slot::F(g(a).f() / g(c).f()),
                    | Inst::FNeg(a) => Slot::F(-g(a).f()),
                    | Inst::FPowI(a, n) => Slot::F(super::tape::int_pow(g(a).f(), *n)),
                    | Inst::FPow(a, c) => Slot::F(g(a).f().powf(g(c).f())),
                    | Inst::FUnary(i, a) => Slot::F(i.apply(g(a).f())),
                    | Inst::FAtan2(a, c) => Slot::F(g(a).f().atan2(g(c).f())),
                    | Inst::FMin(a, c) => Slot::F(fmin(g(a).f(), g(c).f())),
                    | Inst::FMax(a, c) => Slot::F(fmax(g(a).f(), g(c).f())),
                    | Inst::IAdd(a, c) => {
                        match g(a).i().checked_add(g(c).i()) {
                            | Some(x) => Slot::I(x),
                            | None => return Flow::Bail(BAIL_OVERFLOW),
                        }
                    },
                    | Inst::ISub(a, c) => {
                        match g(a).i().checked_sub(g(c).i()) {
                            | Some(x) => Slot::I(x),
                            | None => return Flow::Bail(BAIL_OVERFLOW),
                        }
                    },
                    | Inst::IMul(a, c) => {
                        match g(a).i().checked_mul(g(c).i()) {
                            | Some(x) => Slot::I(x),
                            | None => return Flow::Bail(BAIL_OVERFLOW),
                        }
                    },
                    | Inst::IDiv(a, c) => {
                        match g(a).i().checked_div(g(c).i()) {
                            | Some(x) => Slot::I(x),
                            | None => return Flow::Bail(BAIL_OVERFLOW),
                        }
                    },
                    | Inst::IRem(a, c) => {
                        match g(a).i().checked_rem(g(c).i()) {
                            | Some(x) => Slot::I(x),
                            | None => return Flow::Bail(BAIL_OVERFLOW),
                        }
                    },
                    | Inst::INeg(a) => {
                        match g(a).i().checked_neg() {
                            | Some(x) => Slot::I(x),
                            | None => return Flow::Bail(BAIL_OVERFLOW),
                        }
                    },
                    | Inst::IMin(a, c) => Slot::I(g(a).i().min(g(c).i())),
                    | Inst::IMax(a, c) => Slot::I(g(a).i().max(g(c).i())),
                    | Inst::FCmp(k, a, c) => Slot::B(fcmp(*k, g(a).f(), g(c).f())),
                    | Inst::ICmp(k, a, c) => Slot::B(icmp(*k, g(a).i(), g(c).i())),
                    | Inst::Not(a) => Slot::B(!g(a).b()),
                    | Inst::And(a, c) => Slot::B(g(a).b() && g(c).b()),
                    | Inst::Or(a, c) => Slot::B(g(a).b() || g(c).b()),
                    #[allow(clippy::cast_precision_loss)]
                    | Inst::IToF(a) => Slot::F(g(a).i() as f64),
                    | Inst::FToI(a) => {
                        match f64_to_i64(g(a).f()) {
                            | Some(x) => Slot::I(x),
                            | None => return Flow::Bail(BAIL_CONVERT),
                        }
                    },
                    | Inst::Select(c, a, d) => {
                        if g(c).b() {
                            g(a)
                        } else {
                            g(d)
                        }
                    },
                    | Inst::Len(s) => {
                        match g(s) {
                            | Slot::S(k) => Slot::I(slices.len(k)),
                            | _ => Slot::I(0),
                        }
                    },
                    | Inst::Load(s, i) => {
                        match g(s) {
                            | Slot::S(k) => Slot::F(slices.load(k, g(i).i())),
                            | _ => Slot::F(f64::NAN),
                        }
                    },
                    | Inst::Store(s, i, v) => {
                        if let Slot::S(k) = g(s) {
                            slices.store(k, g(i).i(), g(v).f());
                        }
                        Slot::F(0.0)
                    },
                    | Inst::Call(callee, call_args) => {
                        let values: Vec<Slot> = call_args.iter().map(g).collect();
                        match self.run(*callee, &values, slices, depth.saturating_add(1)) {
                            | Flow::Return(v) => v.first().copied().unwrap_or(Slot::F(f64::NAN)),
                            | Flow::Bail(code) => return Flow::Bail(code),
                        }
                    },
                    | Inst::Eval(function, call_args) => {
                        let values: Vec<f64> = call_args.iter().map(|a| g(a).f()).collect();
                        Slot::F(function(&values))
                    },
                    | Inst::Extern1(function, a) => Slot::F(function(g(a).f())),
                    | Inst::Extern2(function, a, c) => Slot::F(function(g(a).f(), g(c).f())),
                };
                if let Some(slot) = r.get_mut(dst.0 as usize) {
                    *slot = value;
                }
            }
            let read = |regs: &[Reg]| -> Vec<Slot> {
                regs.iter()
                    .map(|x| r.get(x.0 as usize).copied().unwrap_or(Slot::F(f64::NAN)))
                    .collect()
            };
            match &b.term {
                | Term::Jump(t, a) => {
                    incoming = read(a);
                    block = t.0 as usize;
                },
                | Term::Branch(c, (t, ta), (e, ea)) => {
                    let cond = r.get(c.0 as usize).is_some_and(|s| s.b());
                    let (target, a) = if cond { (t, ta) } else { (e, ea) };
                    incoming = read(a);
                    block = target.0 as usize;
                },
                | Term::Return(values) => return Flow::Return(read(values)),
                | Term::Bail(code) => return Flow::Bail(*code),
            }
            if self.steps == 0 {
                return Flow::Bail(BAIL_STEPS);
            }
            self.steps -= 1;
        }
    }
}

// ---------------------------------------------------------------------------
// Lowering from the tape
// ---------------------------------------------------------------------------

/// Lowers a tape to a *batch* function that evaluates it at every row.
///
/// The function's parameters are the input columns (`inputs` slices), the
/// parameter vector (one slice), the output columns (`outputs` slices) and
/// the row count (`I64`): for `row < rows` it stores output `k` of `row` into
/// `out_k[row]`. Missing input elements read NaN, exactly like
/// [`Compiled::call_batch`](super::Compiled::call_batch). Constants and
/// parameters are loop-invariant and evaluated once before the loop.
///
/// Every row gets the bits [`Tape::eval`](super::tape::Tape::eval) gives it:
/// each tape instruction maps to the LIR instruction with the same IEEE
/// semantics. Returns `None` for tapes with `Gather` (its column is chosen
/// per row, which a column-per-parameter signature cannot express).
#[must_use]
pub fn batch_from_tape(tape: &super::tape::Tape) -> Option<Function> {
    use super::tape::Inst as T;
    let outputs = tape.outputs.len();
    let mut sig = vec![Ty::Slice; tape.inputs + 1 + outputs];
    sig.push(Ty::I64);
    let mut f = FunctionBuilder::new("tape_batch", &sig, &[]);
    let column = |k: usize| f.param(k);
    let ins: Vec<Reg> = (0..tape.inputs).map(column).collect();
    let params = f.param(tape.inputs);
    let outs: Vec<Reg> = (0..outputs).map(|k| f.param(tape.inputs + 1 + k)).collect();
    let rows = f.param(tape.inputs + 1 + outputs);
    // Loop-invariant prologue: constants and parameters.
    let mut regs: Vec<Reg> = Vec::with_capacity(tape.insts.len());
    let zero_f = f.ins(Inst::F(0.0));
    let one_f = f.ins(Inst::F(1.0));
    for inst in &tape.insts {
        let r = match *inst {
            | T::Const(c) => f.ins(Inst::F(c)),
            | T::Param(i) => {
                let index = f.ins(Inst::I(i64::from(i)));
                f.ins(Inst::Load(params, index))
            },
            | T::Gather(..) => return None,
            | _ => f.reg(Ty::F64),
        };
        regs.push(r);
    }
    let (head, body, exit) = (f.new_block(), f.new_block(), f.new_block());
    let row = f.block_param(head, Ty::I64);
    let start = f.ins(Inst::I(0));
    f.jump(head, &[start]);
    f.switch_to(head);
    let more = f.ins(Inst::ICmp(super::tape::Cmp::Lt, row, rows));
    f.branch(more, body, &[], exit, &[]);
    f.switch_to(body);
    let reg = |r: u32| regs.get(r as usize).copied().unwrap_or(zero_f);
    for (k, inst) in tape.insts.iter().enumerate() {
        let dst = regs.get(k).copied()?;
        let lowered = match *inst {
            | T::Const(_) | T::Param(_) | T::Gather(..) => continue,
            | T::Input(i) => Inst::Load(*ins.get(i as usize)?, row),
            | T::Add(a, b) => Inst::FAdd(reg(a), reg(b)),
            | T::Sub(a, b) => Inst::FSub(reg(a), reg(b)),
            | T::Mul(a, b) => Inst::FMul(reg(a), reg(b)),
            | T::Div(a, b) => Inst::FDiv(reg(a), reg(b)),
            | T::Neg(a) => Inst::FNeg(reg(a)),
            | T::PowI(a, n) => Inst::FPowI(reg(a), n),
            | T::PowF(a, b) => Inst::FPow(reg(a), reg(b)),
            | T::Unary(i, a) => Inst::FUnary(i, reg(a)),
            | T::Atan2(a, b) => Inst::FAtan2(reg(a), reg(b)),
            | T::Cmp(c, a, b) => {
                // The tape's comparison is the number 1.0 or 0.0.
                let flag = f.ins(Inst::FCmp(c, reg(a), reg(b)));
                Inst::Select(flag, one_f, zero_f)
            },
            | T::Select(c, a, b) => {
                // `c != 0` selects `a`; NaN counts as true.
                let flag = f.ins(Inst::FCmp(super::tape::Cmp::Ne, reg(c), zero_f));
                Inst::Select(flag, reg(a), reg(b))
            },
            | T::Call(function, start, len) => {
                let range = start as usize..(start as usize).saturating_add(len as usize);
                let args = tape.args.get(range)?.iter().map(|&a| reg(a)).collect();
                Inst::Eval(function, args)
            },
            | T::Extern1(function, a) => Inst::Extern1(function, reg(a)),
            | T::Extern2(function, a, b) => Inst::Extern2(function, reg(a), reg(b)),
        };
        f.assign(dst, lowered);
    }
    for (k, &out) in tape.outputs.iter().enumerate() {
        let column = *outs.get(k)?;
        f.ins(Inst::Store(column, row, reg(out)));
    }
    let one = f.ins(Inst::I(1));
    let next = f.ins(Inst::IAdd(row, one));
    f.jump(head, &[next]);
    f.switch_to(exit);
    f.ret(&[]);
    f.finish().ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `fact(n) = if n <= 1 then 1 else n * fact(n - 1)`, recursive.
    fn factorial() -> Module {
        let mut f = FunctionBuilder::new("fact", &[Ty::I64], &[Ty::I64]);
        let n = f.param(0);
        let (base, rec) = (f.new_block(), f.new_block());
        let one = f.ins(Inst::I(1));
        let small = f.ins(Inst::ICmp(Cmp::Le, n, one));
        f.branch(small, base, &[], rec, &[]);
        f.switch_to(base);
        f.ret(&[one]);
        f.switch_to(rec);
        let m = f.ins(Inst::ISub(n, one));
        let r = f.ins_typed(Inst::Call(FuncId(0), vec![m]), Ty::I64);
        let p = f.ins(Inst::IMul(n, r));
        f.ret(&[p]);
        Module {
            functions: vec![f.finish().unwrap_or_else(|e| panic!("{e}"))],
        }
    }

    #[test]
    fn recursion_and_checked_overflow() {
        let m = factorial();
        assert_eq!(m.verify(), Ok(()));
        let out = m.interpret(FuncId(0), &mut [Value::I(20)], 1 << 20);
        assert_eq!(
            out,
            Ok(Outcome::Returned(vec![Value::I(2_432_902_008_176_640_000)]))
        );
        // 21! overflows i64: the call bails instead of wrapping.
        let out = m.interpret(FuncId(0), &mut [Value::I(21)], 1 << 20);
        assert_eq!(out, Ok(Outcome::Bailed(BAIL_OVERFLOW)));
        // Unbounded recursion hits the depth limit.
        let out = m.interpret(FuncId(0), &mut [Value::I(100_000)], u64::MAX);
        assert_eq!(out, Ok(Outcome::Bailed(BAIL_DEPTH)));
    }

    #[test]
    fn slices_are_bounds_checked_and_read_only_slices_are_refused() {
        // y[i] = 2 * x[i] for i in 0..len(y), reading one past the end of x.
        let mut f = FunctionBuilder::new("scale", &[Ty::Slice, Ty::Slice], &[]);
        let (x, y) = (f.param(0), f.param(1));
        let (head, body, exit) = (f.new_block(), f.new_block(), f.new_block());
        let i = f.block_param(head, Ty::I64);
        let zero = f.ins(Inst::I(0));
        f.jump(head, &[zero]);
        f.switch_to(head);
        let n = f.ins(Inst::Len(y));
        let more = f.ins(Inst::ICmp(Cmp::Lt, i, n));
        f.branch(more, body, &[], exit, &[]);
        f.switch_to(body);
        let v = f.ins(Inst::Load(x, i));
        let two = f.ins(Inst::F(2.0));
        let w = f.ins(Inst::FMul(v, two));
        f.ins(Inst::Store(y, i, w));
        let one = f.ins(Inst::I(1));
        let next = f.ins(Inst::IAdd(i, one));
        f.jump(head, &[next]);
        f.switch_to(exit);
        f.ret(&[]);
        let m = Module {
            functions: vec![f.finish().unwrap_or_else(|e| panic!("{e}"))],
        };
        let x = [1.0, 2.0];
        let mut y = [0.0; 3];
        let out = m.interpret(
            FuncId(0),
            &mut [Value::Read(&x), Value::Write(&mut y)],
            1000,
        );
        assert_eq!(out, Ok(Outcome::Returned(vec![])));
        assert_eq!(y[..2], [2.0, 4.0]);
        assert!(y[2].is_nan());
        let mut z = [0.0; 1];
        // The output must be writable; the input may be read-only.
        let refused = m.interpret(FuncId(0), &mut [Value::Read(&x), Value::Read(&[0.0])], 1000);
        assert!(refused.is_err());
        assert_eq!(m.written_params(FuncId(0)), vec![false, true]);
        let ok = m.interpret(
            FuncId(0),
            &mut [Value::Read(&x), Value::Write(&mut z)],
            1000,
        );
        assert!(ok.is_ok());
    }

    #[test]
    fn tape_batches_agree_with_the_scalar_tape() {
        use crate::graph::Graph;
        let mut g = Graph::new();
        assert!(crate::graph::Engine::install(&mut g, &crate::rules::standard()).is_ok());
        let root = g
            .parse("sin(x)*exp(-x^2/2) + x^3 - 2*x/(1 + x^2) + sqrt(x + 1) + erf(y)")
            .unwrap_or_else(|e| panic!("{e}"));
        let x = g.interner_mut().symbol("x");
        let y = g.interner_mut().symbol("y");
        let tape = super::super::tape::lower(&g, &[root], &[x, y], &[])
            .unwrap_or_else(|e| panic!("{e}"))
            .optimise();
        let batch = batch_from_tape(&tape).unwrap_or_else(|| panic!("lowers"));
        let module = Module {
            functions: vec![batch],
        };
        let xs: Vec<f64> = (0..37).map(|i| f64::from(i) * 0.17 - 2.0).collect();
        let ys: Vec<f64> = (0..30).map(|i| f64::from(i) * -0.3).collect();
        let mut out = vec![0.0; 37];
        let outcome = module.interpret(
            FuncId(0),
            &mut [
                Value::Read(&xs),
                Value::Read(&ys),
                Value::Read(&[]),
                Value::Write(&mut out),
                Value::I(37),
            ],
            1 << 24,
        );
        assert_eq!(outcome, Ok(Outcome::Returned(vec![])));
        let mut columns = vec![0.0; 37];
        tape.eval_columns(&[&xs, &ys], &[], &mut [&mut columns]);
        for (row, &got) in out.iter().enumerate() {
            let mut want = [0.0];
            let point = [xs[row], ys.get(row).copied().unwrap_or(f64::NAN)];
            tape.eval(&point, &[], &mut want);
            assert_eq!(got.to_bits(), want[0].to_bits(), "row {row}");
            assert_eq!(
                columns[row].to_bits(),
                want[0].to_bits(),
                "columns, row {row}"
            );
        }
    }

    #[test]
    fn verification_catches_type_errors() {
        let mut f = FunctionBuilder::new("bad", &[Ty::F64], &[Ty::F64]);
        let x = f.param(0);
        let y = f.ins_typed(Inst::IAdd(x, x), Ty::I64);
        f.ret(&[y]);
        let m = Module {
            functions: vec![f.finish().unwrap_or_else(|e| panic!("{e}"))],
        };
        assert!(m.verify().is_err());
    }

    #[test]
    fn registers_start_at_zero_and_jumps_assign_in_parallel() {
        // b1(a, b) swaps: jump b1(b, a) once, then return a - b.
        let mut f = FunctionBuilder::new("swap", &[Ty::F64, Ty::F64], &[Ty::F64]);
        let (x, y) = (f.param(0), f.param(1));
        let (b1, b2) = (f.new_block(), f.new_block());
        let a = f.block_param(b1, Ty::F64);
        let b = f.block_param(b1, Ty::F64);
        let unset = f.reg(Ty::F64);
        f.jump(b1, &[x, y]);
        f.switch_to(b1);
        f.jump(b2, &[]);
        f.switch_to(b2);
        let d = f.ins(Inst::FSub(b, a));
        let e = f.ins(Inst::FAdd(d, unset));
        f.ret(&[e]);
        let m = Module {
            functions: vec![f.finish().unwrap_or_else(|e| panic!("{e}"))],
        };
        let out = m.interpret(FuncId(0), &mut [Value::F(1.0), Value::F(5.0)], 100);
        assert_eq!(out, Ok(Outcome::Returned(vec![Value::F(4.0)])));
    }
}
