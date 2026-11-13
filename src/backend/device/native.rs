//! Native CPU execution of device kernels through LIR (feature `jit`).
//!
//! The CPU devices run a [`KernelIr`] by interpreting it row by row. With the
//! `jit` feature they first try to lower it to a [LIR](crate::backend::lir)
//! function — the row loop, the fold loop and every instruction — and run
//! that as machine code, splitting the rows over the thread pool exactly as
//! the interpreter does (every row is independent, so the split changes
//! nothing).
//!
//! The lowered function is
//!
//! ```text
//! kernel(col₀ … col_{i-1}: Slice, params: Slice, out₀ … out_{k-1}: Slice,
//!        row_lo: I64, rows: I64)
//! ```
//!
//! and computes rows `row_lo .. row_lo + len(out₀)`, writing `out_k[at]` for
//! row `row_lo + at`. The device semantics that differ from LIR's are emulated
//! instruction by instruction: a comparison is the number `1.0`/`0.0`, a
//! `select` takes its first operand for any non-zero or NaN condition, index
//! division and remainder by zero give `0`, `ToIdx` is NaN-to-zero and
//! saturates at the 32-bit limits. Index `+ - *` wrap in the interpreter and
//! are checked in LIR: an overflow *bails* and the launch falls back to the
//! interpreter, so results never differ. The one deliberate difference is the
//! sign of a zero chosen by a `Min`/`Max` fold between `-0.0` and `+0.0`,
//! which the interpreter leaves to `f64::min`.
//!
//! `Gather` (an input column chosen per row) has no LIR form; such kernels are
//! interpreted.

use super::ReduceOp;
use super::kernel::Analysis;
use super::kernel::KernelIr;
use super::kernel::KernelOp;
use super::kernel::Ty as KTy;
use crate::backend::lir::FunctionBuilder;
use crate::backend::lir::Inst as L;
use crate::backend::lir::Module;
use crate::backend::lir::Reg;
use crate::backend::lir::Ty;
use crate::backend::lir_jit::LirModule;
use crate::backend::tape::Cmp;

/// Lowers `ir` (with its analysis) to the LIR kernel described in the module
/// documentation; `None` for kernels LIR cannot express.
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn lower(
    ir: &KernelIr,
    analysis: &Analysis,
) -> Option<Module> {
    let outputs = ir.outputs.len();
    let mut sig = vec![Ty::Slice; ir.inputs + 1 + outputs];
    sig.extend([Ty::I64, Ty::I64]);
    let mut f = FunctionBuilder::new("device_kernel", &sig, &[]);
    let cols: Vec<Reg> = (0..ir.inputs).map(|k| f.param(k)).collect();
    let params = f.param(ir.inputs);
    let outs: Vec<Reg> = (0..outputs).map(|k| f.param(ir.inputs + 1 + k)).collect();
    let row_lo = f.param(ir.inputs + 1 + outputs);
    let rows = f.param(ir.inputs + 2 + outputs);
    // One LIR register per kernel register, of the analysed type.
    let regs: Vec<Reg> = analysis
        .ty
        .iter()
        .map(|t| {
            f.reg(if *t == KTy::Idx {
                Ty::I64
            } else {
                Ty::F64
            })
        })
        .collect();
    let zero = f.ins(L::F(0.0));
    let one = f.ins(L::F(1.0));
    let izero = f.ins(L::I(0));
    let ione = f.ins(L::I(1));
    let count = match outs.first() {
        | Some(&o) => f.ins(L::Len(o)),
        | None => izero,
    };
    let (head, body, exit) = (f.new_block(), f.new_block(), f.new_block());
    let at = f.block_param(head, Ty::I64);
    f.jump(head, &[izero]);
    f.switch_to(head);
    let more = f.ins(L::ICmp(Cmp::Lt, at, count));
    f.branch(more, body, &[], exit, &[]);
    f.switch_to(body);
    let row = f.ins(L::IAdd(row_lo, at));
    let inner = f.reg(Ty::I64);
    let accs: Vec<Reg> = ir
        .fold
        .as_ref()
        .map(|fold| fold.accs.iter().map(|_| f.reg(Ty::F64)).collect())
        .unwrap_or_default();

    // Emits register `d`.
    let emit = |f: &mut FunctionBuilder, d: usize| -> Option<()> {
        let r = |x: u32| regs.get(x as usize).copied();
        let dst = *regs.get(d)?;
        let flag = |f: &mut FunctionBuilder, c: Reg| -> Reg {
            // select_is_true: NaN or non-zero.
            f.ins(L::FCmp(Cmp::Ne, c, zero))
        };
        let inst = match *ir.ops.get(d)? {
            | KernelOp::Const(c) => L::F(c),
            | KernelOp::Input(k) => L::Load(*cols.get(k as usize)?, row),
            | KernelOp::Param(k) => {
                let i = f.ins(L::I(i64::from(k)));
                L::Load(params, i)
            },
            | KernelOp::Add(a, b) => L::FAdd(r(a)?, r(b)?),
            | KernelOp::Sub(a, b) => L::FSub(r(a)?, r(b)?),
            | KernelOp::Mul(a, b) => L::FMul(r(a)?, r(b)?),
            | KernelOp::Div(a, b) => L::FDiv(r(a)?, r(b)?),
            | KernelOp::Neg(a) => L::FNeg(r(a)?),
            | KernelOp::PowI(a, n) => L::FPowI(r(a)?, n),
            | KernelOp::PowF(a, b) => L::FPow(r(a)?, r(b)?),
            | KernelOp::Unary(i, a) => L::FUnary(i, r(a)?),
            | KernelOp::Atan2(a, b) => L::FAtan2(r(a)?, r(b)?),
            | KernelOp::Cmp(c, a, b) => {
                let t = f.ins(L::FCmp(c, r(a)?, r(b)?));
                L::Select(t, one, zero)
            },
            | KernelOp::ICmp(c, a, b) => {
                let t = f.ins(L::ICmp(c, r(a)?, r(b)?));
                L::Select(t, one, zero)
            },
            | KernelOp::Select(c, a, b) => {
                let t = flag(f, r(c)?);
                L::Select(t, r(a)?, r(b)?)
            },
            | KernelOp::Gather(..) => return None,
            | KernelOp::Row => L::Copy(row),
            | KernelOp::Rows => L::Copy(rows),
            | KernelOp::Len(k) => L::Len(*cols.get(k as usize)?),
            | KernelOp::Inner => L::Copy(inner),
            | KernelOp::IConst(c) => L::I(c),
            | KernelOp::IAdd(a, b) => L::IAdd(r(a)?, r(b)?),
            | KernelOp::ISub(a, b) => L::ISub(r(a)?, r(b)?),
            | KernelOp::IMul(a, b) => L::IMul(r(a)?, r(b)?),
            | KernelOp::IDiv(a, b) | KernelOp::IMod(a, b) => {
                // `0` for a zero divisor and for MIN / -1, like checked_div.
                let (x, y) = (r(a)?, r(b)?);
                let zero_div = f.ins(L::ICmp(Cmp::Eq, y, izero));
                let min = f.ins(L::I(i64::MIN));
                let minus = f.ins(L::I(-1));
                let is_min = f.ins(L::ICmp(Cmp::Eq, x, min));
                let is_minus = f.ins(L::ICmp(Cmp::Eq, y, minus));
                let wrap = f.ins(L::And(is_min, is_minus));
                let bad = f.ins(L::Or(zero_div, wrap));
                let safe = f.ins(L::Select(bad, ione, y));
                let q = if matches!(ir.ops.get(d)?, KernelOp::IDiv(..)) {
                    f.ins(L::IDiv(x, safe))
                } else {
                    f.ins(L::IRem(x, safe))
                };
                L::Select(bad, izero, q)
            },
            | KernelOp::IMin(a, b) => L::IMin(r(a)?, r(b)?),
            | KernelOp::IMax(a, b) => L::IMax(r(a)?, r(b)?),
            | KernelOp::ToNum(a) => L::IToF(r(a)?),
            | KernelOp::ToIdx(a) => {
                // NaN -> 0, then clamp to the 32-bit range and truncate.
                let x = r(a)?;
                let nan = f.ins(L::FCmp(Cmp::Ne, x, x));
                let x = f.ins(L::Select(nan, zero, x));
                let lo = f.ins(L::F(f64::from(i32::MIN)));
                let hi = f.ins(L::F(f64::from(i32::MAX)));
                let x = f.ins(L::FMax(x, lo));
                let x = f.ins(L::FMin(x, hi));
                L::FToI(x)
            },
            | KernelOp::Load(k, i) => L::Load(*cols.get(k as usize)?, r(i)?),
            | KernelOp::Acc(k) => L::Copy(*accs.get(k as usize)?),
        };
        f.assign(dst, inst);
        Some(())
    };

    for &d in &analysis.pre {
        emit(&mut f, d)?;
    }
    if let Some(fold) = &ir.fold {
        for (&acc, (op, _)) in accs.iter().zip(&fold.accs) {
            f.assign(acc, L::F(op.identity()));
        }
        let total = *regs.get(fold.count as usize)?;
        let (fhead, fbody, fexit) = (f.new_block(), f.new_block(), f.new_block());
        let k = f.block_param(fhead, Ty::I64);
        f.jump(fhead, &[izero]);
        f.switch_to(fhead);
        let go = f.ins(L::ICmp(Cmp::Lt, k, total));
        f.branch(go, fbody, &[], fexit, &[]);
        f.switch_to(fbody);
        f.assign(inner, L::Copy(k));
        for &d in &analysis.body {
            emit(&mut f, d)?;
        }
        for (&acc, &(op, reg)) in accs.iter().zip(&fold.accs) {
            let v = *regs.get(reg as usize)?;
            match op {
                | ReduceOp::Sum => f.assign(acc, L::FAdd(acc, v)),
                | ReduceOp::Prod => f.assign(acc, L::FMul(acc, v)),
                | ReduceOp::Min | ReduceOp::Max => {
                    // NaN if either is NaN, else the smaller (larger).
                    let c = if op == ReduceOp::Min {
                        Cmp::Lt
                    } else {
                        Cmp::Gt
                    };
                    let better = f.ins(L::FCmp(c, v, acc));
                    let picked = f.ins(L::Select(better, v, acc));
                    let nan_a = f.ins(L::FCmp(Cmp::Ne, acc, acc));
                    let nan_v = f.ins(L::FCmp(Cmp::Ne, v, v));
                    let any = f.ins(L::Or(nan_a, nan_v));
                    let nan = f.ins(L::F(f64::NAN));
                    f.assign(acc, L::Select(any, nan, picked));
                },
            }
        }
        let next = f.ins(L::IAdd(k, ione));
        f.jump(fhead, &[next]);
        f.switch_to(fexit);
        for &d in &analysis.post {
            emit(&mut f, d)?;
        }
    }
    for (&out, &reg) in outs.iter().zip(&ir.outputs) {
        let v = *regs.get(reg as usize)?;
        f.ins(L::Store(out, at, v));
    }
    let next = f.ins(L::IAdd(at, ione));
    f.jump(head, &[next]);
    f.switch_to(exit);
    f.ret(&[]);
    Some(Module {
        functions: vec![f.finish().ok()?],
    })
}

/// A kernel compiled to machine code.
#[derive(Debug)]
pub struct NativeKernel {
    module: LirModule,
    inputs: usize,
}

impl NativeKernel {
    /// Lowers and compiles `ir`; `None` when either step fails.
    #[must_use]
    pub fn compile(
        ir: &KernelIr,
        analysis: &Analysis,
    ) -> Option<Self> {
        let module = lower(ir, analysis)?;
        Some(Self {
            module: LirModule::compile(&module).ok()?,
            inputs: ir.inputs,
        })
    }

    /// Runs rows `row_lo .. row_lo + out[0].len()`; `false` when the code
    /// bailed (the caller then interprets the launch).
    pub fn run(
        &self,
        cols: &[&[f64]],
        params: &[f64],
        out: &mut [&mut [f64]],
        row_lo: usize,
        rows: usize,
    ) -> bool {
        use crate::backend::lir::Value;
        if cols.len() != self.inputs {
            return false;
        }
        let (Ok(lo), Ok(rows)) = (i64::try_from(row_lo), i64::try_from(rows)) else {
            return false;
        };
        let mut args: Vec<Value<'_>> = cols.iter().map(|c| Value::Read(c)).collect();
        args.push(Value::Read(params));
        for o in out.iter_mut() {
            args.push(Value::Write(o));
        }
        args.push(Value::I(lo));
        args.push(Value::I(rows));
        matches!(
            self.module.call(crate::backend::lir::FuncId(0), &mut args),
            Ok(crate::backend::lir::Outcome::Returned(_))
        )
    }
}
