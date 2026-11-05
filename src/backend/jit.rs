//! # Cranelift JIT backend
//!
//! Compiles the [Tape IR](super::tape) to machine code with Cranelift
//! (feature `jit`). Every program becomes one function
//! `fn(inputs: *const f64, params: *const f64, out: *mut f64)`; all of its
//! outputs share their common subexpressions.
//!
//! * Arithmetic, `sqrt`, `abs`, `floor`, `ceil`, comparisons and selects
//!   are native instructions; integer powers are square-and-multiply
//!   chains (the same sequence as the tape interpreter, so results agree
//!   bit for bit); other elementary functions call the platform's math
//!   library through `extern "C"` shims; operators known only by their
//!   [`EvalFn`](crate::graph::op::EvalFn) are called through a shim with
//!   their arguments in a stack slot.
//! * A compiled function owns its code through an [`Arc`]: the machine
//!   code is freed when the last handle is dropped. There is no global
//!   state; [`CraneliftBackend`] keeps a cache keyed by the structural
//!   hash of the optimised tape, so the same function rebuilt in another
//!   graph is compiled once.
//! * [`TieredBackend`] starts with the tape interpreter (no compile
//!   latency) and switches to JIT code compiled on a background thread
//!   once a function has been called often enough.
//!
//! All `unsafe` code lives in the private `abi` module: the calling
//! convention of the generated code and the shims it calls.

use std::collections::HashMap;
use std::sync::Arc;
use std::sync::Mutex;
use std::sync::OnceLock;
use std::sync::atomic::AtomicBool;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;

use cranelift_codegen::ir::AbiParam;
use cranelift_codegen::ir::InstBuilder;
use cranelift_codegen::ir::MemFlagsData;
use cranelift_codegen::ir::StackSlotData;
use cranelift_codegen::ir::StackSlotKind;
use cranelift_codegen::ir::Value;
use cranelift_codegen::ir::condcodes::FloatCC;
use cranelift_codegen::ir::condcodes::IntCC;
use cranelift_codegen::ir::types;
use cranelift_codegen::settings;
use cranelift_codegen::settings::Configurable;
use cranelift_frontend::FunctionBuilder;
use cranelift_frontend::FunctionBuilderContext;
use cranelift_jit::JITBuilder;
use cranelift_jit::JITModule;
use cranelift_module::FuncId;
use cranelift_module::Linkage;
use cranelift_module::Module;

use super::Backend;
use super::BackendError;
use super::Compiled;
use super::CompiledMulti;
use super::Intrinsic;
use super::tape::Cmp;
use super::tape::Inst;
use super::tape::Tape;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::SymbolId;

#[allow(unsafe_code)]
mod abi {
    //! The calling convention of generated code and the shims it calls.

    use crate::graph::op::EvalFn;

    /// Signature of every generated function.
    pub type Entry = unsafe extern "C" fn(*const f64, *const f64, *mut f64);

    /// Calls an operator's scalar semantics from generated code.
    ///
    /// # Safety
    /// `f` must be the address of an [`EvalFn`] and `args` must point to
    /// `len` readable `f64`s; generated code guarantees both.
    pub unsafe extern "C" fn eval_shim(
        f: usize,
        args: *const f64,
        len: usize,
    ) -> f64 {
        // SAFETY: `f` was produced from an `EvalFn` by the code generator.
        let function: EvalFn = unsafe { std::mem::transmute::<usize, EvalFn>(f) };
        // SAFETY: the caller passes a stack slot holding `len` values.
        let slice = unsafe { std::slice::from_raw_parts(args, len) };
        function(slice)
    }

    /// Calls a one-argument C function given by address.
    ///
    /// # Safety
    /// `f` must be the address of an `extern "C" fn(f64) -> f64`.
    pub unsafe extern "C" fn extern1_shim(
        f: usize,
        x: f64,
    ) -> f64 {
        // SAFETY: the address comes from a `Lowering::Extern1`.
        let function: extern "C" fn(f64) -> f64 =
            unsafe { std::mem::transmute::<usize, extern "C" fn(f64) -> f64>(f) };
        function(x)
    }

    /// Calls a two-argument C function given by address.
    ///
    /// # Safety
    /// `f` must be the address of an `extern "C" fn(f64, f64) -> f64`.
    pub unsafe extern "C" fn extern2_shim(
        f: usize,
        x: f64,
        y: f64,
    ) -> f64 {
        // SAFETY: the address comes from a `Lowering::Extern2`.
        let function: extern "C" fn(f64, f64) -> f64 =
            unsafe { std::mem::transmute::<usize, extern "C" fn(f64, f64) -> f64>(f) };
        function(x, y)
    }

    /// Runs generated code.
    ///
    /// # Safety
    /// `entry` must be a finalized function of the live module that
    /// produced it; `inputs`/`params` must hold at least as many values as
    /// the program reads and `out` room for every output.
    pub unsafe fn call(
        entry: Entry,
        inputs: &[f64],
        params: &[f64],
        out: &mut [f64],
    ) {
        // SAFETY: guaranteed by the caller (lengths are checked there).
        unsafe { entry(inputs.as_ptr(), params.as_ptr(), out.as_mut_ptr()) }
    }

    /// Reinterprets a finalized code pointer as an [`Entry`].
    ///
    /// # Safety
    /// `ptr` must point to a function with the [`Entry`] signature.
    pub unsafe fn entry(ptr: *const u8) -> Entry {
        // SAFETY: guaranteed by the caller.
        unsafe { std::mem::transmute::<*const u8, Entry>(ptr) }
    }

    /// Frees the code of a module.
    ///
    /// # Safety
    /// No function of the module may be called afterwards.
    pub unsafe fn free(module: cranelift_jit::JITModule) {
        // SAFETY: guaranteed by the caller (the last owner drops it).
        unsafe { module.free_memory() }
    }

    /// Owner of a module's code; frees it when dropped.
    pub struct Owner(pub Option<cranelift_jit::JITModule>);

    // SAFETY: after finalisation the module is never mutated; its code is
    // read-only and may be executed from any thread. The owner only frees
    // it on drop, when no other handle exists.
    unsafe impl Send for Owner {}
    // SAFETY: see above; shared references never touch the module.
    unsafe impl Sync for Owner {}

    impl Drop for Owner {
        fn drop(&mut self) {
            if let Some(module) = self.0.take() {
                // SAFETY: this is the last owner; no entry point survives it.
                unsafe { free(module) };
            }
        }
    }
}

macro_rules! math_shims {
    ($($name:ident => $f:expr),* $(,)?) => {
        $(
            extern "C" fn $name(x: f64) -> f64 {
                let f: fn(f64) -> f64 = $f;
                f(x)
            }
        )*
    };
}

math_shims! {
    shim_exp => f64::exp,
    shim_ln => f64::ln,
    shim_sin => f64::sin,
    shim_cos => f64::cos,
    shim_tan => f64::tan,
    shim_asin => f64::asin,
    shim_acos => f64::acos,
    shim_atan => f64::atan,
    shim_sinh => f64::sinh,
    shim_cosh => f64::cosh,
    shim_tanh => f64::tanh,
    shim_asinh => f64::asinh,
    shim_acosh => f64::acosh,
    shim_atanh => f64::atanh,
    shim_cbrt => f64::cbrt,
}

extern "C" fn shim_pow(
    x: f64,
    y: f64,
) -> f64 {
    x.powf(y)
}

extern "C" fn shim_atan2(
    y: f64,
    x: f64,
) -> f64 {
    y.atan2(x)
}

fn library(f: Intrinsic) -> Option<(&'static str, extern "C" fn(f64) -> f64)> {
    Some(match f {
        | Intrinsic::Exp => ("rssn_exp", shim_exp),
        | Intrinsic::Ln => ("rssn_ln", shim_ln),
        | Intrinsic::Sin => ("rssn_sin", shim_sin),
        | Intrinsic::Cos => ("rssn_cos", shim_cos),
        | Intrinsic::Tan => ("rssn_tan", shim_tan),
        | Intrinsic::Asin => ("rssn_asin", shim_asin),
        | Intrinsic::Acos => ("rssn_acos", shim_acos),
        | Intrinsic::Atan => ("rssn_atan", shim_atan),
        | Intrinsic::Sinh => ("rssn_sinh", shim_sinh),
        | Intrinsic::Cosh => ("rssn_cosh", shim_cosh),
        | Intrinsic::Tanh => ("rssn_tanh", shim_tanh),
        | Intrinsic::Asinh => ("rssn_asinh", shim_asinh),
        | Intrinsic::Acosh => ("rssn_acosh", shim_acosh),
        | Intrinsic::Atanh => ("rssn_atanh", shim_atanh),
        | Intrinsic::Cbrt => ("rssn_cbrt", shim_cbrt),
        | Intrinsic::Sqrt | Intrinsic::Abs | Intrinsic::Floor | Intrinsic::Ceil => return None,
    })
}

const UNARY_SHIMS: [Intrinsic; 15] = [
    Intrinsic::Exp,
    Intrinsic::Ln,
    Intrinsic::Sin,
    Intrinsic::Cos,
    Intrinsic::Tan,
    Intrinsic::Asin,
    Intrinsic::Acos,
    Intrinsic::Atan,
    Intrinsic::Sinh,
    Intrinsic::Cosh,
    Intrinsic::Tanh,
    Intrinsic::Asinh,
    Intrinsic::Acosh,
    Intrinsic::Atanh,
    Intrinsic::Cbrt,
];

fn error(e: impl std::fmt::Display) -> BackendError {
    BackendError::Codegen(e.to_string())
}

// Cranelift memory offsets are signed 32-bit bytes, not element indices.
// Saturating here silently loads the wrong element for large valid tapes.
fn f64_byte_offset(index: u32) -> Result<i32, BackendError> {
    i32::try_from(index)
        .ok()
        .and_then(|value| value.checked_mul(8))
        .ok_or_else(|| BackendError::Codegen("floating-point load offset exceeds i32".to_owned()))
}

/// A program compiled to machine code.
pub struct JitFunction {
    _owner: Arc<abi::Owner>,
    entry: abi::Entry,
    inputs: usize,
    params: usize,
    outputs: usize,
}

impl JitFunction {
    /// Compiles `tape` (which should already be optimised).
    ///
    /// # Errors
    /// [`BackendError::Codegen`] when Cranelift rejects the program or the
    /// host is unsupported.
    pub fn compile(tape: &Tape) -> Result<Self, BackendError> {
        Self::compile_impl(tape, false)
    }

    fn compile_impl(tape: &Tape, tail_loop: bool) -> Result<Self, BackendError> {
        let mut flags = settings::builder();
        flags
            .set("use_colocated_libcalls", "false")
            .map_err(error)?;
        flags.set("is_pic", "false").map_err(error)?;
        flags.set("opt_level", "speed").map_err(error)?;
        let isa = cranelift_native::builder()
            .map_err(error)?
            .finish(settings::Flags::new(flags))
            .map_err(error)?;
        // `finalize` needs the target's frontend config, and `isa` is moved
        // into the builder below, so take it first.
        let frontend_config = isa.frontend_config();
        let mut builder = JITBuilder::with_isa(isa, cranelift_module::default_libcall_names());
        for f in UNARY_SHIMS {
            if let Some((name, shim)) = library(f) {
                builder.symbol(name, shim as *const u8);
            }
        }
        builder.symbol("rssn_pow", shim_pow as *const u8);
        builder.symbol("rssn_atan2", shim_atan2 as *const u8);
        builder.symbol("rssn_eval", abi::eval_shim as *const u8);
        builder.symbol("rssn_extern1", abi::extern1_shim as *const u8);
        builder.symbol("rssn_extern2", abi::extern2_shim as *const u8);
        let mut module = JITModule::new(builder);
        let pointer = module.target_config().pointer_type();

        // Imported shims.
        let mut sig_unary = module.make_signature();
        sig_unary.params.push(AbiParam::new(types::F64));
        sig_unary.returns.push(AbiParam::new(types::F64));
        let mut sig_binary = sig_unary.clone();
        sig_binary.params.push(AbiParam::new(types::F64));
        let mut sig_eval = module.make_signature();
        sig_eval.params.extend([
            AbiParam::new(pointer),
            AbiParam::new(pointer),
            AbiParam::new(pointer),
        ]);
        sig_eval.returns.push(AbiParam::new(types::F64));
        let mut sig_ext1 = module.make_signature();
        sig_ext1
            .params
            .extend([AbiParam::new(pointer), AbiParam::new(types::F64)]);
        sig_ext1.returns.push(AbiParam::new(types::F64));
        let mut sig_ext2 = sig_ext1.clone();
        sig_ext2.params.push(AbiParam::new(types::F64));
        let mut imports: HashMap<&'static str, FuncId> = HashMap::new();
        for f in UNARY_SHIMS {
            if let Some((name, _)) = library(f) {
                imports.insert(
                    name,
                    module
                        .declare_function(name, Linkage::Import, &sig_unary)
                        .map_err(error)?,
                );
            }
        }
        imports.insert(
            "rssn_pow",
            module
                .declare_function("rssn_pow", Linkage::Import, &sig_binary)
                .map_err(error)?,
        );
        imports.insert(
            "rssn_atan2",
            module
                .declare_function("rssn_atan2", Linkage::Import, &sig_binary)
                .map_err(error)?,
        );
        imports.insert(
            "rssn_eval",
            module
                .declare_function("rssn_eval", Linkage::Import, &sig_eval)
                .map_err(error)?,
        );
        imports.insert(
            "rssn_extern1",
            module
                .declare_function("rssn_extern1", Linkage::Import, &sig_ext1)
                .map_err(error)?,
        );
        imports.insert(
            "rssn_extern2",
            module
                .declare_function("rssn_extern2", Linkage::Import, &sig_ext2)
                .map_err(error)?,
        );

        let mut sig = module.make_signature();
        sig.params.extend([
            AbiParam::new(pointer),
            AbiParam::new(pointer),
            AbiParam::new(pointer),
        ]);
        let id = module
            .declare_function("rssn_program", Linkage::Local, &sig)
            .map_err(error)?;
        let mut ctx = module.make_context();
        ctx.func.signature = sig;
        let mut fctx = FunctionBuilderContext::new();
        {
            let mut b = FunctionBuilder::new(&mut ctx.func, &mut fctx);
            let block = b.create_block();
            b.append_block_params_for_function_params(block);
            b.switch_to_block(block);
            b.seal_block(block);
            let params = b.block_params(block).to_vec();
            let &[in_ptr, par_ptr, out_ptr] = params.as_slice() else {
                return Err(BackendError::Codegen("signature".to_owned()));
            };
            // 0.136 separates the flag *data* (what the builder takes) from
            // the interned `MemFlags` handle the encoder uses.
            let flags = MemFlagsData::trusted();
            let mut refs = HashMap::new();
            for (&name, &fid) in &imports {
                refs.insert(name, module.declare_func_in_func(fid, b.func));
            }
            let call1 = |b: &mut FunctionBuilder<'_>,
                         name: &str,
                         args: &[Value]|
             -> Result<Value, BackendError> {
                let f = *refs
                    .get(name)
                    .ok_or_else(|| BackendError::Codegen(name.to_owned()))?;
                let call = b.ins().call(f, args);
                b.inst_results(call)
                    .first()
                    .copied()
                    .ok_or_else(|| BackendError::Codegen(name.to_owned()))
            };
            // Tail loops use `out[2..]` as a mutable parameter frame. The
            // ordinary Tape lowering below runs unchanged inside this block.
            // The extra final output carries the unconsumed u64 fuel bits.
            let loop_blocks = if tail_loop {
                for i in 0..tape.inputs {
                    let input_offset = f64_byte_offset(u32::try_from(i).map_err(error)?)?;
                    let output_offset = f64_byte_offset(u32::try_from(i + 2).map_err(error)?)?;
                    let value = b.ins().load(types::F64, flags, in_ptr, input_offset);
                    b.ins().store(flags, value, out_ptr, output_offset);
                }
                let initial_fuel = b.ins().load(types::I64, flags, par_ptr, 0);
                let header = b.create_block();
                let more = b.create_block();
                let exit = b.create_block();
                b.append_block_param(header, types::I64);
                b.append_block_param(exit, types::I64);
                b.ins().jump(header, &[initial_fuel.into()]);
                b.switch_to_block(header);
                Some((header, more, exit))
            } else {
                None
            };
            let tape_in_ptr = if tail_loop { b.ins().iadd_imm_s(out_ptr, 16) } else { in_ptr };
            let mut values: Vec<Value> = Vec::with_capacity(tape.insts.len());
            let get = |values: &[Value], r: u32| {
                values
                    .get(r as usize)
                    .copied()
                    .ok_or_else(|| BackendError::Codegen("register".to_owned()))
            };
            for inst in &tape.insts {
                let v = match *inst {
                    | Inst::Const(c) => b.ins().f64const(c),
                    | Inst::Input(i) => {
                        b.ins().load(types::F64, flags, tape_in_ptr, f64_byte_offset(i)?)
                    },
                    | Inst::Param(i) => {
                        b.ins()
                            .load(types::F64, flags, par_ptr, f64_byte_offset(i)?)
                    },
                    | Inst::Add(x, y) => {
                        let (x, y) = (get(&values, x)?, get(&values, y)?);
                        b.ins().fadd(x, y)
                    },
                    | Inst::Sub(x, y) => {
                        let (x, y) = (get(&values, x)?, get(&values, y)?);
                        b.ins().fsub(x, y)
                    },
                    | Inst::Mul(x, y) => {
                        let (x, y) = (get(&values, x)?, get(&values, y)?);
                        b.ins().fmul(x, y)
                    },
                    | Inst::Div(x, y) => {
                        let (x, y) = (get(&values, x)?, get(&values, y)?);
                        b.ins().fdiv(x, y)
                    },
                    | Inst::Neg(x) => {
                        let x = get(&values, x)?;
                        b.ins().fneg(x)
                    },
                    | Inst::PowI(x, n) => {
                        let x = get(&values, x)?;
                        emit_powi(&mut b, x, n)
                    },
                    | Inst::PowF(x, y) => {
                        let (x, y) = (get(&values, x)?, get(&values, y)?);
                        call1(&mut b, "rssn_pow", &[x, y])?
                    },
                    | Inst::Atan2(y, x) => {
                        let (y, x) = (get(&values, y)?, get(&values, x)?);
                        call1(&mut b, "rssn_atan2", &[y, x])?
                    },
                    | Inst::Unary(f, x) => {
                        let x = get(&values, x)?;
                        match f {
                            | Intrinsic::Sqrt => b.ins().sqrt(x),
                            | Intrinsic::Abs => b.ins().fabs(x),
                            | Intrinsic::Floor => b.ins().floor(x),
                            | Intrinsic::Ceil => b.ins().ceil(x),
                            | other => {
                                let (name, _) = library(other)
                                    .ok_or_else(|| BackendError::Codegen("intrinsic".to_owned()))?;
                                call1(&mut b, name, &[x])?
                            },
                        }
                    },
                    | Inst::Cmp(c, x, y) => {
                        let (x, y) = (get(&values, x)?, get(&values, y)?);
                        let cc = match c {
                            | Cmp::Lt => FloatCC::LessThan,
                            | Cmp::Le => FloatCC::LessThanOrEqual,
                            | Cmp::Gt => FloatCC::GreaterThan,
                            | Cmp::Ge => FloatCC::GreaterThanOrEqual,
                            | Cmp::Eq => FloatCC::Equal,
                            | Cmp::Ne => FloatCC::NotEqual,
                        };
                        let flag = b.ins().fcmp(cc, x, y);
                        let (one, zero) = (b.ins().f64const(1.0), b.ins().f64const(0.0));
                        b.ins().select(flag, one, zero)
                    },
                    | Inst::Select(c, x, y) => {
                        let (c, x, y) = (get(&values, c)?, get(&values, x)?, get(&values, y)?);
                        let zero = b.ins().f64const(0.0);
                        let flag = b.ins().fcmp(FloatCC::NotEqual, c, zero);
                        b.ins().select(flag, x, y)
                    },
                    | Inst::Gather(base, index) => {
                        let index = get(&values, index)?;
                        let start = usize::try_from(base).unwrap_or(usize::MAX);
                        let span = tape.inputs.saturating_sub(start);
                        let nan = b.ins().f64const(f64::NAN);
                        if span == 0 {
                            // The tensor is not in the input buffer at all;
                            // reading anything would be a wild load.
                            nan
                        } else {
                            let zero = b.ins().f64const(0.0);
                            let span_f = b
                                .ins()
                                .f64const(f64::from(u32::try_from(span).map_err(error)?));
                            let last = b
                                .ins()
                                .f64const(f64::from(u32::try_from(span - 1).map_err(error)?));
                            let in_range = {
                                let lo = b.ins().fcmp(FloatCC::GreaterThanOrEqual, index, zero);
                                let hi = b.ins().fcmp(FloatCC::LessThan, index, span_f);
                                b.ins().band(lo, hi)
                            };
                            // Clamping before the conversion is what makes
                            // the load below provably in bounds; the `select`
                            // then replaces the value the clamp produced when
                            // the original index was out of range.
                            let floor = b.ins().fmin(index, last);
                            let clamped = b.ins().fmax(floor, zero);
                            let offset = b.ins().fcvt_to_uint_sat(types::I64, clamped);
                            let bytes = b.ins().imul_imm_s(offset, 8);
                            let address = b.ins().iadd_imm_s(
                                bytes,
                                i64::try_from(start).map_err(error)?.wrapping_mul(8),
                            );
                            let pointer = b.ins().iadd(in_ptr, address);
                            let loaded = b.ins().load(types::F64, flags, pointer, 0);
                            b.ins().select(in_range, loaded, nan)
                        }
                    },
                    | Inst::Call(f, start, len) => {
                        let range = start as usize..(start as usize).saturating_add(len as usize);
                        let regs = tape
                            .args
                            .get(range)
                            .ok_or_else(|| BackendError::Codegen("arguments".to_owned()))?;
                        let size =
                            u32::try_from(regs.len().saturating_mul(8).max(8)).map_err(error)?;
                        let slot = b.create_sized_stack_slot(StackSlotData::new(
                            StackSlotKind::ExplicitSlot,
                            size,
                            3,
                        ));
                        for (k, &r) in regs.iter().enumerate() {
                            let v = get(&values, r)?;
                            b.ins().stack_store(
                                types::I64,
                                v,
                                slot,
                                i32::try_from(k.saturating_mul(8)).map_err(error)?,
                            );
                        }
                        let address = b.ins().stack_addr(pointer, slot, 0);
                        let fptr = b
                            .ins()
                            .iconst(pointer, i64::try_from(f as usize).map_err(error)?);
                        let n = b
                            .ins()
                            .iconst(pointer, i64::try_from(regs.len()).map_err(error)?);
                        call1(&mut b, "rssn_eval", &[fptr, address, n])?
                    },
                    | Inst::Extern1(f, x) => {
                        let x = get(&values, x)?;
                        let fptr = b
                            .ins()
                            .iconst(pointer, i64::try_from(f as usize).map_err(error)?);
                        call1(&mut b, "rssn_extern1", &[fptr, x])?
                    },
                    | Inst::Extern2(f, x, y) => {
                        let (x, y) = (get(&values, x)?, get(&values, y)?);
                        let fptr = b
                            .ins()
                            .iconst(pointer, i64::try_from(f as usize).map_err(error)?);
                        call1(&mut b, "rssn_extern2", &[fptr, x, y])?
                    },
                };
                values.push(v);
            }
            for (k, &r) in tape.outputs.iter().enumerate() {
                let v = get(&values, r)?;
                b.ins().store(
                    flags,
                    v,
                    out_ptr,
                    i32::try_from(k.saturating_mul(8)).map_err(error)?,
                );
            }
            if let Some((header, more, exit)) = loop_blocks {
                let fuel = b.block_params(header)[0];
                let is_tail = {
                    let flag = get(&values, tape.outputs[0])?;
                    let zero = b.ins().f64const(0.0);
                    b.ins().fcmp(FloatCC::NotEqual, flag, zero)
                };
                b.ins().brif(is_tail, more, &[], exit, &[fuel.into()]);
                b.switch_to_block(more);
                b.seal_block(more);
                let has_fuel = b.ins().icmp_imm_s(IntCC::NotEqual, fuel, 0);
                let next_fuel = b.ins().iadd_imm_s(fuel, -1);
                b.ins().brif(has_fuel, header, &[next_fuel.into()], exit, &[fuel.into()]);
                b.seal_block(header);
                b.switch_to_block(exit);
                b.seal_block(exit);
                let remaining = b.block_params(exit)[0];
                let fuel_offset = f64_byte_offset(u32::try_from(tape.outputs.len()).map_err(error)?)?;
                b.ins().store(flags, remaining, out_ptr, fuel_offset);
            }
            b.ins().return_(&[]);
            b.finalize(frontend_config);
        }
        module.define_function(id, &mut ctx).map_err(error)?;
        module.clear_context(&mut ctx);
        module.finalize_definitions().map_err(error)?;
        let code = module.get_finalized_function(id);
        // SAFETY (in `abi`): `code` is the finalized `rssn_program` with the
        // entry signature, and the module is kept alive by `_owner`.
        #[allow(unsafe_code)]
        let entry = unsafe { abi::entry(code) };
        Ok(Self {
            _owner: Arc::new(abi::Owner(Some(module))),
            entry,
            inputs: tape.inputs,
            params: if tail_loop { 1 } else { tape.params },
            outputs: tape.outputs.len() + usize::from(tail_loop),
        })
    }

    fn run(
        &self,
        inputs: &[f64],
        params: &[f64],
        out: &mut [f64],
    ) {
        // Pad short argument lists with NaN, as the interpreter does.
        let pad = |values: &[f64], n: usize| -> Option<Vec<f64>> {
            (values.len() < n).then(|| {
                let mut v = values.to_vec();
                v.resize(n, f64::NAN);
                v
            })
        };
        let (padded_in, padded_par) = (
            pad(inputs, self.inputs.max(1)),
            pad(params, self.params.max(1)),
        );
        let inputs = padded_in.as_deref().unwrap_or(inputs);
        let params = padded_par.as_deref().unwrap_or(params);
        let mut scratch;
        let out: &mut [f64] = if out.len() < self.outputs {
            scratch = vec![f64::NAN; self.outputs];
            &mut scratch
        } else {
            out
        };
        // SAFETY (in `abi`): lengths were checked above and `_owner` keeps
        // the code alive for the duration of the call.
        #[allow(unsafe_code)]
        unsafe {
            abi::call(self.entry, inputs, params, out);
        }
    }
}

/// A pure scalar tail-recursive step compiled together with a native loop.
///
/// The tape takes `n` Float inputs and returns `[tail_flag, result, next_0,
/// ..., next_n]`. A nonzero flag requests a tail jump. The generated loop
/// executes at most `fuel` jumps, retaining the final step's outputs when it
/// completes or encounters a jump beyond the budget. Only pure instructions
/// are accepted, so evaluating all Select arms has no externally visible
/// effects. No existing Tape API or interpreter is changed.
pub struct TailLoop(JitFunction);

/// The final step and the remaining number of permitted jumps.
pub struct TailLoopResult {
    /// The final step outputs, including its completion flag.
    pub outputs: Vec<f64>,
    /// Unused tail-jump budget.
    pub remaining: u64,
}

impl TailLoop {
    /// Compile a pure, well-formed multi-output scalar step to native control flow.
    ///
    /// # Errors
    /// Returns a code-generation error for invalid or unsupported tapes or
    /// hosts that cannot generate machine code.
    pub fn compile(tape: &Tape) -> Result<Self, BackendError> {
        let invalid = || BackendError::Codegen("invalid tail-loop tape".to_string());
        let expected = tape.inputs.checked_add(2).ok_or_else(invalid)?;
        if tape.params != 0 || !tape.args.is_empty() || tape.outputs.len() != expected {
            return Err(invalid());
        }
        // Input and result frame offsets must fit Cranelift's signed offset.
        f64_byte_offset(u32::try_from(expected).map_err(error)?)?;
        for (index, inst) in tape.insts.iter().enumerate() {
            let valid = |reg: u32| (reg as usize) < index;
            match *inst {
                Inst::Const(_) => {},
                Inst::Input(i) if (i as usize) < tape.inputs => {
                    f64_byte_offset(i)?;
                },
                Inst::Neg(a) if valid(a) => {},
                Inst::Add(a, b) | Inst::Sub(a, b) | Inst::Mul(a, b) | Inst::Div(a, b)
                    | Inst::Cmp(_, a, b) if valid(a) && valid(b) => {},
                Inst::Select(a, b, c) if valid(a) && valid(b) && valid(c) => {},
                _ => return Err(invalid()),
            }
        }
        if tape.outputs.iter().any(|&output| output as usize >= tape.insts.len()) {
            return Err(invalid());
        }
        JitFunction::compile_impl(tape, true).map(Self)
    }

    /// Run until a result or until the next jump would exceed `fuel`.
    /// Returns `None` for a mismatched input count.
    #[must_use]
    pub fn run(&self, inputs: &[f64], fuel: u64) -> Option<TailLoopResult> {
        if inputs.len() != self.0.inputs { return None; }
        let mut outputs = vec![0.0; self.0.outputs];
        self.0.run(inputs, &[f64::from_bits(fuel)], &mut outputs);
        let remaining = outputs.pop()?.to_bits();
        Some(TailLoopResult { outputs, remaining })
    }
}

impl CompiledMulti for JitFunction {
    fn inputs(&self) -> usize {
        self.inputs
    }

    fn params(&self) -> usize {
        self.params
    }

    fn outputs(&self) -> usize {
        self.outputs
    }

    fn eval(
        &self,
        inputs: &[f64],
        params: &[f64],
        out: &mut [f64],
    ) {
        self.run(inputs, params, out);
    }
}

impl Compiled for JitFunction {
    fn arity(&self) -> usize {
        self.inputs
    }

    fn call(
        &self,
        args: &[f64],
    ) -> f64 {
        let mut out = [f64::NAN];
        self.run(args, &[], &mut out);
        out[0]
    }

    fn call_batch(
        &self,
        columns: &[&[f64]],
        out: &mut [f64],
    ) {
        let mut point = vec![f64::NAN; self.inputs.max(1)];
        let mut value = [f64::NAN];
        for (row, slot) in out.iter_mut().enumerate() {
            for (arg, column) in point.iter_mut().zip(columns) {
                *arg = column.get(row).copied().unwrap_or(f64::NAN);
            }
            self.run(&point, &[], &mut value);
            *slot = value[0];
        }
    }
}

/// `x^n` by square-and-multiply, the same sequence as
/// [`tape::int_pow`](super::tape::int_pow).
fn emit_powi(
    b: &mut FunctionBuilder<'_>,
    x: Value,
    n: i32,
) -> Value {
    let mut e = n.unsigned_abs();
    let mut base = x;
    let mut acc: Option<Value> = None;
    while e > 0 {
        if e & 1 == 1 {
            acc = Some(match acc {
                | Some(a) => b.ins().fmul(a, base),
                | None => base,
            });
        }
        e >>= 1;
        if e > 0 {
            base = b.ins().fmul(base, base);
        }
    }
    let value = acc.unwrap_or_else(|| b.ins().f64const(1.0));
    if n < 0 {
        let one = b.ins().f64const(1.0);
        b.ins().fdiv(one, value)
    } else {
        value
    }
}

/// A shared compiled function as a [`Compiled`] trait object.
struct Handle(Arc<JitFunction>);

impl Compiled for Handle {
    fn arity(&self) -> usize {
        self.0.arity()
    }

    fn call(
        &self,
        args: &[f64],
    ) -> f64 {
        self.0.call(args)
    }

    fn call_batch(
        &self,
        columns: &[&[f64]],
        out: &mut [f64],
    ) {
        self.0.call_batch(columns, out);
    }
}

impl CompiledMulti for Handle {
    fn inputs(&self) -> usize {
        self.0.inputs
    }

    fn params(&self) -> usize {
        self.0.params
    }

    fn outputs(&self) -> usize {
        self.0.outputs
    }

    fn eval(
        &self,
        inputs: &[f64],
        params: &[f64],
        out: &mut [f64],
    ) {
        self.0.run(inputs, params, out);
    }
}

/// The Cranelift backend, with a cache keyed by the structural hash of the
/// optimised tape, mixed with a context key.
///
/// # Why the key has two halves
///
/// Two graphs that lower to the same tape do not necessarily mean the same
/// thing: an operator's *semantics* come from the rule set that installed it,
/// and the same tape compiled for a different CPU may not be valid. The
/// structural hash covers the first half — the program — and
/// [`CraneliftBackend::key_extra`] covers the second: the manifest hash of the
/// build plus the target device. A host that keeps a compiled cache on disk
/// builds the same pair and gets the same key.
#[derive(Default)]
pub struct CraneliftBackend {
    cache: Mutex<HashMap<u64, Arc<JitFunction>>>,
    key_extra: u64,
}

impl CraneliftBackend {
    /// A backend with an empty cache and no context key.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// A backend whose cache keys are mixed with `key_extra`.
    ///
    /// Pass [`target_key`] combined with the manifest hash:
    /// `backend = CraneliftBackend::with_key(manifest.hash_u64() ^ target_key())`.
    #[must_use]
    pub fn with_key(key_extra: u64) -> Self {
        Self {
            cache: Mutex::new(HashMap::new()),
            key_extra,
        }
    }

    /// The context key mixed into every cache key.
    #[must_use]
    pub const fn key_extra(&self) -> u64 {
        self.key_extra
    }

    /// Changes the context key, dropping everything compiled under the old
    /// one.
    pub fn set_key_extra(
        &mut self,
        key_extra: u64,
    ) {
        self.key_extra = key_extra;
        if let Ok(mut cache) = self.cache.lock() {
            cache.clear();
        }
    }

    /// The cache key of `tape` under this backend's context.
    #[must_use]
    pub fn cache_key(
        &self,
        tape: &Tape,
    ) -> u64 {
        mix(tape.structural_hash(), self.key_extra)
    }

    fn get_or_compile(
        &self,
        tape: &Tape,
    ) -> Result<Arc<JitFunction>, BackendError> {
        let key = self.cache_key(tape);
        if let Some(f) = self.cache.lock().ok().and_then(|c| c.get(&key).cloned()) {
            return Ok(f);
        }
        let compiled = Arc::new(JitFunction::compile(tape)?);
        if let Ok(mut cache) = self.cache.lock() {
            cache.insert(key, Arc::clone(&compiled));
        }
        Ok(compiled)
    }

    /// Number of distinct programs compiled so far.
    #[must_use]
    pub fn cached(&self) -> usize {
        self.cache.lock().map_or(0, |c| c.len())
    }
}

/// Mixes two hashes so that neither half can cancel the other.
///
/// A plain `^` would make `(h, k)` and `(h ^ k, 0)` collide; the multiply
/// and rotate below are the finaliser of splitmix64, which is cheap and has
/// no fixed points worth exploiting.
const fn mix(
    left: u64,
    right: u64,
) -> u64 {
    let mut z = left
        .wrapping_mul(0x9E37_79B9_7F4A_7C15)
        .wrapping_add(right.rotate_left(32));
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// A stable identifier for the machine this code will run on.
///
/// The value is a hash of the target triple and the CPU feature set the
/// compiler was configured with, so two builds for different micro-
/// architectures do not share cached machine code.
#[must_use]
pub fn target_key() -> u64 {
    use std::hash::Hash;
    use std::hash::Hasher;
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    std::env::consts::ARCH.hash(&mut hasher);
    std::env::consts::OS.hash(&mut hasher);
    // `cfg!(target_feature = ...)` is the feature set this build was compiled
    // for, which is what the generated code may assume.
    for feature in ["sse2", "sse4.1", "avx", "avx2", "fma", "neon", "sve"] {
        let present = match feature {
            | "sse2" => cfg!(target_feature = "sse2"),
            | "sse4.1" => cfg!(target_feature = "sse4.1"),
            | "avx" => cfg!(target_feature = "avx"),
            | "avx2" => cfg!(target_feature = "avx2"),
            | "fma" => cfg!(target_feature = "fma"),
            | "neon" => cfg!(target_feature = "neon"),
            | "sve" => cfg!(target_feature = "sve"),
            | _ => false,
        };
        feature.hash(&mut hasher);
        present.hash(&mut hasher);
    }
    hasher.finish()
}

impl Backend for CraneliftBackend {
    fn name(&self) -> &'static str {
        "cranelift"
    }

    fn compile(
        &self,
        graph: &Graph,
        root: NodeId,
        inputs: &[SymbolId],
    ) -> Result<Box<dyn Compiled>, BackendError> {
        let tape = super::tape::lower(graph, &[root], inputs, &[])?.optimise();
        Ok(Box::new(Handle(self.get_or_compile(&tape)?)))
    }

    fn compile_multi(
        &self,
        graph: &Graph,
        roots: &[NodeId],
        inputs: &[SymbolId],
        params: &[SymbolId],
    ) -> Result<Box<dyn CompiledMulti>, BackendError> {
        let tape = super::tape::lower(graph, roots, inputs, params)?.optimise();
        Ok(Box::new(Handle(self.get_or_compile(&tape)?)))
    }
}

/// Tiered execution: interpret the tape at first, compile to machine code
/// on a background thread after `threshold` calls, then switch.
#[derive(Clone, Copy, Debug)]
pub struct TieredBackend {
    /// Calls (or batch rows) before compilation starts.
    pub threshold: usize,
}

impl Default for TieredBackend {
    fn default() -> Self {
        Self { threshold: 2000 }
    }
}

struct Tiered {
    tape: Arc<Tape>,
    calls: AtomicUsize,
    started: AtomicBool,
    threshold: usize,
    jit: Arc<OnceLock<JitFunction>>,
}

impl Tiered {
    fn count(
        &self,
        n: usize,
    ) {
        let before = self.calls.fetch_add(n, Ordering::Relaxed);
        if before.saturating_add(n) >= self.threshold && !self.started.swap(true, Ordering::AcqRel)
        {
            let (tape, slot) = (Arc::clone(&self.tape), Arc::clone(&self.jit));
            std::thread::spawn(move || {
                if let Ok(f) = JitFunction::compile(&tape) {
                    let _ = slot.set(f);
                }
            });
        }
    }

    /// Whether the compiled tier is active.
    fn compiled(&self) -> bool {
        self.jit.get().is_some()
    }
}

impl Compiled for Tiered {
    fn arity(&self) -> usize {
        self.tape.inputs
    }

    fn call(
        &self,
        args: &[f64],
    ) -> f64 {
        if let Some(f) = self.jit.get() {
            return f.call(args);
        }
        self.count(1);
        let mut out = [f64::NAN];
        self.tape.eval(args, &[], &mut out);
        out[0]
    }

    fn call_batch(
        &self,
        columns: &[&[f64]],
        out: &mut [f64],
    ) {
        if let Some(f) = self.jit.get() {
            return f.call_batch(columns, out);
        }
        self.count(out.len());
        let mut point = vec![f64::NAN; self.tape.inputs];
        let mut values = Vec::with_capacity(self.tape.insts.len());
        let mut value = [f64::NAN];
        for (row, slot) in out.iter_mut().enumerate() {
            for (arg, column) in point.iter_mut().zip(columns) {
                *arg = column.get(row).copied().unwrap_or(f64::NAN);
            }
            self.tape.run(&point, &[], &mut values, &mut value);
            *slot = value[0];
        }
    }
}

impl Backend for TieredBackend {
    fn name(&self) -> &'static str {
        "tiered"
    }

    fn compile(
        &self,
        graph: &Graph,
        root: NodeId,
        inputs: &[SymbolId],
    ) -> Result<Box<dyn Compiled>, BackendError> {
        let tape = super::tape::lower(graph, &[root], inputs, &[])?.optimise();
        Ok(Box::new(Tiered {
            tape: Arc::new(tape),
            calls: AtomicUsize::new(0),
            started: AtomicBool::new(false),
            threshold: self.threshold,
            jit: Arc::new(OnceLock::new()),
        }))
    }
}

/// Whether a tiered function returned by [`TieredBackend`] has switched to
/// machine code (for diagnostics and tests).
#[must_use]
pub fn is_compiled(f: &dyn std::any::Any) -> bool {
    f.downcast_ref::<Tiered>().is_some_and(Tiered::compiled)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::backend::Interpreter;
    use crate::graph::Engine;
    use crate::rules;

    fn setup(src: &str) -> (Graph, NodeId, Vec<SymbolId>) {
        let mut g = Graph::new();
        assert!(Engine::install(&mut g, &rules::standard()).is_ok());
        let root = g.parse(src).unwrap_or_else(|e| panic!("{e}"));
        let inputs = vec![g.interner_mut().symbol("x"), g.interner_mut().symbol("y")];
        (g, root, inputs)
    }

    fn ulps(
        a: f64,
        b: f64,
    ) -> u64 {
        if a.partial_cmp(&b) == Some(std::cmp::Ordering::Equal) || (a.is_nan() && b.is_nan()) {
            return 0;
        }
        let (ia, ib) = (a.to_bits().cast_signed(), b.to_bits().cast_signed());
        ia.abs_diff(ib)
    }

    #[test]
    fn large_input_and_parameter_offsets_are_checked() {
        // Constructing this tape does not allocate an input or parameter buffer.
        let limit = (i32::MAX / 8) as u32;
        for (input, index, accepted) in [
            (true, limit, true),
            (true, limit + 1, false),
            (false, limit, true),
            (false, limit + 1, false),
        ] {
            let tape = Tape {
                insts: vec![if input {
                    Inst::Input(index)
                } else {
                    Inst::Param(index)
                }],
                args: Vec::new(),
                inputs: if input {
                    index as usize + 1
                } else {
                    0
                },
                params: if input {
                    0
                } else {
                    index as usize + 1
                },
                outputs: vec![0],
            };
            assert_eq!(
                JitFunction::compile(&tape).is_ok(),
                accepted,
                "input={input}, index={index}"
            );
        }
    }

    #[test]
    fn agrees_with_the_interpreter() {
        let backend = CraneliftBackend::new();
        let points = [
            (0.3, -1.2),
            (1.7, 0.4),
            (-2.5, 3.1),
            (0.0, 0.0),
            (1e-3, 7.0),
        ];
        for src in [
            "x^2 + 3*x*y - y/7",
            "sin(x)^2*exp(-y) + (x - y)/(1 + x^2) + pi",
            "sqrt(x^2 + y^2) + x^(3/2)*y^(-1/2) + (x^2 + 1)^(1/3)",
            "atan(x) + atan2(y, x) + tanh(x*y) + ln(1 + x^2) + abs(x - y)",
            "x^7 - 2*x^5*y^2 + y^(-3)",
            "gamma(x^2 + 1) + erf(y)",
            "lt(x, y) + ge(x, 0) + heaviside(x - y)",
            "-x - y - 2*x*y",
            "1/(x - 1)^2",
        ] {
            let (g, root, inputs) = setup(src);
            let reference = Interpreter
                .compile(&g, root, &inputs)
                .unwrap_or_else(|e| panic!("{e}"));
            let jit = backend
                .compile(&g, root, &inputs)
                .unwrap_or_else(|e| panic!("{src}: {e}"));
            for &(x, y) in &points {
                let (want, got) = (reference.call(&[x, y]), jit.call(&[x, y]));
                let close =
                    ulps(want, got) <= 8 || (want - got).abs() <= 1e-14 * want.abs().max(1.0);
                assert!(close, "{src} at ({x}, {y}): interpreter {want}, jit {got}");
            }
        }
    }

    #[test]
    fn tape_and_jit_agree_bitwise() {
        // The JIT executes exactly the tape's operations.
        let backend = CraneliftBackend::new();
        let (g, root, inputs) = setup("x^9*y - 3*x/(y^2 + 1) + sqrt(x*x + 1) - exp(sin(y))");
        let tape = crate::backend::tape::lower(&g, &[root], &inputs, &[])
            .unwrap_or_else(|e| panic!("{e}"))
            .optimise();
        let jit = backend
            .compile(&g, root, &inputs)
            .unwrap_or_else(|e| panic!("{e}"));
        for (x, y) in [(0.1, 0.2), (-3.7, 1.9), (12.0, -0.003)] {
            let mut out = [0.0];
            tape.eval(&[x, y], &[], &mut out);
            assert_eq!(out[0].to_bits(), jit.call(&[x, y]).to_bits());
        }
    }

    #[test]
    fn multiple_outputs_parameters_and_the_cache() {
        let backend = CraneliftBackend::new();
        let mut g = Graph::new();
        assert!(Engine::install(&mut g, &rules::standard()).is_ok());
        let roots: Vec<NodeId> = ["a*x + y", "x*y*a", "sin(a*x + y)"]
            .iter()
            .map(|s| g.parse(s).unwrap_or_else(|e| panic!("{e}")))
            .collect();
        let (x, y, a) = (
            g.interner_mut().symbol("x"),
            g.interner_mut().symbol("y"),
            g.interner_mut().symbol("a"),
        );
        let f = backend
            .compile_multi(&g, &roots, &[x, y], &[a])
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!((f.inputs(), f.params(), f.outputs()), (2, 1, 3));
        let mut out = [0.0; 3];
        f.eval(&[2.0, 3.0], &[0.5], &mut out);
        assert!(
            (out[0] - 4.0).abs() < 1e-15
                && (out[1] - 3.0).abs() < 1e-15
                && (out[2] - 4.0_f64.sin()).abs() < 1e-15
        );
        // Batch evaluation with a hoisted parameter.
        let xs = [1.0, 2.0];
        let ys = [0.0, 1.0];
        let (mut o1, mut o2, mut o3) = ([0.0; 2], [0.0; 2], [0.0; 2]);
        f.eval_batch(&[&xs, &ys], &[2.0], &mut [&mut o1, &mut o2, &mut o3]);
        assert!((o1[1] - 5.0).abs() < 1e-15 && (o2[1] - 4.0).abs() < 1e-15);
        // The same program in a fresh graph hits the cache.
        let before = backend.cached();
        let (g2, root2, inputs2) = setup("x^2 + 3*x*y - y/7");
        let _ = backend.compile(&g2, root2, &inputs2);
        let _ = backend.compile(&g2, root2, &inputs2);
        assert_eq!(backend.cached(), before + 1);
    }

    #[test]
    fn tiered_execution_switches_to_machine_code() {
        let (g, root, inputs) = setup("x^3 - y");
        let tiered = TieredBackend { threshold: 10 };
        let f = tiered
            .compile(&g, root, &inputs)
            .unwrap_or_else(|e| panic!("{e}"));
        for _ in 0..50 {
            assert!((f.call(&[2.0, 1.0]) - 7.0).abs() < 1e-15);
        }
        // Give the background compilation time to finish; results are the
        // same either way.
        std::thread::sleep(std::time::Duration::from_millis(200));
        assert!((f.call(&[3.0, 0.0]) - 27.0).abs() < 1e-15);
    }

    #[test]
    fn kernels_run_on_the_jit() {
        // The numeric quadrature fallback and the rest of the engine work
        // with the JIT as the current backend.
        let rules = [crate::rules::calculus()];
        let (value, error) =
            crate::backend::with_backend(Arc::new(CraneliftBackend::new()), || {
                crate::rules::testing::numeric(
                    &rules,
                    "defint(exp(-x^2)*cos(x^3), x, 0, 2)",
                    &[],
                    1e-12,
                )
            });
        let (reference, _) = crate::rules::testing::numeric(
            &rules,
            "defint(exp(-x^2)*cos(x^3), x, 0, 2)",
            &[],
            1e-12,
        );
        assert!(
            (value - reference).abs() < 1e-12 && error < 1e-10,
            "{value} vs {reference} ± {error}"
        );
    }
}
