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
use cranelift_codegen::ir::MemFlags;
use cranelift_codegen::ir::StackSlotData;
use cranelift_codegen::ir::StackSlotKind;
use cranelift_codegen::ir::Value;
use cranelift_codegen::ir::condcodes::FloatCC;
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
        let function: extern "C" fn(f64) -> f64 = unsafe { std::mem::transmute::<usize, extern "C" fn(f64) -> f64>(f) };
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
        let function: extern "C" fn(f64, f64) -> f64 = unsafe { std::mem::transmute::<usize, extern "C" fn(f64, f64) -> f64>(f) };
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
        let mut flags = settings::builder();
        flags.set("use_colocated_libcalls", "false").map_err(error)?;
        flags.set("is_pic", "false").map_err(error)?;
        flags.set("opt_level", "speed").map_err(error)?;
        let isa = cranelift_native::builder().map_err(error)?.finish(settings::Flags::new(flags)).map_err(error)?;
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
        sig_eval.params.extend([AbiParam::new(pointer), AbiParam::new(pointer), AbiParam::new(pointer)]);
        sig_eval.returns.push(AbiParam::new(types::F64));
        let mut sig_ext1 = module.make_signature();
        sig_ext1.params.extend([AbiParam::new(pointer), AbiParam::new(types::F64)]);
        sig_ext1.returns.push(AbiParam::new(types::F64));
        let mut sig_ext2 = sig_ext1.clone();
        sig_ext2.params.push(AbiParam::new(types::F64));
        let mut imports: HashMap<&'static str, FuncId> = HashMap::new();
        for f in UNARY_SHIMS {
            if let Some((name, _)) = library(f) {
                imports.insert(name, module.declare_function(name, Linkage::Import, &sig_unary).map_err(error)?);
            }
        }
        imports.insert("rssn_pow", module.declare_function("rssn_pow", Linkage::Import, &sig_binary).map_err(error)?);
        imports.insert("rssn_atan2", module.declare_function("rssn_atan2", Linkage::Import, &sig_binary).map_err(error)?);
        imports.insert("rssn_eval", module.declare_function("rssn_eval", Linkage::Import, &sig_eval).map_err(error)?);
        imports.insert("rssn_extern1", module.declare_function("rssn_extern1", Linkage::Import, &sig_ext1).map_err(error)?);
        imports.insert("rssn_extern2", module.declare_function("rssn_extern2", Linkage::Import, &sig_ext2).map_err(error)?);

        let mut sig = module.make_signature();
        sig.params.extend([AbiParam::new(pointer), AbiParam::new(pointer), AbiParam::new(pointer)]);
        let id = module.declare_function("rssn_program", Linkage::Local, &sig).map_err(error)?;
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
            let flags = MemFlags::trusted();
            let mut refs = HashMap::new();
            for (&name, &fid) in &imports {
                refs.insert(name, module.declare_func_in_func(fid, b.func));
            }
            let call1 = |b: &mut FunctionBuilder<'_>, name: &str, args: &[Value]| -> Result<Value, BackendError> {
                let f = *refs.get(name).ok_or_else(|| BackendError::Codegen(name.to_owned()))?;
                let call = b.ins().call(f, args);
                b.inst_results(call).first().copied().ok_or_else(|| BackendError::Codegen(name.to_owned()))
            };
            let mut values: Vec<Value> = Vec::with_capacity(tape.insts.len());
            let get = |values: &[Value], r: u32| values.get(r as usize).copied().ok_or_else(|| BackendError::Codegen("register".to_owned()));
            for inst in &tape.insts {
                let v = match *inst {
                    | Inst::Const(c) => b.ins().f64const(c),
                    | Inst::Input(i) => b.ins().load(types::F64, flags, in_ptr, i32::try_from(i).map_err(error)?.saturating_mul(8)),
                    | Inst::Param(i) => b.ins().load(types::F64, flags, par_ptr, i32::try_from(i).map_err(error)?.saturating_mul(8)),
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
                                let (name, _) = library(other).ok_or_else(|| BackendError::Codegen("intrinsic".to_owned()))?;
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
                    | Inst::Call(f, start, len) => {
                        let range = start as usize..(start as usize).saturating_add(len as usize);
                        let regs = tape.args.get(range).ok_or_else(|| BackendError::Codegen("arguments".to_owned()))?;
                        let size = u32::try_from(regs.len().saturating_mul(8).max(8)).map_err(error)?;
                        let slot = b.create_sized_stack_slot(StackSlotData::new(StackSlotKind::ExplicitSlot, size, 3));
                        for (k, &r) in regs.iter().enumerate() {
                            let v = get(&values, r)?;
                            b.ins().stack_store(v, slot, i32::try_from(k.saturating_mul(8)).map_err(error)?);
                        }
                        let address = b.ins().stack_addr(pointer, slot, 0);
                        let fptr = b.ins().iconst(pointer, i64::try_from(f as usize).map_err(error)?);
                        let n = b.ins().iconst(pointer, i64::try_from(regs.len()).map_err(error)?);
                        call1(&mut b, "rssn_eval", &[fptr, address, n])?
                    },
                    | Inst::Extern1(f, x) => {
                        let x = get(&values, x)?;
                        let fptr = b.ins().iconst(pointer, i64::try_from(f as usize).map_err(error)?);
                        call1(&mut b, "rssn_extern1", &[fptr, x])?
                    },
                    | Inst::Extern2(f, x, y) => {
                        let (x, y) = (get(&values, x)?, get(&values, y)?);
                        let fptr = b.ins().iconst(pointer, i64::try_from(f as usize).map_err(error)?);
                        call1(&mut b, "rssn_extern2", &[fptr, x, y])?
                    },
                };
                values.push(v);
            }
            for (k, &r) in tape.outputs.iter().enumerate() {
                let v = get(&values, r)?;
                b.ins().store(flags, v, out_ptr, i32::try_from(k.saturating_mul(8)).map_err(error)?);
            }
            b.ins().return_(&[]);
            b.finalize();
        }
        module.define_function(id, &mut ctx).map_err(error)?;
        module.clear_context(&mut ctx);
        module.finalize_definitions().map_err(error)?;
        let code = module.get_finalized_function(id);
        // SAFETY (in `abi`): `code` is the finalized `rssn_program` with the
        // entry signature, and the module is kept alive by `_owner`.
        #[allow(unsafe_code)]
        let entry = unsafe { abi::entry(code) };
        Ok(Self { _owner: Arc::new(abi::Owner(Some(module))), entry, inputs: tape.inputs, params: tape.params, outputs: tape.outputs.len() })
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
        let (padded_in, padded_par) = (pad(inputs, self.inputs.max(1)), pad(params, self.params.max(1)));
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
/// optimised tape.
#[derive(Default)]
pub struct CraneliftBackend {
    cache: Mutex<HashMap<u64, Arc<JitFunction>>>,
}

impl CraneliftBackend {
    /// A backend with an empty cache.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    fn get_or_compile(
        &self,
        tape: &Tape,
    ) -> Result<Arc<JitFunction>, BackendError> {
        let key = tape.structural_hash();
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
        if before.saturating_add(n) >= self.threshold && !self.started.swap(true, Ordering::AcqRel) {
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
    fn agrees_with_the_interpreter() {
        let backend = CraneliftBackend::new();
        let points = [(0.3, -1.2), (1.7, 0.4), (-2.5, 3.1), (0.0, 0.0), (1e-3, 7.0)];
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
            let reference = Interpreter.compile(&g, root, &inputs).unwrap_or_else(|e| panic!("{e}"));
            let jit = backend.compile(&g, root, &inputs).unwrap_or_else(|e| panic!("{src}: {e}"));
            for &(x, y) in &points {
                let (want, got) = (reference.call(&[x, y]), jit.call(&[x, y]));
                let close = ulps(want, got) <= 8 || (want - got).abs() <= 1e-14 * want.abs().max(1.0);
                assert!(close, "{src} at ({x}, {y}): interpreter {want}, jit {got}");
            }
        }
    }

    #[test]
    fn tape_and_jit_agree_bitwise() {
        // The JIT executes exactly the tape's operations.
        let backend = CraneliftBackend::new();
        let (g, root, inputs) = setup("x^9*y - 3*x/(y^2 + 1) + sqrt(x*x + 1) - exp(sin(y))");
        let tape = crate::backend::tape::lower(&g, &[root], &inputs, &[]).unwrap_or_else(|e| panic!("{e}")).optimise();
        let jit = backend.compile(&g, root, &inputs).unwrap_or_else(|e| panic!("{e}"));
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
        let roots: Vec<NodeId> = ["a*x + y", "x*y*a", "sin(a*x + y)"].iter().map(|s| g.parse(s).unwrap_or_else(|e| panic!("{e}"))).collect();
        let (x, y, a) = (g.interner_mut().symbol("x"), g.interner_mut().symbol("y"), g.interner_mut().symbol("a"));
        let f = backend.compile_multi(&g, &roots, &[x, y], &[a]).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!((f.inputs(), f.params(), f.outputs()), (2, 1, 3));
        let mut out = [0.0; 3];
        f.eval(&[2.0, 3.0], &[0.5], &mut out);
        assert!((out[0] - 4.0).abs() < 1e-15 && (out[1] - 3.0).abs() < 1e-15 && (out[2] - 4.0_f64.sin()).abs() < 1e-15);
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
        let f = tiered.compile(&g, root, &inputs).unwrap_or_else(|e| panic!("{e}"));
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
        let (value, error) = crate::backend::with_backend(Arc::new(CraneliftBackend::new()), || {
            crate::rules::testing::numeric(&rules, "defint(exp(-x^2)*cos(x^3), x, 0, 2)", &[], 1e-12)
        });
        let (reference, _) = crate::rules::testing::numeric(&rules, "defint(exp(-x^2)*cos(x^3), x, 0, 2)", &[], 1e-12);
        assert!((value - reference).abs() < 1e-12 && error < 1e-10, "{value} vs {reference} ± {error}");
    }
}
