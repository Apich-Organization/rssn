//! # LIR to machine code (Cranelift)
//!
//! Compiles a whole [`lir::Module`](super::lir::Module) into one Cranelift
//! module (feature `jit`). The translation is direct because LIR's registers
//! are Cranelift *variables*: every register becomes one variable (a slice
//! becomes two, address and length), instructions `def_var` their
//! destination, and Cranelift's SSA builder inserts the block parameters.
//!
//! ## Calling convention
//!
//! Every LIR function `f(p₁ … pₙ) -> (r₁ … rₘ)` becomes a native function
//!
//! ```text
//! f(ctx: *mut u32, depth: i64, p₁ … pₙ [, ret: *mut u64]) [-> r₁]
//! ```
//!
//! where a slice parameter expands to `(address, length)`, a boolean is an
//! `i8`, `ret` is present when `m > 1` and the single result is returned
//! directly when `m == 1`. `ctx` holds the bail status: [`Term::Bail`] and
//! every failed check (overflow, conversion, depth) store their code there
//! and return; every call site reloads it and propagates a bail by
//! returning. `depth` counts the call nesting so the depth limit agrees with
//! the interpreter.
//!
//! The host enters through a per-function trampoline
//! `fn(args: *const u64, out: *mut u64) -> u32` that unpacks arguments from
//! 8-byte slots, owns the status word, and returns it. All `unsafe` stays in
//! [`super::jit`]'s `abi` module.
//!
//! ## Agreement with the interpreter
//!
//! The executors agree bit for bit on the same module: float operations are
//! the same IEEE operations in the same order, `FPowI` emits the same
//! square-and-multiply chain, transcendental functions call the same Rust
//! functions through the same shims, `FMin`/`FMax` are IEEE minimum/maximum
//! in both, and every check (overflow, bounds, conversion, depth) has the
//! same outcome. `tests::agrees_with_the_interpreter` checks it.

use std::collections::HashMap;
use std::sync::Arc;

use cranelift_codegen::ir::AbiParam;
use cranelift_codegen::ir::Block as ClBlock;
use cranelift_codegen::ir::InstBuilder;
use cranelift_codegen::ir::MemFlagsData;
use cranelift_codegen::ir::Signature;
use cranelift_codegen::ir::StackSlotData;
use cranelift_codegen::ir::StackSlotKind;
use cranelift_codegen::ir::Type;
use cranelift_codegen::ir::Value as ClValue;
use cranelift_codegen::ir::condcodes::FloatCC;
use cranelift_codegen::ir::condcodes::IntCC;
use cranelift_codegen::ir::types;
use cranelift_codegen::settings;
use cranelift_codegen::settings::Configurable;
use cranelift_frontend::FunctionBuilder as ClBuilder;
use cranelift_frontend::FunctionBuilderContext;
use cranelift_frontend::Variable;
use cranelift_jit::JITBuilder;
use cranelift_jit::JITModule;
use cranelift_module::FuncId as ClFuncId;
use cranelift_module::Linkage;
use cranelift_module::Module as _;

use super::BackendError;
use super::Intrinsic;
use super::jit::abi;
use super::lir::BAIL_CONVERT;
use super::lir::BAIL_DEPTH;
use super::lir::BAIL_OVERFLOW;
use super::lir::FuncId;
use super::lir::Inst;
use super::lir::LirError;
use super::lir::MAX_CALL_DEPTH;
use super::lir::Module;
use super::lir::Outcome;
use super::lir::Reg;
use super::lir::Term;
use super::lir::Ty;
use super::lir::Value;
use super::tape::Cmp;

fn error(e: impl std::fmt::Display) -> BackendError {
    BackendError::Codegen(e.to_string())
}

fn cl_type(ty: Ty) -> Type {
    match ty {
        | Ty::F64 => types::F64,
        | Ty::I64 | Ty::Slice => types::I64,
        | Ty::Bool => types::I8,
    }
}

const fn float_cc(c: Cmp) -> FloatCC {
    match c {
        | Cmp::Lt => FloatCC::LessThan,
        | Cmp::Le => FloatCC::LessThanOrEqual,
        | Cmp::Gt => FloatCC::GreaterThan,
        | Cmp::Ge => FloatCC::GreaterThanOrEqual,
        | Cmp::Eq => FloatCC::Equal,
        | Cmp::Ne => FloatCC::NotEqual,
    }
}

const fn int_cc(c: Cmp) -> IntCC {
    match c {
        | Cmp::Lt => IntCC::SignedLessThan,
        | Cmp::Le => IntCC::SignedLessThanOrEqual,
        | Cmp::Gt => IntCC::SignedGreaterThan,
        | Cmp::Ge => IntCC::SignedGreaterThanOrEqual,
        | Cmp::Eq => IntCC::Equal,
        | Cmp::Ne => IntCC::NotEqual,
    }
}

/// The native parameter types of a function (without `ctx`/`depth`/`ret`).
fn native_params(f: &super::lir::Function) -> Vec<Type> {
    let mut out = Vec::new();
    for &p in f.blocks.first().map_or(&[][..], |e| &e.params[..]) {
        let ty = f.regs.get(p.0 as usize).copied().unwrap_or(Ty::F64);
        out.push(cl_type(ty));
        if ty == Ty::Slice {
            out.push(types::I64);
        }
    }
    out
}

/// A compiled LIR module.
pub struct LirModule {
    _owner: Arc<abi::Owner>,
    entries: Vec<abi::LirEntry>,
    module: Module,
}

impl std::fmt::Debug for LirModule {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        f.debug_struct("LirModule")
            .field("functions", &self.module.functions.len())
            .finish_non_exhaustive()
    }
}

impl LirModule {
    /// Verifies and compiles every function of `module`.
    ///
    /// # Errors
    /// The module does not verify, or Cranelift rejects it.
    pub fn compile(module: &Module) -> Result<Self, BackendError> {
        module.verify()?;
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
        let frontend_config = isa.frontend_config();
        let mut builder = JITBuilder::with_isa(isa, cranelift_module::default_libcall_names());
        super::jit::register_shims(&mut builder);
        let mut jit = JITModule::new(builder);
        let pointer = jit.target_config().pointer_type();
        let imports = super::jit::declare_shims(&mut jit)?;

        // Declare every function first so that calls (and recursion) resolve.
        let mut natives: Vec<(ClFuncId, Signature)> = Vec::new();
        for (index, f) in module.functions.iter().enumerate() {
            let mut sig = jit.make_signature();
            sig.params.push(AbiParam::new(pointer));
            sig.params.push(AbiParam::new(types::I64));
            for t in native_params(f) {
                sig.params.push(AbiParam::new(t));
            }
            match f.results.as_slice() {
                | [] => {},
                | [one] => sig.returns.push(AbiParam::new(cl_type(*one))),
                | _ => sig.params.push(AbiParam::new(pointer)),
            }
            let id = jit
                .declare_function(&format!("lir_{index}"), Linkage::Local, &sig)
                .map_err(error)?;
            natives.push((id, sig));
        }
        let mut trampolines = Vec::new();
        for index in 0..module.functions.len() {
            let mut sig = jit.make_signature();
            sig.params.push(AbiParam::new(pointer));
            sig.params.push(AbiParam::new(pointer));
            sig.returns.push(AbiParam::new(types::I32));
            let id = jit
                .declare_function(&format!("lir_entry_{index}"), Linkage::Local, &sig)
                .map_err(error)?;
            trampolines.push((id, sig));
        }

        let mut fctx = FunctionBuilderContext::new();
        for (index, f) in module.functions.iter().enumerate() {
            let (id, sig) = natives
                .get(index)
                .cloned()
                .ok_or_else(|| error("function table"))?;
            let mut ctx = jit.make_context();
            ctx.func.signature = sig;
            {
                let mut b = ClBuilder::new(&mut ctx.func, &mut fctx);
                let mut refs = HashMap::new();
                for (&name, &fid) in &imports {
                    refs.insert(name, jit.declare_func_in_func(fid, b.func));
                }
                let mut callees = Vec::new();
                for (fid, _) in &natives {
                    callees.push(jit.declare_func_in_func(*fid, b.func));
                }
                let mut cg = Gen {
                    b: &mut b,
                    f,
                    refs: &refs,
                    callees: &callees,
                    pointer,
                    vars: Vec::new(),
                    blocks: Vec::new(),
                    ctx_ptr: None,
                    depth: None,
                    ret_ptr: None,
                    bail: None,
                };
                cg.function()?;
                b.finalize(frontend_config);
            }
            jit.define_function(id, &mut ctx).map_err(error)?;
            jit.clear_context(&mut ctx);
        }
        for (index, f) in module.functions.iter().enumerate() {
            let (id, sig) = trampolines
                .get(index)
                .cloned()
                .ok_or_else(|| error("trampoline table"))?;
            let (native, _) = natives
                .get(index)
                .cloned()
                .ok_or_else(|| error("function table"))?;
            let mut ctx = jit.make_context();
            ctx.func.signature = sig;
            {
                let mut b = ClBuilder::new(&mut ctx.func, &mut fctx);
                let callee = jit.declare_func_in_func(native, b.func);
                trampoline(&mut b, f, callee, pointer);
                b.finalize(frontend_config);
            }
            jit.define_function(id, &mut ctx).map_err(error)?;
            jit.clear_context(&mut ctx);
        }
        jit.finalize_definitions().map_err(error)?;
        let mut entries = Vec::new();
        for (id, _) in &trampolines {
            let ptr = jit.get_finalized_function(*id);
            // SAFETY-relevant step lives in `abi`: the pointer is a trampoline
            // with the `LirEntry` signature, finalized above.
            entries.push(abi::lir_entry(ptr));
        }
        Ok(Self {
            _owner: Arc::new(abi::Owner(Some(jit))),
            entries,
            module: module.clone(),
        })
    }

    /// The module this was compiled from.
    #[must_use]
    pub const fn module(&self) -> &Module {
        &self.module
    }

    /// Calls `func` with `args` in machine code.
    ///
    /// # Errors
    /// The arguments do not match the function (see
    /// [`Module::check_args`]).
    pub fn call(
        &self,
        func: FuncId,
        args: &mut [Value<'_>],
    ) -> Result<Outcome, LirError> {
        let f = self.module.check_args(func, args)?;
        let entry = *self
            .entries
            .get(func.0 as usize)
            .ok_or_else(|| LirError(format!("unknown function f{}", func.0)))?;
        let mut slots: Vec<u64> = Vec::with_capacity(args.len().saturating_mul(2));
        for a in args.iter_mut() {
            match a {
                | Value::F(x) => slots.push(x.to_bits()),
                | Value::I(x) => slots.push(x.cast_unsigned()),
                | Value::B(x) => slots.push(u64::from(*x)),
                | Value::Read(s) => {
                    slots.push(s.as_ptr() as usize as u64);
                    slots.push(s.len() as u64);
                },
                | Value::Write(s) => {
                    slots.push(s.as_mut_ptr() as usize as u64);
                    slots.push(s.len() as u64);
                },
            }
        }
        let mut out = vec![0_u64; f.results.len().max(1)];
        let status = abi::lir_call(entry, &slots, &mut out);
        if status != 0 {
            return Ok(Outcome::Bailed(status));
        }
        Ok(Outcome::Returned(
            f.results
                .iter()
                .zip(&out)
                .map(|(t, &bits)| {
                    match t {
                        | Ty::F64 => Value::F(f64::from_bits(bits)),
                        | Ty::I64 => Value::I(bits.cast_signed()),
                        | Ty::Bool | Ty::Slice => Value::B(bits & 0xff != 0),
                    }
                })
                .collect(),
        ))
    }
}

/// Emits the host trampoline of `f`.
fn trampoline(
    b: &mut ClBuilder<'_>,
    f: &super::lir::Function,
    callee: cranelift_codegen::ir::FuncRef,
    pointer: Type,
) {
    let flags = MemFlagsData::trusted();
    let entry = b.create_block();
    b.append_block_params_for_function_params(entry);
    b.switch_to_block(entry);
    b.seal_block(entry);
    let params = b.block_params(entry).to_vec();
    let (args_ptr, out_ptr) = (params[0], params[1]);
    let status_slot =
        b.create_sized_stack_slot(StackSlotData::new(StackSlotKind::ExplicitSlot, 8, 3));
    let zero32 = b.ins().iconst(types::I32, 0);
    b.ins().stack_store(pointer, zero32, status_slot, 0);
    let ctx = b.ins().stack_addr(pointer, status_slot, 0);
    let depth = b.ins().iconst(types::I64, 0);
    let mut call_args = vec![ctx, depth];
    let mut slot = 0_i32;
    for t in native_params(f) {
        let raw = b.ins().load(types::I64, flags, args_ptr, slot * 8);
        let v = if t == types::F64 {
            b.ins().bitcast(types::F64, MemFlagsData::new(), raw)
        } else if t == types::I8 {
            b.ins().ireduce(types::I8, raw)
        } else {
            raw
        };
        call_args.push(v);
        slot += 1;
    }
    if f.results.len() > 1 {
        call_args.push(out_ptr);
    }
    let call = b.ins().call(callee, &call_args);
    if let [one] = f.results.as_slice() {
        let result = b.inst_results(call)[0];
        let bits = match one {
            | Ty::F64 => b.ins().bitcast(types::I64, MemFlagsData::new(), result),
            | Ty::Bool => b.ins().uextend(types::I64, result),
            | Ty::I64 | Ty::Slice => result,
        };
        b.ins().store(flags, bits, out_ptr, 0);
    }
    let status = b.ins().stack_load(pointer, types::I32, status_slot, 0);
    b.ins().return_(&[status]);
}

/// Code generation state for one function.
struct Gen<'a, 'b> {
    b: &'a mut ClBuilder<'b>,
    f: &'a super::lir::Function,
    refs: &'a HashMap<&'static str, cranelift_codegen::ir::FuncRef>,
    callees: &'a [cranelift_codegen::ir::FuncRef],
    pointer: Type,
    /// Per register: its variable (the address for a slice) and, for a
    /// slice, the length variable.
    vars: Vec<(Variable, Option<Variable>)>,
    blocks: Vec<ClBlock>,
    ctx_ptr: Option<ClValue>,
    depth: Option<ClValue>,
    ret_ptr: Option<ClValue>,
    /// The shared bail block: parameter = code (i32).
    bail: Option<ClBlock>,
}

impl Gen<'_, '_> {
    fn var(
        &self,
        r: Reg,
    ) -> Result<Variable, BackendError> {
        self.vars
            .get(r.0 as usize)
            .map(|v| v.0)
            .ok_or_else(|| error(format!("register r{}", r.0)))
    }

    fn get(
        &mut self,
        r: Reg,
    ) -> Result<ClValue, BackendError> {
        let v = self.var(r)?;
        Ok(self.b.use_var(v))
    }

    fn slice(
        &mut self,
        r: Reg,
    ) -> Result<(ClValue, ClValue), BackendError> {
        let (addr, len) = self
            .vars
            .get(r.0 as usize)
            .copied()
            .ok_or_else(|| error(format!("register r{}", r.0)))?;
        let len = len.ok_or_else(|| error("not a slice"))?;
        Ok((self.b.use_var(addr), self.b.use_var(len)))
    }

    fn ty(
        &self,
        r: Reg,
    ) -> Ty {
        self.f.regs.get(r.0 as usize).copied().unwrap_or(Ty::F64)
    }

    fn zero_of(
        &mut self,
        t: Type,
    ) -> ClValue {
        if t == types::F64 {
            self.b.ins().f64const(0.0)
        } else {
            self.b.ins().iconst(t, 0)
        }
    }

    /// Returns zeros of the function's result types (the value of a bailed
    /// call, which nobody reads).
    fn return_zeros(&mut self) {
        match self.f.results.as_slice() {
            | [one] => {
                let z = self.zero_of(cl_type(*one));
                self.b.ins().return_(&[z]);
            },
            | _ => {
                self.b.ins().return_(&[]);
            },
        }
    }

    fn bail_block(&mut self) -> ClBlock {
        if let Some(block) = self.bail {
            return block;
        }
        let block = self.b.create_block();
        self.b.append_block_param(block, types::I32);
        self.bail = Some(block);
        block
    }

    /// Jumps to the bail block with `code` when `cond` (an `i8`) is set and
    /// continues in a fresh block otherwise.
    fn bail_if(
        &mut self,
        cond: ClValue,
        code: u32,
    ) {
        let bail = self.bail_block();
        let next = self.b.create_block();
        let code = self.b.ins().iconst(types::I32, i64::from(code));
        self.b.ins().brif(cond, bail, &[code.into()], next, &[]);
        self.b.seal_block(next);
        self.b.switch_to_block(next);
    }

    fn call_shim(
        &mut self,
        name: &str,
        args: &[ClValue],
    ) -> Result<ClValue, BackendError> {
        let f = *self.refs.get(name).ok_or_else(|| error(name))?;
        let call = self.b.ins().call(f, args);
        self.b
            .inst_results(call)
            .first()
            .copied()
            .ok_or_else(|| error(name))
    }

    #[allow(clippy::too_many_lines)]
    fn function(&mut self) -> Result<(), BackendError> {
        // Variables for every register.
        for &t in &self.f.regs {
            let main = self.b.declare_var(cl_type(t));
            let len = (t == Ty::Slice).then(|| self.b.declare_var(types::I64));
            self.vars.push((main, len));
        }
        for _ in &self.f.blocks {
            let block = self.b.create_block();
            self.blocks.push(block);
        }
        // The native entry: unpack parameters, zero every register.
        let entry = self.b.create_block();
        self.b.append_block_params_for_function_params(entry);
        self.b.switch_to_block(entry);
        self.b.seal_block(entry);
        let params = self.b.block_params(entry).to_vec();
        self.ctx_ptr = params.first().copied();
        self.depth = params.get(1).copied();
        if self.f.results.len() > 1 {
            self.ret_ptr = params.last().copied();
        }
        for (r, &t) in self.f.regs.clone().iter().enumerate() {
            let (main, len) = self.vars[r];
            let z = self.zero_of(cl_type(t));
            self.b.def_var(main, z);
            if let Some(len) = len {
                let z = self.b.ins().iconst(types::I64, 0);
                self.b.def_var(len, z);
            }
        }
        let mut k = 2_usize;
        let entry_params = self
            .f
            .blocks
            .first()
            .map(|e| e.params.clone())
            .unwrap_or_default();
        for p in entry_params {
            let (main, len) = self.vars[p.0 as usize];
            self.b.def_var(main, params[k]);
            k += 1;
            if let Some(len) = len {
                self.b.def_var(len, params[k]);
                k += 1;
            }
        }
        // Depth limit.
        let depth = self.depth.ok_or_else(|| error("depth"))?;
        let too_deep = self.b.ins().icmp_imm_s(
            IntCC::SignedGreaterThanOrEqual,
            depth,
            i64::from(MAX_CALL_DEPTH),
        );
        self.bail_if(too_deep, BAIL_DEPTH);
        let first = self.blocks[0];
        self.b.ins().jump(first, &[]);

        for (index, block) in self.f.blocks.iter().enumerate() {
            let cl = self.blocks[index];
            self.b.switch_to_block(cl);
            for (dst, inst) in &block.insts {
                self.inst(*dst, inst)?;
            }
            self.term(&block.term)?;
        }
        if let Some(bail) = self.bail {
            self.b.switch_to_block(bail);
            let code = self.b.block_params(bail)[0];
            let ctx = self.ctx_ptr.ok_or_else(|| error("ctx"))?;
            self.b.ins().store(MemFlagsData::trusted(), code, ctx, 0);
            self.return_zeros();
        }
        self.b.seal_all_blocks();
        Ok(())
    }

    fn assign(
        &mut self,
        targets: &[Reg],
        args: &[Reg],
    ) -> Result<(), BackendError> {
        // Parallel assignment: read every argument before writing.
        let mut values = Vec::new();
        for &a in args {
            if self.ty(a) == Ty::Slice {
                let (p, l) = self.slice(a)?;
                values.push((p, Some(l)));
            } else {
                values.push((self.get(a)?, None));
            }
        }
        for (&t, (v, l)) in targets.iter().zip(values) {
            let (main, len) = self.vars[t.0 as usize];
            self.b.def_var(main, v);
            if let (Some(len), Some(l)) = (len, l) {
                self.b.def_var(len, l);
            }
        }
        Ok(())
    }

    fn term(
        &mut self,
        term: &Term,
    ) -> Result<(), BackendError> {
        match term {
            | Term::Jump(t, args) => {
                let params = self.f.blocks[t.0 as usize].params.clone();
                self.assign(&params, args)?;
                let target = self.blocks[t.0 as usize];
                self.b.ins().jump(target, &[]);
            },
            | Term::Branch(c, (t, ta), (e, ea)) => {
                let cond = self.get(*c)?;
                // Each edge gets its own block for its assignments.
                let then_edge = self.b.create_block();
                let else_edge = self.b.create_block();
                self.b.ins().brif(cond, then_edge, &[], else_edge, &[]);
                for (edge, target, args) in [(then_edge, t, ta), (else_edge, e, ea)] {
                    self.b.seal_block(edge);
                    self.b.switch_to_block(edge);
                    let params = self.f.blocks[target.0 as usize].params.clone();
                    self.assign(&params, args)?;
                    let block = self.blocks[target.0 as usize];
                    self.b.ins().jump(block, &[]);
                }
            },
            | Term::Return(values) => {
                match self.f.results.len() {
                    | 0 => {
                        self.b.ins().return_(&[]);
                    },
                    | 1 => {
                        let v = self.get(values[0])?;
                        self.b.ins().return_(&[v]);
                    },
                    | _ => {
                        let ret = self.ret_ptr.ok_or_else(|| error("ret"))?;
                        for (k, &r) in values.iter().enumerate() {
                            let v = self.get(r)?;
                            let bits = match self.ty(r) {
                                | Ty::F64 => {
                                    self.b.ins().bitcast(types::I64, MemFlagsData::new(), v)
                                },
                                | Ty::Bool => self.b.ins().uextend(types::I64, v),
                                | Ty::I64 | Ty::Slice => v,
                            };
                            let offset = i32::try_from(k * 8).map_err(error)?;
                            self.b
                                .ins()
                                .store(MemFlagsData::trusted(), bits, ret, offset);
                        }
                        self.b.ins().return_(&[]);
                    },
                }
            },
            | Term::Bail(code) => {
                let bail = self.bail_block();
                let code = self.b.ins().iconst(types::I32, i64::from(*code));
                self.b.ins().jump(bail, &[code.into()]);
            },
        }
        Ok(())
    }

    #[allow(clippy::too_many_lines)]
    fn inst(
        &mut self,
        dst: Reg,
        inst: &Inst,
    ) -> Result<(), BackendError> {
        let value = match inst {
            | Inst::F(x) => self.b.ins().f64const(*x),
            | Inst::I(x) => self.b.ins().iconst(types::I64, *x),
            | Inst::B(x) => self.b.ins().iconst(types::I8, i64::from(*x)),
            | Inst::Copy(a) => {
                if self.ty(*a) == Ty::Slice {
                    let (p, l) = self.slice(*a)?;
                    let (main, len) = self.vars[dst.0 as usize];
                    self.b.def_var(main, p);
                    if let Some(len) = len {
                        self.b.def_var(len, l);
                    }
                    return Ok(());
                }
                self.get(*a)?
            },
            | Inst::FAdd(a, c) => {
                let (x, y) = (self.get(*a)?, self.get(*c)?);
                self.b.ins().fadd(x, y)
            },
            | Inst::FSub(a, c) => {
                let (x, y) = (self.get(*a)?, self.get(*c)?);
                self.b.ins().fsub(x, y)
            },
            | Inst::FMul(a, c) => {
                let (x, y) = (self.get(*a)?, self.get(*c)?);
                self.b.ins().fmul(x, y)
            },
            | Inst::FDiv(a, c) => {
                let (x, y) = (self.get(*a)?, self.get(*c)?);
                self.b.ins().fdiv(x, y)
            },
            | Inst::FNeg(a) => {
                let x = self.get(*a)?;
                self.b.ins().fneg(x)
            },
            | Inst::FPowI(a, n) => {
                let x = self.get(*a)?;
                super::jit::emit_powi(self.b, x, *n)
            },
            | Inst::FPow(a, c) => {
                let (x, y) = (self.get(*a)?, self.get(*c)?);
                self.call_shim("rssn_pow", &[x, y])?
            },
            | Inst::FAtan2(a, c) => {
                let (x, y) = (self.get(*a)?, self.get(*c)?);
                self.call_shim("rssn_atan2", &[x, y])?
            },
            | Inst::FUnary(i, a) => {
                let x = self.get(*a)?;
                match i {
                    | Intrinsic::Sqrt => self.b.ins().sqrt(x),
                    | Intrinsic::Abs => self.b.ins().fabs(x),
                    | Intrinsic::Floor => self.b.ins().floor(x),
                    | Intrinsic::Ceil => self.b.ins().ceil(x),
                    | other => {
                        let (name, _) =
                            super::jit::library(*other).ok_or_else(|| error("intrinsic"))?;
                        self.call_shim(name, &[x])?
                    },
                }
            },
            | Inst::FMin(a, c) | Inst::FMax(a, c) => {
                let (x, y) = (self.get(*a)?, self.get(*c)?);
                if matches!(inst, Inst::FMin(..)) {
                    self.b.ins().fmin(x, y)
                } else {
                    self.b.ins().fmax(x, y)
                }
            },
            | Inst::IAdd(a, c) | Inst::ISub(a, c) | Inst::IMul(a, c) => {
                let (x, y) = (self.get(*a)?, self.get(*c)?);
                let (v, overflow) = match inst {
                    | Inst::IAdd(..) => self.b.ins().sadd_overflow(x, y),
                    | Inst::ISub(..) => self.b.ins().ssub_overflow(x, y),
                    | _ => self.b.ins().smul_overflow(x, y),
                };
                self.bail_if(overflow, BAIL_OVERFLOW);
                v
            },
            | Inst::IDiv(a, c) | Inst::IRem(a, c) => {
                let (x, y) = (self.get(*a)?, self.get(*c)?);
                let zero = self.b.ins().icmp_imm_s(IntCC::Equal, y, 0);
                let min = self.b.ins().icmp_imm_s(IntCC::Equal, x, i64::MIN);
                let minus_one = self.b.ins().icmp_imm_s(IntCC::Equal, y, -1);
                let wrap = self.b.ins().band(min, minus_one);
                let bad = self.b.ins().bor(zero, wrap);
                self.bail_if(bad, BAIL_OVERFLOW);
                if matches!(inst, Inst::IDiv(..)) {
                    self.b.ins().sdiv(x, y)
                } else {
                    self.b.ins().srem(x, y)
                }
            },
            | Inst::INeg(a) => {
                let x = self.get(*a)?;
                let min = self.b.ins().icmp_imm_s(IntCC::Equal, x, i64::MIN);
                self.bail_if(min, BAIL_OVERFLOW);
                self.b.ins().ineg(x)
            },
            | Inst::IMin(a, c) => {
                let (x, y) = (self.get(*a)?, self.get(*c)?);
                self.b.ins().smin(x, y)
            },
            | Inst::IMax(a, c) => {
                let (x, y) = (self.get(*a)?, self.get(*c)?);
                self.b.ins().smax(x, y)
            },
            | Inst::FCmp(k, a, c) => {
                let (x, y) = (self.get(*a)?, self.get(*c)?);
                self.b.ins().fcmp(float_cc(*k), x, y)
            },
            | Inst::ICmp(k, a, c) => {
                let (x, y) = (self.get(*a)?, self.get(*c)?);
                self.b.ins().icmp(int_cc(*k), x, y)
            },
            | Inst::Not(a) => {
                let x = self.get(*a)?;
                self.b.ins().bxor_imm_u(x, 1)
            },
            | Inst::And(a, c) => {
                let (x, y) = (self.get(*a)?, self.get(*c)?);
                self.b.ins().band(x, y)
            },
            | Inst::Or(a, c) => {
                let (x, y) = (self.get(*a)?, self.get(*c)?);
                self.b.ins().bor(x, y)
            },
            | Inst::IToF(a) => {
                let x = self.get(*a)?;
                self.b.ins().fcvt_from_sint(types::F64, x)
            },
            | Inst::FToI(a) => {
                let x = self.get(*a)?;
                let lo = self.b.ins().f64const(-9_223_372_036_854_775_808.0);
                let hi = self.b.ins().f64const(9_223_372_036_854_775_808.0);
                let ge = self.b.ins().fcmp(FloatCC::GreaterThanOrEqual, x, lo);
                let lt = self.b.ins().fcmp(FloatCC::LessThan, x, hi);
                let ok = self.b.ins().band(ge, lt);
                let bad = self.b.ins().bxor_imm_u(ok, 1);
                self.bail_if(bad, BAIL_CONVERT);
                self.b.ins().fcvt_to_sint(types::I64, x)
            },
            | Inst::Select(c, a, d) => {
                let cond = self.get(*c)?;
                let (x, y) = (self.get(*a)?, self.get(*d)?);
                self.b.ins().select(cond, x, y)
            },
            | Inst::Len(s) => self.slice(*s)?.1,
            | Inst::Load(s, i) => {
                let (ptr, len) = self.slice(*s)?;
                let index = self.get(*i)?;
                let inside = self.b.ins().icmp(IntCC::UnsignedLessThan, index, len);
                let load = self.b.create_block();
                let done = self.b.create_block();
                self.b.append_block_param(done, types::F64);
                let nan = self.b.ins().f64const(f64::NAN);
                self.b.ins().brif(inside, load, &[], done, &[nan.into()]);
                self.b.seal_block(load);
                self.b.switch_to_block(load);
                let offset = self.b.ins().imul_imm_s(index, 8);
                let address = self.b.ins().iadd(ptr, offset);
                let v = self
                    .b
                    .ins()
                    .load(types::F64, MemFlagsData::trusted(), address, 0);
                self.b.ins().jump(done, &[v.into()]);
                self.b.seal_block(done);
                self.b.switch_to_block(done);
                self.b.block_params(done)[0]
            },
            | Inst::Store(s, i, v) => {
                let (ptr, len) = self.slice(*s)?;
                let index = self.get(*i)?;
                let value = self.get(*v)?;
                let inside = self.b.ins().icmp(IntCC::UnsignedLessThan, index, len);
                let store = self.b.create_block();
                let done = self.b.create_block();
                self.b.ins().brif(inside, store, &[], done, &[]);
                self.b.seal_block(store);
                self.b.switch_to_block(store);
                let offset = self.b.ins().imul_imm_s(index, 8);
                let address = self.b.ins().iadd(ptr, offset);
                self.b
                    .ins()
                    .store(MemFlagsData::trusted(), value, address, 0);
                self.b.ins().jump(done, &[]);
                self.b.seal_block(done);
                self.b.switch_to_block(done);
                self.b.ins().f64const(0.0)
            },
            | Inst::Call(callee, args) => {
                let ctx = self.ctx_ptr.ok_or_else(|| error("ctx"))?;
                let depth = self.depth.ok_or_else(|| error("depth"))?;
                let next = self.b.ins().iadd_imm_s(depth, 1);
                let mut values = vec![ctx, next];
                for &a in args {
                    if self.ty(a) == Ty::Slice {
                        let (p, l) = self.slice(a)?;
                        values.push(p);
                        values.push(l);
                    } else {
                        values.push(self.get(a)?);
                    }
                }
                let target = *self
                    .callees
                    .get(callee.0 as usize)
                    .ok_or_else(|| error("callee"))?;
                let call = self.b.ins().call(target, &values);
                let result = self
                    .b
                    .inst_results(call)
                    .first()
                    .copied()
                    .ok_or_else(|| error("call result"))?;
                // A bail inside the callee: propagate it.
                let status = self
                    .b
                    .ins()
                    .load(types::I32, MemFlagsData::trusted(), ctx, 0);
                let failed = self.b.ins().icmp_imm_s(IntCC::NotEqual, status, 0);
                let propagate = self.b.create_block();
                let next_block = self.b.create_block();
                self.b.ins().brif(failed, propagate, &[], next_block, &[]);
                self.b.seal_block(propagate);
                self.b.switch_to_block(propagate);
                self.return_zeros();
                self.b.seal_block(next_block);
                self.b.switch_to_block(next_block);
                result
            },
            | Inst::Eval(function, args) => {
                let n = u32::try_from(args.len().max(1) * 8).map_err(error)?;
                let slot = self.b.create_sized_stack_slot(StackSlotData::new(
                    StackSlotKind::ExplicitSlot,
                    n,
                    3,
                ));
                for (k, &a) in args.iter().enumerate() {
                    let v = self.get(a)?;
                    let offset = i32::try_from(k * 8).map_err(error)?;
                    self.b.ins().stack_store(self.pointer, v, slot, offset);
                }
                let base = self.b.ins().stack_addr(self.pointer, slot, 0);
                let address = self.b.ins().iconst(self.pointer, *function as usize as i64);
                let count = self.b.ins().iconst(self.pointer, args.len() as i64);
                self.call_shim("rssn_eval", &[address, base, count])?
            },
            | Inst::Extern1(function, a) => {
                let x = self.get(*a)?;
                let address = self.b.ins().iconst(self.pointer, *function as usize as i64);
                self.call_shim("rssn_extern1", &[address, x])?
            },
            | Inst::Extern2(function, a, c) => {
                let (x, y) = (self.get(*a)?, self.get(*c)?);
                let address = self.b.ins().iconst(self.pointer, *function as usize as i64);
                self.call_shim("rssn_extern2", &[address, x, y])?
            },
        };
        let main = self.var(dst)?;
        self.b.def_var(main, value);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::super::lir::FunctionBuilder;
    use super::super::lir::Value;
    use super::*;

    fn agree(
        module: &Module,
        func: FuncId,
        args: &[Value<'static>],
    ) {
        let compiled = LirModule::compile(module).unwrap_or_else(|e| panic!("{e}"));
        let copy = |a: &Value<'static>| {
            match a {
                | Value::F(x) => Value::F(*x),
                | Value::I(x) => Value::I(*x),
                | Value::B(x) => Value::B(*x),
                | Value::Read(s) => Value::Read(s),
                | Value::Write(_) => Value::B(false),
            }
        };
        let mut a1: Vec<Value<'_>> = args.iter().map(copy).collect();
        let mut a2: Vec<Value<'_>> = args.iter().map(copy).collect();
        let interpreted = module.interpret(func, &mut a1, u64::MAX);
        let native = compiled.call(func, &mut a2);
        let bits = |o: &Result<Outcome, LirError>| -> String {
            match o {
                | Ok(Outcome::Returned(v)) => {
                    v.iter()
                        .map(|x| {
                            match x {
                                | Value::F(f) => format!("{:x}", f.to_bits()),
                                | other => format!("{other:?}"),
                            }
                        })
                        .collect::<Vec<_>>()
                        .join(",")
                },
                | other => format!("{other:?}"),
            }
        };
        assert_eq!(bits(&interpreted), bits(&native), "{args:?}");
    }

    #[test]
    fn agrees_with_the_interpreter() {
        // Recursive factorial with checked overflow.
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
        // A float loop: x_{k+1} = cos(x_k) + x_k^3 / 7, k < n, min/max mixed in.
        let mut g = FunctionBuilder::new("iterate", &[Ty::F64, Ty::I64], &[Ty::F64]);
        let (x0, count) = (g.param(0), g.param(1));
        let (head, body, exit) = (g.new_block(), g.new_block(), g.new_block());
        let x = g.block_param(head, Ty::F64);
        let k = g.block_param(head, Ty::I64);
        let zero = g.ins(Inst::I(0));
        g.jump(head, &[x0, zero]);
        g.switch_to(head);
        let more = g.ins(Inst::ICmp(Cmp::Lt, k, count));
        g.branch(more, body, &[], exit, &[]);
        g.switch_to(body);
        let c = g.ins(Inst::FUnary(Intrinsic::Cos, x));
        let cube = g.ins(Inst::FPowI(x, 3));
        let seven = g.ins(Inst::F(7.0));
        let q = g.ins(Inst::FDiv(cube, seven));
        let y = g.ins(Inst::FAdd(c, q));
        let lo = g.ins(Inst::F(-3.0));
        let y = g.ins(Inst::FMax(y, lo));
        let one = g.ins(Inst::I(1));
        let k2 = g.ins(Inst::IAdd(k, one));
        g.jump(head, &[y, k2]);
        g.switch_to(exit);
        let t = g.ins(Inst::FToI(x));
        let back = g.ins(Inst::IToF(t));
        let out = g.ins(Inst::FSub(x, back));
        g.ret(&[out]);
        let module = Module {
            functions: vec![
                f.finish().unwrap_or_else(|e| panic!("{e}")),
                g.finish().unwrap_or_else(|e| panic!("{e}")),
            ],
        };
        for n in [0, 1, 5, 20, 21, -3] {
            agree(&module, FuncId(0), &[Value::I(n)]);
        }
        for (x, k) in [(0.5, 10), (1.0, 1000), (-2.0, 3), (f64::NAN, 2)] {
            agree(&module, FuncId(1), &[Value::F(x), Value::I(k)]);
        }
        agree(&module, FuncId(0), &[Value::I(100_000)]);
    }

    #[test]
    fn slices_and_multiple_results() {
        // (sum, max) of a slice, written also into out[0], out[1].
        let mut f = FunctionBuilder::new("stats", &[Ty::Slice, Ty::Slice], &[Ty::F64, Ty::F64]);
        let (xs, out) = (f.param(0), f.param(1));
        let (head, body, exit) = (f.new_block(), f.new_block(), f.new_block());
        let i = f.block_param(head, Ty::I64);
        let s = f.block_param(head, Ty::F64);
        let m = f.block_param(head, Ty::F64);
        let zero = f.ins(Inst::I(0));
        let fz = f.ins(Inst::F(0.0));
        let ninf = f.ins(Inst::F(f64::NEG_INFINITY));
        f.jump(head, &[zero, fz, ninf]);
        f.switch_to(head);
        let n = f.ins(Inst::Len(xs));
        let more = f.ins(Inst::ICmp(Cmp::Lt, i, n));
        f.branch(more, body, &[], exit, &[]);
        f.switch_to(body);
        let x = f.ins(Inst::Load(xs, i));
        let s2 = f.ins(Inst::FAdd(s, x));
        let m2 = f.ins(Inst::FMax(m, x));
        let one = f.ins(Inst::I(1));
        let i2 = f.ins(Inst::IAdd(i, one));
        f.jump(head, &[i2, s2, m2]);
        f.switch_to(exit);
        let z = f.ins(Inst::I(0));
        let o = f.ins(Inst::I(1));
        let far = f.ins(Inst::I(-5));
        f.ins(Inst::Store(out, z, s));
        f.ins(Inst::Store(out, o, m));
        f.ins(Inst::Store(out, far, m));
        f.ret(&[s, m]);
        let module = Module {
            functions: vec![f.finish().unwrap_or_else(|e| panic!("{e}"))],
        };
        let compiled = LirModule::compile(&module).unwrap_or_else(|e| panic!("{e}"));
        let xs = [1.0, 4.0, -2.0];
        let mut out = [0.0; 2];
        let got = compiled.call(FuncId(0), &mut [Value::Read(&xs), Value::Write(&mut out)]);
        assert_eq!(
            got,
            Ok(Outcome::Returned(vec![Value::F(3.0), Value::F(4.0)]))
        );
        assert_eq!(out, [3.0, 4.0]);
        // Read-only output is refused before any code runs.
        let refused = compiled.call(FuncId(0), &mut [Value::Read(&xs), Value::Read(&[0.0; 2])]);
        assert!(refused.is_err());
    }
}
