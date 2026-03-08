//! JIT Engine implementation using Cranelift.

use std::collections::HashMap;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;

#[cfg(feature = "jit")]
use cranelift_codegen::Context;
#[cfg(feature = "jit")]
use cranelift_codegen::ir::AbiParam;
#[cfg(feature = "jit")]
use cranelift_codegen::ir::InstBuilder;
#[cfg(feature = "jit")]
use cranelift_codegen::ir::MemFlags;
#[cfg(feature = "jit")]
use cranelift_codegen::ir::StackSlotData;
#[cfg(feature = "jit")]
use cranelift_codegen::ir::StackSlotKind;
#[cfg(feature = "jit")]
use cranelift_codegen::ir::Type;
#[cfg(feature = "jit")]
use cranelift_codegen::ir::Value;
#[cfg(feature = "jit")]
use cranelift_codegen::ir::condcodes::FloatCC;
#[cfg(feature = "jit")]
use cranelift_codegen::ir::condcodes::IntCC;
#[cfg(feature = "jit")]
use cranelift_codegen::ir::types;
#[cfg(feature = "jit")]
use cranelift_frontend::FunctionBuilder;
#[cfg(feature = "jit")]
use cranelift_frontend::FunctionBuilderContext;
#[cfg(feature = "jit")]
use cranelift_jit::JITBuilder;
#[cfg(feature = "jit")]
use cranelift_jit::JITModule;
#[cfg(feature = "jit")]
use cranelift_module::Linkage;
#[cfg(feature = "jit")]
use cranelift_module::Module;

use crate::jit::instructions::Instruction;
use crate::jit::instructions::JitType;

/// The JIT Engine for compiling and executing code.
/// Represents a region of memory accessible to the JIT sandbox.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct MemoryRegion {
    /// Base address of the memory region.
    pub base: *mut u8,
    /// Size of the memory region.
    pub size: usize,
}

/// The Sandbox Context holding pointers and lengths of memory regions and allowed functions.
#[repr(C)]
#[derive(Debug, Clone)]
pub struct SandboxContext {
    /// Pointer to the array of MemoryRegion structs.
    pub regions: *const MemoryRegion,
    /// Number of memory regions registered.
    pub regions_count: usize,
    /// Pointer to the array of allowed function pointers.
    pub allowed_calls: *const *const u8,
    /// Number of allowed function pointers registered.
    pub allowed_calls_count: usize,
}

/// The JIT Engine for compiling and executing code.
pub struct JitEngine {
    #[cfg(feature = "jit")]
    builder_context: FunctionBuilderContext,
    #[cfg(feature = "jit")]
    ctx: Context,
    #[cfg(feature = "jit")]
    module: JITModule,

    #[cfg(feature = "jit")]
    custom_ops: HashMap<u32, CustomOpData>,

    memory_regions: Vec<MemoryRegion>,
    allowed_calls: Vec<*const u8>,

    function_counter: AtomicUsize,
}

/// Holds the data for a custom operation.
#[cfg(feature = "jit")]
#[derive(Clone)]
pub struct CustomOpData {
    /// A pointer to the custom function.
    pub func_ptr: *const u8,
    /// The number of arguments the custom function takes.
    pub arg_count: usize,
}

impl Default for JitEngine {
    fn default() -> Self {
        Self::new()
    }
}

impl JitEngine {
    /// Creates a new JIT engine.
    ///
    /// # Panics
    ///
    /// Panics if the `JITBuilder` fails to initialize.
    /// Registers a memory region to the sandbox, returning the region ID.
    #[must_use]
    pub fn register_memory_region(
        &mut self,
        base: *mut u8,
        size: usize,
    ) -> u16 {
        let id = self.memory_regions.len() as u16;

        self.memory_regions.push(MemoryRegion { base, size });

        id
    }

    /// Clears all registered memory regions.
    pub fn clear_memory_regions(&mut self) {
        self.memory_regions.clear();
    }

    /// Appends a function pointer to the list of allowed JIT calls.
    pub fn allow_call_target(
        &mut self,
        func_ptr: *const u8,
    ) {
        self.allowed_calls.push(func_ptr);
    }

    /// Builds the `SandboxContext` struct from the current registered lists.
    pub const fn build_sandbox_context(&self) -> SandboxContext {
        SandboxContext {
            regions: self.memory_regions.as_ptr(),
            regions_count: self.memory_regions.len(),
            allowed_calls: self.allowed_calls.as_ptr(),
            allowed_calls_count: self.allowed_calls.len(),
        }
    }

    /// Creates a new JIT engine.
    ///
    /// # Panics
    ///
    /// Panics if the `JITBuilder` fails to initialize.
    pub fn new() -> Self {
        let function_counter = AtomicUsize::new(0);

        #[cfg(feature = "jit")]
        {
            let builder = JITBuilder::new(cranelift_module::default_libcall_names()).unwrap();

            let module = JITModule::new(builder);

            Self {
                builder_context: FunctionBuilderContext::new(),
                ctx: module.make_context(),
                module,
                custom_ops: HashMap::new(),
                memory_regions: Vec::new(),
                allowed_calls: Vec::new(),
                function_counter,
            }
        }

        #[cfg(not(feature = "jit"))]
        {
            Self {
                memory_regions: Vec::new(),
                allowed_calls: Vec::new(),
                function_counter,
            }
        }
    }

    /// Registers a custom operation handler.
    ///
    /// # Safety
    ///
    /// The caller must ensure that `func_ptr` is a valid pointer to a function
    /// that can be safely called with `arg_count` arguments. This pointer must remain
    /// valid for the entire lifetime of the `JitEngine` and any compiled code that uses it.
    pub unsafe fn register_custom_op(
        &mut self,
        opcode: u32,
        func_ptr: *const u8,
        arg_count: usize,
    ) {
        #[cfg(feature = "jit")]
        {
            self.custom_ops
                .insert(opcode, CustomOpData { func_ptr, arg_count });
        }
    }

    /// Compiles the given instructions into a function.
    /// The function takes no arguments and returns an f64.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - The `jit` feature is not enabled.
    /// - The function cannot be declared or defined in the JIT module.
    /// - The module fails to finalize definitions.
    ///
    /// # Safety
    ///
    /// The caller must ensure that the provided instructions are valid and do not cause
    /// undefined behavior when executed (e.g., stack underflow, invalid memory access).
    /// The returned pointer must be used correctly according to the expected signature.
    pub unsafe fn compile(
        &mut self,
        instructions: &[Instruction],
    ) -> Result<*const u8, String> {
        #[cfg(feature = "jit")]
        {
            let id = self.function_counter.fetch_add(1, Ordering::SeqCst);

            let name = format!("jit_fn_{id}");

            // Setup signature: (i64) -> f64
            self.ctx.func.signature.params.clear();

            self.ctx
                .func
                .signature
                .params
                .push(AbiParam::new(types::I64));

            self.ctx.func.signature.returns.clear();

            self.ctx
                .func
                .signature
                .returns
                .push(AbiParam::new(types::F64));

            let mut builder = FunctionBuilder::new(&mut self.ctx.func, &mut self.builder_context);

            // Create blocks for all labels
            let mut blocks = HashMap::new();

            let entry_block = builder.create_block();

            let trap_exit_block = builder.create_block();

            blocks.insert(u32::MAX, entry_block);

            for inst in instructions {
                if let Instruction::Label(id) = inst {
                    if !blocks.contains_key(id) {
                        blocks.insert(*id, builder.create_block());
                    }
                }
            }

            builder.append_block_params_for_function_params(entry_block);

            builder.switch_to_block(entry_block);

            // Stack slot setup
            let stack_slot = builder.create_sized_stack_slot(StackSlotData::new(
                StackSlotKind::ExplicitSlot,
                2048,
                4,
            ));

            let stack_ptr = builder.ins().stack_addr(types::I64, stack_slot, 0);

            let sp_var = builder.declare_var(types::I64);

            let var_tmp = builder.ins().iconst(types::I64, 0);

            builder.def_var(sp_var, var_tmp);

            // Macros for stack manipulation
            macro_rules! pop {
                () => {{
                    let continue_block_pop = builder.create_block();
                    let sp = builder.use_var(sp_var);
                    let cmp_under = builder.ins().icmp_imm(IntCC::SignedLessThanOrEqual, sp, 0);
                    builder
                        .ins()
                        .brif(cmp_under, trap_exit_block, &[], continue_block_pop, &[]);
                    builder.switch_to_block(continue_block_pop);
                    let old_sp = builder.ins().iadd_imm(sp, -8);
                    let addr = builder.ins().iadd(stack_ptr, old_sp);
                    let val = builder.ins().load(types::I64, MemFlags::new(), addr, 0);
                    builder.def_var(sp_var, old_sp);
                    val
                }};
            }

            macro_rules! push {
                ($val:expr) => {{
                    let val = $val;
                    let sp = builder.use_var(sp_var);
                    let cmp_over =
                        builder
                            .ins()
                            .icmp_imm(IntCC::SignedGreaterThanOrEqual, sp, 2048);
                    let continue_block_push = builder.create_block();
                    builder
                        .ins()
                        .brif(cmp_over, trap_exit_block, &[], continue_block_push, &[]);
                    builder.switch_to_block(continue_block_push);

                    let addr = builder.ins().iadd(stack_ptr, sp);
                    builder.ins().store(MemFlags::new(), val, addr, 0);
                    let new_sp = builder.ins().iadd_imm(sp, 8);
                    builder.def_var(sp_var, new_sp);
                }};
            }

            // Loop instructions
            for inst in instructions {
                if let Instruction::Label(id) = inst {
                    let block = blocks[id];
                    if !builder.is_unreachable() {
                        builder.ins().jump(block, &[]);
                    }
                    builder.switch_to_block(block);
                }

                match inst {
                    | Instruction::Label(_) => {},

                    | Instruction::ImmI(val) => {
                        let val_v = builder.ins().iconst(types::I64, *val);
                        push!(val_v);
                    },
                    | Instruction::ImmF(val) => {
                        let val_v = builder.ins().f64const(*val);
                        push!(val_v);
                    },

                    | Instruction::Load(ty) => {
                        let addr_v = pop!();
                        let (load_ty, store_ty) = type_to_clif(*ty);
                        let ctx_ptr = builder.block_params(entry_block)[0];
                        let regions_ptr =
                            builder.ins().load(types::I64, MemFlags::new(), ctx_ptr, 0);
                        let regions_count =
                            builder.ins().load(types::I64, MemFlags::new(), ctx_ptr, 8);
                        let region_id = builder.ins().ushr_imm(addr_v, 48);
                        let offset = builder.ins().band_imm(addr_v, 0x0000FFFFFFFFFFFF);
                        let cmp_region = builder.ins().icmp(
                            IntCC::UnsignedGreaterThanOrEqual,
                            region_id,
                            regions_count,
                        );
                        let check_region_ok = builder.create_block();
                        builder
                            .ins()
                            .brif(cmp_region, trap_exit_block, &[], check_region_ok, &[]);
                        builder.switch_to_block(check_region_ok);
                        let region_idx16 = builder.ins().imul_imm(region_id, 16);
                        let region_struct_ptr = builder.ins().iadd(regions_ptr, region_idx16);
                        let base_ptr =
                            builder
                                .ins()
                                .load(types::I64, MemFlags::new(), region_struct_ptr, 0);
                        let size_val =
                            builder
                                .ins()
                                .load(types::I64, MemFlags::new(), region_struct_ptr, 8);
                        let size_req = builder.ins().iadd_imm(offset, i64::from(load_ty.bytes()));
                        let cmp_size =
                            builder
                                .ins()
                                .icmp(IntCC::UnsignedGreaterThan, size_req, size_val);
                        let check_size_ok = builder.create_block();
                        builder
                            .ins()
                            .brif(cmp_size, trap_exit_block, &[], check_size_ok, &[]);
                        builder.switch_to_block(check_size_ok);
                        let phys_addr = builder.ins().iadd(base_ptr, offset);
                        let val = builder.ins().load(load_ty, MemFlags::new(), phys_addr, 0);
                        let val_64 = cast_to_storage(&mut builder, val, load_ty, store_ty);
                        push!(val_64);
                    },

                    | Instruction::Store(ty) => {
                        let val_raw = pop!();
                        let addr_v = pop!();
                        let (mem_ty, _) = type_to_clif(*ty);
                        let ctx_ptr = builder.block_params(entry_block)[0];
                        let regions_ptr =
                            builder.ins().load(types::I64, MemFlags::new(), ctx_ptr, 0);
                        let regions_count =
                            builder.ins().load(types::I64, MemFlags::new(), ctx_ptr, 8);
                        let region_id = builder.ins().ushr_imm(addr_v, 48);
                        let offset = builder.ins().band_imm(addr_v, 0x0000FFFFFFFFFFFF);
                        let cmp_region = builder.ins().icmp(
                            IntCC::UnsignedGreaterThanOrEqual,
                            region_id,
                            regions_count,
                        );
                        let check_region_ok = builder.create_block();
                        builder
                            .ins()
                            .brif(cmp_region, trap_exit_block, &[], check_region_ok, &[]);
                        builder.switch_to_block(check_region_ok);
                        let region_idx16 = builder.ins().imul_imm(region_id, 16);
                        let region_struct_ptr = builder.ins().iadd(regions_ptr, region_idx16);
                        let base_ptr =
                            builder
                                .ins()
                                .load(types::I64, MemFlags::new(), region_struct_ptr, 0);
                        let size_val =
                            builder
                                .ins()
                                .load(types::I64, MemFlags::new(), region_struct_ptr, 8);
                        let size_req = builder.ins().iadd_imm(offset, i64::from(mem_ty.bytes()));
                        let cmp_size =
                            builder
                                .ins()
                                .icmp(IntCC::UnsignedGreaterThan, size_req, size_val);
                        let check_size_ok = builder.create_block();
                        builder
                            .ins()
                            .brif(cmp_size, trap_exit_block, &[], check_size_ok, &[]);
                        builder.switch_to_block(check_size_ok);
                        let phys_addr = builder.ins().iadd(base_ptr, offset);
                        let val_trunc = cast_from_storage(&mut builder, val_raw, mem_ty);
                        builder
                            .ins()
                            .store(MemFlags::new(), val_trunc, phys_addr, 0);
                    },

                    | Instruction::Custom { opcode, payload: _ } => {
                        if let Some(op_data) = self.custom_ops.get(opcode) {
                            let op_addr = op_data.func_ptr as i64;
                            let arg_cnt = op_data.arg_count;
                            let fn_ptr_v = builder.ins().iconst(types::I64, op_addr);

                            let mut args = Vec::new();
                            for _ in 0..arg_cnt {
                                args.push(pop!());
                            }
                            // args are popped in reverse order (top is last arg)

                            let mut sig = self.module.make_signature();
                            for _ in 0..arg_cnt {
                                sig.params.push(AbiParam::new(types::I64));
                            }
                            sig.returns.push(AbiParam::new(types::I64));

                            let sig_ref = builder.import_signature(sig);

                            // args need to be reversed
                            // CallIndirect expects &[Value].
                            // If vector is [y, x], args are y, x.
                            // This means y is arg0. Incorrect.
                            // So we MUST reverse.
                            args.reverse();

                            let call = builder.ins().call_indirect(sig_ref, fn_ptr_v, &args);
                            let res = builder.inst_results(call)[0];
                            push!(res);
                        }
                    },

                    | Instruction::Add(ty)
                    | Instruction::Sub(ty)
                    | Instruction::Mul(ty)
                    | Instruction::Div(ty) => {
                        let rhs_raw = pop!();
                        let lhs_raw = pop!();

                        let (op_ty, _) = type_to_clif(*ty);
                        let rhs = cast_from_storage(&mut builder, rhs_raw, op_ty);
                        let lhs = cast_from_storage(&mut builder, lhs_raw, op_ty);

                        let res = match inst {
                            | Instruction::Add(_) => {
                                if op_ty.is_int() {
                                    builder.ins().iadd(lhs, rhs)
                                } else {
                                    builder.ins().fadd(lhs, rhs)
                                }
                            },
                            | Instruction::Sub(_) => {
                                if op_ty.is_int() {
                                    builder.ins().isub(lhs, rhs)
                                } else {
                                    builder.ins().fsub(lhs, rhs)
                                }
                            },
                            | Instruction::Mul(_) => {
                                if op_ty.is_int() {
                                    builder.ins().imul(lhs, rhs)
                                } else {
                                    builder.ins().fmul(lhs, rhs)
                                }
                            },
                            | Instruction::Div(_) => {
                                if op_ty.is_int() {
                                    builder.ins().sdiv(lhs, rhs)
                                } else {
                                    builder.ins().fdiv(lhs, rhs)
                                }
                            },
                            | _ => unreachable!(),
                        };

                        let res_chk = cast_to_storage(&mut builder, res, op_ty, types::I64);
                        push!(res_chk);
                    },

                    | Instruction::And | Instruction::Or | Instruction::Xor => {
                        let rhs = pop!();
                        let lhs = pop!();
                        let res = match inst {
                            | Instruction::And => builder.ins().band(lhs, rhs),
                            | Instruction::Or => builder.ins().bor(lhs, rhs),
                            | Instruction::Xor => builder.ins().bxor(lhs, rhs),
                            | _ => unreachable!(),
                        };
                        push!(res);
                    },

                    | Instruction::Not => {
                        let val = pop!();
                        let res = builder.ins().bnot(val);
                        push!(res);
                    },

                    | Instruction::Eq(ty)
                    | Instruction::Ne(ty)
                    | Instruction::Lt(ty)
                    | Instruction::Gt(ty)
                    | Instruction::Le(ty)
                    | Instruction::Ge(ty) => {
                        let rhs_raw = pop!();
                        let lhs_raw = pop!();

                        let (op_ty, _) = type_to_clif(*ty);
                        let rhs = cast_from_storage(&mut builder, rhs_raw, op_ty);
                        let lhs = cast_from_storage(&mut builder, lhs_raw, op_ty);

                        let cond = if op_ty.is_int() {
                            let cc = match inst {
                                | Instruction::Eq(_) => IntCC::Equal,
                                | Instruction::Ne(_) => IntCC::NotEqual,
                                | Instruction::Lt(_) => IntCC::SignedLessThan,
                                | Instruction::Gt(_) => IntCC::SignedGreaterThan,
                                | Instruction::Le(_) => IntCC::SignedLessThanOrEqual,
                                | Instruction::Ge(_) => IntCC::SignedGreaterThanOrEqual,
                                | _ => unreachable!(),
                            };
                            builder.ins().icmp(cc, lhs, rhs)
                        } else {
                            let cc = match inst {
                                | Instruction::Eq(_) => FloatCC::Equal,
                                | Instruction::Ne(_) => FloatCC::NotEqual,
                                | Instruction::Lt(_) => FloatCC::LessThan,
                                | Instruction::Gt(_) => FloatCC::GreaterThan,
                                | Instruction::Le(_) => FloatCC::LessThanOrEqual,
                                | Instruction::Ge(_) => FloatCC::GreaterThanOrEqual,
                                | _ => unreachable!(),
                            };
                            builder.ins().fcmp(cc, lhs, rhs)
                        };

                        let res = builder.ins().uextend(types::I64, cond);
                        push!(res);
                    },

                    | Instruction::Jump(target) => {
                        if let Some(blk) = blocks.get(target) {
                            builder.ins().jump(*blk, &[]);
                        }
                    },

                    | Instruction::BranchIfTrue(target) => {
                        let val = pop!();
                        if let Some(blk) = blocks.get(target) {
                            let cond = builder.ins().icmp_imm(IntCC::NotEqual, val, 0);
                            let continue_block = builder.create_block();
                            builder.ins().brif(cond, *blk, &[], continue_block, &[]);
                            builder.switch_to_block(continue_block);
                        }
                    },

                    | Instruction::BranchIfFalse(target) => {
                        let val = pop!();
                        if let Some(blk) = blocks.get(target) {
                            let cond = builder.ins().icmp_imm(IntCC::Equal, val, 0);
                            let continue_block = builder.create_block();
                            builder.ins().brif(cond, *blk, &[], continue_block, &[]);
                            builder.switch_to_block(continue_block);
                        }
                    },

                    | Instruction::Dup => {
                        let val = pop!();
                        push!(val);
                        push!(val);
                    },

                    | Instruction::Drop => {
                        let _ = pop!();
                    },

                    | Instruction::Swap => {
                        let v1 = pop!();
                        let v2 = pop!();
                        push!(v1);
                        push!(v2);
                    },

                    | Instruction::Call(arg_count) => {
                        let fn_ptr = pop!();
                        let ctx_ptr = builder.block_params(entry_block)[0];
                        let allowed_calls_ptr =
                            builder.ins().load(types::I64, MemFlags::new(), ctx_ptr, 16);
                        let allowed_calls_count =
                            builder.ins().load(types::I64, MemFlags::new(), ctx_ptr, 24);

                        let loop_header = builder.create_block();
                        let loop_body = builder.create_block();
                        let loop_exit = builder.create_block();
                        let trap_block = builder.create_block();

                        let iter_var = builder.declare_var(types::I64);
                        let zero_val = builder.ins().iconst(types::I64, 0);
                        builder.def_var(iter_var, zero_val);
                        builder.ins().jump(loop_header, &[]);

                        builder.switch_to_block(loop_header);
                        let idx = builder.use_var(iter_var);
                        let cmp_done = builder.ins().icmp(
                            IntCC::UnsignedGreaterThanOrEqual,
                            idx,
                            allowed_calls_count,
                        );
                        builder
                            .ins()
                            .brif(cmp_done, trap_block, &[], loop_body, &[]);

                        builder.switch_to_block(loop_body);
                        let idx8 = builder.ins().imul_imm(idx, 8);
                        let ptr_addr = builder.ins().iadd(allowed_calls_ptr, idx8);
                        let allowed_ptr =
                            builder.ins().load(types::I64, MemFlags::new(), ptr_addr, 0);
                        let cmp_match = builder.ins().icmp(IntCC::Equal, fn_ptr, allowed_ptr);

                        let next_idx = builder.ins().iadd_imm(idx, 1);
                        builder.def_var(iter_var, next_idx);
                        builder
                            .ins()
                            .brif(cmp_match, loop_exit, &[], loop_header, &[]);

                        builder.switch_to_block(trap_block);
                        builder.ins().jump(trap_exit_block, &[]);

                        builder.switch_to_block(loop_exit);

                        let mut args = Vec::new();
                        for _ in 0..*arg_count {
                            args.push(pop!());
                        }
                        args.reverse();

                        let mut sig = self.module.make_signature();
                        for _ in 0..*arg_count {
                            sig.params.push(AbiParam::new(types::I64));
                        }
                        sig.returns.push(AbiParam::new(types::I64));

                        let sig_ref = builder.import_signature(sig);
                        let call = builder.ins().call_indirect(sig_ref, fn_ptr, &args);
                        let res = builder.inst_results(call)[0];
                        push!(res);
                    },

                    | Instruction::Return => {
                        let sp = builder.use_var(sp_var);
                        let cmp_under = builder.ins().icmp_imm(IntCC::SignedLessThanOrEqual, sp, 0);
                        let return_ok = builder.create_block();
                        builder
                            .ins()
                            .brif(cmp_under, trap_exit_block, &[], return_ok, &[]);
                        builder.switch_to_block(return_ok);
                        let top_off = builder.ins().iadd_imm(sp, -8);
                        let addr = builder.ins().iadd(stack_ptr, top_off);
                        let val_raw = builder.ins().load(types::I64, MemFlags::new(), addr, 0);
                        let ret_val = builder.ins().bitcast(types::F64, MemFlags::new(), val_raw);
                        builder.ins().return_(&[ret_val]);
                    },
                }
            }

            if !builder.is_unreachable() {
                let zero = builder.ins().f64const(0.0);

                builder.ins().return_(&[zero]);
            }

            builder.switch_to_block(trap_exit_block);

            let nan_val = builder.ins().f64const(f64::NAN);

            builder.ins().return_(&[nan_val]);

            builder.finalize();

            let func_id = self
                .module
                .declare_function(&name, Linkage::Export, &self.ctx.func.signature)
                .map_err(|e| e.to_string())?;

            self.module
                .define_function(func_id, &mut self.ctx)
                .map_err(|e| format!("Compilation error: {e}"))?;

            self.module.clear_context(&mut self.ctx);

            self.module
                .finalize_definitions()
                .map_err(|e| e.to_string())?;

            let code = self.module.get_finalized_function(func_id);

            Ok(code)
        }

        #[cfg(not(feature = "jit"))]
        {
            Err("JIT disabled".to_string())
        }
    }
}

#[cfg(feature = "jit")]
const fn type_to_clif(ty: JitType) -> (Type, Type) {
    match ty {
        | JitType::I8 => (types::I8, types::I64),
        | JitType::I16 => (types::I16, types::I64),
        | JitType::I32 => (types::I32, types::I64),
        | JitType::I64 => (types::I64, types::I64),
        | JitType::F32 => (types::F32, types::I64),
        | JitType::F64 => (types::F64, types::I64),
    }
}

#[cfg(feature = "jit")]
pub(crate) fn cast_to_storage(
    builder: &mut FunctionBuilder<'_>,
    val: Value,
    ty: Type,
    storage_ty: Type,
) -> Value {
    if ty == storage_ty {
        return val;
    }

    if ty.is_int() && storage_ty.is_int() && ty.bits() < storage_ty.bits() {
        return builder.ins().uextend(storage_ty, val);
    }

    if ty.is_float() && storage_ty == types::I64 {
        if ty == types::F32 {
            let bits = builder.ins().bitcast(types::I32, MemFlags::new(), val);

            return builder.ins().uextend(types::I64, bits);
        }

        if ty == types::F64 {
            return builder.ins().bitcast(types::I64, MemFlags::new(), val);
        }
    }

    val
}

#[cfg(feature = "jit")]
fn cast_from_storage(
    builder: &mut FunctionBuilder<'_>,
    val: Value,
    target_ty: Type,
) -> Value {
    if target_ty == types::I64 {
        return val;
    }

    if target_ty.is_int() {
        return builder.ins().ireduce(target_ty, val);
    }

    if target_ty == types::F32 {
        let t = builder.ins().ireduce(types::I32, val);

        return builder.ins().bitcast(types::F32, MemFlags::new(), t);
    }

    if target_ty == types::F64 {
        return builder.ins().bitcast(types::F64, MemFlags::new(), val);
    }

    val
}
