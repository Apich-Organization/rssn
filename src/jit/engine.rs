//! JIT Engine implementation using Cranelift.

use std::collections::HashMap;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;

#[cfg(feature = "jit")]
use cranelift_codegen::ir::condcodes::*;
#[cfg(feature = "jit")]
use cranelift_codegen::ir::*;
#[cfg(feature = "jit")]
use cranelift_codegen::Context;
#[cfg(feature = "jit")]
use cranelift_frontend::FunctionBuilder;
#[cfg(feature = "jit")]
use cranelift_frontend::FunctionBuilderContext;
#[cfg(feature = "jit")]
use cranelift_frontend::Variable;
#[cfg(feature = "jit")]
use cranelift_jit::JITBuilder;
#[cfg(feature = "jit")]
use cranelift_jit::JITModule;
#[cfg(feature = "jit")]
use cranelift_module::DataDescription;
#[cfg(feature = "jit")]
use cranelift_module::Linkage;
#[cfg(feature = "jit")]
use cranelift_module::Module;

use crate::jit::instructions::Instruction;
use crate::jit::instructions::JitType;

/// The JIT Engine for compiling and executing code.

pub struct JitEngine {
    #[cfg(feature = "jit")]
    builder_context:
        FunctionBuilderContext,
    #[cfg(feature = "jit")]
    ctx: Context,
    #[cfg(feature = "jit")]
    module: JITModule,

    #[cfg(feature = "jit")]
    custom_ops:
        HashMap<u32, CustomOpData>,

    function_counter: AtomicUsize,
}

#[cfg(feature = "jit")]
#[derive(Clone)]

pub struct CustomOpData {
    pub func_ptr: *const u8,
    pub arg_count: usize,
}

impl Default for JitEngine {
    fn default() -> Self {

        Self::new()
    }
}

impl JitEngine {
    pub fn new() -> Self {

        let function_counter =
            AtomicUsize::new(0);

        #[cfg(feature = "jit")]
        {

            let builder = JITBuilder::new(cranelift_module::default_libcall_names()).unwrap();

            let module =
                JITModule::new(builder);

            Self {
                builder_context: FunctionBuilderContext::new(),
                ctx: module.make_context(),
                module,
                custom_ops: HashMap::new(),
                function_counter,
            }
        }

        #[cfg(not(feature = "jit"))]
        {

            Self {
                function_counter,
            }
        }
    }

    /// Registers a custom operation handler.

    pub unsafe fn register_custom_op(
        &mut self,
        opcode: u32,
        func_ptr: *const u8,
        arg_count: usize,
    ) {

        #[cfg(feature = "jit")]
        {

            self.custom_ops
                .insert(
                    opcode,
                    CustomOpData {
                        func_ptr,
                        arg_count,
                    },
                );
        }
    }

    /// Compiles the given instructions into a function.
    /// The function takes no arguments and returns an f64.

    pub unsafe fn compile(
        &mut self,
        instructions: &[Instruction],
    ) -> Result<*const u8, String> {

        #[cfg(feature = "jit")]
        {

            let id = self
                .function_counter
                .fetch_add(
                    1,
                    Ordering::SeqCst,
                );

            let name = format!(
                "jit_fn_{}",
                id
            );

            // Setup signature: () -> f64
            self.ctx
                .func
                .signature
                .params
                .clear();

            self.ctx
                .func
                .signature
                .returns
                .clear();

            self.ctx
                .func
                .signature
                .returns
                .push(AbiParam::new(
                    types::F64,
                ));

            let mut builder = FunctionBuilder::new(&mut self.ctx.func, &mut self.builder_context);

            // Create blocks for all labels
            let mut blocks =
                HashMap::new();

            let entry_block =
                builder.create_block();

            blocks.insert(
                u32::MAX,
                entry_block,
            );

            for inst in instructions {

                if let Instruction::Label(id) = inst {
                    if !blocks.contains_key(id) {
                         blocks.insert(*id, builder.create_block());
                    }
                }
            }

            builder.append_block_params_for_function_params(entry_block);

            builder.switch_to_block(
                entry_block,
            );

            // Stack slot setup
            let stack_slot = builder.create_sized_stack_slot(StackSlotData::new(StackSlotKind::ExplicitSlot, 2048, 4));

            let stack_ptr = builder
                .ins()
                .stack_addr(
                    types::I64,
                    stack_slot,
                    0,
                );

            let sp_var = builder
                .declare_var(
                    types::I64,
                );

            let mut var_tmp = builder
                .ins()
                .iconst(types::I64, 0);

            builder.def_var(
                sp_var,
                var_tmp,
            );

            // Macros for stack manipulation
            macro_rules! pop {
                () => {{
                    let sp = builder.use_var(sp_var);
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

                    let sp = builder
                        .use_var(
                            sp_var,
                        );

                    let addr = builder
                        .ins()
                        .iadd(
                            stack_ptr,
                            sp,
                        );

                    builder
                        .ins()
                        .store(
                        MemFlags::new(),
                        val,
                        addr,
                        0,
                    );

                    let new_sp =
                        builder
                            .ins()
                            .iadd_imm(
                                sp, 8,
                            );

                    builder.def_var(
                        sp_var,
                        new_sp,
                    );
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
                    Instruction::Label(_) => {}

                    Instruction::ImmI(val) => {
                         let val_v = builder.ins().iconst(types::I64, *val);
                         push!(val_v);
                    }
                    Instruction::ImmF(val) => {
                         let val_v = builder.ins().f64const(*val);
                         push!(val_v);
                    }

                    Instruction::Load(ty) => {
                        let addr_v = pop!();
                        let (load_ty, store_ty) = type_to_clif(*ty);
                        let val = builder.ins().load(load_ty, MemFlags::new(), addr_v, 0);
                        let val_64 = cast_to_storage(&mut builder, val, load_ty, store_ty);
                        push!(val_64);
                    }

                    Instruction::Store(ty) => {
                        let val_raw = pop!();
                        let addr = pop!();
                        let (mem_ty, _) = type_to_clif(*ty);
                        let val_trunc = cast_from_storage(&mut builder, val_raw, mem_ty);
                        builder.ins().store(MemFlags::new(), val_trunc, addr, 0);
                    }

                    Instruction::Custom { opcode, payload: _ } => {
                         if let Some(op_data) = self.custom_ops.get(opcode) {
                              let op_addr = op_data.func_ptr as i64;
                              let arg_cnt = op_data.arg_count;
                              let fn_ptr_v = builder.ins().iconst(types::I64, op_addr);

                             let mut args = Vec::new();
                             for _ in 0..arg_cnt {
                                  args.push(pop!());
                             }
                             // args are popped in reverse order (top is last arg)
                             // wait, stack: [arg1, arg2], pop -> arg2, pop -> arg1
                             // functions expect (arg1, arg2). So we need to reverse the vector?
                             // No, if I push(1), push(2), stack is [1, 2]. Pop gives 2.
                             // call(1, 2). 
                             // So args are reversed.
                             // But `args` has [2, 1]. `args.reverse()` -> [1, 2].
                             // Logic was correct in previous version.
                             // But here we construct differently.
                             // With pop!, we get them in reverse order.
                             // Let's verify: 
                             // Call(x, y)
                             // Stack: ..., x, y (top)
                             // pop -> y
                             // pop -> x
                             // args = [y, x]
                             // reverse -> [x, y] -> correct.

                             let mut sig = self.module.make_signature();
                             for _ in 0..arg_cnt {
                                 sig.params.push(AbiParam::new(types::I64));
                             }
                             sig.returns.push(AbiParam::new(types::I64));

                             let sig_ref = builder.import_signature(sig);

                             // args need to be reversed?
                             // CallIndirect expects &[Value].
                             // If vector is [y, x], args are y, x.
                             // This means y is arg0. Incorrect.
                             // So we MUST reverse.
                             args.reverse();

                             let call = builder.ins().call_indirect(sig_ref, fn_ptr_v, &args);
                             let res = builder.inst_results(call)[0];
                             push!(res);

                         }
                    }

                    Instruction::Add(ty) | Instruction::Sub(ty) | Instruction::Mul(ty) | Instruction::Div(ty) => {
                         let rhs_raw = pop!();
                         let lhs_raw = pop!();

                         let (op_ty, _) = type_to_clif(*ty);
                         let rhs = cast_from_storage(&mut builder, rhs_raw, op_ty);
                         let lhs = cast_from_storage(&mut builder, lhs_raw, op_ty);

                         let res = match inst {
                             Instruction::Add(_) => if op_ty.is_int() { builder.ins().iadd(lhs, rhs) } else { builder.ins().fadd(lhs, rhs) },
                             Instruction::Sub(_) => if op_ty.is_int() { builder.ins().isub(lhs, rhs) } else { builder.ins().fsub(lhs, rhs) },
                             Instruction::Mul(_) => if op_ty.is_int() { builder.ins().imul(lhs, rhs) } else { builder.ins().fmul(lhs, rhs) },
                             Instruction::Div(_) => if op_ty.is_int() { builder.ins().sdiv(lhs, rhs) } else { builder.ins().fdiv(lhs, rhs) },
                             _ => unreachable!(),
                         };

                         let res_chk = cast_to_storage(&mut builder, res, op_ty, types::I64);
                         push!(res_chk);
                    }

                    Instruction::And | Instruction::Or | Instruction::Xor => {
                         let rhs = pop!();
                         let lhs = pop!();
                         let res = match inst {
                             Instruction::And => builder.ins().band(lhs, rhs),
                             Instruction::Or => builder.ins().bor(lhs, rhs),
                             Instruction::Xor => builder.ins().bxor(lhs, rhs),
                             _ => unreachable!(),
                         };
                         push!(res);
                    }

                    Instruction::Not => {
                         let val = pop!();
                         let res = builder.ins().bnot(val);
                         push!(res);
                    }

                    Instruction::Eq(ty) | Instruction::Ne(ty) | Instruction::Lt(ty) | Instruction::Gt(ty) | Instruction::Le(ty) | Instruction::Ge(ty) => {
                         let rhs_raw = pop!();
                         let lhs_raw = pop!();

                         let (op_ty, _) = type_to_clif(*ty);
                         let rhs = cast_from_storage(&mut builder, rhs_raw, op_ty);
                         let lhs = cast_from_storage(&mut builder, lhs_raw, op_ty);

                         let cond = if op_ty.is_int() {
                              let cc = match inst {
                                   Instruction::Eq(_) => IntCC::Equal,
                                   Instruction::Ne(_) => IntCC::NotEqual,
                                   Instruction::Lt(_) => IntCC::SignedLessThan,
                                   Instruction::Gt(_) => IntCC::SignedGreaterThan,
                                   Instruction::Le(_) => IntCC::SignedLessThanOrEqual,
                                   Instruction::Ge(_) => IntCC::SignedGreaterThanOrEqual,
                                   _ => unreachable!()
                              };
                              builder.ins().icmp(cc, lhs, rhs)
                         } else {
                              let cc = match inst {
                                   Instruction::Eq(_) => FloatCC::Equal,
                                   Instruction::Ne(_) => FloatCC::NotEqual,
                                   Instruction::Lt(_) => FloatCC::LessThan,
                                   Instruction::Gt(_) => FloatCC::GreaterThan,
                                   Instruction::Le(_) => FloatCC::LessThanOrEqual,
                                   Instruction::Ge(_) => FloatCC::GreaterThanOrEqual,
                                   _ => unreachable!()
                              };
                              builder.ins().fcmp(cc, lhs, rhs)
                         };

                         let res = builder.ins().uextend(types::I64, cond);
                         push!(res);
                    }

                    Instruction::Jump(target) => {
                        if let Some(blk) = blocks.get(target) {
                             builder.ins().jump(*blk, &[]);
                        }
                    }

                    Instruction::BranchIfTrue(target) => {
                        let val = pop!();
                        if let Some(blk) = blocks.get(target) {
                             let cond = builder.ins().icmp_imm(IntCC::NotEqual, val, 0);
                             let continue_block = builder.create_block();
                             builder.ins().brif(cond, *blk, &[], continue_block, &[]);
                             builder.switch_to_block(continue_block);
                        }
                    }

                    Instruction::BranchIfFalse(target) => {
                        let val = pop!();
                        if let Some(blk) = blocks.get(target) {
                             let cond = builder.ins().icmp_imm(IntCC::Equal, val, 0);
                             let continue_block = builder.create_block();
                             builder.ins().brif(cond, *blk, &[], continue_block, &[]);
                             builder.switch_to_block(continue_block);
                        }
                    }

                    Instruction::Dup => {
                        let val = pop!();
                        push!(val);
                        push!(val);
                    }

                    Instruction::Drop => {
                         let _ = pop!();
                    }

                    Instruction::Swap => {
                         let v1 = pop!();
                         let v2 = pop!();
                         push!(v1);
                         push!(v2);
                    }

                    Instruction::Call(arg_count) => {
                         let fn_ptr = pop!();
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
                    }

                    Instruction::Return => {
                        let sp = builder.use_var(sp_var);
                        let top_off = builder.ins().iadd_imm(sp, -8);
                        let addr = builder.ins().iadd(stack_ptr, top_off);
                        let val_raw = builder.ins().load(types::I64, MemFlags::new(), addr, 0);
                        let ret_val = builder.ins().bitcast(types::F64, MemFlags::new(), val_raw);
                        builder.ins().return_(&[ret_val]);
                    }
                }
            }

            if !builder.is_unreachable()
            {

                let zero = builder
                    .ins()
                    .f64const(0.0);

                builder
                    .ins()
                    .return_(&[zero]);
            }

            builder.finalize();

            let func_id = self
                .module
                .declare_function(
                    &name,
                    Linkage::Export,
                    &self
                        .ctx
                        .func
                        .signature,
                )
                .map_err(|e| {

                    e.to_string()
                })?;

            self.module
                .define_function(
                    func_id,
                    &mut self.ctx,
                )
                .map_err(|e| {

                    e.to_string()
                })?;

            self.module
                .clear_context(
                    &mut self.ctx,
                );

            self.module
                .finalize_definitions()
                .map_err(|e| {

                    e.to_string()
                })?;

            let code = self
                .module
                .get_finalized_function(
                    func_id,
                );

            Ok(code)
        }

        #[cfg(not(feature = "jit"))]
        {

            Err("JIT disabled"
                .to_string())
        }
    }
}

#[cfg(feature = "jit")]

fn type_to_clif(
    ty: JitType
) -> (Type, Type) {

    match ty {
        | JitType::I8 => {
            (
                types::I8,
                types::I64,
            )
        },
        | JitType::I16 => {
            (
                types::I16,
                types::I64,
            )
        },
        | JitType::I32 => {
            (
                types::I32,
                types::I64,
            )
        },
        | JitType::I64 => {
            (
                types::I64,
                types::I64,
            )
        },
        | JitType::F32 => {
            (
                types::F32,
                types::I64,
            )
        },
        | JitType::F64 => {
            (
                types::F64,
                types::I64,
            )
        },
    }
}

#[cfg(feature = "jit")]

fn cast_to_storage(
    builder: &mut FunctionBuilder<'_>,
    val: Value,
    ty: Type,
    storage_ty: Type,
) -> Value {

    if ty == storage_ty {

        return val;
    }

    if ty.is_int()
        && storage_ty.is_int()
    {

        if ty.bits() < storage_ty.bits()
        {

            return builder
                .ins()
                .uextend(
                    storage_ty,
                    val,
                );
        }
    }

    if ty.is_float()
        && storage_ty == types::I64
    {

        if ty == types::F32 {

            let bits = builder
                .ins()
                .bitcast(
                    types::I32,
                    MemFlags::new(),
                    val,
                );

            return builder
                .ins()
                .uextend(
                    types::I64,
                    bits,
                );
        }

        if ty == types::F64 {

            return builder
                .ins()
                .bitcast(
                    types::I64,
                    MemFlags::new(),
                    val,
                );
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

        return builder
            .ins()
            .ireduce(target_ty, val);
    }

    if target_ty == types::F32 {

        let t = builder
            .ins()
            .ireduce(types::I32, val);

        return builder
            .ins()
            .bitcast(
                types::F32,
                MemFlags::new(),
                t,
            );
    }

    if target_ty == types::F64 {

        return builder
            .ins()
            .bitcast(
                types::F64,
                MemFlags::new(),
                val,
            );
    }

    val
}
