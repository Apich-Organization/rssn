//! JIT Engine implementation using Cranelift.

use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::jit::instructions::{Instruction, JitType};

#[cfg(feature = "jit")]
use cranelift::prelude::*;
#[cfg(feature = "jit")]
use cranelift::jit::{JITBuilder, JITModule};
#[cfg(feature = "jit")]
use cranelift::module::{DataContext, Linkage, Module};

/// The JIT Engine for compiling and executing code.
pub struct JitEngine {
    #[cfg(feature = "jit")]
    builder_context: FunctionBuilderContext,
    #[cfg(feature = "jit")]
    ctx: codegen::Context,
    #[cfg(feature = "jit")]
    module: JITModule,
    #[cfg(feature = "jit")]
    data_ctx: DataContext,
    
    #[cfg(feature = "jit")]
    custom_ops: HashMap<u32, CustomOpData>,

    function_counter: AtomicUsize,
}

#[cfg(feature = "jit")]
#[derive(Clone)]
pub struct CustomOpData {
    pub func_ptr: *const u8,
    pub arg_count: usize,
    // For simplicity, we assume 1 return value pushed to stack if not void, 
    // but JIT stack logic (expressions) generally produces 1 value. 
    // Let's assume custom ops behave like functions: (args...) -> res.
    // If res is void, push 0?
    // Let's enforce 1 return value (I64) for custom ops for now to fit the stack model easily.
}


impl Default for JitEngine {
    fn default() -> Self {
        Self::new()
    }
}

impl JitEngine {
    pub fn new() -> Self {
        let function_counter = AtomicUsize::new(0);
        
        #[cfg(feature = "jit")]
        {
            let builder = JITBuilder::new(cranelift::module::default_libcall_names()).unwrap();
            let module = JITModule::new(builder);
            Self {
                builder_context: FunctionBuilderContext::new(),
                ctx: module.make_context(),
                ctx: module.make_context(),
                module,
                data_ctx: DataContext::new(),
                custom_ops: HashMap::new(),
                function_counter,
            }
        }
        #[cfg(not(feature = "jit"))]
        {
            Self { function_counter }
        }
    }

    /// Registers a custom operation handler.
    /// When `Instruction::Custom { opcode, .. }` is encountered, the JIT will compile a call
    /// to `func_ptr`.
    /// 
    /// # Safety
    /// `func_ptr` must be valid fn(arg1..argN) -> i64.
    pub unsafe fn register_custom_op(&mut self, opcode: u32, func_ptr: *const u8, arg_count: usize) {
        #[cfg(feature = "jit")]
        {
            self.custom_ops.insert(opcode, CustomOpData { func_ptr, arg_count });
        }
    }

    /// Compiles the given instructions into a function.
    /// The function takes no arguments and returns an f64.
    ///
    /// # Safety
    /// Only safe if the instructions are well-formed and memory accesses are valid.
    pub unsafe fn compile(&mut self, instructions: &[Instruction]) -> Result<*const u8, String> {
        #[cfg(feature = "jit")]
        {
            let id = self.function_counter.fetch_add(1, Ordering::SeqCst);
            let name = format!("jit_fn_{}", id);
            
            // Setup signature: () -> f64
            self.ctx.func.signature.params.clear();
            self.ctx.func.signature.returns.clear();
            self.ctx.func.signature.returns.push(AbiParam::new(types::F64));
            
            let mut builder = FunctionBuilder::new(&mut self.ctx.func, &mut self.builder_context);
            
            // Create blocks for all labels
            let mut blocks = HashMap::new();
            let entry_block = builder.create_block();
            blocks.insert(u32::MAX, entry_block); // Special key for entry? No, just start there.
            
            for inst in instructions {
                if let Instruction::Label(id) = inst {
                    if !blocks.contains_key(id) {
                         blocks.insert(*id, builder.create_block());
                    }
                }
            }
            
            builder.append_block_params_for_function_params(entry_block);
            builder.switch_to_block(entry_block);
            
            // We need to manage stack across blocks? 
            // Standard JIT stack approach: Cranelift Values are SSA, not stack slots.
            // If we jump between blocks, we must pass stack values as block arguments OR use a separate in-memory stack.
            // For a simple implementation to support loops, using a real memory stack (alloca) is easiest
            // because mapping stack depth and types across arbitrary jumps in SSA is hard (Phi nodes).
            // Let's use a "shadow stack" allocated on the stack frame.
            
            // Allocate a stack slot for our operand stack.
            // Let's say max stack depth 128.
            let stack_slot = builder.create_stack_slot(StackSlotData::new(StackSlotKind::ExplicitSlot, 1024)); // 128 * 8 bytes
            let stack_ptr = builder.ins().stack_addr(types::I64, stack_slot, 0);
            
            // We need a mutable stack index variable (phi/reg?).
            // For simplicity in this logic-grammar JIT, let's keep `sp` (stack pointer offset) in a Variable?
            // Cranelift Variables are easy for this.
            
            let sp_var = Variable::new(0);
            builder.declare_var(sp_var, types::I64);
            builder.def_var(sp_var, builder.ins().iconst(types::I64, 0));
            
            // Helper to push/pop from memory stack
            // We assume all items are 64-bit for simplicity on the stack storage, even if I8.
            // (Promote everything to 64-bit on stack).
            
            // Note: closures below borrow builder mutably. We can't definte them easily if they capture builder.
            // We'll write logic inline or use macros.
            
            for inst in instructions {
                // If this instruction marks start of a block (Label), switch to it.
                // Note: Fallthrough must be handled. 
                if let Instruction::Label(id) = inst {
                     let block = blocks[id];
                     // If current block not terminated, jump to new block
                     if !builder.is_filled() {
                         builder.ins().jump(block, &[]);
                     }
                     builder.switch_to_block(block);
                     // Note: Variables (sp_var) work across blocks automatically.
                }

                match inst {
                    Instruction::Label(_) => {} // handled above
                    
                    Instruction::ImmI(val) => {
                         let val_v = builder.ins().iconst(types::I64, *val);
                         
                         let sp = builder.use_var(sp_var);
                         builder.ins().store(MemFlags::new(), val_v, stack_ptr, sp);
                         let new_sp = builder.ins().iadd_imm(sp, 8);
                         builder.def_var(sp_var, new_sp);
                    }
                    Instruction::ImmF(val) => {
                         let val_v = builder.ins().f64const(*val);
                         // Store as 64 bits. F64 fits.
                         let sp = builder.use_var(sp_var);
                         builder.ins().store(MemFlags::new(), val_v, stack_ptr, sp);
                         let new_sp = builder.ins().iadd_imm(sp, 8);
                         builder.def_var(sp_var, new_sp);
                    }
                    
                    Instruction::Load(ty) => {
                        // pop addr
                        let sp = builder.use_var(sp_var);
                        let old_sp = builder.ins().isub_imm(sp, 8);
                        let addr_v = builder.ins().load(types::I64, MemFlags::new(), stack_ptr, old_sp);
                        builder.def_var(sp_var, old_sp); // effectively popped
                        
                        let (load_ty, store_ty) = type_to_clif(*ty);
                        let val = builder.ins().load(load_ty, MemFlags::new(), addr_v, 0);
                        
                        // Extend to 64 bit for stack storage if needed
                        let val_64 = cast_to_storage(&mut builder, val, load_ty, store_ty);
                        
                        // Push val
                        let sp = builder.use_var(sp_var);
                        builder.ins().store(MemFlags::new(), val_64, stack_ptr, sp);
                        let new_sp = builder.ins().iadd_imm(sp, 8);
                        builder.def_var(sp_var, new_sp);
                    }
                    
                    Instruction::Store(ty) => {
                        // pop val, pop addr
                        let sp = builder.use_var(sp_var);
                        let sp_1 = builder.ins().isub_imm(sp, 8);
                        let val_raw = builder.ins().load(types::I64, MemFlags::new(), stack_ptr, sp_1);
                        
                        let sp_2 = builder.ins().isub_imm(sp_1, 8);
                        let addr = builder.ins().load(types::I64, MemFlags::new(), stack_ptr, sp_2);
                        builder.def_var(sp_var, sp_2);
                        
                        let (mem_ty, _) = type_to_clif(*ty);
                        // Truncate val_raw to mem_ty
                        let val_trunc = cast_from_storage(&mut builder, val_raw, mem_ty);
                        
                        builder.ins().store(MemFlags::new(), val_trunc, addr, 0);
                    }

                    Instruction::Custom { opcode, payload } => {
                         // Find handler
                         if let Some(op_data) = self.custom_ops.get(opcode) {
                              let op_addr = op_data.func_ptr as i64;
                              let arg_cnt = op_data.arg_count;
                              let fn_ptr_v = builder.ins().iconst(types::I64, op_addr);
                              
                              // We currently don't use payload in the call arguments, 
                              // but we could pass it as a pointer if we allocated it?
                              // For simplicity: Custom instruction = Call registered func with stack args.
                              // Payload is ignored in this simple binding (or user can use ImmI).
                              // Using payload would require data section management.
                              
                              // Pop `arg_cnt`
                             let mut args = Vec::new();
                             let sp = builder.use_var(sp_var);
                             let mut current_sp = sp;
                             
                             for _ in 0..arg_cnt {
                                  current_sp = builder.ins().isub_imm(current_sp, 8);
                                  let arg_val = builder.ins().load(types::I64, MemFlags::new(), stack_ptr, current_sp);
                                  // Assume all args I64
                                  args.push(arg_val);
                             }
                             builder.def_var(sp_var, current_sp);
                             args.reverse();

                             let mut sig = self.module.make_signature();
                             for _ in 0..arg_cnt {
                                 sig.params.push(AbiParam::new(types::I64));
                             }
                             // Assume always returns 1 value I64
                             sig.returns.push(AbiParam::new(types::I64));
                             
                             let sig_ref = builder.import_signature(sig);
                             let call = builder.ins().call_indirect(sig_ref, fn_ptr_v, &args);
                             let res = builder.inst_results(call)[0];
                             
                             // Push result
                             let sp = builder.use_var(sp_var);
                             builder.ins().store(MemFlags::new(), res, stack_ptr, sp);
                             let new_sp = builder.ins().iadd_imm(sp, 8);
                             builder.def_var(sp_var, new_sp);

                         } else {
                             // Opcode not found. Panic? Or No-op.
                             // No-op for safety, but this is a runtime/compile time configuration error.
                         }
                    }

                    // Logic/Math: Pop 2, Op, Push 1
                    // For brevity, we implement a helper logic inline
                    Instruction::Add(ty) | Instruction::Sub(ty) | Instruction::Mul(ty) | Instruction::Div(ty) => {
                         // Pop rhs, lhs
                         let sp = builder.use_var(sp_var);
                         let sp_1 = builder.ins().isub_imm(sp, 8);
                         let rhs_raw = builder.ins().load(types::I64, MemFlags::new(), stack_ptr, sp_1);
                         let sp_2 = builder.ins().isub_imm(sp_1, 8);
                         let lhs_raw = builder.ins().load(types::I64, MemFlags::new(), stack_ptr, sp_2);
                         builder.def_var(sp_var, sp_2);
                         
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
                         
                         let res_chk = cast_to_storage(&mut builder, res, op_ty, types::I64); // or F64
                         
                         // Push
                         let sp = builder.use_var(sp_var);
                         builder.ins().store(MemFlags::new(), res_chk, stack_ptr, sp);
                         let new_sp = builder.ins().iadd_imm(sp, 8);
                         builder.def_var(sp_var, new_sp);
                    }
                    
                    // Bitwise
                    Instruction::And | Instruction::Or | Instruction::Xor => {
                         let sp = builder.use_var(sp_var);
                         let sp_1 = builder.ins().isub_imm(sp, 8);
                         let rhs = builder.ins().load(types::I64, MemFlags::new(), stack_ptr, sp_1);
                         let sp_2 = builder.ins().isub_imm(sp_1, 8);
                         let lhs = builder.ins().load(types::I64, MemFlags::new(), stack_ptr, sp_2);
                         builder.def_var(sp_var, sp_2);
                         
                         let res = match inst {
                             Instruction::And => builder.ins().band(lhs, rhs),
                             Instruction::Or => builder.ins().bor(lhs, rhs),
                             Instruction::Xor => builder.ins().bxor(lhs, rhs),
                             _ => unreachable!(),
                         };
                         
                         let sp = builder.use_var(sp_var);
                         builder.ins().store(MemFlags::new(), res, stack_ptr, sp);
                         let new_sp = builder.ins().iadd_imm(sp, 8);
                         builder.def_var(sp_var, new_sp);
                    }
                    
                    Instruction::Not => {
                         let sp = builder.use_var(sp_var);
                         let sp_1 = builder.ins().isub_imm(sp, 8);
                         let val = builder.ins().load(types::I64, MemFlags::new(), stack_ptr, sp_1);
                         
                         let res = builder.ins().bnot(val);
                         
                         builder.ins().store(MemFlags::new(), res, stack_ptr, sp_1);
                         // sp unchanged
                    }

                    // Comparisons
                    Instruction::Eq(ty) | Instruction::Ne(ty) | Instruction::Lt(ty) | Instruction::Gt(ty) | Instruction::Le(ty) | Instruction::Ge(ty) => {
                         let sp = builder.use_var(sp_var);
                         let sp_1 = builder.ins().isub_imm(sp, 8);
                         let rhs_raw = builder.ins().load(types::I64, MemFlags::new(), stack_ptr, sp_1);
                         let sp_2 = builder.ins().isub_imm(sp_1, 8);
                         let lhs_raw = builder.ins().load(types::I64, MemFlags::new(), stack_ptr, sp_2);
                         builder.def_var(sp_var, sp_2);
                         
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
                         
                         // cond is I8 (boolean) usually 1 or 0. Extend to I64.
                         let res = builder.ins().uextend(types::I64, cond);
                         
                         let sp = builder.use_var(sp_var);
                         builder.ins().store(MemFlags::new(), res, stack_ptr, sp);
                         let new_sp = builder.ins().iadd_imm(sp, 8);
                         builder.def_var(sp_var, new_sp);
                    }

                    Instruction::Jump(target) => {
                        if let Some(blk) = blocks.get(target) {
                             builder.ins().jump(*blk, &[]);
                        }
                    }
                    
                    Instruction::BranchIfTrue(target) => {
                        // pop cond
                        let sp = builder.use_var(sp_var);
                        let sp_1 = builder.ins().isub_imm(sp, 8);
                        let val = builder.ins().load(types::I64, MemFlags::new(), stack_ptr, sp_1);
                        builder.def_var(sp_var, sp_1);
                        
                        if let Some(blk) = blocks.get(target) {
                             // brnz expects integer condition? val is I64.
                             builder.ins().brnz(val, *blk, &[]);
                        }
                        // Fallthrough continues. 
                    }
                    
                     Instruction::BranchIfFalse(target) => {
                        // pop cond
                        let sp = builder.use_var(sp_var);
                        let sp_1 = builder.ins().isub_imm(sp, 8);
                        let val = builder.ins().load(types::I64, MemFlags::new(), stack_ptr, sp_1);
                        builder.def_var(sp_var, sp_1);
                        
                        if let Some(blk) = blocks.get(target) {
                             builder.ins().brz(val, *blk, &[]);
                        }
                    }
                    
                    Instruction::Dup => {
                        let sp = builder.use_var(sp_var);
                        let top = builder.ins().isub_imm(sp, 8);
                        let val = builder.ins().load(types::I64, MemFlags::new(), stack_ptr, top);
                        
                        builder.ins().store(MemFlags::new(), val, stack_ptr, sp);
                        let new_sp = builder.ins().iadd_imm(sp, 8);
                        builder.def_var(sp_var, new_sp);
                    }
                    
                    Instruction::Drop => {
                         let sp = builder.use_var(sp_var);
                         let new_sp = builder.ins().isub_imm(sp, 8);
                         builder.def_var(sp_var, new_sp);
                    }
                    
                    Instruction::Swap => {
                         let sp = builder.use_var(sp_var);
                         let sp_1 = builder.ins().isub_imm(sp, 8);
                         let sp_2 = builder.ins().isub_imm(sp_1, 8);
                         
                         let v1 = builder.ins().load(types::I64, MemFlags::new(), stack_ptr, sp_1);
                         let v2 = builder.ins().load(types::I64, MemFlags::new(), stack_ptr, sp_2);
                         
                         builder.ins().store(MemFlags::new(), v1, stack_ptr, sp_2);
                         builder.ins().store(MemFlags::new(), v2, stack_ptr, sp_1);
                    }

                    Instruction::Call(arg_count) => {
                         // Pop ptr
                         let sp = builder.use_var(sp_var);
                         let sp_ptr = builder.ins().isub_imm(sp, 8);
                         let fn_ptr = builder.ins().load(types::I64, MemFlags::new(), stack_ptr, sp_ptr);
                         
                         // Prepare args.
                         // Pop arg_count items into a vector
                         let mut args = Vec::new();
                         let mut current_sp = sp_ptr;
                         for _ in 0..*arg_count {
                              current_sp = builder.ins().isub_imm(current_sp, 8);
                              let arg_val = builder.ins().load(types::I64, MemFlags::new(), stack_ptr, current_sp);
                              // Assume args are I64 ? Or F64?
                              // A generic call is hard without signature. 
                              // User request is "direct access". 
                              // Let's assume all args are I64 (pointers/ints) or F64 reinterpreted.
                              // We need a specific signature.
                              // For now, let's create a signature: (I64, I64, ...) -> I64.
                              args.push(arg_val);
                         }
                         builder.def_var(sp_var, current_sp);
                         args.reverse(); // Pop order is reverse
                         
                         // Make signature
                         let mut sig = self.module.make_signature();
                         for _ in 0..*arg_count {
                             sig.params.push(AbiParam::new(types::I64));
                         }
                         sig.returns.push(AbiParam::new(types::I64));
                         
                         let sig_ref = builder.import_signature(sig);
                         let call = builder.ins().call_indirect(sig_ref, fn_ptr, &args);
                         let res = builder.inst_results(call)[0];
                         
                         // Push res
                         let sp = builder.use_var(sp_var);
                         builder.ins().store(MemFlags::new(), res, stack_ptr, sp);
                         let new_sp = builder.ins().iadd_imm(sp, 8);
                         builder.def_var(sp_var, new_sp);
                    }

                    Instruction::Return => {
                        let sp = builder.use_var(sp_var);
                        // Check if empty (sp==0)
                        let is_empty = builder.ins().icmp_imm(IntCC::Equal, sp, 0);
                        // We must return F64.
                        // Since basic blocks can't conditional return easily without splitting,
                        // assuming well-formed code (stack has 1 item).
                        
                        let top_off = builder.ins().isub_imm(sp, 8);
                        let val_raw = builder.ins().load(types::I64, MemFlags::new(), stack_ptr, top_off);
                        // Reinterpret as F64 for return
                        let ret_val = builder.ins().bitcast(types::F64, val_raw);
                        
                        builder.ins().return_(&[ret_val]);
                    }
                }
            }
            
             // Fallthrough finalizer
            if !builder.is_filled() {
                 let zero = builder.ins().f64const(0.0);
                 builder.ins().return_(&[zero]);
            }
            
            builder.finalize();
            
             // Declare and define
            let func_id = self.module.declare_function(&name, Linkage::Export, &self.ctx.func.signature)
                .map_err(|e| e.to_string())?;
            
            self.module.define_function(func_id, &mut self.ctx)
                .map_err(|e| e.to_string())?;
            
            self.module.clear_context(&mut self.ctx);
            
            self.module.finalize_definitions()
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
fn type_to_clif(ty: JitType) -> (Type, Type) {
    match ty {
        JitType::I8 => (types::I8, types::I64),
        JitType::I16 => (types::I16, types::I64),
        JitType::I32 => (types::I32, types::I64),
        JitType::I64 => (types::I64, types::I64),
        JitType::F32 => (types::F32, types::I64), // Stored as I64 (f32 extended or bitcast?) F32 usually bitcast to I32 then ext? Or store F32 in 64 slot?
        JitType::F64 => (types::F64, types::I64),
    }
}

#[cfg(feature = "jit")]
fn cast_to_storage(builder: &mut FunctionBuilder, val: Value, ty: Type, storage_ty: Type) -> Value {
    if ty == storage_ty {
        return val;
    }
    if ty.is_int() && storage_ty.is_int() {
        if ty.bits() < storage_ty.bits() {
             return builder.ins().uextend(storage_ty, val);
        }
    }
    if ty.is_float() && storage_ty == types::I64 {
        if ty == types::F32 {
             let bits = builder.ins().bitcast(types::I32, val);
             return builder.ins().uextend(types::I64, bits);
        }
        if ty == types::F64 {
             return builder.ins().bitcast(types::I64, val);
        }
    }
    // Default
    val
}

#[cfg(feature = "jit")]
fn cast_from_storage(builder: &mut FunctionBuilder, val: Value, target_ty: Type) -> Value {
    // val is I64
    if target_ty == types::I64 {
        return val;
    }
    if target_ty.is_int() {
        return builder.ins().ireduce(target_ty, val);
    }
    if target_ty == types::F32 {
        let t = builder.ins().ireduce(types::I32, val);
        return builder.ins().bitcast(types::F32, t);
    }
    if target_ty == types::F64 {
        return builder.ins().bitcast(types::F64, val);
    }
    val
}
