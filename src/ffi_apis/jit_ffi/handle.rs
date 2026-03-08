//! Handle APIs for JIT Engine.

use crate::jit::JitEngine;

/// Creates a new JIT Engine instance.
#[unsafe(no_mangle)]
pub extern "C" fn rssn_jit_create() -> *mut JitEngine {
    Box::into_raw(Box::new(JitEngine::new()))
}

/// Frees a JIT Engine instance.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_jit_free(engine: *mut JitEngine) {
    unsafe {
        if !engine.is_null() {
            let _ = Box::from_raw(engine);
        }
    }
}

/// Executes a JIT-compiled function pointer.
///
/// # Safety
/// The function pointer must be a valid pointer returned by `rssn_jit_compile_*`.
/// It assumes the function signature is `fn() -> f64`.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_jit_execute(
    engine: *mut JitEngine,
    func_ptr: *const u8,
) -> f64 {
    unsafe {
        if engine.is_null() || func_ptr.is_null() {
            return 0.0;
        }

        let sandbox_ctx = (*engine).build_sandbox_context();

        /// # Safety
        ///
        /// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
        /// The caller must ensure:
        /// 1. All pointer arguments are valid and point to initialized memory.
        /// 2. The memory layout of passed structures matches the expected C-ABI layout.
        /// 3. Any pointers returned by this function are managed according to the API's ownership rules.
        let func: unsafe extern "C" fn(*const crate::jit::engine::SandboxContext) -> f64 =
            std::mem::transmute(func_ptr);

        func(&raw const sandbox_ctx)
    }
}

/// Registers a custom instruction handler.
///
/// `opcode`: The exact u32 opcode found in `Instruction::Custom`.
/// `func_ptr`: Pointer to the C function to call. Signature must be `fn(i64, ...) -> i64` where `i64` represents a stack value.
/// `arg_count`: Number of arguments the function expects (popped from stack).
///
/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_jit_register_custom_op(
    engine: *mut JitEngine,
    opcode: u32,
    func_ptr: *const u8,
    arg_count: usize,
) {
    unsafe {
        if engine.is_null() || func_ptr.is_null() {
            return;
        }

        let engine = &mut *engine;

        engine.register_custom_op(opcode, func_ptr, arg_count);
    }
}

/// Registers a memory region to the sandbox, returning the region ID.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_jit_register_memory_region(
    engine: *mut JitEngine,
    base: *mut u8,
    size: usize,
) -> u16 {
    unsafe {
        if engine.is_null() {
            return 0;
        }

        (*engine).register_memory_region(base, size)
    }
}

/// Clears all registered memory regions.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_jit_clear_memory_regions(engine: *mut JitEngine) {
    unsafe {
        if engine.is_null() {
            return;
        }

        (*engine).clear_memory_regions();
    }
}

/// Appends a function pointer to the list of allowed JIT calls.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_jit_allow_call_target(
    engine: *mut JitEngine,
    func_ptr: *const u8,
) {
    unsafe {
        if engine.is_null() || func_ptr.is_null() {
            return;
        }

        (*engine).allow_call_target(func_ptr);
    }
}
