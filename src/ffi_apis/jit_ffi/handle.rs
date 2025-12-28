//! Handle APIs for JIT Engine.

use crate::jit::JitEngine;

/// Creates a new JIT Engine instance.
#[no_mangle]

pub extern "C" fn rssn_jit_create(
) -> *mut JitEngine {

    Box::into_raw(Box::new(
        JitEngine::new(),
    ))
}

/// Frees a JIT Engine instance.
#[no_mangle]

pub unsafe extern "C" fn rssn_jit_free(
    engine: *mut JitEngine
) {

    if !engine.is_null() {

        let _ = Box::from_raw(engine);
    }
}

/// Executes a JIT-compiled function pointer.
///
/// # Safety
/// The function pointer must be a valid pointer returned by `rssn_jit_compile_*`.
/// It assumes the function signature is `fn() -> f64`.
#[no_mangle]

pub unsafe extern "C" fn rssn_jit_execute(
    func_ptr: *const u8
) -> f64 {

    if func_ptr.is_null() {

        return 0.0;
    }

    let func: unsafe extern "C" fn() -> f64 = std::mem::transmute(func_ptr);

    func()
}

/// Registers a custom instruction handler.
///
/// `opcode`: The exact u32 opcode found in `Instruction::Custom`.
/// `func_ptr`: Pointer to the C function to call. Signature must be `fn(i64, ...) -> i64` where `i64` represents a stack value.
/// `arg_count`: Number of arguments the function expects (popped from stack).
#[no_mangle]

pub unsafe extern "C" fn rssn_jit_register_custom_op(
    engine: *mut JitEngine,
    opcode: u32,
    func_ptr: *const u8,
    arg_count: usize,
) {

    if engine.is_null()
        || func_ptr.is_null()
    {

        return;
    }

    let engine = &mut *engine;

    engine.register_custom_op(
        opcode,
        func_ptr,
        arg_count,
    );
}
