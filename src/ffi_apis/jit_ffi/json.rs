//! JSON API for JIT compilation.

use std::ffi::CStr;
use std::ffi::CString;
use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::ffi_api::FfiResult;
use crate::jit::Instruction;
use crate::jit::JitEngine;

#[derive(Deserialize)]
struct JitCompileRequest {
    instructions: Vec<Instruction>,
}

#[derive(Deserialize)]
struct JitSandboxConfig {
    memory_regions: Option<Vec<MemoryRegionConfig>>,
    allowed_calls: Option<Vec<usize>>,
}

#[derive(Deserialize)]
struct MemoryRegionConfig {
    base: usize,
    size: usize,
}

/// Configures the sandbox via JSON.
/// Returns a JSON result containing a boolean indicating success.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_jit_configure_sandbox_json(
    engine: *mut JitEngine,
    json_ptr: *const c_char,
) -> *mut c_char {
    unsafe {
        if engine.is_null() || json_ptr.is_null() {
            return std::ptr::null_mut();
        }

        let engine = &mut *engine;

        let json_str = match CStr::from_ptr(json_ptr).to_str() {
            | Ok(s) => s,
            | Err(_) => return std::ptr::null_mut(),
        };

        let req: JitSandboxConfig = match serde_json::from_str(json_str) {
            | Ok(r) => r,
            | Err(e) => {
                let res: FfiResult<bool, String> = FfiResult {
                    ok: None,
                    err: Some(e.to_string()),
                };

                return CString::new(serde_json::to_string(&res).unwrap())
                    .unwrap()
                    .into_raw();
            },
        };

        engine.clear_memory_regions();

        if let Some(regions) = req.memory_regions {
            for r in regions {
                let _ = engine.register_memory_region(r.base as *mut u8, r.size);
            }
        }

        if let Some(calls) = req.allowed_calls {
            for c in calls {
                engine.allow_call_target(c as *const u8);
            }
        }

        let res: FfiResult<bool, String> = FfiResult {
            ok: Some(true),
            err: None,
        };

        CString::new(serde_json::to_string(&res).unwrap())
            .unwrap()
            .into_raw()
    }
}

/// Compiles a sequence of instructions provided as JSON.
/// Returns a JSON result containing the address (as usize) of the compiled function.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.
///
/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_jit_compile_json(
    engine: *mut JitEngine,
    json_ptr: *const c_char,
) -> *mut c_char {
    unsafe {
        if engine.is_null() || json_ptr.is_null() {
            return std::ptr::null_mut();
        }

        let engine = &mut *engine;

        let json_str = match CStr::from_ptr(json_ptr).to_str() {
            | Ok(s) => s,
            | Err(_) => return std::ptr::null_mut(),
        };

        let req: JitCompileRequest = match serde_json::from_str(json_str) {
            | Ok(r) => r,
            | Err(e) => {
                let res: FfiResult<usize, String> = FfiResult {
                    ok: None,
                    err: Some(e.to_string()),
                };

                return CString::new(serde_json::to_string(&res).unwrap())
                    .unwrap()
                    .into_raw();
            },
        };

        match engine.compile(&req.instructions) {
            | Ok(ptr) => {
                let res: FfiResult<usize, String> = FfiResult {
                    ok: Some(ptr as usize),
                    err: None,
                };

                CString::new(serde_json::to_string(&res).unwrap())
                    .unwrap()
                    .into_raw()
            },
            | Err(e) => {
                let res: FfiResult<usize, String> = FfiResult {
                    ok: None,
                    err: Some(e),
                };

                CString::new(serde_json::to_string(&res).unwrap())
                    .unwrap()
                    .into_raw()
            },
        }
    }
}
