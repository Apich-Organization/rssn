//! Bincode-based FFI API for the Plugin Manager.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::plugins::manager::GLOBAL_PLUGIN_MANAGER;
use crate::symbolic::core::Expr;

#[derive(Deserialize)]

struct PluginExecutionRequest {
    plugin_name: String,
    command: String,
    args: Expr,
}

fn decode<
    T: for<'de> Deserialize<'de>,
>(
    data: *const u8,
    len: usize,
) -> Option<T> {

    if data.is_null() {

        return None;
    }

    let slice = unsafe {

        std::slice::from_raw_parts(
            data, len,
        )
    };

    bincode_next::serde::decode_from_slice(slice, bincode_next::config::standard())
        .ok()
        .map(|(v, _)| v)
}

fn encode<T: Serialize>(
    val: &T
) -> BincodeBuffer {

    match bincode_next::serde::encode_to_vec(val, bincode_next::config::standard()) {
        Ok(bytes) => BincodeBuffer::from_vec(bytes),
        Err(_) => BincodeBuffer::empty(),
    }
}

/// Executes a plugin command via Bincode.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_plugins_execute_bincode(
    data: *const u8,
    len: usize,
) -> BincodeBuffer {

    let req: PluginExecutionRequest =
        match decode(data, len) {
            | Some(r) => r,
            | None => return encode(
                &FfiResult::<
                    Expr,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Decode error"
                            .to_string(
                            ),
                    ),
                },
            ),
        };

    let manager =
        match GLOBAL_PLUGIN_MANAGER
            .read()
        {
            | Ok(m) => m,
            | Err(_) => return encode(
                &FfiResult::<
                    Expr,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Lock poison"
                            .to_string(
                            ),
                    ),
                },
            ),
        };

    match manager.execute_plugin(
        &req.plugin_name,
        &req.command,
        &req.args,
    ) {
        | Ok(result_expr) => {
            encode(&FfiResult::<
                Expr,
                String,
            > {
                ok: Some(result_expr),
                err: None,
            })
        },
        | Err(e) => {
            encode(&FfiResult::<
                Expr,
                String,
            > {
                ok: None,
                err: Some(
                    e.to_string(),
                ),
            })
        },
    }
}
