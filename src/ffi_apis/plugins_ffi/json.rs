//! JSON-based FFI API for the Plugin Manager.

use std::ffi::CStr;
use std::ffi::CString;
use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::ffi_api::FfiResult;
use crate::plugins::manager::GLOBAL_PLUGIN_MANAGER;
use crate::symbolic::core::Expr;

#[derive(Deserialize)]

struct PluginExecutionRequest {
    plugin_name: String,
    command: String,
    args: Expr,
}

/// Executes a plugin command via JSON (args passed as JSON expr).
#[no_mangle]

pub unsafe extern "C" fn rssn_plugins_execute_json(
    json_ptr: *const c_char
) -> *mut c_char {

    if json_ptr.is_null() {

        return std::ptr::null_mut();
    }

    let json_str = match CStr::from_ptr(
        json_ptr,
    )
    .to_str()
    {
        | Ok(s) => s,
        | Err(_) => {
            return std::ptr::null_mut()
        },
    };

    let req: PluginExecutionRequest =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(r) => r,
            | Err(e) => {

                let res: FfiResult<
                    Expr,
                    String,
                > = FfiResult {
                    ok: None,
                    err: Some(
                        e.to_string(),
                    ),
                };

                return CString::new(serde_json::to_string(&res).unwrap()).unwrap().into_raw();
            },
        };

    let manager =
        match GLOBAL_PLUGIN_MANAGER
            .read()
        {
            | Ok(m) => m,
            | Err(_) => {

                let res: FfiResult<
                    Expr,
                    String,
                > = FfiResult {
                    ok: None,
                    err: Some(
                        "Lock poison"
                            .to_string(
                            ),
                    ),
                };

                return CString::new(serde_json::to_string(&res).unwrap()).unwrap().into_raw();
            },
        };

    match manager.execute_plugin(
        &req.plugin_name,
        &req.command,
        &req.args,
    ) {
        | Ok(result_expr) => {

            let res: FfiResult<
                Expr,
                String,
            > = FfiResult {
                ok: Some(result_expr),
                err: None,
            };

            CString::new(
                serde_json::to_string(
                    &res,
                )
                .unwrap(),
            )
            .unwrap()
            .into_raw()
        },
        | Err(e) => {

            let res: FfiResult<
                Expr,
                String,
            > = FfiResult {
                ok: None,
                err: Some(
                    e.to_string(),
                ),
            };

            CString::new(
                serde_json::to_string(
                    &res,
                )
                .unwrap(),
            )
            .unwrap()
            .into_raw()
        },
    }
}
