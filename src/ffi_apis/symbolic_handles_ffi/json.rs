//! JSON-based FFI API for the HandleManager.

use std::os::raw::c_char;

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::handles::HANDLE_MANAGER;

#[derive(Serialize, Deserialize)]

struct HandleInfo {
    handle : usize,
    expression : String,
}

#[derive(Serialize, Deserialize)]
#[allow(dead_code)]

struct HandleListResponse {
    handles : Vec<usize>,
    count : usize,
}

#[derive(Serialize, Deserialize)]

struct HandleStatsResponse {
    count : usize,
    handles : Vec<HandleInfo>,
}

/// Inserts an expression (JSON) into the handle manager.
///
/// Input: JSON-serialized Expr
/// Output: JSON object with "handle" field
#[no_mangle]

pub extern "C" fn rssn_handle_insert_json(json_str : *const c_char) -> *mut c_char {

    let expr : Option<Expr> = from_json_string(json_str);

    let expr = match expr {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    let handle = HANDLE_MANAGER.insert(expr);

    let response = serde_json::json!({ "handle": handle });

    to_json_string(&response)
}

/// Retrieves an expression by handle (JSON).
///
/// Input: JSON object with "handle" field
/// Output: JSON-serialized Expr
#[no_mangle]

pub extern "C" fn rssn_handle_get_json(json_str : *const c_char) -> *mut c_char {

    #[derive(Deserialize)]

    struct Request {
        handle : usize,
    }

    let req : Option<Request> = from_json_string(json_str);

    let req = match req {
        | Some(r) => r,
        | None => return std::ptr::null_mut(),
    };

    match HANDLE_MANAGER.get(req.handle) {
        | Some(arc_expr) => to_json_string(&*arc_expr),
        | None => std::ptr::null_mut(),
    }
}

/// Checks if a handle exists (JSON).
///
/// Input: JSON object with "handle" field
/// Output: JSON object with "exists" boolean field
#[no_mangle]

pub extern "C" fn rssn_handle_exists_json(json_str : *const c_char) -> *mut c_char {

    #[derive(Deserialize)]

    struct Request {
        handle : usize,
    }

    let req : Option<Request> = from_json_string(json_str);

    let req = match req {
        | Some(r) => r,
        | None => return std::ptr::null_mut(),
    };

    let exists = HANDLE_MANAGER.exists(req.handle);

    let response = serde_json::json!({ "exists": exists });

    to_json_string(&response)
}

/// Frees a handle (JSON).
///
/// Input: JSON object with "handle" field
/// Output: JSON object with "freed" boolean field
#[no_mangle]

pub extern "C" fn rssn_handle_free_json(json_str : *const c_char) -> *mut c_char {

    #[derive(Deserialize)]

    struct Request {
        handle : usize,
    }

    let req : Option<Request> = from_json_string(json_str);

    let req = match req {
        | Some(r) => r,
        | None => return std::ptr::null_mut(),
    };

    let freed = HANDLE_MANAGER
        .free(req.handle)
        .is_some();

    let response = serde_json::json!({ "freed": freed });

    to_json_string(&response)
}

/// Returns handle manager statistics (JSON).
///
/// Output: JSON object with "count" and "handles" fields
#[no_mangle]

pub extern "C" fn rssn_handle_stats_json() -> *mut c_char {

    let handles = HANDLE_MANAGER.get_all_handles();

    let mut handle_infos = Vec::new();

    for &handle in &handles {

        if let Some(arc_expr) = HANDLE_MANAGER.get(handle) {

            handle_infos.push(HandleInfo {
                handle,
                expression : format!("{}", arc_expr),
            });
        }
    }

    let response = HandleStatsResponse {
        count : handles.len(),
        handles : handle_infos,
    };

    to_json_string(&response)
}

/// Clears all handles (JSON).
///
/// Output: JSON object with "cleared" boolean field
#[no_mangle]

pub extern "C" fn rssn_handle_clear_json() -> *mut c_char {

    HANDLE_MANAGER.clear();

    let response = serde_json::json!({ "cleared": true });

    to_json_string(&response)
}

/// Clones a handle (JSON).
///
/// Input: JSON object with "handle" field
/// Output: JSON object with "new_handle" field
#[no_mangle]

pub extern "C" fn rssn_handle_clone_json(json_str : *const c_char) -> *mut c_char {

    #[derive(Deserialize)]

    struct Request {
        handle : usize,
    }

    let req : Option<Request> = from_json_string(json_str);

    let req = match req {
        | Some(r) => r,
        | None => return std::ptr::null_mut(),
    };

    match HANDLE_MANAGER.clone_expr(req.handle) {
        | Some(expr) => {

            let new_handle = HANDLE_MANAGER.insert(expr);

            let response = serde_json::json!({ "new_handle": new_handle });

            to_json_string(&response)
        },
        | None => std::ptr::null_mut(),
    }
}
