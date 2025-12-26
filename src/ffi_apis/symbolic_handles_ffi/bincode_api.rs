//! Bincode-based FFI API for the HandleManager.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::handles::HANDLE_MANAGER;

#[derive(Serialize, Deserialize)]

struct HandleResponse {
    handle : usize,
}

#[derive(Serialize, Deserialize)]

struct HandleRequest {
    handle : usize,
}

#[derive(Serialize, Deserialize)]

struct ExistsResponse {
    exists : bool,
}

#[derive(Serialize, Deserialize)]

struct FreedResponse {
    freed : bool,
}

#[derive(Serialize, Deserialize)]

struct HandleListResponse {
    handles : Vec<usize>,
    count : usize,
}

/// Inserts an expression (Bincode) into the handle manager.
#[no_mangle]

pub extern "C" fn rssn_handle_insert_bincode(input : BincodeBuffer) -> BincodeBuffer {

    let expr : Option<Expr> = from_bincode_buffer(&input);

    let expr = match expr {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let handle = HANDLE_MANAGER.insert(expr);

    let response = HandleResponse {
        handle,
    };

    to_bincode_buffer(&response)
}

/// Retrieves an expression by handle (Bincode).
#[no_mangle]

pub extern "C" fn rssn_handle_get_bincode(input : BincodeBuffer) -> BincodeBuffer {

    let req : Option<HandleRequest> = from_bincode_buffer(&input);

    let req = match req {
        | Some(r) => r,
        | None => return BincodeBuffer::empty(),
    };

    match HANDLE_MANAGER.get(req.handle) {
        | Some(arc_expr) => to_bincode_buffer(&*arc_expr),
        | None => BincodeBuffer::empty(),
    }
}

/// Checks if a handle exists (Bincode).
#[no_mangle]

pub extern "C" fn rssn_handle_exists_bincode(input : BincodeBuffer) -> BincodeBuffer {

    let req : Option<HandleRequest> = from_bincode_buffer(&input);

    let req = match req {
        | Some(r) => r,
        | None => return BincodeBuffer::empty(),
    };

    let exists = HANDLE_MANAGER.exists(req.handle);

    let response = ExistsResponse {
        exists,
    };

    to_bincode_buffer(&response)
}

/// Frees a handle (Bincode).
#[no_mangle]

pub extern "C" fn rssn_handle_free_bincode(input : BincodeBuffer) -> BincodeBuffer {

    let req : Option<HandleRequest> = from_bincode_buffer(&input);

    let req = match req {
        | Some(r) => r,
        | None => return BincodeBuffer::empty(),
    };

    let freed = HANDLE_MANAGER
        .free(req.handle)
        .is_some();

    let response = FreedResponse {
        freed,
    };

    to_bincode_buffer(&response)
}

/// Returns all active handles (Bincode).
#[no_mangle]

pub extern "C" fn rssn_handle_get_all_bincode() -> BincodeBuffer {

    let handles = HANDLE_MANAGER.get_all_handles();

    let response = HandleListResponse {
        count : handles.len(),
        handles,
    };

    to_bincode_buffer(&response)
}

/// Clears all handles (Bincode).
#[no_mangle]

pub extern "C" fn rssn_handle_clear_bincode() -> BincodeBuffer {

    HANDLE_MANAGER.clear();

    let response = serde_json::json!({ "cleared": true });

    to_bincode_buffer(&response)
}

/// Clones a handle (Bincode).
#[no_mangle]

pub extern "C" fn rssn_handle_clone_bincode(input : BincodeBuffer) -> BincodeBuffer {

    let req : Option<HandleRequest> = from_bincode_buffer(&input);

    let req = match req {
        | Some(r) => r,
        | None => return BincodeBuffer::empty(),
    };

    match HANDLE_MANAGER.clone_expr(req.handle) {
        | Some(expr) => {

            let new_handle = HANDLE_MANAGER.insert(expr);

            let response = HandleResponse {
                handle : new_handle,
            };

            to_bincode_buffer(&response)
        },
        | None => BincodeBuffer::empty(),
    }
}
