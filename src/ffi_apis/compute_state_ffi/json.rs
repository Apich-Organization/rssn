//! JSON-based FFI API for compute state module.

use std::os::raw::c_char;

use crate::compute::state::State;
use crate::ffi_apis::common::to_c_string;

/// Creates a new State and returns it as a JSON string.
/// The caller must free the returned string using rssn_free_string.
#[no_mangle]

pub extern "C" fn rssn_state_new_json(
) -> *mut c_char {

    let state = State::new();

    match serde_json::to_string(&state)
    {
        | Ok(json) => to_c_string(json),
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

/// Gets the intermediate value from a JSON state string.
/// Returns the value as a plain string (not JSON-encoded).
/// The caller must free the returned string using rssn_free_string.
#[no_mangle]

pub extern "C" fn rssn_state_get_intermediate_value_json(
    json_state : *const c_char
) -> *mut c_char {

    if json_state.is_null() {

        return std::ptr::null_mut();
    }

    unsafe {

        let c_str =
            std::ffi::CStr::from_ptr(
                json_state,
            );

        if let Ok(json_str) =
            c_str.to_str()
        {

            if let Ok(state) =
                serde_json::from_str::<
                    State,
                >(
                    json_str
                )
            {

                return to_c_string(
                    state.intermediate_value,
                );
            }
        }

        std::ptr::null_mut()
    }
}

/// Sets the intermediate value in a JSON state string and returns the updated JSON.
/// The caller must free the returned string using rssn_free_string.
#[no_mangle]

pub extern "C" fn rssn_state_set_intermediate_value_json(
    json_state : *const c_char,
    value : *const c_char,
) -> *mut c_char {

    if json_state.is_null()
        || value.is_null()
    {

        return std::ptr::null_mut();
    }

    unsafe {

        let state_str =
            match std::ffi::CStr::from_ptr(json_state)
                .to_str()
            {
                | Ok(s) => s,
                | Err(_) => return std::ptr::null_mut(),
            };

        let value_str = match std::ffi::CStr::from_ptr(
            value,
        )
        .to_str()
        {
            | Ok(s) => s,
            | Err(_) => return std::ptr::null_mut(),
        };

        let mut state: State =
            match serde_json::from_str(state_str) {
                | Ok(s) => s,
                | Err(_) => return std::ptr::null_mut(),
            };

        state.intermediate_value =
            value_str.to_string();

        match serde_json::to_string(
            &state,
        ) {
            | Ok(json) => {
                to_c_string(json)
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    }
}
