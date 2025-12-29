//! Handle-based FFI API for compute state module.

use std::ffi::CStr;
use std::ffi::CString;
use std::os::raw::c_char;

use crate::compute::state::State;

/// Creates a new State.
/// The caller is responsible for freeing the memory using `rssn_state_free`.
#[no_mangle]

pub extern "C" fn rssn_state_new(
) -> *mut State {

    Box::into_raw(Box::new(
        State::new(),
    ))
}

/// Frees a State.
#[no_mangle]

pub extern "C" fn rssn_state_free(
    state: *mut State
) {

    if state.is_null() {

        return;
    }

    unsafe {

        let _ = Box::from_raw(state);
    }
}

/// Gets the intermediate value from the state.
/// The returned string must be freed by the caller using `rssn_free_string`.
#[no_mangle]

pub extern "C" fn rssn_state_get_intermediate_value(
    state: *const State
) -> *mut c_char {

    if state.is_null() {

        return std::ptr::null_mut();
    }

    unsafe {

        let s = &(*state)
            .intermediate_value;

        match CString::new(s.as_str()) {
            | Ok(c_str) => {
                c_str.into_raw()
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    }
}

/// Sets the intermediate value in the state.
#[no_mangle]

pub extern "C" fn rssn_state_set_intermediate_value(
    state: *mut State,
    value: *const c_char,
) {

    if state.is_null()
        || value.is_null()
    {

        return;
    }

    unsafe {

        let c_str =
            CStr::from_ptr(value);

        if let Ok(s) = c_str.to_str() {

            (*state)
                .intermediate_value =
                s.to_string();
        }
    }
}
