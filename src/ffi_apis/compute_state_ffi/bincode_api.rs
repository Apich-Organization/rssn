//! Bincode-based FFI API for compute state module.

use crate::compute::state::State;
use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;

/// Creates a new State and returns it as a bincode buffer.
/// The caller must free the returned buffer using rssn_free_bincode_buffer.
#[no_mangle]

pub extern "C" fn rssn_state_new_bincode() -> BincodeBuffer {

    let state = State::new();

    to_bincode_buffer(&state)
}

/// Gets the intermediate value from a bincode state buffer.
/// Returns the value as a bincode buffer containing a String.
/// The caller must free the returned buffer using rssn_free_bincode_buffer.
#[no_mangle]

pub extern "C" fn rssn_state_get_intermediate_value_bincode(
    state_buffer : BincodeBuffer
) -> BincodeBuffer {

    let state : Option<State> = from_bincode_buffer(&state_buffer);

    match state {
        | Some(s) => to_bincode_buffer(&s.intermediate_value),
        | None => BincodeBuffer::empty(),
    }
}

/// Sets the intermediate value in a bincode state buffer and returns the updated buffer.
/// The caller must free the returned buffer using rssn_free_bincode_buffer.
#[no_mangle]

pub extern "C" fn rssn_state_set_intermediate_value_bincode(
    state_buffer : BincodeBuffer,
    value_buffer : BincodeBuffer,
) -> BincodeBuffer {

    let state : Option<State> = from_bincode_buffer(&state_buffer);

    let value : Option<String> = from_bincode_buffer(&value_buffer);

    match (state, value) {
        | (Some(mut s), Some(v)) => {

            s.intermediate_value = v;

            to_bincode_buffer(&s)
        },
        | _ => BincodeBuffer::empty(),
    }
}
