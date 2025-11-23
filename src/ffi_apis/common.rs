//! Common FFI utilities shared across all FFI modules.
//!
//! This module provides shared types and functions for the three FFI API versions:
//! - Handle-based (opaque pointers)
//! - JSON-based (string serialization)
//! - Bincode-based (binary serialization)

use std::ffi::CString;
use std::os::raw::c_char;

/// A buffer containing binary data from bincode serialization.
///
/// The caller is responsible for freeing this buffer using `rssn_free_bincode_buffer`.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct BincodeBuffer {
    /// Pointer to the binary data.
    pub data: *mut u8,
    /// Length of the binary data in bytes.
    pub len: usize,
}

impl BincodeBuffer {
    /// Creates a new empty buffer.
    pub fn empty() -> Self {
        Self {
            data: std::ptr::null_mut(),
            len: 0,
        }
    }

    /// Creates a buffer from a Vec<u8>.
    pub fn from_vec(bytes: Vec<u8>) -> Self {
        let len = bytes.len();
        let data = Box::into_raw(bytes.into_boxed_slice()) as *mut u8;
        Self { data, len }
    }

    /// Checks if the buffer is null/empty.
    pub fn is_null(&self) -> bool {
        self.data.is_null() || self.len == 0
    }

    /// Converts the buffer to a slice (unsafe).
    ///
    /// # Safety
    /// The buffer must be valid and not yet freed.
    pub unsafe fn as_slice(&self) -> &[u8] {
        if self.is_null() {
            &[]
        } else {
            std::slice::from_raw_parts(self.data, self.len)
        }
    }
}

/// Frees a string allocated by an FFI function.
///
/// # Safety
/// The string must have been allocated by an FFI function that returns `*mut c_char`.
/// This function should only be called once per string.
#[no_mangle]
pub extern "C" fn rssn_free_string(s: *mut c_char) {
    if !s.is_null() {
        unsafe {
            let _ = CString::from_raw(s);
        }
    }
}

/// Frees a bincode buffer allocated by an FFI function.
///
/// # Safety
/// The buffer must have been allocated by an FFI function that returns `BincodeBuffer`.
/// This function should only be called once per buffer.
#[no_mangle]
pub extern "C" fn rssn_free_bincode_buffer(buffer: BincodeBuffer) {
    if !buffer.is_null() {
        unsafe {
            let _ = Box::from_raw(std::slice::from_raw_parts_mut(buffer.data, buffer.len));
        }
    }
}

/// Helper function to create a C string from a Rust string.
///
/// Returns null on error.
pub fn to_c_string(s: String) -> *mut c_char {
    match CString::new(s) {
        Ok(c_str) => c_str.into_raw(),
        Err(_) => std::ptr::null_mut(),
    }
}

/// Helper function to serialize to JSON and return as C string.
///
/// Returns null on error.
pub fn to_json_string<T: serde::Serialize>(value: &T) -> *mut c_char {
    match serde_json::to_string(value) {
        Ok(json) => to_c_string(json),
        Err(_) => std::ptr::null_mut(),
    }
}

/// Helper function to deserialize from JSON C string.
///
/// Returns None on error.
pub fn from_json_string<T: serde::de::DeserializeOwned>(json: *const c_char) -> Option<T> {
    if json.is_null() {
        return None;
    }
    unsafe {
        let c_str = std::ffi::CStr::from_ptr(json);
        c_str
            .to_str()
            .ok()
            .and_then(|s| serde_json::from_str(s).ok())
    }
}

/// Helper function to serialize to bincode and return as buffer.
///
/// Returns empty buffer on error.
pub fn to_bincode_buffer<T: serde::Serialize>(value: &T) -> BincodeBuffer {
    match bincode::serde::encode_to_vec(value, bincode::config::standard()) {
        Ok(bytes) => BincodeBuffer::from_vec(bytes),
        Err(_) => BincodeBuffer::empty(),
    }
}

/// Helper function to deserialize from bincode buffer.
///
/// Returns None on error.
pub fn from_bincode_buffer<T: serde::de::DeserializeOwned>(buffer: &BincodeBuffer) -> Option<T> {
    if buffer.is_null() {
        return None;
    }
    unsafe {
        let slice = buffer.as_slice();
        bincode::serde::decode_from_slice(slice, bincode::config::standard())
            .ok()
            .map(|(v, _)| v)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bincode_buffer_empty() {
        let buffer = BincodeBuffer::empty();
        assert!(buffer.is_null());
    }

    #[test]
    fn test_bincode_buffer_from_vec() {
        let vec = vec![1, 2, 3, 4];
        let buffer = BincodeBuffer::from_vec(vec);
        assert!(!buffer.is_null());
        assert_eq!(buffer.len, 4);
        unsafe {
            assert_eq!(buffer.as_slice(), &[1, 2, 3, 4]);
        }
        rssn_free_bincode_buffer(buffer);
    }

    #[test]
    fn test_to_c_string() {
        let s = "Hello, World!".to_string();
        let c_str = to_c_string(s);
        assert!(!c_str.is_null());
        rssn_free_string(c_str);
    }
}
