# FFI API Strategy for RSSN

## Overview

All FFI APIs in RSSN must provide **three different versions** to support various use cases and performance requirements:

1. **Handle-based API** - Traditional C-style opaque pointers
2. **JSON-based API** - String-based serialization for easy interop
3. **Bincode-based API** - Binary serialization for performance

## Design Principles

### 1. Handle-based API (Primary)
- **Purpose**: Direct memory access, maximum performance
- **Use case**: Native C/C++/Fortran applications
- **Naming**: `rssn_<module>_<function>`
- **Example**: `rssn_state_new()`, `rssn_state_free()`
- **Memory management**: Caller responsible for freeing

### 2. JSON-based API
- **Purpose**: Language-agnostic, human-readable
- **Use case**: Scripting languages, debugging, web services
- **Naming**: `rssn_<module>_<function>_json`
- **Example**: `rssn_state_new_json()` returns JSON string
- **Memory management**: Always returns allocated strings that must be freed with `rssn_free_string()`

### 3. Bincode-based API
- **Purpose**: Efficient binary serialization
- **Use case**: High-performance cross-language communication
- **Naming**: `rssn_<module>_<function>_bincode`
- **Example**: `rssn_state_new_bincode()` returns binary buffer
- **Memory management**: Returns buffer + length, must be freed with `rssn_free_buffer()`

## Implementation Pattern

For each module (e.g., `state`), create three sub-modules:

```rust
// src/ffi_apis/state_ffi/mod.rs
pub mod handle;  // Handle-based API
pub mod json;    // JSON-based API
pub mod bincode; // Bincode-based API

// Re-export all
pub use handle::*;
pub use json::*;
pub use bincode::*;
```

### Handle-based Implementation

```rust
// src/ffi_apis/state_ffi/handle.rs
use crate::compute::state::State;

#[no_mangle]
pub extern "C" fn rssn_state_new() -> *mut State {
    Box::into_raw(Box::new(State::new()))
}

#[no_mangle]
pub extern "C" fn rssn_state_free(state: *mut State) {
    if !state.is_null() {
        unsafe { let _ = Box::from_raw(state); }
    }
}
```

### JSON-based Implementation

```rust
// src/ffi_apis/state_ffi/json.rs
use crate::compute::state::State;
use std::ffi::CString;
use std::os::raw::c_char;

#[no_mangle]
pub extern "C" fn rssn_state_new_json() -> *mut c_char {
    let state = State::new();
    match serde_json::to_string(&state) {
        Ok(json) => CString::new(json).unwrap().into_raw(),
        Err(_) => std::ptr::null_mut(),
    }
}

#[no_mangle]
pub extern "C" fn rssn_state_from_json(json: *const c_char) -> *mut State {
    if json.is_null() {
        return std::ptr::null_mut();
    }
    unsafe {
        let c_str = std::ffi::CStr::from_ptr(json);
        match c_str.to_str() {
            Ok(s) => match serde_json::from_str::<State>(s) {
                Ok(state) => Box::into_raw(Box::new(state)),
                Err(_) => std::ptr::null_mut(),
            },
            Err(_) => std::ptr::null_mut(),
        }
    }
}

#[no_mangle]
pub extern "C" fn rssn_state_to_json(state: *const State) -> *mut c_char {
    if state.is_null() {
        return std::ptr::null_mut();
    }
    unsafe {
        match serde_json::to_string(&*state) {
            Ok(json) => CString::new(json).unwrap().into_raw(),
            Err(_) => std::ptr::null_mut(),
        }
    }
}
```

### Bincode-based Implementation

```rust
// src/ffi_apis/state_ffi/bincode.rs
use crate::compute::state::State;
use std::os::raw::c_char;

#[repr(C)]
pub struct BincodeBuffer {
    pub data: *mut u8,
    pub len: usize,
}

#[no_mangle]
pub extern "C" fn rssn_state_new_bincode() -> BincodeBuffer {
    let state = State::new();
    match bincode::serialize(&state) {
        Ok(bytes) => {
            let len = bytes.len();
            let data = Box::into_raw(bytes.into_boxed_slice()) as *mut u8;
            BincodeBuffer { data, len }
        },
        Err(_) => BincodeBuffer { data: std::ptr::null_mut(), len: 0 },
    }
}

#[no_mangle]
pub extern "C" fn rssn_state_from_bincode(data: *const u8, len: usize) -> *mut State {
    if data.is_null() || len == 0 {
        return std::ptr::null_mut();
    }
    unsafe {
        let slice = std::slice::from_raw_parts(data, len);
        match bincode::deserialize::<State>(slice) {
            Ok(state) => Box::into_raw(Box::new(state)),
            Err(_) => std::ptr::null_mut(),
        }
    }
}

#[no_mangle]
pub extern "C" fn rssn_state_to_bincode(state: *const State) -> BincodeBuffer {
    if state.is_null() {
        return BincodeBuffer { data: std::ptr::null_mut(), len: 0 };
    }
    unsafe {
        match bincode::serialize(&*state) {
            Ok(bytes) => {
                let len = bytes.len();
                let data = Box::into_raw(bytes.into_boxed_slice()) as *mut u8;
                BincodeBuffer { data, len }
            },
            Err(_) => BincodeBuffer { data: std::ptr::null_mut(), len: 0 },
        }
    }
}

#[no_mangle]
pub extern "C" fn rssn_free_bincode_buffer(buffer: BincodeBuffer) {
    if !buffer.data.is_null() && buffer.len > 0 {
        unsafe {
            let _ = Box::from_raw(std::slice::from_raw_parts_mut(buffer.data, buffer.len));
        }
    }
}
```

## Memory Management Rules

### Handle-based
- Caller allocates via `rssn_<type>_new()`
- Caller frees via `rssn_<type>_free()`

### JSON-based
- All returned strings must be freed with `rssn_free_string()`
- Input strings are borrowed, not owned

### Bincode-based
- All returned buffers must be freed with `rssn_free_bincode_buffer()`
- Input buffers are borrowed, not owned

## Common Utilities

```rust
// src/ffi_apis/common.rs

use std::ffi::CString;
use std::os::raw::c_char;

#[no_mangle]
pub extern "C" fn rssn_free_string(s: *mut c_char) {
    if !s.is_null() {
        unsafe { let _ = CString::from_raw(s); }
    }
}

#[repr(C)]
pub struct BincodeBuffer {
    pub data: *mut u8,
    pub len: usize,
}

#[no_mangle]
pub extern "C" fn rssn_free_bincode_buffer(buffer: BincodeBuffer) {
    if !buffer.data.is_null() && buffer.len > 0 {
        unsafe {
            let _ = Box::from_raw(std::slice::from_raw_parts_mut(buffer.data, buffer.len));
        }
    }
}
```

## Migration Plan

### Phase 1: Core Infrastructure (Current)
- [x] Add `bincode` dependency to Cargo.toml
- [ ] Create `src/ffi_apis/common.rs` with shared utilities
- [ ] Update existing FFI modules to use new structure

### Phase 2: Refactor Existing FFI
- [ ] Refactor `constant_ffi` to three-version structure
- [ ] Refactor `compute_cache_ffi` to three-version structure
- [ ] Refactor `compute_state_ffi` to three-version structure

### Phase 3: New FFI Modules
- [ ] Implement three-version FFI for all compute modules
- [ ] Implement three-version FFI for symbolic modules
- [ ] Implement three-version FFI for numerical modules

### Phase 4: Testing & Documentation
- [ ] Create C header files for all three APIs
- [ ] Create example programs in C/C++/Fortran
- [ ] Add comprehensive FFI tests
- [ ] Document FFI usage patterns

## Feature Flags

Add Cargo features for optional bincode support:

```toml
[features]
default = ["ffi_api"]
ffi_api = []
ffi_json = ["ffi_api", "serde_json"]
ffi_bincode = ["ffi_api", "bincode"]
ffi_full = ["ffi_json", "ffi_bincode"]
```

## Error Handling

All FFI functions should:
1. Check for null pointers
2. Return null/empty on error (or use error codes)
3. Never panic across FFI boundary
4. Log errors internally if possible

## Testing Strategy

For each module, create:
1. Unit tests in Rust
2. Integration tests calling FFI from Rust
3. Example C programs demonstrating usage
4. Benchmark comparisons between the three versions

## Performance Considerations

- **Handle-based**: Fastest, zero serialization overhead
- **JSON-based**: Slowest, human-readable, ~10-100x slower
- **Bincode-based**: Fast, binary, ~2-5x slower than handle-based

## Security Considerations

- Always validate buffer lengths
- Check for null pointers before dereferencing
- Use `CStr::from_ptr()` safely
- Avoid buffer overflows in bincode operations
- Consider adding checksums for bincode buffers
