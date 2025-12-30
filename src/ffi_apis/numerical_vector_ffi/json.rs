//! JSON-based FFI API for numerical vector operations.

use std::ffi::CStr;
use std::ffi::CString;
use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::vector;

#[derive(Deserialize)]

struct VecInput {
    v: Vec<f64>,
}

#[derive(Deserialize)]

struct TwoVecInput {
    v1: Vec<f64>,
    v2: Vec<f64>,
}

#[derive(Deserialize)]

struct VecScalarInput {
    v: Vec<f64>,
    s: f64,
}

#[derive(Deserialize)]

struct VecNormInput {
    v: Vec<f64>,
    p: f64,
}

#[derive(Deserialize)]

struct VecEpsilonInput {
    v1: Vec<f64>,
    v2: Vec<f64>,
    epsilon: f64,
}

#[derive(Deserialize)]

struct LerpInput {
    v1: Vec<f64>,
    v2: Vec<f64>,
    t: f64,
}

/// JSON FFI for `vec_add`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_vec_add_json(
    json_ptr: *const c_char
) -> *mut c_char { unsafe {

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

    let input: TwoVecInput =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(v) => v,
            | Err(e) => {
                return CString::new(
                    format!(
                        "{{\"err\": \
                         \"{e}\"}}"
                    ),
                )
                .unwrap()
                .into_raw()
            },
        };

    let res = match vector::vec_add(
        &input.v1,
        &input.v2,
    ) {
        | Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None::<String>,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        },
    };

    CString::new(
        serde_json::to_string(&res)
            .unwrap(),
    )
    .unwrap()
    .into_raw()
}}

/// JSON FFI for `vec_sub`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_vec_sub_json(
    json_ptr: *const c_char
) -> *mut c_char { unsafe {

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

    let input: TwoVecInput =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(v) => v,
            | Err(e) => {
                return CString::new(
                    format!(
                        "{{\"err\": \
                         \"{e}\"}}"
                    ),
                )
                .unwrap()
                .into_raw()
            },
        };

    let res = match vector::vec_sub(
        &input.v1,
        &input.v2,
    ) {
        | Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None::<String>,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        },
    };

    CString::new(
        serde_json::to_string(&res)
            .unwrap(),
    )
    .unwrap()
    .into_raw()
}}

/// JSON FFI for `scalar_mul`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_vec_scalar_mul_json(
    json_ptr: *const c_char
) -> *mut c_char { unsafe {

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

    let input: VecScalarInput =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(v) => v,
            | Err(e) => {
                return CString::new(
                    format!(
                        "{{\"err\": \
                         \"{e}\"}}"
                    ),
                )
                .unwrap()
                .into_raw()
            },
        };

    let v = vector::scalar_mul(
        &input.v,
        input.s,
    );

    let res = FfiResult {
        ok: Some(v),
        err: None::<String>,
    };

    CString::new(
        serde_json::to_string(&res)
            .unwrap(),
    )
    .unwrap()
    .into_raw()
}}

/// JSON FFI for `dot_product`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_vec_dot_product_json(
    json_ptr: *const c_char
) -> *mut c_char { unsafe {

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

    let input: TwoVecInput =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(v) => v,
            | Err(e) => {
                return CString::new(
                    format!(
                        "{{\"err\": \
                         \"{e}\"}}"
                    ),
                )
                .unwrap()
                .into_raw()
            },
        };

    let res = match vector::dot_product(
        &input.v1,
        &input.v2,
    ) {
        | Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None::<String>,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        },
    };

    CString::new(
        serde_json::to_string(&res)
            .unwrap(),
    )
    .unwrap()
    .into_raw()
}}

/// JSON FFI for norm ($`L_2`$).
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_vec_norm_json(
    json_ptr: *const c_char
) -> *mut c_char { unsafe {

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

    let input: VecInput =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(v) => v,
            | Err(e) => {
                return CString::new(
                    format!(
                        "{{\"err\": \
                         \"{e}\"}}"
                    ),
                )
                .unwrap()
                .into_raw()
            },
        };

    let v = vector::norm(&input.v);

    let res = FfiResult {
        ok: Some(v),
        err: None::<String>,
    };

    CString::new(
        serde_json::to_string(&res)
            .unwrap(),
    )
    .unwrap()
    .into_raw()
}}

/// JSON FFI for `lp_norm`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_vec_lp_norm_json(
    json_ptr: *const c_char
) -> *mut c_char { unsafe {

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

    let input: VecNormInput =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(v) => v,
            | Err(e) => {
                return CString::new(
                    format!(
                        "{{\"err\": \
                         \"{e}\"}}"
                    ),
                )
                .unwrap()
                .into_raw()
            },
        };

    let v = vector::lp_norm(
        &input.v,
        input.p,
    );

    let res = FfiResult {
        ok: Some(v),
        err: None::<String>,
    };

    CString::new(
        serde_json::to_string(&res)
            .unwrap(),
    )
    .unwrap()
    .into_raw()
}}

/// JSON FFI for normalize.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_vec_normalize_json(
    json_ptr: *const c_char
) -> *mut c_char { unsafe {

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

    let input: VecInput =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(v) => v,
            | Err(e) => {
                return CString::new(
                    format!(
                        "{{\"err\": \
                         \"{e}\"}}"
                    ),
                )
                .unwrap()
                .into_raw()
            },
        };

    let res = match vector::normalize(
        &input.v,
    ) {
        | Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None::<String>,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        },
    };

    CString::new(
        serde_json::to_string(&res)
            .unwrap(),
    )
    .unwrap()
    .into_raw()
}}

/// JSON FFI for `cross_product`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_vec_cross_product_json(
    json_ptr: *const c_char
) -> *mut c_char { unsafe {

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

    let input: TwoVecInput =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(v) => v,
            | Err(e) => {
                return CString::new(
                    format!(
                        "{{\"err\": \
                         \"{e}\"}}"
                    ),
                )
                .unwrap()
                .into_raw()
            },
        };

    let res =
        match vector::cross_product(
            &input.v1,
            &input.v2,
        ) {
            | Ok(v) => {
                FfiResult {
                    ok: Some(v),
                    err: None::<String>,
                }
            },
            | Err(e) => {
                FfiResult {
                    ok: None,
                    err: Some(e),
                }
            },
        };

    CString::new(
        serde_json::to_string(&res)
            .unwrap(),
    )
    .unwrap()
    .into_raw()
}}

/// JSON FFI for distance.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_vec_distance_json(
    json_ptr: *const c_char
) -> *mut c_char { unsafe {

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

    let input: TwoVecInput =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(v) => v,
            | Err(e) => {
                return CString::new(
                    format!(
                        "{{\"err\": \
                         \"{e}\"}}"
                    ),
                )
                .unwrap()
                .into_raw()
            },
        };

    let res = match vector::distance(
        &input.v1,
        &input.v2,
    ) {
        | Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None::<String>,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        },
    };

    CString::new(
        serde_json::to_string(&res)
            .unwrap(),
    )
    .unwrap()
    .into_raw()
}}

/// JSON FFI for angle.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_vec_angle_json(
    json_ptr: *const c_char
) -> *mut c_char { unsafe {

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

    let input: TwoVecInput =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(v) => v,
            | Err(e) => {
                return CString::new(
                    format!(
                        "{{\"err\": \
                         \"{e}\"}}"
                    ),
                )
                .unwrap()
                .into_raw()
            },
        };

    let res = match vector::angle(
        &input.v1,
        &input.v2,
    ) {
        | Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None::<String>,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        },
    };

    CString::new(
        serde_json::to_string(&res)
            .unwrap(),
    )
    .unwrap()
    .into_raw()
}}

/// JSON FFI for project.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_vec_project_json(
    json_ptr: *const c_char
) -> *mut c_char { unsafe {

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

    let input: TwoVecInput =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(v) => v,
            | Err(e) => {
                return CString::new(
                    format!(
                        "{{\"err\": \
                         \"{e}\"}}"
                    ),
                )
                .unwrap()
                .into_raw()
            },
        };

    let res = match vector::project(
        &input.v1,
        &input.v2,
    ) {
        | Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None::<String>,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        },
    };

    CString::new(
        serde_json::to_string(&res)
            .unwrap(),
    )
    .unwrap()
    .into_raw()
}}

/// JSON FFI for reflect.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_vec_reflect_json(
    json_ptr: *const c_char
) -> *mut c_char { unsafe {

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

    let input: TwoVecInput =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(v) => v,
            | Err(e) => {
                return CString::new(
                    format!(
                        "{{\"err\": \
                         \"{e}\"}}"
                    ),
                )
                .unwrap()
                .into_raw()
            },
        };

    let res = match vector::reflect(
        &input.v1,
        &input.v2,
    ) {
        | Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None::<String>,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        },
    };

    CString::new(
        serde_json::to_string(&res)
            .unwrap(),
    )
    .unwrap()
    .into_raw()
}}

/// JSON FFI for lerp.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_vec_lerp_json(
    json_ptr: *const c_char
) -> *mut c_char { unsafe {

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

    let input: LerpInput =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(v) => v,
            | Err(e) => {
                return CString::new(
                    format!(
                        "{{\"err\": \
                         \"{e}\"}}"
                    ),
                )
                .unwrap()
                .into_raw()
            },
        };

    let res = match vector::lerp(
        &input.v1,
        &input.v2,
        input.t,
    ) {
        | Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None::<String>,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        },
    };

    CString::new(
        serde_json::to_string(&res)
            .unwrap(),
    )
    .unwrap()
    .into_raw()
}}

/// JSON FFI for `is_orthogonal`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_vec_is_orthogonal_json(
    json_ptr: *const c_char
) -> *mut c_char { unsafe {

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

    let input: VecEpsilonInput =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(v) => v,
            | Err(e) => {
                return CString::new(
                    format!(
                        "{{\"err\": \
                         \"{e}\"}}"
                    ),
                )
                .unwrap()
                .into_raw()
            },
        };

    let res =
        match vector::is_orthogonal(
            &input.v1,
            &input.v2,
            input.epsilon,
        ) {
            | Ok(v) => {
                FfiResult {
                    ok: Some(v),
                    err: None::<String>,
                }
            },
            | Err(e) => {
                FfiResult {
                    ok: None,
                    err: Some(e),
                }
            },
        };

    CString::new(
        serde_json::to_string(&res)
            .unwrap(),
    )
    .unwrap()
    .into_raw()
}}

/// JSON FFI for `is_parallel`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_vec_is_parallel_json(
    json_ptr: *const c_char
) -> *mut c_char { unsafe {

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

    let input: VecEpsilonInput =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(v) => v,
            | Err(e) => {
                return CString::new(
                    format!(
                        "{{\"err\": \
                         \"{e}\"}}"
                    ),
                )
                .unwrap()
                .into_raw()
            },
        };

    let res = match vector::is_parallel(
        &input.v1,
        &input.v2,
        input.epsilon,
    ) {
        | Ok(v) => {
            FfiResult {
                ok: Some(v),
                err: None::<String>,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        },
    };

    CString::new(
        serde_json::to_string(&res)
            .unwrap(),
    )
    .unwrap()
    .into_raw()
}}

/// JSON FFI for `cosine_similarity`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_vec_cosine_similarity_json(
    json_ptr: *const c_char
) -> *mut c_char { unsafe {

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

    let input: TwoVecInput =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(v) => v,
            | Err(e) => {
                return CString::new(
                    format!(
                        "{{\"err\": \
                         \"{e}\"}}"
                    ),
                )
                .unwrap()
                .into_raw()
            },
        };

    let res =
        match vector::cosine_similarity(
            &input.v1,
            &input.v2,
        ) {
            | Ok(v) => {
                FfiResult {
                    ok: Some(v),
                    err: None::<String>,
                }
            },
            | Err(e) => {
                FfiResult {
                    ok: None,
                    err: Some(e),
                }
            },
        };

    CString::new(
        serde_json::to_string(&res)
            .unwrap(),
    )
    .unwrap()
    .into_raw()
}}
