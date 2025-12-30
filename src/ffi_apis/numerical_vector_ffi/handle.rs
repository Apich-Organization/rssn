//! Raw pointer-based FFI API for numerical vector operations.


use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::vector;

/// Creates a new numerical vector from a raw array of doubles.
/// The caller is responsible for freeing the returned pointer using `rssn_num_vec_free`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_vec_create(
    data: *const f64,
    len: usize,
) -> *mut Vec<f64> { unsafe {

    if data.is_null() {

        update_last_error(
            "Null pointer passed to \
             rssn_num_vec_create"
                .to_string(),
        );

        return std::ptr::null_mut();
    }

    let v =  {

        std::slice::from_raw_parts(
            data, len,
        )
    }
    .to_vec();

    Box::into_raw(Box::new(v))
}}

/// Frees a numerical vector allocated by the library.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_vec_free(
    v: *mut Vec<f64>
) {

    if !v.is_null() {

        unsafe {

            let _ = Box::from_raw(v);
        }
    }
}

/// Returns the length of a numerical vector.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub const unsafe extern "C" fn rssn_num_vec_len(
    v: *const Vec<f64>
) -> usize {

    if v.is_null() {

        return 0;
    }

    unsafe {

        (*v).len()
    }
}

/// Returns a pointer to the underlying data of a numerical vector.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub const unsafe extern "C" fn rssn_num_vec_data(
    v: *const Vec<f64>
) -> *const f64 {

    if v.is_null() {

        return std::ptr::null();
    }

    unsafe {

        (*v).as_ptr()
    }
}

/// Computes the sum of two vectors.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_vec_add(
    v1: *const Vec<f64>,
    v2: *const Vec<f64>,
) -> *mut Vec<f64> { unsafe {

    if v1.is_null() || v2.is_null() {

        return std::ptr::null_mut();
    }

    let res = vector::vec_add(
         {

            &*v1
        },
         {

            &*v2
        },
    );

    match res {
        | Ok(v) => {
            Box::into_raw(Box::new(v))
        },
        | Err(e) => {

            update_last_error(e);

            std::ptr::null_mut()
        },
    }
}}

/// Computes the difference of two vectors.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_vec_sub(
    v1: *const Vec<f64>,
    v2: *const Vec<f64>,
) -> *mut Vec<f64> { unsafe {

    if v1.is_null() || v2.is_null() {

        return std::ptr::null_mut();
    }

    let res = vector::vec_sub(
         {

            &*v1
        },
         {

            &*v2
        },
    );

    match res {
        | Ok(v) => {
            Box::into_raw(Box::new(v))
        },
        | Err(e) => {

            update_last_error(e);

            std::ptr::null_mut()
        },
    }
}}

/// Multiplies a vector by a scalar.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_vec_scalar_mul(
    v: *const Vec<f64>,
    s: f64,
) -> *mut Vec<f64> {

    if v.is_null() {

        return std::ptr::null_mut();
    }

    let res = vector::scalar_mul(
        unsafe {

            &*v
        },
        s,
    );

    Box::into_raw(Box::new(res))
}

/// Computes the dot product of two vectors.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_vec_dot_product(
    v1: *const Vec<f64>,
    v2: *const Vec<f64>,
    result: *mut f64,
) -> i32 { unsafe {

    if v1.is_null()
        || v2.is_null()
        || result.is_null()
    {

        return -1;
    }

    match vector::dot_product(
         {

            &*v1
        },
         {

            &*v2
        },
    ) {
        | Ok(val) => {

             {

                *result = val;
            }

            0
        },
        | Err(e) => {

            update_last_error(e);

            -1
        },
    }
}}

/// Computes the L2 norm of a vector.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_vec_norm(
    v: *const Vec<f64>,
    result: *mut f64,
) -> i32 {

    if v.is_null() || result.is_null() {

        return -1;
    }

    unsafe {

        *result = vector::norm(&*v);
    }

    0
}

/// Computes the Lp norm of a vector.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_vec_lp_norm(
    v: *const Vec<f64>,
    p: f64,
    result: *mut f64,
) -> i32 {

    if v.is_null() || result.is_null() {

        return -1;
    }

    unsafe {

        *result =
            vector::lp_norm(&*v, p);
    }

    0
}

/// Normalizes a vector.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_vec_normalize(
    v: *const Vec<f64>
) -> *mut Vec<f64> { unsafe {

    if v.is_null() {

        return std::ptr::null_mut();
    }

    match vector::normalize( {

        &*v
    }) {
        | Ok(res) => {
            Box::into_raw(Box::new(res))
        },
        | Err(e) => {

            update_last_error(e);

            std::ptr::null_mut()
        },
    }
}}

/// Computes the cross product of two 3D vectors.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_vec_cross_product(
    v1: *const Vec<f64>,
    v2: *const Vec<f64>,
) -> *mut Vec<f64> { unsafe {

    if v1.is_null() || v2.is_null() {

        return std::ptr::null_mut();
    }

    match vector::cross_product(
         {

            &*v1
        },
         {

            &*v2
        },
    ) {
        | Ok(res) => {
            Box::into_raw(Box::new(res))
        },
        | Err(e) => {

            update_last_error(e);

            std::ptr::null_mut()
        },
    }
}}

/// Computes the angle between two vectors.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_vec_angle(
    v1: *const Vec<f64>,
    v2: *const Vec<f64>,
    result: *mut f64,
) -> i32 { unsafe {

    if v1.is_null()
        || v2.is_null()
        || result.is_null()
    {

        return -1;
    }

    match vector::angle(
         {

            &*v1
        },
         {

            &*v2
        },
    ) {
        | Ok(val) => {

             {

                *result = val;
            }

            0
        },
        | Err(e) => {

            update_last_error(e);

            -1
        },
    }
}}

/// Projects v1 onto v2.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_vec_project(
    v1: *const Vec<f64>,
    v2: *const Vec<f64>,
) -> *mut Vec<f64> { unsafe {

    if v1.is_null() || v2.is_null() {

        return std::ptr::null_mut();
    }

    match vector::project(
         {

            &*v1
        },
         {

            &*v2
        },
    ) {
        | Ok(res) => {
            Box::into_raw(Box::new(res))
        },
        | Err(e) => {

            update_last_error(e);

            std::ptr::null_mut()
        },
    }
}}

/// Reflects v about n.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_vec_reflect(
    v: *const Vec<f64>,
    n: *const Vec<f64>,
) -> *mut Vec<f64> { unsafe {

    if v.is_null() || n.is_null() {

        return std::ptr::null_mut();
    }

    match vector::reflect(
         {

            &*v
        },
         {

            &*n
        },
    ) {
        | Ok(res) => {
            Box::into_raw(Box::new(res))
        },
        | Err(e) => {

            update_last_error(e);

            std::ptr::null_mut()
        },
    }
}}

// Correction for reflect: vector::reflect(&*v, &*n)
