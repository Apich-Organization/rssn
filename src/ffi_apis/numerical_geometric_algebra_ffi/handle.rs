//! Handle-based FFI API for numerical geometric algebra operations.

use std::ptr;

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::geometric_algebra::Multivector3D;

/// Creates a new `Multivector3D`.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_ga_create(
    s: f64,
    v1: f64,
    v2: f64,
    v3: f64,
    b12: f64,
    b23: f64,
    b31: f64,
    pss: f64,
) -> *mut Multivector3D {

    let mv = Multivector3D::new(
        s, v1, v2, v3, b12, b23, b31,
        pss,
    );

    Box::into_raw(Box::new(mv))
}

/// Frees a `Multivector3D`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_ga_free(
    mv: *mut Multivector3D
) {

    if !mv.is_null() {

        unsafe {

            let _ = Box::from_raw(mv);
        }
    }
}

/// Gets components of a `Multivector3D`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_ga_get_components(
    mv: *const Multivector3D,
    s: *mut f64,
    v1: *mut f64,
    v2: *mut f64,
    v3: *mut f64,
    b12: *mut f64,
    b23: *mut f64,
    b31: *mut f64,
    pss: *mut f64,
) -> i32 {

    unsafe {

        if mv.is_null() {

            update_last_error("Null pointer passed to rssn_num_ga_get_components".to_string());

            return -1;
        }

        let m = {

            &*mv
        };

        if !s.is_null() {

            *s = m.s;
        }

        if !v1.is_null() {

            *v1 = m.v1;
        }

        if !v2.is_null() {

            *v2 = m.v2;
        }

        if !v3.is_null() {

            *v3 = m.v3;
        }

        if !b12.is_null() {

            *b12 = m.b12;
        }

        if !b23.is_null() {

            *b23 = m.b23;
        }

        if !b31.is_null() {

            *b31 = m.b31;
        }

        if !pss.is_null() {

            *pss = m.pss;
        }

        0
    }
}

/// Performs multivector addition.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_ga_add(
    mv1: *const Multivector3D,
    mv2: *const Multivector3D,
) -> *mut Multivector3D {

    if mv1.is_null() || mv2.is_null() {

        return ptr::null_mut();
    }

    let a = unsafe {

        &*mv1
    };

    let b = unsafe {

        &*mv2
    };

    Box::into_raw(Box::new(*a + *b))
}

/// Performs multivector subtraction.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_ga_sub(
    mv1: *const Multivector3D,
    mv2: *const Multivector3D,
) -> *mut Multivector3D {

    if mv1.is_null() || mv2.is_null() {

        return ptr::null_mut();
    }

    let a = unsafe {

        &*mv1
    };

    let b = unsafe {

        &*mv2
    };

    Box::into_raw(Box::new(*a - *b))
}

/// Performs geometric product.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_ga_mul(
    mv1: *const Multivector3D,
    mv2: *const Multivector3D,
) -> *mut Multivector3D {

    if mv1.is_null() || mv2.is_null() {

        return ptr::null_mut();
    }

    let a = unsafe {

        &*mv1
    };

    let b = unsafe {

        &*mv2
    };

    Box::into_raw(Box::new(*a * *b))
}

/// Performs outer product.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_ga_wedge(
    mv1: *const Multivector3D,
    mv2: *const Multivector3D,
) -> *mut Multivector3D {

    if mv1.is_null() || mv2.is_null() {

        return ptr::null_mut();
    }

    let a = unsafe {

        &*mv1
    };

    let b = unsafe {

        &*mv2
    };

    Box::into_raw(Box::new(
        a.wedge(*b),
    ))
}

/// Performs inner product.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_ga_dot(
    mv1: *const Multivector3D,
    mv2: *const Multivector3D,
) -> *mut Multivector3D {

    if mv1.is_null() || mv2.is_null() {

        return ptr::null_mut();
    }

    let a = unsafe {

        &*mv1
    };

    let b = unsafe {

        &*mv2
    };

    Box::into_raw(Box::new(a.dot(*b)))
}

/// Returns the reverse of a `Multivector3D`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_ga_reverse(
    mv: *const Multivector3D
) -> *mut Multivector3D {

    if mv.is_null() {

        return ptr::null_mut();
    }

    let a = unsafe {

        &*mv
    };

    Box::into_raw(Box::new(
        a.reverse(),
    ))
}

/// Returns the norm of a `Multivector3D`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_ga_norm(
    mv: *const Multivector3D
) -> f64 {

    if mv.is_null() {

        return 0.0;
    }

    let a = unsafe {

        &*mv
    };

    a.norm()
}

/// Returns the inverse of a `Multivector3D`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_ga_inv(
    mv: *const Multivector3D
) -> *mut Multivector3D {

    unsafe {

        if mv.is_null() {

            return ptr::null_mut();
        }

        let a = {

            &*mv
        };

        match a.inv() {
            | Some(res) => {
                Box::into_raw(Box::new(
                    res,
                ))
            },
            | None => {

                update_last_error(
                    "Multivector is \
                     not invertible"
                        .to_string(),
                );

                ptr::null_mut()
            },
        }
    }
}
