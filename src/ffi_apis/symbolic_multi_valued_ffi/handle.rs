use crate::symbolic::core::Expr;
use crate::symbolic::multi_valued::abs;
use crate::symbolic::multi_valued::arg;
use crate::symbolic::multi_valued::general_arccos;
use crate::symbolic::multi_valued::general_arcsin;
use crate::symbolic::multi_valued::general_arctan;
use crate::symbolic::multi_valued::general_log;
use crate::symbolic::multi_valued::general_nth_root;
use crate::symbolic::multi_valued::general_power;
use crate::symbolic::multi_valued::general_sqrt;

/// Computes general multi-valued logarithm (Handle)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_general_log_handle(
    z: *const Expr,
    k: *const Expr,
) -> *mut Expr {

    let z_ref = unsafe {

        if z.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*z
    };

    let k_ref = unsafe {

        if k.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*k
    };

    let result =
        general_log(z_ref, k_ref);

    Box::into_raw(Box::new(result))
}

/// Computes general multi-valued square root (Handle)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_general_sqrt_handle(
    z: *const Expr,
    k: *const Expr,
) -> *mut Expr {

    let z_ref = unsafe {

        if z.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*z
    };

    let k_ref = unsafe {

        if k.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*k
    };

    let result =
        general_sqrt(z_ref, k_ref);

    Box::into_raw(Box::new(result))
}

/// Computes general multi-valued power (Handle)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_general_power_handle(
    z: *const Expr,
    w: *const Expr,
    k: *const Expr,
) -> *mut Expr {

    let z_ref = unsafe {

        if z.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*z
    };

    let w_ref = unsafe {

        if w.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*w
    };

    let k_ref = unsafe {

        if k.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*k
    };

    let result = general_power(
        z_ref, w_ref, k_ref,
    );

    Box::into_raw(Box::new(result))
}

/// Computes general multi-valued n-th root (Handle)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_general_nth_root_handle(
    z: *const Expr,
    n: *const Expr,
    k: *const Expr,
) -> *mut Expr {

    let z_ref = unsafe {

        if z.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*z
    };

    let n_ref = unsafe {

        if n.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*n
    };

    let k_ref = unsafe {

        if k.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*k
    };

    let result = general_nth_root(
        z_ref, n_ref, k_ref,
    );

    Box::into_raw(Box::new(result))
}

/// Computes general multi-valued arcsin (Handle)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_general_arcsin_handle(
    z: *const Expr,
    k: *const Expr,
) -> *mut Expr {

    let z_ref = unsafe {

        if z.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*z
    };

    let k_ref = unsafe {

        if k.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*k
    };

    let result =
        general_arcsin(z_ref, k_ref);

    Box::into_raw(Box::new(result))
}

/// Computes general multi-valued arccos (Handle)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_general_arccos_handle(
    z: *const Expr,
    k: *const Expr,
    s: *const Expr,
) -> *mut Expr {

    let z_ref = unsafe {

        if z.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*z
    };

    let k_ref = unsafe {

        if k.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*k
    };

    let s_ref = unsafe {

        if s.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*s
    };

    let result = general_arccos(
        z_ref, k_ref, s_ref,
    );

    Box::into_raw(Box::new(result))
}

/// Computes general multi-valued arctan (Handle)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_general_arctan_handle(
    z: *const Expr,
    k: *const Expr,
) -> *mut Expr {

    let z_ref = unsafe {

        if z.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*z
    };

    let k_ref = unsafe {

        if k.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*k
    };

    let result =
        general_arctan(z_ref, k_ref);

    Box::into_raw(Box::new(result))
}

/// Computes argument (angle) of complex number (Handle)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_arg_handle(
    z: *const Expr
) -> *mut Expr {

    let z_ref = unsafe {

        if z.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*z
    };

    let result = arg(z_ref);

    Box::into_raw(Box::new(result))
}

/// Computes absolute value (magnitude) of complex number (Handle)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_abs_handle(
    z: *const Expr
) -> *mut Expr {

    let z_ref = unsafe {

        if z.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*z
    };

    let result = abs(z_ref);

    Box::into_raw(Box::new(result))
}
