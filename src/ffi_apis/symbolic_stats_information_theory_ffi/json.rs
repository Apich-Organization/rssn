use std::os::raw::c_char;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_json_string;
use crate::symbolic::core::Expr;
use crate::symbolic::stats_information_theory;

/// Computes the Shannon entropy of a probability distribution.

///

/// Takes a JSON string representing a `Vec<Expr>` (probabilities).

/// Returns a JSON string representing the `Expr` of the entropy.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_shannon_entropy(
    probs_json: *const c_char
) -> *mut c_char {

    let probs: Option<Vec<Expr>> =
        from_json_string(probs_json);

    if let Some(p) = probs {

        let res = stats_information_theory::shannon_entropy(&p);

        to_json_string(&res)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the Kullback-Leibler divergence between two probability distributions.

///

/// Takes JSON strings representing two `Vec<Expr>` (probability distributions `p` and `q`).

/// Returns a JSON string representing the `Expr` of the KL divergence.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_kl_divergence(
    p_probs_json: *const c_char,
    q_probs_json: *const c_char,
) -> *mut c_char {

    let p_probs: Option<Vec<Expr>> =
        from_json_string(p_probs_json);

    let q_probs: Option<Vec<Expr>> =
        from_json_string(q_probs_json);

    if let (Some(p), Some(q)) =
        (p_probs, q_probs)
    {

        match stats_information_theory::kl_divergence(&p, &q) {
            | Ok(res) => to_json_string(&res),
            | Err(_) => std::ptr::null_mut(),
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the cross-entropy between two probability distributions.

///

/// Takes two JSON strings representing `Vec<Expr>` (probability distributions `p` and `q`).

/// Returns a JSON string representing the `Expr` of the cross-entropy.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_cross_entropy(
    p_probs_json: *const c_char,
    q_probs_json: *const c_char,
) -> *mut c_char {

    let p_probs: Option<Vec<Expr>> =
        from_json_string(p_probs_json);

    let q_probs: Option<Vec<Expr>> =
        from_json_string(q_probs_json);

    if let (Some(p), Some(q)) =
        (p_probs, q_probs)
    {

        match stats_information_theory::cross_entropy(&p, &q) {
            | Ok(res) => to_json_string(&res),
            | Err(_) => std::ptr::null_mut(),
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the Gini impurity of a probability distribution.

///

/// Takes a JSON string representing a `Vec<Expr>` (probabilities).

/// Returns a JSON string representing the `Expr` of the Gini impurity.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_gini_impurity(
    probs_json: *const c_char
) -> *mut c_char {

    let probs: Option<Vec<Expr>> =
        from_json_string(probs_json);

    if let Some(p) = probs {

        let res = stats_information_theory::gini_impurity(&p);

        to_json_string(&res)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the joint entropy of a joint probability distribution.

///

/// Takes a JSON string representing an `Expr` (joint probability distribution).

/// Returns a JSON string representing the `Expr` of the joint entropy.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_joint_entropy(
    joint_probs_json: *const c_char
) -> *mut c_char {

    let joint: Option<Expr> =
        from_json_string(
            joint_probs_json,
        );

    if let Some(j) = joint {

        match stats_information_theory::joint_entropy(&j) {
            | Ok(res) => to_json_string(&res),
            | Err(_) => std::ptr::null_mut(),
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the conditional entropy of a joint probability distribution.

///

/// Takes a JSON string representing an `Expr` (joint probability distribution).

/// Returns a JSON string representing the `Expr` of the conditional entropy.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_conditional_entropy(
    joint_probs_json: *const c_char
) -> *mut c_char {

    let joint: Option<Expr> =
        from_json_string(
            joint_probs_json,
        );

    if let Some(j) = joint {

        match stats_information_theory::conditional_entropy(&j) {
            | Ok(res) => to_json_string(&res),
            | Err(_) => std::ptr::null_mut(),
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the mutual information between two random variables from their joint probability distribution.

///

/// Takes a JSON string representing an `Expr` (joint probability distribution).

/// Returns a JSON string representing the `Expr` of the mutual information.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_mutual_information(
    joint_probs_json: *const c_char
) -> *mut c_char {

    let joint: Option<Expr> =
        from_json_string(
            joint_probs_json,
        );

    if let Some(j) = joint {

        match stats_information_theory::mutual_information(&j) {
            | Ok(res) => to_json_string(&res),
            | Err(_) => std::ptr::null_mut(),
        }
    } else {

        std::ptr::null_mut()
    }
}
