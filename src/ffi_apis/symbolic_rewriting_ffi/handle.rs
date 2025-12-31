//! Handle-based FFI API for term rewriting systems.

use std::os::raw::c_char;

use crate::ffi_apis::common::to_c_string;
use crate::symbolic::core::Expr;
use crate::symbolic::rewriting::RewriteRule;
use crate::symbolic::rewriting::apply_rules_to_normal_form;
use crate::symbolic::rewriting::knuth_bendix;

/// Creates a new rewrite rule from lhs and rhs expressions.
///
/// # Safety
/// The caller must ensure `lhs` and `rhs` are valid Expr pointers.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_rewrite_rule_new(
    lhs: *const Expr,
    rhs: *const Expr,
) -> *mut RewriteRule {

    unsafe {

        if lhs.is_null()
            || rhs.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        let rule = RewriteRule {
            lhs: (*lhs).clone(),
            rhs: (*rhs).clone(),
        };

        Box::into_raw(Box::new(rule))
    }
}

/// Frees a rewrite rule.
///
/// # Safety
/// The caller must ensure `rule` was created by this module and hasn't been freed yet.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_rewrite_rule_free(
    rule: *mut RewriteRule
) {

    unsafe {

        if !rule.is_null() {

            let _ = Box::from_raw(rule);
        }
    }
}

/// Gets the LHS of a rewrite rule.
///
/// Returns a new owned Expr pointer that must be freed by the caller.
///
/// # Safety
/// The caller must ensure `rule` is a valid `RewriteRule` pointer.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_rewrite_rule_get_lhs(
    rule: *const RewriteRule
) -> *mut Expr {

    unsafe {

        if rule.is_null() {

            return std::ptr::null_mut(
            );
        }

        Box::into_raw(Box::new(
            (*rule).lhs.clone(),
        ))
    }
}

/// Gets the RHS of a rewrite rule.
///
/// Returns a new owned Expr pointer that must be freed by the caller.
///
/// # Safety
/// The caller must ensure `rule` is a valid `RewriteRule` pointer.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_rewrite_rule_get_rhs(
    rule: *const RewriteRule
) -> *mut Expr {

    unsafe {

        if rule.is_null() {

            return std::ptr::null_mut(
            );
        }

        Box::into_raw(Box::new(
            (*rule).rhs.clone(),
        ))
    }
}

/// Applies a set of rewrite rules to an expression until a normal form is reached.
///
/// # Safety
/// The caller must ensure `expr` is a valid Expr pointer and `rules` is a valid array.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_apply_rules_to_normal_form(
    expr: *const Expr,
    rules: *const *const RewriteRule,
    rules_len: usize,
) -> *mut Expr {

    unsafe {

        if expr.is_null()
            || (rules_len > 0
                && rules.is_null())
        {

            return std::ptr::null_mut(
            );
        }

        let expr_ref = &*expr;

        // Convert rules array
        let mut rules_vec =
            Vec::with_capacity(
                rules_len,
            );

        if rules_len > 0 {

            let rules_slice =
            std::slice::from_raw_parts(
                rules,
                rules_len,
            );

            for &rule_ptr in rules_slice
            {

                if rule_ptr.is_null() {

                    return std::ptr::null_mut();
                }

                rules_vec.push(
                    (*rule_ptr).clone(),
                );
            }
        }

        let result =
            apply_rules_to_normal_form(
                expr_ref,
                &rules_vec,
            );

        Box::into_raw(Box::new(result))
    }
}

/// Applies the Knuth-Bendix completion algorithm to a set of equations.
///
/// Returns a pointer to a Vec<RewriteRule> on success, or null on failure.
///
/// # Safety
/// The caller must ensure `equations` is a valid array of Expr pointers.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_knuth_bendix(
    equations: *const *const Expr,
    equations_len: usize,
) -> *mut Vec<RewriteRule> {

    unsafe {

        if equations_len > 0
            && equations.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        // Convert equations array
        let mut equations_vec =
            Vec::with_capacity(
                equations_len,
            );

        if equations_len > 0 {

            let equations_slice =
            std::slice::from_raw_parts(
                equations,
                equations_len,
            );

            for &eq_ptr in
                equations_slice
            {

                if eq_ptr.is_null() {

                    return std::ptr::null_mut();
                }

                equations_vec.push(
                    (*eq_ptr).clone(),
                );
            }
        }

        match knuth_bendix(
            &equations_vec,
        ) {
            | Ok(rules) => {
                Box::into_raw(Box::new(
                    rules,
                ))
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    }
}

/// Gets the length of a rules vector.
///
/// # Safety
/// The caller must ensure `rules` is a valid Vec<RewriteRule> pointer.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub const unsafe extern "C" fn rssn_rules_vec_len(
    rules: *const Vec<RewriteRule>
) -> usize {

    unsafe {

        if rules.is_null() {

            0
        } else {

            (*rules).len()
        }
    }
}

/// Gets a rule from a rules vector by index.
///
/// Returns a new owned `RewriteRule` pointer that must be freed by the caller.
///
/// # Safety
/// The caller must ensure `rules` is a valid Vec<RewriteRule> pointer.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_rules_vec_get(
    rules: *const Vec<RewriteRule>,
    index: usize,
) -> *mut RewriteRule {

    unsafe {

        if rules.is_null() {

            return std::ptr::null_mut(
            );
        }

        let rules_ref = &*rules;

        if index >= rules_ref.len() {

            return std::ptr::null_mut(
            );
        }

        Box::into_raw(Box::new(
            rules_ref[index].clone(),
        ))
    }
}

/// Frees a rules vector.
///
/// # Safety
/// The caller must ensure `rules` was created by this module and hasn't been freed yet.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_rules_vec_free(
    rules: *mut Vec<RewriteRule>
) {

    unsafe {

        if !rules.is_null() {

            let _ =
                Box::from_raw(rules);
        }
    }
}

/// Converts a rewrite rule to a string representation.
///
/// The returned string must be freed using `rssn_free_string`.
///
/// # Safety
/// The caller must free the returned string.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_rewrite_rule_to_string(
    rule: *const RewriteRule
) -> *mut c_char {

    unsafe {

        if rule.is_null() {

            return std::ptr::null_mut(
            );
        }

        let rule_ref = &*rule;

        let rule_str = format!(
            "{} -> {}",
            rule_ref.lhs, rule_ref.rhs
        );

        to_c_string(rule_str)
    }
}
