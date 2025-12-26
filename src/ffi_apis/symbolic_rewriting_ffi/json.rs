//! JSON-based FFI API for term rewriting systems.

use std::os::raw::c_char;

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::rewriting::apply_rules_to_normal_form;
use crate::symbolic::rewriting::knuth_bendix;
use crate::symbolic::rewriting::RewriteRule;

#[derive(Serialize, Deserialize)]

struct ApplyRulesInput {
    expr: Expr,
    rules: Vec<RewriteRule>,
}

/// Applies rewrite rules to an expression (JSON).
///
/// Input: JSON object with "expr" and "rules" fields
/// Output: JSON-serialized Expr (the normal form)
#[no_mangle]

pub extern "C" fn rssn_apply_rules_to_normal_form_json(
    json_str: *const c_char
) -> *mut c_char {

    let input: Option<ApplyRulesInput> =
        from_json_string(json_str);

    let input = match input {
        | Some(i) => i,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let result =
        apply_rules_to_normal_form(
            &input.expr,
            &input.rules,
        );

    to_json_string(&result)
}

/// Applies the Knuth-Bendix completion algorithm (JSON).
///
/// Input: JSON array of equations (Expr::Eq)
/// Output: JSON array of RewriteRule objects
#[no_mangle]

pub extern "C" fn rssn_knuth_bendix_json(
    json_str: *const c_char
) -> *mut c_char {

    let equations: Option<Vec<Expr>> =
        from_json_string(json_str);

    let equations = match equations {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut()
        },
    };

    match knuth_bendix(&equations) {
        | Ok(rules) => {
            to_json_string(&rules)
        },
        | Err(err) => {

            let error_response = serde_json::json!({ "error": err });

            to_json_string(
                &error_response,
            )
        },
    }
}

/// Creates a rewrite rule from JSON.
///
/// Input: JSON object with "lhs" and "rhs" fields (both Expr)
/// Output: JSON-serialized RewriteRule
#[no_mangle]

pub extern "C" fn rssn_rewrite_rule_new_json(
    json_str: *const c_char
) -> *mut c_char {

    #[derive(Deserialize)]

    struct RuleInput {
        lhs: Expr,
        rhs: Expr,
    }

    let input: Option<RuleInput> =
        from_json_string(json_str);

    let input = match input {
        | Some(i) => i,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let rule = RewriteRule {
        lhs: input.lhs,
        rhs: input.rhs,
    };

    to_json_string(&rule)
}

/// Converts a rewrite rule to a human-readable string (JSON).
///
/// Input: JSON-serialized RewriteRule
/// Output: JSON object with "string" field
#[no_mangle]

pub extern "C" fn rssn_rewrite_rule_to_string_json(
    json_str: *const c_char
) -> *mut c_char {

    let rule: Option<RewriteRule> =
        from_json_string(json_str);

    let rule = match rule {
        | Some(r) => r,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let rule_str = format!(
        "{} -> {}",
        rule.lhs, rule.rhs
    );

    let response = serde_json::json!({ "string": rule_str });

    to_json_string(&response)
}
