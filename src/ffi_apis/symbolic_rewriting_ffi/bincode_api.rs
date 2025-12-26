//! Bincode-based FFI API for term rewriting systems.

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::rewriting::{
    apply_rules_to_normal_form,
    knuth_bendix,
    RewriteRule,
};
use serde::{
    Deserialize,
    Serialize,
};

#[derive(Serialize, Deserialize)]

struct ApplyRulesInput {
    expr: Expr,
    rules: Vec<RewriteRule>,
}

#[derive(Serialize, Deserialize)]

struct RuleInput {
    lhs: Expr,
    rhs: Expr,
}

#[derive(Serialize, Deserialize)]

struct StringResponse {
    string: String,
}

#[derive(Serialize, Deserialize)]

struct ErrorResponse {
    error: String,
}

/// Applies rewrite rules to an expression (Bincode).
#[no_mangle]

pub extern "C" fn rssn_apply_rules_to_normal_form_bincode(
    input: BincodeBuffer
) -> BincodeBuffer {

    let input_data: Option<
        ApplyRulesInput,
    > = from_bincode_buffer(&input);

    let input_data = match input_data {
        | Some(i) => i,
        | None => {
            return BincodeBuffer::empty(
            )
        },
    };

    let result =
        apply_rules_to_normal_form(
            &input_data.expr,
            &input_data.rules,
        );

    to_bincode_buffer(&result)
}

/// Applies the Knuth-Bendix completion algorithm (Bincode).
#[no_mangle]

pub extern "C" fn rssn_knuth_bendix_bincode(
    input: BincodeBuffer
) -> BincodeBuffer {

    let equations: Option<Vec<Expr>> =
        from_bincode_buffer(&input);

    let equations = match equations {
        | Some(e) => e,
        | None => {
            return BincodeBuffer::empty(
            )
        },
    };

    match knuth_bendix(&equations) {
        | Ok(rules) => {
            to_bincode_buffer(&rules)
        },
        | Err(err) => {

            let error_response =
                ErrorResponse {
                    error: err,
                };

            to_bincode_buffer(
                &error_response,
            )
        },
    }
}

/// Creates a rewrite rule from Bincode.
#[no_mangle]

pub extern "C" fn rssn_rewrite_rule_new_bincode(
    input: BincodeBuffer
) -> BincodeBuffer {

    let input_data: Option<RuleInput> =
        from_bincode_buffer(&input);

    let input_data = match input_data {
        | Some(i) => i,
        | None => {
            return BincodeBuffer::empty(
            )
        },
    };

    let rule = RewriteRule {
        lhs: input_data.lhs,
        rhs: input_data.rhs,
    };

    to_bincode_buffer(&rule)
}

/// Converts a rewrite rule to a human-readable string (Bincode).
#[no_mangle]

pub extern "C" fn rssn_rewrite_rule_to_string_bincode(
    input: BincodeBuffer
) -> BincodeBuffer {

    let rule: Option<RewriteRule> =
        from_bincode_buffer(&input);

    let rule = match rule {
        | Some(r) => r,
        | None => {
            return BincodeBuffer::empty(
            )
        },
    };

    let rule_str = format!(
        "{} -> {}",
        rule.lhs, rule.rhs
    );

    let response = StringResponse {
        string: rule_str,
    };

    to_bincode_buffer(&response)
}
