//! Configurations and `rssn_compute`.

use std::os::raw::c_char;

use super::RSSN_NO_TERM;
use super::RssnStatus;
use super::compute_error;
use super::fail;
use super::guard;
use super::term::RssnSession;
use super::term::borrow_session;
use super::term::resolve;
use super::text;
use crate::api::Config;
use crate::graph::Facts;

/// What to compute with and what to aim for. Opaque.
pub struct RssnConfig(pub(crate) Config);

/// The result of `rssn_compute`.
#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct RssnAnswer {
    /// The answer term: a closed form, or a float literal for a numeric
    /// target.
    pub term: u32,
    /// Whether `value` and `error` are meaningful.
    pub has_value: bool,
    /// The numeric value, when known.
    pub value: f64,
    /// Estimated absolute error of `value`.
    pub error: f64,
    /// Whether the request was fully reduced; `false` means the term still
    /// contains unevaluated operators.
    pub reduced: bool,
    /// Outer iterations the engine performed.
    pub iterations: usize,
}

/// A configuration with every standard rule set and a symbolic target.
/// Free it with `rssn_config_free`.
#[unsafe(no_mangle)]
pub extern "C" fn rssn_config_new() -> *mut RssnConfig {
    std::panic::catch_unwind(|| Box::into_raw(Box::new(RssnConfig(Config::new()))))
        .unwrap_or(std::ptr::null_mut())
}

/// Releases a configuration. Null is ignored.
///
/// # Safety
/// `config` must be null or come from `rssn_config_new`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_config_free(config: *mut RssnConfig) {
    if !config.is_null() {
        // SAFETY: allocated by `rssn_config_new`.
        drop(unsafe { Box::from_raw(config) });
    }
}

/// Applies `change` to the configuration behind `config`.
fn edit(
    config: *mut RssnConfig,
    change: impl FnOnce(Config) -> Result<Config, RssnStatus>,
) -> RssnStatus {
    guard(|| {
        // SAFETY: the caller passes null or a live configuration.
        let Some(slot) = (unsafe { config.as_mut() }) else {
            return fail(RssnStatus::NullArgument, "null configuration");
        };
        match change(slot.0.clone()) {
            | Ok(next) => {
                slot.0 = next;
                RssnStatus::Ok
            },
            | Err(e) => e,
        }
    })
}

/// Asks for a closed-form answer (the default).
///
/// # Safety
/// `config` must be a live configuration.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_config_symbolic(config: *mut RssnConfig) -> RssnStatus {
    edit(config, |c| Ok(c.symbolic()))
}

/// Asks for a numeric answer within the absolute `tolerance`.
///
/// # Safety
/// `config` must be a live configuration.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_config_numeric(
    config: *mut RssnConfig,
    tolerance: f64,
) -> RssnStatus {
    edit(config, |c| Ok(c.numeric(tolerance)))
}

/// Binds the symbol `name` to `value` for numeric evaluation.
///
/// # Safety
/// `config` must be a live configuration; `name` a NUL-terminated string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_config_bind(
    config: *mut RssnConfig,
    name: *const c_char,
    value: f64,
) -> RssnStatus {
    // SAFETY: caller contract.
    edit(config, |c| Ok(c.bind(unsafe { text(name) }?, value)))
}

/// Declares facts about the symbol `name`. `facts` is a bit set: 1 real,
/// 2 nonzero, 4 nonnegative, 8 nonpositive, 16 integer (so 7 is
/// "positive").
///
/// # Safety
/// `config` must be a live configuration; `name` a NUL-terminated string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_config_assume(
    config: *mut RssnConfig,
    name: *const c_char,
    facts: u8,
) -> RssnStatus {
    // SAFETY: caller contract.
    edit(config, |c| Ok(c.assume(unsafe { text(name) }?, Facts::from_bits(facts))))
}

/// Limits the search: at most `max_iterations` outer iterations and
/// `max_nodes` graph nodes. Zero keeps the current value.
///
/// # Safety
/// `config` must be a live configuration.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_config_budget(
    config: *mut RssnConfig,
    max_iterations: usize,
    max_nodes: usize,
) -> RssnStatus {
    edit(config, |c| {
        let mut budget = c.current_budget();
        if max_iterations > 0 {
            budget.max_iterations = max_iterations;
        }
        if max_nodes > 0 {
            budget.max_nodes = max_nodes;
        }
        Ok(c.budget(budget))
    })
}

/// Reduces `term` to the target of `config` (null means the defaults) and
/// stores the result in `*out`.
///
/// # Safety
/// `session` must be a live session, `config` null or a live
/// configuration, and `out` a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_compute(
    session: *const RssnSession,
    term: u32,
    config: *const RssnConfig,
    out: *mut RssnAnswer,
) -> RssnStatus {
    guard(|| {
        if out.is_null() {
            return fail(RssnStatus::NullArgument, "null output pointer");
        }
        // SAFETY: caller contract.
        let request = match unsafe { borrow_session(session) }.and_then(|s| resolve(s, term)) {
            | Ok(t) => t,
            | Err(e) => return e,
        };
        let default;
        // SAFETY: caller contract.
        let config = match unsafe { config.as_ref() } {
            | Some(c) => &c.0,
            | None => {
                default = Config::new();
                &default
            },
        };
        let answer = match request.session().compute(request, config) {
            | Ok(a) => a,
            | Err(e) => return compute_error(&e),
        };
        let result = RssnAnswer {
            term: answer.term.id(),
            has_value: answer.value.is_some(),
            value: answer.value.unwrap_or(f64::NAN),
            error: answer.error.unwrap_or(f64::NAN),
            reduced: answer.reduced,
            iterations: answer.report.iterations,
        };
        // SAFETY: non-null and valid by contract.
        unsafe { *out = result };
        RssnStatus::Ok
    })
}

/// Shorthand: parses `source`, computes it with the defaults and returns
/// the answer term, or `RSSN_NO_TERM` (see `rssn_last_error`).
///
/// # Safety
/// `session` must be a live session; `source` a NUL-terminated string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_simplify(
    session: *const RssnSession,
    source: *const c_char,
) -> u32 {
    // SAFETY: forwarded caller contract.
    let request = unsafe { super::rssn_term_parse(session, source) };
    if request == RSSN_NO_TERM {
        return RSSN_NO_TERM;
    }
    let mut answer = RssnAnswer {
        term: RSSN_NO_TERM,
        has_value: false,
        value: f64::NAN,
        error: f64::NAN,
        reduced: false,
        iterations: 0,
    };
    // SAFETY: forwarded caller contract; `answer` is a local.
    match unsafe { rssn_compute(session, request, std::ptr::null(), &raw mut answer) } {
        | RssnStatus::Ok => answer.term,
        | _ => RSSN_NO_TERM,
    }
}
