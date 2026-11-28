//! Configurations and `rssn_compute`.

use std::ffi::CString;
use std::mem::MaybeUninit;
use std::os::raw::c_char;
use std::ptr;

use serde::Deserialize;
use serde::Serialize;

use super::RSSN_NO_TERM;
use super::RssnStatus;
use super::compute_error;
use super::fail;
use super::guard;
use super::term::RssnSession;
use super::term::borrow_session;
use super::term::resolve;
use super::text;
use super::to_c;
use crate::api::Answer;
use crate::api::Config;
use crate::graph::DimensionMode;
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

/// Per-request metadata returned by `rssn_compute_all_with_metadata`.
///
/// Each non-null pointer is an owned string. Release it with `rssn_string_free`.
#[repr(C)]
#[derive(Debug)]
pub struct RssnAnswerMetadata {
    /// Inferred dimension of a successful answer, or null if unavailable.
    pub dimension: *mut c_char,
    /// Stable diagnostic codes separated by line feeds, or null if empty.
    pub diagnostic_codes: *mut c_char,
    /// Error message for this request, or null on success.
    pub message: *mut c_char,
}

#[derive(Default)]
struct OwnedAnswerMetadata {
    dimension: Option<CString>,
    diagnostic_codes: Option<CString>,
    message: Option<CString>,
}

impl OwnedAnswerMetadata {
    fn computed(answer: &Answer<'_>) -> Self {
        Self {
            dimension: answer
                .dimension
                .as_ref()
                .and_then(|dimension| CString::new(dimension.to_string()).ok()),
            diagnostic_codes: (!answer.diagnostics.is_empty())
                .then(|| {
                    answer
                        .diagnostics
                        .iter()
                        .map(crate::graph::DimError::code)
                        .collect::<Vec<_>>()
                        .join("\n")
                })
                .and_then(|codes| CString::new(codes).ok()),
            message: None,
        }
    }

    fn failed(
        message: &str,
        code: Option<&str>,
    ) -> Self {
        Self {
            dimension: None,
            diagnostic_codes: code.and_then(|code| CString::new(code).ok()),
            message: CString::new(message.replace('\0', " ")).ok(),
        }
    }

    fn into_ffi(self) -> RssnAnswerMetadata {
        RssnAnswerMetadata {
            dimension: self.dimension.map_or(ptr::null_mut(), CString::into_raw),
            diagnostic_codes: self
                .diagnostic_codes
                .map_or(ptr::null_mut(), CString::into_raw),
            message: self.message.map_or(ptr::null_mut(), CString::into_raw),
        }
    }
}

fn raw_answer(answer: &Answer<'_>) -> RssnAnswer {
    RssnAnswer {
        term: answer.term.id(),
        has_value: answer.value.is_some(),
        value: answer.value.unwrap_or(f64::NAN),
        error: answer.error.unwrap_or(f64::NAN),
        reduced: answer.reduced,
        iterations: answer.report.iterations,
    }
}

/// The answer of one computation, as the record both `rssn_compute_json` and
/// [`rssn_compute_bincode`] carry.
///
/// The JSON text and the `bincode-next` bytes are two encodings of this one
/// struct, so the two paths cannot drift. The field names are the JSON keys;
/// the byte form is
/// [`encode_binary`](crate::kernel::codec::encode_binary) — the
/// same configuration plugin ABI messages and content-address digests use —
/// and
/// [`decode_binary`](crate::kernel::codec::decode_binary) reads
/// it back.
#[derive(Clone, Debug, Serialize, Deserialize, bincode_next::Encode, bincode_next::Decode)]
pub struct RssnComputeRecord {
    /// The answer term as infix text: what `rssn_term_to_string` returns for
    /// it.
    pub term: String,
    /// The numeric value, when one is known.
    pub value: Option<f64>,
    /// Estimated absolute error of `value`, when one is known.
    pub error: Option<f64>,
    /// Whether the request was fully reduced; `false` means the answer still
    /// contains unevaluated operators.
    pub reduced: bool,
    /// Outer iterations the engine performed.
    ///
    /// `u64` rather than `usize` so the encoding is the same on 32- and
    /// 64-bit hosts.
    pub iterations: u64,
    /// The inferred dimension as text, when one was inferred.
    pub dimension: Option<String>,
    /// Advisory dimensional findings, as text.
    pub diagnostics: Vec<String>,
    /// Stable machine-readable kinds for findings.
    pub diagnostic_codes: Vec<String>,
}

/// The record both wire formats carry.
fn record(answer: &Answer<'_>) -> RssnComputeRecord {
    super::remember_dimensions(&answer.diagnostics);
    RssnComputeRecord {
        term: answer.term.to_string(),
        value: answer.value,
        error: answer.error,
        reduced: answer.reduced,
        iterations: u64::try_from(answer.report.iterations).unwrap_or(u64::MAX),
        dimension: answer
            .dimension
            .as_ref()
            .map(std::string::ToString::to_string),
        diagnostics: answer
            .diagnostics
            .iter()
            .map(std::string::ToString::to_string)
            .collect(),
        diagnostic_codes: answer
            .diagnostics
            .iter()
            .map(|d| d.code().to_string())
            .collect(),
    }
}

/// One request's entry of a [`RssnComputeBatchRecord`], as
/// `rssn_compute_all_bincode` writes it.
#[derive(Clone, Debug, Serialize, Deserialize, bincode_next::Encode, bincode_next::Decode)]
pub struct RssnBatchEntry {
    /// The status of this request, as the numeric value of [`RssnStatus`]:
    /// the same code a C host reads from the `statuses` array of
    /// `rssn_compute_all`.
    pub status: u32,
    /// The answer, present exactly when `status` is `RssnStatus::Ok`.
    pub answer: Option<RssnComputeRecord>,
}

/// The answer of `rssn_compute_all_bincode`: one [`RssnBatchEntry`] per input
/// handle, in input order.
#[derive(Clone, Debug, Serialize, Deserialize, bincode_next::Encode, bincode_next::Decode)]
pub struct RssnComputeBatchRecord {
    /// The entries, in input order: entry `k` belongs to `terms[k]`.
    pub entries: Vec<RssnBatchEntry>,
}

/// Runs one request under `config`; null means the defaults.
///
/// # Safety
/// `session` must be a live session and `config` null or a live
/// configuration.
unsafe fn run_request<'s>(
    session: *const RssnSession,
    term: u32,
    config: *const RssnConfig,
) -> Result<Answer<'s>, RssnStatus> {
    // SAFETY: caller contract.
    let request = match unsafe { borrow_session(session) }.and_then(|s| resolve(s, term)) {
        | Ok(t) => t,
        | Err(status) => return Err(status),
    };
    let default;
    // SAFETY: caller contract.
    let config = match unsafe { config.as_ref() } {
        | Some(config) => &config.0,
        | None => {
            default = Config::new();
            &default
        },
    };
    let result = request.session().compute(request, config);
    // The dimension codes are part of the answer a C host may need without
    // parsing JSON: record them here, where every compute entry point meets.
    match &result {
        | Ok(answer) => super::remember_dimensions(&answer.diagnostics),
        | Err(crate::api::ComputeError::Dimension(error)) => {
            super::remember_dimension_codes([error.code()]);
        },
        | Err(_) => super::clear_dimensions(),
    }
    result.map_err(|error| compute_error(&error))
}

/// A configuration with every standard rule set and a symbolic target.
/// Free it with `rssn_config_free`. Returns null if the built-in rule sets
/// fail to install; the reason is then in `rssn_last_error`.
#[unsafe(no_mangle)]
pub extern "C" fn rssn_config_new() -> *mut RssnConfig {
    let mut config = std::ptr::null_mut();
    let _ = guard(|| {
        config = Box::into_raw(Box::new(RssnConfig(Config::new())));
        RssnStatus::Ok
    });
    config
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
    edit(config, |c| {
        Ok(c.assume(unsafe { text(name) }?, Facts::from_bits(facts)))
    })
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

/// Sets the dimension checking mode: 0 off, 1 advisory (the default), 2
/// strict. Any other value is rejected.
///
/// # Safety
/// `config` must be a live configuration.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_config_dimension_mode(
    config: *mut RssnConfig,
    mode: u8,
) -> RssnStatus {
    edit(config, |c| {
        let mode = match mode {
            | 0 => DimensionMode::Off,
            | 1 => DimensionMode::Advisory,
            | 2 => DimensionMode::Strict,
            | _ => {
                return Err(fail(
                    RssnStatus::InvalidArgument,
                    "dimension mode must be 0, 1 or 2",
                ));
            },
        };
        Ok(c.dimensions(mode))
    })
}

/// Declares the dimension of `symbol` as `text` (a unit expression, a
/// base-dimension product, or `"any"`).
///
/// # Safety
/// `config` must be a live configuration; both strings NUL-terminated.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_config_dimension(
    config: *mut RssnConfig,
    symbol: *const c_char,
    text: *const c_char,
) -> RssnStatus {
    // SAFETY: caller contract.
    edit(config, |c| {
        Ok(c.dimension(unsafe { super::text(symbol) }?, unsafe {
            super::text(text)
        }?))
    })
}

/// Registers a new base dimension `name` (for example `"Currency"`).
///
/// # Safety
/// `config` must be a live configuration; `name` NUL-terminated.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_config_define_base_dimension(
    config: *mut RssnConfig,
    name: *const c_char,
) -> RssnStatus {
    // SAFETY: caller contract.
    edit(config, |c| {
        Ok(c.define_base_dimension(unsafe { text(name) }?))
    })
}

/// Disables hard dimensional errors even in strict mode (the explicit
/// escape hatch). Advisory diagnostics are still reported.
///
/// # Safety
/// `config` must be a live configuration.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_config_bypass_dimensions(config: *mut RssnConfig) -> RssnStatus {
    edit(config, |c| Ok(c.bypass_dimension_checks()))
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
        let answer = match unsafe { run_request(session, term, config) } {
            | Ok(answer) => answer,
            | Err(status) => return status,
        };
        // SAFETY: non-null and valid by contract.
        unsafe { *out = raw_answer(&answer) };
        RssnStatus::Ok
    })
}

/// Computes `term` and writes a JSON object to `*out`.
///
/// The object is [`RssnComputeRecord`] serialised as JSON: `term`, `value`,
/// `error`, `reduced`, `iterations`, `dimension` (a dimension string or
/// `null`) and `diagnostics` (an array of strings). Free the string with
/// `rssn_string_free`. [`rssn_compute_bincode`] carries the same record in
/// the crate's binary wire format.
///
/// # Safety
/// `session` must be a live session, `config` null or a live
/// configuration, and `out` a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_compute_json(
    session: *const RssnSession,
    term: u32,
    config: *const RssnConfig,
    out: *mut *mut c_char,
) -> RssnStatus {
    guard(|| {
        if out.is_null() {
            return fail(RssnStatus::NullArgument, "null output pointer");
        }
        // SAFETY: caller contract.
        let answer = match unsafe { run_request(session, term, config) } {
            | Ok(answer) => answer,
            | Err(status) => return status,
        };
        let document = match serde_json::to_value(record(&answer)) {
            | Ok(document) => document,
            | Err(error) => {
                return fail(
                    RssnStatus::Panic,
                    &format!("the answer could not be encoded: {error}"),
                );
            },
        };
        let text = to_c(document.to_string());
        if text.is_null() {
            return fail(
                RssnStatus::Panic,
                "the answer is not representable as a C string",
            );
        }
        // SAFETY: `out` is non-null by the check above.
        unsafe { *out = text };
        RssnStatus::Ok
    })
}

/// Computes `term` and writes its inferred dimension as text (`""` when no
/// dimension was inferred or dimensional checking is off). Free the string
/// with `rssn_string_free`.
///
/// # Safety
/// `session` must be a live session, `config` null or a live
/// configuration, and `out` a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_dimension(
    session: *const RssnSession,
    term: u32,
    config: *const RssnConfig,
    out: *mut *mut c_char,
) -> RssnStatus {
    guard(|| {
        if out.is_null() {
            return fail(RssnStatus::NullArgument, "null output pointer");
        }
        // SAFETY: caller contract.
        let answer = match unsafe { run_request(session, term, config) } {
            | Ok(answer) => answer,
            | Err(status) => return status,
        };
        let text = answer
            .dimension
            .map_or_else(String::new, |dimension| dimension.to_string());
        let text = to_c(text);
        if text.is_null() {
            return fail(
                RssnStatus::Panic,
                "the dimension is not representable as a C string",
            );
        }
        // SAFETY: `out` is non-null by the check above.
        unsafe { *out = text };
        RssnStatus::Ok
    })
}

/// Computes `term` and writes the answer as a `bincode-next` byte string.
///
/// The bytes are [`RssnComputeRecord`] — the same record
/// [`rssn_compute_json`] encodes as JSON — in the deterministic
/// `bincode-next` format of
/// [`encode_binary`](crate::kernel::codec::encode_binary), which
/// is also the format of the kernel manifest and the plugin ABI. `*out_len`
/// receives the byte count and `*out` the buffer; release it with
/// [`rssn_bytes_free`](super::rssn_bytes_free). Neither pointer is written
/// when the call fails.
///
/// Decode the buffer using the field order and types of [`RssnComputeRecord`].
///
/// # Safety
/// `session` must be a live session, `config` null or a live configuration,
/// and `out` and `out_len` valid pointers.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_compute_bincode(
    session: *const RssnSession,
    term: u32,
    config: *const RssnConfig,
    out: *mut *mut u8,
    out_len: *mut usize,
) -> RssnStatus {
    guard(|| {
        if out.is_null() || out_len.is_null() {
            return fail(RssnStatus::NullArgument, "null output pointer");
        }
        // SAFETY: caller contract.
        let answer = match unsafe { run_request(session, term, config) } {
            | Ok(answer) => answer,
            | Err(status) => return status,
        };
        let bytes = match crate::kernel::codec::encode_binary(&record(&answer)) {
            | Ok(bytes) => bytes,
            | Err(error) => {
                return fail(
                    RssnStatus::Panic,
                    &format!("the answer could not be encoded: {error}"),
                );
            },
        };
        let len = bytes.len();
        let data = Box::into_raw(bytes.into_boxed_slice()).cast::<u8>();
        // SAFETY: both pointers are non-null by the check above and valid by
        // the caller's contract.
        unsafe {
            *out = data;
            *out_len = len;
        }
        RssnStatus::Ok
    })
}

/// Computes `n` requests at once (in parallel inside the library) with one configuration.
///
/// `terms` holds the `n` request handles; `answers` and `statuses` are caller-provided arrays of
/// `n` entries that receive, for each request, its `RssnAnswer` and its `RssnStatus` (an entry
/// whose status is not `Ok` has `term == RSSN_NO_TERM`). The answers equal `n` calls of
/// `rssn_compute`. The return value is `Ok` when the call itself was well formed, even if some
/// requests failed; the message of the first failing request is kept for `rssn_last_error`.
/// `config` may be null for the defaults.
///
/// # Safety
/// `session` must be a live session, `config` null or a live configuration, and `terms`,
/// `answers`, `statuses` valid for `n` elements (they may be null only when `n` is 0).
/// The arrays must not overlap. Output arrays may be uninitialized and all
/// their elements are initialized on an `Ok` return.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_compute_all(
    session: *const RssnSession,
    terms: *const u32,
    n: usize,
    config: *const RssnConfig,
    answers: *mut RssnAnswer,
    statuses: *mut RssnStatus,
) -> RssnStatus {
    // SAFETY: the caller supplies the arrays required by compute_batch.
    unsafe { compute_batch(session, terms, n, config, answers, statuses, None) }
}

/// Computes requests in parallel and writes one metadata record per request.
///
/// The answer and status arrays follow the same contract as `rssn_compute_all`.
/// Metadata stays in request order, including when some requests fail. Each
/// non-null string belongs to the caller and is released with `rssn_string_free`.
///
/// # Safety
/// The session, configuration and arrays must satisfy `rssn_compute_all`'s
/// contract. The metadata array must be writable for `n` elements and disjoint
/// from the other arrays. It may be null when `n` is zero. Every element is
/// initialized on an `Ok` return.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_compute_all_with_metadata(
    session: *const RssnSession,
    terms: *const u32,
    n: usize,
    config: *const RssnConfig,
    answers: *mut RssnAnswer,
    statuses: *mut RssnStatus,
    metadata: *mut RssnAnswerMetadata,
) -> RssnStatus {
    // SAFETY: the caller supplies the arrays required by compute_batch.
    unsafe { compute_batch(session, terms, n, config, answers, statuses, Some(metadata)) }
}

#[allow(
    clippy::too_many_lines,
    reason = "keeps handle resolution, parallel dispatch and ordered writes together"
)]
unsafe fn compute_batch(
    session: *const RssnSession,
    terms: *const u32,
    n: usize,
    config: *const RssnConfig,
    answers: *mut RssnAnswer,
    statuses: *mut RssnStatus,
    metadata: Option<*mut RssnAnswerMetadata>,
) -> RssnStatus {
    guard(|| {
        if n > 0
            && (terms.is_null()
                || answers.is_null()
                || statuses.is_null()
                || metadata.is_some_and(<*mut RssnAnswerMetadata>::is_null))
        {
            return fail(RssnStatus::NullArgument, "null array argument");
        }
        // SAFETY: caller contract.
        let session = match unsafe { borrow_session(session) } {
            | Ok(s) => s,
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
        let empty = RssnAnswer {
            term: RSSN_NO_TERM,
            has_value: false,
            value: f64::NAN,
            error: f64::NAN,
            reduced: false,
            iterations: 0,
        };
        // SAFETY: the caller provides separate arrays valid for `n` elements.
        // MaybeUninit preserves the outputs' layout and allows unwritten contents.
        let (ids, answers, statuses) = unsafe {
            if n == 0 {
                (&[][..], &mut [][..], &mut [][..])
            } else {
                (
                    std::slice::from_raw_parts(terms, n),
                    std::slice::from_raw_parts_mut(answers.cast::<MaybeUninit<RssnAnswer>>(), n),
                    std::slice::from_raw_parts_mut(statuses.cast::<MaybeUninit<RssnStatus>>(), n),
                )
            }
        };
        let mut live = Vec::new();
        let mut slots = Vec::new();
        let mut reports = metadata.map(|_| {
            let mut reports = Vec::with_capacity(n);
            reports.resize_with(n, OwnedAnswerMetadata::default);
            reports
        });
        // `resolve` already records the missing-handle message.
        let mut first_missing = n;
        for (index, ((answer, status), &id)) in
            answers.iter_mut().zip(&mut *statuses).zip(ids).enumerate()
        {
            answer.write(empty);
            let current = match resolve(session, id) {
                | Ok(term) => {
                    live.push(term);
                    slots.push(index);
                    // An internal failure to return a result still leaves the slot initialized.
                    RssnStatus::Panic
                },
                | Err(code) => {
                    first_missing = first_missing.min(index);
                    if let Some(report) = reports.as_mut().and_then(|rows| rows.get_mut(index)) {
                        *report = OwnedAnswerMetadata::failed(&super::last_error_text(), None);
                    }
                    code
                },
            };
            status.write(current);
        }
        let mut first_error: Option<String> = None;
        for (&k, result) in slots.iter().zip(session.compute_all(&live, config)) {
            let report = reports.as_mut().and_then(|rows| rows.get_mut(k));
            match result {
                | Ok(answer) => {
                    super::remember_dimensions(&answer.diagnostics);
                    statuses[k].write(RssnStatus::Ok);
                    answers[k].write(raw_answer(&answer));
                    if let Some(report) = report {
                        *report = OwnedAnswerMetadata::computed(&answer);
                    }
                },
                | Err(e) => {
                    let dimension_code = if let crate::api::ComputeError::Dimension(error) = &e {
                        super::remember_dimension_codes([error.code()]);
                        Some(error.code())
                    } else {
                        super::clear_dimensions();
                        None
                    };
                    statuses[k].write(super::status_of(&e));
                    if let Some(report) = report {
                        let message = e.to_string();
                        *report = OwnedAnswerMetadata::failed(&message, dimension_code);
                        if k < first_missing {
                            first_error.get_or_insert(message);
                        }
                    } else if k < first_missing {
                        first_error.get_or_insert_with(|| e.to_string());
                    }
                },
            }
        }
        if let Some(message) = first_error {
            let _ = fail(RssnStatus::Ok, &message);
        }
        if let (Some(metadata), Some(reports)) = (metadata, reports) {
            for (index, report) in reports.into_iter().enumerate() {
                // SAFETY: the caller provides n writable metadata slots. No
                // allocations or fallible work remain after taking ownership
                // of the strings from the Rust records.
                unsafe { metadata.add(index).write(report.into_ffi()) };
            }
        }
        RssnStatus::Ok
    })
}

/// Computes `n` requests at once and writes the batch answer as a
/// `bincode-next` byte string.
///
/// The bytes are [`RssnComputeBatchRecord`]: one [`RssnBatchEntry`] per input
/// handle, in input order, each carrying the [`RssnComputeRecord`] that
/// [`rssn_compute_bincode`] would write for a request that succeeded. As with
/// [`rssn_compute_all`], the return value is `Ok` when the call itself was
/// well formed even if some requests failed: the per-request [`RssnStatus`] is
/// inside the record, and the message of the first failing request is kept for
/// `rssn_last_error`. `config` may be null for the defaults.
///
/// `*out_len` receives the byte count and `*out` the buffer; release it with
/// [`rssn_bytes_free`](super::rssn_bytes_free). Neither output pointer is
/// written when the call fails.
///
/// The bytes are opaque to C: a host decodes them from a schema generated from
/// [`RssnComputeBatchRecord`].
///
/// # Safety
/// `session` must be a live session, `config` null or a live configuration,
/// and `terms` valid for `n` handles (it may be null only when `n` is 0).
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_compute_all_bincode(
    session: *const RssnSession,
    terms: *const u32,
    n: usize,
    config: *const RssnConfig,
    out: *mut *mut u8,
    out_len: *mut usize,
) -> RssnStatus {
    guard(|| {
        if out.is_null() || out_len.is_null() {
            return fail(RssnStatus::NullArgument, "null output pointer");
        }
        if n > 0 && terms.is_null() {
            return fail(RssnStatus::NullArgument, "null array argument");
        }
        // SAFETY: caller contract.
        let session = match unsafe { borrow_session(session) } {
            | Ok(session) => session,
            | Err(status) => return status,
        };
        let default;
        // SAFETY: caller contract.
        let config = match unsafe { config.as_ref() } {
            | Some(config) => &config.0,
            | None => {
                default = Config::new();
                &default
            },
        };
        // SAFETY: `terms` holds `n` handles by contract (null only when 0).
        let ids: &[u32] = unsafe {
            if n == 0 {
                &[]
            } else {
                std::slice::from_raw_parts(terms, n)
            }
        };
        // In input order. `None` marks a handle that resolved, whose result
        // therefore comes from `compute_all` below; `Some` is the status of a
        // handle that did not resolve (its message is already in
        // `rssn_last_error`).
        let mut statuses: Vec<Option<RssnStatus>> = vec![None; n];
        let mut live = Vec::new();
        let mut first_missing = n;
        for (index, (status, &id)) in statuses.iter_mut().zip(ids).enumerate() {
            match resolve(session, id) {
                | Ok(term) => live.push(term),
                | Err(error) => {
                    first_missing = first_missing.min(index);
                    *status = Some(error);
                },
            }
        }
        // `compute_all` answers the resolved handles in their order, which is
        // the order of the `None` slots, so one pass over `statuses` zips the
        // results back into input order.
        let mut results = session.compute_all(&live, config).into_iter();
        let mut first_error: Option<String> = None;
        let mut entries: Vec<RssnBatchEntry> = Vec::with_capacity(n);
        for (index, status) in statuses.into_iter().enumerate() {
            let (status, answer) = match status {
                | Some(status) => (status, None),
                | None => {
                    match results.next() {
                        | Some(Ok(answer)) => (RssnStatus::Ok, Some(record(&answer))),
                        | Some(Err(error)) => {
                            if let crate::api::ComputeError::Dimension(dimension) = &error {
                                super::remember_dimension_codes([dimension.code()]);
                            } else {
                                super::clear_dimensions();
                            }
                            if index < first_missing {
                                first_error.get_or_insert_with(|| error.to_string());
                            }
                            (super::status_of(&error), None)
                        },
                        // Every resolved handle has a result. This cannot happen;
                        // a wrong status is still cheaper than unwinding a C
                        // caller.
                        | None => (RssnStatus::Panic, None),
                    }
                },
            };
            entries.push(RssnBatchEntry {
                status: status as u32,
                answer,
            });
        }
        if let Some(message) = first_error {
            let _ = fail(RssnStatus::Ok, &message);
        }
        let batch = RssnComputeBatchRecord { entries };
        let bytes = match crate::kernel::codec::encode_binary(&batch) {
            | Ok(bytes) => bytes,
            | Err(error) => {
                return fail(
                    RssnStatus::Panic,
                    &format!("the batch could not be encoded: {error}"),
                );
            },
        };
        let len = bytes.len();
        let data = Box::into_raw(bytes.into_boxed_slice()).cast::<u8>();
        // SAFETY: both pointers are non-null by the check above and valid by
        // the caller's contract.
        unsafe {
            *out = data;
            *out_len = len;
        }
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
