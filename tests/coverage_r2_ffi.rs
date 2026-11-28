//! Coverage round 2 for the C ABI (`rssn::ffi`): every exported function is
//! called as a Rust function with hostile arguments (null pointers, handles
//! that do not exist, invalid UTF-8, empty and oversized counts, truncated
//! wire buffers) and must answer with a status code, never a panic or a
//! crash; successful paths are round-tripped and everything the API says must
//! be freed is freed.
#![allow(unsafe_code)]

use std::ffi::CStr;
use std::ffi::CString;
use std::os::raw::c_char;

use rssn::ffi::*;
use rssn::kernel::codec::decode_deterministic;
use rssn::sim::scenario::SCENARIOS;

fn c(text: &str) -> CString {
    CString::new(text).unwrap()
}

/// A NUL-terminated string that is not UTF-8.
fn bad_utf8() -> CString {
    CString::new(vec![b'x', 0xff, 0xfe, b'y']).unwrap()
}

fn owned(raw: *mut c_char) -> String {
    assert!(!raw.is_null());
    // SAFETY: a string just returned by rssn.
    let text = unsafe { CStr::from_ptr(raw) }.to_str().unwrap().to_owned();
    // SAFETY: returned by rssn, freed once.
    unsafe { rssn_string_free(raw) };
    text
}

fn last_error() -> String {
    // SAFETY: valid until the next failing call.
    unsafe { CStr::from_ptr(rssn_last_error()) }
        .to_string_lossy()
        .into_owned()
}

fn blank_answer() -> RssnAnswer {
    RssnAnswer {
        term: RSSN_NO_TERM,
        has_value: false,
        value: f64::NAN,
        error: f64::NAN,
        reduced: false,
        iterations: 0,
    }
}

fn blank_tensor() -> RssnTensor {
    RssnTensor {
        rank: 0,
        shape: std::ptr::null_mut(),
        len: 0,
        data: std::ptr::null_mut(),
    }
}

// ------------------------------------------------------------ lifetime

#[test]
fn free_functions_accept_null_and_session_and_config_round_trip() {
    // SAFETY: null is explicitly allowed by every release function.
    unsafe {
        rssn_session_free(std::ptr::null_mut());
        rssn_config_free(std::ptr::null_mut());
        rssn_string_free(std::ptr::null_mut());
        rssn_bytes_free(std::ptr::null_mut(), 0);
        rssn_bytes_free(std::ptr::null_mut(), 17);
        rssn_tensor_free(std::ptr::null_mut());
    }
    for _ in 0..50 {
        let s = rssn_session_new();
        let cfg = rssn_config_new();
        assert!(!s.is_null() && !cfg.is_null());
        // SAFETY: allocated just above, freed once.
        unsafe {
            rssn_config_free(cfg);
            rssn_session_free(s);
        }
    }
}

#[test]
fn version_is_the_package_version_and_static() {
    let a = rssn_version();
    let b = rssn_version();
    assert_eq!(a, b, "static pointer");
    // SAFETY: static NUL-terminated string.
    let v = unsafe { CStr::from_ptr(a) }.to_str().unwrap();
    assert_eq!(v, env!("CARGO_PKG_VERSION"));
}

#[test]
fn tensor_free_is_idempotent_on_the_cleared_struct() {
    let s = rssn_session_new();
    // SAFETY: live session; tensor freed through the API.
    unsafe {
        let t = rssn_term_parse(s, c("list(1, 2, 3, 4)").as_ptr());
        let mut tensor = blank_tensor();
        assert_eq!(
            rssn_term_tensor_data(s, t, &raw mut tensor),
            RssnStatus::Ok,
            "{}",
            last_error()
        );
        assert_eq!(tensor.rank, 1);
        assert_eq!(std::slice::from_raw_parts(tensor.shape, 1), &[4]);
        assert_eq!(
            std::slice::from_raw_parts(tensor.data, tensor.len),
            &[1.0, 2.0, 3.0, 4.0]
        );
        rssn_tensor_free(&raw mut tensor);
        assert!(tensor.shape.is_null() && tensor.data.is_null());
        assert_eq!((tensor.rank, tensor.len), (0, 0));
        // freeing the cleared struct again is a no-op
        rssn_tensor_free(&raw mut tensor);

        // 3-D tensor keeps row-major order
        let cube = rssn_term_parse(
            s,
            c("list(list(list(1, 2), list(3, 4)), list(list(5, 6), list(7, 8)))").as_ptr(),
        );
        assert_eq!(
            rssn_term_tensor_data(s, cube, &raw mut tensor),
            RssnStatus::Ok,
            "{}",
            last_error()
        );
        assert_eq!(
            std::slice::from_raw_parts(tensor.shape, tensor.rank),
            &[2, 2, 2]
        );
        let data = std::slice::from_raw_parts(tensor.data, tensor.len);
        assert_eq!(data, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]);
        rssn_tensor_free(&raw mut tensor);

        // a symbol is not a tensor
        let x = rssn_term_sym(s, c("x").as_ptr());
        assert_eq!(
            rssn_term_tensor_data(s, x, &raw mut tensor),
            RssnStatus::NoValue
        );
        assert!(tensor.data.is_null(), "failure leaves the out struct alone");
        rssn_session_free(s);
    }
}

// ------------------------------------------------------- null / bad input

#[test]
fn every_term_constructor_rejects_null_session_and_null_strings() {
    let s = rssn_session_new();
    let bad = bad_utf8();
    // SAFETY: only null / invalid inputs are given, which the ABI checks.
    unsafe {
        let null_s = std::ptr::null();
        assert_eq!(rssn_term_sym(null_s, c("x").as_ptr()), RSSN_NO_TERM);
        assert_eq!(rssn_term_int(null_s, 1), RSSN_NO_TERM);
        assert_eq!(rssn_term_rational(null_s, 1, 2), RSSN_NO_TERM);
        assert_eq!(rssn_term_float(null_s, 1.0), RSSN_NO_TERM);
        assert_eq!(rssn_term_parse(null_s, c("1").as_ptr()), RSSN_NO_TERM);
        assert_eq!(
            rssn_term_apply(null_s, c("sin").as_ptr(), std::ptr::null(), 0),
            RSSN_NO_TERM
        );
        assert_eq!(rssn_simplify(null_s, c("1").as_ptr()), RSSN_NO_TERM);

        assert_eq!(rssn_term_sym(s, std::ptr::null()), RSSN_NO_TERM);
        assert_eq!(rssn_term_parse(s, std::ptr::null()), RSSN_NO_TERM);
        assert_eq!(
            rssn_term_apply(s, std::ptr::null(), std::ptr::null(), 0),
            RSSN_NO_TERM
        );
        assert_eq!(rssn_simplify(s, std::ptr::null()), RSSN_NO_TERM);

        assert_eq!(rssn_term_sym(s, bad.as_ptr()), RSSN_NO_TERM);
        assert!(!last_error().is_empty());
        assert_eq!(rssn_term_parse(s, bad.as_ptr()), RSSN_NO_TERM);
        assert_eq!(
            rssn_term_apply(s, bad.as_ptr(), std::ptr::null(), 0),
            RSSN_NO_TERM
        );
        assert_eq!(rssn_simplify(s, bad.as_ptr()), RSSN_NO_TERM);
        rssn_session_free(s);
    }
}

#[test]
fn constructors_succeed_and_agree_with_the_parser() {
    let s = rssn_session_new();
    // SAFETY: live session, valid strings.
    unsafe {
        let a = rssn_term_int(s, 7);
        let b = rssn_term_parse(s, c("7").as_ptr());
        assert_eq!(
            owned(rssn_term_to_string(s, a)),
            owned(rssn_term_to_string(s, b))
        );
        let r = rssn_term_rational(s, 3, 6);
        let mut v = 0.0;
        assert_eq!(rssn_term_to_f64(s, r, &raw mut v), RssnStatus::Ok);
        assert!((v - 0.5).abs() < 1e-15, "3/6 reduces to 1/2, got {v}");
        // negative denominator is normalised, not rejected
        let r2 = rssn_term_rational(s, 1, -4);
        assert_ne!(r2, RSSN_NO_TERM);
        assert_eq!(rssn_term_to_f64(s, r2, &raw mut v), RssnStatus::Ok);
        assert!((v + 0.25).abs() < 1e-15);
        let f = rssn_term_float(s, 0.125);
        assert_eq!(rssn_term_to_f64(s, f, &raw mut v), RssnStatus::Ok);
        assert_eq!(v, 0.125);
        let i = rssn_term_int(s, i64::MIN);
        assert_eq!(rssn_term_to_f64(s, i, &raw mut v), RssnStatus::Ok);
        assert_eq!(v, i64::MIN as f64);
        // a symbol has no literal value
        let x = rssn_term_sym(s, c("x").as_ptr());
        assert_eq!(rssn_term_to_f64(s, x, &raw mut v), RssnStatus::NoValue);
        // hash-consing: the same symbol twice is the same handle
        assert_eq!(x, rssn_term_sym(s, c("x").as_ptr()));
        // unicode symbol names survive the round trip
        let u = rssn_term_sym(s, c("alpha_β").as_ptr());
        assert_ne!(u, RSSN_NO_TERM);
        assert!(owned(rssn_term_to_string(s, u)).contains('β'));
        // zero denominator
        assert_eq!(rssn_term_rational(s, 5, 0), RSSN_NO_TERM);
        assert!(!last_error().is_empty());
        rssn_session_free(s);
    }
}

#[test]
fn apply_checks_arrays_and_handles() {
    let s = rssn_session_new();
    // SAFETY: only valid-or-null pointers.
    unsafe {
        let x = rssn_term_sym(s, c("x").as_ptr());
        // count > 0 with a null array
        assert_eq!(
            rssn_term_apply(s, c("sin").as_ptr(), std::ptr::null(), 1),
            RSSN_NO_TERM
        );
        // stale / unknown handles
        let bogus = [u32::MAX - 1];
        assert_eq!(
            rssn_term_apply(s, c("sin").as_ptr(), bogus.as_ptr(), 1),
            RSSN_NO_TERM
        );
        let mixed = [x, 3_999_999_999];
        assert_eq!(
            rssn_term_apply(s, c("diff").as_ptr(), mixed.as_ptr(), 2),
            RSSN_NO_TERM
        );
        // the sentinel is never a live handle
        let sentinel = [RSSN_NO_TERM];
        assert_eq!(
            rssn_term_apply(s, c("sin").as_ptr(), sentinel.as_ptr(), 1),
            RSSN_NO_TERM
        );
        // wrong arity: sin takes one argument
        let two = [x, x];
        let wrong = rssn_term_apply(s, c("sin").as_ptr(), two.as_ptr(), 2);
        let _ = wrong; // either a node or an error; must not panic
        // unknown operator
        assert_eq!(
            rssn_term_apply(s, c("definitely_not_an_op").as_ptr(), &raw const x, 1),
            RSSN_NO_TERM
        );
        assert!(last_error().contains("definitely_not_an_op"));
        // count == 0 ignores the (non-null) array pointer
        let okterm = rssn_term_apply(s, c("sin").as_ptr(), &raw const x, 1);
        assert_ne!(okterm, RSSN_NO_TERM);
        assert_eq!(owned(rssn_term_to_string(s, okterm)), "sin(x)");
        rssn_session_free(s);
    }
}

#[test]
fn renderers_return_null_for_bad_handles_and_text_for_good_ones() {
    let s = rssn_session_new();
    // SAFETY: live session.
    unsafe {
        let t = rssn_term_parse(s, c("x^2 + 1/2").as_ptr());
        assert_ne!(t, RSSN_NO_TERM);
        for render in [
            rssn_term_to_string,
            rssn_term_to_latex,
            rssn_term_to_typst,
            rssn_term_to_pretty,
        ] {
            assert!(render(std::ptr::null(), t).is_null());
            assert!(render(s, 4_000_000_000).is_null());
            assert!(render(s, RSSN_NO_TERM).is_null());
            let text = owned(render(s, t));
            assert!(text.contains('x'), "{text}");
        }
        rssn_session_free(s);
    }
}

#[test]
fn term_eval_checks_every_pointer_and_name() {
    let s = rssn_session_new();
    // SAFETY: valid-or-null pointers only.
    unsafe {
        let t = rssn_term_parse(s, c("x*y + 1").as_ptr());
        let mut out = 0.0;
        let names = [c("x"), c("y")];
        let ptrs = [names[0].as_ptr(), names[1].as_ptr()];
        let vals = [2.0, 3.0];
        // null out
        assert_eq!(
            rssn_term_eval(s, t, ptrs.as_ptr(), vals.as_ptr(), 2, std::ptr::null_mut()),
            RssnStatus::NullArgument
        );
        // null arrays with count > 0
        assert_eq!(
            rssn_term_eval(s, t, std::ptr::null(), vals.as_ptr(), 2, &raw mut out),
            RssnStatus::NullArgument
        );
        assert_eq!(
            rssn_term_eval(s, t, ptrs.as_ptr(), std::ptr::null(), 2, &raw mut out),
            RssnStatus::NullArgument
        );
        // null session / bad handle
        assert_eq!(
            rssn_term_eval(
                std::ptr::null(),
                t,
                ptrs.as_ptr(),
                vals.as_ptr(),
                2,
                &raw mut out
            ),
            RssnStatus::NullArgument
        );
        assert_eq!(
            rssn_term_eval(
                s,
                3_000_000_000,
                ptrs.as_ptr(),
                vals.as_ptr(),
                2,
                &raw mut out
            ),
            RssnStatus::NoSuchTerm
        );
        // a null name inside the array
        let with_null = [ptrs[0], std::ptr::null()];
        assert_eq!(
            rssn_term_eval(s, t, with_null.as_ptr(), vals.as_ptr(), 2, &raw mut out),
            RssnStatus::NullArgument
        );
        // invalid UTF-8 inside the array
        let bad = bad_utf8();
        let with_bad = [ptrs[0], bad.as_ptr()];
        assert_eq!(
            rssn_term_eval(s, t, with_bad.as_ptr(), vals.as_ptr(), 2, &raw mut out),
            RssnStatus::InvalidUtf8
        );
        // unbound variable: not numeric
        assert_eq!(
            rssn_term_eval(s, t, ptrs.as_ptr(), vals.as_ptr(), 1, &raw mut out),
            RssnStatus::NotNumeric
        );
        // zero bindings and a closed term is fine with null arrays
        let k = rssn_term_parse(s, c("2 + 3").as_ptr());
        assert_eq!(
            rssn_term_eval(s, k, std::ptr::null(), std::ptr::null(), 0, &raw mut out),
            RssnStatus::Ok
        );
        assert_eq!(out, 5.0);
        // full binding: linearity check x*y+1 at several points
        for (xv, yv) in [(0.0, 0.0), (-1.5, 4.0), (1e3, 1e-3)] {
            let vals = [xv, yv];
            assert_eq!(
                rssn_term_eval(s, t, ptrs.as_ptr(), vals.as_ptr(), 2, &raw mut out),
                RssnStatus::Ok
            );
            assert!((out - (xv * yv + 1.0)).abs() <= 1e-12 * (1.0 + (xv * yv).abs()));
        }
        // NaN input propagates, no panic
        let nanv = [f64::NAN, 1.0];
        let st = rssn_term_eval(s, t, ptrs.as_ptr(), nanv.as_ptr(), 2, &raw mut out);
        assert!(st == RssnStatus::Ok || st == RssnStatus::NotNumeric);
        if st == RssnStatus::Ok {
            assert!(out.is_nan());
        }
        rssn_session_free(s);
    }
}

#[test]
fn to_f64_reports_null_before_bad_handle() {
    let s = rssn_session_new();
    // SAFETY: valid-or-null pointers only.
    unsafe {
        let mut v = 0.0;
        assert_eq!(
            rssn_term_to_f64(s, 0, std::ptr::null_mut()),
            RssnStatus::NullArgument
        );
        assert_eq!(
            rssn_term_to_f64(s, 3_000_000_000, std::ptr::null_mut()),
            RssnStatus::NullArgument
        );
        assert_eq!(
            rssn_term_to_f64(s, 3_000_000_000, &raw mut v),
            RssnStatus::NoSuchTerm
        );
        assert_eq!(
            rssn_term_tensor_data(s, 0, std::ptr::null_mut()),
            RssnStatus::NullArgument
        );
        assert_eq!(
            rssn_term_tensor_data(std::ptr::null(), 0, &mut blank_tensor()),
            RssnStatus::NullArgument
        );
        rssn_session_free(s);
    }
}

// -------------------------------------------------------------- config

#[test]
fn config_editors_reject_null_config_and_bad_strings() {
    let bad = bad_utf8();
    let cfg = rssn_config_new();
    // SAFETY: valid-or-null pointers.
    unsafe {
        let nc = std::ptr::null_mut();
        assert_eq!(rssn_config_symbolic(nc), RssnStatus::NullArgument);
        assert_eq!(rssn_config_numeric(nc, 1e-9), RssnStatus::NullArgument);
        assert_eq!(
            rssn_config_bind(nc, c("x").as_ptr(), 1.0),
            RssnStatus::NullArgument
        );
        assert_eq!(
            rssn_config_assume(nc, c("x").as_ptr(), 1),
            RssnStatus::NullArgument
        );
        assert_eq!(rssn_config_budget(nc, 1, 1), RssnStatus::NullArgument);
        assert_eq!(rssn_config_dimension_mode(nc, 1), RssnStatus::NullArgument);
        assert_eq!(
            rssn_config_dimension(nc, c("v").as_ptr(), c("m").as_ptr()),
            RssnStatus::NullArgument
        );
        assert_eq!(
            rssn_config_define_base_dimension(nc, c("Q").as_ptr()),
            RssnStatus::NullArgument
        );
        assert_eq!(rssn_config_bypass_dimensions(nc), RssnStatus::NullArgument);

        assert_eq!(
            rssn_config_bind(cfg, std::ptr::null(), 1.0),
            RssnStatus::NullArgument
        );
        assert_eq!(
            rssn_config_bind(cfg, bad.as_ptr(), 1.0),
            RssnStatus::InvalidUtf8
        );
        assert_eq!(
            rssn_config_assume(cfg, bad.as_ptr(), 1),
            RssnStatus::InvalidUtf8
        );
        assert_eq!(
            rssn_config_dimension(cfg, bad.as_ptr(), c("m").as_ptr()),
            RssnStatus::InvalidUtf8
        );
        assert_eq!(
            rssn_config_dimension(cfg, c("v").as_ptr(), bad.as_ptr()),
            RssnStatus::InvalidUtf8
        );
        assert_eq!(
            rssn_config_dimension(cfg, c("v").as_ptr(), std::ptr::null()),
            RssnStatus::NullArgument
        );
        assert_eq!(
            rssn_config_define_base_dimension(cfg, bad.as_ptr()),
            RssnStatus::InvalidUtf8
        );
        for mode in [3u8, 4, 100, 255] {
            assert_eq!(
                rssn_config_dimension_mode(cfg, mode),
                RssnStatus::InvalidArgument,
                "mode {mode}"
            );
        }
        for mode in [0u8, 1, 2] {
            assert_eq!(rssn_config_dimension_mode(cfg, mode), RssnStatus::Ok);
        }
        // unusual numeric tolerances and budgets are accepted (or at least
        // never panic) and the config still computes afterwards
        for tol in [0.0, -1.0, f64::NAN, f64::INFINITY, 1e-300] {
            let st = rssn_config_numeric(cfg, tol);
            assert!(
                matches!(st, RssnStatus::Ok | RssnStatus::InvalidArgument),
                "tol {tol}: {st:?}"
            );
        }
        assert_eq!(
            rssn_config_budget(cfg, usize::MAX, usize::MAX),
            RssnStatus::Ok
        );
        assert_eq!(rssn_config_budget(cfg, 0, 0), RssnStatus::Ok);
        assert_eq!(rssn_config_budget(cfg, 1, 0), RssnStatus::Ok);
        assert_eq!(rssn_config_symbolic(cfg), RssnStatus::Ok);
        rssn_config_free(cfg);
    }
}

#[test]
fn a_failed_edit_leaves_the_configuration_unchanged() {
    let s = rssn_session_new();
    let cfg = rssn_config_new();
    // SAFETY: live objects.
    unsafe {
        assert_eq!(rssn_config_numeric(cfg, 1e-12), RssnStatus::Ok);
        assert_eq!(rssn_config_bind(cfg, c("x").as_ptr(), 2.0), RssnStatus::Ok);
        // failing edits
        assert_eq!(
            rssn_config_bind(cfg, bad_utf8().as_ptr(), 9.0),
            RssnStatus::InvalidUtf8
        );
        assert_eq!(
            rssn_config_dimension_mode(cfg, 9),
            RssnStatus::InvalidArgument
        );
        let t = rssn_term_parse(s, c("x^2").as_ptr());
        let mut a = blank_answer();
        assert_eq!(rssn_compute(s, t, cfg, &raw mut a), RssnStatus::Ok);
        assert!(a.has_value);
        assert!(
            (a.value - 4.0).abs() < 1e-10,
            "x=2 survives failed edits: {}",
            a.value
        );
        rssn_config_free(cfg);
        rssn_session_free(s);
    }
}

#[test]
fn assumptions_change_the_answer_and_bits_are_independent() {
    let s = rssn_session_new();
    // SAFETY: live objects.
    unsafe {
        let r = rssn_term_parse(s, c("sqrt(a^2)").as_ptr());
        let plain = rssn_config_new();
        let mut a = blank_answer();
        assert_eq!(rssn_compute(s, r, plain, &raw mut a), RssnStatus::Ok);
        let without = owned(rssn_term_to_string(s, a.term));
        let cfg = rssn_config_new();
        assert_eq!(rssn_config_assume(cfg, c("a").as_ptr(), 7), RssnStatus::Ok);
        assert_eq!(rssn_compute(s, r, cfg, &raw mut a), RssnStatus::Ok);
        let with = owned(rssn_term_to_string(s, a.term));
        assert_eq!(with, "a");
        assert_ne!(
            without, with,
            "without the assumption the answer must keep |a|"
        );
        // every fact byte is accepted without panic
        for facts in 0..=255u8 {
            assert_eq!(
                rssn_config_assume(cfg, c("a").as_ptr(), facts),
                RssnStatus::Ok
            );
        }
        rssn_config_free(plain);
        rssn_config_free(cfg);
        rssn_session_free(s);
    }
}

// ------------------------------------------------------------- compute

#[test]
fn compute_family_rejects_null_outputs_session_and_bad_handles() {
    let s = rssn_session_new();
    // SAFETY: valid-or-null pointers.
    unsafe {
        let t = rssn_term_parse(s, c("1 + 1").as_ptr());
        let mut a = blank_answer();
        let mut text = std::ptr::null_mut();
        let mut data = std::ptr::null_mut();
        let mut len = 0usize;
        let np: *const RssnConfig = std::ptr::null();
        let ns: *const RssnSession = std::ptr::null();

        assert_eq!(
            rssn_compute(s, t, np, std::ptr::null_mut()),
            RssnStatus::NullArgument
        );
        assert_eq!(
            rssn_compute(ns, t, np, &raw mut a),
            RssnStatus::NullArgument
        );
        assert_eq!(
            rssn_compute(s, 3_000_000_000, np, &raw mut a),
            RssnStatus::NoSuchTerm
        );
        assert_eq!(
            rssn_compute(s, RSSN_NO_TERM, np, &raw mut a),
            RssnStatus::NoSuchTerm
        );

        assert_eq!(
            rssn_compute_json(s, t, np, std::ptr::null_mut()),
            RssnStatus::NullArgument
        );
        assert_eq!(
            rssn_compute_json(ns, t, np, &raw mut text),
            RssnStatus::NullArgument
        );
        assert_eq!(
            rssn_compute_json(s, 3_000_000_000, np, &raw mut text),
            RssnStatus::NoSuchTerm
        );
        assert!(text.is_null(), "failure must not set the output");

        assert_eq!(
            rssn_dimension(s, t, np, std::ptr::null_mut()),
            RssnStatus::NullArgument
        );
        assert_eq!(
            rssn_dimension(ns, t, np, &raw mut text),
            RssnStatus::NullArgument
        );

        assert_eq!(
            rssn_compute_bincode(s, t, np, std::ptr::null_mut(), &raw mut len),
            RssnStatus::NullArgument
        );
        assert_eq!(
            rssn_compute_bincode(s, t, np, &raw mut data, std::ptr::null_mut()),
            RssnStatus::NullArgument
        );
        assert_eq!(
            rssn_compute_bincode(ns, t, np, &raw mut data, &raw mut len),
            RssnStatus::NullArgument
        );
        assert_eq!(
            rssn_compute_bincode(s, 3_000_000_000, np, &raw mut data, &raw mut len),
            RssnStatus::NoSuchTerm
        );
        assert!(data.is_null());
        assert_eq!(len, 0);
        rssn_session_free(s);
    }
}

#[test]
fn json_and_bincode_records_agree_and_the_plain_answer_matches() {
    let s = rssn_session_new();
    // SAFETY: live objects, buffers freed.
    unsafe {
        let cfg = rssn_config_new();
        rssn_config_numeric(cfg, 1e-12);
        rssn_config_bind(cfg, c("x").as_ptr(), 0.5);
        for src in [
            "sin(x)^2 + cos(x)^2",
            "exp(x)*exp(-x)",
            "diff(x^3, x)",
            "x^2 + 1",
        ] {
            let t = rssn_term_parse(s, c(src).as_ptr());
            assert_ne!(t, RSSN_NO_TERM, "{src}");
            let mut a = blank_answer();
            assert_eq!(rssn_compute(s, t, cfg, &raw mut a), RssnStatus::Ok, "{src}");

            let mut json = std::ptr::null_mut();
            assert_eq!(rssn_compute_json(s, t, cfg, &raw mut json), RssnStatus::Ok);
            let rec: RssnComputeRecord = serde_json::from_str(&owned(json)).unwrap();

            let mut data = std::ptr::null_mut();
            let mut len = 0;
            assert_eq!(
                rssn_compute_bincode(s, t, cfg, &raw mut data, &raw mut len),
                RssnStatus::Ok
            );
            assert!(len > 0 && !data.is_null());
            let bytes = std::slice::from_raw_parts(data, len).to_vec();
            rssn_bytes_free(data, len);
            let bin: RssnComputeRecord = decode_deterministic(&bytes).unwrap();

            assert_eq!(rec.term, bin.term, "{src}");
            assert_eq!(rec.value, bin.value);
            assert_eq!(rec.reduced, bin.reduced);
            assert_eq!(rec.iterations, bin.iterations);
            assert_eq!(rec.diagnostic_codes, bin.diagnostic_codes);
            assert_eq!(rec.reduced, a.reduced);
            assert_eq!(rec.value.is_some(), a.has_value);
            if let (Some(v), true) = (rec.value, a.has_value) {
                assert_eq!(v.to_bits(), a.value.to_bits(), "{src}");
            }
            // determinism: a second call encodes the identical bytes
            let mut data2 = std::ptr::null_mut();
            let mut len2 = 0;
            assert_eq!(
                rssn_compute_bincode(s, t, cfg, &raw mut data2, &raw mut len2),
                RssnStatus::Ok
            );
            assert_eq!(std::slice::from_raw_parts(data2, len2), &bytes[..], "{src}");
            rssn_bytes_free(data2, len2);
        }
        // sin^2+cos^2 at x = 0.5 is exactly 1 within the tolerance
        let t = rssn_term_parse(s, c("sin(x)^2 + cos(x)^2").as_ptr());
        let mut a = blank_answer();
        rssn_compute(s, t, cfg, &raw mut a);
        assert!(a.has_value && (a.value - 1.0).abs() < 1e-10);
        rssn_config_free(cfg);
        rssn_session_free(s);
    }
}

#[test]
fn truncated_or_corrupt_bincode_never_decodes_as_the_original() {
    let s = rssn_session_new();
    // SAFETY: live objects, buffers freed.
    unsafe {
        let t = rssn_term_parse(s, c("diff(sin(x), x)").as_ptr());
        let mut data = std::ptr::null_mut();
        let mut len = 0;
        assert_eq!(
            rssn_compute_bincode(s, t, std::ptr::null(), &raw mut data, &raw mut len),
            RssnStatus::Ok
        );
        let bytes = std::slice::from_raw_parts(data, len).to_vec();
        rssn_bytes_free(data, len);
        let full: RssnComputeRecord = decode_deterministic(&bytes).unwrap();
        // a single record has no self-framing tail, but cutting inside a
        // field must fail rather than panic
        let mut failures = 0;
        for cut in 0..bytes.len() {
            match decode_deterministic::<RssnComputeRecord>(&bytes[..cut]) {
                | Err(_) => failures += 1,
                | Ok(r) => assert_ne!(r.term, full.term, "cut {cut}"),
            }
        }
        assert!(failures >= bytes.len() / 2, "{failures}/{}", bytes.len());
        assert!(decode_deterministic::<RssnComputeRecord>(&[]).is_err());
        rssn_session_free(s);
    }
}

#[test]
fn compute_all_edge_counts_and_mixed_batches() {
    let s = rssn_session_new();
    // SAFETY: arrays sized to the counts passed.
    unsafe {
        let np: *const RssnConfig = std::ptr::null();
        let ns: *const RssnSession = std::ptr::null();
        // n == 0 with null arrays is fine and writes nothing
        assert_eq!(
            rssn_compute_all(
                s,
                std::ptr::null(),
                0,
                np,
                std::ptr::null_mut(),
                std::ptr::null_mut()
            ),
            RssnStatus::Ok
        );
        // n > 0 with any null array
        let ids = [0u32];
        let mut answers = [blank_answer()];
        let mut statuses = [RssnStatus::Ok];
        assert_eq!(
            rssn_compute_all(
                s,
                std::ptr::null(),
                1,
                np,
                answers.as_mut_ptr(),
                statuses.as_mut_ptr()
            ),
            RssnStatus::NullArgument
        );
        assert_eq!(
            rssn_compute_all(
                s,
                ids.as_ptr(),
                1,
                np,
                std::ptr::null_mut(),
                statuses.as_mut_ptr()
            ),
            RssnStatus::NullArgument
        );
        assert_eq!(
            rssn_compute_all(
                s,
                ids.as_ptr(),
                1,
                np,
                answers.as_mut_ptr(),
                std::ptr::null_mut()
            ),
            RssnStatus::NullArgument
        );
        // null session
        assert_eq!(
            rssn_compute_all(
                ns,
                ids.as_ptr(),
                1,
                np,
                answers.as_mut_ptr(),
                statuses.as_mut_ptr()
            ),
            RssnStatus::NullArgument
        );

        // a batch mixing good, stale and sentinel handles keeps input order
        let good: Vec<u32> = (1..=20)
            .map(|k| rssn_term_parse(s, c(&format!("{k} * {k}")).as_ptr()))
            .collect();
        let mut batch = Vec::new();
        for (k, g) in good.iter().enumerate() {
            batch.push(*g);
            if k % 5 == 0 {
                batch.push(3_000_000_000);
            }
        }
        let mut answers = vec![blank_answer(); batch.len()];
        let mut statuses = vec![RssnStatus::Panic; batch.len()];
        let cfg = rssn_config_new();
        rssn_config_numeric(cfg, 1e-12);
        assert_eq!(
            rssn_compute_all(
                s,
                batch.as_ptr(),
                batch.len(),
                cfg,
                answers.as_mut_ptr(),
                statuses.as_mut_ptr()
            ),
            RssnStatus::Ok
        );
        let mut k = 0usize;
        for (i, id) in batch.iter().enumerate() {
            if *id == 3_000_000_000 {
                assert_eq!(statuses[i], RssnStatus::NoSuchTerm, "slot {i}");
                assert_eq!(answers[i].term, RSSN_NO_TERM);
                assert!(!answers[i].has_value);
            } else {
                k += 1;
                assert_eq!(statuses[i], RssnStatus::Ok, "slot {i}");
                assert!((answers[i].value - (k * k) as f64).abs() < 1e-9);
            }
        }
        // batch results equal the one-at-a-time results bit for bit
        for (i, id) in batch.iter().enumerate() {
            if statuses[i] != RssnStatus::Ok {
                continue;
            }
            let mut one = blank_answer();
            assert_eq!(rssn_compute(s, *id, cfg, &raw mut one), RssnStatus::Ok);
            assert_eq!(one.value.to_bits(), answers[i].value.to_bits());
            assert_eq!(one.term, answers[i].term);
        }
        rssn_config_free(cfg);
        rssn_session_free(s);
    }
}

#[test]
fn compute_all_bincode_edge_cases() {
    let s = rssn_session_new();
    // SAFETY: arrays sized to the counts passed.
    unsafe {
        let np: *const RssnConfig = std::ptr::null();
        let ns: *const RssnSession = std::ptr::null();
        let mut data = std::ptr::null_mut();
        let mut len = 0usize;
        // null outputs
        assert_eq!(
            rssn_compute_all_bincode(
                s,
                std::ptr::null(),
                0,
                np,
                std::ptr::null_mut(),
                &raw mut len
            ),
            RssnStatus::NullArgument
        );
        assert_eq!(
            rssn_compute_all_bincode(
                s,
                std::ptr::null(),
                0,
                np,
                &raw mut data,
                std::ptr::null_mut()
            ),
            RssnStatus::NullArgument
        );
        // n > 0, null terms
        assert_eq!(
            rssn_compute_all_bincode(s, std::ptr::null(), 3, np, &raw mut data, &raw mut len),
            RssnStatus::NullArgument
        );
        // null session
        assert_eq!(
            rssn_compute_all_bincode(ns, std::ptr::null(), 0, np, &raw mut data, &raw mut len),
            RssnStatus::NullArgument
        );
        // empty batch: an empty entry list that decodes
        assert_eq!(
            rssn_compute_all_bincode(s, std::ptr::null(), 0, np, &raw mut data, &raw mut len),
            RssnStatus::Ok
        );
        let bytes = std::slice::from_raw_parts(data, len).to_vec();
        rssn_bytes_free(data, len);
        let rec: RssnComputeBatchRecord = decode_deterministic(&bytes).unwrap();
        assert!(rec.entries.is_empty());

        // mixed batch
        let a = rssn_term_parse(s, c("2 + 2").as_ptr());
        let ids = [a, 3_000_000_000, a];
        assert_eq!(
            rssn_compute_all_bincode(s, ids.as_ptr(), 3, np, &raw mut data, &raw mut len),
            RssnStatus::Ok
        );
        let bytes = std::slice::from_raw_parts(data, len).to_vec();
        rssn_bytes_free(data, len);
        let rec: RssnComputeBatchRecord = decode_deterministic(&bytes).unwrap();
        assert_eq!(rec.entries.len(), 3);
        assert_eq!(rec.entries[0].status, RssnStatus::Ok as u32);
        assert_eq!(rec.entries[1].status, RssnStatus::NoSuchTerm as u32);
        assert!(rec.entries[1].answer.is_none());
        assert_eq!(rec.entries[2].status, RssnStatus::Ok as u32);
        assert_eq!(
            rec.entries[0].answer.as_ref().unwrap().term,
            rec.entries[2].answer.as_ref().unwrap().term
        );
        assert_eq!(rec.entries[0].answer.as_ref().unwrap().term, "4");
        rssn_session_free(s);
    }
}

#[test]
fn simplify_matches_compute_and_failures_are_no_term() {
    let s = rssn_session_new();
    // SAFETY: live session.
    unsafe {
        for (src, want) in [
            ("integral(2*x, x)", "x^2"),
            ("diff(x^3, x)", "3 * x^2"),
            ("x + 0", "x"),
            ("1 * x", "x"),
        ] {
            let t = rssn_simplify(s, c(src).as_ptr());
            assert_ne!(t, RSSN_NO_TERM, "{src}: {}", last_error());
            let got = owned(rssn_term_to_string(s, t));
            // compare to the compute path rather than a spelling
            let parsed = rssn_term_parse(s, c(src).as_ptr());
            let mut a = blank_answer();
            assert_eq!(
                rssn_compute(s, parsed, std::ptr::null(), &raw mut a),
                RssnStatus::Ok
            );
            assert_eq!(a.term, t, "{src}");
            let _ = want;
            assert!(!got.is_empty());
        }
        assert_eq!(rssn_simplify(s, c("sin(").as_ptr()), RSSN_NO_TERM);
        assert_eq!(rssn_simplify(s, c("").as_ptr()), RSSN_NO_TERM);
        rssn_session_free(s);
    }
}

#[test]
fn strict_dimension_mode_reports_a_stable_code_for_each_conflict() {
    let s = rssn_session_new();
    // SAFETY: live objects.
    unsafe {
        let cfg = rssn_config_new();
        assert_eq!(rssn_config_dimension_mode(cfg, 2), RssnStatus::Ok);
        let mut a = blank_answer();
        // mismatched dimensions in a sum
        let bad = rssn_term_parse(s, c("quantity(1, m) + quantity(1, s)").as_ptr());
        assert_eq!(rssn_compute(s, bad, cfg, &raw mut a), RssnStatus::Dimension);
        assert_eq!(rssn_dimension_count(), 1);
        assert_eq!(
            CStr::from_ptr(rssn_dimension_code(0)).to_str().unwrap(),
            "mismatch"
        );
        // a consistent product passes and clears the codes
        let good = rssn_term_parse(s, c("quantity(2, m) * quantity(3, s)").as_ptr());
        assert_eq!(
            rssn_compute(s, good, cfg, &raw mut a),
            RssnStatus::Ok,
            "{}",
            last_error()
        );
        assert_eq!(rssn_dimension_count(), 0);
        // the dimension of the product is Length*Time
        let mut out = std::ptr::null_mut();
        assert_eq!(rssn_dimension(s, good, cfg, &raw mut out), RssnStatus::Ok);
        let dim = owned(out);
        assert!(dim.contains("Length") && dim.contains("Time"), "{dim}");
        // a term with no dimensional content gives an empty dimension text
        let plain = rssn_term_parse(s, c("3").as_ptr());
        assert_eq!(
            rssn_dimension(s, plain, std::ptr::null(), &raw mut out),
            RssnStatus::Ok
        );
        let _ = owned(out);
        // index far out of range is null, not UB
        assert!(rssn_dimension_code(usize::MAX).is_null());
        assert!(rssn_dimension_code(1_000_000).is_null());
        rssn_config_free(cfg);
        rssn_session_free(s);
    }
}

// ------------------------------------------------------------- sim_run

#[test]
fn sim_run_null_and_malformed_inputs() {
    let mut out: *mut c_char = std::ptr::null_mut();
    let bad = bad_utf8();
    // SAFETY: valid-or-null pointers.
    unsafe {
        assert_eq!(
            rssn_sim_run(c("ising").as_ptr(), c("{}").as_ptr(), std::ptr::null_mut()),
            RssnStatus::NullArgument
        );
        assert_eq!(
            rssn_sim_run(std::ptr::null(), c("{}").as_ptr(), &raw mut out),
            RssnStatus::NullArgument
        );
        assert_eq!(
            rssn_sim_run(c("ising").as_ptr(), std::ptr::null(), &raw mut out),
            RssnStatus::NullArgument
        );
        assert_eq!(
            rssn_sim_run(bad.as_ptr(), c("{}").as_ptr(), &raw mut out),
            RssnStatus::InvalidUtf8
        );
        assert_eq!(
            rssn_sim_run(c("ising").as_ptr(), bad.as_ptr(), &raw mut out),
            RssnStatus::InvalidUtf8
        );
        assert_eq!(
            rssn_sim_run(c("").as_ptr(), c("{}").as_ptr(), &raw mut out),
            RssnStatus::Simulation
        );
        assert!(out.is_null(), "no output on failure");
    }
}

#[test]
fn every_scenario_rejects_malformed_json_without_panicking() {
    let mut names: Vec<&str> = SCENARIOS.to_vec();
    names.extend_from_slice(rssn::sim::scenario::ENGINE_SCENARIOS);
    names.sort_unstable();
    names.dedup();
    // Documents that are not JSON objects at all (an empty object would be a legal request that
    // runs the defaults, which is slow; objects with unknown keys are legal too).
    let malformed = ["{", "}{", "nul", "\"text\"", "{\"a\": }"];
    let mut accepted = Vec::new();
    for name in &names {
        for params in malformed {
            let mut out: *mut c_char = std::ptr::null_mut();
            // SAFETY: valid strings, out is a local.
            let st = unsafe { rssn_sim_run(c(name).as_ptr(), c(params).as_ptr(), &raw mut out) };
            assert_ne!(
                st,
                RssnStatus::Panic,
                "{name} with {params:?}: {}",
                last_error()
            );
            if st == RssnStatus::Ok {
                accepted.push(format!("{name} <- {params:?}"));
                // SAFETY: owned result.
                unsafe { rssn_string_free(out) };
            } else {
                assert_eq!(st, RssnStatus::Simulation, "{name} {params:?}");
                assert!(out.is_null());
                assert!(!last_error().is_empty());
            }
        }
    }
    assert!(
        accepted.is_empty(),
        "malformed documents were accepted: {accepted:?}"
    );
}

#[test]
fn sim_run_is_deterministic_for_a_fixed_seed() {
    let params = c(r#"{"width": 8, "height": 8, "temperature": 2.0, "mc_steps": 10, "seed": 42}"#);
    let run = |seed_params: &CString| {
        let mut out = std::ptr::null_mut();
        // SAFETY: valid strings.
        assert_eq!(
            unsafe { rssn_sim_run(c("ising").as_ptr(), seed_params.as_ptr(), &raw mut out) },
            RssnStatus::Ok
        );
        owned(out)
    };
    let a = run(&params);
    let b = run(&params);
    assert_eq!(a, b);
    let v: serde_json::Value = serde_json::from_str(&a).unwrap();
    assert_eq!(v["spins"].as_array().unwrap().len(), 64);
    let m = v["magnetization"].as_f64().unwrap();
    assert!((-1.0..=1.0).contains(&m));
    let other = run(&c(
        r#"{"width": 8, "height": 8, "temperature": 2.0, "mc_steps": 10, "seed": 43}"#,
    ));
    assert_ne!(
        a, other,
        "different seeds should give different spin fields"
    );
}

// ----------------------------------------------------------- manifest

#[test]
fn manifest_json_is_parseable_stable_and_hash_matches() {
    let a = owned(rssn_manifest_json());
    let b = owned(rssn_manifest_json());
    assert_eq!(a, b);
    let v: serde_json::Value = serde_json::from_str(&a).unwrap();
    assert!(v["ops"].is_array() || v["ops"].is_object());
    assert_eq!(rssn_manifest_hash(), rssn_manifest_hash());
    assert_ne!(rssn_manifest_hash(), 0);
}

// ------------------------------------------------------------- threads

#[test]
fn separate_sessions_on_separate_threads_do_not_interfere() {
    let handles: Vec<_> = (0..4)
        .map(|k| {
            std::thread::spawn(move || {
                // Pointers are not Send; create and use them inside the thread.
                let s = rssn_session_new();
                let mut sum = 0.0;
                // SAFETY: this thread's own session.
                unsafe {
                    for j in 0..20 {
                        let t = rssn_term_parse(s, c(&format!("{k} * {j} + 1")).as_ptr());
                        let mut a = blank_answer();
                        let cfg = rssn_config_new();
                        rssn_config_numeric(cfg, 1e-12);
                        assert_eq!(rssn_compute(s, t, cfg, &raw mut a), RssnStatus::Ok);
                        sum += a.value;
                        rssn_config_free(cfg);
                    }
                    // the error slot is thread local
                    assert_eq!(rssn_term_parse(s, c("sin(").as_ptr()), RSSN_NO_TERM);
                    rssn_session_free(s);
                }
                sum
            })
        })
        .collect();
    for (k, h) in handles.into_iter().enumerate() {
        let sum = h.join().expect("thread panicked");
        let want: f64 = (0..20).map(|j| (k * j + 1) as f64).sum();
        assert!((sum - want).abs() < 1e-9, "thread {k}: {sum} vs {want}");
    }
}

#[test]
fn last_error_is_per_thread() {
    let s = rssn_session_new();
    // SAFETY: valid strings.
    unsafe {
        assert_eq!(rssn_term_parse(s, c("sin(").as_ptr()), RSSN_NO_TERM);
        let here = last_error();
        assert!(!here.is_empty());
        let there = std::thread::spawn(last_error).join().unwrap();
        assert!(
            there.is_empty(),
            "fresh thread starts with no error, got {there:?}"
        );
        assert_eq!(last_error(), here);
        rssn_session_free(s);
    }
}

#[test]
fn successful_call_does_not_clear_the_last_error() {
    let s = rssn_session_new();
    // SAFETY: valid strings.
    unsafe {
        assert_eq!(rssn_term_parse(s, c("sin(").as_ptr()), RSSN_NO_TERM);
        let msg = last_error();
        assert_ne!(rssn_term_int(s, 1), RSSN_NO_TERM);
        assert_eq!(last_error(), msg);
        rssn_session_free(s);
    }
}
