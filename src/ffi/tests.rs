use std::ffi::CStr;
use std::ffi::CString;

use super::*;

fn c(text: &str) -> CString {
    CString::new(text).unwrap()
}

/// Takes ownership of an rssn string.
fn owned(raw: *mut c_char) -> String {
    assert!(!raw.is_null(), "{}", last_error());
    // SAFETY: a string just returned by rssn.
    let text = unsafe { CStr::from_ptr(raw) }.to_str().unwrap().to_owned();
    // SAFETY: returned by rssn, freed once.
    unsafe { rssn_string_free(raw) };
    text
}

fn last_error() -> String {
    // SAFETY: rssn keeps the message alive until the next failure.
    unsafe { CStr::from_ptr(rssn_last_error()) }.to_string_lossy().into_owned()
}

#[test]
fn compute_round_trip() {
    let s = rssn_session_new();
    unsafe {
        let x = rssn_term_sym(s, c("x").as_ptr());
        let sin = rssn_term_apply(s, c("sin").as_ptr(), &x, 1);
        let request = rssn_term_apply(s, c("diff").as_ptr(), [sin, x].as_ptr(), 2);
        assert_ne!(request, RSSN_NO_TERM);
        let mut answer = std::mem::zeroed::<RssnAnswer>();
        assert_eq!(rssn_compute(s, request, std::ptr::null(), &mut answer), RssnStatus::Ok);
        assert!(answer.reduced);
        assert_eq!(owned(rssn_term_to_string(s, answer.term)), "cos(x)");
        assert_eq!(owned(rssn_term_to_latex(s, answer.term)), "\\cos(x)");

        // Numeric target with a binding.
        let config = rssn_config_new();
        assert_eq!(rssn_config_numeric(config, 1e-12), RssnStatus::Ok);
        assert_eq!(rssn_config_bind(config, c("x").as_ptr(), 0.0), RssnStatus::Ok);
        assert_eq!(rssn_config_budget(config, 0, 0), RssnStatus::Ok);
        assert_eq!(rssn_compute(s, request, config, &mut answer), RssnStatus::Ok);
        assert!(answer.has_value);
        assert!((answer.value - 1.0).abs() < 1e-12);
        let mut v = 0.0;
        assert_eq!(rssn_term_to_f64(s, answer.term, &mut v), RssnStatus::Ok);
        assert!((v - 1.0).abs() < 1e-12);
        rssn_config_free(config);

        // Assumptions reach conditional identities.
        let config = rssn_config_new();
        assert_eq!(rssn_config_assume(config, c("a").as_ptr(), 7), RssnStatus::Ok);
        let r = rssn_term_parse(s, c("sqrt(a^2)").as_ptr());
        assert_eq!(rssn_compute(s, r, config, &mut answer), RssnStatus::Ok);
        assert_eq!(owned(rssn_term_to_string(s, answer.term)), "a");
        rssn_config_free(config);

        let simplified = rssn_simplify(s, c("integral(2*x, x)").as_ptr());
        assert_eq!(owned(rssn_term_to_string(s, simplified)), "x^2");

        let mut out = 0.0;
        let names = [c("x"), c("y")];
        let name_ptrs = [names[0].as_ptr(), names[1].as_ptr()];
        let t = rssn_term_parse(s, c("x*y + 1").as_ptr());
        assert_eq!(rssn_term_eval(s, t, name_ptrs.as_ptr(), [2.0, 3.0].as_ptr(), 2, &mut out), RssnStatus::Ok);
        assert!((out - 7.0).abs() < 1e-15);
        rssn_session_free(s);
    }
}

#[test]
fn tensors_and_errors() {
    let s = rssn_session_new();
    unsafe {
        let m = rssn_term_parse(s, c("list(list(1, 2, 3), list(4, 5, 6))").as_ptr());
        let mut t = std::mem::zeroed::<RssnTensor>();
        assert_eq!(rssn_term_tensor_data(s, m, &mut t), RssnStatus::Ok, "{}", last_error());
        assert_eq!(std::slice::from_raw_parts(t.shape, t.rank), &[2, 3]);
        assert_eq!(std::slice::from_raw_parts(t.data, t.len), &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        rssn_tensor_free(&mut t);
        assert!(t.data.is_null());

        let ragged = rssn_term_parse(s, c("list(list(1, 2), list(3))").as_ptr());
        assert_eq!(rssn_term_tensor_data(s, ragged, &mut t), RssnStatus::NoValue);

        assert_eq!(rssn_term_parse(s, c("sin(").as_ptr()), RSSN_NO_TERM);
        assert!(!last_error().is_empty());
        assert_eq!(rssn_term_apply(s, c("no_such_op").as_ptr(), std::ptr::null(), 0), RSSN_NO_TERM);
        assert!(last_error().contains("no_such_op"));
        let mut v = 0.0;
        assert_eq!(rssn_term_to_f64(s, 4_000_000_000, &mut v), RssnStatus::NoSuchTerm);
        assert_eq!(rssn_term_to_f64(std::ptr::null(), 0, &mut v), RssnStatus::NullArgument);
        assert_eq!(rssn_term_rational(s, 1, 0), RSSN_NO_TERM);
        rssn_session_free(s);
    }
}

#[test]
fn simulation_by_name() {
    let mut out = std::ptr::null_mut();
    let params = c(r#"{"width": 6, "height": 6, "temperature": 2.0, "mc_steps": 5, "seed": 1}"#);
    unsafe {
        assert_eq!(rssn_sim_run(c("ising").as_ptr(), params.as_ptr(), &mut out), RssnStatus::Ok);
        assert!(owned(out).contains("magnetization"));
        assert_eq!(rssn_sim_run(c("nope").as_ptr(), params.as_ptr(), &mut out), RssnStatus::Simulation);
    }
    assert!(last_error().contains("unknown scenario"));
    // SAFETY: a static NUL-terminated string.
    let version = unsafe { CStr::from_ptr(rssn_version()) };
    assert_eq!(version.to_str().unwrap(), env!("CARGO_PKG_VERSION"));
}
