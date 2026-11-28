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
    unsafe { CStr::from_ptr(rssn_last_error()) }
        .to_string_lossy()
        .into_owned()
}

/// A status-returning entry point failed with `expected` and left a message.
fn assert_status(
    status: RssnStatus,
    expected: RssnStatus,
) {
    assert_eq!(status, expected, "{}", last_error());
    assert_ne!(last_error(), "", "a failing call must leave a message");
}

/// A call succeeded (a success is allowed to leave an old message alone).
fn assert_ok(status: RssnStatus) {
    assert_eq!(status, RssnStatus::Ok, "{}", last_error());
}

/// A handle-returning entry point failed and left a message.
fn assert_no_term(id: u32) {
    assert_eq!(id, RSSN_NO_TERM, "{}", last_error());
    assert_ne!(last_error(), "", "a failing call must leave a message");
}

/// The bincode bytes of one answer and the record they decode to.
unsafe fn bincode_record(
    session: *const RssnSession,
    term: u32,
    config: *const RssnConfig,
) -> (Vec<u8>, RssnComputeRecord) {
    let mut data = std::ptr::null_mut();
    let mut len = 0_usize;
    // SAFETY: the caller passes a live session and configuration.
    assert_eq!(
        unsafe { rssn_compute_bincode(session, term, config, &raw mut data, &raw mut len) },
        RssnStatus::Ok,
        "{}",
        last_error()
    );
    // SAFETY: the call succeeded, so `data` holds `len` bytes.
    let bytes = unsafe { std::slice::from_raw_parts(data, len) }.to_vec();
    // SAFETY: returned by rssn and freed once.
    unsafe { rssn_bytes_free(data, len) };
    let record = crate::kernel::codec::decode_deterministic(&bytes)
        .unwrap_or_else(|e| panic!("the bytes rssn wrote must decode: {e}"));
    (bytes, record)
}

/// The same answer from the JSON entry point, decoded.
unsafe fn json_record(
    session: *const RssnSession,
    term: u32,
    config: *const RssnConfig,
) -> serde_json::Value {
    let mut raw = std::ptr::null_mut();
    // SAFETY: the caller passes a live session and configuration.
    assert_eq!(
        unsafe { rssn_compute_json(session, term, config, &raw mut raw) },
        RssnStatus::Ok,
        "{}",
        last_error()
    );
    serde_json::from_str(&owned(raw)).unwrap_or_else(|e| panic!("rssn wrote invalid JSON: {e}"))
}

#[test]
fn list_schema_keeps_coherent_float_dimensions_without_assigning_units_to_bare_zero() {
    use crate::kernel::signature::DimExpr;
    use crate::kernel::signature::TypeExpr;
    use crate::kernel::value::DynValue;
    use crate::kernel::value::Value;

    let measured = |f, dim| {
        Value::Dyn(DynValue {
            value: Box::new(Value::Float(f)),
            dim,
        })
    };
    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let same_length = DimExpr::Mul(vec![length.clone(), DimExpr::Dimless]);
    let a = measured(2.0, length.clone());
    let b = measured(-0.0, same_length);
    let t = measured(0.0, time);
    let list_schema = |data| {
        let Value::Tensor(tensor) = super::value::list_value(data) else {
            unreachable!()
        };
        tensor.elem
    };
    assert_eq!(
        list_schema(vec![a.clone(), b]),
        TypeExpr::Float.with_dim(length)
    );
    assert_eq!(
        list_schema(vec![a.clone(), Value::Float(0.0)]),
        crate::kernel::infer::term()
    );
    assert_eq!(list_schema(vec![a, t]), crate::kernel::infer::term());
    assert_eq!(
        list_schema(vec![Value::int(0), Value::int(1)]),
        TypeExpr::Int
    );
}

#[test]
fn compute_round_trip() {
    let s = rssn_session_new();
    unsafe {
        let x = rssn_term_sym(s, c("x").as_ptr());
        let sin = rssn_term_apply(s, c("sin").as_ptr(), &raw const x, 1);
        let request = rssn_term_apply(s, c("diff").as_ptr(), [sin, x].as_ptr(), 2);
        assert_ne!(request, RSSN_NO_TERM);
        let mut answer = std::mem::zeroed::<RssnAnswer>();
        assert_eq!(
            rssn_compute(s, request, std::ptr::null(), &raw mut answer),
            RssnStatus::Ok
        );
        assert!(answer.reduced);
        assert_eq!(owned(rssn_term_to_string(s, answer.term)), "cos(x)");
        assert_eq!(owned(rssn_term_to_latex(s, answer.term)), "\\cos(x)");

        // Numeric target with a binding.
        let config = rssn_config_new();
        assert_eq!(rssn_config_numeric(config, 1e-12), RssnStatus::Ok);
        assert_eq!(
            rssn_config_bind(config, c("x").as_ptr(), 0.0),
            RssnStatus::Ok
        );
        assert_eq!(rssn_config_budget(config, 0, 0), RssnStatus::Ok);
        assert_eq!(
            rssn_compute(s, request, config, &raw mut answer),
            RssnStatus::Ok
        );
        assert!(answer.has_value);
        assert!((answer.value - 1.0).abs() < 1e-12);
        let mut v = 0.0;
        assert_eq!(rssn_term_to_f64(s, answer.term, &raw mut v), RssnStatus::Ok);
        assert!((v - 1.0).abs() < 1e-12);
        rssn_config_free(config);

        // Assumptions reach conditional identities.
        let config = rssn_config_new();
        assert_eq!(
            rssn_config_assume(config, c("a").as_ptr(), 7),
            RssnStatus::Ok
        );
        let r = rssn_term_parse(s, c("sqrt(a^2)").as_ptr());
        assert_eq!(rssn_compute(s, r, config, &raw mut answer), RssnStatus::Ok);
        assert_eq!(owned(rssn_term_to_string(s, answer.term)), "a");
        rssn_config_free(config);

        let simplified = rssn_simplify(s, c("integral(2*x, x)").as_ptr());
        assert_eq!(owned(rssn_term_to_string(s, simplified)), "x^2");

        let mut out = 0.0;
        let names = [c("x"), c("y")];
        let name_ptrs = [names[0].as_ptr(), names[1].as_ptr()];
        let t = rssn_term_parse(s, c("x*y + 1").as_ptr());
        assert_eq!(
            rssn_term_eval(
                s,
                t,
                name_ptrs.as_ptr(),
                [2.0, 3.0].as_ptr(),
                2,
                &raw mut out
            ),
            RssnStatus::Ok
        );
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
        assert_eq!(
            rssn_term_tensor_data(s, m, &raw mut t),
            RssnStatus::Ok,
            "{}",
            last_error()
        );
        assert_eq!(std::slice::from_raw_parts(t.shape, t.rank), &[2, 3]);
        assert_eq!(
            std::slice::from_raw_parts(t.data, t.len),
            &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        );
        rssn_tensor_free(&raw mut t);
        assert!(t.data.is_null());

        let ragged = rssn_term_parse(s, c("list(list(1, 2), list(3))").as_ptr());
        assert_eq!(
            rssn_term_tensor_data(s, ragged, &raw mut t),
            RssnStatus::NoValue
        );

        assert_eq!(rssn_term_parse(s, c("sin(").as_ptr()), RSSN_NO_TERM);
        assert_ne!(last_error(), "");
        assert_eq!(
            rssn_term_apply(s, c("no_such_op").as_ptr(), std::ptr::null(), 0),
            RSSN_NO_TERM
        );
        assert!(last_error().contains("no_such_op"));
        let mut v = 0.0;
        assert_eq!(
            rssn_term_to_f64(s, 4_000_000_000, &raw mut v),
            RssnStatus::NoSuchTerm
        );
        assert_eq!(
            rssn_term_to_f64(std::ptr::null(), 0, &raw mut v),
            RssnStatus::NullArgument
        );
        assert_eq!(rssn_term_rational(s, 1, 0), RSSN_NO_TERM);
        rssn_session_free(s);
    }
}

#[test]
fn simulation_by_name() {
    let mut out = std::ptr::null_mut();
    let params = c(r#"{"width": 6, "height": 6, "temperature": 2.0, "mc_steps": 5, "seed": 1}"#);
    unsafe {
        assert_eq!(
            rssn_sim_run(c("ising").as_ptr(), params.as_ptr(), &raw mut out),
            RssnStatus::Ok
        );
        assert!(owned(out).contains("magnetization"));
        assert_eq!(
            rssn_sim_run(c("nope").as_ptr(), params.as_ptr(), &raw mut out),
            RssnStatus::Simulation
        );
    }
    assert!(last_error().contains("unknown scenario"));
    // SAFETY: a static NUL-terminated string.
    let version = unsafe { CStr::from_ptr(rssn_version()) };
    assert_eq!(version.to_str().unwrap(), env!("CARGO_PKG_VERSION"));
}

#[test]
fn batch_metadata_stays_with_its_request() {
    use std::mem::MaybeUninit;
    use std::ptr;

    let session = rssn_session_new();
    // SAFETY: the session is live and the provided arrays have three slots.
    unsafe {
        let mixed = rssn_term_parse(session, c("quantity(3, m) + quantity(4, s)").as_ptr());
        let length = rssn_term_parse(session, c("quantity(2, m)").as_ptr());
        let terms = [mixed, u32::MAX, length];
        let mut answers = MaybeUninit::<[RssnAnswer; 3]>::uninit();
        let mut statuses = MaybeUninit::<[RssnStatus; 3]>::uninit();
        let mut metadata = MaybeUninit::<[RssnAnswerMetadata; 3]>::uninit();
        assert_status(
            rssn_compute_all_with_metadata(
                session,
                terms.as_ptr(),
                terms.len(),
                ptr::null(),
                answers.as_mut_ptr().cast(),
                statuses.as_mut_ptr().cast(),
                ptr::null_mut(),
            ),
            RssnStatus::NullArgument,
        );
        assert_ok(rssn_compute_all_with_metadata(
            session,
            terms.as_ptr(),
            terms.len(),
            ptr::null(),
            answers.as_mut_ptr().cast(),
            statuses.as_mut_ptr().cast(),
            metadata.as_mut_ptr().cast(),
        ));
        assert_eq!(
            statuses.assume_init(),
            [RssnStatus::Ok, RssnStatus::NoSuchTerm, RssnStatus::Ok]
        );
        let [mixed, missing, length] = metadata.assume_init();
        let mut json = ptr::null_mut();
        assert_ok(rssn_compute_json(
            session,
            terms[0],
            ptr::null(),
            &raw mut json,
        ));
        let json: serde_json::Value = serde_json::from_str(&owned(json)).unwrap();
        assert_eq!(owned(mixed.dimension), json["dimension"].as_str().unwrap());
        assert!(
            owned(mixed.diagnostic_codes)
                .lines()
                .any(|code| code == "mismatch")
        );
        assert!(mixed.message.is_null());
        assert!(missing.dimension.is_null() && missing.diagnostic_codes.is_null());
        assert_ne!(owned(missing.message), "");
        assert_eq!(owned(length.dimension), "Length");
        assert!(length.diagnostic_codes.is_null() && length.message.is_null());
        rssn_session_free(session);
    }
}

#[test]
fn compute_all_matches_compute_and_reports_per_request_codes() {
    let s = rssn_session_new();
    unsafe {
        let x = rssn_term_sym(s, c("x").as_ptr());
        let sin = rssn_term_apply(s, c("sin").as_ptr(), &raw const x, 1);
        let d1 = rssn_term_apply(s, c("diff").as_ptr(), [sin, x].as_ptr(), 2);
        let cos = rssn_term_apply(s, c("cos").as_ptr(), &raw const x, 1);
        let d2 = rssn_term_apply(s, c("diff").as_ptr(), [cos, x].as_ptr(), 2);
        let handles = [d1, 999_999, d2];
        // C callers may reserve output buffers without initializing any entries.
        let mut answers = std::mem::MaybeUninit::<[RssnAnswer; 3]>::uninit();
        let mut statuses = std::mem::MaybeUninit::<[RssnStatus; 3]>::uninit();
        assert_eq!(
            rssn_compute_all(
                s,
                handles.as_ptr(),
                3,
                std::ptr::null(),
                answers.as_mut_ptr().cast::<RssnAnswer>(),
                statuses.as_mut_ptr().cast::<RssnStatus>()
            ),
            RssnStatus::Ok
        );
        // An Ok batch initializes both arrays, including any failed request.
        let answers = answers.assume_init();
        let statuses = statuses.assume_init();
        assert_eq!(
            statuses,
            [RssnStatus::Ok, RssnStatus::NoSuchTerm, RssnStatus::Ok]
        );
        assert_eq!(answers[1].term, RSSN_NO_TERM);
        for k in [0, 2] {
            let mut single = std::mem::zeroed::<RssnAnswer>();
            assert_eq!(
                rssn_compute(s, handles[k], std::ptr::null(), &raw mut single),
                RssnStatus::Ok
            );
            assert_eq!(answers[k].term, single.term);
            assert_eq!(answers[k].reduced, single.reduced);
        }
        // null arrays are only accepted for n = 0
        assert_eq!(
            rssn_compute_all(
                s,
                std::ptr::null(),
                0,
                std::ptr::null(),
                std::ptr::null_mut(),
                std::ptr::null_mut()
            ),
            RssnStatus::Ok
        );
        assert_eq!(
            rssn_compute_all(
                s,
                std::ptr::null(),
                1,
                std::ptr::null(),
                std::ptr::null_mut(),
                std::ptr::null_mut()
            ),
            RssnStatus::NullArgument
        );
        rssn_session_free(s);
    }
}

#[test]
fn dimension_c_api_and_json() {
    let s = rssn_session_new();
    unsafe {
        let config = rssn_config_new();
        let quantity = rssn_term_parse(s, c("quantity(2, m)").as_ptr());
        assert_ne!(quantity, RSSN_NO_TERM);
        let mut out = std::ptr::null_mut();
        assert_eq!(
            rssn_dimension(s, quantity, config, &raw mut out),
            RssnStatus::Ok
        );
        assert_eq!(owned(out), "Length");

        assert_eq!(
            rssn_compute_json(s, quantity, config, &raw mut out),
            RssnStatus::Ok
        );
        let document: serde_json::Value = serde_json::from_str(&owned(out)).unwrap();
        assert_eq!(document["dimension"], "Length");
        assert_eq!(document["term"], "quantity(2, m)");
        assert!(document["diagnostics"].is_array());

        // A declared symbol dimension reaches the C API.
        assert_eq!(
            rssn_config_dimension(config, c("v").as_ptr(), c("m/s").as_ptr()),
            RssnStatus::Ok
        );
        let v = rssn_term_parse(s, c("v").as_ptr());
        assert_eq!(rssn_dimension(s, v, config, &raw mut out), RssnStatus::Ok);
        assert_eq!(owned(out), "Length/Time");

        // Custom base dimension.
        assert_eq!(
            rssn_config_define_base_dimension(config, c("Currency").as_ptr()),
            RssnStatus::Ok
        );
        assert_eq!(
            rssn_config_dimension(config, c("price").as_ptr(), c("Currency").as_ptr()),
            RssnStatus::Ok
        );
        let price = rssn_term_parse(s, c("price").as_ptr());
        assert_eq!(
            rssn_dimension(s, price, config, &raw mut out),
            RssnStatus::Ok
        );
        assert_eq!(owned(out), "Currency");

        // Strict mode rejects a conflict; bypass turns it back into an answer.
        assert_eq!(rssn_config_dimension_mode(config, 2), RssnStatus::Ok);
        let bad = rssn_term_parse(s, c("quantity(1, m) + quantity(1, s)").as_ptr());
        let mut answer = std::mem::zeroed::<RssnAnswer>();
        assert_eq!(
            rssn_compute(s, bad, config, &raw mut answer),
            RssnStatus::Dimension
        );
        // The stable codes of the findings are exposed without parsing JSON.
        assert_eq!(rssn_dimension_count(), 1);
        let code = rssn_dimension_code(0);
        assert!(!code.is_null());
        // SAFETY: rssn keeps the code alive until the next computation, and
        // this whole test body is inside one `unsafe` block.
        assert_eq!(CStr::from_ptr(code).to_str().unwrap(), "mismatch");
        assert!(rssn_dimension_code(1).is_null(), "out of range is null");
        // An advisory run keeps the answer and reports its finding.
        assert_eq!(rssn_config_dimension_mode(config, 1), RssnStatus::Ok);
        assert_eq!(
            rssn_compute(s, bad, config, &raw mut answer),
            RssnStatus::Ok
        );
        assert!(rssn_dimension_count() >= 1);
        let code = rssn_dimension_code(0);
        assert!(!code.is_null());
        // SAFETY: as above.
        assert_eq!(CStr::from_ptr(code).to_str().unwrap(), "mismatch");
        // A computation without findings clears them, so a host must read
        // immediately after the call it cares about.
        let clean = rssn_term_parse(s, c("quantity(2, m)").as_ptr());
        assert_eq!(
            rssn_compute(s, clean, config, &raw mut answer),
            RssnStatus::Ok
        );
        assert_eq!(rssn_dimension_count(), 0);
        assert!(rssn_dimension_code(0).is_null());
        assert_eq!(rssn_config_bypass_dimensions(config), RssnStatus::Ok);
        assert_eq!(
            rssn_compute(s, bad, config, &raw mut answer),
            RssnStatus::Ok
        );
        assert_eq!(
            rssn_config_dimension_mode(config, 3),
            RssnStatus::InvalidArgument
        );
        rssn_config_free(config);
        rssn_session_free(s);
    }
}

/// The error matrix the module documentation promises: every failure path has
/// a precise [`RssnStatus`] and leaves a message in `rssn_last_error`; output
/// and array pointers are checked before handles are resolved; and a
/// successful call does not clear the previous message.
#[test]
fn status_matrix_maps_failures_and_sets_last_error() {
    let s = rssn_session_new();
    assert!(!s.is_null());
    let config = rssn_config_new();
    assert!(!config.is_null());
    unsafe {
        let x = rssn_term_sym(s, c("x").as_ptr());
        assert_ne!(x, RSSN_NO_TERM);
        let mut answer = std::mem::zeroed::<RssnAnswer>();
        let mut json = std::ptr::null_mut();
        let mut value = 0.0_f64;
        let mut tensor = std::mem::zeroed::<RssnTensor>();
        let mut bytes = std::ptr::null_mut();
        let mut bytes_len = usize::MAX;

        // Output and array pointers are validated before handles.
        assert_status(
            rssn_compute(s, x, std::ptr::null(), std::ptr::null_mut()),
            RssnStatus::NullArgument,
        );
        assert_status(
            rssn_compute_json(s, x, std::ptr::null(), std::ptr::null_mut()),
            RssnStatus::NullArgument,
        );
        assert_status(
            rssn_dimension(s, x, std::ptr::null(), std::ptr::null_mut()),
            RssnStatus::NullArgument,
        );
        assert_status(
            rssn_term_to_f64(s, 999_999, std::ptr::null_mut()),
            RssnStatus::NullArgument,
        );
        assert_status(
            rssn_term_tensor_data(s, x, std::ptr::null_mut()),
            RssnStatus::NullArgument,
        );
        assert_status(
            rssn_compute_bincode(
                s,
                x,
                std::ptr::null(),
                std::ptr::null_mut(),
                &raw mut bytes_len,
            ),
            RssnStatus::NullArgument,
        );
        assert_status(
            rssn_compute_bincode(s, x, std::ptr::null(), &raw mut bytes, std::ptr::null_mut()),
            RssnStatus::NullArgument,
        );
        assert_eq!(
            bytes_len,
            usize::MAX,
            "a failed call must not write the length"
        );
        assert!(bytes.is_null(), "a failed call must not write the buffer");

        // Null sessions and configurations.
        assert_status(
            rssn_compute(std::ptr::null(), 0, std::ptr::null(), &raw mut answer),
            RssnStatus::NullArgument,
        );
        assert_status(
            rssn_compute_bincode(
                std::ptr::null(),
                0,
                std::ptr::null(),
                &raw mut bytes,
                &raw mut bytes_len,
            ),
            RssnStatus::NullArgument,
        );
        assert_status(
            rssn_term_to_f64(std::ptr::null(), 0, &raw mut value),
            RssnStatus::NullArgument,
        );
        assert_status(
            rssn_config_symbolic(std::ptr::null_mut()),
            RssnStatus::NullArgument,
        );
        assert_status(
            rssn_config_numeric(std::ptr::null_mut(), 1e-9),
            RssnStatus::NullArgument,
        );
        assert_status(
            rssn_config_budget(std::ptr::null_mut(), 1, 1),
            RssnStatus::NullArgument,
        );
        assert_status(
            rssn_config_bypass_dimensions(std::ptr::null_mut()),
            RssnStatus::NullArgument,
        );
        // The null configuration is what the caller must fix first, even when
        // the mode is also invalid.
        assert_status(
            rssn_config_dimension_mode(std::ptr::null_mut(), 3),
            RssnStatus::NullArgument,
        );
        assert_status(
            rssn_config_dimension_mode(config, 3),
            RssnStatus::InvalidArgument,
        );
        assert_ok(rssn_config_dimension_mode(config, 1));

        // Unknown handles.
        assert_status(
            rssn_compute(s, 999_999, std::ptr::null(), &raw mut answer),
            RssnStatus::NoSuchTerm,
        );
        assert_status(
            rssn_compute_json(s, 999_999, std::ptr::null(), &raw mut json),
            RssnStatus::NoSuchTerm,
        );
        assert_status(
            rssn_dimension(s, 999_999, std::ptr::null(), &raw mut json),
            RssnStatus::NoSuchTerm,
        );
        assert_status(
            rssn_term_to_f64(s, 999_999, &raw mut value),
            RssnStatus::NoSuchTerm,
        );
        assert_status(
            rssn_term_tensor_data(s, 999_999, &raw mut tensor),
            RssnStatus::NoSuchTerm,
        );
        assert_status(
            rssn_compute_bincode(
                s,
                999_999,
                std::ptr::null(),
                &raw mut bytes,
                &raw mut bytes_len,
            ),
            RssnStatus::NoSuchTerm,
        );
        let rendered = rssn_term_to_string(s, 999_999);
        assert!(rendered.is_null());
        assert_ne!(last_error(), "", "a renderer reports its reason only here");

        // Malformed text and unknown operators.
        assert_no_term(rssn_term_parse(s, c("sin(").as_ptr()));
        assert_no_term(rssn_simplify(s, c("sin(").as_ptr()));
        assert_no_term(rssn_term_apply(
            s,
            c("frobnicate").as_ptr(),
            std::ptr::null(),
            0,
        ));
        assert!(last_error().contains("frobnicate"));
        assert_no_term(rssn_term_apply(s, c("sin").as_ptr(), std::ptr::null(), 1));
        assert_no_term(rssn_term_apply(
            s,
            c("sin").as_ptr(),
            std::ptr::null(),
            usize::MAX,
        ));
        let bad_args = [999_999_u32];
        assert_no_term(rssn_term_apply(s, c("sin").as_ptr(), bad_args.as_ptr(), 1));

        // Bad sizes, null arrays and non-numeric terms.
        assert_status(
            rssn_term_eval(s, x, std::ptr::null(), std::ptr::null(), 1, &raw mut value),
            RssnStatus::NullArgument,
        );
        let name = c("x");
        let names = [name.as_ptr(), std::ptr::null()];
        let bindings = [1.0_f64, 2.0];
        assert_status(
            rssn_term_eval(s, x, names.as_ptr(), bindings.as_ptr(), 2, &raw mut value),
            RssnStatus::NullArgument,
        );
        let symbolic = rssn_term_parse(s, c("sin(x)").as_ptr());
        assert_status(
            rssn_term_eval(
                s,
                symbolic,
                std::ptr::null(),
                std::ptr::null(),
                0,
                &raw mut value,
            ),
            RssnStatus::NotNumeric,
        );
        assert_status(
            rssn_term_to_f64(s, symbolic, &raw mut value),
            RssnStatus::NoValue,
        );
        assert_status(
            rssn_compute_all(
                s,
                std::ptr::null(),
                1,
                std::ptr::null(),
                std::ptr::null_mut(),
                std::ptr::null_mut(),
            ),
            RssnStatus::NullArgument,
        );
        assert_status(
            rssn_compute_all(
                s,
                std::ptr::null(),
                0,
                std::ptr::null(),
                std::ptr::null_mut(),
                std::ptr::null_mut(),
            ),
            RssnStatus::Ok,
        );

        // A numeric target without a binding is NotNumeric.
        let no_value = rssn_term_parse(s, c("x + 1").as_ptr());
        assert_ok(rssn_config_numeric(config, 1e-9));
        assert_status(
            rssn_compute(s, no_value, config, &raw mut answer),
            RssnStatus::NotNumeric,
        );
        assert_ok(rssn_config_symbolic(config));

        // Null names and strings that are not UTF-8.
        assert_no_term(rssn_term_sym(s, std::ptr::null()));
        assert_no_term(rssn_term_parse(s, std::ptr::null()));
        assert_status(
            rssn_config_bind(config, std::ptr::null(), 1.0),
            RssnStatus::NullArgument,
        );
        assert_status(
            rssn_config_assume(config, std::ptr::null(), 1),
            RssnStatus::NullArgument,
        );
        assert_status(
            rssn_config_dimension(config, std::ptr::null(), c("m").as_ptr()),
            RssnStatus::NullArgument,
        );
        assert_status(
            rssn_config_dimension(config, c("y").as_ptr(), std::ptr::null()),
            RssnStatus::NullArgument,
        );
        assert_status(
            rssn_config_define_base_dimension(config, std::ptr::null()),
            RssnStatus::NullArgument,
        );
        let invalid_utf8 = CString::new(vec![0xFF_u8]).unwrap_or_default();
        assert_no_term(rssn_term_sym(s, invalid_utf8.as_ptr()));
        assert_no_term(rssn_term_parse(s, invalid_utf8.as_ptr()));
        assert_status(
            rssn_config_bind(config, invalid_utf8.as_ptr(), 1.0),
            RssnStatus::InvalidUtf8,
        );

        // A simulation error (bad JSON, unknown scenario) is `Simulation`.
        let mut sim = std::ptr::null_mut();
        assert_status(
            rssn_sim_run(c("ising").as_ptr(), c("{").as_ptr(), &raw mut sim),
            RssnStatus::Simulation,
        );
        assert_status(
            rssn_sim_run(c("nope").as_ptr(), c("{}").as_ptr(), &raw mut sim),
            RssnStatus::Simulation,
        );
        assert_status(
            rssn_sim_run(c("ising").as_ptr(), std::ptr::null(), &raw mut sim),
            RssnStatus::NullArgument,
        );
        assert_status(
            rssn_sim_run(c("ising").as_ptr(), c("{}").as_ptr(), std::ptr::null_mut()),
            RssnStatus::NullArgument,
        );

        // A success keeps the previous failure message: rssn_last_error is
        // only written by a failing call, never cleared by a successful one.
        assert_no_term(rssn_term_parse(s, c("sin(").as_ptr()));
        let previous = last_error();
        let _ = rssn_manifest_hash();
        assert_eq!(last_error(), previous);

        rssn_config_free(config);
        rssn_session_free(s);
    }
}

/// `rssn_last_error` is per thread and starts empty.
#[test]
fn last_error_is_thread_local_and_starts_empty() {
    assert_eq!(last_error(), "");
    let spawned = std::thread::spawn(last_error);
    assert_eq!(spawned.join().unwrap_or_default(), "");
}

/// The bincode path carries the same record as the JSON path, and bytes that
/// are not that record are an error rather than a panic.
#[test]
fn compute_bincode_matches_the_json_path() {
    use crate::kernel::codec::decode_deterministic;

    let s = rssn_session_new();
    let config = rssn_config_new();
    assert!(!s.is_null());
    assert!(!config.is_null());
    unsafe {
        assert_ok(rssn_config_dimension(
            config,
            c("v").as_ptr(),
            c("m/s").as_ptr(),
        ));
        let term = rssn_term_parse(s, c("v").as_ptr());
        let (bytes, record) = bincode_record(s, term, config);
        assert_ne!(bytes, Vec::<u8>::new());
        let document = json_record(s, term, config);
        assert_eq!(record.term, document["term"].as_str().unwrap_or_default());
        assert_eq!(
            record.dimension,
            document["dimension"].as_str().map(str::to_owned)
        );
        assert_eq!(record.value, document["value"].as_f64());
        assert_eq!(
            record.reduced,
            document["reduced"].as_bool().unwrap_or(false)
        );
        assert_eq!(Some(record.iterations), document["iterations"].as_u64());
        let diagnostics: Vec<String> =
            document["diagnostics"]
                .as_array()
                .map_or_default(|findings| {
                    findings
                        .iter()
                        .filter_map(serde_json::Value::as_str)
                        .map(str::to_owned)
                        .collect()
                });
        assert_eq!(record.diagnostics, diagnostics);

        // A numeric target fills the Option fields.
        assert_ok(rssn_config_numeric(config, 1e-9));
        assert_ok(rssn_config_bind(config, c("v").as_ptr(), 3.0));
        let numeric = rssn_term_parse(s, c("v * 2").as_ptr());
        let (_, record) = bincode_record(s, numeric, config);
        assert_eq!(record.value, Some(6.0));
        assert!(record.error.is_some());
        let document = json_record(s, numeric, config);
        assert_eq!(record.value, document["value"].as_f64());

        // Bad bytes: empty, truncated, garbage and JSON are all errors.
        assert!(decode_deterministic::<RssnComputeRecord>(&[]).is_err());
        assert!(decode_deterministic::<RssnComputeRecord>(&[0xFF]).is_err());
        assert!(decode_deterministic::<RssnComputeRecord>(&bytes[..bytes.len() - 1]).is_err());
        assert!(decode_deterministic::<RssnComputeRecord>(b"{\"term\":\"x\"}").is_err());

        // Null is ignored by the free; no double free happens.
        rssn_bytes_free(std::ptr::null_mut(), 0);

        rssn_config_free(config);
        rssn_session_free(s);
    }
}

/// A dimension text that cannot be parsed is a declaration: it is a
/// diagnostic in the default advisory mode and a hard `Dimension` error in
/// strict mode, never a silently accepted dimension.
#[test]
fn a_bad_dimension_declaration_is_diagnosed_and_strict_rejects_it() {
    let s = rssn_session_new();
    let config = rssn_config_new();
    assert!(!s.is_null());
    assert!(!config.is_null());
    unsafe {
        assert_ok(rssn_config_dimension(
            config,
            c("z").as_ptr(),
            c("not_a_unit").as_ptr(),
        ));
        let term = rssn_term_parse(s, c("z + z").as_ptr());
        let mut answer = std::mem::zeroed::<RssnAnswer>();
        assert_ok(rssn_compute(s, term, config, &raw mut answer));
        let document = json_record(s, term, config);
        assert!(
            document["diagnostics"]
                .as_array()
                .is_some_and(|findings| !findings.is_empty()),
            "an unresolvable declaration must be a diagnostic: {document}"
        );
        assert_eq!(
            document["diagnostic_codes"].as_array().map(|codes| {
                codes
                    .iter()
                    .filter_map(serde_json::Value::as_str)
                    .map(str::to_owned)
                    .collect::<Vec<_>>()
            }),
            Some(vec!["unknown_unit".to_string()])
        );

        // Strict mode makes the declaration error blocking.
        assert_ok(rssn_config_dimension_mode(config, 2));
        assert_status(
            rssn_compute(s, term, config, &raw mut answer),
            RssnStatus::Dimension,
        );
        assert_ne!(last_error(), "");

        rssn_config_free(config);
        rssn_session_free(s);
    }
}

#[test]
fn diagnostic_codes_populated_in_compute_record() {
    let s = rssn_session_new();
    let config = rssn_config_new();
    assert!(!s.is_null());
    assert!(!config.is_null());
    unsafe {
        assert_ok(rssn_config_dimension(
            config,
            c("z").as_ptr(),
            c("not_a_unit").as_ptr(),
        ));
        let term = rssn_term_parse(s, c("z + z").as_ptr());
        let (_, record) = bincode_record(s, term, config);
        assert_ne!(record.diagnostic_codes, [] as [std::string::String; 0]);
        assert_eq!(record.diagnostic_codes, vec!["unknown_unit".to_string()]);
        for code in &record.diagnostic_codes {
            assert!(
                code.chars().all(|c| c.is_ascii_lowercase() || c == '_'),
                "`{code}` is not snake_case"
            );
        }

        let document = json_record(s, term, config);
        assert_eq!(
            document["diagnostic_codes"].as_array().map(|codes| {
                codes
                    .iter()
                    .filter_map(serde_json::Value::as_str)
                    .map(str::to_owned)
                    .collect::<Vec<_>>()
            }),
            Some(vec!["unknown_unit".to_string()])
        );

        rssn_config_free(config);
        rssn_session_free(s);
    }
}

/// The bytes of a batch answer and the record they decode to.
unsafe fn batch_record(
    session: *const RssnSession,
    terms: &[u32],
    config: *const RssnConfig,
) -> (Vec<u8>, RssnComputeBatchRecord) {
    let mut data = std::ptr::null_mut();
    let mut len = 0_usize;
    // SAFETY: the caller passes a live session, handles and configuration;
    // `terms.as_ptr()` is valid for `terms.len()` values (possibly empty).
    assert_eq!(
        unsafe {
            rssn_compute_all_bincode(
                session,
                terms.as_ptr(),
                terms.len(),
                config,
                &raw mut data,
                &raw mut len,
            )
        },
        RssnStatus::Ok,
        "{}",
        last_error()
    );
    // SAFETY: the call succeeded, so `data` holds `len` bytes.
    let bytes = unsafe { std::slice::from_raw_parts(data, len) }.to_vec();
    // SAFETY: returned by rssn and freed once.
    unsafe { rssn_bytes_free(data, len) };
    let record = crate::kernel::codec::decode_deterministic(&bytes)
        .unwrap_or_else(|e| panic!("the batch bytes rssn wrote must decode: {e}"));
    (bytes, record)
}

#[test]
fn compute_batches_keep_the_first_failure_message_in_input_order() {
    let session = rssn_session_new();
    let config = rssn_config_new();
    assert!(!session.is_null());
    assert!(!config.is_null());
    unsafe {
        assert_ok(rssn_config_numeric(config, 1e-9));
        let unbound = rssn_term_parse(session, c("y + 1").as_ptr());
        assert_ne!(unbound, RSSN_NO_TERM);
        let mut single = std::mem::MaybeUninit::<RssnAnswer>::uninit();
        assert_status(
            rssn_compute(session, RSSN_NO_TERM, config, single.as_mut_ptr()),
            RssnStatus::NoSuchTerm,
        );
        let missing = last_error();
        assert_status(
            rssn_compute(session, unbound, config, single.as_mut_ptr()),
            RssnStatus::NotNumeric,
        );
        let not_numeric = last_error();
        assert_ne!(missing, not_numeric);

        let mut terms = [RSSN_NO_TERM, unbound];
        let mut expected_statuses = [RssnStatus::NoSuchTerm, RssnStatus::NotNumeric];
        for expected_message in [missing, not_numeric] {
            let mut answers = [const { std::mem::MaybeUninit::<RssnAnswer>::uninit() }; 2];
            let mut statuses = [RssnStatus::Panic; 2];
            assert_ok(rssn_compute_all(
                session,
                terms.as_ptr(),
                terms.len(),
                config,
                answers.as_mut_ptr().cast::<RssnAnswer>(),
                statuses.as_mut_ptr(),
            ));
            assert_eq!(statuses, expected_statuses);
            assert_eq!(last_error(), expected_message, "C array batch");

            let (_, batch) = batch_record(session, &terms, config);
            assert_eq!(batch.entries.len(), terms.len());
            for (entry, status) in batch.entries.iter().zip(expected_statuses) {
                assert_eq!(entry.status, status as u32);
                assert!(entry.answer.is_none());
            }
            assert_eq!(last_error(), expected_message, "bincode batch");
            terms.reverse();
            expected_statuses.reverse();
        }
        rssn_config_free(config);
        rssn_session_free(session);
    }
}

/// The batch bincode path carries one entry per input handle, in input order,
/// with the same per-request statuses as `rssn_compute_all` and the same
/// answers as the JSON and single bincode paths; bad bytes and bad arguments
/// are errors, not panics.
#[test]
fn compute_all_bincode_matches_json_and_reports_per_request_status() {
    use crate::kernel::codec::decode_deterministic;

    let s = rssn_session_new();
    let config = rssn_config_new();
    assert!(!s.is_null());
    assert!(!config.is_null());
    unsafe {
        assert_ok(rssn_config_dimension(
            config,
            c("v").as_ptr(),
            c("m/s").as_ptr(),
        ));
        let v = rssn_term_parse(s, c("v").as_ptr());
        let sum = rssn_term_parse(s, c("v + 1").as_ptr());
        let unbound = rssn_term_parse(s, c("y + 1").as_ptr());
        let handles = [v, 999_999, sum];

        let (bytes, batch) = batch_record(s, &handles, config);
        assert_ne!(bytes, Vec::<u8>::new());
        assert_eq!(batch.entries.len(), 3);
        assert_eq!(batch.entries[0].status, RssnStatus::Ok as u32);
        assert_eq!(batch.entries[1].status, RssnStatus::NoSuchTerm as u32);
        assert!(batch.entries[1].answer.is_none());
        assert_eq!(batch.entries[2].status, RssnStatus::Ok as u32);
        assert_ne!(last_error(), "", "a failing slot leaves its message");

        // The C-array batch reports the same codes and answers.
        let mut answers = [std::mem::zeroed::<RssnAnswer>(); 3];
        let mut statuses = [RssnStatus::Panic; 3];
        assert_ok(rssn_compute_all(
            s,
            handles.as_ptr(),
            3,
            config,
            answers.as_mut_ptr(),
            statuses.as_mut_ptr(),
        ));
        for k in [0, 2] {
            assert_eq!(statuses[k], RssnStatus::Ok);
            let document = json_record(s, handles[k], config);
            let (_, single) = bincode_record(s, handles[k], config);
            let entry = batch.entries[k]
                .answer
                .clone()
                .unwrap_or_else(|| panic!("entry {k} must be Ok"));
            assert_eq!(entry.term, document["term"].as_str().unwrap_or_default());
            assert_eq!(entry.term, single.term);
            assert_eq!(
                entry.dimension,
                document["dimension"].as_str().map(str::to_owned)
            );
            assert_eq!(entry.value, document["value"].as_f64());
            assert_eq!(
                entry.reduced,
                document["reduced"].as_bool().unwrap_or(false)
            );
            assert_eq!(entry.term, owned(rssn_term_to_string(s, answers[k].term)));
        }

        // A per-request computation failure travels inside the record: the
        // numeric phase cannot evaluate the unbound term.
        assert_ok(rssn_config_numeric(config, 1e-9));
        assert_ok(rssn_config_bind(config, c("v").as_ptr(), 3.0));
        let numeric = [sum, unbound];
        let (_, batch) = batch_record(s, &numeric, config);
        assert_eq!(batch.entries[0].status, RssnStatus::Ok as u32);
        assert_eq!(
            batch.entries[0]
                .answer
                .as_ref()
                .and_then(|entry| entry.value),
            Some(4.0)
        );
        assert_eq!(batch.entries[1].status, RssnStatus::NotNumeric as u32);
        assert!(batch.entries[1].answer.is_none());
        assert_ne!(last_error(), "");

        // An empty batch is well formed and encodes as an empty entry list.
        let (empty_bytes, empty) = batch_record(s, &[], config);
        assert!(empty.entries.is_empty());
        assert_ne!(empty_bytes, Vec::<u8>::new());

        // Bad arguments: precise codes, nothing written.
        let mut data = std::ptr::null_mut();
        let mut len = usize::MAX;
        assert_status(
            rssn_compute_all_bincode(
                s,
                handles.as_ptr(),
                3,
                config,
                std::ptr::null_mut(),
                &raw mut len,
            ),
            RssnStatus::NullArgument,
        );
        assert_eq!(len, usize::MAX);
        assert_status(
            rssn_compute_all_bincode(
                s,
                handles.as_ptr(),
                3,
                config,
                &raw mut data,
                std::ptr::null_mut(),
            ),
            RssnStatus::NullArgument,
        );
        assert!(data.is_null());
        assert_status(
            rssn_compute_all_bincode(s, std::ptr::null(), 1, config, &raw mut data, &raw mut len),
            RssnStatus::NullArgument,
        );
        assert_status(
            rssn_compute_all_bincode(
                std::ptr::null(),
                std::ptr::null(),
                0,
                config,
                &raw mut data,
                &raw mut len,
            ),
            RssnStatus::NullArgument,
        );
        assert!(data.is_null());
        assert_eq!(len, usize::MAX);

        // Bad bytes are errors, never panics. Note that a bincode payload is
        // not self-framing: a byte string cut at a sequence boundary can still
        // decode as a *shorter* batch, so a host must treat `*out_len` as the
        // exact length (and check a digest for untrusted storage) rather than
        // relying on the decoder to notice truncation. The cases below are
        // structurally invalid, not merely short.
        assert!(decode_deterministic::<RssnComputeBatchRecord>(&[]).is_err());
        assert!(decode_deterministic::<RssnComputeBatchRecord>(&[0xFF]).is_err());

        rssn_config_free(config);
        rssn_session_free(s);
    }
}

// ---------------------------------------------------------------------------
// host callbacks
// ---------------------------------------------------------------------------

/// The callback the tests below register: `x -> x*x - 2`, the residual whose
/// positive root is `sqrt(2)`.
extern "C" fn square_minus_two(
    _user: *mut std::ffi::c_void,
    session: *const RssnSession,
    args: *const u32,
    count: usize,
    out: *mut u32,
) -> RssnStatus {
    if session.is_null() || args.is_null() || out.is_null() || count != 1 {
        return fail(
            RssnStatus::InvalidArgument,
            "square_minus_two wants exactly one argument",
        );
    }
    let mut x = 0.0;
    // SAFETY: the callback contract: `count` live argument cells, one session.
    let status = unsafe { rssn_value_to_f64(session, *args, &raw mut x) };
    if status != RssnStatus::Ok {
        return status;
    }
    // SAFETY: as above; `out` is a valid pointer by the check above.
    unsafe { rssn_value_new_float(session, x.mul_add(x, -2.0), out) }
}

/// A callback that always fails, to exercise the error path.
extern "C" fn refuses(
    _user: *mut std::ffi::c_void,
    _session: *const RssnSession,
    _args: *const u32,
    _count: usize,
    _out: *mut u32,
) -> RssnStatus {
    let message = c("the host refuses to evaluate");
    // SAFETY: a NUL-terminated string that outlives the call.
    unsafe { rssn_callback_fail(message.as_ptr()) };
    RssnStatus::InvalidArgument
}

/// Stores `value` in a fresh cell of `session`.
unsafe fn cell_of(
    session: *const RssnSession,
    value: f64,
) -> u32 {
    let mut cell = 0u32;
    // SAFETY: the caller passes a live session.
    assert_ok(unsafe { rssn_value_new_float(session, value, &raw mut cell) });
    cell
}

/// Stores `value` in a fresh cell of `session`.
unsafe fn int_cell_of(
    session: *const RssnSession,
    value: i64,
) -> u32 {
    let mut cell = 0u32;
    // SAFETY: the caller passes a live session.
    assert_ok(unsafe { rssn_value_new_int(session, value, &raw mut cell) });
    cell
}

#[test]
fn dense_value_cells_validate_rank_zero_and_extent_products() {
    let session = rssn_session_new();
    assert!(!session.is_null());
    // SAFETY: the session, extent storage and out pointers stay live. No data
    // pointer is needed for an empty tensor.
    unsafe {
        let empty_tensor = |shape: &[usize], out: *mut u32| {
            rssn_value_new_tensor(
                session,
                std::ptr::null(),
                0,
                shape.as_ptr(),
                shape.len(),
                out,
            )
        };
        let mut cell = 0;
        assert_eq!(
            empty_tensor(&[], &raw mut cell),
            RssnStatus::InvalidArgument
        );
        assert_eq!(
            empty_tensor(&[usize::MAX, 2], &raw mut cell),
            RssnStatus::InvalidArgument
        );
        assert_eq!(cell, 0);
        assert_ok(empty_tensor(&[2, 0, 3], &raw mut cell));
        let mut tensor = RssnTensor {
            shape: std::ptr::null_mut(),
            rank: 0,
            len: 0,
            data: std::ptr::null_mut(),
        };
        assert_ok(rssn_value_to_tensor(session, cell, &raw mut tensor));
        assert_eq!(tensor.len, 0);
        assert_eq!(
            std::slice::from_raw_parts(tensor.shape, tensor.rank),
            &[2, 0, 3]
        );
        rssn_tensor_free(&raw mut tensor);
        rssn_value_free(session, cell);
        rssn_session_free(session);
    }
}

#[test]
fn decimal_integer_cells_preserve_sign_and_every_digit() {
    let session = rssn_session_new();
    assert!(!session.is_null());
    let digits = "-1267650600228229401496703205383";
    let text = c(digits);
    let invalid = c("1e7");
    // SAFETY: the session, input C strings and output storage remain live.
    unsafe {
        let mut cell = 0;
        assert_ok(rssn_value_new_int_decimal(
            session,
            text.as_ptr(),
            &raw mut cell,
        ));
        let mut result = std::ptr::null_mut();
        assert_ok(rssn_value_to_text(session, cell, &raw mut result));
        assert_eq!(owned(result), digits);
        rssn_value_free(session, cell);
        assert_status(
            rssn_value_new_int_decimal(session, invalid.as_ptr(), &raw mut cell),
            RssnStatus::InvalidArgument,
        );
        assert_eq!(rssn_value_count(session), 0);
        rssn_session_free(session);
    }
}

#[test]
fn cloned_cells_have_independent_lifetimes() {
    let session = rssn_session_new();
    assert!(!session.is_null());
    // SAFETY: the session and all output pointers are live for the calls.
    unsafe {
        let original = int_cell_of(session, 42);
        let mut copy = 0;
        assert_ok(rssn_value_clone(session, original, &raw mut copy));
        assert_ne!(copy, original);
        rssn_value_free(session, original);
        let mut number = 0.0;
        assert_ok(rssn_value_to_f64(session, copy, &raw mut number));
        assert_eq!(number.to_bits(), 42.0_f64.to_bits());
        let mut missing = 0;
        assert_eq!(
            rssn_value_clone(session, original, &raw mut missing),
            RssnStatus::NoValue
        );
        assert_eq!(missing, 0);
        rssn_value_free(session, copy);
        assert_eq!(rssn_value_count(session), 0);
        rssn_session_free(session);
    }
}

#[test]
fn a_host_callback_reaches_a_higher_order_operator() {
    let session = rssn_session_new();
    assert!(!session.is_null(), "{}", last_error());

    // The signature comes from the operator's own declaration:
    // `num_brent(a, b, tolerance, iterations, f)`.
    let mut callback = 0u32;
    let op = c("num_brent");
    assert_ok(unsafe {
        rssn_callback_new_for(
            session,
            op.as_ptr(),
            4,
            square_minus_two,
            std::ptr::null_mut(),
            &raw mut callback,
        )
    });
    assert_ne!(callback, 0);
    assert_eq!(unsafe { rssn_callback_count(session) }, 1);

    // SAFETY: a live session; the cells stay alive for the call.
    let cells = unsafe {
        [
            cell_of(session, 1.0),
            cell_of(session, 2.0),
            cell_of(session, 1e-12),
            int_cell_of(session, 64),
            callback,
        ]
    };
    let mut answer = 0u32;
    // SAFETY: `cells` holds five live ids and `answer` is a valid pointer.
    assert_ok(unsafe {
        rssn_value_call(
            session,
            op.as_ptr(),
            cells.as_ptr(),
            cells.len(),
            &raw mut answer,
        )
    });

    // `num_brent` answers a record; its first field is the root.
    let mut root_cell = 0u32;
    // SAFETY: a live record cell and a valid out pointer.
    assert_ok(unsafe { rssn_value_at(session, answer, 0, &raw mut root_cell) });
    let mut root = 0.0;
    // SAFETY: as above; `_to_f64` looks through the value's dimension wrapper.
    assert_ok(unsafe { rssn_value_to_f64(session, root_cell, &raw mut root) });
    assert!(
        (root - std::f64::consts::SQRT_2).abs() < 1e-9,
        "root = {root}"
    );

    // Retiring the callback empties the bank.
    unsafe { rssn_callback_free(session, callback) };
    assert_eq!(unsafe { rssn_callback_count(session) }, 0);
    // SAFETY: the session is still live.
    unsafe { rssn_session_free(session) };
}

#[test]
fn an_explicit_signature_is_enough_to_register_a_callback() {
    let session = rssn_session_new();
    assert!(!session.is_null(), "{}", last_error());
    let signature = c(r#"{"params": ["real"], "ret": "real"}"#);
    let mut callback = 0u32;
    // SAFETY: a live session, a NUL-terminated signature, a valid out pointer.
    assert_ok(unsafe {
        rssn_callback_new(
            session,
            signature.as_ptr(),
            square_minus_two,
            std::ptr::null_mut(),
            &raw mut callback,
        )
    });
    assert_ne!(callback, 0);

    // A malformed signature is a status, never a panic.
    let bad = c(r#"{"params": "real"}"#);
    let mut ignored = 0u32;
    // SAFETY: as above.
    assert_status(
        unsafe {
            rssn_callback_new(
                session,
                bad.as_ptr(),
                square_minus_two,
                std::ptr::null_mut(),
                &raw mut ignored,
            )
        },
        RssnStatus::InvalidArgument,
    );
    // SAFETY: the session is still live.
    unsafe { rssn_session_free(session) };
}

#[test]
fn callback_signatures_are_listed_and_failures_propagate() {
    let session = rssn_session_new();
    assert!(!session.is_null(), "{}", last_error());
    let op = c("num_brent");
    let mut raw = std::ptr::null_mut();
    // SAFETY: a live session, a NUL-terminated operator name, a valid pointer.
    assert_ok(unsafe { rssn_callback_signatures(session, op.as_ptr(), 4, &raw mut raw) });
    let listed: serde_json::Value = serde_json::from_str(&owned(raw)).expect("the list is JSON");
    let first = &listed[0];
    assert_eq!(first["params"][0], serde_json::json!("float"));
    assert_eq!(first["ret"], serde_json::json!("float"));

    // An operator with no function at that index lists nothing.
    let integrate = c("num_integrate");
    // SAFETY: as above.
    assert_ok(unsafe { rssn_callback_signatures(session, integrate.as_ptr(), 0, &raw mut raw) });
    assert_eq!(owned(raw), "[]");

    // A callback that fails aborts the operator with the host's message.
    let mut callback = 0u32;
    // SAFETY: as above.
    assert_ok(unsafe {
        rssn_callback_new_for(
            session,
            op.as_ptr(),
            4,
            refuses,
            std::ptr::null_mut(),
            &raw mut callback,
        )
    });
    // SAFETY: a live session; the cells stay alive for the call.
    let cells = unsafe {
        [
            cell_of(session, 1.0),
            cell_of(session, 2.0),
            cell_of(session, 1e-12),
            int_cell_of(session, 8),
            callback,
        ]
    };
    let mut answer = 0u32;
    // SAFETY: five live ids and a valid out pointer.
    assert_status(
        unsafe {
            rssn_value_call(
                session,
                op.as_ptr(),
                cells.as_ptr(),
                cells.len(),
                &raw mut answer,
            )
        },
        RssnStatus::InvalidArgument,
    );
    assert!(
        last_error().contains("refuses"),
        "the host message crosses: {}",
        last_error()
    );
    // SAFETY: the session is still live.
    unsafe { rssn_session_free(session) };
}

#[test]
fn coherent_dimension_value_failures_leave_output_and_cell_bank_unchanged() {
    use super::value::MAX_COHERENT_DIMENSION_JSON_BYTES;

    const SENTINEL: u32 = u32::MAX;
    let session = rssn_session_new();
    let mut input = 0;
    let mut integer = 0;
    // SAFETY: session and all pointers are valid for each call.
    unsafe {
        assert_ok(rssn_value_new_float(session, -0.0, &raw mut input));
        assert_ok(rssn_value_new_int(session, 0, &raw mut integer));
        let count = rssn_value_count(session);
        let check = |cell, bytes: &[u8], expected| {
            let mut out = SENTINEL;
            assert_status(
                rssn_value_with_coherent_dimension(
                    session,
                    cell,
                    bytes.as_ptr(),
                    bytes.len(),
                    &raw mut out,
                ),
                expected,
            );
            assert_eq!(out, SENTINEL);
            assert_eq!(rssn_value_count(session), count);
        };
        for bytes in [
            &b""[..],
            br#""s""#,
            br#""m/N/s""#,
            br#""any""#,
            br#""non_metric""#,
            br#""container""#,
            br#"{"var":"D"}"#,
            br#"{"pow":[{"base":"Length"},{"symbol":"n"}]}"#,
            br#"{"pow":[{"base":"Length"},{"rat":{"numer":1,"denom":0}}]}"#,
            b"{\"base\":\"Length\"}\0",
        ] {
            check(input, bytes, RssnStatus::InvalidArgument);
        }
        check(input, &[0xff], RssnStatus::InvalidUtf8);
        check(
            integer,
            br#"{"base":"Length"}"#,
            RssnStatus::InvalidArgument,
        );
        check(SENTINEL, br#"{"base":"Length"}"#, RssnStatus::NoValue);

        let mut boundary = br#"{"base":"Length"}"#.to_vec();
        boundary.resize(MAX_COHERENT_DIMENSION_JSON_BYTES + 1, b' ');
        check(input, &boundary, RssnStatus::InvalidArgument);
        let mut output = SENTINEL;
        assert_status(
            rssn_value_with_coherent_dimension(
                session,
                input,
                std::ptr::null(),
                0,
                &raw mut output,
            ),
            RssnStatus::NullArgument,
        );
        assert_eq!(output, SENTINEL);
        assert_status(
            rssn_value_with_coherent_dimension(
                session,
                input,
                boundary.as_ptr(),
                MAX_COHERENT_DIMENSION_JSON_BYTES,
                std::ptr::null_mut(),
            ),
            RssnStatus::NullArgument,
        );
        assert_eq!(rssn_value_count(session), count);

        assert_ok(rssn_value_with_coherent_dimension(
            session,
            input,
            boundary.as_ptr(),
            MAX_COHERENT_DIMENSION_JSON_BYTES,
            &raw mut output,
        ));
        assert_ne!(output, input);
        assert_ne!(output, SENTINEL);
        assert_eq!(rssn_value_count(session), count + 1);
        let mut kind = u8::MAX;
        assert_ok(rssn_value_kind(session, output, &raw mut kind));
        assert_eq!(kind, 2);
        let mut magnitude = 1.0;
        assert_ok(rssn_value_to_f64(session, output, &raw mut magnitude));
        assert_f64_is_negative_zero(magnitude);
        let mut wire = std::ptr::null_mut();
        assert_ok(rssn_value_json(session, output, &raw mut wire));
        let json: serde_json::Value = serde_json::from_str(&owned(wire)).unwrap();
        assert_eq!(json["dyn"]["dim"], serde_json::json!({"base": "Length"}));
        rssn_value_free(session, output);
        assert_eq!(rssn_value_count(session), count);
        rssn_session_free(session);
    }
}

fn assert_f64_is_negative_zero(value: f64) {
    assert_eq!(value.to_bits(), (-0.0_f64).to_bits());
}

#[test]
fn coherent_dimension_value_preserves_rank_shape_and_rejects_retagging() {
    let session = rssn_session_new();
    let length = br#"{"base":"Length"}"#;
    let time = br#"{"base":"Time"}"#;
    // SAFETY: session and all pointers are valid for each call.
    unsafe {
        for (shape, magnitude) in [
            (&[][..], &[1.25][..]),
            (&[0, 3][..], &[][..]),
            (&[3, 0][..], &[][..]),
            (&[2, 2][..], &[1.0, -0.0, 3.0, 4.0][..]),
        ] {
            let mut plain = 0;
            assert_ok(rssn_value_new_tensor(
                session,
                magnitude.as_ptr(),
                magnitude.len(),
                shape.as_ptr(),
                shape.len(),
                &raw mut plain,
            ));
            let mut measured = 0;
            assert_ok(rssn_value_with_coherent_dimension(
                session,
                plain,
                length.as_ptr(),
                length.len(),
                &raw mut measured,
            ));
            let mut wire = std::ptr::null_mut();
            assert_ok(rssn_value_json(session, measured, &raw mut wire));
            let json: serde_json::Value = serde_json::from_str(&owned(wire)).unwrap();
            let tensor = &json["tensor"];
            assert_eq!(tensor["shape"], serde_json::json!(shape));
            assert_eq!(
                tensor["elem"],
                serde_json::json!({"dim": ["float", {"base": "Length"}]})
            );
            assert_eq!(tensor["data"].as_array().unwrap().len(), magnitude.len());
            let count = rssn_value_count(session);
            let mut conflicted = u32::MAX;
            assert_status(
                rssn_value_with_coherent_dimension(
                    session,
                    measured,
                    time.as_ptr(),
                    time.len(),
                    &raw mut conflicted,
                ),
                RssnStatus::InvalidArgument,
            );
            assert_eq!(conflicted, u32::MAX);
            assert_eq!(rssn_value_count(session), count);
            if let Some(&first) = magnitude.first() {
                let mut element = 0;
                assert_ok(rssn_value_at(session, measured, 0, &raw mut element));
                let mut element_wire = std::ptr::null_mut();
                assert_ok(rssn_value_json(session, element, &raw mut element_wire));
                let element_json: serde_json::Value =
                    serde_json::from_str(&owned(element_wire)).unwrap();
                assert_eq!(
                    element_json["dyn"]["dim"],
                    serde_json::json!({"base":"Length"})
                );
                assert_eq!(element_json["dyn"]["value"]["float"].as_f64(), Some(first));
                rssn_value_free(session, element);
            }
            rssn_value_free(session, measured);
            rssn_value_free(session, plain);
            assert_eq!(rssn_value_count(session), 0);
        }
        rssn_session_free(session);
    }
}
