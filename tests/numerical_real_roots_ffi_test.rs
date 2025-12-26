use rssn::ffi_apis::numerical_real_roots_ffi::*;
use rssn::prelude::numerical::*;
use std::ffi::{CStr, CString};

#[test]

fn test_find_roots_handle_ffi() {

    unsafe {

        // x^2 - 1 = 0 => roots -1, 1
        let coeffs = vec![1.0, 0.0, -1.0];

        let roots_ptr = handle::rssn_real_roots_find_roots(coeffs.as_ptr(), coeffs.len(), 1e-9);

        assert!(!roots_ptr.is_null());

        let len = handle::rssn_real_roots_get_vec_len(roots_ptr);

        assert_eq!(len, 2);

        let mut out = vec![0.0; len];

        handle::rssn_real_roots_get_vec_data(roots_ptr, out.as_mut_ptr());

        // Sorted: -1.0, 1.0
        assert!((out[0] - (-1.0)).abs() < 1e-9);

        assert!((out[1] - 1.0).abs() < 1e-9);

        handle::rssn_real_roots_free_vec(roots_ptr);
    }
}

#[test]

fn test_find_roots_json_ffi() {

    unsafe {

        // x^2 - 1 = 0
        let json_input = r#"{"coeffs": [1.0, 0.0, -1.0], "tolerance": 1e-9}"#;

        let c_json = CString::new(json_input).unwrap();

        let res_ptr = json::rssn_real_roots_find_roots_json(c_json.as_ptr());

        let res_str = CStr::from_ptr(res_ptr).to_str().unwrap();

        let v: serde_json::Value = serde_json::from_str(res_str).unwrap();

        // FfiResult with Option<Vec<f64>>
        let ok = &v["ok"];

        assert!(ok.is_array());

        let roots = ok.as_array().unwrap();

        assert_eq!(roots.len(), 2);

        assert!((roots[0].as_f64().unwrap() - (-1.0)).abs() < 1e-9);

        let _ = CString::from_raw(res_ptr);
    }
}
