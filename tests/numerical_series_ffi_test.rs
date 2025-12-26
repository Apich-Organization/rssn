use assert_approx_eq::assert_approx_eq;
use rssn::ffi_apis::common::{rssn_free_bincode_buffer, rssn_free_string};
use rssn::ffi_apis::numerical_series_ffi::*;
use rssn::symbolic::core::Expr;
use std::ffi::{CStr, CString};

#[test]

fn test_numerical_series_handle_ffi() {

    unsafe {

        let x = Expr::new_variable("x");

        let f = Expr::new_pow(
            x,
            Expr::new_constant(2.0),
        );

        let var_name = CString::new("x").unwrap();

        let coeffs_ptr = handle::rssn_numerical_taylor_coefficients(
            &f,
            var_name.as_ptr(),
            0.0,
            2,
        );

        assert!(!coeffs_ptr.is_null());

        let coeffs = &*coeffs_ptr;

        assert_approx_eq!(coeffs[2], 1.0, 1e-10f64);

        let val = handle::rssn_numerical_evaluate_power_series(coeffs_ptr, 0.0, 1.0);

        assert_approx_eq!(val, 1.0, 1e-10f64);

        // Freeing the Vec pointer
        let _ = Box::from_raw(coeffs_ptr);
    }
}

#[test]

fn test_numerical_sum_series_json_ffi() {

    unsafe {

        let n = Expr::new_variable("n");

        let f = n;

        let f_json = serde_json::to_string(&f).unwrap();

        let json_input = format!(
            r#"{{"expr": {}, "var": "n", "start": 1, "end": 10}}"#,
            f_json
        );

        let c_json = CString::new(json_input).unwrap();

        let res_ptr = json::rssn_numerical_sum_series_json(c_json.as_ptr());

        assert!(!res_ptr.is_null());

        let res_str = CStr::from_ptr(res_ptr)
            .to_str()
            .unwrap();

        let v: serde_json::Value =
            serde_json::from_str(res_str).expect("Failed to parse result JSON");

        assert_approx_eq!(
            v["ok"]
                .as_f64()
                .unwrap(),
            55.0,
            1e-10f64
        );

        rssn_free_string(res_ptr);
    }
}

#[test]

fn test_numerical_taylor_bincode_ffi() {

    unsafe {

        use rssn::ffi_apis::common::{from_bincode_buffer, to_bincode_buffer};
        use serde::{Deserialize, Serialize};

        #[derive(Serialize)]

        struct TaylorInput {
            expr: Expr,
            var: String,
            at_point: f64,
            order: usize,
        }

        let x = Expr::new_variable("x");

        let f = Expr::new_pow(
            x,
            Expr::new_constant(2.0),
        );

        let input = TaylorInput {
            expr: f,
            var: "x".to_string(),
            at_point: 0.0,
            order: 2,
        };

        let buffer = to_bincode_buffer(&input);

        let res_buffer = bincode_api::rssn_numerical_taylor_coefficients_bincode(buffer);

        assert!(!res_buffer.is_null());

        #[derive(Deserialize)]

        struct FfiResult<T, E> {
            ok: Option<T>,
            #[allow(dead_code)]
            err: Option<E>,
        }

        let res: FfiResult<Vec<f64>, String> =
            from_bincode_buffer(&res_buffer).expect("Failed to decode bincode result");

        assert!(res.err.is_none());

        assert_approx_eq!(
            res.ok.unwrap()[2],
            1.0,
            1e-10f64
        );

        rssn_free_bincode_buffer(res_buffer);

        rssn_free_bincode_buffer(buffer);
    }
}
