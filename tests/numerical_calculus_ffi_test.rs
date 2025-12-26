use assert_approx_eq::assert_approx_eq;
use rssn::ffi_apis::common::{rssn_free_bincode_buffer, rssn_free_string};
use rssn::ffi_apis::numerical_calculus_ffi::*;
use rssn::symbolic::core::Expr;
use std::ffi::{CStr, CString};

#[test]

fn test_numerical_calculus_handle_ffi() {

    unsafe {

        let x = Expr::new_variable("x");

        let two = Expr::new_constant(2.0);

        let f = Expr::new_pow(x, two);

        let mut res = 0.0;

        let var_name = CString::new("x").unwrap();

        let status =
            handle::rssn_num_calculus_partial_derivative(&f, var_name.as_ptr(), 2.0, &mut res);

        assert_eq!(status, 0);

        assert_approx_eq!(res, 4.0, 1e-5f64);
    }
}

#[test]

fn test_numerical_gradient_json_ffi() {

    unsafe {

        let x = Expr::new_variable("x");

        let y = Expr::new_variable("y");

        let f = Expr::new_add(Expr::new_pow(x, Expr::new_constant(2.0)), y);

        let f_json = serde_json::to_string(&f).unwrap();

        let json_input = format!(
            r#"{{"expr": {}, "vars": ["x", "y"], "point": [2.0, 3.0]}}"#,
            f_json
        );

        let c_json = CString::new(json_input).unwrap();

        let res_ptr = json::rssn_numerical_gradient_json(c_json.as_ptr());

        assert!(!res_ptr.is_null());

        let res_str = CStr::from_ptr(res_ptr).to_str().unwrap();

        let v: serde_json::Value =
            serde_json::from_str(res_str).expect("Failed to parse result JSON");

        let ok = v["ok"].as_array().expect("Result 'ok' should be an array");

        assert_approx_eq!(ok[0].as_f64().unwrap(), 4.0, 1e-5f64);

        assert_approx_eq!(ok[1].as_f64().unwrap(), 1.0, 1e-5f64);

        rssn_free_string(res_ptr);
    }
}

#[test]

fn test_numerical_hessian_bincode_ffi() {

    unsafe {

        use rssn::ffi_apis::common::{from_bincode_buffer, to_bincode_buffer};
        use serde::{Deserialize, Serialize};

        #[derive(Serialize)]

        struct HessianInput {
            expr: Expr,
            vars: Vec<String>,
            point: Vec<f64>,
        }

        let x = Expr::new_variable("x");

        let f = Expr::new_pow(x, Expr::new_constant(2.0));

        let input = HessianInput {
            expr: f,
            vars: vec!["x".to_string()],
            point: vec![2.0],
        };

        let buffer = to_bincode_buffer(&input);

        let res_buffer = bincode_api::rssn_numerical_hessian_bincode(buffer);

        assert!(!res_buffer.is_null());

        #[derive(Deserialize)]

        struct FfiResult<T, E> {
            ok: Option<T>,
            #[allow(dead_code)]
            err: Option<E>,
        }

        let res: FfiResult<Vec<Vec<f64>>, String> =
            from_bincode_buffer(&res_buffer).expect("Failed to decode bincode result");

        assert!(res.err.is_none());

        assert_approx_eq!(res.ok.unwrap()[0][0], 2.0, 1e-4f64);

        rssn_free_bincode_buffer(res_buffer);

        rssn_free_bincode_buffer(buffer);
    }
}
