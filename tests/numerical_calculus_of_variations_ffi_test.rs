use assert_approx_eq::assert_approx_eq;
use rssn::ffi_apis::common::{rssn_free_bincode_buffer, rssn_free_string};
use rssn::ffi_apis::numerical_calculus_of_variations_ffi::{bincode_api, handle, json};
use rssn::symbolic::core::Expr;
use std::ffi::{CStr, CString};

#[test]

fn test_cov_handle_ffi() {

    unsafe {

        let t = Expr::new_variable("t");

        let y_dot = Expr::new_variable("y_dot");

        let lagrangian = Expr::new_mul(
            Expr::new_constant(0.5),
            Expr::new_pow(y_dot, Expr::new_constant(2.0)),
        );

        let path = Expr::new_mul(Expr::new_constant(2.0), t);

        let t_var = CString::new("t").unwrap();

        let y_var = CString::new("y").unwrap();

        let yd_var = CString::new("y_dot").unwrap();

        let mut result = 0.0;

        let status = handle::rssn_num_cov_evaluate_action(
            &lagrangian,
            &path,
            t_var.as_ptr(),
            y_var.as_ptr(),
            yd_var.as_ptr(),
            0.0,
            1.0,
            &mut result,
        );

        assert_eq!(status, 0);

        assert_approx_eq!(result, 2.0, 1e-5);
    }
}

#[test]

fn test_cov_json_ffi() {

    unsafe {

        let t = Expr::new_variable("t");

        let y_dot = Expr::new_variable("y_dot");

        let lagrangian = Expr::new_mul(
            Expr::new_constant(0.5),
            Expr::new_pow(y_dot, Expr::new_constant(2.0)),
        );

        let path = t;

        let json_input = format!(
            r#"{{"lagrangian": {}, "path": {}, "t_var": "t", "path_var": "y", "path_dot_var": "y_dot", "t_range": [0.0, 1.0]}}"#,
            serde_json::to_string(&lagrangian).unwrap(),
            serde_json::to_string(&path).unwrap()
        );

        let c_json = CString::new(json_input).unwrap();

        let res_ptr = json::rssn_num_cov_evaluate_action_json(c_json.as_ptr());

        assert!(!res_ptr.is_null());

        let res_str = CStr::from_ptr(res_ptr).to_str().unwrap();

        let v: serde_json::Value = serde_json::from_str(res_str).unwrap();

        assert_approx_eq!(v["ok"].as_f64().unwrap(), 0.5, 1e-5);

        rssn_free_string(res_ptr);
    }
}

#[test]

fn test_cov_bincode_ffi() {

    unsafe {

        use rssn::ffi_apis::common::{from_bincode_buffer, to_bincode_buffer};
        use serde::{Deserialize, Serialize};

        #[derive(Serialize)]

        struct ActionInput {
            lagrangian: Expr,
            path: Expr,
            t_var: String,
            path_var: String,
            path_dot_var: String,
            t_range: (f64, f64),
        }

        let input = ActionInput {
            lagrangian: Expr::new_mul(
                Expr::new_constant(0.5),
                Expr::new_pow(Expr::new_variable("y_dot"), Expr::new_constant(2.0)),
            ),
            path: Expr::new_variable("t"),
            t_var: "t".to_string(),
            path_var: "y".to_string(),
            path_dot_var: "y_dot".to_string(),
            t_range: (0.0, 1.0),
        };

        let buffer = to_bincode_buffer(&input);

        let res_buffer = bincode_api::rssn_num_cov_evaluate_action_bincode(buffer);

        assert!(!res_buffer.is_null());

        #[derive(Deserialize)]

        struct FfiResult<T, E> {
            ok: Option<T>,
            #[allow(dead_code)]
            err: Option<E>,
        }

        let res: FfiResult<f64, String> = from_bincode_buffer(&res_buffer).unwrap();

        assert_approx_eq!(res.ok.unwrap(), 0.5, 1e-5);

        rssn_free_bincode_buffer(res_buffer);

        rssn_free_bincode_buffer(buffer);
    }
}
