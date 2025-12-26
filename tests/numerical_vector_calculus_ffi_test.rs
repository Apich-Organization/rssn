use std::ffi::CStr;
use std::ffi::CString;

use assert_approx_eq::assert_approx_eq;
use rssn::ffi_apis::common::rssn_free_bincode_buffer;
use rssn::ffi_apis::common::rssn_free_string;
use rssn::ffi_apis::numerical_vector_calculus_ffi::bincode_api;
use rssn::ffi_apis::numerical_vector_calculus_ffi::handle;
use rssn::ffi_apis::numerical_vector_calculus_ffi::json;
use rssn::symbolic::core::Expr;

#[test]

fn test_numerical_vector_calculus_handle_ffi() {

    unsafe {

        let x = Expr::new_variable("x");

        let y = Expr::new_variable("y");

        let funcs = vec![
            &Expr::new_pow(
                x.clone(),
                Expr::new_constant(2.0),
            ) as *const Expr,
            &Expr::new_pow(
                y.clone(),
                Expr::new_constant(2.0),
            ) as *const Expr,
        ];

        let var_x = CString::new("x").unwrap();

        let var_y = CString::new("y").unwrap();

        let vars = vec![
            var_x.as_ptr(),
            var_y.as_ptr(),
        ];

        let point = vec![1.0, 2.0];

        let mut result = 0.0;

        // Divergence
        let status = handle::rssn_num_vector_calculus_divergence(
            funcs.as_ptr(),
            2,
            vars.as_ptr(),
            point.as_ptr(),
            2,
            &mut result,
        );

        assert_eq!(status, 0);

        assert_approx_eq!(result, 6.0, 1e-5);

        // Laplacian
        let f = Expr::new_pow(
            x.clone(),
            Expr::new_constant(2.0),
        );

        let mut lap_res = 0.0;

        let status2 = handle::rssn_num_vector_calculus_laplacian(
            &f,
            vars.as_ptr(),
            point.as_ptr(),
            2,
            &mut lap_res,
        );

        assert_eq!(status2, 0);

        assert_approx_eq!(lap_res, 2.0, 1e-5);
    }
}

#[test]

fn test_numerical_vector_calculus_json_ffi() {

    unsafe {

        let x = Expr::new_variable("x");

        let y = Expr::new_variable("y");

        let f1 = Expr::new_pow(
            x.clone(),
            Expr::new_constant(2.0),
        );

        let f2 = Expr::new_pow(
            y.clone(),
            Expr::new_constant(2.0),
        );

        let json_input = format!(
            r#"{{"funcs": [{}, {}], "vars": ["x", "y"], "point": [1.0, 2.0]}}"#,
            serde_json::to_string(&f1).unwrap(),
            serde_json::to_string(&f2).unwrap()
        );

        let c_json = CString::new(json_input).unwrap();

        let res_ptr = json::rssn_num_vector_calculus_divergence_json(c_json.as_ptr());

        assert!(!res_ptr.is_null());

        let res_str = CStr::from_ptr(res_ptr)
            .to_str()
            .unwrap();

        let v : serde_json::Value = serde_json::from_str(res_str).unwrap();

        assert_approx_eq!(
            v["ok"]
                .as_f64()
                .unwrap(),
            6.0,
            1e-5
        );

        rssn_free_string(res_ptr);
    }
}

#[test]

fn test_numerical_vector_calculus_bincode_ffi() {

    unsafe {

        use rssn::ffi_apis::common::from_bincode_buffer;
        use rssn::ffi_apis::common::to_bincode_buffer;
        use serde::Deserialize;
        use serde::Serialize;

        #[derive(Serialize)]

        struct LaplacianInput {
            f : Expr,
            vars : Vec<String>,
            point : Vec<f64>,
        }

        let x = Expr::new_variable("x");

        let input = LaplacianInput {
            f : Expr::new_pow(
                x,
                Expr::new_constant(2.0),
            ),
            vars : vec!["x".to_string()],
            point : vec![1.0],
        };

        let buffer = to_bincode_buffer(&input);

        let res_buffer = bincode_api::rssn_num_vector_calculus_laplacian_bincode(buffer);

        assert!(!res_buffer.is_null());

        #[derive(Deserialize)]

        struct FfiResult<T, E> {
            ok : Option<T>,
            #[allow(dead_code)]
            err : Option<E>,
        }

        let res : FfiResult<f64, String> = from_bincode_buffer(&res_buffer).unwrap();

        assert_approx_eq!(
            res.ok.unwrap(),
            2.0,
            1e-5
        );

        rssn_free_bincode_buffer(res_buffer);

        rssn_free_bincode_buffer(buffer);
    }
}
