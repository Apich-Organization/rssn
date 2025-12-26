use std::ffi::CStr;
use std::ffi::CString;

use assert_approx_eq::assert_approx_eq;
use rssn::ffi_apis::common::rssn_free_bincode_buffer;
use rssn::ffi_apis::common::rssn_free_string;
use rssn::ffi_apis::numerical_integrate_ffi::*;
use rssn::numerical::integrate::QuadratureMethod;
use rssn::symbolic::core::Expr;

#[test]

fn test_numerical_quadrature_handle_ffi(
) {

    unsafe {

        // Integrate x^2 from 0 to 1
        let x = Expr::new_variable("x");

        let two =
            Expr::new_constant(2.0);

        let expr =
            Expr::new_pow(x, two);

        let var =
            CString::new("x").unwrap();

        let mut result: f64 = 0.0;

        // Method 1: Simpson (method=1)
        let status = handle::rssn_numerical_quadrature(
            &expr,
            var.as_ptr(),
            0.0,
            1.0,
            100,
            1, // Simpson
            &mut result,
        );

        assert_eq!(status, 0);

        assert_approx_eq!(
            result,
            1.0 / 3.0,
            1e-10f64
        );

        // Method 2: Gauss-Legendre (method=4)
        let status_gl = handle::rssn_numerical_quadrature(
            &expr,
            var.as_ptr(),
            0.0,
            1.0,
            0, // n_steps not used for GL n=5
            4, // GaussLegendre
            &mut result,
        );

        assert_eq!(status_gl, 0);

        assert_approx_eq!(
            result,
            1.0 / 3.0,
            1e-10f64
        );
    }
}

#[test]

fn test_numerical_quadrature_json_ffi()
{

    unsafe {

        // Integrate x^2 from 0 to 1 using JSON
        let x = Expr::new_variable("x");

        let two =
            Expr::new_constant(2.0);

        let expr =
            Expr::new_pow(x, two);

        let expr_json =
            serde_json::to_string(
                &expr,
            )
            .unwrap();

        let json_input = format!(
            r#"{{"expr": {}, "var": "x", "a": 0.0, "b": 1.0, "n_steps": 100, "method": "Simpson"}}"#,
            expr_json
        );

        let c_json =
            CString::new(json_input)
                .unwrap();

        let res_ptr = json::rssn_numerical_quadrature_json(c_json.as_ptr());

        assert!(!res_ptr.is_null());

        let res_str =
            CStr::from_ptr(res_ptr)
                .to_str()
                .unwrap();

        let v: serde_json::Value =
            serde_json::from_str(
                res_str,
            )
            .expect(
                "Failed to parse \
                 result JSON",
            );

        if let Some(err) = v.get("err")
        {

            if !err.is_null() {

                panic!(
                    "FFI error: {}",
                    err
                );
            }
        }

        let ok = v["ok"]
            .as_f64()
            .expect(
                "Result 'ok' should \
                 be f64",
            );

        assert_approx_eq!(
            ok,
            1.0 / 3.0,
            1e-10f64
        );

        rssn_free_string(res_ptr);
    }
}

#[test]

fn test_numerical_quadrature_bincode_ffi(
) {

    unsafe {

        use rssn::ffi_apis::common::from_bincode_buffer;
        use rssn::ffi_apis::common::to_bincode_buffer;
        use serde::Deserialize;
        use serde::Serialize;

        #[derive(Serialize)]

        struct QuadratureInput {
            expr: Expr,
            var: String,
            a: f64,
            b: f64,
            n_steps: usize,
            method: QuadratureMethod,
        }

        let x = Expr::new_variable("x");

        let two =
            Expr::new_constant(2.0);

        let expr =
            Expr::new_pow(x, two);

        let input = QuadratureInput {
            expr,
            var : "x".to_string(),
            a : 0.0,
            b : 2.0,
            n_steps : 100,
            method : QuadratureMethod::Romberg,
        };

        let buffer =
            to_bincode_buffer(&input);

        let res_buffer = bincode_api::rssn_numerical_quadrature_bincode(buffer);

        assert!(!res_buffer.is_null());

        #[derive(Deserialize)]

        struct FfiResult<T, E> {
            #[allow(dead_code)]
            ok: Option<T>,
            #[allow(dead_code)]
            err: Option<E>,
        }

        let res: FfiResult<
            f64,
            String,
        > = from_bincode_buffer(
            &res_buffer,
        )
        .expect(
            "Failed to decode bincode \
             result",
        );

        assert!(res.err.is_none());

        // Integral of x^2 from 0 to 2 is 8/3
        assert_approx_eq!(
            res.ok.unwrap(),
            8.0 / 3.0,
            1e-8f64
        );

        rssn_free_bincode_buffer(
            res_buffer,
        );

        rssn_free_bincode_buffer(
            buffer,
        );
    }
}
