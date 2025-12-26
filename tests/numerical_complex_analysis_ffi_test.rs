use std::ffi::CStr;
use std::ffi::CString;

use num_complex::Complex;
use rssn::ffi_apis::common::rssn_free_bincode_buffer;
use rssn::ffi_apis::common::rssn_free_string;
use rssn::ffi_apis::numerical_complex_analysis_ffi::bincode_api;
use rssn::ffi_apis::numerical_complex_analysis_ffi::handle;
use rssn::ffi_apis::numerical_complex_analysis_ffi::json;
use rssn::symbolic::core::Expr;

#[test]

fn test_complex_handle_ffi() {

    unsafe {

        let z = Expr::Variable(
            "z".to_string(),
        );

        let expr = Expr::new_pow(
            z,
            Expr::Constant(2.0),
        );

        let z_name =
            CString::new("z").unwrap();

        let var_names =
            [z_name.as_ptr()];

        let var_re = [2.0];

        let var_im = [0.0];

        let mut res_re = 0.0;

        let mut res_im = 0.0;

        let status = handle::rssn_num_complex_eval(
            &expr,
            var_names.as_ptr(),
            var_re.as_ptr(),
            var_im.as_ptr(),
            1,
            &mut res_re,
            &mut res_im,
        );

        if status != 0 {

            let err =
                CStr::from_ptr(rssn::ffi_apis::ffi_api::rssn_get_last_error()).to_string_lossy();

            panic!(
                "FFI call failed with \
                 status {}: {}",
                status, err
            );
        }

        assert_eq!(res_re, 4.0);

        assert_eq!(res_im, 0.0);
    }
}

#[test]

fn test_complex_json_ffi() {

    unsafe {

        let json_input = r#"{
            "expr": {"Variable": "z"},
            "vars": {"z": [0.0, 1.0]}
        }"#;

        let c_json =
            CString::new(json_input)
                .unwrap();

        let res_ptr = json::rssn_num_complex_eval_json(c_json.as_ptr());

        assert!(!res_ptr.is_null());

        let res_str =
            CStr::from_ptr(res_ptr)
                .to_str()
                .unwrap();

        let v: serde_json::Value =
            serde_json::from_str(
                res_str,
            )
            .unwrap();

        if v["ok"].is_null() {

            panic!(
                "FFI JSON call \
                 failed: {}",
                v["err"]
            );
        }

        let res = v["ok"]
            .as_array()
            .unwrap();

        assert_eq!(
            res[0]
                .as_f64()
                .unwrap(),
            0.0
        );

        assert_eq!(
            res[1]
                .as_f64()
                .unwrap(),
            1.0
        );

        rssn_free_string(res_ptr);
    }
}

#[test]

fn test_complex_bincode_ffi() {

    unsafe {

        use std::collections::HashMap;

        use rssn::ffi_apis::common::from_bincode_buffer;
        use rssn::ffi_apis::common::to_bincode_buffer;
        use serde::Deserialize;
        use serde::Serialize;

        #[derive(Serialize)]

        struct EvalInput {
            expr: Expr,
            vars: HashMap<
                String,
                Complex<f64>,
            >,
        }

        let mut vars = HashMap::new();

        vars.insert(
            "z".to_string(),
            Complex::new(0.0, 1.0),
        );

        let input = EvalInput {
            expr: Expr::Variable(
                "z".to_string(),
            ),
            vars,
        };

        let buffer =
            to_bincode_buffer(&input);

        let res_buffer = bincode_api::rssn_num_complex_eval_bincode(buffer);

        assert!(!res_buffer.is_null());

        #[derive(Deserialize)]

        struct FfiResult<T, E> {
            ok: Option<T>,
            #[allow(dead_code)]
            err: Option<E>,
        }

        let res: FfiResult<
            Complex<f64>,
            String,
        > = from_bincode_buffer(
            &res_buffer,
        )
        .unwrap();

        assert_eq!(
            res.ok.unwrap(),
            Complex::new(0.0, 1.0)
        );

        rssn_free_bincode_buffer(
            res_buffer,
        );

        rssn_free_bincode_buffer(
            buffer,
        );
    }
}
