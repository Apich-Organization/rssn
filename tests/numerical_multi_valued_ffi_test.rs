use rssn::ffi_apis::common::{
    from_bincode_buffer,
    rssn_free_bincode_buffer,
    rssn_free_string,
    to_bincode_buffer,
    BincodeBuffer,
};
use rssn::ffi_apis::ffi_api::FfiResult;
use rssn::ffi_apis::numerical_multi_valued_ffi::{
    bincode_api,
    handle,
    json,
};
use rssn::symbolic::core::Expr;
use serde::{
    Deserialize,
    Serialize,
};
use std::ffi::{
    CStr,
    CString,
};

#[test]

fn test_mv_handle_ffi() {

    unsafe {

        let z = Expr::Variable("z".to_string());

        let f = Expr::new_sub(
            Expr::new_pow(
                z.clone(),
                Expr::Constant(2.0),
            ),
            Expr::Constant(1.0),
        );

        let f_prime = Expr::new_mul(
            Expr::Constant(2.0),
            z.clone(),
        );

        let mut res_re = 0.0;

        let mut res_im = 0.0;

        // Root near 2.0 should be 1.0
        let status = handle::rssn_num_mv_newton_method_complex(
            &f,
            &f_prime,
            2.0,
            0.0,
            1e-6,
            100,
            &mut res_re,
            &mut res_im,
        );

        assert_eq!(status, 0);

        assert!((res_re - 1.0).abs() < 1e-5);

        assert!(res_im.abs() < 1e-5);
    }
}

#[test]

fn test_mv_json_ffi() {

    unsafe {

        let z = Expr::Variable("z".to_string());

        let f = Expr::new_sub(
            Expr::new_pow(
                z.clone(),
                Expr::Constant(2.0),
            ),
            Expr::Constant(1.0),
        );

        let f_prime = Expr::new_mul(
            Expr::Constant(2.0),
            z.clone(),
        );

        #[derive(Serialize)]

        struct NewtonInput {
            f: Expr,
            f_prime: Expr,
            start_re: f64,
            start_im: f64,
            tolerance: f64,
            max_iter: usize,
        }

        let input = NewtonInput {
            f,
            f_prime,
            start_re: 2.0,
            start_im: 0.0,
            tolerance: 1e-6,
            max_iter: 100,
        };

        let json_str = serde_json::to_string(&input).unwrap();

        let c_json = CString::new(json_str).unwrap();

        let res_ptr = json::rssn_num_mv_newton_method_complex_json(c_json.as_ptr());

        assert!(!res_ptr.is_null());

        let res_str = CStr::from_ptr(res_ptr)
            .to_str()
            .unwrap();

        let v: serde_json::Value = serde_json::from_str(res_str).unwrap();

        // Check result
        let res_obj = v["ok"]
            .as_object()
            .unwrap();

        let re = res_obj["re"]
            .as_f64()
            .unwrap();

        assert!((re - 1.0).abs() < 1e-5);

        rssn_free_string(res_ptr);
    }
}

#[test]

fn test_mv_bincode_ffi() {

    unsafe {

        let z = Expr::Variable("z".to_string());

        let f = Expr::new_sub(
            Expr::new_pow(
                z.clone(),
                Expr::Constant(2.0),
            ),
            Expr::Constant(1.0),
        );

        let f_prime = Expr::new_mul(
            Expr::Constant(2.0),
            z.clone(),
        );

        #[derive(Serialize)]

        struct NewtonInput {
            f: Expr,
            f_prime: Expr,
            start_re: f64,
            start_im: f64,
            tolerance: f64,
            max_iter: usize,
        }

        #[derive(Deserialize)]

        struct ComplexResult {
            re: f64,
            im: f64,
        }

        let input = NewtonInput {
            f,
            f_prime,
            start_re: 2.0,
            start_im: 0.0,
            tolerance: 1e-6,
            max_iter: 100,
        };

        let buffer = to_bincode_buffer(&input);

        let res_buffer = bincode_api::rssn_num_mv_newton_method_complex_bincode(buffer);

        assert!(!res_buffer.is_null());

        let res: FfiResult<ComplexResult, String> = from_bincode_buffer(&res_buffer).unwrap();

        let root = res.ok.unwrap();

        assert!((root.re - 1.0).abs() < 1e-5);

        rssn_free_bincode_buffer(res_buffer);

        rssn_free_bincode_buffer(buffer);
    }
}

#[test]

fn test_mv_handle_others() {

    unsafe {

        let mut res_re = 0.0;

        let mut res_im = 0.0;

        // log(1, 0) = 0
        handle::rssn_num_mv_complex_log_k(
            1.0,
            0.0,
            0,
            &mut res_re,
            &mut res_im,
        );

        assert!(res_re.abs() < 1e-9);

        assert!(res_im.abs() < 1e-9);

        // sqrt(1, 1) = -1
        handle::rssn_num_mv_complex_sqrt_k(
            1.0,
            0.0,
            1,
            &mut res_re,
            &mut res_im,
        );

        assert!((res_re + 1.0).abs() < 1e-9);

        assert!(res_im.abs() < 1e-9);
    }
}

#[test]

fn test_mv_json_others() {

    unsafe {

        #[derive(Serialize)]

        struct LogSqrtInput {
            re: f64,
            im: f64,
            k: i32,
        }

        // log(1, 0) = 0
        let input = LogSqrtInput {
            re: 1.0,
            im: 0.0,
            k: 0,
        };

        let json_str = serde_json::to_string(&input).unwrap();

        let c_json = CString::new(json_str).unwrap();

        let res_ptr = json::rssn_num_mv_complex_log_k_json(c_json.as_ptr());

        assert!(!res_ptr.is_null());

        let res_str = CStr::from_ptr(res_ptr)
            .to_str()
            .unwrap();

        let v: serde_json::Value = serde_json::from_str(res_str).unwrap();

        let res_obj = v["ok"]
            .as_object()
            .unwrap();

        let re = res_obj["re"]
            .as_f64()
            .unwrap();

        let im = res_obj["im"]
            .as_f64()
            .unwrap();

        assert!(re.abs() < 1e-9);

        assert!(im.abs() < 1e-9);

        rssn_free_string(res_ptr);
    }
}

#[test]

fn test_mv_bincode_others() {

    unsafe {

        #[derive(Serialize)]

        struct LogSqrtInput {
            re: f64,
            im: f64,
            k: i32,
        }

        #[derive(Deserialize)]

        struct ComplexResult {
            re: f64,
            im: f64,
        }

        // sqrt(1, 0) = 1
        let input = LogSqrtInput {
            re: 1.0,
            im: 0.0,
            k: 0,
        };

        let buffer = to_bincode_buffer(&input);

        let res_buffer = bincode_api::rssn_num_mv_complex_sqrt_k_bincode(buffer);

        assert!(!res_buffer.is_null());

        let res: FfiResult<ComplexResult, String> = from_bincode_buffer(&res_buffer).unwrap();

        let root = res.ok.unwrap();

        assert!((root.re - 1.0).abs() < 1e-9);

        assert!(root.im.abs() < 1e-9);

        rssn_free_bincode_buffer(res_buffer);

        rssn_free_bincode_buffer(buffer);
    }
}
