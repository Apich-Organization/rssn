use std::ffi::CStr;
use std::ffi::CString;

use rssn::ffi_apis::common::from_bincode_buffer;
use rssn::ffi_apis::common::rssn_free_bincode_buffer;
use rssn::ffi_apis::common::rssn_free_string;
use rssn::ffi_apis::common::to_bincode_buffer;
use rssn::ffi_apis::ffi_api::FfiResult;
use rssn::ffi_apis::numerical_special_ffi::bincode_api;
use rssn::ffi_apis::numerical_special_ffi::handle;
use rssn::ffi_apis::numerical_special_ffi::json;
use serde::Serialize;

#[test]

fn test_special_handle_ffi() {

    // Gamma
    let g =
        handle::rssn_num_special_gamma(
            5.0,
        );

    assert!((g - 24.0).abs() < 1e-10);

    // Ln Gamma
    let lg = handle::rssn_num_special_ln_gamma(5.0);

    assert!(
        (lg - 24_f64.ln()).abs()
            < 1e-10
    );

    // Digamma
    let dg = handle::rssn_num_special_digamma(1.0);

    assert!(
        (dg - (-0.5772156649015329))
            .abs()
            < 1e-10
    );

    // Beta
    let b =
        handle::rssn_num_special_beta(
            2.0, 3.0,
        );

    assert!(b > 0.0);

    // Erf
    let e =
        handle::rssn_num_special_erf(
            0.0,
        );

    assert!(e.abs() < 1e-10);

    let e2 =
        handle::rssn_num_special_erf(
            10.0,
        );

    assert!((e2 - 1.0).abs() < 1e-10);

    // Erfc
    let ec =
        handle::rssn_num_special_erfc(
            0.0,
        );

    assert!((ec - 1.0).abs() < 1e-10);

    // Bessel J0
    let j0 = handle::rssn_num_special_bessel_j0(0.0);

    assert!((j0 - 1.0).abs() < 1e-10);

    // Bessel J1
    let j1 = handle::rssn_num_special_bessel_j1(0.0);

    assert!(j1.abs() < 1e-10);

    // Legendre
    let p2 = handle::rssn_num_special_legendre_p(2, 0.5);

    let expected =
        (3.0 * 0.5 * 0.5 - 1.0) / 2.0;

    assert!(
        (p2 - expected).abs() < 1e-10
    );

    // Chebyshev T
    let t1 = handle::rssn_num_special_chebyshev_t(1, 0.5);

    assert!((t1 - 0.5).abs() < 1e-10);

    // Hermite
    let h1 = handle::rssn_num_special_hermite_h(1, 1.0);

    assert!((h1 - 2.0).abs() < 1e-10);

    // Factorial
    let f5 = handle::rssn_num_special_factorial(5);

    assert!((f5 - 120.0).abs() < 1e-10);

    // Double Factorial
    let df5 = handle::rssn_num_special_double_factorial(5);

    assert!((df5 - 15.0).abs() < 1e-10);

    // Binomial
    let c52 = handle::rssn_num_special_binomial(5, 2);

    assert!((c52 - 10.0).abs() < 1e-10);

    // Zeta
    let z2 =
        handle::rssn_num_special_zeta(
            2.0,
        );

    assert!(
        (z2 - std::f64::consts::PI
            .powi(2)
            / 6.0)
            .abs()
            < 1e-3
    );

    // Sinc
    let s0 =
        handle::rssn_num_special_sinc(
            0.0,
        );

    assert!((s0 - 1.0).abs() < 1e-10);

    // Sigmoid
    let sig0 = handle::rssn_num_special_sigmoid(0.0);

    assert!((sig0 - 0.5).abs() < 1e-10);

    // Softplus
    let sp0 = handle::rssn_num_special_softplus(0.0);

    assert!(
        (sp0 - 2_f64.ln()).abs()
            < 1e-10
    );

    // Logit
    let logit_half =
        handle::rssn_num_special_logit(
            0.5,
        );

    assert!(logit_half.abs() < 1e-10);
}

#[derive(Serialize)]

struct SingleInput {
    x : f64,
}

#[derive(Serialize)]

struct TwoInput {
    a : f64,
    b : f64,
}

#[derive(Serialize)]

struct PolyInput {
    n : u32,
    x : f64,
}

#[derive(Serialize)]

struct IntInput {
    n : u64,
}

#[derive(Serialize)]

struct BinomialInput {
    n : u64,
    k : u64,
}

#[test]

fn test_special_json_ffi() {

    unsafe {

        // Gamma
        let input = SingleInput {
            x : 5.0,
        };

        let json_str =
            serde_json::to_string(
                &input,
            )
            .unwrap();

        let c_json =
            CString::new(json_str)
                .unwrap();

        let res_ptr = json::rssn_num_special_gamma_json(c_json.as_ptr());

        let res_str =
            CStr::from_ptr(res_ptr)
                .to_str()
                .unwrap();

        let v : serde_json::Value =
            serde_json::from_str(
                res_str,
            )
            .unwrap();

        assert!(
            (v["ok"]
                .as_f64()
                .unwrap()
                - 24.0)
                .abs()
                < 1e-10
        );

        rssn_free_string(res_ptr);

        // Erf
        let input = SingleInput {
            x : 0.0,
        };

        let json_str =
            serde_json::to_string(
                &input,
            )
            .unwrap();

        let c_json =
            CString::new(json_str)
                .unwrap();

        let res_ptr = json::rssn_num_special_erf_json(c_json.as_ptr());

        let res_str =
            CStr::from_ptr(res_ptr)
                .to_str()
                .unwrap();

        let v : serde_json::Value =
            serde_json::from_str(
                res_str,
            )
            .unwrap();

        assert!(
            v["ok"]
                .as_f64()
                .unwrap()
                .abs()
                < 1e-10
        );

        rssn_free_string(res_ptr);

        // Factorial
        let input = IntInput {
            n : 5,
        };

        let json_str =
            serde_json::to_string(
                &input,
            )
            .unwrap();

        let c_json =
            CString::new(json_str)
                .unwrap();

        let res_ptr = json::rssn_num_special_factorial_json(c_json.as_ptr());

        let res_str =
            CStr::from_ptr(res_ptr)
                .to_str()
                .unwrap();

        let v : serde_json::Value =
            serde_json::from_str(
                res_str,
            )
            .unwrap();

        assert!(
            (v["ok"]
                .as_f64()
                .unwrap()
                - 120.0)
                .abs()
                < 1e-10
        );

        rssn_free_string(res_ptr);

        // Binomial
        let input = BinomialInput {
            n : 5,
            k : 2,
        };

        let json_str =
            serde_json::to_string(
                &input,
            )
            .unwrap();

        let c_json =
            CString::new(json_str)
                .unwrap();

        let res_ptr = json::rssn_num_special_binomial_json(c_json.as_ptr());

        let res_str =
            CStr::from_ptr(res_ptr)
                .to_str()
                .unwrap();

        let v : serde_json::Value =
            serde_json::from_str(
                res_str,
            )
            .unwrap();

        assert!(
            (v["ok"]
                .as_f64()
                .unwrap()
                - 10.0)
                .abs()
                < 1e-10
        );

        rssn_free_string(res_ptr);

        // Legendre
        let input = PolyInput {
            n : 2,
            x : 0.5,
        };

        let json_str =
            serde_json::to_string(
                &input,
            )
            .unwrap();

        let c_json =
            CString::new(json_str)
                .unwrap();

        let res_ptr = json::rssn_num_special_legendre_p_json(c_json.as_ptr());

        let res_str =
            CStr::from_ptr(res_ptr)
                .to_str()
                .unwrap();

        let v : serde_json::Value =
            serde_json::from_str(
                res_str,
            )
            .unwrap();

        let expected =
            (3.0 * 0.5 * 0.5 - 1.0)
                / 2.0;

        assert!(
            (v["ok"]
                .as_f64()
                .unwrap()
                - expected)
                .abs()
                < 1e-10
        );

        rssn_free_string(res_ptr);
    }
}

#[test]

fn test_special_bincode_ffi() {

    unsafe {

        #[derive(Serialize)]

        struct SingleInputB {
            x : f64,
        }

        // Gamma
        let input = SingleInputB {
            x : 5.0,
        };

        let buffer =
            to_bincode_buffer(&input);

        let res_buffer = bincode_api::rssn_num_special_gamma_bincode(buffer);

        let res : FfiResult<
            f64,
            String,
        > = from_bincode_buffer(
            &res_buffer,
        )
        .unwrap();

        assert!(
            (res.ok.unwrap() - 24.0)
                .abs()
                < 1e-10
        );

        rssn_free_bincode_buffer(
            res_buffer,
        );

        rssn_free_bincode_buffer(
            buffer,
        );

        // Erf
        let input = SingleInputB {
            x : 0.0,
        };

        let buffer =
            to_bincode_buffer(&input);

        let res_buffer = bincode_api::rssn_num_special_erf_bincode(buffer);

        let res : FfiResult<
            f64,
            String,
        > = from_bincode_buffer(
            &res_buffer,
        )
        .unwrap();

        assert!(
            res.ok
                .unwrap()
                .abs()
                < 1e-10
        );

        rssn_free_bincode_buffer(
            res_buffer,
        );

        rssn_free_bincode_buffer(
            buffer,
        );

        // Factorial
        #[derive(Serialize)]

        struct IntInputB {
            n : u64,
        }

        let input = IntInputB {
            n : 5,
        };

        let buffer =
            to_bincode_buffer(&input);

        let res_buffer = bincode_api::rssn_num_special_factorial_bincode(buffer);

        let res : FfiResult<
            f64,
            String,
        > = from_bincode_buffer(
            &res_buffer,
        )
        .unwrap();

        assert!(
            (res.ok.unwrap() - 120.0)
                .abs()
                < 1e-10
        );

        rssn_free_bincode_buffer(
            res_buffer,
        );

        rssn_free_bincode_buffer(
            buffer,
        );
    }
}
