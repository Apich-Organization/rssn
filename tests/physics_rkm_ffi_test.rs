//! FFI tests for physics RKM module.

use rssn::ffi_apis::ffi_api::FfiResult;
use std::ffi::CString;

#[test]

fn test_lorenz_handle_ffi() {

    unsafe {

        let matrix_ptr =
            rssn::ffi_apis::physics_rkm_ffi::handle::rssn_physics_rkm_simulate_lorenz();

        assert!(!matrix_ptr.is_null());

        let matrix = &*matrix_ptr;

        assert!(matrix.rows() > 0);

        assert_eq!(matrix.cols(), 4); // time, x, y, z
                                      // Clean up
        let _ =
            Box::from_raw(matrix_ptr);
    }
}

#[test]

fn test_lorenz_json_ffi() {

    let input = r#"{
        "sigma": 10.0,
        "rho": 28.0,
        "beta": 2.6666,
        "y0": [1.0, 1.0, 1.0],
        "t_span": [0.0, 1.0],
        "dt_initial": 0.01,
        "tol": [1e-6, 1e-6]
    }"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let res_ptr =
            rssn::ffi_apis::physics_rkm_ffi::json::rssn_physics_rkm_lorenz_json(c_input.as_ptr());

        assert!(!res_ptr.is_null());

        let res_str =
            std::ffi::CStr::from_ptr(
                res_ptr,
            )
            .to_string_lossy();

        let _res: serde_json::Value =
            serde_json::from_str(
                &res_str,
            )
            .unwrap();

        // Check if OK
        assert!(
            res_str.contains("\"ok\":")
        );
    }
}

#[test]

fn test_damped_oscillator_json_ffi() {

    let input = r#"{
        "omega": 1.0,
        "zeta": 0.15,
        "y0": [1.0, 0.0],
        "t_span": [0.0, 1.0],
        "dt": 0.1
    }"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let res_ptr =
            rssn::ffi_apis::physics_rkm_ffi::json::rssn_physics_rkm_damped_oscillator_json(
                c_input.as_ptr(),
            );

        assert!(!res_ptr.is_null());

        let res_str =
            std::ffi::CStr::from_ptr(
                res_ptr,
            )
            .to_string_lossy();

        assert!(
            res_str.contains("\"ok\":")
        );
    }
}
