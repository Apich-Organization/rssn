//! FFI tests for the physics CNM module.

use rssn::ffi_apis::ffi_api::FfiResult;
use std::ffi::{
    CStr,
    CString,
};

#[test]

fn test_cnm_heat_1d_handle_ffi() {

    let initial = [
        1.0, 1.0, 1.0, 1.0, 1.0,
    ];

    let mut out_size = 0;

    unsafe {

        let res_ptr = rssn::ffi_apis::physics_cnm_ffi::handle::rssn_physics_cnm_solve_heat_1d(
            initial.as_ptr(),
            5,
            0.1,
            0.01,
            0.1,
            5,
            &mut out_size,
        );

        assert!(!res_ptr.is_null());

        assert_eq!(out_size, 5);

        rssn::ffi_apis::physics_cnm_ffi::handle::rssn_free_f64_cnm_array(res_ptr, out_size);
    }
}

#[test]

fn test_cnm_heat_2d_json_ffi() {

    let input = r#"{
        "initial_condition": [0,0,0,0,0,0,1,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0],
        "config": {
            "nx": 5,
            "ny": 5,
            "dx": 0.2,
            "dy": 0.2,
            "dt": 0.01,
            "d_coeff": 0.1,
            "steps": 2
        }
    }"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let res_ptr = rssn::ffi_apis::physics_cnm_ffi::json::rssn_physics_cnm_solve_heat_2d_json(
            c_input.as_ptr(),
        );

        assert!(!res_ptr.is_null());

        let res_str =
            CStr::from_ptr(res_ptr)
                .to_string_lossy();

        assert!(
            res_str.contains("\"ok\":")
        );

        rssn::ffi_apis::ffi_api::free_string(res_ptr);
    }
}
