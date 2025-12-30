//! FFI tests for the physics SM module.

use std::ffi::CStr;
use std::ffi::CString;

#[test]

fn test_sm_simulate_1d_handle_ffi() {

    unsafe {

        let matrix_ptr =
            rssn::ffi_apis::physics_sm_ffi::handle::rssn_physics_sm_simulate_1d_advection();

        assert!(!matrix_ptr.is_null());

        let matrix = &*matrix_ptr;

        assert_eq!(matrix.rows(), 1);

        assert_eq!(matrix.cols(), 128);

        rssn::ffi_apis::numerical_matrix_ffi::handle::rssn_num_matrix_free(matrix_ptr as *mut rssn::ffi_apis::numerical_matrix_ffi::handle::RssnMatrixHandle);
    }
}

#[test]

fn test_sm_solve_2d_json_ffi() {

    let input = r#"{
        "initial_condition": [1.0, 0.0, 0.0, 0.0],
        "config": {
            "width": 2,
            "height": 2,
            "dx": 1.0,
            "dy": 1.0,
            "c": [0.1, 0.1],
            "d": 0.01,
            "dt": 0.1,
            "steps": 1
        }
    }"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let res_ptr = rssn::ffi_apis::physics_sm_ffi::json::rssn_physics_sm_solve_advection_2d_json(
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
