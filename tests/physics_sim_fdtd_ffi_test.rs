//! FFI tests for the physics sim FDTD electrodynamics module.

use std::ffi::CStr;
use std::ffi::CString;

#[test]

fn test_fdtd_handle_ffi() {

    unsafe {

        let matrix_ptr = rssn::ffi_apis::physics_sim_fdtd_ffi::handle::rssn_physics_sim_fdtd_run_2d(
            40, 40, 50, 20, 20, 0.1,
        );

        assert!(!matrix_ptr.is_null());

        let matrix = &*matrix_ptr;

        assert_eq!(matrix.rows(), 40);

        assert_eq!(matrix.cols(), 40);

        rssn::ffi_apis::numerical_matrix_ffi::handle::rssn_num_matrix_free(matrix_ptr);
    }
}

#[test]

fn test_fdtd_json_ffi() {

    let input = r#"{
        "width": 30,
        "height": 30,
        "time_steps": 20,
        "source_pos": [15, 15],
        "source_freq": 0.1
    }"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let res_ptr = rssn::ffi_apis::physics_sim_fdtd_ffi::json::rssn_physics_sim_fdtd_run_json(
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
