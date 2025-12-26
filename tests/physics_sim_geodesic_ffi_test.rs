//! FFI tests for the physics sim geodesic relativity module.

use std::ffi::{CStr, CString};

#[test]
fn test_geodesic_handle_ffi() {
    unsafe {
        let matrix_ptr =
            rssn::ffi_apis::physics_sim_geodesic_ffi::handle::rssn_physics_sim_geodesic_run(
                1.0, 10.0, 0.0, 0.0, 0.035, 100.0, 0.1,
            );
        assert!(!matrix_ptr.is_null());
        let matrix = &*matrix_ptr;
        assert!(matrix.rows() > 0);
        assert_eq!(matrix.cols(), 2);
        rssn::ffi_apis::numerical_matrix_ffi::handle::rssn_num_matrix_free(matrix_ptr);
    }
}

#[test]
fn test_geodesic_json_ffi() {
    let input = r#"{
        "black_hole_mass": 1.0,
        "initial_state": [10.0, 0.0, 0.0, 0.035],
        "proper_time_end": 50.0,
        "initial_dt": 0.1
    }"#;
    let c_input = CString::new(input).unwrap();
    unsafe {
        let res_ptr =
            rssn::ffi_apis::physics_sim_geodesic_ffi::json::rssn_physics_sim_geodesic_run_json(
                c_input.as_ptr(),
            );
        assert!(!res_ptr.is_null());
        let res_str = CStr::from_ptr(res_ptr).to_string_lossy();
        assert!(res_str.contains("\"ok\":"));
        rssn::ffi_apis::ffi_api::free_string(res_ptr);
    }
}
