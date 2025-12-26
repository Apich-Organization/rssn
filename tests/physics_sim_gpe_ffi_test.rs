//! FFI tests for the physics sim GPE superfluidity module.

use std::ffi::{CStr, CString};

#[test]
fn test_gpe_handle_ffi() {
    unsafe {
        let matrix_ptr = rssn::ffi_apis::physics_sim_gpe_ffi::handle::rssn_physics_sim_gpe_run_ground_state_finder(
            16, 16, 10.0, 10.0, 0.1, 10, 50.0, 1.0
        );
        assert!(!matrix_ptr.is_null());
        let matrix = &*matrix_ptr;
        assert_eq!(matrix.rows(), 16);
        assert_eq!(matrix.cols(), 16);
        rssn::ffi_apis::numerical_matrix_ffi::handle::rssn_num_matrix_free(matrix_ptr);
    }
}

#[test]
fn test_gpe_json_ffi() {
    let input = r#"{
        "nx": 16,
        "ny": 16,
        "lx": 10.0,
        "ly": 10.0,
        "d_tau": 0.1,
        "time_steps": 5,
        "g": 50.0,
        "trap_strength": 1.0
    }"#;
    let c_input = CString::new(input).unwrap();
    unsafe {
        let res_ptr = rssn::ffi_apis::physics_sim_gpe_ffi::json::rssn_physics_sim_gpe_run_json(c_input.as_ptr());
        assert!(!res_ptr.is_null());
        let res_str = CStr::from_ptr(res_ptr).to_string_lossy();
        assert!(res_str.contains("\"ok\":"));
        rssn::ffi_apis::ffi_api::free_string(res_ptr);
    }
}
