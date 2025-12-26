use std::ffi::{CStr, CString};

#[test]
fn test_em_simulate_oscillator_handle_ffi() {
    unsafe {
        let matrix_ptr =
            rssn::ffi_apis::physics_em_ffi::handle::rssn_physics_em_simulate_oscillator_forward();
        assert!(!matrix_ptr.is_null());
        let matrix = &*matrix_ptr;
        assert!(matrix.rows() > 0);
        assert_eq!(matrix.cols(), 3); // time, y0, y1
        rssn::ffi_apis::numerical_matrix_ffi::handle::rssn_num_matrix_free(matrix_ptr);
    }
}

#[test]
fn test_em_solve_json_ffi() {
    let input = r#"{
        "system_type": "oscillator",
        "params": {
            "omega": 1.0,
            "zeta": 0.1
        },
        "y0": [1.0, 0.0],
        "t_span": [0.0, 1.0],
        "dt": 0.1,
        "method": "heun"
    }"#;
    let c_input = CString::new(input).unwrap();
    unsafe {
        let res_ptr =
            rssn::ffi_apis::physics_em_ffi::json::rssn_physics_em_solve_json(c_input.as_ptr());
        assert!(!res_ptr.is_null());
        let res_str = CStr::from_ptr(res_ptr).to_string_lossy();
        assert!(res_str.contains("\"ok\":"));
        rssn::ffi_apis::ffi_api::free_string(res_ptr);
    }
}
