//! FFI tests for the physics sim linear elasticity module.

use std::ffi::CStr;
use std::ffi::CString;

#[test]

fn test_linear_elasticity_handle_ffi() {

    unsafe {

        let matrix_ptr = rssn::ffi_apis::physics_sim_linear_elasticity_ffi::handle::rssn_physics_sim_linear_elasticity_simulate_cantilever();

        assert!(!matrix_ptr.is_null());

        let matrix = &*matrix_ptr;

        assert!(matrix.rows() > 0);

        assert_eq!(matrix.cols(), 2);

        rssn::ffi_apis::numerical_matrix_ffi::handle::rssn_num_matrix_free(matrix_ptr);
    }
}

#[test]

fn test_linear_elasticity_json_ffi() {

    let input = r#"{
        "nodes": [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        "elements": [[0, 1, 2, 3]],
        "youngs_modulus": 1e7,
        "poissons_ratio": 0.3,
        "fixed_nodes": [0, 3],
        "loads": [[1, 1000.0, 0.0], [2, 1000.0, 0.0]]
    }"#;

    let c_input = CString::new(input).unwrap();

    unsafe {

        let res_ptr = rssn::ffi_apis::physics_sim_linear_elasticity_ffi::json::rssn_physics_sim_linear_elasticity_run_json(c_input.as_ptr());

        assert!(!res_ptr.is_null());

        let res_str = CStr::from_ptr(res_ptr).to_string_lossy();

        assert!(res_str.contains("\"ok\":"));

        rssn::ffi_apis::ffi_api::free_string(res_ptr);
    }
}
