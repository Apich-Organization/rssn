//! FFI tests for physics FEM module.

use std::ffi::CString;

#[test]

fn test_fem_1d_json_ffi() {

    let input = r#"{
        "n_elements": 10,
        "domain_length": 1.0
    }"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let res_ptr = rssn::ffi_apis::physics_fem_ffi::json::rssn_physics_fem_solve_poisson_1d_json(
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

        rssn::ffi_apis::ffi_api::free_string(res_ptr);
    }
}

#[test]

fn test_fem_1d_handle_ffi() {

    let mut out_size = 0;

    unsafe {

        let res_ptr = rssn::ffi_apis::physics_fem_ffi::handle::rssn_physics_fem_solve_poisson_1d(
            10,
            1.0,
            &mut out_size,
        );

        assert!(!res_ptr.is_null());

        assert_eq!(out_size, 11);

        rssn::ffi_apis::physics_fem_ffi::handle::rssn_free_f64_array(res_ptr, out_size);
    }
}
