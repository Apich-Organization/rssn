//! FFI tests for the physics MTM (Multigrid) module.

use std::ffi::{
    CStr,
    CString,
};

#[test]

fn test_mtm_1d_handle_ffi() {

    let f = vec![-2.0; 31]; // 2^5-1
    let mut out_size = 0;

    unsafe {

        let res_ptr = rssn::ffi_apis::physics_mtm_ffi::handle::rssn_physics_mtm_solve_poisson_1d(
            31,
            f.as_ptr(),
            10,
            &mut out_size,
        );

        assert!(!res_ptr.is_null());

        assert_eq!(out_size, 33);

        rssn::ffi_apis::physics_mtm_ffi::handle::rssn_free_f64_mtm_array(res_ptr, out_size);
    }
}

#[test]

fn test_mtm_2d_json_ffi() {

    let n = 9; // 2^3+1
    let f = vec![0.0; n * n];

    let input = format!(
        r#"{{
        "n": {},
        "f": {:?},
        "num_cycles": 2
    }}"#,
        n, f
    );

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let res_ptr = rssn::ffi_apis::physics_mtm_ffi::json::rssn_physics_mtm_solve_poisson_2d_json(
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
