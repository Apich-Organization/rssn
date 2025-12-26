//! FFI tests for the physics BEM module.

use std::ffi::CStr;
use std::ffi::CString;

use rssn::ffi_apis::ffi_api::FfiResult;

#[test]

fn test_bem_2d_json_ffi() {

    let input = r#"{
        "points": [[0,0], [1,0], [1,1], [0,1]],
        "bcs": [
            {"Flux": 0.0},
            {"Potential": 0.0},
            {"Flux": 0.0},
            {"Potential": 100.0}
        ]
    }"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let res_ptr = rssn::ffi_apis::physics_bem_ffi::json::rssn_physics_bem_solve_laplace_2d_json(
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

#[test]

fn test_bem_2d_handle_ffi() {

    let points_x = [0.0, 1.0, 1.0, 0.0];

    let points_y = [0.0, 0.0, 1.0, 1.0];

    let bcs_type = [1, 0, 1, 0]; // 1=Flux, 0=Potential
    let bcs_value =
        [0.0, 0.0, 0.0, 100.0];

    let mut out_u = [0.0; 4];

    let mut out_q = [0.0; 4];

    unsafe {

        let res = rssn::ffi_apis::physics_bem_ffi::handle::rssn_physics_bem_solve_laplace_2d(
            points_x.as_ptr(),
            points_y.as_ptr(),
            bcs_type.as_ptr(),
            bcs_value.as_ptr(),
            4,
            out_u.as_mut_ptr(),
            out_q.as_mut_ptr(),
        );

        assert_eq!(res, 0);

        assert_eq!(out_u[1], 0.0);

        assert_eq!(out_u[3], 100.0);
    }
}
