//! FFI tests for the physics sim Navier-Stokes module.

use std::ffi::{
    CStr,
    CString,
};

#[test]

fn test_navier_stokes_handle_ffi() {

    unsafe {

        let handles = rssn::ffi_apis::physics_sim_navier_stokes_ffi::handle::rssn_physics_sim_navier_stokes_run_lid_driven_cavity(
            9, 9, 100.0, 0.01, 10, 1.0
        );

        assert!(!handles.u.is_null());

        assert!(!handles.v.is_null());

        assert!(!handles.p.is_null());

        assert_eq!(
            (*handles.u).rows(),
            9
        );

        assert_eq!(
            (*handles.u).cols(),
            9
        );

        rssn::ffi_apis::physics_sim_navier_stokes_ffi::handle::rssn_physics_sim_navier_stokes_free_results(handles);
    }
}

#[test]

fn test_navier_stokes_json_ffi() {

    let input = r#"{
        "nx": 9,
        "ny": 9,
        "re": 100.0,
        "dt": 0.01,
        "n_iter": 5,
        "lid_velocity": 1.0
    }"#;

    let c_input = CString::new(input).unwrap();

    unsafe {

        let res_ptr = rssn::ffi_apis::physics_sim_navier_stokes_ffi::json::rssn_physics_sim_navier_stokes_run_json(c_input.as_ptr());

        assert!(!res_ptr.is_null());

        let res_str = CStr::from_ptr(res_ptr).to_string_lossy();

        assert!(res_str.contains("\"ok\":"));

        rssn::ffi_apis::ffi_api::free_string(res_ptr);
    }
}
