//! FFI tests for physics FDM module.

use rssn::physics::physics_fdm::FdmGrid;
use std::ffi::CString;

#[test]
fn test_fdm_handle_ffi() {
    unsafe {
        let grid_ptr = rssn::ffi_apis::physics_fdm_ffi::handle::rssn_physics_fdm_simulate_heat_2d();
        assert!(!grid_ptr.is_null());
        let grid = &*grid_ptr;
        assert!(grid.len() > 0);

        // Clean up
        rssn::ffi_apis::physics_fdm_ffi::handle::rssn_physics_fdm_grid_free(grid_ptr);
    }
}

#[test]
fn test_fdm_heat_json_ffi() {
    let input = r#"{
        "width": 20,
        "height": 20,
        "alpha": 0.01,
        "dx": 1.0,
        "dy": 1.0,
        "dt": 0.1,
        "steps": 10,
        "initial_temp": 100.0
    }"#;
    let c_input = CString::new(input).unwrap();
    unsafe {
        let res_ptr =
            rssn::ffi_apis::physics_fdm_ffi::json::rssn_physics_fdm_heat_json(c_input.as_ptr());
        assert!(!res_ptr.is_null());
        let res_str = std::ffi::CStr::from_ptr(res_ptr).to_string_lossy();
        assert!(res_str.contains("\"ok\":"));
        // Clean up
        rssn::ffi_apis::ffi_api::free_string(res_ptr);
    }
}

#[test]
fn test_fdm_burgers_json_ffi() {
    let input = r#"{
        "initial_u": [1.0, 1.0, 0.0, 0.0],
        "dx": 1.0,
        "nu": 0.1,
        "dt": 0.01,
        "steps": 5
    }"#;
    let c_input = CString::new(input).unwrap();
    unsafe {
        let res_ptr =
            rssn::ffi_apis::physics_fdm_ffi::json::rssn_physics_fdm_burgers_json(c_input.as_ptr());
        assert!(!res_ptr.is_null());
        let res_str = std::ffi::CStr::from_ptr(res_ptr).to_string_lossy();
        assert!(res_str.contains("\"ok\":"));
        rssn::ffi_apis::ffi_api::free_string(res_ptr);
    }
}
