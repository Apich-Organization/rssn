//! FFI tests for the physics sim Ising statistical module.

use std::ffi::{CStr, CString};

#[test]
fn test_ising_handle_ffi() {
    unsafe {
        let handle = rssn::ffi_apis::physics_sim_ising_ffi::handle::rssn_physics_sim_ising_run(
            16, 16, 0.5, 10
        );
        
        assert!(!handle.grid.is_null());
        assert_eq!((*handle.grid).rows(), 16);
        assert_eq!((*handle.grid).cols(), 16);
        assert!(handle.magnetization >= 0.0 && handle.magnetization <= 1.0);
        
        rssn::ffi_apis::physics_sim_ising_ffi::handle::rssn_physics_sim_ising_free_result(handle);
    }
}

#[test]
fn test_ising_json_ffi() {
    let input = r#"{
        "width": 10,
        "height": 10,
        "temperature": 1.0,
        "mc_steps": 10
    }"#;
    let c_input = CString::new(input).unwrap();
    unsafe {
        let res_ptr = rssn::ffi_apis::physics_sim_ising_ffi::json::rssn_physics_sim_ising_run_json(c_input.as_ptr());
        assert!(!res_ptr.is_null());
        let res_str = CStr::from_ptr(res_ptr).to_string_lossy();
        assert!(res_str.contains("\"ok\":"));
        rssn::ffi_apis::ffi_api::free_string(res_ptr);
    }
}
