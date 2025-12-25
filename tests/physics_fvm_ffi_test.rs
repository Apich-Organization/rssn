//! FFI tests for physics FVM module.

use std::ffi::CString;
use rssn::ffi_apis::ffi_api::FfiResult;

#[test]
fn test_fvm_handle_ffi() {
    unsafe {
        let mesh_ptr = rssn::ffi_apis::physics_fvm_ffi::handle::rssn_physics_fvm_mesh_new(50, 1.0);
        assert!(!mesh_ptr.is_null());
        let data_ptr = rssn::ffi_apis::physics_fvm_ffi::handle::rssn_physics_fvm_mesh_data(mesh_ptr);
        assert!(!data_ptr.is_null());
        rssn::ffi_apis::physics_fvm_ffi::handle::rssn_physics_fvm_mesh_free(mesh_ptr);
    }
}

#[test]
fn test_fvm_advection_json_ffi() {
    let input = r#"{
        "num_cells": 20,
        "domain_size": 1.0,
        "velocity": 1.0,
        "dt": 0.01,
        "steps": 10,
        "initial_values": [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    }"#;
    let c_input = CString::new(input).unwrap();
    unsafe {
        let res_ptr = rssn::ffi_apis::physics_fvm_ffi::json::rssn_physics_fvm_advection_json(c_input.as_ptr());
        assert!(!res_ptr.is_null());
        let res_str = std::ffi::CStr::from_ptr(res_ptr).to_string_lossy();
        assert!(res_str.contains("\"ok\":"));
        rssn::ffi_apis::ffi_api::free_string(res_ptr);
    }
}
