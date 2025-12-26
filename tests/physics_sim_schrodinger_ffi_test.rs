//! FFI tests for the physics sim Schrodinger quantum module.

use std::ffi::{CStr, CString};
use num_complex::Complex;

#[test]
fn test_schrodinger_handle_ffi() {
    unsafe {
        let nx = 16;
        let ny = 16;
        let n = nx * ny;
        let potential = vec![0.0; n];
        let psi_re = vec![1.0; n];
        let psi_im = vec![0.0; n];
        
        let matrix_ptr = rssn::ffi_apis::physics_sim_schrodinger_ffi::handle::rssn_physics_sim_schrodinger_run_2d(
            nx, ny, 10.0, 10.0, 0.1, 10, 1.0, 1.0,
            potential.as_ptr(), psi_re.as_ptr(), psi_im.as_ptr()
        );
        
        assert!(!matrix_ptr.is_null());
        let matrix = &*matrix_ptr;
        assert_eq!(matrix.rows(), ny);
        assert_eq!(matrix.cols(), nx);
        
        rssn::ffi_apis::numerical_matrix_ffi::handle::rssn_num_matrix_free(matrix_ptr);
    }
}

#[test]
fn test_schrodinger_json_ffi() {
    let nx = 8;
    let ny = 8;
    let n = nx * ny;
    let potential = vec![0.0; n];
    let psi_re = vec![1.0; n];
    let psi_im = vec![0.0; n];

    let input = format!(r#"{{
        "params": {{
            "nx": {},
            "ny": {},
            "lx": 10.0,
            "ly": 10.0,
            "dt": 0.1,
            "time_steps": 5,
            "hbar": 1.0,
            "mass": 1.0,
            "potential": {:?}
        }},
        "initial_psi_re": {:?},
        "initial_psi_im": {:?}
    }}"#, nx, ny, potential, psi_re, psi_im);
    
    let c_input = CString::new(input).unwrap();
    unsafe {
        let res_ptr = rssn::ffi_apis::physics_sim_schrodinger_ffi::json::rssn_physics_sim_schrodinger_run_json(c_input.as_ptr());
        assert!(!res_ptr.is_null());
        let res_str = CStr::from_ptr(res_ptr).to_string_lossy();
        assert!(res_str.contains("\"ok\":"));
        rssn::ffi_apis::ffi_api::free_string(res_ptr);
    }
}
