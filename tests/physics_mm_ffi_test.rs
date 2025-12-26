//! FFI tests for the physics MM module.

use std::ffi::{CStr, CString};

#[test]
fn test_mm_simulate_dam_break_handle_ffi() {
    unsafe {
        let matrix_ptr = rssn::ffi_apis::physics_mm_ffi::handle::rssn_physics_mm_simulate_dam_break();
        assert!(!matrix_ptr.is_null());
        let matrix = &*matrix_ptr;
        assert_eq!(matrix.rows(), 200);
        assert_eq!(matrix.cols(), 2);
        rssn::ffi_apis::numerical_matrix_ffi::handle::rssn_num_matrix_free(matrix_ptr);
    }
}

#[test]
fn test_mm_sph_handle_life_cycle() {
    unsafe {
        let system_ptr = rssn::ffi_apis::physics_mm_ffi::handle::rssn_physics_mm_sph_create(0.1, 4.0, 4.0);
        assert!(!system_ptr.is_null());
        
        rssn::ffi_apis::physics_mm_ffi::handle::rssn_physics_mm_sph_add_particle(system_ptr, 1.0, 1.0, 0.0, 0.0, 1.0);
        rssn::ffi_apis::physics_mm_ffi::handle::rssn_physics_mm_sph_add_particle(system_ptr, 1.05, 1.0, 0.0, 0.0, 1.0);
        
        assert_eq!(rssn::ffi_apis::physics_mm_ffi::handle::rssn_physics_mm_sph_get_particle_count(system_ptr), 2);
        
        rssn::ffi_apis::physics_mm_ffi::handle::rssn_physics_mm_sph_update(system_ptr, 0.005);
        
        let pos_ptr = rssn::ffi_apis::physics_mm_ffi::handle::rssn_physics_mm_sph_get_positions(system_ptr);
        assert!(!pos_ptr.is_null());
        assert_eq!((*pos_ptr).rows(), 2);
        
        rssn::ffi_apis::numerical_matrix_ffi::handle::rssn_num_matrix_free(pos_ptr);
        rssn::ffi_apis::physics_mm_ffi::handle::rssn_physics_mm_sph_free(system_ptr);
    }
}
