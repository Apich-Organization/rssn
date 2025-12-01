use crate::symbolic::core::Expr;
use crate::symbolic::differential_geometry::*;
use crate::symbolic::vector::Vector;
use std::ffi::CStr;
use std::os::raw::{c_char, c_int};

unsafe fn parse_c_str_array(arr: *const *const c_char, len: usize) -> Option<Vec<String>> {
    if arr.is_null() {
        return None;
    }
    let mut vars = Vec::with_capacity(len);
    for i in 0..len {
        let ptr = *arr.add(i);
        if ptr.is_null() {
            return None;
        }
        let c_str = CStr::from_ptr(ptr);
        match c_str.to_str() {
            Ok(s) => vars.push(s.to_string()),
            Err(_) => return None,
        }
    }
    Some(vars)
}

/// Computes the exterior derivative of a differential form (Handle)
#[no_mangle]
pub extern "C" fn rssn_exterior_derivative_handle(
    form_ptr: *const DifferentialForm,
    vars_ptr: *const *const c_char,
    vars_len: c_int,
) -> *mut DifferentialForm {
    if form_ptr.is_null() {
        return std::ptr::null_mut();
    }
    unsafe {
        let form = &*form_ptr;
        let vars_strings = match parse_c_str_array(vars_ptr, vars_len as usize) {
            Some(v) => v,
            None => return std::ptr::null_mut(),
        };
        let vars_refs: Vec<&str> = vars_strings.iter().map(|s| s.as_str()).collect();
        
        let result = exterior_derivative(form, &vars_refs);
        Box::into_raw(Box::new(result))
    }
}

/// Computes the wedge product of two differential forms (Handle)
#[no_mangle]
pub extern "C" fn rssn_wedge_product_handle(
    form1_ptr: *const DifferentialForm,
    form2_ptr: *const DifferentialForm,
) -> *mut DifferentialForm {
    if form1_ptr.is_null() || form2_ptr.is_null() {
        return std::ptr::null_mut();
    }
    unsafe {
        let form1 = &*form1_ptr;
        let form2 = &*form2_ptr;
        
        let result = wedge_product(form1, form2);
        Box::into_raw(Box::new(result))
    }
}

/// Computes the boundary of a domain (Handle)
#[no_mangle]
pub extern "C" fn rssn_boundary_handle(domain_ptr: *const Expr) -> *mut Expr {
    if domain_ptr.is_null() {
        return std::ptr::null_mut();
    }
    unsafe {
        let domain = &*domain_ptr;
        let result = boundary(domain);
        Box::into_raw(Box::new(result))
    }
}

/// Represents the generalized Stokes' theorem (Handle)
#[no_mangle]
pub extern "C" fn rssn_generalized_stokes_theorem_handle(
    omega_ptr: *const DifferentialForm,
    manifold_ptr: *const Expr,
    vars_ptr: *const *const c_char,
    vars_len: c_int,
) -> *mut Expr {
    if omega_ptr.is_null() || manifold_ptr.is_null() {
        return std::ptr::null_mut();
    }
    unsafe {
        let omega = &*omega_ptr;
        let manifold = &*manifold_ptr;
        let vars_strings = match parse_c_str_array(vars_ptr, vars_len as usize) {
            Some(v) => v,
            None => return std::ptr::null_mut(),
        };
        let vars_refs: Vec<&str> = vars_strings.iter().map(|s| s.as_str()).collect();
        
        let result = generalized_stokes_theorem(omega, manifold, &vars_refs);
        Box::into_raw(Box::new(result))
    }
}

/// Represents Gauss's theorem (Handle)
#[no_mangle]
pub extern "C" fn rssn_gauss_theorem_handle(
    vector_field_ptr: *const Vector,
    volume_ptr: *const Expr,
) -> *mut Expr {
    if vector_field_ptr.is_null() || volume_ptr.is_null() {
        return std::ptr::null_mut();
    }
    unsafe {
        let vector_field = &*vector_field_ptr;
        let volume = &*volume_ptr;
        
        let result = gauss_theorem(vector_field, volume);
        Box::into_raw(Box::new(result))
    }
}

/// Represents Stokes' theorem (Handle)
#[no_mangle]
pub extern "C" fn rssn_stokes_theorem_handle(
    vector_field_ptr: *const Vector,
    surface_ptr: *const Expr,
) -> *mut Expr {
    if vector_field_ptr.is_null() || surface_ptr.is_null() {
        return std::ptr::null_mut();
    }
    unsafe {
        let vector_field = &*vector_field_ptr;
        let surface = &*surface_ptr;
        
        let result = stokes_theorem(vector_field, surface);
        Box::into_raw(Box::new(result))
    }
}

/// Represents Green's theorem (Handle)
#[no_mangle]
pub extern "C" fn rssn_greens_theorem_handle(
    p_ptr: *const Expr,
    q_ptr: *const Expr,
    domain_ptr: *const Expr,
) -> *mut Expr {
    if p_ptr.is_null() || q_ptr.is_null() || domain_ptr.is_null() {
        return std::ptr::null_mut();
    }
    unsafe {
        let p = &*p_ptr;
        let q = &*q_ptr;
        let domain = &*domain_ptr;
        
        let result = greens_theorem(p, q, domain);
        Box::into_raw(Box::new(result))
    }
}

/// Frees a DifferentialForm handle
#[no_mangle]
pub extern "C" fn rssn_free_differential_form_handle(ptr: *mut DifferentialForm) {
    if !ptr.is_null() {
        unsafe {
            let _ = Box::from_raw(ptr);
        }
    }
}
