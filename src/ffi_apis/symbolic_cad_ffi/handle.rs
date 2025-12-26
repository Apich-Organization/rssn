use crate::symbolic::cad::{cad, Cad};
use crate::symbolic::core::{Expr, SparsePolynomial};
use std::ffi::{CStr, CString};
use std::os::raw::c_char;

/// Computes CAD for a set of polynomials (Handle).
///
/// Expects an array of Expr handles (which must be SparsePolynomial variants)
/// and an array of variable name strings.
#[no_mangle]

pub extern "C" fn rssn_cad_handle(
    polys: *const *const Expr,
    polys_count: usize,
    vars: *const *const c_char,
    vars_count: usize,
) -> *mut Cad {

    if polys.is_null() || vars.is_null() {

        return std::ptr::null_mut();
    }

    let polys_slice = unsafe {

        std::slice::from_raw_parts(polys, polys_count)
    };

    let mut sparse_polys = Vec::new();

    for &p_ptr in polys_slice {

        if p_ptr.is_null() {

            return std::ptr::null_mut();
        }

        if let Expr::SparsePolynomial(sp) = unsafe {

            &*p_ptr
        } {

            sparse_polys.push(sp.clone());
        } else {

            return std::ptr::null_mut();
        }
    }

    let vars_slice = unsafe {

        std::slice::from_raw_parts(vars, vars_count)
    };

    let mut vars_vec = Vec::new();

    for &v_ptr in vars_slice {

        if v_ptr.is_null() {

            return std::ptr::null_mut();
        }

        let c_str = unsafe {

            CStr::from_ptr(v_ptr)
        };

        if let Ok(s) = c_str.to_str() {

            vars_vec.push(s);
        } else {

            return std::ptr::null_mut();
        }
    }

    match cad(
        &sparse_polys,
        &vars_vec,
    ) {
        Ok(c) => Box::into_raw(Box::new(c)),
        Err(_) => std::ptr::null_mut(),
    }
}

/// Frees a CAD handle.
#[no_mangle]

pub extern "C" fn rssn_free_cad_handle(ptr: *mut Cad) {

    if !ptr.is_null() {

        unsafe {

            let _ = Box::from_raw(ptr);
        }
    }
}

/// Gets the number of cells in a CAD.
#[no_mangle]

pub extern "C" fn rssn_cad_get_cell_count(ptr: *const Cad) -> usize {

    if ptr.is_null() {

        return 0;
    }

    unsafe {

        (*ptr).cells.len()
    }
}
