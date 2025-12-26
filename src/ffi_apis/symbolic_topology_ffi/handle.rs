use std::ffi::CStr;
use std::os::raw::c_char;
use std::os::raw::c_double;
use std::os::raw::c_int;

use crate::symbolic::core::Expr;
use crate::symbolic::topology::*;

// --- Simplex ---

/// Creates a new Simplex (Handle)
#[no_mangle]

pub extern "C" fn rssn_simplex_create(
    vertices_ptr : *const usize,
    len : usize,
) -> *mut Simplex {

    if vertices_ptr.is_null() {

        return std::ptr::null_mut();
    }

    unsafe {

        let vertices = std::slice::from_raw_parts(vertices_ptr, len);

        let simplex = Simplex::new(vertices);

        Box::into_raw(Box::new(simplex))
    }
}

/// Frees a Simplex handle
#[no_mangle]

pub extern "C" fn rssn_simplex_free(ptr : *mut Simplex) {

    if !ptr.is_null() {

        unsafe {

            let _ = Box::from_raw(ptr);
        }
    }
}

/// Gets the dimension of a Simplex
#[no_mangle]

pub extern "C" fn rssn_simplex_dimension(ptr : *const Simplex) -> usize {

    if ptr.is_null() {

        return 0;
    }

    unsafe {

        (*ptr).dimension()
    }
}

// --- SimplicialComplex ---

/// Creates a new SimplicialComplex (Handle)
#[no_mangle]

pub extern "C" fn rssn_simplicial_complex_create() -> *mut SimplicialComplex {

    let complex = SimplicialComplex::new();

    Box::into_raw(Box::new(complex))
}

/// Frees a SimplicialComplex handle
#[no_mangle]

pub extern "C" fn rssn_simplicial_complex_free(ptr : *mut SimplicialComplex) {

    if !ptr.is_null() {

        unsafe {

            let _ = Box::from_raw(ptr);
        }
    }
}

/// Adds a simplex to a SimplicialComplex
#[no_mangle]

pub extern "C" fn rssn_simplicial_complex_add_simplex(
    complex_ptr : *mut SimplicialComplex,
    vertices_ptr : *const usize,
    len : usize,
) {

    if complex_ptr.is_null() || vertices_ptr.is_null() {

        return;
    }

    unsafe {

        let complex = &mut *complex_ptr;

        let vertices = std::slice::from_raw_parts(vertices_ptr, len);

        complex.add_simplex(vertices);
    }
}

/// Gets the dimension of a SimplicialComplex
#[no_mangle]

pub extern "C" fn rssn_simplicial_complex_dimension(ptr : *const SimplicialComplex) -> c_int {

    if ptr.is_null() {

        return -1;
    }

    unsafe {

        match (*ptr).dimension() {
            | Some(d) => d as c_int,
            | None => -1,
        }
    }
}

/// Computes the Euler characteristic
#[no_mangle]

pub extern "C" fn rssn_simplicial_complex_euler_characteristic(
    ptr : *const SimplicialComplex
) -> isize {

    if ptr.is_null() {

        return 0;
    }

    unsafe {

        (*ptr).compute_euler_characteristic()
    }
}

/// Gets the symbolic boundary matrix for dimension k
#[no_mangle]

pub extern "C" fn rssn_simplicial_complex_get_symbolic_boundary_matrix(
    complex_ptr : *const SimplicialComplex,
    k : usize,
) -> *mut Expr {

    if complex_ptr.is_null() {

        return std::ptr::null_mut();
    }

    unsafe {

        match (*complex_ptr).get_symbolic_boundary_matrix(k) {
            | Some(expr) => Box::into_raw(Box::new(expr)),
            | None => std::ptr::null_mut(),
        }
    }
}

// --- SymbolicChain ---

/// Creates a new SymbolicChain (Handle)
#[no_mangle]

pub extern "C" fn rssn_symbolic_chain_create(dimension : usize) -> *mut SymbolicChain {

    let chain = SymbolicChain::new(dimension);

    Box::into_raw(Box::new(chain))
}

/// Frees a SymbolicChain handle
#[no_mangle]

pub extern "C" fn rssn_symbolic_chain_free(ptr : *mut SymbolicChain) {

    if !ptr.is_null() {

        unsafe {

            let _ = Box::from_raw(ptr);
        }
    }
}

/// Adds a term to a SymbolicChain
#[no_mangle]

pub extern "C" fn rssn_symbolic_chain_add_term(
    chain_ptr : *mut SymbolicChain,
    simplex_ptr : *const Simplex,
    coeff_ptr : *const Expr,
) -> bool {

    if chain_ptr.is_null() || simplex_ptr.is_null() || coeff_ptr.is_null() {

        return false;
    }

    unsafe {

        let chain = &mut *chain_ptr;

        let simplex = &*simplex_ptr;

        let coeff = &*coeff_ptr;

        match chain.add_term(
            simplex.clone(),
            coeff.clone(),
        ) {
            | Ok(_) => true,
            | Err(_) => false,
        }
    }
}

/// Applies the symbolic boundary operator to a SymbolicChain
#[no_mangle]

pub extern "C" fn rssn_simplicial_complex_apply_symbolic_boundary_operator(
    complex_ptr : *const SimplicialComplex,
    chain_ptr : *const SymbolicChain,
) -> *mut SymbolicChain {

    if complex_ptr.is_null() || chain_ptr.is_null() {

        return std::ptr::null_mut();
    }

    unsafe {

        let complex = &*complex_ptr;

        let chain = &*chain_ptr;

        match complex.apply_symbolic_boundary_operator(chain) {
            | Some(result) => Box::into_raw(Box::new(result)),
            | None => std::ptr::null_mut(),
        }
    }
}

// --- Generators ---

/// Creates a grid complex
#[no_mangle]

pub extern "C" fn rssn_create_grid_complex(
    width : usize,
    height : usize,
) -> *mut SimplicialComplex {

    let complex = create_grid_complex(width, height);

    Box::into_raw(Box::new(complex))
}

/// Creates a torus complex
#[no_mangle]

pub extern "C" fn rssn_create_torus_complex(
    m : usize,
    n : usize,
) -> *mut SimplicialComplex {

    let complex = create_torus_complex(m, n);

    Box::into_raw(Box::new(complex))
}
