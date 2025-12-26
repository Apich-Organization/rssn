use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::topology::*;

/// Creates a new Simplex (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_simplex_create(
    vertices_json: *const std::os::raw::c_char,
) -> *mut std::os::raw::c_char {

    let vertices: Option<Vec<usize>> = from_json_string(vertices_json);

    if let Some(v) = vertices {

        let simplex = Simplex::new(&v);

        to_json_string(&simplex)
    } else {

        std::ptr::null_mut()
    }
}

/// Gets the dimension of a Simplex (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_simplex_dimension(
    simplex_json: *const std::os::raw::c_char,
) -> *mut std::os::raw::c_char {

    let simplex: Option<Simplex> = from_json_string(simplex_json);

    if let Some(s) = simplex {

        to_json_string(&s.dimension())
    } else {

        std::ptr::null_mut()
    }
}

/// Creates a new SimplicialComplex (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_simplicial_complex_create() -> *mut std::os::raw::c_char {

    let complex = SimplicialComplex::new();

    to_json_string(&complex)
}

/// Adds a simplex to a SimplicialComplex (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_simplicial_complex_add_simplex(
    complex_json: *const std::os::raw::c_char,
    vertices_json: *const std::os::raw::c_char,
) -> *mut std::os::raw::c_char {

    let complex: Option<SimplicialComplex> = from_json_string(complex_json);

    let vertices: Option<Vec<usize>> = from_json_string(vertices_json);

    if let (Some(mut c), Some(v)) = (complex, vertices) {

        c.add_simplex(&v);

        to_json_string(&c)
    } else {

        std::ptr::null_mut()
    }
}

/// Gets the symbolic boundary matrix for dimension k (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_simplicial_complex_get_symbolic_boundary_matrix(
    complex_json: *const std::os::raw::c_char,
    k: usize,
) -> *mut std::os::raw::c_char {

    let complex: Option<SimplicialComplex> = from_json_string(complex_json);

    if let Some(c) = complex {

        match c.get_symbolic_boundary_matrix(k) {
            Some(expr) => to_json_string(&expr),
            None => std::ptr::null_mut(),
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Creates a new SymbolicChain (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_symbolic_chain_create(dimension: usize) -> *mut std::os::raw::c_char {

    let chain = SymbolicChain::new(dimension);

    to_json_string(&chain)
}

/// Adds a term to a SymbolicChain (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_symbolic_chain_add_term(
    chain_json: *const std::os::raw::c_char,
    simplex_json: *const std::os::raw::c_char,
    coeff_json: *const std::os::raw::c_char,
) -> *mut std::os::raw::c_char {

    let chain: Option<SymbolicChain> = from_json_string(chain_json);

    let simplex: Option<Simplex> = from_json_string(simplex_json);

    let coeff: Option<Expr> = from_json_string(coeff_json);

    if let (Some(mut c), Some(s), Some(coeff)) = (chain, simplex, coeff) {

        match c.add_term(s, coeff) {
            Ok(_) => to_json_string(&c),
            Err(_) => std::ptr::null_mut(),
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Applies the symbolic boundary operator to a SymbolicChain (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_simplicial_complex_apply_symbolic_boundary_operator(
    complex_json: *const std::os::raw::c_char,
    chain_json: *const std::os::raw::c_char,
) -> *mut std::os::raw::c_char {

    let complex: Option<SimplicialComplex> = from_json_string(complex_json);

    let chain: Option<SymbolicChain> = from_json_string(chain_json);

    if let (Some(c), Some(ch)) = (complex, chain) {

        match c.apply_symbolic_boundary_operator(&ch) {
            Some(result) => to_json_string(&result),
            None => std::ptr::null_mut(),
        }
    } else {

        std::ptr::null_mut()
    }
}
