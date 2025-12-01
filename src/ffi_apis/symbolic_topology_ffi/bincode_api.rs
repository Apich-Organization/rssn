use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::topology::*;

/// Creates a new Simplex (Bincode)
#[no_mangle]
pub extern "C" fn rssn_bincode_simplex_create(
    vertices_buf: BincodeBuffer,
) -> BincodeBuffer {
    let vertices: Option<Vec<usize>> = from_bincode_buffer(&vertices_buf);
    if let Some(v) = vertices {
        let simplex = Simplex::new(&v);
        to_bincode_buffer(&simplex)
    } else {
        BincodeBuffer::empty()
    }
}

/// Gets the dimension of a Simplex (Bincode)
#[no_mangle]
pub extern "C" fn rssn_bincode_simplex_dimension(
    simplex_buf: BincodeBuffer,
) -> BincodeBuffer {
    let simplex: Option<Simplex> = from_bincode_buffer(&simplex_buf);
    if let Some(s) = simplex {
        to_bincode_buffer(&s.dimension())
    } else {
        BincodeBuffer::empty()
    }
}

/// Creates a new SimplicialComplex (Bincode)
#[no_mangle]
pub extern "C" fn rssn_bincode_simplicial_complex_create() -> BincodeBuffer {
    let complex = SimplicialComplex::new();
    to_bincode_buffer(&complex)
}

/// Adds a simplex to a SimplicialComplex (Bincode)
#[no_mangle]
pub extern "C" fn rssn_bincode_simplicial_complex_add_simplex(
    complex_buf: BincodeBuffer,
    vertices_buf: BincodeBuffer,
) -> BincodeBuffer {
    let complex: Option<SimplicialComplex> = from_bincode_buffer(&complex_buf);
    let vertices: Option<Vec<usize>> = from_bincode_buffer(&vertices_buf);
    
    if let (Some(mut c), Some(v)) = (complex, vertices) {
        c.add_simplex(&v);
        to_bincode_buffer(&c)
    } else {
        BincodeBuffer::empty()
    }
}

/// Gets the symbolic boundary matrix for dimension k (Bincode)
#[no_mangle]
pub extern "C" fn rssn_bincode_simplicial_complex_get_symbolic_boundary_matrix(
    complex_buf: BincodeBuffer,
    k: usize,
) -> BincodeBuffer {
    let complex: Option<SimplicialComplex> = from_bincode_buffer(&complex_buf);
    
    if let Some(c) = complex {
        match c.get_symbolic_boundary_matrix(k) {
            Some(expr) => to_bincode_buffer(&expr),
            None => BincodeBuffer::empty(),
        }
    } else {
        BincodeBuffer::empty()
    }
}

/// Creates a new SymbolicChain (Bincode)
#[no_mangle]
pub extern "C" fn rssn_bincode_symbolic_chain_create(
    dimension: usize,
) -> BincodeBuffer {
    let chain = SymbolicChain::new(dimension);
    to_bincode_buffer(&chain)
}

/// Adds a term to a SymbolicChain (Bincode)
#[no_mangle]
pub extern "C" fn rssn_bincode_symbolic_chain_add_term(
    chain_buf: BincodeBuffer,
    simplex_buf: BincodeBuffer,
    coeff_buf: BincodeBuffer,
) -> BincodeBuffer {
    let chain: Option<SymbolicChain> = from_bincode_buffer(&chain_buf);
    let simplex: Option<Simplex> = from_bincode_buffer(&simplex_buf);
    let coeff: Option<Expr> = from_bincode_buffer(&coeff_buf);
    
    if let (Some(mut c), Some(s), Some(coeff)) = (chain, simplex, coeff) {
        match c.add_term(s, coeff) {
            Ok(_) => to_bincode_buffer(&c),
            Err(_) => BincodeBuffer::empty(),
        }
    } else {
        BincodeBuffer::empty()
    }
}

/// Applies the symbolic boundary operator to a SymbolicChain (Bincode)
#[no_mangle]
pub extern "C" fn rssn_bincode_simplicial_complex_apply_symbolic_boundary_operator(
    complex_buf: BincodeBuffer,
    chain_buf: BincodeBuffer,
) -> BincodeBuffer {
    let complex: Option<SimplicialComplex> = from_bincode_buffer(&complex_buf);
    let chain: Option<SymbolicChain> = from_bincode_buffer(&chain_buf);
    
    if let (Some(c), Some(ch)) = (complex, chain) {
        match c.apply_symbolic_boundary_operator(&ch) {
            Some(result) => to_bincode_buffer(&result),
            None => BincodeBuffer::empty(),
        }
    } else {
        BincodeBuffer::empty()
    }
}
