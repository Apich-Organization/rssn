use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::geometric_algebra::Multivector;

/// Creates a new scalar multivector (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_multivector_scalar(
    p : u32,
    q : u32,
    r : u32,
    value_buf : BincodeBuffer,
) -> BincodeBuffer {

    let value : Option<Expr> = from_bincode_buffer(&value_buf);

    if let Some(val) = value {

        let mv = Multivector::scalar((p, q, r), val);

        to_bincode_buffer(&mv)
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes geometric product (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_multivector_geometric_product(
    a_buf : BincodeBuffer,
    b_buf : BincodeBuffer,
) -> BincodeBuffer {

    let a : Option<Multivector> = from_bincode_buffer(&a_buf);

    let b : Option<Multivector> = from_bincode_buffer(&b_buf);

    if let (Some(mv_a), Some(mv_b)) = (a, b) {

        let result = mv_a.geometric_product(&mv_b);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes outer product (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_multivector_outer_product(
    a_buf : BincodeBuffer,
    b_buf : BincodeBuffer,
) -> BincodeBuffer {

    let a : Option<Multivector> = from_bincode_buffer(&a_buf);

    let b : Option<Multivector> = from_bincode_buffer(&b_buf);

    if let (Some(mv_a), Some(mv_b)) = (a, b) {

        let result = mv_a.outer_product(&mv_b);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes inner product (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_multivector_inner_product(
    a_buf : BincodeBuffer,
    b_buf : BincodeBuffer,
) -> BincodeBuffer {

    let a : Option<Multivector> = from_bincode_buffer(&a_buf);

    let b : Option<Multivector> = from_bincode_buffer(&b_buf);

    if let (Some(mv_a), Some(mv_b)) = (a, b) {

        let result = mv_a.inner_product(&mv_b);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes reverse (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_multivector_reverse(mv_buf : BincodeBuffer) -> BincodeBuffer {

    let mv : Option<Multivector> = from_bincode_buffer(&mv_buf);

    if let Some(multivector) = mv {

        let result = multivector.reverse();

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes grade projection (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_multivector_grade_projection(
    mv_buf : BincodeBuffer,
    grade : u32,
) -> BincodeBuffer {

    let mv : Option<Multivector> = from_bincode_buffer(&mv_buf);

    if let Some(multivector) = mv {

        let result = multivector.grade_projection(grade);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes magnitude (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_multivector_magnitude(mv_buf : BincodeBuffer) -> BincodeBuffer {

    let mv : Option<Multivector> = from_bincode_buffer(&mv_buf);

    if let Some(multivector) = mv {

        let result = multivector.magnitude();

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}
