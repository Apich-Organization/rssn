use crate::ffi_apis::common::{BincodeBuffer, from_bincode_buffer, to_bincode_buffer};
use crate::symbolic::tensor::Tensor;

/// Performs tensor addition.

///

/// Takes two bincode-serialized `Tensor` objects as input,

/// and returns a bincode-serialized `Tensor` representing their sum.

#[no_mangle]

pub extern "C" fn rssn_bincode_tensor_add(
    t1_buf: BincodeBuffer,
    t2_buf: BincodeBuffer,
) -> BincodeBuffer {

    let t1: Option<Tensor> =
        from_bincode_buffer(&t1_buf);

    let t2: Option<Tensor> =
        from_bincode_buffer(&t2_buf);

    if let (
        Some(tensor1),
        Some(tensor2),
    ) = (t1, t2)
    {

        match tensor1.add(&tensor2) {
            | Ok(result) => {
                to_bincode_buffer(
                    &result,
                )
            },
            | Err(_) => {
                BincodeBuffer::empty()
            },
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Performs scalar multiplication on a tensor.

///

/// Takes a bincode-serialized `Tensor` and a bincode-serialized `Expr` (scalar).

/// Returns a bincode-serialized `Tensor` representing the result.

#[no_mangle]

pub extern "C" fn rssn_bincode_tensor_scalar_mul(
    t_buf: BincodeBuffer,
    scalar_buf: BincodeBuffer,
) -> BincodeBuffer {

    let t: Option<Tensor> =
        from_bincode_buffer(&t_buf);

    let scalar: Option<
        crate::symbolic::core::Expr,
    > = from_bincode_buffer(
        &scalar_buf,
    );

    if let (Some(tensor), Some(s)) =
        (t, scalar)
    {

        match tensor.scalar_mul(&s) {
            | Ok(result) => {
                to_bincode_buffer(
                    &result,
                )
            },
            | Err(_) => {
                BincodeBuffer::empty()
            },
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the outer product of two tensors.

///

/// Takes two bincode-serialized `Tensor` objects as input,

/// and returns a bincode-serialized `Tensor` representing their outer product.

#[no_mangle]

pub extern "C" fn rssn_bincode_tensor_outer_product(
    t1_buf: BincodeBuffer,
    t2_buf: BincodeBuffer,
) -> BincodeBuffer {

    let t1: Option<Tensor> =
        from_bincode_buffer(&t1_buf);

    let t2: Option<Tensor> =
        from_bincode_buffer(&t2_buf);

    if let (
        Some(tensor1),
        Some(tensor2),
    ) = (t1, t2)
    {

        match tensor1
            .outer_product(&tensor2)
        {
            | Ok(result) => {
                to_bincode_buffer(
                    &result,
                )
            },
            | Err(_) => {
                BincodeBuffer::empty()
            },
        }
    } else {

        BincodeBuffer::empty()
    }
}
