use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::symbolic::tensor::Tensor;

/// Performs tensor addition.

///

/// Takes two bincode-serialized `Tensor` objects as input,

/// and returns a bincode-serialized `Tensor` representing their sum.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_tensor_add(
    t1_buf: BincodeBuffer,
    t2_buf: BincodeBuffer,
) -> BincodeBuffer {

    let t1: Option<Tensor> =
        from_bincode_buffer(&t1_buf);

    let t2: Option<Tensor> =
        from_bincode_buffer(&t2_buf);

    match (t1, t2) {
        | (
            Some(tensor1),
            Some(tensor2),
        ) => {

            match tensor1.add(&tensor2)
            {
                | Ok(result) => {
                    to_bincode_buffer(
                        &result,
                    )
                },
                | Err(_) => {
                    BincodeBuffer::empty(
                    )
                },
            }
        },
        | _ => BincodeBuffer::empty(),
    }
}

/// Performs scalar multiplication on a tensor.

///

/// Takes a bincode-serialized `Tensor` and a bincode-serialized `Expr` (scalar).

/// Returns a bincode-serialized `Tensor` representing the result.

#[unsafe(no_mangle)]

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

    match (t, scalar) {
        | (Some(tensor), Some(s)) => {
            match tensor.scalar_mul(&s)
            {
                | Ok(result) => {
                    to_bincode_buffer(
                        &result,
                    )
                },
                | Err(_) => {
                    BincodeBuffer::empty(
                    )
                },
            }
        },
        | _ => BincodeBuffer::empty(),
    }
}

/// Computes the outer product of two tensors.

///

/// Takes two bincode-serialized `Tensor` objects as input,

/// and returns a bincode-serialized `Tensor` representing their outer product.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_tensor_outer_product(
    t1_buf: BincodeBuffer,
    t2_buf: BincodeBuffer,
) -> BincodeBuffer {

    let t1: Option<Tensor> =
        from_bincode_buffer(&t1_buf);

    let t2: Option<Tensor> =
        from_bincode_buffer(&t2_buf);

    match (t1, t2) {
        | (
            Some(tensor1),
            Some(tensor2),
        ) => {

            match tensor1
                .outer_product(&tensor2)
            {
                | Ok(result) => {
                    to_bincode_buffer(
                        &result,
                    )
                },
                | Err(_) => {
                    BincodeBuffer::empty(
                    )
                },
            }
        },
        | _ => BincodeBuffer::empty(),
    }
}
