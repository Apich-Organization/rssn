//! Bincode-based FFI API for symbolic statistics functions.

use crate::ffi_apis::common::{BincodeBuffer, from_bincode_buffer, to_bincode_buffer};
use crate::symbolic::core::Expr;
use crate::symbolic::stats;

/// Computes the symbolic mean of a set of expressions using Bincode.
#[no_mangle]

pub extern "C" fn rssn_bincode_mean(
    data_buf: BincodeBuffer
) -> BincodeBuffer {

    let data: Option<Vec<Expr>> =
        from_bincode_buffer(&data_buf);

    if let Some(d) = data {

        to_bincode_buffer(&stats::mean(
            &d,
        ))
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the symbolic variance of a set of expressions using Bincode.
#[no_mangle]

pub extern "C" fn rssn_bincode_variance(
    data_buf: BincodeBuffer
) -> BincodeBuffer {

    let data: Option<Vec<Expr>> =
        from_bincode_buffer(&data_buf);

    if let Some(d) = data {

        to_bincode_buffer(
            &stats::variance(&d),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the symbolic standard deviation of a set of expressions using Bincode.
#[no_mangle]

pub extern "C" fn rssn_bincode_std_dev(
    data_buf: BincodeBuffer
) -> BincodeBuffer {

    let data: Option<Vec<Expr>> =
        from_bincode_buffer(&data_buf);

    if let Some(d) = data {

        to_bincode_buffer(
            &stats::std_dev(&d),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the symbolic covariance of two sets of expressions using Bincode.
#[no_mangle]

pub extern "C" fn rssn_bincode_covariance(
    data1_buf: BincodeBuffer,
    data2_buf: BincodeBuffer,
) -> BincodeBuffer {

    let data1: Option<Vec<Expr>> =
        from_bincode_buffer(&data1_buf);

    let data2: Option<Vec<Expr>> =
        from_bincode_buffer(&data2_buf);

    if let (Some(d1), Some(d2)) =
        (data1, data2)
    {

        to_bincode_buffer(
            &stats::covariance(
                &d1, &d2,
            ),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the symbolic Pearson correlation coefficient using Bincode.
#[no_mangle]

pub extern "C" fn rssn_bincode_correlation(
    data1_buf: BincodeBuffer,
    data2_buf: BincodeBuffer,
) -> BincodeBuffer {

    let data1: Option<Vec<Expr>> =
        from_bincode_buffer(&data1_buf);

    let data2: Option<Vec<Expr>> =
        from_bincode_buffer(&data2_buf);

    if let (Some(d1), Some(d2)) =
        (data1, data2)
    {

        to_bincode_buffer(
            &stats::correlation(
                &d1, &d2,
            ),
        )
    } else {

        BincodeBuffer::empty()
    }
}
