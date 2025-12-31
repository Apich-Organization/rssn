//! Bincode-based FFI API for symbolic integral transforms.
//!
//! This module provides binary serialization-based FFI functions for Fourier, Laplace, and Z-transforms,
//! offering efficient binary data interchange for high-performance applications.

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::symbolic::core::Expr;
use crate::symbolic::transforms;

// --- Fourier Transform ---

/// Computes the Fourier transform of an expression.

///

/// Takes bincode-serialized `Expr` (expression), `String` (input variable), and `String` (output variable).

/// Returns a bincode-serialized `Expr` representing the Fourier transform.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_fourier_transform(
    expr_buf: BincodeBuffer,
    in_var_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let in_var: Option<String> =
        from_bincode_buffer(
            &in_var_buf,
        );

    let out_var: Option<String> =
        from_bincode_buffer(
            &out_var_buf,
        );

    match (
        expr,
        in_var,
        out_var,
    ) { (
        Some(e),
        Some(iv),
        Some(ov),
    ) => {

        to_bincode_buffer(&transforms::fourier_transform(&e, &iv, &ov))
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Computes the inverse Fourier transform of an expression.

///

/// Takes bincode-serialized `Expr` (expression), `String` (input variable), and `String` (output variable).

/// Returns a bincode-serialized `Expr` representing the inverse Fourier transform.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_inverse_fourier_transform(
    expr_buf: BincodeBuffer,
    in_var_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let in_var: Option<String> =
        from_bincode_buffer(
            &in_var_buf,
        );

    let out_var: Option<String> =
        from_bincode_buffer(
            &out_var_buf,
        );

    match (
        expr,
        in_var,
        out_var,
    ) { (
        Some(e),
        Some(iv),
        Some(ov),
    ) => {

        to_bincode_buffer(&transforms::inverse_fourier_transform(&e, &iv, &ov))
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Applies the time shift property of the Fourier transform.

///

/// Takes bincode-serialized `Expr` (frequency domain expression), `Expr` (time shift amount),

/// and `String` (output variable).

/// Returns a bincode-serialized `Expr` representing the transformed expression.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_fourier_time_shift(
    f_omega_buf: BincodeBuffer,
    a_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> =
        from_bincode_buffer(
            &f_omega_buf,
        );

    let a: Option<Expr> =
        from_bincode_buffer(&a_buf);

    let out_var: Option<String> =
        from_bincode_buffer(
            &out_var_buf,
        );

    match (f, a, out_var)
    { (
        Some(f),
        Some(a),
        Some(ov),
    ) => {

        to_bincode_buffer(&transforms::fourier_time_shift(&f, &a, &ov))
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Applies the frequency shift property of the Fourier transform.

///

/// Takes bincode-serialized `Expr` (frequency domain expression), `Expr` (frequency shift amount),

/// and `String` (output variable).

/// Returns a bincode-serialized `Expr` representing the transformed expression.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_fourier_frequency_shift(
    f_omega_buf: BincodeBuffer,
    a_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> =
        from_bincode_buffer(
            &f_omega_buf,
        );

    let a: Option<Expr> =
        from_bincode_buffer(&a_buf);

    let out_var: Option<String> =
        from_bincode_buffer(
            &out_var_buf,
        );

    match (f, a, out_var)
    { (
        Some(f),
        Some(a),
        Some(ov),
    ) => {

        to_bincode_buffer(&transforms::fourier_frequency_shift(&f, &a, &ov))
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Applies the scaling property of the Fourier transform.

///

/// Takes bincode-serialized `Expr` (frequency domain expression), `Expr` (scaling factor `a`),

/// and `String` (output variable).

/// Returns a bincode-serialized `Expr` representing the transformed expression.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_fourier_scaling(
    f_omega_buf: BincodeBuffer,
    a_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> =
        from_bincode_buffer(
            &f_omega_buf,
        );

    let a: Option<Expr> =
        from_bincode_buffer(&a_buf);

    let out_var: Option<String> =
        from_bincode_buffer(
            &out_var_buf,
        );

    match (f, a, out_var)
    { (
        Some(f),
        Some(a),
        Some(ov),
    ) => {

        to_bincode_buffer(&transforms::fourier_scaling(&f, &a, &ov))
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Applies the differentiation property of the Fourier transform.

///

/// Takes bincode-serialized `Expr` (frequency domain expression) and `String` (output variable).

/// Returns a bincode-serialized `Expr` representing the transformed expression.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_fourier_differentiation(
    f_omega_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> =
        from_bincode_buffer(
            &f_omega_buf,
        );

    let out_var: Option<String> =
        from_bincode_buffer(
            &out_var_buf,
        );

    match (f, out_var)
    { (Some(f), Some(ov)) => {

        to_bincode_buffer(&transforms::fourier_differentiation(&f, &ov))
    } _ => {

        BincodeBuffer::empty()
    }}
}

// --- Laplace Transform ---

/// Computes the Laplace transform of an expression.

///

/// Takes bincode-serialized `Expr` (expression), `String` (input variable), and `String` (output variable).

/// Returns a bincode-serialized `Expr` representing the Laplace transform.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_laplace_transform(
    expr_buf: BincodeBuffer,
    in_var_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let in_var: Option<String> =
        from_bincode_buffer(
            &in_var_buf,
        );

    let out_var: Option<String> =
        from_bincode_buffer(
            &out_var_buf,
        );

    match (
        expr,
        in_var,
        out_var,
    ) { (
        Some(e),
        Some(iv),
        Some(ov),
    ) => {

        to_bincode_buffer(&transforms::laplace_transform(&e, &iv, &ov))
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Computes the inverse Laplace transform of an expression.

///

/// Takes bincode-serialized `Expr` (expression), `String` (input variable), and `String` (output variable).

/// Returns a bincode-serialized `Expr` representing the inverse Laplace transform.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_inverse_laplace_transform(
    expr_buf: BincodeBuffer,
    in_var_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let in_var: Option<String> =
        from_bincode_buffer(
            &in_var_buf,
        );

    let out_var: Option<String> =
        from_bincode_buffer(
            &out_var_buf,
        );

    match (
        expr,
        in_var,
        out_var,
    ) { (
        Some(e),
        Some(iv),
        Some(ov),
    ) => {

        to_bincode_buffer(&transforms::inverse_laplace_transform(&e, &iv, &ov))
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Applies the time shift property of the Laplace transform.

///

/// Takes bincode-serialized `Expr` (s-domain expression), `Expr` (time shift amount `a`),

/// and `String` (output variable).

/// Returns a bincode-serialized `Expr` representing the transformed expression.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_laplace_time_shift(
    f_s_buf: BincodeBuffer,
    a_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> =
        from_bincode_buffer(&f_s_buf);

    let a: Option<Expr> =
        from_bincode_buffer(&a_buf);

    let out_var: Option<String> =
        from_bincode_buffer(
            &out_var_buf,
        );

    match (f, a, out_var)
    { (
        Some(f),
        Some(a),
        Some(ov),
    ) => {

        to_bincode_buffer(&transforms::laplace_time_shift(&f, &a, &ov))
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Applies the frequency shift property of the Laplace transform.

///

/// Takes bincode-serialized `Expr` (s-domain expression), `Expr` (frequency shift amount `a`),

/// and `String` (output variable).

/// Returns a bincode-serialized `Expr` representing the transformed expression.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_laplace_frequency_shift(
    f_s_buf: BincodeBuffer,
    a_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> =
        from_bincode_buffer(&f_s_buf);

    let a: Option<Expr> =
        from_bincode_buffer(&a_buf);

    let out_var: Option<String> =
        from_bincode_buffer(
            &out_var_buf,
        );

    match (f, a, out_var)
    { (
        Some(f),
        Some(a),
        Some(ov),
    ) => {

        to_bincode_buffer(&transforms::laplace_frequency_shift(&f, &a, &ov))
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Applies the scaling property of the Laplace transform.

///

/// Takes bincode-serialized `Expr` (s-domain expression), `Expr` (scaling factor `a`),

/// and `String` (output variable).

/// Returns a bincode-serialized `Expr` representing the transformed expression.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_laplace_scaling(
    f_s_buf: BincodeBuffer,
    a_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> =
        from_bincode_buffer(&f_s_buf);

    let a: Option<Expr> =
        from_bincode_buffer(&a_buf);

    let out_var: Option<String> =
        from_bincode_buffer(
            &out_var_buf,
        );

    match (f, a, out_var)
    { (
        Some(f),
        Some(a),
        Some(ov),
    ) => {

        to_bincode_buffer(&transforms::laplace_scaling(&f, &a, &ov))
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Applies the differentiation property of the Laplace transform.

///

/// Takes bincode-serialized `Expr` (s-domain expression), `String` (output variable),

/// and `Expr` (`f(0)` - initial condition).

/// Returns a bincode-serialized `Expr` representing the transformed expression.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_laplace_differentiation(
    f_s_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
    f_zero_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> =
        from_bincode_buffer(&f_s_buf);

    let out_var: Option<String> =
        from_bincode_buffer(
            &out_var_buf,
        );

    let f_zero: Option<Expr> =
        from_bincode_buffer(
            &f_zero_buf,
        );

    match (f, out_var, f_zero)
    { (
        Some(f),
        Some(ov),
        Some(fz),
    ) => {

        to_bincode_buffer(&transforms::laplace_differentiation(&f, &ov, &fz))
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Applies the integration property of the Laplace transform.

///

/// Takes bincode-serialized `Expr` (s-domain expression) and `String` (output variable).

/// Returns a bincode-serialized `Expr` representing the transformed expression.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_laplace_integration(
    f_s_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> =
        from_bincode_buffer(&f_s_buf);

    let out_var: Option<String> =
        from_bincode_buffer(
            &out_var_buf,
        );

    match (f, out_var)
    { (Some(f), Some(ov)) => {

        to_bincode_buffer(&transforms::laplace_integration(&f, &ov))
    } _ => {

        BincodeBuffer::empty()
    }}
}

// --- Z-Transform ---

/// Computes the Z-transform of an expression.

///

/// Takes bincode-serialized `Expr` (expression), `String` (input variable), and `String` (output variable).

/// Returns a bincode-serialized `Expr` representing the Z-transform.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_z_transform(
    expr_buf: BincodeBuffer,
    in_var_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let in_var: Option<String> =
        from_bincode_buffer(
            &in_var_buf,
        );

    let out_var: Option<String> =
        from_bincode_buffer(
            &out_var_buf,
        );

    match (
        expr,
        in_var,
        out_var,
    ) {
        | (
            Some(e),
            Some(iv),
            Some(ov),
        ) => to_bincode_buffer(
            &transforms::z_transform(
                &e, &iv, &ov,
            ),
        ),
        | _ => BincodeBuffer::empty(),
    }
}

/// Computes the inverse Z-transform of an expression.

///

/// Takes bincode-serialized `Expr` (expression), `String` (input variable), and `String` (output variable).

/// Returns a bincode-serialized `Expr` representing the inverse Z-transform.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_inverse_z_transform(
    expr_buf: BincodeBuffer,
    in_var_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let in_var: Option<String> =
        from_bincode_buffer(
            &in_var_buf,
        );

    let out_var: Option<String> =
        from_bincode_buffer(
            &out_var_buf,
        );

    match (
        expr,
        in_var,
        out_var,
    ) { (
        Some(e),
        Some(iv),
        Some(ov),
    ) => {

        to_bincode_buffer(&transforms::inverse_z_transform(&e, &iv, &ov))
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Applies the time shift property of the Z-transform.

///

/// Takes bincode-serialized `Expr` (z-domain expression), `Expr` (time shift amount `k`),

/// and `String` (output variable).

/// Returns a bincode-serialized `Expr` representing the transformed expression.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_z_time_shift(
    f_z_buf: BincodeBuffer,
    k_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> =
        from_bincode_buffer(&f_z_buf);

    let k: Option<Expr> =
        from_bincode_buffer(&k_buf);

    let out_var: Option<String> =
        from_bincode_buffer(
            &out_var_buf,
        );

    match (f, k, out_var) {
        | (
            Some(f),
            Some(k),
            Some(ov),
        ) => to_bincode_buffer(
            &transforms::z_time_shift(
                &f, &k, &ov,
            ),
        ),
        | _ => BincodeBuffer::empty(),
    }
}

/// Applies the scaling property of the Z-transform.

///

/// Takes bincode-serialized `Expr` (z-domain expression), `Expr` (scaling factor `a`),

/// and `String` (output variable).

/// Returns a bincode-serialized `Expr` representing the transformed expression.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_z_scaling(
    f_z_buf: BincodeBuffer,
    a_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> =
        from_bincode_buffer(&f_z_buf);

    let a: Option<Expr> =
        from_bincode_buffer(&a_buf);

    let out_var: Option<String> =
        from_bincode_buffer(
            &out_var_buf,
        );

    match (f, a, out_var) {
        | (
            Some(f),
            Some(a),
            Some(ov),
        ) => {
            to_bincode_buffer(
                &transforms::z_scaling(
                    &f, &a, &ov,
                ),
            )
        },
        | _ => BincodeBuffer::empty(),
    }
}

/// Applies the differentiation property of the Z-transform.

///

/// Takes bincode-serialized `Expr` (z-domain expression) and `String` (output variable).

/// Returns a bincode-serialized `Expr` representing the transformed expression.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_z_differentiation(
    f_z_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> =
        from_bincode_buffer(&f_z_buf);

    let out_var: Option<String> =
        from_bincode_buffer(
            &out_var_buf,
        );

    match (f, out_var)
    { (Some(f), Some(ov)) => {

        to_bincode_buffer(&transforms::z_differentiation(&f, &ov))
    } _ => {

        BincodeBuffer::empty()
    }}
}

// --- Utils ---

/// Computes the convolution of two functions using the Fourier transform property.

///

/// Takes bincode-serialized `Expr` representing two functions (`f` and `g`),

/// and `String`s for input and output variables.

/// Returns a bincode-serialized `Expr` representing the convolution.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_convolution_fourier(
    f_buf: BincodeBuffer,
    g_buf: BincodeBuffer,
    in_var_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> =
        from_bincode_buffer(&f_buf);

    let g: Option<Expr> =
        from_bincode_buffer(&g_buf);

    let in_var: Option<String> =
        from_bincode_buffer(
            &in_var_buf,
        );

    let out_var: Option<String> =
        from_bincode_buffer(
            &out_var_buf,
        );

    match (
        f,
        g,
        in_var,
        out_var,
    ) { (
        Some(f),
        Some(g),
        Some(iv),
        Some(ov),
    ) => {

        to_bincode_buffer(&transforms::convolution_fourier(&f, &g, &iv, &ov))
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Computes the convolution of two functions using the Laplace transform property.

///

/// Takes bincode-serialized `Expr` representing two functions (`f` and `g`),

/// and `String`s for input and output variables.

/// Returns a bincode-serialized `Expr` representing the convolution.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_convolution_laplace(
    f_buf: BincodeBuffer,
    g_buf: BincodeBuffer,
    in_var_buf: BincodeBuffer,
    out_var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let f: Option<Expr> =
        from_bincode_buffer(&f_buf);

    let g: Option<Expr> =
        from_bincode_buffer(&g_buf);

    let in_var: Option<String> =
        from_bincode_buffer(
            &in_var_buf,
        );

    let out_var: Option<String> =
        from_bincode_buffer(
            &out_var_buf,
        );

    match (
        f,
        g,
        in_var,
        out_var,
    ) { (
        Some(f),
        Some(g),
        Some(iv),
        Some(ov),
    ) => {

        to_bincode_buffer(&transforms::convolution_laplace(&f, &g, &iv, &ov))
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Computes the partial fraction decomposition of an expression.

///

/// Takes bincode-serialized `Expr` (expression) and `String` (variable).

/// Returns a bincode-serialized `Expr` representing the partial fraction decomposition.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_partial_fraction_decomposition(
    expr_buf: BincodeBuffer,
    var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let var: Option<String> =
        from_bincode_buffer(&var_buf);

    if let (Some(expr), Some(var)) =
        (expr, var)
    {

        if let Some(result) = transforms::partial_fraction_decomposition(&expr, &var) {

            return to_bincode_buffer(&result);
        }
    }

    BincodeBuffer::empty()
}
