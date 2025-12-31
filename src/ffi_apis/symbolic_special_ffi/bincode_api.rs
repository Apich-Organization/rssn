//! Bincode-based FFI API for numerical special functions.
//!
//! This module provides binary serialization-based FFI functions for various special
//! mathematical functions, offering efficient binary data interchange for high-performance
//! applications.

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::symbolic::special;

// ============================================================================
// Gamma and Related Functions
// ============================================================================

/// Computes Γ(x) via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_gamma_numerical(
    val_buf: BincodeBuffer
) -> BincodeBuffer {

    let val: Option<f64> =
        from_bincode_buffer(&val_buf);

    if let Some(v) = val {

        to_bincode_buffer(
            &special::gamma_numerical(
                v,
            ),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes ln(Γ(x)) via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_ln_gamma_numerical(
    val_buf: BincodeBuffer
) -> BincodeBuffer {

    let val: Option<f64> =
        from_bincode_buffer(&val_buf);

    if let Some(v) = val {

        to_bincode_buffer(&special::ln_gamma_numerical(v))
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes ψ(x) via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_digamma_numerical(
    val_buf: BincodeBuffer
) -> BincodeBuffer {

    let val: Option<f64> =
        from_bincode_buffer(&val_buf);

    if let Some(v) = val {

        to_bincode_buffer(
            &special::digamma_numerical(
                v,
            ),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes B(a, b) via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_beta_numerical(
    a_buf: BincodeBuffer,
    b_buf: BincodeBuffer,
) -> BincodeBuffer {

    let a: Option<f64> =
        from_bincode_buffer(&a_buf);

    let b: Option<f64> =
        from_bincode_buffer(&b_buf);

    if let (Some(val_a), Some(val_b)) =
        (a, b)
    {

        to_bincode_buffer(
            &special::beta_numerical(
                val_a, val_b,
            ),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes ln(B(a, b)) via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_ln_beta_numerical(
    a_buf: BincodeBuffer,
    b_buf: BincodeBuffer,
) -> BincodeBuffer {

    let a: Option<f64> =
        from_bincode_buffer(&a_buf);

    let b: Option<f64> =
        from_bincode_buffer(&b_buf);

    if let (Some(val_a), Some(val_b)) =
        (a, b)
    {

        to_bincode_buffer(
            &special::ln_beta_numerical(
                val_a, val_b,
            ),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes Iₓ(a, b) via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_regularized_incomplete_beta(
    a_buf: BincodeBuffer,
    b_buf: BincodeBuffer,
    x_buf: BincodeBuffer,
) -> BincodeBuffer {

    let a: Option<f64> =
        from_bincode_buffer(&a_buf);

    let b: Option<f64> =
        from_bincode_buffer(&b_buf);

    let x: Option<f64> =
        from_bincode_buffer(&x_buf);

    if let (
        Some(va),
        Some(vb),
        Some(vx),
    ) = (a, b, x)
    {

        to_bincode_buffer(&special::regularized_incomplete_beta(va, vb, vx))
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes P(a, x) via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_regularized_gamma_p(
    a_buf: BincodeBuffer,
    x_buf: BincodeBuffer,
) -> BincodeBuffer {

    let a: Option<f64> =
        from_bincode_buffer(&a_buf);

    let x: Option<f64> =
        from_bincode_buffer(&x_buf);

    if let (Some(va), Some(vx)) = (a, x)
    {

        to_bincode_buffer(&special::regularized_gamma_p(va, vx))
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes Q(a, x) via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_regularized_gamma_q(
    a_buf: BincodeBuffer,
    x_buf: BincodeBuffer,
) -> BincodeBuffer {

    let a: Option<f64> =
        from_bincode_buffer(&a_buf);

    let x: Option<f64> =
        from_bincode_buffer(&x_buf);

    if let (Some(va), Some(vx)) = (a, x)
    {

        to_bincode_buffer(&special::regularized_gamma_q(va, vx))
    } else {

        BincodeBuffer::empty()
    }
}

// ============================================================================
// Error Functions
// ============================================================================

/// Computes erf(x) via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_erf_numerical(
    val_buf: BincodeBuffer
) -> BincodeBuffer {

    let val: Option<f64> =
        from_bincode_buffer(&val_buf);

    if let Some(v) = val {

        to_bincode_buffer(
            &special::erf_numerical(v),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes erfc(x) via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_erfc_numerical(
    val_buf: BincodeBuffer
) -> BincodeBuffer {

    let val: Option<f64> =
        from_bincode_buffer(&val_buf);

    if let Some(v) = val {

        to_bincode_buffer(
            &special::erfc_numerical(v),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes erf⁻¹(x) via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_inverse_erf(
    val_buf: BincodeBuffer
) -> BincodeBuffer {

    let val: Option<f64> =
        from_bincode_buffer(&val_buf);

    if let Some(v) = val {

        to_bincode_buffer(
            &special::inverse_erf(v),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes erfc⁻¹(x) via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_inverse_erfc(
    val_buf: BincodeBuffer
) -> BincodeBuffer {

    let val: Option<f64> =
        from_bincode_buffer(&val_buf);

    if let Some(v) = val {

        to_bincode_buffer(
            &special::inverse_erfc(v),
        )
    } else {

        BincodeBuffer::empty()
    }
}

// ============================================================================
// Combinatorial Functions
// ============================================================================

/// Computes n! via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_factorial(
    n_buf: BincodeBuffer
) -> BincodeBuffer {

    let n: Option<u64> =
        from_bincode_buffer(&n_buf);

    if let Some(v) = n {

        to_bincode_buffer(
            &special::factorial(v),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes n!! via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_double_factorial(
    n_buf: BincodeBuffer
) -> BincodeBuffer {

    let n: Option<u64> =
        from_bincode_buffer(&n_buf);

    if let Some(v) = n {

        to_bincode_buffer(
            &special::double_factorial(
                v,
            ),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes C(n, k) via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_binomial(
    n_buf: BincodeBuffer,
    k_buf: BincodeBuffer,
) -> BincodeBuffer {

    let n: Option<u64> =
        from_bincode_buffer(&n_buf);

    let k: Option<u64> =
        from_bincode_buffer(&k_buf);

    if let (Some(vn), Some(vk)) = (n, k)
    {

        to_bincode_buffer(
            &special::binomial(vn, vk),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes (x)ₙ via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_rising_factorial(
    x_buf: BincodeBuffer,
    n_buf: BincodeBuffer,
) -> BincodeBuffer {

    let x: Option<f64> =
        from_bincode_buffer(&x_buf);

    let n: Option<u32> =
        from_bincode_buffer(&n_buf);

    if let (Some(vx), Some(vn)) = (x, n)
    {

        to_bincode_buffer(
            &special::rising_factorial(
                vx, vn,
            ),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes (x)₍ₙ₎ via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_falling_factorial(
    x_buf: BincodeBuffer,
    n_buf: BincodeBuffer,
) -> BincodeBuffer {

    let x: Option<f64> =
        from_bincode_buffer(&x_buf);

    let n: Option<u32> =
        from_bincode_buffer(&n_buf);

    if let (Some(vx), Some(vn)) = (x, n)
    {

        to_bincode_buffer(
            &special::falling_factorial(
                vx, vn,
            ),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes ln(n!) via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_ln_factorial(
    n_buf: BincodeBuffer
) -> BincodeBuffer {

    let n: Option<u64> =
        from_bincode_buffer(&n_buf);

    if let Some(v) = n {

        to_bincode_buffer(
            &special::ln_factorial(v),
        )
    } else {

        BincodeBuffer::empty()
    }
}

// ============================================================================
// Bessel Functions
// ============================================================================

/// Computes J₀(x) via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_bessel_j0(
    val_buf: BincodeBuffer
) -> BincodeBuffer {

    let val: Option<f64> =
        from_bincode_buffer(&val_buf);

    if let Some(v) = val {

        to_bincode_buffer(
            &special::bessel_j0(v),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes J₁(x) via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_bessel_j1(
    val_buf: BincodeBuffer
) -> BincodeBuffer {

    let val: Option<f64> =
        from_bincode_buffer(&val_buf);

    if let Some(v) = val {

        to_bincode_buffer(
            &special::bessel_j1(v),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes Y₀(x) via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_bessel_y0(
    val_buf: BincodeBuffer
) -> BincodeBuffer {

    let val: Option<f64> =
        from_bincode_buffer(&val_buf);

    if let Some(v) = val {

        to_bincode_buffer(
            &special::bessel_y0(v),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes Y₁(x) via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_bessel_y1(
    val_buf: BincodeBuffer
) -> BincodeBuffer {

    let val: Option<f64> =
        from_bincode_buffer(&val_buf);

    if let Some(v) = val {

        to_bincode_buffer(
            &special::bessel_y1(v),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes I₀(x) via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_bessel_i0(
    val_buf: BincodeBuffer
) -> BincodeBuffer {

    let val: Option<f64> =
        from_bincode_buffer(&val_buf);

    if let Some(v) = val {

        to_bincode_buffer(
            &special::bessel_i0(v),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes I₁(x) via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_bessel_i1(
    val_buf: BincodeBuffer
) -> BincodeBuffer {

    let val: Option<f64> =
        from_bincode_buffer(&val_buf);

    if let Some(v) = val {

        to_bincode_buffer(
            &special::bessel_i1(v),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes K₀(x) via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_bessel_k0(
    val_buf: BincodeBuffer
) -> BincodeBuffer {

    let val: Option<f64> =
        from_bincode_buffer(&val_buf);

    if let Some(v) = val {

        to_bincode_buffer(
            &special::bessel_k0(v),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes K₁(x) via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_bessel_k1(
    val_buf: BincodeBuffer
) -> BincodeBuffer {

    let val: Option<f64> =
        from_bincode_buffer(&val_buf);

    if let Some(v) = val {

        to_bincode_buffer(
            &special::bessel_k1(v),
        )
    } else {

        BincodeBuffer::empty()
    }
}

// ============================================================================
// Other Special Functions
// ============================================================================

/// Computes sinc(x) via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_sinc(
    val_buf: BincodeBuffer
) -> BincodeBuffer {

    let val: Option<f64> =
        from_bincode_buffer(&val_buf);

    if let Some(v) = val {

        to_bincode_buffer(
            &special::sinc(v),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes ζ(s) via Bincode interface.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_zeta_numerical(
    val_buf: BincodeBuffer
) -> BincodeBuffer {

    let val: Option<f64> =
        from_bincode_buffer(&val_buf);

    if let Some(v) = val {

        to_bincode_buffer(
            &special::zeta(v),
        )
    } else {

        BincodeBuffer::empty()
    }
}
