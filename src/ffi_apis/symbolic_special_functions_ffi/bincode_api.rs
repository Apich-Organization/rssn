use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::symbolic::core::Expr;
use crate::symbolic::special_functions;

/// Computes the gamma function Γ(z).

///

/// Takes a bincode-serialized `Expr` representing `z` as input,

/// and returns a bincode-serialized `Expr` representing Γ(z).

#[no_mangle]

pub extern "C" fn rssn_bincode_gamma(
    arg_buf: BincodeBuffer
) -> BincodeBuffer {

    let arg: Option<Expr> =
        from_bincode_buffer(&arg_buf);

    if let Some(a) = arg {

        to_bincode_buffer(
            &special_functions::gamma(
                a,
            ),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the log-gamma function ln(Γ(z)).

///

/// Takes a bincode-serialized `Expr` representing `z` as input,

/// and returns a bincode-serialized `Expr` representing ln(Γ(z)).

#[no_mangle]

pub extern "C" fn rssn_bincode_ln_gamma(
    arg_buf: BincodeBuffer
) -> BincodeBuffer {

    let arg: Option<Expr> =
        from_bincode_buffer(&arg_buf);

    if let Some(a) = arg {

        to_bincode_buffer(&special_functions::ln_gamma(a))
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the beta function B(a, b).

///

/// Takes bincode-serialized `Expr` representing `a` and `b` as inputs,

/// and returns a bincode-serialized `Expr` representing B(a, b).

#[no_mangle]

pub extern "C" fn rssn_bincode_beta(
    a_buf: BincodeBuffer,
    b_buf: BincodeBuffer,
) -> BincodeBuffer {

    let a: Option<Expr> =
        from_bincode_buffer(&a_buf);

    let b: Option<Expr> =
        from_bincode_buffer(&b_buf);

    if let (Some(val_a), Some(val_b)) =
        (a, b)
    {

        to_bincode_buffer(
            &special_functions::beta(
                val_a, val_b,
            ),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the error function erf(z).

///

/// Takes a bincode-serialized `Expr` representing `z` as input,

/// and returns a bincode-serialized `Expr` representing erf(z).

#[no_mangle]

pub extern "C" fn rssn_bincode_erf(
    arg_buf: BincodeBuffer
) -> BincodeBuffer {

    let arg: Option<Expr> =
        from_bincode_buffer(&arg_buf);

    if let Some(a) = arg {

        to_bincode_buffer(
            &special_functions::erf(a),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the complementary error function erfc(z).

///

/// Takes a bincode-serialized `Expr` representing `z` as input,

/// and returns a bincode-serialized `Expr` representing erfc(z).

#[no_mangle]

pub extern "C" fn rssn_bincode_erfc(
    arg_buf: BincodeBuffer
) -> BincodeBuffer {

    let arg: Option<Expr> =
        from_bincode_buffer(&arg_buf);

    if let Some(a) = arg {

        to_bincode_buffer(
            &special_functions::erfc(a),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the imaginary error function erfi(z).

///

/// Takes a bincode-serialized `Expr` representing `z` as input,

/// and returns a bincode-serialized `Expr` representing erfi(z).

#[no_mangle]

pub extern "C" fn rssn_bincode_erfi(
    arg_buf: BincodeBuffer
) -> BincodeBuffer {

    let arg: Option<Expr> =
        from_bincode_buffer(&arg_buf);

    if let Some(a) = arg {

        to_bincode_buffer(
            &special_functions::erfi(a),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the Riemann zeta function ζ(s).

///

/// Takes a bincode-serialized `Expr` representing `s` as input,

/// and returns a bincode-serialized `Expr` representing ζ(s).

#[no_mangle]

pub extern "C" fn rssn_bincode_zeta(
    arg_buf: BincodeBuffer
) -> BincodeBuffer {

    let arg: Option<Expr> =
        from_bincode_buffer(&arg_buf);

    if let Some(a) = arg {

        to_bincode_buffer(
            &special_functions::zeta(a),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the digamma function ψ(z).

///

/// Takes a bincode-serialized `Expr` representing `z` as input,

/// and returns a bincode-serialized `Expr` representing ψ(z).

#[no_mangle]

pub extern "C" fn rssn_bincode_digamma(
    arg_buf: BincodeBuffer
) -> BincodeBuffer {

    let arg: Option<Expr> =
        from_bincode_buffer(&arg_buf);

    if let Some(a) = arg {

        to_bincode_buffer(
            &special_functions::digamma(
                a,
            ),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the polygamma function ψ⁽ⁿ⁾(z).

///

/// Takes bincode-serialized `Expr` representing `n` and `z` as inputs,

/// and returns a bincode-serialized `Expr` representing ψ⁽ⁿ⁾(z).

#[no_mangle]

pub extern "C" fn rssn_bincode_polygamma(
    n_buf: BincodeBuffer,
    z_buf: BincodeBuffer,
) -> BincodeBuffer {

    let n: Option<Expr> =
        from_bincode_buffer(&n_buf);

    let z: Option<Expr> =
        from_bincode_buffer(&z_buf);

    if let (Some(n), Some(z)) = (n, z) {

        to_bincode_buffer(&special_functions::polygamma(n, z))
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the Bessel function of the first kind `J_n(x)`.

///

/// Takes bincode-serialized `Expr` representing `n` and `x` as inputs,

/// and returns a bincode-serialized `Expr` representing `J_n(x)`.

#[no_mangle]

pub extern "C" fn rssn_bincode_bessel_j(
    order_buf: BincodeBuffer,
    arg_buf: BincodeBuffer,
) -> BincodeBuffer {

    let order: Option<Expr> =
        from_bincode_buffer(&order_buf);

    let arg: Option<Expr> =
        from_bincode_buffer(&arg_buf);

    if let (Some(o), Some(a)) =
        (order, arg)
    {

        to_bincode_buffer(&special_functions::bessel_j(o, a))
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the Bessel function of the second kind `Y_n(x)`.

///

/// Takes bincode-serialized `Expr` representing `n` and `x` as inputs,

/// and returns a bincode-serialized `Expr` representing `Y_n(x)`.

#[no_mangle]

pub extern "C" fn rssn_bincode_bessel_y(
    order_buf: BincodeBuffer,
    arg_buf: BincodeBuffer,
) -> BincodeBuffer {

    let order: Option<Expr> =
        from_bincode_buffer(&order_buf);

    let arg: Option<Expr> =
        from_bincode_buffer(&arg_buf);

    if let (Some(o), Some(a)) =
        (order, arg)
    {

        to_bincode_buffer(&special_functions::bessel_y(o, a))
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the modified Bessel function of the first kind `I_n(x)`.

///

/// Takes bincode-serialized `Expr` representing `n` and `x` as inputs,

/// and returns a bincode-serialized `Expr` representing `I_n(x)`.

#[no_mangle]

pub extern "C" fn rssn_bincode_bessel_i(
    order_buf: BincodeBuffer,
    arg_buf: BincodeBuffer,
) -> BincodeBuffer {

    let order: Option<Expr> =
        from_bincode_buffer(&order_buf);

    let arg: Option<Expr> =
        from_bincode_buffer(&arg_buf);

    if let (Some(o), Some(a)) =
        (order, arg)
    {

        to_bincode_buffer(&special_functions::bessel_i(o, a))
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the modified Bessel function of the second kind `K_n(x)`.

///

/// Takes bincode-serialized `Expr` representing `n` and `x` as inputs,

/// and returns a bincode-serialized `Expr` representing `K_n(x)`.

#[no_mangle]

pub extern "C" fn rssn_bincode_bessel_k(
    order_buf: BincodeBuffer,
    arg_buf: BincodeBuffer,
) -> BincodeBuffer {

    let order: Option<Expr> =
        from_bincode_buffer(&order_buf);

    let arg: Option<Expr> =
        from_bincode_buffer(&arg_buf);

    if let (Some(o), Some(a)) =
        (order, arg)
    {

        to_bincode_buffer(&special_functions::bessel_k(o, a))
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the Legendre polynomial `P_n(x)`.

///

/// Takes bincode-serialized `Expr` representing `n` and `x` as inputs,

/// and returns a bincode-serialized `Expr` representing `P_n(x)`.

#[no_mangle]

pub extern "C" fn rssn_bincode_legendre_p(
    degree_buf: BincodeBuffer,
    arg_buf: BincodeBuffer,
) -> BincodeBuffer {

    let degree: Option<Expr> =
        from_bincode_buffer(
            &degree_buf,
        );

    let arg: Option<Expr> =
        from_bincode_buffer(&arg_buf);

    if let (Some(d), Some(a)) =
        (degree, arg)
    {

        to_bincode_buffer(&special_functions::legendre_p(d, a))
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the Laguerre polynomial `L_n(x)`.

///

/// Takes bincode-serialized `Expr` representing `n` and `x` as inputs,

/// and returns a bincode-serialized `Expr` representing `L_n(x)`.

#[no_mangle]

pub extern "C" fn rssn_bincode_laguerre_l(
    degree_buf: BincodeBuffer,
    arg_buf: BincodeBuffer,
) -> BincodeBuffer {

    let degree: Option<Expr> =
        from_bincode_buffer(
            &degree_buf,
        );

    let arg: Option<Expr> =
        from_bincode_buffer(&arg_buf);

    if let (Some(d), Some(a)) =
        (degree, arg)
    {

        to_bincode_buffer(&special_functions::laguerre_l(d, a))
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the generalized Laguerre polynomial `L_n^α(x)`.

///

/// Takes bincode-serialized `Expr` representing `n`, `alpha`, and `x` as inputs,

/// and returns a bincode-serialized `Expr` representing `L_n^α(x)`.

#[no_mangle]

pub extern "C" fn rssn_bincode_generalized_laguerre(
    n_buf: BincodeBuffer,
    alpha_buf: BincodeBuffer,
    x_buf: BincodeBuffer,
) -> BincodeBuffer {

    let n: Option<Expr> =
        from_bincode_buffer(&n_buf);

    let alpha: Option<Expr> =
        from_bincode_buffer(&alpha_buf);

    let x: Option<Expr> =
        from_bincode_buffer(&x_buf);

    if let (
        Some(n),
        Some(alpha),
        Some(x),
    ) = (n, alpha, x)
    {

        to_bincode_buffer(&special_functions::generalized_laguerre(n, alpha, x))
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the Hermite polynomial `H_n(x)`.

///

/// Takes bincode-serialized `Expr` representing `n` and `x` as inputs,

/// and returns a bincode-serialized `Expr` representing `H_n(x)`.

#[no_mangle]

pub extern "C" fn rssn_bincode_hermite_h(
    degree_buf: BincodeBuffer,
    arg_buf: BincodeBuffer,
) -> BincodeBuffer {

    let degree: Option<Expr> =
        from_bincode_buffer(
            &degree_buf,
        );

    let arg: Option<Expr> =
        from_bincode_buffer(&arg_buf);

    if let (Some(d), Some(a)) =
        (degree, arg)
    {

        to_bincode_buffer(&special_functions::hermite_h(d, a))
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the Chebyshev polynomial of the first kind `T_n(x)`.

///

/// Takes bincode-serialized `Expr` representing `n` and `x` as inputs,

/// and returns a bincode-serialized `Expr` representing `T_n(x)`.

#[no_mangle]

pub extern "C" fn rssn_bincode_chebyshev_t(
    n_buf: BincodeBuffer,
    x_buf: BincodeBuffer,
) -> BincodeBuffer {

    let n: Option<Expr> =
        from_bincode_buffer(&n_buf);

    let x: Option<Expr> =
        from_bincode_buffer(&x_buf);

    if let (Some(n), Some(x)) = (n, x) {

        to_bincode_buffer(&special_functions::chebyshev_t(n, x))
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the Chebyshev polynomial of the second kind `U_n(x)`.

///

/// Takes bincode-serialized `Expr` representing `n` and `x` as inputs,

/// and returns a bincode-serialized `Expr` representing `U_n(x)`.

#[no_mangle]

pub extern "C" fn rssn_bincode_chebyshev_u(
    n_buf: BincodeBuffer,
    x_buf: BincodeBuffer,
) -> BincodeBuffer {

    let n: Option<Expr> =
        from_bincode_buffer(&n_buf);

    let x: Option<Expr> =
        from_bincode_buffer(&x_buf);

    if let (Some(n), Some(x)) = (n, x) {

        to_bincode_buffer(&special_functions::chebyshev_u(n, x))
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the Bessel differential equation.

///

/// Takes bincode-serialized `Expr` representing `y`, `x`, and `n` as inputs,

/// and returns a bincode-serialized `Expr` representing the equation.

#[no_mangle]

pub extern "C" fn rssn_bincode_bessel_differential_equation(
    y_buf: BincodeBuffer,
    x_buf: BincodeBuffer,
    n_buf: BincodeBuffer,
) -> BincodeBuffer {

    let y: Option<Expr> =
        from_bincode_buffer(&y_buf);

    let x: Option<Expr> =
        from_bincode_buffer(&x_buf);

    let n: Option<Expr> =
        from_bincode_buffer(&n_buf);

    if let (Some(y), Some(x), Some(n)) =
        (y, x, n)
    {

        to_bincode_buffer(&special_functions::bessel_differential_equation(&y, &x, &n))
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the Legendre differential equation.

///

/// Takes bincode-serialized `Expr` representing `y`, `x`, and `n` as inputs,

/// and returns a bincode-serialized `Expr` representing the equation.

#[no_mangle]

pub extern "C" fn rssn_bincode_legendre_differential_equation(
    y_buf: BincodeBuffer,
    x_buf: BincodeBuffer,
    n_buf: BincodeBuffer,
) -> BincodeBuffer {

    let y: Option<Expr> =
        from_bincode_buffer(&y_buf);

    let x: Option<Expr> =
        from_bincode_buffer(&x_buf);

    let n: Option<Expr> =
        from_bincode_buffer(&n_buf);

    if let (Some(y), Some(x), Some(n)) =
        (y, x, n)
    {

        to_bincode_buffer(&special_functions::legendre_differential_equation(&y, &x, &n))
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the Laguerre differential equation.

///

/// Takes bincode-serialized `Expr` representing `y`, `x`, and `n` as inputs,

/// and returns a bincode-serialized `Expr` representing the equation.

#[no_mangle]

pub extern "C" fn rssn_bincode_laguerre_differential_equation(
    y_buf: BincodeBuffer,
    x_buf: BincodeBuffer,
    n_buf: BincodeBuffer,
) -> BincodeBuffer {

    let y: Option<Expr> =
        from_bincode_buffer(&y_buf);

    let x: Option<Expr> =
        from_bincode_buffer(&x_buf);

    let n: Option<Expr> =
        from_bincode_buffer(&n_buf);

    if let (Some(y), Some(x), Some(n)) =
        (y, x, n)
    {

        to_bincode_buffer(&special_functions::laguerre_differential_equation(&y, &x, &n))
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the Hermite differential equation.

///

/// Takes bincode-serialized `Expr` representing `y`, `x`, and `n` as inputs,

/// and returns a bincode-serialized `Expr` representing the equation.

#[no_mangle]

pub extern "C" fn rssn_bincode_hermite_differential_equation(
    y_buf: BincodeBuffer,
    x_buf: BincodeBuffer,
    n_buf: BincodeBuffer,
) -> BincodeBuffer {

    let y: Option<Expr> =
        from_bincode_buffer(&y_buf);

    let x: Option<Expr> =
        from_bincode_buffer(&x_buf);

    let n: Option<Expr> =
        from_bincode_buffer(&n_buf);

    if let (Some(y), Some(x), Some(n)) =
        (y, x, n)
    {

        to_bincode_buffer(&special_functions::hermite_differential_equation(&y, &x, &n))
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the Chebyshev differential equation.

///

/// Takes bincode-serialized `Expr` representing `y`, `x`, and `n` as inputs,

/// and returns a bincode-serialized `Expr` representing the equation.

#[no_mangle]

pub extern "C" fn rssn_bincode_chebyshev_differential_equation(
    y_buf: BincodeBuffer,
    x_buf: BincodeBuffer,
    n_buf: BincodeBuffer,
) -> BincodeBuffer {

    let y: Option<Expr> =
        from_bincode_buffer(&y_buf);

    let x: Option<Expr> =
        from_bincode_buffer(&x_buf);

    let n: Option<Expr> =
        from_bincode_buffer(&n_buf);

    if let (Some(y), Some(x), Some(n)) =
        (y, x, n)
    {

        to_bincode_buffer(&special_functions::chebyshev_differential_equation(&y, &x, &n))
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the Rodrigues' formula for Legendre polynomials.

///

/// Takes bincode-serialized `Expr` representing `n` and `x` as inputs,

/// and returns a bincode-serialized `Expr` representing the formula.

#[no_mangle]

pub extern "C" fn rssn_bincode_legendre_rodrigues_formula(
    n_buf: BincodeBuffer,
    x_buf: BincodeBuffer,
) -> BincodeBuffer {

    let n: Option<Expr> =
        from_bincode_buffer(&n_buf);

    let x: Option<Expr> =
        from_bincode_buffer(&x_buf);

    if let (Some(n), Some(x)) = (n, x) {

        to_bincode_buffer(&special_functions::legendre_rodrigues_formula(&n, &x))
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the Rodrigues' formula for Hermite polynomials.

///

/// Takes bincode-serialized `Expr` representing `n` and `x` as inputs,

/// and returns a bincode-serialized `Expr` representing the formula.

#[no_mangle]

pub extern "C" fn rssn_bincode_hermite_rodrigues_formula(
    n_buf: BincodeBuffer,
    x_buf: BincodeBuffer,
) -> BincodeBuffer {

    let n: Option<Expr> =
        from_bincode_buffer(&n_buf);

    let x: Option<Expr> =
        from_bincode_buffer(&x_buf);

    if let (Some(n), Some(x)) = (n, x) {

        to_bincode_buffer(&special_functions::hermite_rodrigues_formula(&n, &x))
    } else {

        BincodeBuffer::empty()
    }
}
