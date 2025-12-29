//! JSON-based FFI API for cryptographic operations.
//!
//! This module provides JSON string-based FFI functions for elliptic curve cryptography.

use std::os::raw::c_char;
use std::str::FromStr;

use num_bigint::BigInt;

use crate::ffi_apis::common::*;
use crate::symbolic::cryptography::ecdsa_sign;
use crate::symbolic::cryptography::ecdsa_verify;
use crate::symbolic::cryptography::generate_keypair;
use crate::symbolic::cryptography::generate_shared_secret;
use crate::symbolic::cryptography::point_compress;
use crate::symbolic::cryptography::point_decompress;
use crate::symbolic::cryptography::CurvePoint;
use crate::symbolic::cryptography::EcdsaSignature;
use crate::symbolic::cryptography::EllipticCurve;
use crate::symbolic::finite_field::PrimeField;
use crate::symbolic::finite_field::PrimeFieldElement;

/// Helper to parse BigInt from string or JSON string.

fn parse_bigint(
    s: Option<String>
) -> Option<BigInt> {

    s.and_then(|str_val| {

        BigInt::from_str(&str_val).ok()
    })
}

/// Creates a new elliptic curve.
/// Arguments: a (str), b (str), modulus (str)
#[no_mangle]

pub unsafe extern "C" fn rssn_json_elliptic_curve_new(
    a_json: *const c_char,
    b_json: *const c_char,
    modulus_json: *const c_char,
) -> *mut c_char {

    let a_str: Option<String> =
        from_json_string(a_json);

    let b_str: Option<String> =
        from_json_string(b_json);

    let mod_str: Option<String> =
        from_json_string(modulus_json);

    if let (Some(a), Some(b), Some(m)) = (
        parse_bigint(a_str),
        parse_bigint(b_str),
        parse_bigint(mod_str),
    ) {

        to_json_string(
            &EllipticCurve::new(
                a, b, m,
            ),
        )
    } else {

        std::ptr::null_mut()
    }
}

/// Creates an affine curve point.
/// Arguments: x (str), y (str), modulus (str)
#[no_mangle]

pub unsafe extern "C" fn rssn_json_curve_point_affine(
    x_json: *const c_char,
    y_json: *const c_char,
    modulus_json: *const c_char,
) -> *mut c_char {

    let x_str: Option<String> =
        from_json_string(x_json);

    let y_str: Option<String> =
        from_json_string(y_json);

    let mod_str: Option<String> =
        from_json_string(modulus_json);

    if let (Some(x), Some(y), Some(m)) = (
        parse_bigint(x_str),
        parse_bigint(y_str),
        parse_bigint(mod_str),
    ) {

        let field = PrimeField::new(m);

        let point = CurvePoint::Affine {
            x : PrimeFieldElement::new(x, field.clone()),
            y : PrimeFieldElement::new(y, field),
        };

        to_json_string(&point)
    } else {

        std::ptr::null_mut()
    }
}

/// Creates a point at infinity.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_curve_point_infinity(
) -> *mut c_char {

    to_json_string(
        &CurvePoint::Infinity,
    )
}

/// Checks if a point is on the curve.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_curve_is_on_curve(
    curve_json: *const c_char,
    point_json: *const c_char,
) -> *mut c_char {

    let curve: Option<EllipticCurve> =
        from_json_string(curve_json);

    let point: Option<CurvePoint> =
        from_json_string(point_json);

    if let (Some(c), Some(p)) =
        (curve, point)
    {

        to_json_string(
            &c.is_on_curve(&p),
        )
    } else {

        std::ptr::null_mut()
    }
}

/// Negates a point.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_curve_negate(
    curve_json: *const c_char,
    point_json: *const c_char,
) -> *mut c_char {

    let curve: Option<EllipticCurve> =
        from_json_string(curve_json);

    let point: Option<CurvePoint> =
        from_json_string(point_json);

    if let (Some(c), Some(p)) =
        (curve, point)
    {

        to_json_string(&c.negate(&p))
    } else {

        std::ptr::null_mut()
    }
}

/// Doubles a point.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_curve_double(
    curve_json: *const c_char,
    point_json: *const c_char,
) -> *mut c_char {

    let curve: Option<EllipticCurve> =
        from_json_string(curve_json);

    let point: Option<CurvePoint> =
        from_json_string(point_json);

    if let (Some(c), Some(p)) =
        (curve, point)
    {

        to_json_string(&c.double(&p))
    } else {

        std::ptr::null_mut()
    }
}

/// Adds two points.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_curve_add(
    curve_json: *const c_char,
    p1_json: *const c_char,
    p2_json: *const c_char,
) -> *mut c_char {

    let curve: Option<EllipticCurve> =
        from_json_string(curve_json);

    let p1: Option<CurvePoint> =
        from_json_string(p1_json);

    let p2: Option<CurvePoint> =
        from_json_string(p2_json);

    if let (
        Some(c),
        Some(p1),
        Some(p2),
    ) = (curve, p1, p2)
    {

        to_json_string(&c.add(&p1, &p2))
    } else {

        std::ptr::null_mut()
    }
}

/// Scalar multiplication.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_curve_scalar_mult(
    curve_json: *const c_char,
    k_json: *const c_char,
    p_json: *const c_char,
) -> *mut c_char {

    let curve: Option<EllipticCurve> =
        from_json_string(curve_json);

    let k_str: Option<String> =
        from_json_string(k_json);

    let p: Option<CurvePoint> =
        from_json_string(p_json);

    if let (Some(c), Some(k), Some(p)) = (
        curve,
        parse_bigint(k_str),
        p,
    ) {

        to_json_string(
            &c.scalar_mult(&k, &p),
        )
    } else {

        std::ptr::null_mut()
    }
}

/// Generates a key pair.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_generate_keypair(
    curve_json: *const c_char,
    generator_json: *const c_char,
) -> *mut c_char {

    let curve: Option<EllipticCurve> =
        from_json_string(curve_json);

    let gen: Option<CurvePoint> =
        from_json_string(
            generator_json,
        );

    if let (Some(c), Some(g)) =
        (curve, gen)
    {

        to_json_string(
            &generate_keypair(&c, &g),
        )
    } else {

        std::ptr::null_mut()
    }
}

/// Generates a shared secret.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_generate_shared_secret(
    curve_json: *const c_char,
    private_key_json: *const c_char,
    other_public_key_json : *const c_char,
) -> *mut c_char {

    let curve: Option<EllipticCurve> =
        from_json_string(curve_json);

    let pk_str: Option<String> =
        from_json_string(
            private_key_json,
        );

    let other_pub: Option<CurvePoint> =
        from_json_string(
            other_public_key_json,
        );

    if let (
        Some(c),
        Some(pk),
        Some(opub),
    ) = (
        curve,
        parse_bigint(pk_str),
        other_pub,
    ) {

        to_json_string(
            &generate_shared_secret(
                &c, &pk, &opub,
            ),
        )
    } else {

        std::ptr::null_mut()
    }
}

/// Signs a message.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_ecdsa_sign(
    message_hash_json: *const c_char,
    private_key_json: *const c_char,
    curve_json: *const c_char,
    generator_json: *const c_char,
    order_json: *const c_char,
) -> *mut c_char {

    let hash_str: Option<String> =
        from_json_string(
            message_hash_json,
        );

    let pk_str: Option<String> =
        from_json_string(
            private_key_json,
        );

    let curve: Option<EllipticCurve> =
        from_json_string(curve_json);

    let gen: Option<CurvePoint> =
        from_json_string(
            generator_json,
        );

    let order_str: Option<String> =
        from_json_string(order_json);

    if let (
        Some(h),
        Some(pk),
        Some(c),
        Some(g),
        Some(o),
    ) = (
        parse_bigint(hash_str),
        parse_bigint(pk_str),
        curve,
        gen,
        parse_bigint(order_str),
    ) {

        if let Some(sig) = ecdsa_sign(
            &h, &pk, &c, &g, &o,
        ) {

            to_json_string(&sig)
        } else {

            std::ptr::null_mut()
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Verifies a signature.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_ecdsa_verify(
    message_hash_json: *const c_char,
    signature_json: *const c_char,
    public_key_json: *const c_char,
    curve_json: *const c_char,
    generator_json: *const c_char,
    order_json: *const c_char,
) -> *mut c_char {

    let hash_str: Option<String> =
        from_json_string(
            message_hash_json,
        );

    let sig: Option<EcdsaSignature> =
        from_json_string(
            signature_json,
        );

    let pub_key: Option<CurvePoint> =
        from_json_string(
            public_key_json,
        );

    let curve: Option<EllipticCurve> =
        from_json_string(curve_json);

    let gen: Option<CurvePoint> =
        from_json_string(
            generator_json,
        );

    let order_str: Option<String> =
        from_json_string(order_json);

    if let (
        Some(h),
        Some(sig),
        Some(pk),
        Some(c),
        Some(g),
        Some(o),
    ) = (
        parse_bigint(hash_str),
        sig,
        pub_key,
        curve,
        gen,
        parse_bigint(order_str),
    ) {

        to_json_string(&ecdsa_verify(
            &h, &sig, &pk, &c, &g, &o,
        ))
    } else {

        std::ptr::null_mut()
    }
}

/// Compresses a point.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_point_compress(
    point_json: *const c_char
) -> *mut c_char {

    let point: Option<CurvePoint> =
        from_json_string(point_json);

    if let Some(p) = point {

        if let Some((x, is_odd)) =
            point_compress(&p)
        {

            let obj = serde_json::json!({
                "x": x.to_string(),
                "is_odd": is_odd
            });

            to_json_string(&obj)
        } else {

            std::ptr::null_mut()
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Decompresses a point.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_point_decompress(
    x_json: *const c_char,
    is_odd_json: *const c_char,
    curve_json: *const c_char,
) -> *mut c_char {

    let x_str: Option<String> =
        from_json_string(x_json);

    let is_odd: Option<bool> =
        from_json_string(is_odd_json);

    let curve: Option<EllipticCurve> =
        from_json_string(curve_json);

    if let (
        Some(x),
        Some(io),
        Some(c),
    ) = (
        parse_bigint(x_str),
        is_odd,
        curve,
    ) {

        if let Some(p) =
            point_decompress(x, io, &c)
        {

            to_json_string(&p)
        } else {

            std::ptr::null_mut()
        }
    } else {

        std::ptr::null_mut()
    }
}
