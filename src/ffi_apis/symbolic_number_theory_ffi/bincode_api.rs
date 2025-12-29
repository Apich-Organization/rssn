use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::symbolic::core::Expr;
use crate::symbolic::number_theory::chinese_remainder;
use crate::symbolic::number_theory::extended_gcd;
use crate::symbolic::number_theory::is_prime;
use crate::symbolic::number_theory::solve_diophantine;

#[derive(serde::Deserialize)]

struct Congruence {
    remainder: Expr,
    modulus: Expr,
}

/// Solves a Diophantine equation.

///

/// Takes bincode-serialized `Expr` representing the equation and `Vec<String>`

/// representing variables. Returns a bincode-serialized `Vec<Expr>` of solutions.

#[no_mangle]

pub extern "C" fn rssn_bincode_solve_diophantine(
    equation_buf: BincodeBuffer,
    vars_buf: BincodeBuffer,
) -> BincodeBuffer {

    let equation: Option<Expr> =
        from_bincode_buffer(
            &equation_buf,
        );

    let vars: Option<Vec<String>> =
        from_bincode_buffer(&vars_buf);

    if let (Some(eq), Some(v)) =
        (equation, vars)
    {

        let v_str: Vec<&str> = v
            .iter()
            .map(std::string::String::as_str)
            .collect();

        match solve_diophantine(
            &eq,
            &v_str,
        ) {
            | Ok(solutions) => {
                to_bincode_buffer(
                    &solutions,
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

/// Computes the extended greatest common divisor (GCD).

///

/// Takes bincode-serialized `Expr` representing two numbers (`a` and `b`),

/// and returns a bincode-serialized tuple `(g, x, y)` where `g` is the GCD

/// and `x`, `y` are coefficients such that `ax + by = g`.

#[no_mangle]

pub extern "C" fn rssn_bincode_extended_gcd(
    a_buf: BincodeBuffer,
    b_buf: BincodeBuffer,
) -> BincodeBuffer {

    let a: Option<Expr> =
        from_bincode_buffer(&a_buf);

    let b: Option<Expr> =
        from_bincode_buffer(&b_buf);

    if let (
        Some(a_expr),
        Some(b_expr),
    ) = (a, b)
    {

        let (g, x, y) = extended_gcd(
            &a_expr,
            &b_expr,
        );

        to_bincode_buffer(&(g, x, y))
    } else {

        BincodeBuffer::empty()
    }
}

/// Checks if a number is prime.

///

/// Takes a bincode-serialized `Expr` representing a number,

/// and returns a bincode-serialized boolean indicating whether the number is prime.

#[no_mangle]

pub extern "C" fn rssn_bincode_is_prime(
    n_buf: BincodeBuffer
) -> BincodeBuffer {

    let n: Option<Expr> =
        from_bincode_buffer(&n_buf);

    if let Some(n_expr) = n {

        let result = is_prime(&n_expr);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Solves the Chinese Remainder Theorem.

///

/// Takes a bincode-serialized vector of congruences (`(remainder, modulus)`),

/// and returns a bincode-serialized `Expr` representing the solution.

#[no_mangle]

pub extern "C" fn rssn_bincode_chinese_remainder(
    congruences_buf: BincodeBuffer
) -> BincodeBuffer {

    let congruences_input: Option<
        Vec<Congruence>,
    > = from_bincode_buffer(
        &congruences_buf,
    );

    if let Some(input) =
        congruences_input
    {

        let congruences: Vec<(
            Expr,
            Expr,
        )> = input
            .into_iter()
            .map(|c| {

                (
                    c.remainder,
                    c.modulus,
                )
            })
            .collect();

        match chinese_remainder(
            &congruences,
        ) {
            | Some(result) => {
                to_bincode_buffer(
                    &result,
                )
            },
            | None => {
                BincodeBuffer::empty()
            },
        }
    } else {

        BincodeBuffer::empty()
    }
}
