//! Bincode-based FFI API for symbolic integral equations.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::{BincodeBuffer, from_bincode_buffer, to_bincode_buffer};
use crate::symbolic::core::Expr;
use crate::symbolic::integral_equations::solve_airfoil_equation;
use crate::symbolic::integral_equations::FredholmEquation;
use crate::symbolic::integral_equations::VolterraEquation;

#[derive(Serialize, Deserialize)]

struct FredholmNeumannInput {
    equation: FredholmEquation,
    iterations: usize,
}

#[derive(Serialize, Deserialize)]

struct FredholmSeparableInput {
    equation: FredholmEquation,
    a_funcs: Vec<Expr>,
    b_funcs: Vec<Expr>,
}

#[derive(Serialize, Deserialize)]

struct VolterraSuccessiveInput {
    equation: VolterraEquation,
    iterations: usize,
}

#[derive(Serialize, Deserialize)]

struct AirfoilInput {
    f_x: Expr,
    var_x: String,
    var_t: String,
}

/// Solves a Fredholm equation using the Neumann series method (Bincode).
#[unsafe(no_mangle)]

pub extern "C" fn rssn_fredholm_solve_neumann_bincode(
    input_ptr: *const u8,
    input_len: usize,
) -> BincodeBuffer {

    let input_buffer = BincodeBuffer {
        data: input_ptr.cast_mut(),
        len: input_len,
    };

    let input: Option<
        FredholmNeumannInput,
    > = from_bincode_buffer(
        &input_buffer,
    );

    let input = match input {
        | Some(i) => i,
        | None => {
            return BincodeBuffer::empty(
            )
        },
    };

    let result = input
        .equation
        .solve_neumann_series(
            input.iterations,
        );

    to_bincode_buffer(&result)
}

/// Solves a Fredholm equation with a separable kernel (Bincode).
#[unsafe(no_mangle)]

pub extern "C" fn rssn_fredholm_solve_separable_bincode(
    input_ptr: *const u8,
    input_len: usize,
) -> BincodeBuffer {

    let input_buffer = BincodeBuffer {
        data: input_ptr.cast_mut(),
        len: input_len,
    };

    let input: Option<
        FredholmSeparableInput,
    > = from_bincode_buffer(
        &input_buffer,
    );

    let input = match input {
        | Some(i) => i,
        | None => {
            return BincodeBuffer::empty(
            )
        },
    };

    match input
        .equation
        .solve_separable_kernel(
            input.a_funcs,
            input.b_funcs,
        ) {
        | Ok(result) => {
            to_bincode_buffer(&result)
        },
        | Err(_) => {
            BincodeBuffer::empty()
        },
    }
}

/// Solves a Volterra equation using successive approximations (Bincode).
#[unsafe(no_mangle)]

pub extern "C" fn rssn_volterra_solve_successive_bincode(
    input_ptr: *const u8,
    input_len: usize,
) -> BincodeBuffer {

    let input_buffer = BincodeBuffer {
        data: input_ptr.cast_mut(),
        len: input_len,
    };

    let input: Option<
        VolterraSuccessiveInput,
    > = from_bincode_buffer(
        &input_buffer,
    );

    let input = match input {
        | Some(i) => i,
        | None => {
            return BincodeBuffer::empty(
            )
        },
    };

    let result = input
        .equation
        .solve_successive_approximations(input.iterations);

    to_bincode_buffer(&result)
}

/// Solves a Volterra equation by differentiation (Bincode).
#[unsafe(no_mangle)]

pub extern "C" fn rssn_volterra_solve_by_differentiation_bincode(
    input_ptr: *const u8,
    input_len: usize,
) -> BincodeBuffer {

    let input_buffer = BincodeBuffer {
        data: input_ptr.cast_mut(),
        len: input_len,
    };

    let equation: Option<
        VolterraEquation,
    > = from_bincode_buffer(
        &input_buffer,
    );

    let equation = match equation {
        | Some(e) => e,
        | None => {
            return BincodeBuffer::empty(
            )
        },
    };

    match equation
        .solve_by_differentiation()
    {
        | Ok(result) => {
            to_bincode_buffer(&result)
        },
        | Err(_) => {
            BincodeBuffer::empty()
        },
    }
}

/// Solves the airfoil singular integral equation (Bincode).
#[unsafe(no_mangle)]

pub extern "C" fn rssn_solve_airfoil_equation_bincode(
    input_ptr: *const u8,
    input_len: usize,
) -> BincodeBuffer {

    let input_buffer = BincodeBuffer {
        data: input_ptr.cast_mut(),
        len: input_len,
    };

    let input: Option<AirfoilInput> =
        from_bincode_buffer(
            &input_buffer,
        );

    let input = match input {
        | Some(i) => i,
        | None => {
            return BincodeBuffer::empty(
            )
        },
    };

    let result = solve_airfoil_equation(
        &input.f_x,
        &input.var_x,
        &input.var_t,
    );

    to_bincode_buffer(&result)
}
