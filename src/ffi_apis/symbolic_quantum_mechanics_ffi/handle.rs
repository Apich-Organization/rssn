use crate::symbolic::core::Expr;
use crate::symbolic::quantum_mechanics::{
    self,
    Bra,
    Ket,
    Operator,
};
use std::os::raw::c_char;

/// Creates a new Ket from a state expression.
#[no_mangle]

pub unsafe extern "C" fn rssn_ket_new(
    state: *const Expr
) -> *mut Ket {

    if state.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(Ket {
        state: (*state).clone(),
    }))
}

/// Frees a Ket.
#[no_mangle]

pub unsafe extern "C" fn rssn_ket_free(
    ket: *mut Ket
) {

    if !ket.is_null() {

        drop(Box::from_raw(ket));
    }
}

/// Creates a new Bra from a state expression.
#[no_mangle]

pub unsafe extern "C" fn rssn_bra_new(
    state: *const Expr
) -> *mut Bra {

    if state.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(Bra {
        state: (*state).clone(),
    }))
}

/// Frees a Bra.
#[no_mangle]

pub unsafe extern "C" fn rssn_bra_free(
    bra: *mut Bra
) {

    if !bra.is_null() {

        drop(Box::from_raw(bra));
    }
}

/// Creates a new Operator from an expression.
#[no_mangle]

pub unsafe extern "C" fn rssn_operator_new(
    op: *const Expr
) -> *mut Operator {

    if op.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(Operator {
        op: (*op).clone(),
    }))
}

/// Frees an Operator.
#[no_mangle]

pub unsafe extern "C" fn rssn_operator_free(
    op_ptr: *mut Operator
) {

    if !op_ptr.is_null() {

        drop(Box::from_raw(
            op_ptr,
        ));
    }
}

/// Computes the inner product <Bra|Ket>.
#[no_mangle]

pub unsafe extern "C" fn rssn_bra_ket(
    bra: *const Bra,
    ket: *const Ket,
) -> *mut Expr {

    if bra.is_null() || ket.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        quantum_mechanics::bra_ket(
            &*bra, &*ket,
        ),
    ))
}

/// Computes the commutator [A, B] acting on a Ket.
#[no_mangle]

pub unsafe extern "C" fn rssn_commutator(
    a: *const Operator,
    b: *const Operator,
    ket: *const Ket,
) -> *mut Expr {

    if a.is_null()
        || b.is_null()
        || ket.is_null()
    {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        quantum_mechanics::commutator(
            &*a, &*b, &*ket,
        ),
    ))
}

/// Computes the expectation value <A>.
#[no_mangle]

pub unsafe extern "C" fn rssn_expectation_value(
    op: *const Operator,
    psi: *const Ket,
) -> *mut Expr {

    if op.is_null() || psi.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        quantum_mechanics::expectation_value(&*op, &*psi),
    ))
}

/// Computes the uncertainty ΔA.
#[no_mangle]

pub unsafe extern "C" fn rssn_uncertainty(
    op: *const Operator,
    psi: *const Ket,
) -> *mut Expr {

    if op.is_null() || psi.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        quantum_mechanics::uncertainty(
            &*op, &*psi,
        ),
    ))
}

/// Computes the probability density |ψ(x)|^2.
#[no_mangle]

pub unsafe extern "C" fn rssn_probability_density(
    psi: *const Ket
) -> *mut Expr {

    if psi.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        quantum_mechanics::probability_density(&*psi),
    ))
}

/// Hamiltonian for a free particle.
#[no_mangle]

pub unsafe extern "C" fn rssn_hamiltonian_free_particle(
    m: *const Expr
) -> *mut Operator {

    if m.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        quantum_mechanics::hamiltonian_free_particle(&*m),
    ))
}

/// Hamiltonian for a harmonic oscillator.
#[no_mangle]

pub unsafe extern "C" fn rssn_hamiltonian_harmonic_oscillator(
    m: *const Expr,
    omega: *const Expr,
) -> *mut Operator {

    if m.is_null() || omega.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        quantum_mechanics::hamiltonian_harmonic_oscillator(
            &*m, &*omega,
        ),
    ))
}

/// Pauli matrices σ_x, σ_y, σ_z.
#[no_mangle]

pub unsafe extern "C" fn rssn_pauli_matrices(
    sigma_x: *mut *mut Expr,
    sigma_y: *mut *mut Expr,
    sigma_z: *mut *mut Expr,
) {

    let (sx, sy, sz) = quantum_mechanics::pauli_matrices();

    if !sigma_x.is_null() {

        *sigma_x =
            Box::into_raw(Box::new(sx));
    }

    if !sigma_y.is_null() {

        *sigma_y =
            Box::into_raw(Box::new(sy));
    }

    if !sigma_z.is_null() {

        *sigma_z =
            Box::into_raw(Box::new(sz));
    }
}

/// Spin operator S = hbar/2 * σ.
#[no_mangle]

pub unsafe extern "C" fn rssn_spin_operator(
    pauli: *const Expr
) -> *mut Expr {

    if pauli.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        quantum_mechanics::spin_operator(&*pauli),
    ))
}

/// Time-dependent Schrödinger equation.
#[no_mangle]

pub unsafe extern "C" fn rssn_time_dependent_schrodinger_equation(
    hamiltonian: *const Operator,
    wave_function: *const Ket,
) -> *mut Expr {

    if hamiltonian.is_null()
        || wave_function.is_null()
    {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        quantum_mechanics::time_dependent_schrodinger_equation(
            &*hamiltonian,
            &*wave_function,
        ),
    ))
}

/// First-order energy correction.
#[no_mangle]

pub unsafe extern "C" fn rssn_first_order_energy_correction(
    perturbation: *const Operator,
    unperturbed_state: *const Ket,
) -> *mut Expr {

    if perturbation.is_null()
        || unperturbed_state.is_null()
    {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        quantum_mechanics::first_order_energy_correction(
            &*perturbation,
            &*unperturbed_state,
        ),
    ))
}

/// Dirac equation for a free particle.
#[no_mangle]

pub unsafe extern "C" fn rssn_dirac_equation(
    psi: *const Expr,
    m: *const Expr,
) -> *mut Expr {

    if psi.is_null() || m.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        quantum_mechanics::dirac_equation(&*psi, &*m),
    ))
}

/// Klein-Gordon equation.
#[no_mangle]

pub unsafe extern "C" fn rssn_klein_gordon_equation(
    psi: *const Expr,
    m: *const Expr,
) -> *mut Expr {

    if psi.is_null() || m.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        quantum_mechanics::klein_gordon_equation(
            &*psi, &*m,
        ),
    ))
}

/// Scattering amplitude.
#[no_mangle]

pub unsafe extern "C" fn rssn_scattering_amplitude(
    initial_state: *const Ket,
    final_state: *const Ket,
    potential: *const Operator,
) -> *mut Expr {

    if initial_state.is_null()
        || final_state.is_null()
        || potential.is_null()
    {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        quantum_mechanics::scattering_amplitude(
            &*initial_state,
            &*final_state,
            &*potential,
        ),
    ))
}
