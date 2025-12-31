use crate::symbolic::core::Expr;
use crate::symbolic::quantum_mechanics::Bra;
use crate::symbolic::quantum_mechanics::Ket;
use crate::symbolic::quantum_mechanics::Operator;
use crate::symbolic::quantum_mechanics::{
    self,
};

/// Creates a new Ket from a state expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_ket_new(
    state: *const Expr
) -> *mut Ket {

    unsafe {

        if state.is_null() {

            return std::ptr::null_mut(
            );
        }

        Box::into_raw(Box::new(Ket {
            state: (*state).clone(),
        }))
    }
}

/// Frees a Ket.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_ket_free(
    ket: *mut Ket
) {

    unsafe {

        if !ket.is_null() {

            drop(Box::from_raw(ket));
        }
    }
}

/// Creates a new Bra from a state expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bra_new(
    state: *const Expr
) -> *mut Bra {

    unsafe {

        if state.is_null() {

            return std::ptr::null_mut(
            );
        }

        Box::into_raw(Box::new(Bra {
            state: (*state).clone(),
        }))
    }
}

/// Frees a Bra.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bra_free(
    bra: *mut Bra
) {

    unsafe {

        if !bra.is_null() {

            drop(Box::from_raw(bra));
        }
    }
}

/// Creates a new Operator from an expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_operator_new(
    op: *const Expr
) -> *mut Operator {

    unsafe {

        if op.is_null() {

            return std::ptr::null_mut(
            );
        }

        Box::into_raw(Box::new(
            Operator {
                op: (*op).clone(),
            },
        ))
    }
}

/// Frees an Operator.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_operator_free(
    op_ptr: *mut Operator
) {

    unsafe {

        if !op_ptr.is_null() {

            drop(Box::from_raw(
                op_ptr,
            ));
        }
    }
}

/// Computes the inner product <Bra|Ket>.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bra_ket(
    bra: *const Bra,
    ket: *const Ket,
) -> *mut Expr {

    unsafe {

        if bra.is_null()
            || ket.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        Box::into_raw(Box::new(
            quantum_mechanics::bra_ket(
                &*bra, &*ket,
            ),
        ))
    }
}

/// Computes the commutator [A, B] acting on a Ket.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_commutator(
    a: *const Operator,
    b: *const Operator,
    ket: *const Ket,
) -> *mut Expr {

    unsafe {

        if a.is_null()
            || b.is_null()
            || ket.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        Box::into_raw(Box::new(
        quantum_mechanics::commutator(
            &*a, &*b, &*ket,
        ),
    ))
    }
}

/// Computes the expectation value `<A>`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_expectation_value(
    op: *const Operator,
    psi: *const Ket,
) -> *mut Expr {

    unsafe {

        if op.is_null() || psi.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        Box::into_raw(Box::new(
        quantum_mechanics::expectation_value(&*op, &*psi),
    ))
    }
}

/// Computes the uncertainty ΔA.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_uncertainty(
    op: *const Operator,
    psi: *const Ket,
) -> *mut Expr {

    unsafe {

        if op.is_null() || psi.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        Box::into_raw(Box::new(
        quantum_mechanics::uncertainty(
            &*op, &*psi,
        ),
    ))
    }
}

/// Computes the probability density |ψ(x)|^2.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_probability_density(
    psi: *const Ket
) -> *mut Expr {

    unsafe {

        if psi.is_null() {

            return std::ptr::null_mut(
            );
        }

        Box::into_raw(Box::new(
        quantum_mechanics::probability_density(&*psi),
    ))
    }
}

/// Hamiltonian for a free particle.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_hamiltonian_free_particle(
    m: *const Expr
) -> *mut Operator {

    unsafe {

        if m.is_null() {

            return std::ptr::null_mut(
            );
        }

        Box::into_raw(Box::new(
        quantum_mechanics::hamiltonian_free_particle(&*m),
    ))
    }
}

/// Hamiltonian for a harmonic oscillator.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_hamiltonian_harmonic_oscillator(
    m: *const Expr,
    omega: *const Expr,
) -> *mut Operator {

    unsafe {

        if m.is_null()
            || omega.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        Box::into_raw(Box::new(
        quantum_mechanics::hamiltonian_harmonic_oscillator(&*m, &*omega),
    ))
    }
}

/// Pauli matrices `σ_x`, `σ_y`, `σ_z`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_pauli_matrices(
    sigma_x: *mut *mut Expr,
    sigma_y: *mut *mut Expr,
    sigma_z: *mut *mut Expr,
) {

    unsafe {

        let (sx, sy, sz) = quantum_mechanics::pauli_matrices();

        if !sigma_x.is_null() {

            *sigma_x = Box::into_raw(
                Box::new(sx),
            );
        }

        if !sigma_y.is_null() {

            *sigma_y = Box::into_raw(
                Box::new(sy),
            );
        }

        if !sigma_z.is_null() {

            *sigma_z = Box::into_raw(
                Box::new(sz),
            );
        }
    }
}

/// Spin operator S = hbar/2 * σ.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_spin_operator(
    pauli: *const Expr
) -> *mut Expr {

    unsafe {

        if pauli.is_null() {

            return std::ptr::null_mut(
            );
        }

        Box::into_raw(Box::new(
        quantum_mechanics::spin_operator(&*pauli),
    ))
    }
}

/// Time-dependent Schrödinger equation.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_time_dependent_schrodinger_equation(
    hamiltonian: *const Operator,
    wave_function: *const Ket,
) -> *mut Expr {

    unsafe {

        if hamiltonian.is_null()
            || wave_function.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        Box::into_raw(Box::new(
        quantum_mechanics::time_dependent_schrodinger_equation(
            &*hamiltonian,
            &*wave_function,
        ),
    ))
    }
}

/// First-order energy correction.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_first_order_energy_correction(
    perturbation: *const Operator,
    unperturbed_state: *const Ket,
) -> *mut Expr {

    unsafe {

        if perturbation.is_null()
            || unperturbed_state
                .is_null()
        {

            return std::ptr::null_mut(
            );
        }

        Box::into_raw(Box::new(
        quantum_mechanics::first_order_energy_correction(
            &*perturbation,
            &*unperturbed_state,
        ),
    ))
    }
}

/// Dirac equation for a free particle.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_dirac_equation(
    psi: *const Expr,
    m: *const Expr,
) -> *mut Expr {

    unsafe {

        if psi.is_null() || m.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        Box::into_raw(Box::new(
        quantum_mechanics::dirac_equation(&*psi, &*m),
    ))
    }
}

/// Klein-Gordon equation.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_klein_gordon_equation(
    psi: *const Expr,
    m: *const Expr,
) -> *mut Expr {

    unsafe {

        if psi.is_null() || m.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        Box::into_raw(Box::new(
        quantum_mechanics::klein_gordon_equation(&*psi, &*m),
    ))
    }
}

/// Scattering amplitude.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_scattering_amplitude(
    initial_state: *const Ket,
    final_state: *const Ket,
    potential: *const Operator,
) -> *mut Expr {

    unsafe {

        if initial_state.is_null()
            || final_state.is_null()
            || potential.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        Box::into_raw(Box::new(
        quantum_mechanics::scattering_amplitude(
            &*initial_state,
            &*final_state,
            &*potential,
        ),
    ))
    }
}
