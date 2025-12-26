use crate::symbolic::core::Expr;
use crate::symbolic::quantum_field_theory;
use std::os::raw::c_char;

/// Computes the Dirac adjoint of a fermion field.
#[no_mangle]

pub unsafe extern "C" fn rssn_dirac_adjoint(psi: *const Expr) -> *mut Expr {

    if psi.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        quantum_field_theory::dirac_adjoint(&*psi),
    ))
}

/// Computes the Feynman slash notation.
#[no_mangle]

pub unsafe extern "C" fn rssn_feynman_slash(v_mu: *const Expr) -> *mut Expr {

    if v_mu.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        quantum_field_theory::feynman_slash(&*v_mu),
    ))
}

/// Lagrangian density for a free real scalar field.
#[no_mangle]

pub unsafe extern "C" fn rssn_scalar_field_lagrangian(
    phi: *const Expr,
    m: *const Expr,
) -> *mut Expr {

    if phi.is_null() || m.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        quantum_field_theory::scalar_field_lagrangian(&*phi, &*m),
    ))
}

/// Lagrangian density for QED.
#[no_mangle]

pub unsafe extern "C" fn rssn_qed_lagrangian(
    psi_bar: *const Expr,
    psi: *const Expr,
    a_mu: *const Expr,
    m: *const Expr,
    e: *const Expr,
) -> *mut Expr {

    if psi_bar.is_null() || psi.is_null() || a_mu.is_null() || m.is_null() || e.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        quantum_field_theory::qed_lagrangian(
            &*psi_bar, &*psi, &*a_mu, &*m, &*e,
        ),
    ))
}

/// Lagrangian density for QCD.
#[no_mangle]

pub unsafe extern "C" fn rssn_qcd_lagrangian(
    psi_bar: *const Expr,
    psi: *const Expr,
    g_mu: *const Expr,
    m: *const Expr,
    gs: *const Expr,
) -> *mut Expr {

    if psi_bar.is_null() || psi.is_null() || g_mu.is_null() || m.is_null() || gs.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        quantum_field_theory::qcd_lagrangian(
            &*psi_bar, &*psi, &*g_mu, &*m, &*gs,
        ),
    ))
}

/// Computes a propagator for a particle in QFT.
#[no_mangle]

pub unsafe extern "C" fn rssn_qft_propagator(
    p: *const Expr,
    m: *const Expr,
    is_fermion: bool,
) -> *mut Expr {

    if p.is_null() || m.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        quantum_field_theory::propagator(&*p, &*m, is_fermion),
    ))
}

/// Scattering cross-section.
#[no_mangle]

pub unsafe extern "C" fn rssn_qft_scattering_cross_section(
    matrix_element: *const Expr,
    flux: *const Expr,
    phase_space: *const Expr,
) -> *mut Expr {

    if matrix_element.is_null() || flux.is_null() || phase_space.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        quantum_field_theory::scattering_cross_section(
            &*matrix_element,
            &*flux,
            &*phase_space,
        ),
    ))
}

/// Feynman propagator in position space.
#[no_mangle]

pub unsafe extern "C" fn rssn_feynman_propagator_position_space(
    x: *const Expr,
    y: *const Expr,
    m: *const Expr,
) -> *mut Expr {

    if x.is_null() || y.is_null() || m.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        quantum_field_theory::feynman_propagator_position_space(&*x, &*y, &*m),
    ))
}
