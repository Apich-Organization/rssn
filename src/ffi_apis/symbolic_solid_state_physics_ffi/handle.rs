//! Handle-based FFI API for solid-state physics functions.

use crate::symbolic::core::Expr;
use crate::symbolic::solid_state_physics::{CrystalLattice, density_of_states_3d, fermi_energy_3d, drude_conductivity, hall_coefficient};
use crate::symbolic::vector::Vector;

/// Creates a new `CrystalLattice`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_crystal_lattice_new(
    a1: *const Vector,
    a2: *const Vector,
    a3: *const Vector,
) -> *mut CrystalLattice { unsafe {

    if a1.is_null()
        || a2.is_null()
        || a3.is_null()
    {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        CrystalLattice::new(
            (*a1).clone(),
            (*a2).clone(),
            (*a3).clone(),
        ),
    ))
}}

/// Frees a `CrystalLattice`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_crystal_lattice_free(
    ptr: *mut CrystalLattice
) { unsafe {

    if !ptr.is_null() {

        let _ = Box::from_raw(ptr);
    }
}}

/// Computes the volume of the unit cell.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_crystal_lattice_volume(
    ptr: *const CrystalLattice
) -> *mut Expr { unsafe {

    if ptr.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        (*ptr).volume(),
    ))
}}

/// Computes reciprocal lattice vectors.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_crystal_lattice_reciprocal_vectors(
    ptr: *const CrystalLattice,
    b1: *mut *mut Vector,
    b2: *mut *mut Vector,
    b3: *mut *mut Vector,
) { unsafe {

    if ptr.is_null()
        || b1.is_null()
        || b2.is_null()
        || b3.is_null()
    {

        return;
    }

    let (v1, v2, v3) = (*ptr)
        .reciprocal_lattice_vectors();

    *b1 = Box::into_raw(Box::new(v1));

    *b2 = Box::into_raw(Box::new(v2));

    *b3 = Box::into_raw(Box::new(v3));
}}

/// Computes the density of states for a 3D electron gas.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_density_of_states_3d(
    energy: *const Expr,
    effective_mass: *const Expr,
    volume: *const Expr,
) -> *mut Expr { unsafe {

    if energy.is_null()
        || effective_mass.is_null()
        || volume.is_null()
    {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        density_of_states_3d(
            &*energy,
            &*effective_mass,
            &*volume,
        ),
    ))
}}

/// Computes Fermi energy for a 3D electron gas.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_fermi_energy_3d(
    concentration: *const Expr,
    effective_mass: *const Expr,
) -> *mut Expr { unsafe {

    if concentration.is_null()
        || effective_mass.is_null()
    {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        fermi_energy_3d(
            &*concentration,
            &*effective_mass,
        ),
    ))
}}

/// Computes Drude conductivity.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_drude_conductivity(
    n: *const Expr,
    e_charge: *const Expr,
    tau: *const Expr,
    m_star: *const Expr,
) -> *mut Expr { unsafe {

    if n.is_null()
        || e_charge.is_null()
        || tau.is_null()
        || m_star.is_null()
    {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        drude_conductivity(
            &*n,
            &*e_charge,
            &*tau,
            &*m_star,
        ),
    ))
}}

/// Computes Hall coefficient.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_hall_coefficient(
    n: *const Expr,
    q: *const Expr,
) -> *mut Expr { unsafe {

    if n.is_null() || q.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        hall_coefficient(&*n, &*q),
    ))
}}
