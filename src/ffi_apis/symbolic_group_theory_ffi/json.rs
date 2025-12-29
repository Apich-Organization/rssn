use std::collections::HashMap;
use std::os::raw::c_char;

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::group_theory::*;

#[no_mangle]

/// Creates a group from a JSON-encoded description.
///
/// The input string encodes a [`Group`] (e.g., its elements and multiplication law),
/// which is deserialized and returned in canonical internal form.
///
/// # Arguments
///
/// * `json_str` - C string pointer containing JSON for a `Group` description.
///
/// # Returns
///
/// A C string pointer containing JSON-encoded [`Group`], or null if the input
/// cannot be deserialized.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw C string pointer and
/// returns ownership of a heap-allocated C string that must be freed by the caller.
pub unsafe extern "C" fn rssn_json_group_create(
    json_str: *const c_char
) -> *mut c_char {

    let group : Group = match from_json_string(json_str) {
        | Some(g) => g,
        | None => return std::ptr::null_mut(),
    };

    to_json_string(&group)
}

#[no_mangle]

/// Multiplies two group elements using a JSON-encoded group and elements.
///
/// # Arguments
///
/// * `group_json` - C string pointer with JSON-encoded [`Group`].
/// * `a_json` - C string pointer with JSON-encoded [`GroupElement`] for the left factor.
/// * `b_json` - C string pointer with JSON-encoded [`GroupElement`] for the right factor.
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `GroupElement` representing the
/// product, or null if deserialization fails or the multiplication is undefined.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw C string pointers and
/// returns ownership of a heap-allocated C string that must be freed by the caller.
pub unsafe extern "C" fn rssn_json_group_multiply(
    group_json: *const c_char,
    a_json: *const c_char,
    b_json: *const c_char,
) -> *mut c_char {

    let group : Group = match from_json_string(group_json) {
        | Some(g) => g,
        | None => return std::ptr::null_mut(),
    };

    let a : GroupElement = match from_json_string(a_json) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    let b : GroupElement = match from_json_string(b_json) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    let result = group.multiply(&a, &b);

    to_json_string(&result)
}

#[no_mangle]

/// Computes the inverse of a group element using a JSON-encoded group.
///
/// # Arguments
///
/// * `group_json` - C string pointer with JSON-encoded [`Group`].
/// * `a_json` - C string pointer with JSON-encoded [`GroupElement`] whose inverse is sought.
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `GroupElement` for the inverse, or
/// null if deserialization fails or the element has no inverse in the group.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw C string pointers and
/// returns ownership of a heap-allocated C string that must be freed by the caller.
pub unsafe extern "C" fn rssn_json_group_inverse(
    group_json: *const c_char,
    a_json: *const c_char,
) -> *mut c_char {

    let group : Group = match from_json_string(group_json) {
        | Some(g) => g,
        | None => return std::ptr::null_mut(),
    };

    let a : GroupElement = match from_json_string(a_json) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    let result = group.inverse(&a);

    to_json_string(&result)
}

#[no_mangle]

/// Tests whether a group is Abelian using a JSON-encoded group.
///
/// A group is Abelian if its operation is commutative for all pairs of elements.
///
/// # Arguments
///
/// * `group_json` - C string pointer with JSON-encoded [`Group`].
///
/// # Returns
///
/// `true` if the group is Abelian, `false` otherwise or if deserialization fails.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw C string pointer; the
/// caller must ensure it points to a valid JSON string.
pub unsafe extern "C" fn rssn_json_group_is_abelian(
    group_json: *const c_char
) -> bool {

    let group: Group =
        match from_json_string(
            group_json,
        ) {
            | Some(g) => g,
            | None => return false,
        };

    group.is_abelian()
}

#[no_mangle]

/// Computes the order of a group element using a JSON-encoded group.
///
/// The order is the smallest positive integer `n` such that `a^n` is the identity,
/// if such an integer exists.
///
/// # Arguments
///
/// * `group_json` - C string pointer with JSON-encoded [`Group`].
/// * `a_json` - C string pointer with JSON-encoded [`GroupElement`].
///
/// # Returns
///
/// The order of the element as a `usize`, or `0` if the order is undefined or
/// deserialization fails.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw C string pointers; the
/// caller must ensure they point to valid JSON strings.
pub unsafe extern "C" fn rssn_json_group_element_order(
    group_json: *const c_char,
    a_json: *const c_char,
) -> usize {

    let group: Group =
        match from_json_string(
            group_json,
        ) {
            | Some(g) => g,
            | None => return 0,
        };

    let a: GroupElement =
        match from_json_string(a_json) {
            | Some(e) => e,
            | None => return 0,
        };

    group
        .element_order(&a)
        .unwrap_or(0)
}

#[no_mangle]

/// Computes the conjugacy classes of a group using a JSON-encoded group.
///
/// Conjugacy classes partition the group into sets of elements related by
/// conjugation and are returned as JSON-encoded collections.
///
/// # Arguments
///
/// * `group_json` - C string pointer with JSON-encoded [`Group`].
///
/// # Returns
///
/// A C string pointer containing JSON-encoded conjugacy classes (typically as
/// lists of `GroupElement` values), or null if deserialization fails.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw C string pointer and
/// returns ownership of a heap-allocated C string that must be freed by the caller.
pub unsafe extern "C" fn rssn_json_group_conjugacy_classes(
    group_json: *const c_char
) -> *mut c_char {

    let group : Group = match from_json_string(group_json) {
        | Some(g) => g,
        | None => return std::ptr::null_mut(),
    };

    let classes =
        group.conjugacy_classes();

    to_json_string(&classes)
}

#[no_mangle]

/// Computes the center of a group using a JSON-encoded group.
///
/// The center consists of elements that commute with every element of the group.
///
/// # Arguments
///
/// * `group_json` - C string pointer with JSON-encoded [`Group`].
///
/// # Returns
///
/// A C string pointer containing JSON-encoded center elements, or null if
/// deserialization fails.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw C string pointer and
/// returns ownership of a heap-allocated C string that must be freed by the caller.
pub unsafe extern "C" fn rssn_json_group_center(
    group_json: *const c_char
) -> *mut c_char {

    let group : Group = match from_json_string(group_json) {
        | Some(g) => g,
        | None => return std::ptr::null_mut(),
    };

    let center = group.center();

    to_json_string(&center)
}

#[no_mangle]

/// Creates a group representation from a JSON-encoded description.
///
/// A representation assigns matrices or linear operators to group elements,
/// defining a homomorphism into a general linear group.
///
/// # Arguments
///
/// * `json_str` - C string pointer containing JSON for a [`Representation`].
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `Representation`, or null if the
/// input cannot be deserialized.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw C string pointer and
/// returns ownership of a heap-allocated C string that must be freed by the caller.
pub unsafe extern "C" fn rssn_json_representation_create(
    json_str: *const c_char
) -> *mut c_char {

    let rep : Representation = match from_json_string(json_str) {
        | Some(r) => r,
        | None => return std::ptr::null_mut(),
    };

    to_json_string(&rep)
}

#[no_mangle]

/// Checks whether a representation is valid for a given group using JSON-encoded inputs.
///
/// A representation is valid if it respects the group operation, i.e., defines a
/// group homomorphism.
///
/// # Arguments
///
/// * `rep_json` - C string pointer with JSON-encoded [`Representation`].
/// * `group_json` - C string pointer with JSON-encoded [`Group`].
///
/// # Returns
///
/// `true` if the representation is valid for the group, `false` otherwise or if
/// deserialization fails.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw C string pointers; the
/// caller must ensure they point to valid JSON strings.
pub unsafe extern "C" fn rssn_json_representation_is_valid(
    rep_json: *const c_char,
    group_json: *const c_char,
) -> bool {

    let rep: Representation =
        match from_json_string(rep_json)
        {
            | Some(r) => r,
            | None => return false,
        };

    let group: Group =
        match from_json_string(
            group_json,
        ) {
            | Some(g) => g,
            | None => return false,
        };

    rep.is_valid(&group)
}

#[no_mangle]

/// Computes the character of a representation using a JSON-encoded input.
///
/// The character is the trace of each representation matrix and is returned as a
/// JSON-encoded mapping from group elements to scalar values.
///
/// # Arguments
///
/// * `rep_json` - C string pointer with JSON-encoded [`Representation`].
///
/// # Returns
///
/// A C string pointer containing JSON-encoded character values, or null if
/// deserialization fails.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw C string pointer and
/// returns ownership of a heap-allocated C string that must be freed by the caller.
pub unsafe extern "C" fn rssn_json_character(
    rep_json: *const c_char
) -> *mut c_char {

    let rep : Representation = match from_json_string(rep_json) {
        | Some(r) => r,
        | None => return std::ptr::null_mut(),
    };

    let chars = character(&rep);

    to_json_string(&chars)
}
