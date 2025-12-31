use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::symbolic::group_theory::Group;
use crate::symbolic::group_theory::GroupElement;
use crate::symbolic::group_theory::Representation;
use crate::symbolic::group_theory::character;

#[unsafe(no_mangle)]

/// Creates a group from a bincode-encoded description.
///
/// The input buffer encodes a [`Group`] (e.g., its underlying set and operation),
/// which is deserialized and returned in canonical internal form.
///
/// # Arguments
///
/// * `buf` - `BincodeBuffer` containing a serialized `Group` description.
///
/// # Returns
///
/// A `BincodeBuffer` containing the canonicalized `Group`, or an empty buffer if
/// deserialization fails.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point; callers
/// must treat the returned buffer as opaque and only pass it to compatible APIs.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_group_create(
    buf: BincodeBuffer
) -> BincodeBuffer {

    let group : Group = match from_bincode_buffer(&buf) {
        | Some(g) => g,
        | None => return BincodeBuffer::empty(),
    };

    to_bincode_buffer(&group)
}

#[unsafe(no_mangle)]

/// Multiplies two group elements using a bincode-encoded group.
///
/// Given a `Group` and two `GroupElement` values, this applies the group
/// operation to compute their product.
///
/// # Arguments
///
/// * `group_buf` - `BincodeBuffer` encoding a [`Group`].
/// * `a_buf` - `BincodeBuffer` encoding a [`GroupElement`] for the left factor.
/// * `b_buf` - `BincodeBuffer` encoding a [`GroupElement`] for the right factor.
///
/// # Returns
///
/// A `BincodeBuffer` encoding the resulting `GroupElement`, or an empty buffer if
/// any input fails to deserialize.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point; callers
/// must treat the returned buffer as opaque and only pass it to compatible APIs.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_group_multiply(
    group_buf: BincodeBuffer,
    a_buf: BincodeBuffer,
    b_buf: BincodeBuffer,
) -> BincodeBuffer {

    let group : Group = match from_bincode_buffer(&group_buf) {
        | Some(g) => g,
        | None => return BincodeBuffer::empty(),
    };

    let a : GroupElement = match from_bincode_buffer(&a_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let b : GroupElement = match from_bincode_buffer(&b_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let result = group.multiply(&a, &b);

    to_bincode_buffer(&result)
}

#[unsafe(no_mangle)]

/// Computes the inverse of a group element using a bincode-encoded group.
///
/// # Arguments
///
/// * `group_buf` - `BincodeBuffer` encoding a [`Group`].
/// * `a_buf` - `BincodeBuffer` encoding a [`GroupElement`] whose inverse is sought.
///
/// # Returns
///
/// A `BincodeBuffer` encoding the inverse `GroupElement`, or an empty buffer if
/// deserialization fails.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point; callers
/// must treat the returned buffer as opaque and only pass it to compatible APIs.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_group_inverse(
    group_buf: BincodeBuffer,
    a_buf: BincodeBuffer,
) -> BincodeBuffer {

    let group : Group = match from_bincode_buffer(&group_buf) {
        | Some(g) => g,
        | None => return BincodeBuffer::empty(),
    };

    let a : GroupElement = match from_bincode_buffer(&a_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let result = group.inverse(&a);

    to_bincode_buffer(&result)
}

#[unsafe(no_mangle)]

/// Tests whether a bincode-encoded group is Abelian.
///
/// A group is Abelian if its operation is commutative, i.e., \(ab = ba\) for all
/// elements \(a, b\).
///
/// # Arguments
///
/// * `group_buf` - `BincodeBuffer` encoding a [`Group`].
///
/// # Returns
///
/// `true` if the group is Abelian, `false` if it is non-Abelian or deserialization
/// fails.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point; callers
/// must ensure the buffer is a valid encoding of a `Group`.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_group_is_abelian(
    group_buf: BincodeBuffer
) -> bool {

    let group: Group =
        match from_bincode_buffer(
            &group_buf,
        ) {
            | Some(g) => g,
            | None => return false,
        };

    group.is_abelian()
}

#[unsafe(no_mangle)]

/// Computes the order of a group element using a bincode-encoded group.
///
/// The order of an element is the smallest positive integer \(n\) such that
/// \(a^n = e\), where \(e\) is the identity element.
///
/// # Arguments
///
/// * `group_buf` - `BincodeBuffer` encoding a [`Group`].
/// * `a_buf` - `BincodeBuffer` encoding a [`GroupElement`].
///
/// # Returns
///
/// The order of the element as a `usize`, or `0` if the order is undefined or
/// deserialization fails.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point; callers
/// must ensure the buffers encode a compatible group and element.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_group_element_order(
    group_buf: BincodeBuffer,
    a_buf: BincodeBuffer,
) -> usize {

    let group: Group =
        match from_bincode_buffer(
            &group_buf,
        ) {
            | Some(g) => g,
            | None => return 0,
        };

    let a: GroupElement =
        match from_bincode_buffer(
            &a_buf,
        ) {
            | Some(e) => e,
            | None => return 0,
        };

    group
        .element_order(&a)
        .unwrap_or(0)
}

#[unsafe(no_mangle)]

/// Computes the conjugacy classes of a bincode-encoded group.
///
/// Conjugacy classes partition the group into sets of elements related by
/// conjugation \(gag^{-1}\).
///
/// # Arguments
///
/// * `group_buf` - `BincodeBuffer` encoding a [`Group`].
///
/// # Returns
///
/// A `BincodeBuffer` encoding the conjugacy classes (typically as a collection of
/// `Vec<GroupElement>`), or an empty buffer if deserialization fails.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point; callers
/// must treat the returned buffer as opaque and only pass it to compatible APIs.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_group_conjugacy_classes(
    group_buf: BincodeBuffer
) -> BincodeBuffer {

    let group : Group = match from_bincode_buffer(&group_buf) {
        | Some(g) => g,
        | None => return BincodeBuffer::empty(),
    };

    let classes =
        group.conjugacy_classes();

    to_bincode_buffer(&classes)
}

#[unsafe(no_mangle)]

/// Computes the center of a bincode-encoded group.
///
/// The center is the subset of elements that commute with every group element,
/// i.e., \(Z(G) = \{z \in G \mid zg = gz, \forall g \in G\}\).
///
/// # Arguments
///
/// * `group_buf` - `BincodeBuffer` encoding a [`Group`].
///
/// # Returns
///
/// A `BincodeBuffer` encoding the center (typically as a collection of
/// `GroupElement` values), or an empty buffer if deserialization fails.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point; callers
/// must treat the returned buffer as opaque and only pass it to compatible APIs.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_group_center(
    group_buf: BincodeBuffer
) -> BincodeBuffer {

    let group : Group = match from_bincode_buffer(&group_buf) {
        | Some(g) => g,
        | None => return BincodeBuffer::empty(),
    };

    let center = group.center();

    to_bincode_buffer(&center)
}

#[unsafe(no_mangle)]

/// Creates a group representation from a bincode-encoded description.
///
/// A representation assigns linear operators to each group element, typically as
/// matrices acting on a vector space.
///
/// # Arguments
///
/// * `buf` - `BincodeBuffer` containing a serialized [`Representation`].
///
/// # Returns
///
/// A `BincodeBuffer` containing the canonicalized `Representation`, or an empty
/// buffer if deserialization fails.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point; callers
/// must treat the returned buffer as opaque and only pass it to compatible APIs.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_representation_create(
    buf: BincodeBuffer
) -> BincodeBuffer {

    let rep : Representation = match from_bincode_buffer(&buf) {
        | Some(r) => r,
        | None => return BincodeBuffer::empty(),
    };

    to_bincode_buffer(&rep)
}

#[unsafe(no_mangle)]

/// Checks whether a representation is valid for a given group using bincode serialization.
///
/// A representation is valid if it respects the group operation, i.e., the image
/// of the group under the representation is a homomorphism.
///
/// # Arguments
///
/// * `rep_buf` - `BincodeBuffer` encoding a [`Representation`].
/// * `group_buf` - `BincodeBuffer` encoding a [`Group`].
///
/// # Returns
///
/// `true` if the representation is valid for the group, `false` otherwise or if
/// deserialization fails.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point; callers
/// must ensure the buffers encode compatible objects.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_representation_is_valid(
    rep_buf: BincodeBuffer,
    group_buf: BincodeBuffer,
) -> bool {

    let rep: Representation =
        match from_bincode_buffer(
            &rep_buf,
        ) {
            | Some(r) => r,
            | None => return false,
        };

    let group: Group =
        match from_bincode_buffer(
            &group_buf,
        ) {
            | Some(g) => g,
            | None => return false,
        };

    rep.is_valid(&group)
}

#[unsafe(no_mangle)]

/// Computes the character of a representation using bincode serialization.
///
/// The character is the trace of each representation matrix over the group,
/// viewed as a class function on the group.
///
/// # Arguments
///
/// * `rep_buf` - `BincodeBuffer` encoding a [`Representation`].
///
/// # Returns
///
/// A `BincodeBuffer` encoding the character values (e.g., as a vector indexed by
/// group elements or conjugacy classes), or an empty buffer if deserialization fails.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point; callers
/// must treat the returned buffer as opaque and only pass it to compatible APIs.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_character(
    rep_buf: BincodeBuffer
) -> BincodeBuffer {

    let rep : Representation = match from_bincode_buffer(&rep_buf) {
        | Some(r) => r,
        | None => return BincodeBuffer::empty(),
    };

    let chars = character(&rep);

    to_bincode_buffer(&chars)
}
