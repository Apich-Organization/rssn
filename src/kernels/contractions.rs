//! Rectangular contractions over caller-defined scalar arithmetic.

/// Failure of a checked rectangular contraction.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ContractionError<E> {
    /// The result dimensions cannot be represented or allocated.
    ShapeOverflow,
    /// The caller's scalar arithmetic or indexed access failed.
    Element(E),
}

/// Compute a rectangular product in row-major, first-product accumulation order.
///
/// `accumulate(row, inner, column, previous)` performs one scalar product and,
/// when present, adds it to `previous`. This supports exact, complex and physical
/// host scalars without importing their runtime representation into the kernel.
/// Zero-length contractions obtain each additive identity from `zero`.
///
/// # Errors
/// Returns [`ContractionError::ShapeOverflow`] when the output allocation size
/// overflows, or [`ContractionError::Element`] for the first scalar failure.
pub fn contract_rectangular<T, E>(
    rows: usize,
    inner: usize,
    cols: usize,
    mut zero: impl FnMut() -> T,
    mut accumulate: impl FnMut(usize, usize, usize, Option<T>) -> Result<T, E>,
) -> Result<Vec<T>, ContractionError<E>> {
    let size = rows
        .checked_mul(cols)
        .ok_or(ContractionError::ShapeOverflow)?;
    if size
        .checked_mul(std::mem::size_of::<T>())
        .is_none_or(|bytes| bytes > isize::MAX as usize)
    {
        return Err(ContractionError::ShapeOverflow);
    }
    let mut output = Vec::with_capacity(size);
    for row in 0..rows {
        for column in 0..cols {
            let mut sum = None;
            for index in 0..inner {
                sum = Some(accumulate(row, index, column, sum).map_err(ContractionError::Element)?);
            }
            output.push(sum.unwrap_or_else(&mut zero));
        }
    }
    Ok(output)
}
