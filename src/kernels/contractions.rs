//! Rectangular contractions over caller-defined scalar arithmetic.

/// Failure of a checked rectangular contraction.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ContractionError<E> {
    /// The result dimensions cannot be represented or allocated.
    ShapeOverflow,
    /// The caller's scalar arithmetic, indexed access, or cell consumer failed.
    Element(E),
}

/// Compute a rectangular product in row-major, first-product accumulation order.
///
/// `accumulate(row, index, column, previous)` performs one scalar product and,
/// when present, adds it to `previous`. This supports exact, complex and physical
/// host scalars without importing their runtime representation into the kernel.
/// Zero-length contractions obtain each additive identity from `zero`.
///
/// # Errors
/// Returns [`ContractionError::ShapeOverflow`] when the output size is invalid
/// or the allocation fails, or [`ContractionError::Element`] for the first scalar failure.
pub fn contract_rectangular<T, E>(
    rows: usize,
    inner: usize,
    cols: usize,
    zero: impl FnMut() -> T,
    accumulate: impl FnMut(usize, usize, usize, Option<T>) -> Result<T, E>,
) -> Result<Vec<T>, ContractionError<E>> {
    let mut output = Vec::new();
    contract_rectangular_append(rows, inner, cols, &mut output, zero, accumulate)?;
    Ok(output)
}

/// Appends a rectangular contraction to a caller-owned vector in row-major order.
///
/// Existing elements are preserved. Clear the vector first to reuse its allocation for a
/// new result. The callback and identity order is the same as for [`contract_rectangular`].
///
/// # Errors
/// Invalid output size, allocation failure, or the first scalar error returns an error.
/// The vector's original contents and length are restored on error or unwinding.
pub fn contract_rectangular_append<T, E>(
    rows: usize,
    inner: usize,
    cols: usize,
    output: &mut Vec<T>,
    zero: impl FnMut() -> T,
    accumulate: impl FnMut(usize, usize, usize, Option<T>) -> Result<T, E>,
) -> Result<(), ContractionError<E>> {
    let size = rows
        .checked_mul(cols)
        .ok_or(ContractionError::ShapeOverflow)?;
    if size == 0 {
        return Ok(());
    }
    output
        .try_reserve_exact(size)
        .map_err(|_| ContractionError::ShapeOverflow)?;
    let mut append = Appended {
        keep: output.len(),
        output,
    };
    contract_cells(rows, inner, cols, zero, accumulate, |_, _, value| {
        append.output.push(value);
        Ok(())
    })?;
    append.keep = append.output.len();
    Ok(())
}

/// Sends each completed cell to `consume(row, column, value)` without collecting the output.
///
/// Cells are visited in row-major order and use the same scalar callback order as
/// [`contract_rectangular`]. The kernel retains only the current scalar accumulator.
///
/// # Errors
/// Returns [`ContractionError::ShapeOverflow`] if the number of output cells is
/// unrepresentable, or [`ContractionError::Element`] for the first scalar or consumer
/// error. Earlier calls to the consumer remain effective.
pub fn contract_rectangular_for_each<T, E>(
    rows: usize,
    inner: usize,
    cols: usize,
    zero: impl FnMut() -> T,
    accumulate: impl FnMut(usize, usize, usize, Option<T>) -> Result<T, E>,
    consume: impl FnMut(usize, usize, T) -> Result<(), E>,
) -> Result<(), ContractionError<E>> {
    let size = rows
        .checked_mul(cols)
        .ok_or(ContractionError::ShapeOverflow)?;
    if size == 0 {
        return Ok(());
    }
    contract_cells(rows, inner, cols, zero, accumulate, consume)
}

struct Appended<'a, T> {
    output: &'a mut Vec<T>,
    keep: usize,
}

impl<T> Drop for Appended<'_, T> {
    fn drop(&mut self) {
        self.output.truncate(self.keep);
    }
}

fn contract_cells<T, E>(
    rows: usize,
    inner: usize,
    cols: usize,
    mut zero: impl FnMut() -> T,
    mut accumulate: impl FnMut(usize, usize, usize, Option<T>) -> Result<T, E>,
    mut consume: impl FnMut(usize, usize, T) -> Result<(), E>,
) -> Result<(), ContractionError<E>> {
    for row in 0..rows {
        for column in 0..cols {
            let mut sum = None;
            for index in 0..inner {
                sum = Some(accumulate(row, index, column, sum).map_err(ContractionError::Element)?);
            }
            consume(row, column, sum.unwrap_or_else(&mut zero))
                .map_err(ContractionError::Element)?;
        }
    }
    Ok(())
}
