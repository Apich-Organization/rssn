//! Validated `NumPy` float-array codec for RSSN and its host applications.

const NPY_MAGIC: &[u8] = b"\x93NUMPY";
const MAX_HEADER_LEN: usize = 10 * 1024;
const PARALLEL_MATRIX_VALUES: usize = 262_144;
const PARALLEL_FORTRAN_VALUES: usize = 524_288;
const MATRIX_ROWS_PER_JOB: usize = 16;
const FORTRAN_VALUES_PER_JOB: usize = 16_384;
const DECODE_ALLOCATION_ERROR: &str = "E0903: cannot allocate the decoded NPY data";

/// Check that the element count and row-major suffix strides fit `usize`.
/// Empty arrays still require representable strides for downstream tensor views.
#[must_use]
pub fn checked_element_count(shape: &[usize]) -> Option<usize> {
    let count = shape
        .iter()
        .try_fold(1usize, |n, extent| n.checked_mul(*extent))?;
    shape
        .iter()
        .rev()
        .try_fold(1usize, |n, extent| n.checked_mul(*extent))?;
    Some(count)
}

/// Encode a dimensionless row-major Float array as NPY 1.0.
///
/// # Errors
/// Rejects inconsistent or unrepresentable shape and header sizes.
pub fn encode(
    shape: &[usize],
    values: &[f64],
) -> Result<Vec<u8>, String> {
    let count = checked_element_count(shape)
        .ok_or_else(|| "E0902: NPY tensor shape or strides overflow".to_string())?;
    if values.len() != count {
        return Err("E0902: NPY tensor storage does not match its shape".to_string());
    }
    if shape.len() > MAX_HEADER_LEN / 2 {
        return Err("E0903: the NPY header exceeds 10 KiB".to_string());
    }
    let shape = if shape.is_empty() {
        "()".to_string()
    } else {
        let rank = shape.len();
        let mut shape = shape
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(", ");
        if rank == 1 {
            shape.push(',');
        }
        format!("({shape})")
    };
    let mut header = format!("{{'descr': '<f8', 'fortran_order': False, 'shape': {shape}, }}");
    // magic(6) + version(2) + header length(2) + header must be a multiple of 64.
    let unpadded_len = header.len().saturating_add(11);
    let padding = 64_usize.saturating_sub(unpadded_len % 64) % 64;
    header.push_str(&" ".repeat(padding));
    header.push('\n');
    if header.len() > MAX_HEADER_LEN {
        return Err("E0903: the NPY header exceeds 10 KiB".to_string());
    }
    let header_len = u16::try_from(header.len())
        .map_err(|_| "E0903: the NPY header is too large".to_string())?;
    let capacity = std::mem::size_of_val(values)
        .checked_add(header.len())
        .and_then(|length| length.checked_add(10))
        .ok_or_else(|| "E0903: the NPY file is too large".to_string())?;
    let mut bytes = Vec::new();
    bytes
        .try_reserve_exact(capacity)
        .map_err(|_| "E0903: cannot allocate the NPY data".to_string())?;
    bytes.extend_from_slice(NPY_MAGIC);
    bytes.push(1); // major
    bytes.push(0); // minor
    bytes.extend_from_slice(&header_len.to_le_bytes());
    bytes.extend_from_slice(header.as_bytes());
    for value in values {
        bytes.extend_from_slice(&value.to_le_bytes());
    }
    Ok(bytes)
}

// NPY uses a Python-literal dictionary, not JSON. Parse the supported
// primitive header grammar directly; never execute or substring-search it.
struct NpyHeader<'a> {
    remaining: &'a [u8],
}
impl<'a> NpyHeader<'a> {
    fn invalid() -> String {
        "E0903: malformed NPY header dictionary".to_string()
    }

    fn take(
        &mut self,
        byte: u8,
    ) -> bool {
        self.remaining = self.remaining.trim_ascii_start();
        if let Some(rest) = self.remaining.strip_prefix(&[byte]) {
            self.remaining = rest;
            true
        } else {
            false
        }
    }

    fn string(&mut self) -> Result<&'a str, String> {
        let Some((&quote @ (b'\'' | b'"'), rest)) = self.remaining.trim_ascii_start().split_first()
        else {
            return Err(Self::invalid());
        };
        for (index, &byte) in rest.iter().enumerate() {
            if byte == quote {
                let (content, closing) = rest.split_at_checked(index).ok_or_else(Self::invalid)?;
                let (_, remaining) = closing.split_first().ok_or_else(Self::invalid)?;
                self.remaining = remaining;
                return std::str::from_utf8(content).map_err(|_| Self::invalid());
            }
            if byte == b'\\' || byte.is_ascii_control() {
                return Err(Self::invalid());
            }
        }
        Err(Self::invalid())
    }

    fn boolean(&mut self) -> Result<bool, String> {
        let trimmed = self.remaining.trim_ascii_start();
        if let Some(rest) = trimmed.strip_prefix(b"True") {
            self.remaining = rest;
            Ok(true)
        } else if let Some(rest) = trimmed.strip_prefix(b"False") {
            self.remaining = rest;
            Ok(false)
        } else {
            Err(Self::invalid())
        }
    }

    fn shape(&mut self) -> Result<Vec<usize>, String> {
        if !self.take(b'(') {
            return Err(Self::invalid());
        }
        let mut shape = Vec::new();
        if self.take(b')') {
            return Ok(shape);
        }
        loop {
            let trimmed = self.remaining.trim_ascii_start();
            let end = trimmed
                .iter()
                .position(|byte| !byte.is_ascii_digit())
                .unwrap_or(trimmed.len());
            if end == 0 {
                return Err(Self::invalid());
            }
            let (digits, rest) = trimmed.split_at_checked(end).ok_or_else(Self::invalid)?;
            self.remaining = rest;
            let extent = std::str::from_utf8(digits)
                .map_err(|_| Self::invalid())?
                .parse::<usize>()
                .map_err(|_| "E0903: NPY shape entry overflows".to_string())?;
            shape.push(extent);
            if self.take(b',') {
                if self.take(b')') {
                    return Ok(shape);
                }
            } else if self.take(b')') && shape.len() > 1 {
                return Ok(shape);
            } else {
                return Err(Self::invalid());
            }
        }
    }

    #[allow(
        clippy::single_call_fn,
        reason = "The Python literal grammar is independent of NPY framing and array order."
    )]
    fn parse(header: &'a str) -> Result<(&'a str, bool, Vec<usize>), String> {
        let mut parser = Self {
            remaining: header.as_bytes(),
        };
        if !parser.take(b'{') {
            return Err(Self::invalid());
        }
        let (mut descr, mut order, mut shape) = (None, None, None);
        if !parser.take(b'}') {
            loop {
                let key = parser.string()?;
                if !parser.take(b':') {
                    return Err(Self::invalid());
                }
                match key {
                    | "descr" if descr.is_none() => descr = Some(parser.string()?),
                    | "fortran_order" if order.is_none() => order = Some(parser.boolean()?),
                    | "shape" if shape.is_none() => shape = Some(parser.shape()?),
                    | _ => return Err(Self::invalid()),
                }
                if parser.take(b',') {
                    if parser.take(b'}') {
                        break;
                    }
                } else if parser.take(b'}') {
                    break;
                } else {
                    return Err(Self::invalid());
                }
            }
        }
        if !parser.remaining.trim_ascii().is_empty() {
            return Err(Self::invalid());
        }
        Ok((
            descr.ok_or_else(Self::invalid)?,
            order.ok_or_else(Self::invalid)?,
            shape.ok_or_else(Self::invalid)?,
        ))
    }
}

/// Decodes a `.npy` file (§8.11).
///
/// Fortran-order input is converted to row-major storage. Large transpositions
/// can share work across the active RSSN CPU pool.
///
/// # Errors
///
/// `E0903` for malformed/truncated data, `E0904` for an unsupported dtype. The message carries the code, because the language-level
/// result is `Result<Tensor<…>, DataError>`.
pub fn decode(bytes: &[u8]) -> Result<(Vec<usize>, Vec<f64>), String> {
    decode_with_limit(bytes, usize::MAX)
}

pub(crate) fn decode_with_limit(
    bytes: &[u8],
    element_limit: usize,
) -> Result<(Vec<usize>, Vec<f64>), String> {
    let Some(rest) = bytes.strip_prefix(NPY_MAGIC) else {
        return Err("E0903: not a NumPy `.npy` file (bad magic)".to_string());
    };
    let Some((&[major, minor], rest)) = rest.split_first_chunk::<2>() else {
        return Err("E0903: truncated NPY version".to_string());
    };
    let truncated = || "E0903: truncated NPY header".to_string();
    let (header_len, rest) = match (major, minor) {
        | (1, 0) => {
            let (bytes, rest) = rest.split_first_chunk::<2>().ok_or_else(truncated)?;
            (usize::from(u16::from_le_bytes(*bytes)), rest)
        },
        | (2 | 3, 0) => {
            let (bytes, rest) = rest.split_first_chunk::<4>().ok_or_else(truncated)?;
            let length = usize::try_from(u32::from_le_bytes(*bytes))
                .map_err(|_| "E0903: the NPY header is too large".to_string())?;
            (length, rest)
        },
        | _ => return Err(format!("E0904: unsupported NPY version {major}.{minor}")),
    };
    if header_len > MAX_HEADER_LEN {
        return Err("E0903: the NPY header exceeds 10 KiB".to_string());
    }
    let (header, data) = rest.split_at_checked(header_len).ok_or_else(truncated)?;
    let header = std::str::from_utf8(header)
        .map_err(|_| "E0903: the NPY header is not UTF-8".to_string())?;

    let (descr, fortran, shape) = NpyHeader::parse(header)?;
    let count = checked_element_count(&shape)
        .ok_or_else(|| "E0903: NPY shape or strides overflow".to_string())?;
    if count > element_limit {
        return Err(format!(
            "E0903: NPY array exceeds the {element_limit} value limit"
        ));
    }
    let reorder = fortran && shape.len() > 1 && count > 0;
    let values = match descr {
        | "<f8" | "|f8" => decode_payload(data, count, &shape, reorder, f64::from_le_bytes)?,
        | ">f8" => decode_payload(data, count, &shape, reorder, f64::from_be_bytes)?,
        | "=f8" => decode_payload(data, count, &shape, reorder, f64::from_ne_bytes)?,
        | "<f4" | "|f4" => {
            decode_payload(data, count, &shape, reorder, |bytes| {
                f64::from(f32::from_le_bytes(bytes))
            })?
        },
        | ">f4" => {
            decode_payload(data, count, &shape, reorder, |bytes| {
                f64::from(f32::from_be_bytes(bytes))
            })?
        },
        | "=f4" => {
            decode_payload(data, count, &shape, reorder, |bytes| {
                f64::from(f32::from_ne_bytes(bytes))
            })?
        },
        | _ => {
            return Err(format!(
                "E0904: unsupported NPY dtype `{descr}` (only float32/float64)"
            ));
        },
    };
    Ok((shape, values))
}

struct FortranAxis {
    stride: usize,
    rewind: usize,
    last: usize,
    position: usize,
}

const fn payload_matches(
    chunks: usize,
    count: usize,
    trailing: &[u8],
) -> bool {
    chunks == count && trailing.is_empty()
}

fn decode_payload<const WIDTH: usize>(
    data: &[u8],
    count: usize,
    shape: &[usize],
    reorder: bool,
    convert: impl Fn([u8; WIDTH]) -> f64 + Sync,
) -> Result<Vec<f64>, String> {
    // Small, low-rank arrays use less metadata with direct indexing. Larger
    // Fortran arrays amortize the coordinate updates over a row.
    if reorder && (count > 8 || shape.len() > 3) {
        let (chunks, trailing) = data.as_chunks::<WIDTH>();
        if payload_matches(chunks.len(), count, trailing) {
            let mut values = Vec::new();
            values
                .try_reserve_exact(count)
                .map_err(|_| DECODE_ALLOCATION_ERROR.to_string())?;
            return decode_fortran(chunks, values, shape, convert);
        }
    }
    decode_values(data, count, shape, reorder, convert)
}

#[inline(never)]
#[expect(
    clippy::single_call_fn,
    reason = "preserves the optimized contiguous copy and the direct tiny-array decoder"
)]
fn decode_values<const WIDTH: usize>(
    data: &[u8],
    count: usize,
    shape: &[usize],
    reorder: bool,
    convert: impl Fn([u8; WIDTH]) -> f64,
) -> Result<Vec<f64>, String> {
    let (chunks, trailing) = data.as_chunks::<WIDTH>();
    if !payload_matches(chunks.len(), count, trailing) {
        let overflow = || "E0903: NPY shape overflows".to_string();
        let expected = count.checked_mul(WIDTH).ok_or_else(overflow)?;
        return Err(format!(
            "E0903: NPY data is {} byte(s), the shape needs {expected}",
            data.len()
        ));
    }

    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| DECODE_ALLOCATION_ERROR.to_string())?;
    if !reorder {
        values.extend(chunks.iter().copied().map(convert));
        return Ok(values);
    }

    let overflow = || "E0903: NPY shape overflows".to_string();
    let mut strides = Vec::with_capacity(shape.len());
    let mut stride = 1_usize;
    for &extent in shape {
        strides.push(stride);
        stride = stride.checked_mul(extent).ok_or_else(overflow)?;
    }
    for index in 0..count {
        let mut remainder = index;
        let mut source = 0_usize;
        for (&extent, &stride) in shape.iter().zip(&strides).rev() {
            // A positive count guarantees that all axes are nonzero. Each coordinate is
            // bounded by its axis and the sum forms an index smaller than the element count.
            let coordinate = remainder.checked_rem(extent).ok_or_else(overflow)?;
            remainder = remainder.checked_div(extent).ok_or_else(overflow)?;
            source = coordinate
                .checked_mul(stride)
                .and_then(|offset| source.checked_add(offset))
                .ok_or_else(overflow)?;
        }
        let bytes = chunks.get(source).ok_or_else(overflow)?;
        values.push(convert(*bytes));
    }
    Ok(values)
}

#[inline(never)]
#[expect(
    clippy::single_call_fn,
    reason = "keeps the row-major decoder independent of the Fortran layout conversion"
)]
fn decode_fortran<const WIDTH: usize>(
    chunks: &[[u8; WIDTH]],
    mut values: Vec<f64>,
    shape: &[usize],
    convert: impl Fn([u8; WIDTH]) -> f64 + Sync,
) -> Result<Vec<f64>, String> {
    let outer_shape = shape
        .iter()
        .rposition(|&extent| extent > 1)
        .and_then(|axis| shape.get(..axis))
        .unwrap_or_default();
    let overflow = || "E0903: NPY shape overflows".to_string();
    let outer_extents = outer_shape.iter().copied().filter(|&extent| extent > 1);
    let mut axes = Vec::with_capacity(outer_extents.clone().count());
    let mut row_count = 1_usize;
    for extent in outer_extents {
        let next = row_count.checked_mul(extent).ok_or_else(overflow)?;
        axes.push(FortranAxis {
            stride: row_count,
            rewind: next.checked_sub(row_count).ok_or_else(overflow)?,
            last: extent.saturating_sub(1),
            position: 0,
        });
        row_count = next;
    }

    // With two varying axes, each output row starts at its row index
    // in the column-major input. Singleton axes do not change the mapping.
    if axes.len() == 1
        && chunks.len() >= PARALLEL_MATRIX_VALUES
        && crate::execution::should_split(row_count.div_ceil(MATRIX_ROWS_PER_JOB), 2)
    {
        let width = chunks.len().checked_div(row_count).ok_or_else(overflow)?;
        values.resize(chunks.len(), 0.0);
        crate::kernels::par::for_chunks_mut(
            &mut values,
            width.saturating_mul(MATRIX_ROWS_PER_JOB),
            2,
            |job, out| {
                for (r, row) in out.chunks_exact_mut(width).enumerate() {
                    let source_row = job.saturating_mul(MATRIX_ROWS_PER_JOB).saturating_add(r);
                    let source = chunks.iter().skip(source_row).step_by(row_count);
                    for (value, &bytes) in row.iter_mut().zip(source) {
                        *value = convert(bytes);
                    }
                }
            },
        );
        return Ok(values);
    }

    if axes.len() > 1 && chunks.len() >= PARALLEL_FORTRAN_VALUES {
        let width = chunks.len().checked_div(row_count).ok_or_else(overflow)?;
        let rows_per_job = FORTRAN_VALUES_PER_JOB.div_ceil(width).min(row_count);
        if crate::execution::should_split(row_count.div_ceil(rows_per_job), 2) {
            values.resize(chunks.len(), 0.0);
            crate::kernels::par::for_chunks_mut(
                &mut values,
                width.saturating_mul(rows_per_job),
                2,
                |job, out| {
                    decode_fortran_rows(
                        chunks,
                        out,
                        &axes,
                        row_count,
                        width,
                        job.saturating_mul(rows_per_job),
                        &convert,
                    );
                },
            );
            return Ok(values);
        }
    }

    // The last varying axis is a row in C order. Its stride in the source is
    // the product of the preceding extents, which is also the number of rows.
    // Each prefix has a distinct source offset below that stride.
    let mut source = 0_usize;
    for _ in 0..row_count {
        let start = chunks.get(source..).ok_or_else(overflow)?;
        values.extend(start.iter().step_by(row_count).copied().map(&convert));
        for axis in axes.iter_mut().rev() {
            if axis.position < axis.last {
                axis.position = axis.position.saturating_add(1);
                source = source.checked_add(axis.stride).ok_or_else(overflow)?;
                break;
            }
            axis.position = 0;
            source = source.checked_sub(axis.rewind).ok_or_else(overflow)?;
        }
    }
    Ok(values)
}

#[expect(
    clippy::arithmetic_side_effects,
    reason = "validated positive dimensions bound coordinates and offsets by the row count"
)]
#[cfg_attr(
    not(test),
    expect(
        clippy::single_call_fn,
        reason = "isolates the independently testable row mapping"
    )
)]
fn decode_fortran_rows<const WIDTH: usize>(
    chunks: &[[u8; WIDTH]],
    out: &mut [f64],
    axes: &[FortranAxis],
    row_count: usize,
    width: usize,
    mut logical_row: usize,
    convert: &impl Fn([u8; WIDTH]) -> f64,
) {
    // Each varying axis consumes at least one bit of the validated element
    // count. This bounds the worker's coordinates without a heap allocation.
    let mut coordinates = [0; usize::BITS as usize];
    let mut source = 0;
    for (axis, coordinate) in axes.iter().zip(&mut coordinates).rev() {
        let extent = axis.last + 1;
        *coordinate = logical_row % extent;
        logical_row /= extent;
        source += *coordinate * axis.stride;
    }

    for row in out.chunks_exact_mut(width) {
        for (value, &bytes) in row
            .iter_mut()
            .zip(chunks.iter().skip(source).step_by(row_count))
        {
            *value = convert(bytes);
        }
        for (axis, coordinate) in axes.iter().zip(&mut coordinates).rev() {
            if *coordinate < axis.last {
                *coordinate += 1;
                source += axis.stride;
                break;
            }
            *coordinate = 0;
            source -= axis.rewind;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fortran_f64_file(
        shape: &[usize],
        source_bits: &[u64],
    ) -> Vec<u8> {
        use std::fmt::Write as _;

        let mut header = String::from("{'shape': (");
        for extent in shape {
            write!(header, "{extent}, ").unwrap();
        }
        header.push_str("), 'descr': '<f8', 'fortran_order': True, }");
        let padding = (64 - ((header.len() + 13) % 64)) % 64;
        header.push_str(&" ".repeat(padding));
        header.push('\n');
        let mut file = b"\x93NUMPY\x03\x00".to_vec();
        file.extend_from_slice(&u32::try_from(header.len()).unwrap().to_le_bytes());
        file.extend_from_slice(header.as_bytes());
        file.extend(source_bits.iter().flat_map(|bits| bits.to_le_bytes()));
        file
    }

    #[test]
    fn fortran_layout_matches_independent_coordinates_and_preserves_bits() {
        for shape in [
            &[][..],
            &[1, 1, 1],
            &[2, 2],
            &[2, 3, 4],
            &[2, 1, 3, 1, 4, 1],
            &[1, 11, 1],
            &[11, 1, 1],
            &[1, 512, 1, 513, 1],
            &[67, 73, 109],
            &[29, 1, 37, 46, 1, 11],
            &[1, 1, 11],
            &[3, 0, 2],
        ] {
            let count = shape.iter().product::<usize>();
            let expected = (0..count)
                .map(|index| {
                    match index {
                        | 0 => (-0.0_f64).to_bits(),
                        | 1 => 0x7ff8_0000_0000_1234,
                        | 2 => f64::NEG_INFINITY.to_bits(),
                        | _ => (f64::from(u32::try_from(index).unwrap()) + 0.5).to_bits(),
                    }
                })
                .collect::<Vec<_>>();
            let source_bits = (0..expected.len())
                .map(|mut offset| {
                    let mut row_major = 0;
                    for &extent in shape {
                        row_major = row_major * extent + offset % extent;
                        offset /= extent;
                    }
                    expected[row_major]
                })
                .collect::<Vec<_>>();
            let (decoded_shape, values) = decode(&fortran_f64_file(shape, &source_bits)).unwrap();
            assert_eq!(decoded_shape, shape);
            assert_eq!(
                values.into_iter().map(f64::to_bits).collect::<Vec<_>>(),
                expected,
                "{shape:?}"
            );
        }
    }

    #[test]
    fn fortran_row_jobs_start_at_mixed_radix_offsets() {
        let rows = 3 * 5 * 7;
        let width = 9;
        let axes = [
            FortranAxis {
                stride: 1,
                rewind: 2,
                last: 2,
                position: 0,
            },
            FortranAxis {
                stride: 3,
                rewind: 12,
                last: 4,
                position: 0,
            },
            FortranAxis {
                stride: 15,
                rewind: 90,
                last: 6,
                position: 0,
            },
        ];
        let payload = (0..rows * width)
            .map(|source| (source as f32).to_be_bytes())
            .collect::<Vec<_>>();
        for first in [0, 1, 6, 7, 34, 35, rows - 4] {
            let mut values = [f64::NAN; 4 * 9];
            decode_fortran_rows(&payload, &mut values, &axes, rows, width, first, &|bytes| {
                f64::from(f32::from_be_bytes(bytes))
            });
            for (offset, value) in values.into_iter().enumerate() {
                let row = first + offset / width;
                let start = row / 35 + (row / 7 % 5) * 3 + (row % 7) * 15;
                assert_eq!(value, (start + (offset % width) * rows) as f64);
            }
        }
    }

    #[test]
    fn fortran_reordering_preserves_framing_and_empty_stride_errors() {
        let empty = [usize::MAX, 0, 2];
        assert_eq!(
            decode(&fortran_f64_file(&empty, &[])).unwrap(),
            (empty.to_vec(), vec![])
        );
        assert_eq!(
            decode(&fortran_f64_file(&[0, usize::MAX, 2], &[])).unwrap_err(),
            "E0903: NPY shape or strides overflow"
        );
        assert_eq!(
            decode(&fortran_f64_file(&[2, 1, 3, 4], &[])).unwrap_err(),
            "E0903: NPY data is 0 byte(s), the shape needs 192"
        );
    }

    #[test]
    fn big_endian_fortran_matrix_decodes_in_row_major_order() {
        let header = b"{'descr': '>f4', 'fortran_order': True, 'shape': (2, 3), }\n";
        let mut bytes = NPY_MAGIC.to_vec();
        bytes.extend_from_slice(&[2, 0]);
        bytes.extend_from_slice(&u32::try_from(header.len()).unwrap().to_le_bytes());
        bytes.extend_from_slice(header);
        for value in [1.0_f32, 4.0, 2.0, 5.0, 3.0, 6.0] {
            bytes.extend_from_slice(&value.to_be_bytes());
        }
        assert_eq!(
            decode(&bytes).unwrap(),
            (vec![2, 3], vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        );
        bytes.push(0);
        assert!(decode(&bytes).unwrap_err().starts_with("E0903: NPY data"));
    }

    #[test]
    fn encoder_rejects_headers_its_decoder_would_reject() {
        assert!(
            encode(&vec![1; MAX_HEADER_LEN / 2], &[0.0])
                .unwrap_err()
                .contains("header exceeds 10 KiB")
        );
        let shape = vec![1; 1000];
        let bytes = encode(&shape, &[3.5]).unwrap();
        assert_eq!(decode(&bytes).unwrap(), (shape, vec![3.5]));
    }
}
