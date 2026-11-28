//! Typed column storage shared by every frame operation.

/// The logical type of a column's values.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DType {
    /// Signed 64-bit integers.
    Int,
    /// Double-precision floating point numbers.
    Float,
    /// Booleans.
    Bool,
    /// UTF-8 strings.
    Str,
    /// Raw byte sequences.
    Bytes,
    /// Variable-length sequences of one homogeneous element type.
    List,
}

/// A homogeneous column that may carry missing values.
///
/// Nulls are expressed with [`Column::Nullable`]: the inner column keeps a
/// placeholder value at every missing position and the mask marks which
/// positions are absent (`true` means missing). Nesting a nullable column
/// inside another one is allowed but unusual.
#[derive(Debug, Clone, PartialEq)]
pub enum Column {
    /// Signed 64-bit integers.
    Int(Vec<i64>),
    /// Double-precision floats.
    Float(Vec<f64>),
    /// Booleans.
    Bool(Vec<bool>),
    /// UTF-8 strings.
    Str(Vec<String>),
    /// Raw byte sequences.
    Bytes(Vec<Vec<u8>>),
    /// A variable-length sequence of values.
    ///
    /// `values` holds every element in row order and `offsets` has one more
    /// entry than there are rows: row `r` owns
    /// `values[offsets[r]..offsets[r + 1]]`. An empty row is an empty slice, so
    /// a null row needs the [`Column::Nullable`] wrapper.
    List {
        /// The flattened elements.
        values: Box<Self>,
        /// Row boundaries into `values`.
        offsets: Vec<usize>,
    },
    /// An inner column plus a mask where `true` marks a missing value.
    Nullable(Box<Self>, Vec<bool>),
}

impl Column {
    /// Returns the number of slots in the column, including missing ones.
    #[must_use]
    pub fn len(&self) -> usize {
        match self {
            | Self::Int(values) => values.len(),
            | Self::Float(values) => values.len(),
            | Self::Bool(values) => values.len(),
            | Self::Str(values) => values.len(),
            | Self::Bytes(values) => values.len(),
            | Self::List { offsets, .. } => offsets.len().saturating_sub(1),
            | Self::Nullable(_, mask) => mask.len(),
        }
    }

    /// Returns `true` when the column holds no slots at all.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the logical type, looking through any nullable wrapper.
    #[must_use]
    pub fn dtype(&self) -> DType {
        match self {
            | Self::Int(_) => DType::Int,
            | Self::Float(_) => DType::Float,
            | Self::Bool(_) => DType::Bool,
            | Self::Str(_) => DType::Str,
            | Self::Bytes(_) => DType::Bytes,
            | Self::List { .. } => DType::List,
            | Self::Nullable(inner, _) => inner.dtype(),
        }
    }

    /// Returns a numeric view of the value at `index`.
    ///
    /// Integers and booleans (`false` is `0.0`, `true` is `1.0`) are widened
    /// to `f64`; strings and raw bytes have no numeric view. Missing positions
    /// and out-of-range indices return `None`.
    #[must_use]
    pub fn get_f64(
        &self,
        index: usize,
    ) -> Option<f64> {
        match self {
            | Self::Int(values) => values.get(index).map(|value| *value as f64),
            | Self::Float(values) => values.get(index).copied(),
            | Self::Bool(values) => {
                values
                    .get(index)
                    .map(|value| if *value { 1.0 } else { 0.0 })
            },
            | Self::Str(_) | Self::Bytes(_) | Self::List { .. } => None,
            | Self::Nullable(inner, mask) => {
                if mask.get(index).copied().unwrap_or(true) {
                    None
                } else {
                    inner.get_f64(index)
                }
            },
        }
    }

    /// Returns `true` when the slot at `index` is missing.
    ///
    /// Out-of-range indices count as missing.
    #[must_use]
    pub fn is_null(
        &self,
        index: usize,
    ) -> bool {
        match self {
            | Self::Nullable(_, mask) => mask.get(index).copied().unwrap_or(true),
            | _ => false,
        }
    }

    /// Returns the short name of the logical type (`"int"`, `"float"`,
    /// `"bool"`, `"str"`, `"bytes"` or `"list"`).
    #[must_use]
    pub fn name_of(&self) -> &'static str {
        match self.dtype() {
            | DType::Int => "int",
            | DType::Float => "float",
            | DType::Bool => "bool",
            | DType::Str => "str",
            | DType::Bytes => "bytes",
            | DType::List => "list",
        }
    }

    /// Returns `true` when the column carries a null mask.
    #[must_use]
    pub fn is_nullable(&self) -> bool {
        matches!(self, Self::Nullable(_, _))
    }

    /// Returns the null mask, when the column has one.
    #[must_use]
    pub fn null_mask(&self) -> Option<&[bool]> {
        match self {
            | Self::Nullable(_, mask) => Some(mask),
            | _ => None,
        }
    }

    /// Builds a nullable column, checking that the mask matches the values.
    ///
    /// # Errors
    ///
    /// Returns an error when the mask length differs from the inner column.
    pub fn new_nullable(
        inner: Self,
        mask: Vec<bool>,
    ) -> Result<Self, String> {
        if inner.len() != mask.len() {
            return Err(format!(
                "null mask has {} entries but the column has {} values",
                mask.len(),
                inner.len()
            ));
        }
        Ok(Self::Nullable(Box::new(inner), mask))
    }

    /// Returns a copy of the slots in `start..end`.
    ///
    /// Out-of-range bounds are clamped to the column length.
    #[must_use]
    pub fn slice(
        &self,
        start: usize,
        end: usize,
    ) -> Self {
        let len = self.len();
        let start = start.min(len);
        let end = end.min(len).max(start);
        match self {
            | Self::Int(values) => Self::Int(values[start..end].to_vec()),
            | Self::Float(values) => Self::Float(values[start..end].to_vec()),
            | Self::Bool(values) => Self::Bool(values[start..end].to_vec()),
            | Self::Str(values) => Self::Str(values[start..end].to_vec()),
            | Self::Bytes(values) => Self::Bytes(values[start..end].to_vec()),
            | Self::List { values, offsets } => {
                if offsets.len() < 2 {
                    return self.clone();
                }
                let base = offsets[start];
                let tail = offsets[end];
                let shifted: Vec<usize> = offsets[start..=end]
                    .iter()
                    .map(|offset| offset.saturating_sub(base))
                    .collect();
                Self::List {
                    values: Box::new(values.slice(base, tail)),
                    offsets: shifted,
                }
            },
            | Self::Nullable(inner, mask) => {
                Self::Nullable(Box::new(inner.slice(start, end)), mask[start..end].to_vec())
            },
        }
    }

    /// Builds a list column from a flattened element column and row offsets.
    ///
    /// # Errors
    ///
    /// Returns an error when `offsets` is not a non-decreasing run that starts
    /// at zero and ends at `values.len()`.
    pub fn new_list(
        values: Self,
        offsets: Vec<usize>,
    ) -> Result<Self, String> {
        if offsets.is_empty() {
            if values.is_empty() {
                return Ok(Self::List {
                    values: Box::new(values),
                    offsets: vec![0],
                });
            }
            return Err("a list column needs at least one offset".to_string());
        }
        if offsets[0] != 0 {
            return Err("the first list offset must be zero".to_string());
        }
        for pair in offsets.windows(2) {
            if pair[1] < pair[0] {
                return Err("list offsets must not decrease".to_string());
            }
        }
        let last = offsets.last().copied().unwrap_or(0);
        if last != values.len() {
            return Err(format!(
                "the last list offset is {last} but there are {} elements",
                values.len()
            ));
        }
        Ok(Self::List {
            values: Box::new(values),
            offsets,
        })
    }

    /// Returns the elements of row `index`, or `None` when the index is out of
    /// range or the slot is missing.
    #[must_use]
    pub fn list_row(
        &self,
        index: usize,
    ) -> Option<Self> {
        let (start, end) = self.list_bounds(index)?;
        let Self::List { values, .. } = self else {
            return None;
        };
        Some(values.slice(start, end))
    }

    /// Returns the element range of row `index`.
    #[must_use]
    pub fn list_bounds(
        &self,
        index: usize,
    ) -> Option<(usize, usize)> {
        match self {
            | Self::List { offsets, .. } => {
                let start = *offsets.get(index)?;
                let end = *offsets.get(index.checked_add(1)?)?;
                Some((start, end))
            },
            | Self::Nullable(inner, mask) => {
                if mask.get(index).copied().unwrap_or(true) {
                    None
                } else {
                    inner.list_bounds(index)
                }
            },
            | _ => None,
        }
    }

    /// Returns the flattened element column of a list, looking through nulls.
    #[must_use]
    pub fn list_values(&self) -> Option<&Self> {
        match self {
            | Self::List { values, .. } => Some(values),
            | Self::Nullable(inner, _) => inner.list_values(),
            | _ => None,
        }
    }
}

/// Encodes a byte slice as lowercase hexadecimal.
pub fn encode_hex(bytes: &[u8]) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let mut out = String::with_capacity(bytes.len().saturating_mul(2));
    for &byte in bytes {
        out.push(char::from(HEX[usize::from(byte >> 4)]));
        out.push(char::from(HEX[usize::from(byte & 0x0f)]));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numeric_view_and_nulls() {
        let column = Column::Nullable(
            Box::new(Column::Int(vec![1, 2, 3])),
            vec![false, true, false],
        );
        assert_eq!(column.len(), 3);
        assert!(!column.is_empty());
        assert_eq!(column.dtype(), DType::Int);
        assert_eq!(column.name_of(), "int");
        assert!(column.is_nullable());
        assert_eq!(column.get_f64(0), Some(1.0));
        assert_eq!(column.get_f64(1), None);
        assert!(column.is_null(1));
        assert!(!column.is_null(0));
        assert!(column.is_null(99));
        assert_eq!(column.null_mask(), Some(&[false, true, false][..]));
    }

    #[test]
    fn bool_string_and_bytes_views() {
        let booleans = Column::Bool(vec![true, false]);
        assert_eq!(booleans.get_f64(0), Some(1.0));
        assert_eq!(booleans.get_f64(1), Some(0.0));
        let strings = Column::Str(vec!["a".to_string()]);
        assert_eq!(strings.get_f64(0), None);
        assert_eq!(strings.name_of(), "str");
        let bytes = Column::Bytes(vec![vec![0x00, 0xff], vec![0xab]]);
        assert_eq!(bytes.len(), 2);
        assert_eq!(bytes.dtype(), DType::Bytes);
        assert_eq!(bytes.name_of(), "bytes");
        assert_eq!(bytes.get_f64(0), None);
        assert_eq!(encode_hex(&[0x00, 0x0a, 0xff]), "000aff");
        assert!(Column::Int(Vec::new()).is_empty());
    }

    #[test]
    fn nullable_length_must_match() {
        assert!(Column::new_nullable(Column::Bool(vec![true]), Vec::new()).is_err());
        assert!(Column::new_nullable(Column::Bool(vec![true]), vec![false]).is_ok());
    }

    #[test]
    fn lists_expose_rows_and_slice() {
        let list =
            Column::new_list(Column::Int(vec![1, 2, 3, 4]), vec![0, 2, 2, 4]).expect("list column");
        assert_eq!(list.len(), 3);
        assert_eq!(list.dtype(), DType::List);
        assert_eq!(list.name_of(), "list");
        assert_eq!(list.list_bounds(0), Some((0, 2)));
        assert_eq!(list.list_bounds(1), Some((2, 2)));
        assert_eq!(list.list_bounds(2), Some((2, 4)));
        assert_eq!(list.list_bounds(3), None);
        assert_eq!(list.list_row(0), Some(Column::Int(vec![1, 2])));
        assert_eq!(list.list_row(1), Some(Column::Int(Vec::new())));
        assert_eq!(list.get_f64(0), None);

        let tail = list.slice(1, 3);
        assert_eq!(tail.len(), 2);
        assert_eq!(tail.list_bounds(0), Some((0, 0)));
        assert_eq!(tail.list_bounds(1), Some((0, 2)));
        assert_eq!(tail.list_row(1), Some(Column::Int(vec![3, 4])));

        let nested = Column::Nullable(Box::new(list), vec![false, true, false]);
        assert_eq!(nested.len(), 3);
        assert_eq!(nested.list_bounds(1), None);
        assert_eq!(nested.list_bounds(0), Some((0, 2)));
        assert_eq!(nested.slice(0, 2).len(), 2);
    }

    #[test]
    fn list_offsets_are_validated() {
        assert!(Column::new_list(Column::Int(vec![1]), vec![1, 1]).is_err());
        assert!(Column::new_list(Column::Int(vec![1]), vec![0, 0]).is_err());
        assert!(Column::new_list(Column::Int(vec![1]), vec![0, 2]).is_err());
        assert!(Column::new_list(Column::Int(vec![1, 2]), vec![0, 2, 1]).is_err());
        assert!(Column::new_list(Column::Int(Vec::new()), Vec::new()).is_ok());
        assert_eq!(
            Column::new_list(Column::Int(Vec::new()), Vec::new())
                .expect("empty list")
                .len(),
            0
        );
    }
}
