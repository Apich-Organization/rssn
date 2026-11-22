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
            | Self::Nullable(inner, _) => inner.dtype(),
        }
    }

    /// Returns a numeric view of the value at `index`.
    ///
    /// Integers and booleans (`false` is `0.0`, `true` is `1.0`) are widened
    /// to `f64`; strings have no numeric view. Missing positions and
    /// out-of-range indices return `None`.
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
            | Self::Str(_) => None,
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
    /// `"bool"` or `"str"`).
    #[must_use]
    pub fn name_of(&self) -> &'static str {
        match self.dtype() {
            | DType::Int => "int",
            | DType::Float => "float",
            | DType::Bool => "bool",
            | DType::Str => "str",
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
    fn bool_and_string_views() {
        let booleans = Column::Bool(vec![true, false]);
        assert_eq!(booleans.get_f64(0), Some(1.0));
        assert_eq!(booleans.get_f64(1), Some(0.0));
        let strings = Column::Str(vec!["a".to_string()]);
        assert_eq!(strings.get_f64(0), None);
        assert_eq!(strings.name_of(), "str");
        assert!(Column::Int(Vec::new()).is_empty());
    }

    #[test]
    fn nullable_length_must_match() {
        assert!(Column::new_nullable(Column::Bool(vec![true]), Vec::new()).is_err());
        assert!(Column::new_nullable(Column::Bool(vec![true]), vec![false]).is_ok());
    }
}
