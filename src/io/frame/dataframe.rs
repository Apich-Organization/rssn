//! Data frames: named, equal-length columns with selection and aggregation.

use std::cmp::Ordering;
use std::collections::HashMap;

use super::column::Column;
use super::column::DType;
use super::column::encode_hex;

/// An aggregation applied by [`DataFrame::group_by_agg`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Agg {
    /// Count of non-missing values.
    Count,
    /// Sum of numeric values.
    Sum,
    /// Arithmetic mean of numeric values.
    Mean,
    /// Smallest value.
    Min,
    /// Largest value.
    Max,
    /// First non-missing value.
    First,
}

impl Agg {
    /// Returns the label used in the aggregated column name.
    fn label(self) -> &'static str {
        match self {
            | Self::Count => "count",
            | Self::Sum => "sum",
            | Self::Mean => "mean",
            | Self::Min => "min",
            | Self::Max => "max",
            | Self::First => "first",
        }
    }
}

/// A collection of named, equal-length columns.
#[derive(Debug, Clone, PartialEq)]
pub struct DataFrame {
    columns: Vec<(String, Column)>,
}

/// A single aggregated value, before it is stored back into a column.
#[derive(Debug, Clone)]
enum Scalar {
    Int(i64),
    Float(f64),
    Bool(bool),
    Str(String),
    Bytes(Vec<u8>),
    List(Vec<Scalar>),
}

impl DataFrame {
    /// Builds a frame from named columns.
    ///
    /// # Errors
    ///
    /// Returns an error when names are empty or duplicated, or when the
    /// columns do not all have the same length.
    pub fn new(columns: Vec<(String, Column)>) -> Result<Self, String> {
        let mut seen: Vec<&str> = Vec::with_capacity(columns.len());
        let mut length: Option<usize> = None;
        for (name, column) in &columns {
            if name.is_empty() {
                return Err("column names must not be empty".to_string());
            }
            if seen.contains(&name.as_str()) {
                return Err(format!("duplicate column name `{name}`"));
            }
            seen.push(name);
            match length {
                | None => length = Some(column.len()),
                | Some(expected) if expected != column.len() => {
                    return Err(format!(
                        "column `{name}` has {} values but the frame has {expected} rows",
                        column.len()
                    ));
                },
                | Some(_) => {},
            }
        }
        Ok(Self { columns })
    }

    /// Returns the number of rows.
    #[must_use]
    pub fn n_rows(&self) -> usize {
        self.columns.first().map_or(0, |(_, column)| column.len())
    }

    /// Returns the number of columns.
    #[must_use]
    pub fn n_cols(&self) -> usize {
        self.columns.len()
    }

    /// Returns the column names in storage order.
    #[must_use]
    pub fn column_names(&self) -> Vec<&str> {
        self.columns.iter().map(|(name, _)| name.as_str()).collect()
    }

    /// Looks up a column by name.
    #[must_use]
    pub fn column(
        &self,
        name: &str,
    ) -> Option<&Column> {
        self.columns
            .iter()
            .find(|(existing, _)| existing == name)
            .map(|(_, column)| column)
    }

    /// Returns the column at `index`, used by the CSV writer.
    pub(crate) fn column_at(
        &self,
        index: usize,
    ) -> Option<&Column> {
        self.columns.get(index).map(|(_, column)| column)
    }

    /// Returns a frame with the requested columns, in the requested order.
    ///
    /// # Errors
    ///
    /// Returns an error when a name is unknown or repeated.
    pub fn select(
        &self,
        names: &[&str],
    ) -> Result<Self, String> {
        let mut selected = Vec::with_capacity(names.len());
        for name in names {
            let column = self
                .column(name)
                .ok_or_else(|| format!("unknown column `{name}`"))?;
            selected.push(((*name).to_string(), column.clone()));
        }
        Self::new(selected)
    }

    /// Returns the first `n` rows (or fewer when the frame is shorter).
    #[must_use]
    pub fn head(
        &self,
        n: usize,
    ) -> Self {
        let end = n.min(self.n_rows());
        self.take_rows(&(0..end).collect::<Vec<_>>())
    }

    /// Returns `len` rows starting at `start`.
    ///
    /// # Errors
    ///
    /// Returns an error when the range does not fit inside the frame.
    pub fn slice(
        &self,
        start: usize,
        len: usize,
    ) -> Result<Self, String> {
        let end = start
            .checked_add(len)
            .ok_or_else(|| "slice range overflows".to_string())?;
        if end > self.n_rows() {
            return Err(format!(
                "slice {start}..{end} is out of bounds for {} rows",
                self.n_rows()
            ));
        }
        Ok(self.take_rows(&(start..end).collect::<Vec<_>>()))
    }

    /// Keeps the rows selected by `mask`.
    ///
    /// # Errors
    ///
    /// Returns an error when the mask length differs from the row count.
    pub fn filter(
        &self,
        mask: &[bool],
    ) -> Result<Self, String> {
        if mask.len() != self.n_rows() {
            return Err(format!(
                "filter mask has {} entries but the frame has {} rows",
                mask.len(),
                self.n_rows()
            ));
        }
        let rows: Vec<usize> = mask
            .iter()
            .enumerate()
            .filter_map(|(index, keep)| if *keep { Some(index) } else { None })
            .collect();
        Ok(self.take_rows(&rows))
    }

    /// Returns the frame sorted by one column; missing values sort last.
    ///
    /// # Errors
    ///
    /// Returns an error when the column does not exist.
    pub fn sort_by(
        &self,
        name: &str,
        descending: bool,
    ) -> Result<Self, String> {
        let column = self
            .column(name)
            .ok_or_else(|| format!("unknown column `{name}`"))?;
        let mut rows: Vec<usize> = (0..self.n_rows()).collect();
        rows.sort_by(|&left, &right| compare_for_sort(column, left, right, descending));
        Ok(self.take_rows(&rows))
    }

    /// Returns the frame with every row that contains a missing value removed.
    ///
    /// The surviving columns lose their null layer.
    #[must_use]
    pub fn drop_nulls(&self) -> Self {
        let rows: Vec<usize> = (0..self.n_rows())
            .filter(|&row| !self.columns.iter().any(|(_, column)| column.is_null(row)))
            .collect();
        let columns = self
            .columns
            .iter()
            .map(|(name, column)| {
                let compacted = match column {
                    | Column::Nullable(inner, _) => take(inner, &rows),
                    | other => take(other, &rows),
                };
                (name.clone(), compacted)
            })
            .collect();
        Self { columns }
    }

    /// Returns a frame with `column` stored under `name`.
    ///
    /// An existing column with the same name is replaced. Adding the first
    /// column to an empty frame fixes the row count.
    ///
    /// # Errors
    ///
    /// Returns an error when the column length does not match the frame.
    pub fn with_column(
        &self,
        name: String,
        column: Column,
    ) -> Result<Self, String> {
        if name.is_empty() {
            return Err("column names must not be empty".to_string());
        }
        if !self.columns.is_empty() && column.len() != self.n_rows() {
            return Err(format!(
                "column `{name}` has {} values but the frame has {} rows",
                column.len(),
                self.n_rows()
            ));
        }
        let mut columns = self.columns.clone();
        if let Some(slot) = columns.iter_mut().find(|(existing, _)| *existing == name) {
            slot.1 = column;
        } else {
            columns.push((name, column));
        }
        Self::new(columns)
    }

    /// Stacks the rows of `other` below this frame.
    ///
    /// # Errors
    ///
    /// Returns an error when the schemas (names, types and nullability) do
    /// not match.
    pub fn concat_rows(
        &self,
        other: &Self,
    ) -> Result<Self, String> {
        if self.n_cols() != other.n_cols() {
            return Err(format!(
                "cannot concatenate {} columns with {} columns",
                self.n_cols(),
                other.n_cols()
            ));
        }
        for ((left_name, left), (right_name, right)) in self.columns.iter().zip(&other.columns) {
            if left_name != right_name {
                return Err(format!(
                    "cannot concatenate column `{left_name}` with `{right_name}`"
                ));
            }
            if left.dtype() != right.dtype() {
                return Err(format!("column `{left_name}` has mismatched types"));
            }
            if left.is_nullable() != right.is_nullable() {
                return Err(format!("column `{left_name}` has mismatched nullability"));
            }
        }
        let mut columns = Vec::with_capacity(self.n_cols());
        for ((name, left), (_, right)) in self.columns.iter().zip(&other.columns) {
            columns.push((name.clone(), concat_columns(left, right)?));
        }
        Self::new(columns)
    }

    /// Groups rows by the key columns and aggregates `value`.
    ///
    /// The result holds the key columns followed by a column named
    /// `"{agg}_{value}"`.
    ///
    /// # Errors
    ///
    /// Returns an error when a key or value column is unknown, when the
    /// aggregation is numeric but the value column is not, or when integer
    /// sums overflow.
    pub fn group_by_agg(
        &self,
        keys: &[&str],
        value: &str,
        agg: Agg,
    ) -> Result<Self, String> {
        let value_column = self
            .column(value)
            .ok_or_else(|| format!("unknown column `{value}`"))?;
        let mut key_columns: Vec<(&str, &Column)> = Vec::with_capacity(keys.len());
        for key in keys {
            let column = self
                .column(key)
                .ok_or_else(|| format!("unknown column `{key}`"))?;
            key_columns.push((key, column));
        }
        let mut groups: Vec<Vec<usize>> = Vec::new();
        let mut index: HashMap<String, usize> = HashMap::new();
        for row in 0..self.n_rows() {
            let token = group_key(&key_columns, row);
            if let Some(group) = index.get(&token) {
                if let Some(rows) = groups.get_mut(*group) {
                    rows.push(row);
                }
            } else {
                index.insert(token, groups.len());
                groups.push(vec![row]);
            }
        }
        let first_rows: Vec<usize> = groups
            .iter()
            .map(|rows| rows.first().copied().unwrap_or(0))
            .collect();
        let mut columns: Vec<(String, Column)> = Vec::with_capacity(keys.len() + 1);
        for (key, column) in &key_columns {
            columns.push(((*key).to_string(), take(column, &first_rows)));
        }
        columns.push((
            format!("{}_{}", agg.label(), value),
            aggregate(agg, value_column, &groups)?,
        ));
        Self::new(columns)
    }

    /// Joins with `other` on equal values of the `on` column.
    ///
    /// Rows whose key is missing never match. The output keeps the left
    /// columns, then the right columns except the key; clashing right names
    /// gain a `_right` suffix.
    ///
    /// # Errors
    ///
    /// Returns an error when the key is missing from either frame or the two
    /// key columns have different types.
    pub fn inner_join(
        &self,
        other: &Self,
        on: &str,
    ) -> Result<Self, String> {
        let left_key = self
            .column(on)
            .ok_or_else(|| format!("unknown column `{on}` in the left frame"))?;
        let right_key = other
            .column(on)
            .ok_or_else(|| format!("unknown column `{on}` in the right frame"))?;
        if left_key.dtype() != right_key.dtype() {
            return Err(format!("join column `{on}` has different types"));
        }
        let mut lookup: HashMap<String, Vec<usize>> = HashMap::new();
        for row in 0..other.n_rows() {
            if right_key.is_null(row) {
                continue;
            }
            lookup
                .entry(value_token(right_key, row))
                .or_default()
                .push(row);
        }
        let mut left_rows = Vec::new();
        let mut right_rows = Vec::new();
        for row in 0..self.n_rows() {
            if left_key.is_null(row) {
                continue;
            }
            if let Some(matches) = lookup.get(&value_token(left_key, row)) {
                for &matched in matches {
                    left_rows.push(row);
                    right_rows.push(matched);
                }
            }
        }
        let mut columns: Vec<(String, Column)> = Vec::with_capacity(self.n_cols() + other.n_cols());
        for (name, column) in &self.columns {
            columns.push((name.clone(), take(column, &left_rows)));
        }
        let mut used: Vec<String> = columns.iter().map(|(name, _)| name.clone()).collect();
        for (name, column) in &other.columns {
            if name == on {
                continue;
            }
            let mut output_name = name.clone();
            if used.contains(&output_name) {
                output_name = format!("{name}_right");
            }
            let mut suffix = 2usize;
            while used.contains(&output_name) {
                output_name = format!("{name}_right{suffix}");
                suffix = suffix.saturating_add(1);
            }
            used.push(output_name.clone());
            columns.push((output_name, take(column, &right_rows)));
        }
        Self::new(columns)
    }

    /// Returns the named numeric columns as a row-major matrix.
    ///
    /// # Errors
    ///
    /// Returns an error when a name is unknown, the column is not numeric, or
    /// a value is missing.
    pub fn to_matrix(
        &self,
        names: &[&str],
    ) -> Result<Vec<Vec<f64>>, String> {
        let mut selected: Vec<(&str, &Column)> = Vec::with_capacity(names.len());
        for name in names {
            let column = self
                .column(name)
                .ok_or_else(|| format!("unknown column `{name}`"))?;
            if !matches!(column.dtype(), DType::Int | DType::Float) {
                return Err(format!("column `{name}` is not numeric"));
            }
            selected.push((name, column));
        }
        let mut rows = Vec::with_capacity(self.n_rows());
        for row in 0..self.n_rows() {
            let mut values = Vec::with_capacity(selected.len());
            for (name, column) in &selected {
                let value = column
                    .get_f64(row)
                    .ok_or_else(|| format!("column `{name}` has a missing value at row {row}"))?;
                values.push(value);
            }
            rows.push(values);
        }
        Ok(rows)
    }

    /// Materialises the rows at `rows` as a new frame.
    fn take_rows(
        &self,
        rows: &[usize],
    ) -> Self {
        Self {
            columns: self
                .columns
                .iter()
                .map(|(name, column)| (name.clone(), take(column, rows)))
                .collect(),
        }
    }
}

/// Copies the selected rows out of a column, preserving null masks.
fn take(
    column: &Column,
    rows: &[usize],
) -> Column {
    match column {
        | Column::Int(values) => {
            Column::Int(
                rows.iter()
                    .map(|&row| values.get(row).copied().unwrap_or(0))
                    .collect(),
            )
        },
        | Column::Float(values) => {
            Column::Float(
                rows.iter()
                    .map(|&row| values.get(row).copied().unwrap_or(0.0))
                    .collect(),
            )
        },
        | Column::Bool(values) => {
            Column::Bool(
                rows.iter()
                    .map(|&row| values.get(row).copied().unwrap_or(false))
                    .collect(),
            )
        },
        | Column::Str(values) => {
            Column::Str(
                rows.iter()
                    .map(|&row| values.get(row).cloned().unwrap_or_default())
                    .collect(),
            )
        },
        | Column::Bytes(values) => {
            Column::Bytes(
                rows.iter()
                    .map(|&row| values.get(row).cloned().unwrap_or_default())
                    .collect(),
            )
        },
        | Column::List { values, offsets } => {
            let mut flat = Vec::new();
            let mut taken = Vec::with_capacity(rows.len() + 1);
            taken.push(0usize);
            for &row in rows {
                let start = offsets.get(row).copied().unwrap_or(0);
                let end = offsets
                    .get(row + 1)
                    .copied()
                    .unwrap_or(start)
                    .min(values.len());
                let start = start.min(end);
                flat.extend(start..end);
                taken.push(flat.len());
            }
            Column::List {
                values: Box::new(take(values, &flat)),
                offsets: taken,
            }
        },
        | Column::Nullable(inner, mask) => {
            let taken_mask = rows
                .iter()
                .map(|&row| mask.get(row).copied().unwrap_or(true))
                .collect();
            Column::Nullable(Box::new(take(inner, rows)), taken_mask)
        },
    }
}

/// Concatenates two columns of the same variant.
fn concat_columns(
    left: &Column,
    right: &Column,
) -> Result<Column, String> {
    match (left, right) {
        | (Column::Int(first), Column::Int(second)) => {
            let mut values = first.clone();
            values.extend_from_slice(second);
            Ok(Column::Int(values))
        },
        | (Column::Float(first), Column::Float(second)) => {
            let mut values = first.clone();
            values.extend_from_slice(second);
            Ok(Column::Float(values))
        },
        | (Column::Bool(first), Column::Bool(second)) => {
            let mut values = first.clone();
            values.extend_from_slice(second);
            Ok(Column::Bool(values))
        },
        | (Column::Str(first), Column::Str(second)) => {
            let mut values = first.clone();
            values.extend_from_slice(second);
            Ok(Column::Str(values))
        },
        | (Column::Bytes(first), Column::Bytes(second)) => {
            let mut values = first.clone();
            values.extend_from_slice(second);
            Ok(Column::Bytes(values))
        },
        | (
            Column::List {
                values: first,
                offsets: first_offsets,
            },
            Column::List {
                values: second,
                offsets: second_offsets,
            },
        ) => {
            let values = concat_columns(first, second)?;
            let base = first.len();
            let mut offsets = if first_offsets.is_empty() {
                vec![0]
            } else {
                first_offsets.clone()
            };
            offsets.extend(second_offsets.iter().skip(1).map(|offset| offset + base));
            Ok(Column::List {
                values: Box::new(values),
                offsets,
            })
        },
        | (Column::Nullable(first, first_mask), Column::Nullable(second, second_mask)) => {
            let inner = concat_columns(first, second)?;
            let mut mask = first_mask.clone();
            mask.extend_from_slice(second_mask);
            Ok(Column::Nullable(Box::new(inner), mask))
        },
        | _ => Err("column types do not match".to_string()),
    }
}

/// Compares two slots of a column, placing missing values last.
fn compare_for_sort(
    column: &Column,
    left: usize,
    right: usize,
    descending: bool,
) -> Ordering {
    match (column.is_null(left), column.is_null(right)) {
        | (true, true) => Ordering::Equal,
        | (true, false) => Ordering::Greater,
        | (false, true) => Ordering::Less,
        | (false, false) => {
            let order = compare_values(column, left, right);
            if descending {
                order.reverse()
            } else {
                order
            }
        },
    }
}

/// Compares two non-missing slots.
fn compare_values(
    column: &Column,
    left: usize,
    right: usize,
) -> Ordering {
    match column {
        | Column::Int(values) => compare_pairs(values.get(left), values.get(right)),
        | Column::Float(values) => {
            match (values.get(left), values.get(right)) {
                | (Some(first), Some(second)) => first.total_cmp(second),
                | _ => Ordering::Equal,
            }
        },
        | Column::Bool(values) => compare_pairs(values.get(left), values.get(right)),
        | Column::Str(values) => compare_pairs(values.get(left), values.get(right)),
        | Column::Bytes(values) => compare_pairs(values.get(left), values.get(right)),
        | Column::List { values, offsets } => {
            let left_start = offsets.get(left).copied().unwrap_or(0);
            let left_end = offsets.get(left + 1).copied().unwrap_or(left_start);
            let right_start = offsets.get(right).copied().unwrap_or(0);
            let right_end = offsets.get(right + 1).copied().unwrap_or(right_start);
            let mut first = left_start;
            let mut second = right_start;
            while first < left_end && second < right_end {
                let order = compare_values(values, first, second);
                if order != Ordering::Equal {
                    return order;
                }
                first += 1;
                second += 1;
            }
            (left_end.saturating_sub(left_start)).cmp(&right_end.saturating_sub(right_start))
        },
        | Column::Nullable(inner, _) => compare_values(inner, left, right),
    }
}

/// Compares two optional values, treating a missing side as equal.
fn compare_pairs<T: Ord>(
    left: Option<&T>,
    right: Option<&T>,
) -> Ordering {
    match (left, right) {
        | (Some(first), Some(second)) => first.cmp(second),
        | _ => Ordering::Equal,
    }
}

/// Encodes the key columns of one row into a collision-resistant token.
fn group_key(
    columns: &[(&str, &Column)],
    row: usize,
) -> String {
    let mut key = String::new();
    for (_, column) in columns {
        let token = if column.is_null(row) {
            "n".to_string()
        } else {
            value_token(column, row)
        };
        key.push_str(&token.len().to_string());
        key.push(':');
        key.push_str(&token);
    }
    key
}

/// Encodes one non-missing value with a type tag.
fn value_token(
    column: &Column,
    row: usize,
) -> String {
    match column {
        | Column::Int(values) => format!("i{}", values.get(row).copied().unwrap_or(0)),
        | Column::Float(values) => {
            format!(
                "f{:016x}",
                values.get(row).copied().unwrap_or(0.0).to_bits()
            )
        },
        | Column::Bool(values) => format!("b{}", values.get(row).copied().unwrap_or(false)),
        | Column::Str(values) => {
            format!("s{}", values.get(row).map_or("", String::as_str))
        },
        | Column::Bytes(values) => {
            format!(
                "x{}",
                values.get(row).map_or_else(String::new, |v| encode_hex(v))
            )
        },
        | Column::List { values, offsets } => {
            let start = offsets.get(row).copied().unwrap_or(0);
            let end = offsets.get(row + 1).copied().unwrap_or(start);
            let mut token = String::from("l[");
            for index in start..end {
                token.push_str(&value_token(values, index));
                token.push(',');
            }
            token.push(']');
            token
        },
        | Column::Nullable(inner, mask) => {
            if mask.get(row).copied().unwrap_or(true) {
                "n".to_string()
            } else {
                value_token(inner, row)
            }
        },
    }
}

/// Aggregates one group into a single scalar.
fn aggregate_group(
    agg: Agg,
    column: &Column,
    rows: &[usize],
) -> Result<Option<Scalar>, String> {
    match agg {
        | Agg::Count => {
            Ok(Some(Scalar::Int(
                rows.iter().filter(|&&row| !column.is_null(row)).count() as i64,
            )))
        },
        | Agg::Sum => {
            match column.dtype() {
                | DType::Int => {
                    let mut total = 0i64;
                    let mut seen = false;
                    for &row in rows {
                        if let Some(value) = int_at(column, row) {
                            total = total
                                .checked_add(value)
                                .ok_or_else(|| "integer overflow while summing".to_string())?;
                            seen = true;
                        }
                    }
                    Ok(if seen {
                        Some(Scalar::Int(total))
                    } else {
                        None
                    })
                },
                | DType::Float => {
                    let mut total = 0.0f64;
                    let mut seen = false;
                    for &row in rows {
                        if let Some(value) = column.get_f64(row) {
                            total += value;
                            seen = true;
                        }
                    }
                    Ok(if seen {
                        Some(Scalar::Float(total))
                    } else {
                        None
                    })
                },
                | _ => Err("`sum` requires a numeric column".to_string()),
            }
        },
        | Agg::Mean => {
            let mut total = 0.0f64;
            let mut count = 0usize;
            for &row in rows {
                if let Some(value) = column.get_f64(row) {
                    total += value;
                    count += 1;
                }
            }
            Ok(if count == 0 {
                None
            } else {
                Some(Scalar::Float(total / count as f64))
            })
        },
        | Agg::Min | Agg::Max => {
            let mut best: Option<usize> = None;
            for &row in rows {
                if column.is_null(row) {
                    continue;
                }
                best = match best {
                    | None => Some(row),
                    | Some(current) => {
                        let order = compare_values(column, row, current);
                        let replace = if agg == Agg::Min {
                            order == Ordering::Less
                        } else {
                            order == Ordering::Greater
                        };
                        if replace {
                            Some(row)
                        } else {
                            Some(current)
                        }
                    },
                };
            }
            match best {
                | None => Ok(None),
                | Some(row) => Ok(scalar_at(column, row)),
            }
        },
        | Agg::First => {
            for &row in rows {
                if !column.is_null(row) {
                    return Ok(scalar_at(column, row));
                }
            }
            Ok(None)
        },
    }
}

/// Aggregates every group into a column.
fn aggregate(
    agg: Agg,
    column: &Column,
    groups: &[Vec<usize>],
) -> Result<Column, String> {
    let source = column.dtype();
    if matches!(agg, Agg::Sum | Agg::Mean) && !matches!(source, DType::Int | DType::Float) {
        return Err(format!(
            "`{}` requires a numeric column but `{}` was given",
            agg.label(),
            column.name_of()
        ));
    }
    let target = match agg {
        | Agg::Count => DType::Int,
        | Agg::Sum if source == DType::Int => DType::Int,
        | Agg::Sum | Agg::Mean => DType::Float,
        | Agg::Min | Agg::Max | Agg::First => source,
    };
    let mut scalars = Vec::with_capacity(groups.len());
    for rows in groups {
        scalars.push(aggregate_group(agg, column, rows)?);
    }
    materialize(target, &scalars)
}

/// Reads a non-missing integer slot, looking through nullable wrappers.
fn int_at(
    column: &Column,
    row: usize,
) -> Option<i64> {
    match column {
        | Column::Int(values) => values.get(row).copied(),
        | Column::Nullable(inner, mask) => {
            if mask.get(row).copied().unwrap_or(true) {
                None
            } else {
                int_at(inner, row)
            }
        },
        | _ => None,
    }
}

/// Reads a non-missing slot as a scalar.
fn scalar_at(
    column: &Column,
    row: usize,
) -> Option<Scalar> {
    match column {
        | Column::Int(values) => values.get(row).map(|value| Scalar::Int(*value)),
        | Column::Float(values) => values.get(row).map(|value| Scalar::Float(*value)),
        | Column::Bool(values) => values.get(row).map(|value| Scalar::Bool(*value)),
        | Column::Str(values) => values.get(row).map(|value| Scalar::Str(value.clone())),
        | Column::Bytes(values) => values.get(row).map(|value| Scalar::Bytes(value.clone())),
        | Column::List { values, offsets } => {
            let start = offsets.get(row).copied()?;
            let end = offsets.get(row + 1).copied()?;
            let mut items = Vec::with_capacity(end.saturating_sub(start));
            for index in start..end {
                items.push(scalar_at(values, index)?);
            }
            Some(Scalar::List(items))
        },
        | Column::Nullable(inner, mask) => {
            if mask.get(row).copied().unwrap_or(true) {
                None
            } else {
                scalar_at(inner, row)
            }
        },
    }
}

/// Turns aggregated scalars into a typed column, adding a mask when needed.
fn materialize(
    dtype: DType,
    scalars: &[Option<Scalar>],
) -> Result<Column, String> {
    if dtype == DType::List {
        return materialize_lists(scalars);
    }
    let any_missing = scalars.iter().any(Option::is_none);
    let mut mask = Vec::with_capacity(scalars.len());
    let mut ints = Vec::with_capacity(scalars.len());
    let mut floats = Vec::with_capacity(scalars.len());
    let mut booleans = Vec::with_capacity(scalars.len());
    let mut strings = Vec::with_capacity(scalars.len());
    let mut bytes_vec = Vec::with_capacity(scalars.len());
    for scalar in scalars {
        mask.push(scalar.is_none());
        match dtype {
            | DType::Int => {
                ints.push(match scalar {
                    | Some(Scalar::Int(value)) => *value,
                    | Some(Scalar::Float(value)) => *value as i64,
                    | _ => 0,
                })
            },
            | DType::Float => {
                floats.push(match scalar {
                    | Some(Scalar::Float(value)) => *value,
                    | Some(Scalar::Int(value)) => *value as f64,
                    | _ => 0.0,
                })
            },
            | DType::Bool => {
                booleans.push(match scalar {
                    | Some(Scalar::Bool(value)) => *value,
                    | _ => false,
                })
            },
            | DType::Str => {
                strings.push(match scalar {
                    | Some(Scalar::Str(value)) => value.clone(),
                    | _ => String::new(),
                })
            },
            | DType::Bytes => {
                bytes_vec.push(match scalar {
                    | Some(Scalar::Bytes(value)) => value.clone(),
                    | _ => Vec::new(),
                })
            },
            | DType::List => return materialize_lists(scalars),
        }
    }
    let base = match dtype {
        | DType::Int => Column::Int(ints),
        | DType::Float => Column::Float(floats),
        | DType::Bool => Column::Bool(booleans),
        | DType::Str => Column::Str(strings),
        | DType::Bytes => Column::Bytes(bytes_vec),
        | DType::List => return materialize_lists(scalars),
    };
    if any_missing {
        Ok(Column::Nullable(Box::new(base), mask))
    } else {
        Ok(base)
    }
}

/// Builds a list column out of aggregated scalar lists.
fn materialize_lists(scalars: &[Option<Scalar>]) -> Result<Column, String> {
    let mut offsets = Vec::with_capacity(scalars.len() + 1);
    offsets.push(0usize);
    let mut flat: Vec<Option<Scalar>> = Vec::new();
    let mut mask = Vec::with_capacity(scalars.len());
    for scalar in scalars {
        match scalar {
            | Some(Scalar::List(items)) => {
                mask.push(false);
                flat.extend(items.iter().cloned().map(Some));
            },
            | Some(_) => {
                return Err("a list column cannot be built from non-list values".to_string());
            },
            | None => mask.push(true),
        }
        offsets.push(flat.len());
    }
    let element = flat
        .iter()
        .flatten()
        .next()
        .map_or(DType::Int, scalar_dtype);
    let values = materialize(element, &flat)?;
    let list = Column::List {
        values: Box::new(values),
        offsets,
    };
    if mask.iter().any(|&missing| missing) {
        Ok(Column::Nullable(Box::new(list), mask))
    } else {
        Ok(list)
    }
}

/// The logical type a scalar is stored as.
fn scalar_dtype(scalar: &Scalar) -> DType {
    match scalar {
        | Scalar::Int(_) => DType::Int,
        | Scalar::Float(_) => DType::Float,
        | Scalar::Bool(_) => DType::Bool,
        | Scalar::Str(_) => DType::Str,
        | Scalar::Bytes(_) => DType::Bytes,
        | Scalar::List(_) => DType::List,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> DataFrame {
        DataFrame::new(vec![
            ("id".to_string(), Column::Int(vec![3, 1, 2, 4])),
            (
                "name".to_string(),
                Column::Str(vec![
                    "c".to_string(),
                    "a".to_string(),
                    "b".to_string(),
                    "d".to_string(),
                ]),
            ),
            (
                "score".to_string(),
                Column::Float(vec![30.0, 10.0, 20.0, 40.0]),
            ),
        ])
        .expect("sample frame")
    }

    #[test]
    fn filter_sort_select_chain() {
        let frame = sample();
        let filtered = frame
            .filter(&[true, true, false, true])
            .expect("filter works");
        assert_eq!(filtered.n_rows(), 3);
        let sorted = filtered.sort_by("score", false).expect("sort works");
        assert_eq!(sorted.column("id"), Some(&Column::Int(vec![1, 3, 4])));
        let selected = sorted.select(&["name", "id"]).expect("select works");
        assert_eq!(selected.column_names(), vec!["name", "id"]);
        assert_eq!(
            selected.column("name"),
            Some(&Column::Str(vec![
                "a".to_string(),
                "c".to_string(),
                "d".to_string()
            ]))
        );
    }

    #[test]
    fn sort_puts_nulls_last_in_both_directions() {
        let frame = DataFrame::new(vec![(
            "value".to_string(),
            Column::Nullable(
                Box::new(Column::Int(vec![2, 0, 1])),
                vec![false, true, false],
            ),
        )])
        .expect("frame");
        let ascending = frame.sort_by("value", false).expect("sort");
        let ascending_mask = ascending
            .column("value")
            .and_then(Column::null_mask)
            .expect("mask");
        assert_eq!(ascending_mask, &[false, false, true]);
        assert_eq!(
            ascending
                .column("value")
                .and_then(|column| column.get_f64(0)),
            Some(1.0)
        );
        assert_eq!(
            ascending
                .column("value")
                .and_then(|column| column.get_f64(1)),
            Some(2.0)
        );
        let descending = frame.sort_by("value", true).expect("sort");
        let mask = descending
            .column("value")
            .and_then(Column::null_mask)
            .expect("mask");
        assert_eq!(mask, &[false, false, true]);
    }

    #[test]
    fn group_by_sum_and_count() {
        let frame = DataFrame::new(vec![
            (
                "team".to_string(),
                Column::Str(vec![
                    "a".to_string(),
                    "b".to_string(),
                    "a".to_string(),
                    "b".to_string(),
                ]),
            ),
            ("points".to_string(), Column::Int(vec![1, 2, 3, 4])),
        ])
        .expect("frame");
        let summed = frame
            .group_by_agg(&["team"], "points", Agg::Sum)
            .expect("group");
        assert_eq!(summed.column_names(), vec!["team", "sum_points"]);
        assert_eq!(
            summed.column("team"),
            Some(&Column::Str(vec!["a".to_string(), "b".to_string()]))
        );
        assert_eq!(summed.column("sum_points"), Some(&Column::Int(vec![4, 6])));
        let counted = frame
            .group_by_agg(&["team"], "points", Agg::Count)
            .expect("group");
        assert_eq!(
            counted.column("count_points"),
            Some(&Column::Int(vec![2, 2]))
        );
    }

    #[test]
    fn group_by_mean_skips_nulls() {
        let frame = DataFrame::new(vec![
            (
                "key".to_string(),
                Column::Str(vec!["x".to_string(), "x".to_string()]),
            ),
            (
                "value".to_string(),
                Column::Nullable(Box::new(Column::Float(vec![1.0, 3.0])), vec![false, true]),
            ),
        ])
        .expect("frame");
        let grouped = frame
            .group_by_agg(&["key"], "value", Agg::Mean)
            .expect("group");
        assert_eq!(
            grouped.column("mean_value"),
            Some(&Column::Float(vec![1.0]))
        );
    }

    #[test]
    fn inner_join_expands_duplicates() {
        let left = DataFrame::new(vec![
            ("key".to_string(), Column::Int(vec![1, 2])),
            (
                "left".to_string(),
                Column::Str(vec!["l1".to_string(), "l2".to_string()]),
            ),
        ])
        .expect("left");
        let right = DataFrame::new(vec![
            ("key".to_string(), Column::Int(vec![1, 1, 3])),
            ("right".to_string(), Column::Int(vec![10, 11, 12])),
        ])
        .expect("right");
        let joined = left.inner_join(&right, "key").expect("join");
        assert_eq!(joined.n_rows(), 2);
        assert_eq!(joined.column_names(), vec!["key", "left", "right"]);
        assert_eq!(joined.column("right"), Some(&Column::Int(vec![10, 11])));
        assert_eq!(
            joined.column("left"),
            Some(&Column::Str(vec!["l1".to_string(), "l1".to_string()]))
        );
    }

    #[test]
    fn inner_join_renames_clashing_columns() {
        let left = DataFrame::new(vec![
            ("key".to_string(), Column::Int(vec![1])),
            ("value".to_string(), Column::Int(vec![10])),
        ])
        .expect("left");
        let right = DataFrame::new(vec![
            ("key".to_string(), Column::Int(vec![1])),
            ("value".to_string(), Column::Int(vec![20])),
        ])
        .expect("right");
        let joined = left.inner_join(&right, "key").expect("join");
        assert_eq!(joined.column_names(), vec!["key", "value", "value_right"]);
        assert_eq!(joined.column("value_right"), Some(&Column::Int(vec![20])));
    }

    #[test]
    fn concat_rows_requires_matching_schema() {
        let frame = sample();
        let doubled = frame.concat_rows(&frame).expect("concat");
        assert_eq!(doubled.n_rows(), 8);
        let other = DataFrame::new(vec![("id".to_string(), Column::Int(vec![1]))]).expect("frame");
        assert!(frame.concat_rows(&other).is_err());
    }

    #[test]
    fn drop_nulls_removes_rows_and_masks() {
        let frame = DataFrame::new(vec![
            (
                "a".to_string(),
                Column::Nullable(Box::new(Column::Int(vec![1, 2])), vec![false, true]),
            ),
            ("b".to_string(), Column::Bool(vec![true, false])),
        ])
        .expect("frame");
        let dropped = frame.drop_nulls();
        assert_eq!(dropped.n_rows(), 1);
        assert_eq!(dropped.column("a"), Some(&Column::Int(vec![1])));
    }

    #[test]
    fn slice_and_head_are_bounds_checked() {
        let frame = sample();
        assert_eq!(frame.head(2).n_rows(), 2);
        assert_eq!(frame.slice(2, 2).expect("slice").n_rows(), 2);
        assert!(frame.slice(3, 2).is_err());
        assert!(frame.slice(usize::MAX, 1).is_err());
        assert!(frame.filter(&[true]).is_err());
    }

    #[test]
    fn with_column_appends_and_replaces() {
        let frame = sample();
        let extended = frame
            .with_column(
                "extra".to_string(),
                Column::Bool(vec![true, false, true, false]),
            )
            .expect("with_column");
        assert_eq!(extended.n_cols(), 4);
        let replaced = frame
            .with_column("id".to_string(), Column::Int(vec![0, 0, 0, 0]))
            .expect("replace");
        assert_eq!(replaced.column("id"), Some(&Column::Int(vec![0, 0, 0, 0])));
        assert!(
            frame
                .with_column("bad".to_string(), Column::Int(vec![1]))
                .is_err()
        );
    }

    #[test]
    fn to_matrix_requires_numbers() {
        let frame = sample();
        let matrix = frame.to_matrix(&["id", "score"]).expect("matrix");
        assert_eq!(matrix.len(), 4);
        assert_eq!(matrix[0], vec![3.0, 30.0]);
        assert!(frame.to_matrix(&["name"]).is_err());
    }

    #[test]
    fn new_rejects_bad_frames() {
        assert!(DataFrame::new(vec![("a".to_string(), Column::Int(vec![1]))]).is_ok());
        assert!(
            DataFrame::new(vec![
                ("a".to_string(), Column::Int(vec![1])),
                ("a".to_string(), Column::Int(vec![2])),
            ])
            .is_err()
        );
        assert!(
            DataFrame::new(vec![
                ("a".to_string(), Column::Int(vec![1])),
                ("b".to_string(), Column::Int(vec![1, 2])),
            ])
            .is_err()
        );
        assert!(
            DataFrame::new(Vec::new())
                .expect("empty")
                .column_names()
                .is_empty()
        );
    }

    #[test]
    fn bytes_column_sort_group_join_and_concat() {
        let frame = DataFrame::new(vec![
            (
                "k".to_string(),
                Column::Bytes(vec![vec![0xff, 0x01], vec![0x00, 0x02], vec![0xff, 0x01]]),
            ),
            ("v".to_string(), Column::Int(vec![10, 20, 30])),
        ])
        .expect("frame");
        let sorted = frame.sort_by("k", false).expect("sort");
        assert_eq!(
            sorted.column("k"),
            Some(&Column::Bytes(vec![
                vec![0x00, 0x02],
                vec![0xff, 0x01],
                vec![0xff, 0x01],
            ]))
        );
        let grouped = frame.group_by_agg(&["k"], "v", Agg::Sum).expect("group");
        assert_eq!(grouped.n_rows(), 2);
        assert_eq!(grouped.column("sum_v"), Some(&Column::Int(vec![40, 20])));
        let first_bytes = frame.group_by_agg(&["k"], "k", Agg::First).expect("first");
        assert_eq!(
            first_bytes.column("first_k"),
            Some(&Column::Bytes(vec![vec![0xff, 0x01], vec![0x00, 0x02]]))
        );
        let right = DataFrame::new(vec![
            ("k".to_string(), Column::Bytes(vec![vec![0x00, 0x02]])),
            ("tag".to_string(), Column::Str(vec!["hit".to_string()])),
        ])
        .expect("right");
        let joined = frame.inner_join(&right, "k").expect("join");
        assert_eq!(joined.n_rows(), 1);
        assert_eq!(joined.column("v"), Some(&Column::Int(vec![20])));
        let doubled = frame.concat_rows(&frame).expect("concat");
        assert_eq!(doubled.n_rows(), 6);
    }
}
