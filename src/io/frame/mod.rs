//! Typed tabular data: columns, data frames, typed CSV and JSON.
//!
//! Unlike the numeric array layer in [`crate::io`], this module keeps a
//! logical type per column: `Int`, `Float`, `Bool`, `Str`, `Bytes` or `List`,
//! optionally wrapped in a null mask. A [`DataFrame`] is a list of equal-length
//! named columns with selection, filtering, sorting, grouping, joining and
//! matrix conversion. A list column holds a flattened element column plus row
//! offsets, so a null list, an empty list and a null element stay distinct.
//!
//! CSV input follows RFC 4180 quoting and infers the column type from its
//! non-missing values (preserving `Nullable(Int)` by default, or widening to
//! `Nullable(Float)` when [`CsvOptions::preserve_nullable_int`] is `false`);
//! JSON is read as either an array of row objects or an object of equal-length
//! arrays, and written with exact column order preserved.
//!
//! # Limitations
//!
//! - Types are limited to 64-bit integers, doubles, booleans, UTF-8 strings,
//!   raw byte sequences and one level of list nesting: no dates, decimals,
//!   categoricals or records.
//! - Nulls are a single mask per column; there is no per-value type union.
//! - CSV inference is per column and all-or-nothing: one non-numeric token
//!   demotes the whole column to `Str`.
//! - Joins and groupings are in-memory hash operations with no spill to disk,
//!   and `group_by_agg` produces one output column per call.
//!
//! [`DataFrame`]: crate::io::frame::DataFrame
//! [`CsvOptions::preserve_nullable_int`]: crate::io::frame::CsvOptions::preserve_nullable_int

mod column;
mod csv;
mod dataframe;
mod json;

pub use column::Column;
pub use column::DType;
pub use csv::CsvOptions;
pub use csv::read_csv_file;
pub use csv::read_csv_str;
pub use csv::write_csv;
pub use csv::write_csv_file;
pub use dataframe::Agg;
pub use dataframe::DataFrame;
pub use json::read_json_file;
pub use json::read_json_str;
pub use json::write_json;
pub use json::write_json_file;
