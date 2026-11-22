//! Typed tabular data: columns, data frames, typed CSV and JSON.
//!
//! Unlike the numeric array layer in [`crate::io`], this module keeps a
//! logical type per column: `Int`, `Float`, `Bool` or `Str`, optionally
//! wrapped in a null mask. A [`DataFrame`] is a list of equal-length named
//! columns with selection, filtering, sorting, grouping, joining and matrix
//! conversion.
//!
//! CSV input follows RFC 4180 quoting and infers the column type from its
//! non-missing values; JSON is read as either an array of row objects or an
//! object of equal-length arrays.
//!
//! # Limitations
//!
//! - Types are limited to 64-bit integers, doubles, booleans and strings: no
//!   dates, decimals, categoricals, lists or nested structures.
//! - Nulls are a single mask per column; there is no per-value type union.
//! - CSV inference is per column and all-or-nothing: one non-numeric token
//!   demotes the whole column to `Str`, and a column of integers with a
//!   missing value is widened to `Float`.
//! - Joins and groupings are in-memory hash operations with no spill to disk,
//!   and `group_by_agg` produces one output column per call.
//! - JSON objects lose their key order on the way out because `serde_json`
//!   stores objects in a sorted map by default.

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
