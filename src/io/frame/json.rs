//! JSON reading and writing for typed frames.

use std::collections::HashSet;
use std::path::Path;

use serde_json::Map;
use serde_json::Number;
use serde_json::Value;

use super::column::Column;
use super::column::DType;
use super::column::encode_hex;
use super::dataframe::DataFrame;

/// Reads a frame from JSON text.
///
/// Two layouts are accepted: an array of objects (one object per row) and an
/// object whose values are equal-length arrays (one key per column). Missing
/// keys and JSON `null` become missing values.
///
/// # Errors
///
/// Returns an error for invalid JSON, mismatched array lengths, or values
/// that are neither scalars nor strings.
pub fn read_json_str(text: &str) -> Result<DataFrame, String> {
    let value: Value = serde_json::from_str(text).map_err(|error| error.to_string())?;
    match value {
        | Value::Array(rows) => from_rows(&rows),
        | Value::Object(object) => from_columns(&object),
        | _ => Err("json must be an array of objects or an object of arrays".to_string()),
    }
}

/// Writes a frame as JSON text: an array with one object per row, preserving
/// the exact column order of `df`.
///
/// # Errors
///
/// Returns an error when serialisation fails.
pub fn write_json(
    df: &DataFrame,
    pretty: bool,
) -> Result<String, String> {
    let n_rows = df.n_rows();
    if n_rows == 0 {
        return Ok("[]".to_string());
    }
    let names = df.column_names();
    let mut encoded_cols: Vec<(String, &Column)> = Vec::with_capacity(names.len());
    for name in &names {
        let column = df
            .column(name)
            .ok_or_else(|| format!("unknown column `{name}`"))?;
        let key = serde_json::to_string(name).map_err(|error| error.to_string())?;
        encoded_cols.push((key, column));
    }
    let mut out = String::new();
    out.push('[');
    if pretty {
        out.push('\n');
    }
    for row in 0..n_rows {
        if row > 0 {
            out.push(',');
            if pretty {
                out.push('\n');
            }
        }
        if pretty {
            out.push_str("  {");
            if !encoded_cols.is_empty() {
                out.push('\n');
            }
        } else {
            out.push('{');
        }
        for (col_idx, (key, column)) in encoded_cols.iter().enumerate() {
            if col_idx > 0 {
                out.push(',');
                if pretty {
                    out.push('\n');
                }
            }
            if pretty {
                out.push_str("    ");
                out.push_str(key);
                out.push_str(": ");
            } else {
                out.push_str(key);
                out.push(':');
            }
            write_json_cell(&mut out, column, row)?;
        }
        if pretty {
            if !encoded_cols.is_empty() {
                out.push('\n');
                out.push_str("  ");
            }
        }
        out.push('}');
    }
    if pretty {
        out.push('\n');
    }
    out.push(']');
    Ok(out)
}

/// Reads a frame from a JSON file.
///
/// # Errors
///
/// Returns an error when the file cannot be read or its contents are invalid.
pub fn read_json_file(path: &Path) -> Result<DataFrame, String> {
    let text = std::fs::read_to_string(path).map_err(|error| error.to_string())?;
    read_json_str(&text)
}

/// Writes a frame to a JSON file.
///
/// # Errors
///
/// Returns an error when the file cannot be written.
pub fn write_json_file(
    path: &Path,
    df: &DataFrame,
) -> Result<(), String> {
    let text = write_json(df, true)?;
    std::fs::write(path, text).map_err(|error| error.to_string())
}

/// Builds a frame from an array of row objects.
fn from_rows(rows: &[Value]) -> Result<DataFrame, String> {
    let mut names: Vec<String> = Vec::new();
    let mut seen: HashSet<&str> = HashSet::new();
    for row in rows {
        let object = row
            .as_object()
            .ok_or_else(|| "json array elements must be objects".to_string())?;
        for key in object.keys() {
            if seen.insert(key.as_str()) {
                names.push(key.clone());
            }
        }
    }
    let mut columns = Vec::with_capacity(names.len());
    for name in &names {
        let values: Vec<Option<&Value>> = rows
            .iter()
            .map(|row| {
                row.as_object()
                    .and_then(|object| object.get(name))
                    .filter(|value| !value.is_null())
            })
            .collect();
        columns.push((name.clone(), json_column(&values)?));
    }
    DataFrame::new(columns)
}

/// Builds a frame from an object of equal-length arrays.
fn from_columns(object: &Map<String, Value>) -> Result<DataFrame, String> {
    let mut names = Vec::with_capacity(object.len());
    let mut arrays: Vec<&Vec<Value>> = Vec::with_capacity(object.len());
    let mut length: Option<usize> = None;
    for (name, value) in object {
        let array = value
            .as_array()
            .ok_or_else(|| format!("json field `{name}` must be an array"))?;
        match length {
            | None => length = Some(array.len()),
            | Some(expected) if expected != array.len() => {
                return Err("json arrays have different lengths".to_string());
            },
            | Some(_) => {},
        }
        names.push(name.clone());
        arrays.push(array);
    }
    let mut columns = Vec::with_capacity(names.len());
    for (name, array) in names.into_iter().zip(arrays) {
        let values: Vec<Option<&Value>> = array
            .iter()
            .map(|value| {
                if value.is_null() {
                    None
                } else {
                    Some(value)
                }
            })
            .collect();
        columns.push((name, json_column(&values)?));
    }
    DataFrame::new(columns)
}

/// Infers a column type from JSON scalars and materialises it.
fn json_column(values: &[Option<&Value>]) -> Result<Column, String> {
    let has_null = values.iter().any(Option::is_none);
    let present: Vec<&Value> = values.iter().filter_map(|value| *value).collect();
    if !present.is_empty() && present.iter().all(|value| value.is_array()) {
        return json_list_column(values);
    }
    let dtype = if present.is_empty() {
        DType::Str
    } else if present
        .iter()
        .all(|value| value.is_number() && value.as_i64().is_some())
    {
        DType::Int
    } else if present.iter().all(|value| value.is_number()) {
        DType::Float
    } else if present.iter().all(|value| value.is_boolean()) {
        DType::Bool
    } else if present
        .iter()
        .all(|value| value.is_number() || value.is_boolean() || value.is_string())
    {
        DType::Str
    } else {
        return Err("json values must be numbers, booleans or strings".to_string());
    };
    let mut mask = Vec::with_capacity(values.len());
    let mut ints = Vec::with_capacity(values.len());
    let mut floats = Vec::with_capacity(values.len());
    let mut booleans = Vec::with_capacity(values.len());
    let mut strings = Vec::with_capacity(values.len());
    let mut bytes_vec = Vec::with_capacity(values.len());
    for value in values {
        mask.push(value.is_none());
        match dtype {
            | DType::Int => {
                ints.push(match value {
                    | Some(value) => {
                        value
                            .as_i64()
                            .ok_or_else(|| "expected an integer json value".to_string())?
                    },
                    | None => 0,
                });
            },
            | DType::Float => {
                floats.push(match value {
                    | Some(value) => {
                        value
                            .as_f64()
                            .ok_or_else(|| "expected a numeric json value".to_string())?
                    },
                    | None => 0.0,
                });
            },
            | DType::Bool => {
                booleans.push(match value {
                    | Some(value) => {
                        value
                            .as_bool()
                            .ok_or_else(|| "expected a boolean json value".to_string())?
                    },
                    | None => false,
                });
            },
            | DType::Str => {
                strings.push(match value {
                    | Some(value) => {
                        stringify(value).ok_or_else(|| {
                            "unsupported json value in a string column".to_string()
                        })?
                    },
                    | None => String::new(),
                });
            },
            | DType::Bytes => {
                bytes_vec.push(match value {
                    | Some(value) => {
                        stringify(value)
                            .ok_or_else(|| "unsupported json value in a bytes column".to_string())?
                            .into_bytes()
                    },
                    | None => Vec::new(),
                });
            },
            | DType::List => return Err("a json list column needs array values".to_string()),
        }
    }
    let base = match dtype {
        | DType::Int => Column::Int(ints),
        | DType::Float => Column::Float(floats),
        | DType::Bool => Column::Bool(booleans),
        | DType::Str => Column::Str(strings),
        | DType::Bytes => Column::Bytes(bytes_vec),
        | DType::List => return Err("a json list column needs array values".to_string()),
    };
    if has_null {
        Ok(Column::Nullable(Box::new(base), mask))
    } else {
        Ok(base)
    }
}

/// Builds a list column from JSON arrays, recursing for nested arrays.
fn json_list_column(values: &[Option<&Value>]) -> Result<Column, String> {
    let mut offsets = Vec::with_capacity(values.len() + 1);
    offsets.push(0usize);
    let mut flat: Vec<Option<&Value>> = Vec::new();
    let mut mask = Vec::with_capacity(values.len());
    for value in values {
        match value {
            | Some(Value::Array(items)) => {
                mask.push(false);
                flat.extend(items.iter().map(|item| {
                    if item.is_null() {
                        None
                    } else {
                        Some(item)
                    }
                }));
            },
            | Some(_) => return Err("json list column values must be arrays".to_string()),
            | None => mask.push(true),
        }
        offsets.push(flat.len());
    }
    let inner = json_column(&flat)?;
    let list = Column::List {
        values: Box::new(inner),
        offsets,
    };
    if mask.iter().any(|&missing| missing) {
        Ok(Column::Nullable(Box::new(list), mask))
    } else {
        Ok(list)
    }
}

/// Renders a JSON scalar as text.
fn stringify(value: &Value) -> Option<String> {
    match value {
        | Value::String(text) => Some(text.clone()),
        | Value::Number(number) => Some(number.to_string()),
        | Value::Bool(flag) => Some(flag.to_string()),
        | _ => None,
    }
}

/// Appends one frame cell as a JSON value to `out`.
fn write_json_cell(
    out: &mut String,
    column: &Column,
    row: usize,
) -> Result<(), String> {
    if column.is_null(row) {
        out.push_str("null");
        return Ok(());
    }
    match column {
        | Column::Int(values) => {
            let value = values
                .get(row)
                .ok_or_else(|| format!("row {row} is out of bounds"))?;
            out.push_str(&value.to_string());
            Ok(())
        },
        | Column::Float(values) => {
            let value = values
                .get(row)
                .copied()
                .ok_or_else(|| format!("row {row} is out of bounds"))?;
            match Number::from_f64(value) {
                | Some(number) => out.push_str(&number.to_string()),
                | None => out.push_str("null"),
            }
            Ok(())
        },
        | Column::Bool(values) => {
            let value = values
                .get(row)
                .copied()
                .ok_or_else(|| format!("row {row} is out of bounds"))?;
            out.push_str(if value { "true" } else { "false" });
            Ok(())
        },
        | Column::Str(values) => {
            let value = values
                .get(row)
                .ok_or_else(|| format!("row {row} is out of bounds"))?;
            let encoded = serde_json::to_string(value).map_err(|error| error.to_string())?;
            out.push_str(&encoded);
            Ok(())
        },
        | Column::Bytes(values) => {
            let value = values
                .get(row)
                .ok_or_else(|| format!("row {row} is out of bounds"))?;
            let encoded =
                serde_json::to_string(&encode_hex(value)).map_err(|error| error.to_string())?;
            out.push_str(&encoded);
            Ok(())
        },
        | Column::Nullable(inner, _) => write_json_cell(out, inner, row),
        | Column::List { values, offsets } => {
            let start = offsets.get(row).copied().unwrap_or(0);
            let end = offsets.get(row + 1).copied().unwrap_or(start);
            out.push('[');
            for index in start..end {
                if index > start {
                    out.push(',');
                }
                write_json_cell(out, values, index)?;
            }
            out.push(']');
            Ok(())
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_all_types() {
        let frame = DataFrame::new(vec![
            ("id".to_string(), Column::Int(vec![1, 2, 3])),
            ("score".to_string(), Column::Float(vec![1.5, 2.5, 3.5])),
            ("ok".to_string(), Column::Bool(vec![true, false, true])),
            (
                "name".to_string(),
                Column::Str(vec!["a".to_string(), "b".to_string(), "c".to_string()]),
            ),
        ])
        .expect("frame");
        let text = write_json(&frame, true).expect("write");
        let parsed = read_json_str(&text).expect("read");
        assert_eq!(parsed.n_rows(), 3);
        for name in frame.column_names() {
            assert_eq!(parsed.column(name), frame.column(name), "column {name}");
        }
    }

    #[test]
    fn array_of_objects_with_missing_keys() {
        let text = r#"[{"a":1,"b":"x"},{"a":2},{"b":"y"}]"#;
        let frame = read_json_str(text).expect("read");
        assert_eq!(frame.n_rows(), 3);
        let a = frame.column("a").expect("a");
        assert!(a.is_null(2));
        assert_eq!(a.get_f64(0), Some(1.0));
        let b = frame.column("b").expect("b");
        assert!(b.is_null(1));
    }

    #[test]
    fn object_of_arrays() {
        let text = r#"{"a":[1,2],"b":[true,false]}"#;
        let frame = read_json_str(text).expect("read");
        assert_eq!(frame.n_rows(), 2);
        assert_eq!(frame.column("a"), Some(&Column::Int(vec![1, 2])));
        assert_eq!(frame.column("b"), Some(&Column::Bool(vec![true, false])));
    }

    #[test]
    fn mixed_scalars_become_strings() {
        let text = r#"[{"v":1},{"v":"two"}]"#;
        let frame = read_json_str(text).expect("read");
        assert_eq!(
            frame.column("v"),
            Some(&Column::Str(vec!["1".to_string(), "two".to_string()]))
        );
    }

    #[test]
    fn nullable_column_round_trips() {
        let frame = DataFrame::new(vec![(
            "v".to_string(),
            Column::Nullable(
                Box::new(Column::Float(vec![1.0, 0.0, 3.0])),
                vec![false, true, false],
            ),
        )])
        .expect("frame");
        let text = write_json(&frame, false).expect("write");
        let parsed = read_json_str(&text).expect("read");
        let column = parsed.column("v").expect("v");
        assert!(column.is_null(1));
        assert_eq!(column.get_f64(0), Some(1.0));
        assert_eq!(column.get_f64(2), Some(3.0));
    }

    #[test]
    fn invalid_documents_are_rejected() {
        assert!(read_json_str("not json").is_err());
        assert!(read_json_str("42").is_err());
        assert!(read_json_str(r#"[{"a":{"b":1}}]"#).is_err());
        assert!(read_json_str(r#"{"a":[1],"b":[1,2]}"#).is_err());
    }

    #[test]
    fn file_round_trip() {
        let path = std::env::temp_dir().join(format!(
            "rssn_frame_json_{}_{}.json",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map_or(0, |duration| duration.as_nanos())
        ));
        let frame =
            DataFrame::new(vec![("n".to_string(), Column::Int(vec![7, 8]))]).expect("frame");
        write_json_file(&path, &frame).expect("write");
        let read = read_json_file(&path).expect("read");
        let _ = std::fs::remove_file(&path);
        assert_eq!(read.column("n"), frame.column("n"));
    }

    #[test]
    fn write_json_preserves_column_order() {
        let frame = DataFrame::new(vec![
            ("zeta".to_string(), Column::Int(vec![1, 2])),
            (
                "alpha".to_string(),
                Column::Str(vec!["a\"b".to_string(), "c".to_string()]),
            ),
            (
                "mu".to_string(),
                Column::Bytes(vec![vec![0x0a, 0xff], Vec::new()]),
            ),
        ])
        .expect("frame");
        let compact = write_json(&frame, false).expect("compact");
        assert_eq!(
            compact,
            r#"[{"zeta":1,"alpha":"a\"b","mu":"0aff"},{"zeta":2,"alpha":"c","mu":""}]"#
        );
        let pretty = write_json(&frame, true).expect("pretty");
        let pos_zeta = pretty.find("\"zeta\"").expect("zeta");
        let pos_alpha = pretty.find("\"alpha\"").expect("alpha");
        let pos_mu = pretty.find("\"mu\"").expect("mu");
        assert!(pos_zeta < pos_alpha && pos_alpha < pos_mu);
        let reparsed = read_json_str(&pretty).expect("reparse");
        assert_eq!(reparsed.n_rows(), 2);
    }
}
