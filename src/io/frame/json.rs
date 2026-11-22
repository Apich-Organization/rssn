//! JSON reading and writing for typed frames.

use std::collections::HashSet;
use std::path::Path;

use serde_json::Map;
use serde_json::Number;
use serde_json::Value;

use super::column::Column;
use super::column::DType;
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

/// Writes a frame as JSON text: an array with one object per row.
///
/// # Errors
///
/// Returns an error when serialisation fails.
pub fn write_json(
    df: &DataFrame,
    pretty: bool,
) -> Result<String, String> {
    let names = df.column_names();
    let mut rows = Vec::with_capacity(df.n_rows());
    for row in 0..df.n_rows() {
        let mut object = Map::new();
        for name in &names {
            let column = df
                .column(name)
                .ok_or_else(|| format!("unknown column `{name}`"))?;
            object.insert((*name).to_string(), json_value(column, row)?);
        }
        rows.push(Value::Object(object));
    }
    let value = Value::Array(rows);
    if pretty {
        serde_json::to_string_pretty(&value)
    } else {
        serde_json::to_string(&value)
    }
    .map_err(|error| error.to_string())
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
                })
            },
            | DType::Float => {
                floats.push(match value {
                    | Some(value) => {
                        value
                            .as_f64()
                            .ok_or_else(|| "expected a numeric json value".to_string())?
                    },
                    | None => 0.0,
                })
            },
            | DType::Bool => {
                booleans.push(match value {
                    | Some(value) => {
                        value
                            .as_bool()
                            .ok_or_else(|| "expected a boolean json value".to_string())?
                    },
                    | None => false,
                })
            },
            | DType::Str => {
                strings.push(match value {
                    | Some(value) => {
                        stringify(value).ok_or_else(|| {
                            "unsupported json value in a string column".to_string()
                        })?
                    },
                    | None => String::new(),
                })
            },
        }
    }
    let base = match dtype {
        | DType::Int => Column::Int(ints),
        | DType::Float => Column::Float(floats),
        | DType::Bool => Column::Bool(booleans),
        | DType::Str => Column::Str(strings),
    };
    if has_null {
        Ok(Column::Nullable(Box::new(base), mask))
    } else {
        Ok(base)
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

/// Converts one frame cell into a JSON value.
fn json_value(
    column: &Column,
    row: usize,
) -> Result<Value, String> {
    if column.is_null(row) {
        return Ok(Value::Null);
    }
    match column {
        | Column::Int(values) => {
            values
                .get(row)
                .map(|value| Value::from(*value))
                .ok_or_else(|| format!("row {row} is out of bounds"))
        },
        | Column::Float(values) => {
            let value = values
                .get(row)
                .copied()
                .ok_or_else(|| format!("row {row} is out of bounds"))?;
            Ok(Number::from_f64(value).map_or(Value::Null, Value::Number))
        },
        | Column::Bool(values) => {
            values
                .get(row)
                .map(|value| Value::Bool(*value))
                .ok_or_else(|| format!("row {row} is out of bounds"))
        },
        | Column::Str(values) => {
            values
                .get(row)
                .map(|value| Value::String(value.clone()))
                .ok_or_else(|| format!("row {row} is out of bounds"))
        },
        | Column::Nullable(inner, _) => json_value(inner, row),
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
        assert!(read_json_str(r#"[{"a":[1]}]"#).is_err());
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
}
