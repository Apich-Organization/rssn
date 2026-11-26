//! RFC 4180 CSV parsing and writing for typed frames.

use std::path::Path;

use super::column::Column;
use super::column::DType;
use super::column::encode_hex;
use super::dataframe::DataFrame;

/// Options controlling CSV parsing and writing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CsvOptions {
    /// Field separator, normally `b','`.
    pub delimiter: u8,
    /// Whether the first record holds the column names.
    pub has_header: bool,
    /// Whether to infer `Int`/`Float`/`Bool`/`Str` per column.
    ///
    /// With inference disabled every column is read as `Str`.
    pub infer_types: bool,
    /// Unquoted fields equal to one of these tokens are read as missing.
    pub null_values: Vec<String>,
    /// Whether an integer column containing missing values stays `Nullable(Int)`.
    ///
    /// When `false`, integer columns with missing values are widened to
    /// `Nullable(Float)`.
    pub preserve_nullable_int: bool,
}

impl Default for CsvOptions {
    fn default() -> Self {
        Self {
            delimiter: b',',
            has_header: true,
            infer_types: true,
            null_values: vec![
                String::new(),
                "null".to_string(),
                "NULL".to_string(),
                "NA".to_string(),
            ],
            preserve_nullable_int: true,
        }
    }
}

/// Reads a frame from CSV text.
///
/// Quoting follows RFC 4180: fields may be quoted with `"`, quotes are
/// escaped by doubling, and quoted fields may contain the delimiter and line
/// breaks. A quoted empty field is an empty string; an unquoted token listed
/// in [`CsvOptions::null_values`] is missing.
///
/// # Errors
///
/// Returns an error for unterminated quotes, rows whose field count differs
/// from the header, or duplicate column names.
pub fn read_csv_str(
    text: &str,
    opts: &CsvOptions,
) -> Result<DataFrame, String> {
    let records = parse_records(text, opts.delimiter)?;
    if records.is_empty() {
        return DataFrame::new(Vec::new());
    }
    let (names, rows): (Vec<String>, &[Vec<(String, bool)>]) = if opts.has_header {
        (
            records.first().map_or_else(Vec::new, |record| {
                record.iter().map(|(value, _)| value.clone()).collect()
            }),
            records.get(1..).unwrap_or_default(),
        )
    } else {
        let width = records.first().map_or(0, Vec::len);
        (
            (0..width).map(|index| format!("col{index}")).collect(),
            &records[..],
        )
    };
    let width = names.len();
    for (index, record) in rows.iter().enumerate() {
        if record.len() != width {
            return Err(format!(
                "row {} has {} fields but the header has {width}",
                index + 1,
                record.len()
            ));
        }
    }
    let mut columns = Vec::with_capacity(width);
    for column in 0..width {
        let tokens: Vec<Option<String>> = rows
            .iter()
            .map(|record| {
                record.get(column).and_then(|(value, quoted)| {
                    if !quoted && opts.null_values.iter().any(|null| null == value) {
                        None
                    } else {
                        Some(value.clone())
                    }
                })
            })
            .collect();
        let name = names
            .get(column)
            .cloned()
            .ok_or_else(|| "missing column name".to_string())?;
        columns.push((
            name,
            build_column(&tokens, opts.infer_types, opts.preserve_nullable_int)?,
        ));
    }
    DataFrame::new(columns)
}

/// Writes a frame as CSV text.
///
/// Fields are quoted when they contain the delimiter, a quote or a line
/// break; empty strings are quoted so they stay distinct from missing values.
#[must_use]
pub fn write_csv(
    df: &DataFrame,
    opts: &CsvOptions,
) -> String {
    let mut output = String::new();
    if opts.has_header {
        for (index, name) in df.column_names().iter().enumerate() {
            if index > 0 {
                output.push(char::from(opts.delimiter));
            }
            output.push_str(&quote(name, opts.delimiter));
        }
        output.push('\n');
    }
    for row in 0..df.n_rows() {
        for column in 0..df.n_cols() {
            if column > 0 {
                output.push(char::from(opts.delimiter));
            }
            let value = df.column_at(column).and_then(|values| cell(values, row));
            if let Some(text) = value {
                output.push_str(&quote(&text, opts.delimiter));
            }
        }
        output.push('\n');
    }
    output
}

/// Reads a frame from a CSV file.
///
/// # Errors
///
/// Returns an error when the file cannot be read or the contents are invalid.
pub fn read_csv_file(
    path: &Path,
    opts: &CsvOptions,
) -> Result<DataFrame, String> {
    let text = std::fs::read_to_string(path).map_err(|error| error.to_string())?;
    read_csv_str(&text, opts)
}

/// Writes a frame to a CSV file.
///
/// # Errors
///
/// Returns an error when the file cannot be written.
pub fn write_csv_file(
    path: &Path,
    df: &DataFrame,
    opts: &CsvOptions,
) -> Result<(), String> {
    std::fs::write(path, write_csv(df, opts)).map_err(|error| error.to_string())
}

/// One parsed field: its text and whether it was quoted.
type Field = (String, bool);

/// Splits CSV text into records of fields.
fn parse_records(
    text: &str,
    delimiter: u8,
) -> Result<Vec<Vec<Field>>, String> {
    let bytes = text.as_bytes();
    let mut records: Vec<Vec<Field>> = Vec::new();
    let mut record: Vec<Field> = Vec::new();
    let mut field: Vec<u8> = Vec::new();
    let mut quoted = false;
    let mut in_quotes = false;
    let mut field_started = false;
    let mut index = 0usize;
    while index < bytes.len() {
        let byte = bytes.get(index).copied().unwrap_or(0);
        if in_quotes {
            if byte == b'"' {
                if bytes.get(index + 1).copied() == Some(b'"') {
                    field.push(b'"');
                    index += 2;
                    continue;
                }
                in_quotes = false;
                index += 1;
                continue;
            }
            field.push(byte);
            index += 1;
            continue;
        }
        match byte {
            | b'"' if field.is_empty() => {
                in_quotes = true;
                quoted = true;
                field_started = true;
            },
            | b'"' => {
                field.push(b'"');
                field_started = true;
            },
            | _ if byte == delimiter => {
                record.push((finish_field(&field)?, quoted));
                field.clear();
                quoted = false;
                field_started = false;
            },
            | b'\n' => {
                record.push((finish_field(&field)?, quoted));
                field.clear();
                quoted = false;
                field_started = false;
                records.push(std::mem::take(&mut record));
            },
            | b'\r' => {
                record.push((finish_field(&field)?, quoted));
                field.clear();
                quoted = false;
                field_started = false;
                records.push(std::mem::take(&mut record));
                if bytes.get(index + 1).copied() == Some(b'\n') {
                    index += 1;
                }
            },
            | _ => {
                field.push(byte);
                field_started = true;
            },
        }
        index += 1;
    }
    if in_quotes {
        return Err("unterminated quoted csv field".to_string());
    }
    if field_started || !field.is_empty() || !record.is_empty() {
        record.push((finish_field(&field)?, quoted));
        records.push(record);
    }
    records.retain(|row| {
        !(row.len() == 1
            && row
                .first()
                .is_some_and(|(text, quoted)| text.is_empty() && !*quoted))
    });
    Ok(records)
}

/// Converts the raw bytes of a field into text.
fn finish_field(field: &[u8]) -> Result<String, String> {
    String::from_utf8(field.to_vec()).map_err(|_| "csv field is not valid utf-8".to_string())
}

/// Infers a column type and materialises the column.
fn build_column(
    tokens: &[Option<String>],
    infer_types: bool,
    preserve_nullable_int: bool,
) -> Result<Column, String> {
    let has_null = tokens.iter().any(Option::is_none);
    let present: Vec<&str> = tokens.iter().filter_map(|token| token.as_deref()).collect();
    let dtype = if !infer_types || present.is_empty() {
        DType::Str
    } else if present.iter().all(|value| value.parse::<i64>().is_ok()) {
        if has_null && !preserve_nullable_int {
            DType::Float
        } else {
            DType::Int
        }
    } else if present.iter().all(|value| value.parse::<f64>().is_ok()) {
        DType::Float
    } else if present.iter().all(|value| parse_bool(value).is_some()) {
        DType::Bool
    } else {
        DType::Str
    };
    let mut mask = Vec::with_capacity(tokens.len());
    let mut ints = Vec::with_capacity(tokens.len());
    let mut floats = Vec::with_capacity(tokens.len());
    let mut booleans = Vec::with_capacity(tokens.len());
    let mut strings = Vec::with_capacity(tokens.len());
    let mut bytes_vec = Vec::with_capacity(tokens.len());
    for token in tokens {
        mask.push(token.is_none());
        match dtype {
            | DType::Int => {
                ints.push(match token {
                    | Some(value) => {
                        value
                            .parse::<i64>()
                            .map_err(|_| format!("`{value}` is not an integer"))?
                    },
                    | None => 0,
                })
            },
            | DType::Float => {
                floats.push(match token {
                    | Some(value) => {
                        value
                            .parse::<f64>()
                            .map_err(|_| format!("`{value}` is not a number"))?
                    },
                    | None => 0.0,
                })
            },
            | DType::Bool => {
                booleans.push(match token {
                    | Some(value) => {
                        parse_bool(value).ok_or_else(|| format!("`{value}` is not a boolean"))?
                    },
                    | None => false,
                })
            },
            | DType::Str => strings.push(token.clone().unwrap_or_default()),
            | DType::Bytes => {
                bytes_vec.push(
                    token
                        .as_deref()
                        .map_or_else(Vec::new, |s| s.as_bytes().to_vec()),
                )
            },
            | DType::List => {
                return Err("csv cannot infer a list column; use JSON or Parquet".to_string());
            },
        }
    }
    let base = match dtype {
        | DType::Int => Column::Int(ints),
        | DType::Float => Column::Float(floats),
        | DType::Bool => Column::Bool(booleans),
        | DType::Str => Column::Str(strings),
        | DType::Bytes => Column::Bytes(bytes_vec),
        | DType::List => {
            return Err("csv cannot infer a list column; use JSON or Parquet".to_string());
        },
    };
    if has_null {
        Ok(Column::Nullable(Box::new(base), mask))
    } else {
        Ok(base)
    }
}

/// Parses the two spellings accepted for booleans.
fn parse_bool(text: &str) -> Option<bool> {
    if text.eq_ignore_ascii_case("true") {
        Some(true)
    } else if text.eq_ignore_ascii_case("false") {
        Some(false)
    } else {
        None
    }
}

/// Formats one non-missing cell.
fn cell(
    column: &Column,
    row: usize,
) -> Option<String> {
    if column.is_null(row) {
        return None;
    }
    match column {
        | Column::Int(values) => values.get(row).map(ToString::to_string),
        | Column::Float(values) => values.get(row).map(ToString::to_string),
        | Column::Bool(values) => {
            values
                .get(row)
                .map(|value| if *value { "true" } else { "false" }.to_string())
        },
        | Column::Str(values) => values.get(row).cloned(),
        | Column::Bytes(values) => values.get(row).map(|bytes| encode_hex(bytes)),
        | Column::List { values, offsets } => {
            let start = offsets.get(row).copied().unwrap_or(0);
            let end = offsets.get(row + 1).copied().unwrap_or(start);
            let mut text = String::from("[");
            for index in start..end {
                if index > start {
                    text.push(',');
                }
                text.push_str(&json_scalar(values, index));
            }
            text.push(']');
            Some(text)
        },
        | Column::Nullable(inner, _) => cell(inner, row),
    }
}

/// Renders one element as a JSON literal, used for list cells in a CSV.
fn json_scalar(
    column: &Column,
    row: usize,
) -> String {
    match column {
        | Column::Int(values) => {
            values
                .get(row)
                .map_or_else(|| "null".to_string(), |v| v.to_string())
        },
        | Column::Float(values) => {
            values.get(row).map_or_else(
                || "null".to_string(),
                |value| {
                    serde_json::Number::from_f64(*value)
                        .map_or_else(|| "null".to_string(), |n| n.to_string())
                },
            )
        },
        | Column::Bool(values) => {
            values.get(row).map_or_else(
                || "null".to_string(),
                |v| if *v { "true" } else { "false" }.to_string(),
            )
        },
        | Column::Str(values) => {
            values.get(row).map_or_else(
                || "null".to_string(),
                |value| serde_json::to_string(value).unwrap_or_else(|_| "\"\"".to_string()),
            )
        },
        | Column::Bytes(values) => {
            values.get(row).map_or_else(
                || "null".to_string(),
                |value| {
                    serde_json::to_string(&encode_hex(value)).unwrap_or_else(|_| "\"\"".to_string())
                },
            )
        },
        | Column::List { values, offsets } => {
            let start = offsets.get(row).copied().unwrap_or(0);
            let end = offsets.get(row + 1).copied().unwrap_or(start);
            let mut text = String::from("[");
            for index in start..end {
                if index > start {
                    text.push(',');
                }
                text.push_str(&json_scalar(values, index));
            }
            text.push(']');
            text
        },
        | Column::Nullable(inner, mask) => {
            if mask.get(row).copied().unwrap_or(true) {
                "null".to_string()
            } else {
                json_scalar(inner, row)
            }
        },
    }
}

/// Quotes a field when required and doubles embedded quotes.
fn quote(
    value: &str,
    delimiter: u8,
) -> String {
    let needs_quotes = value.is_empty()
        || value
            .bytes()
            .any(|byte| byte == delimiter || byte == b'"' || byte == b'\n' || byte == b'\r');
    if needs_quotes {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quoted_fields_round_trip() {
        let frame = DataFrame::new(vec![
            (
                "text".to_string(),
                Column::Str(vec![
                    "plain".to_string(),
                    "a,b".to_string(),
                    "say \"hi\"".to_string(),
                    "two\nlines".to_string(),
                ]),
            ),
            ("count".to_string(), Column::Int(vec![1, 2, 3, 4])),
        ])
        .expect("frame");
        let text = write_csv(&frame, &CsvOptions::default());
        let parsed = read_csv_str(&text, &CsvOptions::default()).expect("parse");
        assert_eq!(parsed.n_rows(), 4);
        assert_eq!(parsed.column("text"), frame.column("text"));
        assert_eq!(parsed.column("count"), frame.column("count"));
    }

    #[test]
    fn blank_in_integer_column_preserves_nullable_int_by_default() {
        let text = "n,label\n1,a\n,b\n3,c\n";
        let frame = read_csv_str(text, &CsvOptions::default()).expect("parse");
        let column = frame.column("n").expect("column");
        assert_eq!(column.dtype(), DType::Int);
        assert!(column.is_nullable());
        assert!(column.is_null(1));
        assert_eq!(column.get_f64(0), Some(1.0));
        assert_eq!(column.get_f64(2), Some(3.0));
        assert_eq!(column.len(), 3);

        let widen_opts = CsvOptions {
            preserve_nullable_int: false,
            ..CsvOptions::default()
        };
        let widened = read_csv_str(text, &widen_opts).expect("parse widened");
        let wcol = widened.column("n").expect("column");
        assert_eq!(wcol.dtype(), DType::Float);
        assert!(wcol.is_nullable());
        assert!(wcol.is_null(1));
        assert_eq!(wcol.get_f64(0), Some(1.0));
        assert_eq!(wcol.get_f64(2), Some(3.0));
    }

    #[test]
    fn bytes_column_writes_hex_in_csv() {
        let frame = DataFrame::new(vec![(
            "raw".to_string(),
            Column::Bytes(vec![vec![0x00, 0xab, 0xff], Vec::new()]),
        )])
        .expect("frame");
        let text = write_csv(&frame, &CsvOptions::default());
        assert_eq!(text, "raw\n00abff\n\"\"\n");
    }

    #[test]
    fn inference_covers_int_float_bool_str() {
        let text = "i,f,b,s\n1,1.5,true,x\n2,2.5,false,y\n";
        let frame = read_csv_str(text, &CsvOptions::default()).expect("parse");
        assert_eq!(frame.column("i"), Some(&Column::Int(vec![1, 2])));
        assert_eq!(frame.column("f"), Some(&Column::Float(vec![1.5, 2.5])));
        assert_eq!(frame.column("b"), Some(&Column::Bool(vec![true, false])));
        assert_eq!(
            frame.column("s"),
            Some(&Column::Str(vec!["x".to_string(), "y".to_string()]))
        );
    }

    #[test]
    fn hostile_row_with_extra_fields_is_rejected() {
        let text = "a,b\n1,2,3\n";
        assert!(read_csv_str(text, &CsvOptions::default()).is_err());
        let short = "a,b\n1\n";
        assert!(read_csv_str(short, &CsvOptions::default()).is_err());
    }

    #[test]
    fn unterminated_quote_is_rejected() {
        assert!(read_csv_str("a\n\"open\n", &CsvOptions::default()).is_err());
    }

    #[test]
    fn list_cells_are_written_as_json_arrays() {
        let list = Column::new_list(
            Column::Nullable(
                Box::new(Column::Int(vec![1, 0, 3])),
                vec![false, true, false],
            ),
            vec![0, 2, 3],
        )
        .expect("list");
        let frame = DataFrame::new(vec![("xs".to_string(), list)]).expect("frame");
        let text = write_csv(&frame, &CsvOptions::default());
        // The comma inside the JSON array is quoted like any other cell.
        assert_eq!(
            text,
            "xs
\"[1,null]\"
[3]
"
        );
    }

    #[test]
    fn custom_delimiter_and_no_header() {
        let options = CsvOptions {
            delimiter: b';',
            has_header: false,
            infer_types: true,
            null_values: vec![String::new()],
            preserve_nullable_int: true,
        };
        let frame = read_csv_str("1;2\n3;4\n", &options).expect("parse");
        assert_eq!(frame.column_names(), vec!["col0", "col1"]);
        assert_eq!(frame.column("col1"), Some(&Column::Int(vec![2, 4])));
        let text = write_csv(&frame, &options);
        assert!(text.starts_with("1;2\n"));
    }

    #[test]
    fn quoted_empty_string_is_not_null() {
        let frame = read_csv_str("a\n\"\"\n", &CsvOptions::default()).expect("parse");
        assert_eq!(frame.column("a"), Some(&Column::Str(vec![String::new()])));
    }

    #[test]
    fn file_round_trip() {
        let path = std::env::temp_dir().join(format!(
            "rssn_frame_csv_{}_{}.csv",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map_or(0, |duration| duration.as_nanos())
        ));
        let frame = DataFrame::new(vec![(
            "value".to_string(),
            Column::Float(vec![1.0, 2.5, -3.25]),
        )])
        .expect("frame");
        write_csv_file(&path, &frame, &CsvOptions::default()).expect("write");
        let read = read_csv_file(&path, &CsvOptions::default()).expect("read");
        let _ = std::fs::remove_file(&path);
        assert_eq!(read.column("value"), frame.column("value"));
    }
}
