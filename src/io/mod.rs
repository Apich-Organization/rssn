//! # Array file I/O
//!
//! Reading and writing dense 2-D arrays as `.npy`, CSV and JSON. Simulation
//! models use this to dump fields; nothing here knows about the expression
//! graph. The sibling modules add standalone SVG plotting ([`plot`]), 8-bit
//! PGM/PPM raster export ([`raster`]) and term rendering for LaTeX
//! ([`latex`]), Typst ([`typst`]) and a Unicode pretty printer ([`pretty`]).
//!
//! # Limitations
//!
//! Prioritised gaps of the layer as a whole; each submodule documents its own
//! details.
//!
//! - **Data formats.** CSV, JSON and `.npy` store dense `f64` 2-D arrays only:
//!   no integer/complex/string/structured dtypes, no headers or metadata, no
//!   compression, and no streaming (whole files are built in memory). CSV is
//!   unquoted (a value must not contain a comma) and strict on input; JSON
//!   uses `serde_json`'s default float formatting, which is not guaranteed
//!   bit-exact for every `f64` on re-read (`float_roundtrip` is not enabled).
//! - **`.npy` feature gating.** `.npy` support needs the `npy` feature, which
//!   is off by default. Without it [`write_npy_file`] and [`read_npy_file`]
//!   still compile but always return an error naming the missing feature.
//! - **Plots are SVG** (text, annotated, standalone files): 2-D coverage is
//!   functions, series, expression terms, vector fields, scaled quivers,
//!   contour lines, filled contour bands and streamlines; 3-D coverage is
//!   surfaces, parametric curves, polylines and vector fields through one
//!   fixed oblique projection. Missing: log axes, subplots, histograms/box
//!   plots, error bars, chart-like band/area plots, contour labels and any
//!   interactive output.
//! - **Raster export is 8-bit netpbm only** ([`raster`]: binary PGM P5,
//!   8/16-bit, and PPM P6): [`plot::plot_scalar_field_pgm`] and
//!   [`plot::plot_scalar_field_ppm`] rasterise a sampled field, but there is no
//!   PNG/JPEG/other compressed format (that needs an encoder dependency), no
//!   alpha/16-bit PPM, and no resampling, text or annotation in the raster
//!   path.
//! - **LaTeX and Typst** emit math fragments, not documents: no delimiters,
//!   preamble or styling commands, a fixed function/Greek/operator table, and
//!   no `cases`/`aligned` equivalent. LaTeX escapes literal strings for
//!   `\text{...}` and Typst escapes string literals, but symbol names are
//!   assumed identifier-like. Lowercase `omicron` and variant letters
//!   (`varepsilon`, `varphi`, ...) are not distinguished.
//! - **Markup walker** (shared by LaTeX/Typst) special-cases core arithmetic,
//!   lists/matrices, derivatives, binders, limits and relations; other
//!   operators fall back to a plain call. Assumptions, domains, piecewise
//!   definitions, sets/intervals and vectors get no dedicated notation.
//! - **Pretty printer** is a monospace Unicode grid: wide/combining glyphs
//!   misalign, Unicode super/subscripts cover digits (plus `-` and `n` above)
//!   and no wrapping, colours or ANSI styles exist.
//! - **No new dependencies** are allowed for this layer, so anything needing
//!   an encoder, a font or a typesetter stays a documented gap instead.

pub mod latex;
pub(crate) mod markup;
pub mod plot;
pub mod pretty;
pub mod raster;
pub mod typst;

use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Write;
use std::path::Path;

use ndarray::Array2;
#[cfg(feature = "npy")]
use ndarray_npy::read_npy;
#[cfg(feature = "npy")]
use ndarray_npy::write_npy;
use serde_json;

/// Writes a 2D `ndarray::Array` to a `.npy` file.
///
/// # Arguments
/// * `filename` - The path to the `.npy` file.
/// * `arr` - The array to write.
///
/// # Errors
///
/// This function will return an error if the file cannot be written to.
///
/// # Panics
/// Panics if the write fails.
#[cfg(feature = "npy")]
pub fn write_npy_file<P: AsRef<Path>>(
    filename: P,
    arr: &Array2<f64>,
) -> Result<(), String> {
    write_npy(filename, arr).map_err(|e| e.to_string())
}

/// Writes a 2D array to a `.npy` file.
///
/// # Errors
/// Always fails: the crate was built without the `npy` feature.
#[cfg(not(feature = "npy"))]
pub fn write_npy_file<P: AsRef<Path>>(
    _filename: P,
    _arr: &Array2<f64>,
) -> Result<(), String> {
    Err("Feature 'npy' is required \
         for .npy support"
        .to_string())
}

/// Reads a 2D `ndarray::Array` from a `.npy` file.
///
/// # Arguments
/// * `filename` - The path to the `.npy` file.
///
/// # Returns
/// The read array as an `ndarray::Array2<f64>`.
///
/// # Errors
///
/// This function will return an error if the file cannot be read.
///
/// # Panics
/// Panics if the read fails.
#[cfg(feature = "npy")]
pub fn read_npy_file<P: AsRef<Path>>(filename: P) -> Result<Array2<f64>, String> {
    read_npy(filename).map_err(|e| e.to_string())
}

/// Reads a 2D array from a `.npy` file.
///
/// # Errors
/// Always fails: the crate was built without the `npy` feature.
#[cfg(not(feature = "npy"))]
pub fn read_npy_file<P: AsRef<Path>>(_filename: P) -> Result<Array2<f64>, String> {
    Err("Feature 'npy' is required \
         for .npy support"
        .to_string())
}

/// Writes a 2D `ndarray::Array` to a CSV file.
///
/// # Errors
///
/// This function will return an error if the file cannot be written to.
pub fn write_csv_file<P: AsRef<Path>>(
    filename: P,
    arr: &Array2<f64>,
) -> Result<(), String> {
    let mut file = File::create(filename).map_err(|e| e.to_string())?;

    for row in arr.outer_iter() {
        let line = row
            .iter()
            .map(|&v| v.to_string())
            .collect::<Vec<_>>()
            .join(",");

        writeln!(file, "{line}").map_err(|e| e.to_string())?;
    }

    Ok(())
}

/// Reads a 2D `ndarray::Array` from a CSV file.
///
/// # Errors
///
/// This function will return an error if the file cannot be read or parsed.
pub fn read_csv_file<P: AsRef<Path>>(filename: P) -> Result<Array2<f64>, String> {
    let file = File::open(filename).map_err(|e| e.to_string())?;

    let reader = BufReader::new(file);

    let mut data = Vec::new();

    let mut rows = 0;

    let mut cols = 0;

    for line_res in reader.lines() {
        let line = line_res.map_err(|e| e.to_string())?;

        let line = line.trim();

        if line.is_empty() {
            continue;
        }

        let row_data: Vec<f64> = line
            .split(',')
            .map(|s| s.trim().parse::<f64>().map_err(|e| e.to_string()))
            .collect::<Result<Vec<_>, String>>()?;

        if rows == 0 {
            cols = row_data.len();
        } else if row_data.len() != cols {
            return Err("Inconsistent column \
                 count in CSV"
                .to_string());
        }

        data.extend(row_data);

        rows += 1;
    }

    if rows == 0 {
        return Ok(Array2::zeros((0, 0)));
    }

    Array2::from_shape_vec((rows, cols), data).map_err(|e| e.to_string())
}

/// Writes a 2D `ndarray::Array` to a JSON file.
///
/// # Errors
///
/// This function will return an error if the file cannot be written to.
pub fn write_json_file<P: AsRef<Path>>(
    filename: P,
    arr: &Array2<f64>,
) -> Result<(), String> {
    let file = File::create(filename).map_err(|e| e.to_string())?;

    serde_json::to_writer_pretty(file, arr).map_err(|e| e.to_string())
}

/// Reads a 2D `ndarray::Array` from a JSON file.
///
/// # Errors
///
/// This function will return an error if the file cannot be read or parsed.
pub fn read_json_file<P: AsRef<Path>>(filename: P) -> Result<Array2<f64>, String> {
    let file = File::open(filename).map_err(|e| e.to_string())?;

    let reader = BufReader::new(file);

    serde_json::from_reader(reader).map_err(|e| e.to_string())
}
