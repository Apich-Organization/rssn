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
//! - **Data formats.** The file helpers use dense rank-two `f64` arrays.
//!   CSV is unquoted and strict on input. JSON uses the default `serde_json`
//!   float parser, which can change the last bit on re-read. The `npy` codec
//!   reads arbitrary-rank `f32` and `f64` arrays and writes little-endian `f64`.
//!   It supports both byte orders and both C and Fortran input layouts,
//!   including rank zero, empty arrays and binary64 NaN payloads.
//! - **`.npy` feature gating.** The default build enables the `npy` feature.
//!   A minimal build without it retains [`write_npy_file`] and [`read_npy_file`].
//!   These functions return an error naming the missing feature.
//! - **Plots are SVG** (text, annotated, standalone files): 2-D coverage is
//!   functions, series, expression terms, vector fields, scaled quivers,
//!   contour lines, filled contour bands and streamlines; 3-D coverage is
//!   surfaces, parametric curves, polylines and vector fields through one
//!   fixed oblique projection. Missing: log axes, subplots, histograms/box
//!   plots, error bars, chart-like band/area plots, contour labels and any
//!   interactive output.
//! - **Raster export** ([`raster`]) supports binary PGM P5,
//!   8/16-bit, and PPM P6: [`plot::plot_scalar_field_pgm`] and
//!   [`plot::plot_scalar_field_ppm`] rasterise a sampled field, but there is no
//!   PNG/JPEG/other compressed format (that needs an encoder dependency), no
//!   alpha or 16-bit PPM, and no resampling, text or annotation in the raster
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
//! - **No PNG/JPEG** in this layer, so anything needing those encoders stays
//!   a documented gap; GIF animation I/O is in [`gif`] and SVG parsing/I/O is
//!   in [`svg`].

pub mod colormap;
/// Typed tabular data: columns, data frames, typed CSV and JSON.
#[allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    clippy::single_call_fn,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::redundant_else,
    clippy::if_not_else,
    clippy::suboptimal_flops,
    clippy::doc_markdown,
    clippy::needless_pass_by_value,
    clippy::too_many_lines,
    clippy::too_many_arguments,
    clippy::similar_names,
    clippy::manual_let_else,
    clippy::option_if_let_else,
    clippy::collapsible_if,
    clippy::unnecessary_wraps,
    clippy::many_single_char_names,
    clippy::manual_clamp,
    clippy::needless_range_loop,
    clippy::type_complexity,
    clippy::module_name_repetitions,
    clippy::must_use_candidate,
    clippy::missing_const_for_fn,
    clippy::use_self
)]
pub mod frame;
pub mod gif;
/// HDF5 files: groups, attributes and contiguous or chunked datasets.
#[allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    clippy::single_call_fn,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::redundant_else,
    clippy::if_not_else,
    clippy::suboptimal_flops,
    clippy::doc_markdown,
    clippy::needless_pass_by_value,
    clippy::too_many_lines,
    clippy::too_many_arguments,
    clippy::similar_names,
    clippy::manual_let_else,
    clippy::option_if_let_else,
    clippy::collapsible_if,
    clippy::unnecessary_wraps,
    clippy::many_single_char_names,
    clippy::manual_clamp,
    clippy::needless_range_loop,
    clippy::type_complexity,
    clippy::module_name_repetitions,
    clippy::must_use_candidate,
    clippy::missing_const_for_fn,
    clippy::use_self
)]
pub mod hdf5;
pub mod latex;
pub(crate) mod markup;
/// NetCDF classic and 64-bit offset files.
#[allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    clippy::single_call_fn,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::redundant_else,
    clippy::if_not_else,
    clippy::suboptimal_flops,
    clippy::doc_markdown,
    clippy::needless_pass_by_value,
    clippy::too_many_lines,
    clippy::too_many_arguments,
    clippy::similar_names,
    clippy::manual_let_else,
    clippy::option_if_let_else,
    clippy::collapsible_if,
    clippy::unnecessary_wraps,
    clippy::many_single_char_names,
    clippy::manual_clamp,
    clippy::needless_range_loop,
    clippy::type_complexity,
    clippy::module_name_repetitions,
    clippy::must_use_candidate,
    clippy::missing_const_for_fn,
    clippy::use_self
)]
pub mod netcdf;
pub mod npy;
/// Apache Parquet files with a flat schema.
#[allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    clippy::single_call_fn,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::redundant_else,
    clippy::if_not_else,
    clippy::suboptimal_flops,
    clippy::doc_markdown,
    clippy::needless_pass_by_value,
    clippy::too_many_lines,
    clippy::too_many_arguments,
    clippy::similar_names,
    clippy::manual_let_else,
    clippy::option_if_let_else,
    clippy::collapsible_if,
    clippy::unnecessary_wraps,
    clippy::many_single_char_names,
    clippy::manual_clamp,
    clippy::needless_range_loop,
    clippy::type_complexity,
    clippy::module_name_repetitions,
    clippy::must_use_candidate,
    clippy::missing_const_for_fn,
    clippy::use_self
)]
pub mod parquet;
pub mod plot;
pub mod png;
pub mod pretty;
pub mod raster;
pub mod sim_vis;
pub mod svg;
pub mod typst;
pub mod vtk;
pub mod wav;

pub use colormap::Colormap;
pub use gif::GifAnimation;
pub use gif::GifDisposal;
pub use gif::GifFrame;
pub use gif::GifPlotOptions;
pub use gif::decode_gif;
pub use gif::decode_gif_composited;
pub use gif::encode_gif;
pub use gif::encode_particles_gif;
pub use gif::encode_scalar_field_gif;
pub use gif::encode_series_animation_gif;
pub use gif::read_gif_file;
pub use gif::write_gif_file;
pub use gif::write_particles_gif;
pub use gif::write_scalar_field_gif;
pub use gif::write_series_animation_gif;
pub use png::PngColorType;
pub use png::PngImage;
pub use png::decode_png;
pub use png::encode_png;
pub use png::encode_png_with_metadata;
pub use png::read_png_file;
pub use png::write_png_file;
pub use svg::SvgDocument;
pub use svg::SvgElement;
pub use svg::parse_svg;
pub use svg::read_svg_file;
pub use svg::write_svg_file;
pub use vtk::VtkAttributeData;
pub use vtk::VtkDataset;
pub use vtk::VtkFile;
pub use vtk::VtkFormat;
pub use vtk::encode_vtk;
pub use vtk::parse_vtk;
pub use vtk::read_vtk;
pub use vtk::write_vtk;
pub use wav::WavData;
pub use wav::WavSpec;
pub use wav::decode_wav;
pub use wav::encode_wav;
pub use wav::read_wav_file;
pub use wav::write_signal_wav;
pub use wav::write_wav_file;

use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::BufWriter;
use std::io::Write;
use std::path::Path;

use ndarray::Array2;
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
#[cfg(feature = "npy")]
pub fn write_npy_file<P: AsRef<Path>>(
    filename: P,
    arr: &Array2<f64>,
) -> Result<(), String> {
    let data;
    let values = if let Some(contiguous) = arr.as_slice() {
        contiguous
    } else {
        data = arr.iter().copied().collect::<Vec<_>>();
        &data
    };
    let bytes = npy::encode(arr.shape(), values)?;
    std::fs::write(filename, bytes).map_err(|e| e.to_string())
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
#[cfg(feature = "npy")]
pub fn read_npy_file<P: AsRef<Path>>(filename: P) -> Result<Array2<f64>, String> {
    let bytes = std::fs::read(filename).map_err(|e| e.to_string())?;
    let (shape, values) = npy::decode(&bytes)?;
    let [rows, columns] = shape.as_slice() else {
        return Err(format!(
            "expected a rank-two NPY array, found rank {}",
            shape.len()
        ));
    };
    Array2::from_shape_vec((*rows, *columns), values).map_err(|e| e.to_string())
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
    let mut file = BufWriter::new(File::create(filename).map_err(|e| e.to_string())?);

    for row in arr.outer_iter() {
        let mut values = row.iter();
        if let Some(first) = values.next() {
            write!(file, "{first}").map_err(|e| e.to_string())?;
            for value in values {
                write!(file, ",{value}").map_err(|e| e.to_string())?;
            }
        }
        file.write_all(b"\n").map_err(|e| e.to_string())?;
    }

    file.flush().map_err(|e| e.to_string())
}

/// Reads a 2D `ndarray::Array` from a CSV file.
///
/// # Errors
///
/// This function will return an error if the file cannot be read or parsed.
pub fn read_csv_file<P: AsRef<Path>>(filename: P) -> Result<Array2<f64>, String> {
    let file = File::open(filename).map_err(|e| e.to_string())?;

    let mut reader = BufReader::new(file);
    let mut line = String::new();
    let mut data = Vec::new();
    let mut rows = 0;
    let mut cols = 0;

    while reader.read_line(&mut line).map_err(|e| e.to_string())? != 0 {
        let fields = line.trim();
        if !fields.is_empty() {
            let row_start = data.len();
            for s in fields.split(',') {
                data.push(s.trim().parse::<f64>().map_err(|e| e.to_string())?);
            }
            let row_len = data.len().saturating_sub(row_start);
            if rows == 0 {
                cols = row_len;
            } else if row_len != cols {
                return Err("Inconsistent column count in CSV".to_string());
            }
            rows += 1;
        }
        line.clear();
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
    let mut file = BufWriter::new(File::create(filename).map_err(|e| e.to_string())?);
    serde_json::to_writer_pretty(&mut file, arr).map_err(|e| e.to_string())?;
    file.flush().map_err(|e| e.to_string())
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
