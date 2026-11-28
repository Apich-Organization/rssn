//! Ergonomic post-processing and visualization helpers for RSSN simulations.
//!
//! Bridges [`crate::sim::structured::Grid`], particle snapshots, and
//! [`crate::sim::scenario::run`] JSON outputs directly to SVG, PNG, animated GIF,
//! VTK, and compressed HDF5 archives.

use std::path::Path;
use std::path::PathBuf;

use ndarray::Array2;
use serde_json::Value;

use super::colormap::Colormap;
use super::gif::GifPlotOptions;
use super::gif::write_particles_gif;
use super::gif::write_scalar_field_gif;
use super::hdf5::AttrValue;
use super::hdf5::Attribute;
use super::hdf5::Dataset;
use super::hdf5::DatasetWriteOptions;
use super::hdf5::H5Filter;
use super::hdf5::H5Type;
use super::hdf5::H5WriteOptions;
use super::hdf5::write_hdf5_with_options;
use super::plot::plot_heatmap_2d;
use super::plot::plot_scalar_field_png;
use super::plot::plot_series_2d;
use super::vtk::write_vtk_polydata_points;
use super::vtk::write_vtk_structured_points_2d;
use super::vtk::write_vtk_structured_points_3d;
use crate::sim::structured::Grid;

fn grid_to_2d_array(
    grid: &Grid,
    field: &[f64],
) -> Result<Array2<f64>, String> {
    match grid.dim {
        | 2 => {
            let (nx, ny) = (grid.n[0], grid.n[1]);
            if field.len() != nx * ny {
                return Err(format!(
                    "expected {} values for 2D grid {nx}x{ny}, got {}",
                    nx * ny,
                    field.len()
                ));
            }
            Array2::from_shape_vec((nx, ny), field.to_vec()).map_err(|e| e.to_string())
        },
        | 3 => {
            let (nx, ny, nz) = (grid.n[0], grid.n[1], grid.n[2]);
            if field.len() != nx * ny * nz {
                return Err(format!(
                    "expected {} values for 3D grid {nx}x{ny}x{nz}, got {}",
                    nx * ny * nz,
                    field.len()
                ));
            }
            let mid = nz / 2;
            Ok(Array2::from_shape_fn((nx, ny), |(i, j)| {
                field[i * (ny * nz) + j * nz + mid]
            }))
        },
        | other => Err(format!("expected 2D or 3D grid, got {other}D")),
    }
}

/// Plots a [`Grid`] scalar field to SVG (1D curve, 2D heatmap, or 3D mid-plane slice).
///
/// # Errors
/// Fails when dimensions/length mismatch or writing fails.
pub fn plot_grid_field_svg<P: AsRef<Path>>(
    grid: &Grid,
    field: &[f64],
    path: P,
) -> Result<(), String> {
    match grid.dim {
        | 1 => {
            let n = grid.n[0];
            if field.len() != n {
                return Err(format!(
                    "expected {n} values for 1D grid, got {}",
                    field.len()
                ));
            }
            let pts: Vec<(f64, f64)> = (0..n).map(|i| (grid.coord(i)[0], field[i])).collect();
            plot_series_2d(&[("field".to_owned(), pts)], path.as_ref())
        },
        | 2 | 3 => {
            let arr = grid_to_2d_array(grid, field)?;
            plot_heatmap_2d(&arr, path.as_ref())
        },
        | d => Err(format!("unsupported grid dimension {d}")),
    }
}

/// Plots a 2D or 3D-midplane [`Grid`] scalar field as a PNG image.
///
/// # Errors
/// Fails for 1D grids, mismatched lengths, or writing errors.
pub fn plot_grid_field_png<P: AsRef<Path>>(
    grid: &Grid,
    field: &[f64],
    colormap: Colormap,
    path: P,
) -> Result<(), String> {
    let arr = grid_to_2d_array(grid, field)?;
    plot_scalar_field_png(&arr, colormap, path.as_ref())
}

/// Exports a 2D or 3D [`Grid`] scalar field to VTK Legacy ASCII (`STRUCTURED_POINTS`).
///
/// # Errors
/// Fails for 1D grids, mismatched lengths, or writing errors.
pub fn export_grid_field_vtk<P: AsRef<Path>>(
    grid: &Grid,
    field: &[f64],
    field_name: &str,
    path: P,
) -> Result<(), String> {
    match grid.dim {
        | 2 => {
            let arr = grid_to_2d_array(grid, field)?;
            write_vtk_structured_points_2d(
                path,
                &arr,
                field_name,
                [grid.h(0), grid.h(1)],
                [grid.lo[0], grid.lo[1]],
            )
        },
        | 3 => {
            write_vtk_structured_points_3d(
                path,
                grid.n,
                field,
                field_name,
                [grid.h(0), grid.h(1), grid.h(2)],
                grid.lo,
            )
        },
        | d => Err(format!("VTK grid export supports 2D/3D grids, got {d}D")),
    }
}

/// Animates a sequence of 2D or 3D-midplane [`Grid`] snapshots into an animated GIF.
///
/// # Errors
/// Fails when `snapshots` is empty, dimensions mismatch, or GIF encoding fails.
pub fn animate_grid_fields_gif<P: AsRef<Path>>(
    grid: &Grid,
    snapshots: &[Vec<f64>],
    opts: &GifPlotOptions,
    path: P,
) -> Result<(), String> {
    let mut arrays = Vec::with_capacity(snapshots.len());
    for snap in snapshots {
        arrays.push(grid_to_2d_array(grid, snap)?);
    }
    write_scalar_field_gif(path, &arrays, opts)
}

/// Exports [`Grid`] time-series snapshots to a chunked, deflate-compressed HDF5 archive.
///
/// # Errors
/// Fails when `snapshots` is empty, lengths mismatch, or HDF5 writing fails.
pub fn export_grid_snapshots_hdf5<P: AsRef<Path>>(
    grid: &Grid,
    snapshots: &[Vec<f64>],
    times: &[f64],
    path: P,
) -> Result<(), String> {
    if snapshots.is_empty() || snapshots.len() != times.len() {
        return Err("snapshots and times must be non-empty and equal length".to_owned());
    }
    let grid_len: usize = grid.n[..grid.dim].iter().product();
    let mut flat = Vec::with_capacity(snapshots.len() * grid_len);
    for snap in snapshots {
        if snap.len() != grid_len {
            return Err(format!(
                "snapshot length {} does not match grid size {grid_len}",
                snap.len()
            ));
        }
        flat.extend_from_slice(snap);
    }
    let mut shape = vec![snapshots.len()];
    shape.extend_from_slice(&grid.n[..grid.dim]);
    let mut chunk_dims = vec![1];
    chunk_dims.extend_from_slice(&grid.n[..grid.dim]);

    let datasets = vec![
        Dataset::new_numeric("/snapshots", shape, H5Type::F64, flat),
        Dataset::new_numeric("/times", vec![times.len()], H5Type::F64, times.to_vec()),
    ];
    let attrs = vec![(
        "/".to_owned(),
        vec![
            Attribute {
                name: "dim".to_owned(),
                value: AttrValue::Int(grid.dim as i64),
            },
            Attribute {
                name: "lo".to_owned(),
                value: AttrValue::Floats(grid.lo[..grid.dim].to_vec()),
            },
            Attribute {
                name: "hi".to_owned(),
                value: AttrValue::Floats(grid.hi[..grid.dim].to_vec()),
            },
        ],
    )];
    let options = H5WriteOptions {
        dense_groups: false,
        dataset_options: vec![(
            "/snapshots".to_owned(),
            DatasetWriteOptions {
                chunk_dims: Some(chunk_dims),
                filters: vec![H5Filter::Shuffle, H5Filter::Deflate(6)],
            },
        )],
        soft_links: Vec::new(),
    };
    write_hdf5_with_options(path.as_ref(), &datasets, &attrs, &options)
}

/// Visualizes 2D particle positions as an SVG scatter/trajectory plot.
///
/// # Errors
/// Fails when `positions` is empty or writing fails.
pub fn visualize_particles_svg<P: AsRef<Path>>(
    positions: &[[f64; 2]],
    path: P,
) -> Result<(), String> {
    if positions.is_empty() {
        return Err("particle positions cannot be empty".to_owned());
    }
    let pts: Vec<(f64, f64)> = positions.iter().map(|p| (p[0], p[1])).collect();
    plot_series_2d(&[("particles".to_owned(), pts)], path.as_ref())
}

/// Animates 2D particle trajectories into a GIF (re-exported convenience wrapper).
///
/// # Errors
/// As [`write_particles_gif`].
pub fn animate_particles_gif<P: AsRef<Path>>(
    frames: &[Vec<[f64; 2]>],
    width: u16,
    height: u16,
    delay_cs: u16,
    trail_len: usize,
    path: P,
) -> Result<(), String> {
    write_particles_gif(path, frames, width, height, delay_cs, trail_len)
}

/// Exports 3D particle positions to VTK `POLYDATA`.
///
/// # Errors
/// As [`write_vtk_polydata_points`].
pub fn export_particles_vtk<P: AsRef<Path>>(
    positions: &[[f64; 3]],
    scalars: Option<(&str, &[f64])>,
    path: P,
) -> Result<(), String> {
    write_vtk_polydata_points(path, positions, scalars)
}

/// Writes visualizations for a JSON simulation result.
///
/// The top-level object returned by [`crate::sim::scenario::run`] may contain 1D numeric
/// vectors, 2D matrices, and 2D time-series snapshots. These produce SVG, PNG, and GIF files,
/// respectively.
///
/// Returns the list of generated file paths.
///
/// # Errors
/// Fails if `out_dir` cannot be created or any visualization write fails.
pub fn visualize_scenario_json(
    result_json: &Value,
    out_dir: &Path,
    prefix: &str,
) -> Result<Vec<PathBuf>, String> {
    std::fs::create_dir_all(out_dir)
        .map_err(|e| format!("cannot create {}: {e}", out_dir.display()))?;
    let mut written = Vec::new();
    let Some(obj) = result_json.as_object() else {
        return Ok(written);
    };

    for (key, val) in obj {
        if let Some(arr) = val.as_array() {
            // Check 1D numeric array
            if !arr.is_empty() && arr.iter().all(Value::is_number) {
                let pts: Vec<(f64, f64)> = arr
                    .iter()
                    .enumerate()
                    .filter_map(|(i, v)| v.as_f64().map(|y| (i as f64, y)))
                    .collect();
                if !pts.is_empty() {
                    let p = out_dir.join(format!("{prefix}_{key}.svg"));
                    plot_series_2d(&[(key.clone(), pts)], &p)?;
                    written.push(p);
                }
                continue;
            }
            // Check 2D numeric matrix
            if !arr.is_empty()
                && arr.iter().all(|row| {
                    row.as_array()
                        .is_some_and(|r| !r.is_empty() && r.iter().all(Value::is_number))
                })
            {
                let rows = arr.len();
                let cols = arr[0].as_array().map_or(0, Vec::len);
                if cols > 0 && arr.iter().all(|r| r.as_array().map_or(0, Vec::len) == cols) {
                    let flat: Vec<f64> = arr
                        .iter()
                        .filter_map(Value::as_array)
                        .flatten()
                        .filter_map(Value::as_f64)
                        .collect();
                    if flat.len() == rows * cols {
                        let m = Array2::from_shape_vec((rows, cols), flat)
                            .map_err(|e| e.to_string())?;
                        let png_p = out_dir.join(format!("{prefix}_{key}.png"));
                        plot_scalar_field_png(&m, Colormap::Viridis, &png_p)?;
                        written.push(png_p);
                    }
                }
                continue;
            }
            // Check 3D array (time-series of 2D matrices) -> animated GIF
            if !arr.is_empty()
                && arr.iter().all(|frame| {
                    frame.as_array().is_some_and(|rows| {
                        !rows.is_empty()
                            && rows.iter().all(|r| {
                                r.as_array().is_some_and(|c| {
                                    !c.is_empty() && c.iter().all(Value::is_number)
                                })
                            })
                    })
                })
            {
                let mut frames = Vec::with_capacity(arr.len());
                for r_arr in arr.iter().filter_map(Value::as_array) {
                    let rows = r_arr.len();
                    let cols = r_arr[0].as_array().map_or(0, Vec::len);
                    let flat: Vec<f64> = r_arr
                        .iter()
                        .filter_map(Value::as_array)
                        .flatten()
                        .filter_map(Value::as_f64)
                        .collect();
                    if flat.len() == rows * cols {
                        frames.push(
                            Array2::from_shape_vec((rows, cols), flat)
                                .map_err(|e| e.to_string())?,
                        );
                    }
                }
                if !frames.is_empty() {
                    let gif_p = out_dir.join(format!("{prefix}_{key}.gif"));
                    write_scalar_field_gif(&gif_p, &frames, &GifPlotOptions::default())?;
                    written.push(gif_p);
                }
            }
        }
    }

    Ok(written)
}
