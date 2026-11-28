//! Simulation post-processing and visualization example.
//!
//! Run with `cargo run --locked --example sim_postprocess_vis`.
//!
//! Solves a 2D heat diffusion problem on a structured grid and exports:
//! - Chunked, deflate-compressed HDF5 snapshot archive
//! - VTK Legacy ASCII & Binary structured-points fields
//! - PNG colormap heatmap
//! - Animated GIF of the time evolution
//! - SVG line, histogram, and subplot summaries
//! - WAV sonification of the center temperature time series

use std::f64::consts::PI;

use ndarray::Array2;
use rssn::io::Colormap;
use rssn::io::GifPlotOptions;
use rssn::io::plot::AxisScale;
use rssn::io::plot::plot_histogram;
use rssn::io::plot::plot_scalar_field_png;
use rssn::io::plot::plot_series_2d_scaled;
use rssn::io::plot::render_errorbar_svg;
use rssn::io::plot::render_histogram_svg;
use rssn::io::plot::render_subplots_svg;
use rssn::io::sim_vis::animate_grid_fields_gif;
use rssn::io::sim_vis::export_grid_field_vtk;
use rssn::io::sim_vis::export_grid_snapshots_hdf5;
use rssn::io::sim_vis::plot_grid_field_svg;
use rssn::io::vtk::VtkAttributeData;
use rssn::io::vtk::VtkDataset;
use rssn::io::vtk::VtkFile;
use rssn::io::vtk::VtkFormat;
use rssn::io::vtk::write_vtk;
use rssn::io::write_signal_wav;
use rssn::sim::boundary::FaceBcs;
use rssn::sim::structured::Grid;

fn main() -> Result<(), String> {
    let out_dir = std::env::temp_dir().join("rssn-sim-postprocess-vis");
    std::fs::create_dir_all(&out_dir)
        .map_err(|e| format!("cannot create {}: {e}", out_dir.display()))?;

    let nx = 24;
    let ny = 24;
    let grid = Grid::new(
        &[nx, ny],
        &[0.0, 0.0],
        &[1.0, 1.0],
        &FaceBcs::dirichlet_zero(),
    )?;

    let mut snapshots = Vec::new();
    let mut times = Vec::new();
    let mut center_history = Vec::new();

    for step in 0..8 {
        let t = (step as f64) * 0.05;
        times.push(t);
        let decay = (-2.0 * PI * PI * 0.25 * t).exp();
        let mut field = Vec::with_capacity(nx * ny);
        for i in 0..nx {
            let x = (i as f64) / ((nx - 1) as f64);
            for j in 0..ny {
                let y = (j as f64) / ((ny - 1) as f64);
                let val = (PI * x).sin() * (PI * y).sin() * decay;
                field.push(val);
            }
        }
        center_history.push(field[(nx / 2) * ny + (ny / 2)]);
        snapshots.push(field);
    }

    // 1. HDF5 chunked + deflate snapshot archive
    export_grid_snapshots_hdf5(&grid, &snapshots, &times, out_dir.join("heat_snapshots.h5"))?;

    // 2. VTK ASCII & Binary exports
    export_grid_field_vtk(
        &grid,
        &snapshots[0],
        "temperature",
        out_dir.join("heat_t0_ascii.vtk"),
    )?;
    write_vtk(
        out_dir.join("heat_t0_binary.vtk"),
        &VtkFile {
            title: "2D Heat Diffusion Binary".to_owned(),
            format: VtkFormat::Binary,
            dataset: VtkDataset::StructuredPoints {
                dims: [nx, ny, 1],
                origin: [0.0, 0.0, 0.0],
                spacing: [grid.h(0), grid.h(1), 1.0],
            },
            point_data: vec![VtkAttributeData::Scalars {
                name: "temperature".to_owned(),
                values: snapshots[0].clone(),
            }],
            cell_data: Vec::new(),
        },
    )?;

    // 3. SVG & PNG static plots
    plot_grid_field_svg(&grid, &snapshots[0], out_dir.join("heat_t0.svg"))?;
    let arr0 = Array2::from_shape_vec((nx, ny), snapshots[0].clone()).map_err(|e| e.to_string())?;
    plot_scalar_field_png(&arr0, Colormap::Inferno, &out_dir.join("heat_t0.png"))?;

    // 4. Histogram, log-scale decay curve, and subplots
    plot_histogram(&snapshots[0], 10, &out_dir.join("heat_hist.svg"))?;
    let decay_pts: Vec<(f64, f64)> = times
        .iter()
        .copied()
        .zip(center_history.iter().copied())
        .collect();
    plot_series_2d_scaled(
        &[("center_T".to_owned(), decay_pts)],
        AxisScale::Linear,
        AxisScale::Log10,
        &out_dir.join("heat_decay_log.svg"),
    )?;
    let p1 = render_histogram_svg(&snapshots[0], 8)?;
    let p2 = render_errorbar_svg(&times, &center_history, &vec![0.02; times.len()])?;
    let subplots = render_subplots_svg(&[&p1, &p2], 1, 2)?;
    std::fs::write(out_dir.join("heat_summary.svg"), subplots).map_err(|e| e.to_string())?;

    // 5. Animated GIF of time evolution
    animate_grid_fields_gif(
        &grid,
        &snapshots,
        &GifPlotOptions {
            colormap: Colormap::Viridis,
            scale: 4,
            delay_cs: 15,
            loop_count: 0,
            fixed_range: Some((0.0, 1.0)),
        },
        out_dir.join("heat_evolution.gif"),
    )?;

    // 6. WAV sonification of the center temperature
    let audio: Vec<f64> = (0..4000)
        .map(|i| {
            let s = (i as f64) / 8000.0;
            let env = (-3.0 * s).exp();
            env * (2.0 * PI * 440.0 * s).sin()
        })
        .collect();
    write_signal_wav(out_dir.join("heat_sonification.wav"), &audio, 8000)?;

    println!(
        "Wrote simulation post-processing artifacts to {}",
        out_dir.display()
    );
    Ok(())
}
