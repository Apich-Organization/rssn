//! Media, VTK, Plot, and Simulation Visualization integration tests.

use std::path::PathBuf;

use ndarray::Array2;
use rssn::io::plot::*;
use rssn::io::sim_vis::*;
use rssn::io::vtk::*;
use rssn::io::*;
use rssn::sim::boundary::FaceBcs;
use rssn::sim::structured::Grid;

#[test]
fn png_gif_wav_vtk_and_sim_vis_round_trip() {
    let dir = std::env::temp_dir().join(format!("rssn_media_vis_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();

    // 1. PNG via plot_scalar_field_png
    let field = Array2::from_shape_fn((12, 16), |(r, c)| ((r * c) as f64) / 100.0);
    let png_path = dir.join("field.png");
    plot_scalar_field_png(&field, Colormap::Plasma, &png_path).unwrap();
    let img = read_png_file(&png_path).unwrap();
    assert_eq!((img.width, img.height), (16, 12));

    // 2. Dynamic GIFs (scalar field, 1D series, particles)
    let frames: Vec<Array2<f64>> = (0..4)
        .map(|t| Array2::from_shape_fn((8, 8), |(r, c)| ((r + c + t) as f64) * 0.1))
        .collect();
    let gif_path = dir.join("anim.gif");
    write_scalar_field_gif(
        &gif_path,
        &frames,
        &GifPlotOptions {
            scale: 2,
            ..Default::default()
        },
    )
    .unwrap();
    let anim = read_gif_file(&gif_path).unwrap();
    assert_eq!(anim.frames.len(), 4);

    let series_frames: Vec<Vec<(f64, f64)>> = (0..3)
        .map(|t| {
            (0..20)
                .map(|i| (i as f64, ((i + t) as f64).sin()))
                .collect()
        })
        .collect();
    let sgif = encode_series_animation_gif(&series_frames, 120, 80, 10).unwrap();
    assert_eq!(decode_gif(&sgif).unwrap().frames.len(), 3);

    let pframes: Vec<Vec<[f64; 2]>> = (0..3)
        .map(|t| vec![[t as f64 * 0.2, 0.5], [-0.3, t as f64 * 0.1]])
        .collect();
    let pgif = encode_particles_gif(&pframes, 100, 100, 10, 2).unwrap();
    assert_eq!(decode_gif(&pgif).unwrap().frames.len(), 3);

    // 3. WAV audio
    let wav_path = dir.join("tone.wav");
    let signal: Vec<f64> = (0..1600)
        .map(|i| (2.0 * std::f64::consts::PI * 440.0 * (i as f64) / 16000.0).sin())
        .collect();
    write_signal_wav(&wav_path, &signal, 16_000).unwrap();
    let wav = read_wav_file(&wav_path).unwrap();
    assert_eq!(wav.spec.sample_rate, 16_000);
    assert_eq!(wav.samples.len(), 1600);

    // 4. VTK export + read back
    let vtk_path = dir.join("field.vtk");
    write_vtk_structured_points_2d(&vtk_path, &field, "u", [0.1, 0.1], [0.0, 0.0]).unwrap();
    let vtk_back = read_vtk_structured_points_2d(&vtk_path).unwrap();
    assert_eq!(vtk_back.dim(), (12, 16));

    // 5. New SVG plots (histogram, errorbar, log10, subplots)
    let h_svg = render_histogram_svg(&[1.0, 1.2, 2.3, 2.5, 2.7, 3.9], 4).unwrap();
    let e_svg = render_errorbar_svg(&[1.0, 2.0, 3.0], &[2.0, 4.0, 6.0], &[0.2, 0.3, 0.1]).unwrap();
    let l_svg = render_series_2d_scaled_svg(
        &[(
            "exp".to_owned(),
            vec![(1.0, 10.0), (2.0, 100.0), (3.0, 1000.0)],
        )],
        AxisScale::Linear,
        AxisScale::Log10,
    )
    .unwrap();
    let sub_svg = render_subplots_svg(&[&h_svg, &e_svg, &l_svg], 1, 3).unwrap();
    assert!(sub_svg.contains("<g transform="));

    // 6. SIM visualization (`Grid` + `visualize_scenario_json`)
    let grid = Grid::new(
        &[8, 8],
        &[0.0, 0.0],
        &[1.0, 1.0],
        &FaceBcs::dirichlet_zero(),
    )
    .unwrap();
    let snap0: Vec<f64> = (0..64).map(|i| (i as f64) * 0.01).collect();
    let snap1: Vec<f64> = (0..64).map(|i| (i as f64) * 0.02).collect();
    plot_grid_field_svg(&grid, &snap0, dir.join("grid.svg")).unwrap();
    plot_grid_field_png(&grid, &snap0, Colormap::Turbo, dir.join("grid.png")).unwrap();
    export_grid_field_vtk(&grid, &snap0, "T", dir.join("grid.vtk")).unwrap();
    export_grid_snapshots_hdf5(&grid, &[snap0, snap1], &[0.0, 0.1], dir.join("grid.h5")).unwrap();

    let scenario_json = serde_json::json!({
        "energy": [1.0, 0.9, 0.85, 0.82],
        "final_field": [[0.0, 0.5], [0.5, 1.0]]
    });
    let generated = visualize_scenario_json(&scenario_json, &dir, "sim").unwrap();
    assert_eq!(generated.len(), 2);

    let _ = std::fs::remove_dir_all(&dir);

    if let Ok(vdir) = std::env::var("RSSN_VERIFY_DIR") {
        let vdir = PathBuf::from(vdir);
        std::fs::create_dir_all(&vdir).unwrap();
        plot_scalar_field_png(&field, Colormap::Viridis, &vdir.join("rssn_field.png")).unwrap();
        write_scalar_field_gif(
            &vdir.join("rssn_anim.gif"),
            &frames,
            &GifPlotOptions {
                scale: 2,
                ..Default::default()
            },
        )
        .unwrap();
        write_signal_wav(&vdir.join("rssn_tone.wav"), &signal, 16_000).unwrap();
        write_vtk(
            &vdir.join("rssn_bin.vtk"),
            &VtkFile {
                title: "Verify Binary".to_owned(),
                format: VtkFormat::Binary,
                dataset: VtkDataset::StructuredPoints {
                    dims: [3, 2, 1],
                    origin: [0.0, 0.0, 0.0],
                    spacing: [1.0, 1.0, 1.0],
                },
                point_data: vec![VtkAttributeData::Scalars {
                    name: "u".to_owned(),
                    values: vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
                }],
                cell_data: Vec::new(),
            },
        )
        .unwrap();
    }
}

#[test]
fn golden_media_and_vtk_fixtures_decode_strictly() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let media = root.join("media");
    let vtk_dir = root.join("vtk");

    // Golden PNGs
    assert_eq!(
        read_png_file(media.join("golden_rgb8.png"))
            .unwrap()
            .color_type,
        PngColorType::Rgb8
    );
    assert_eq!(
        read_png_file(media.join("golden_rgba8.png"))
            .unwrap()
            .color_type,
        PngColorType::Rgba8
    );
    assert_eq!(
        read_png_file(media.join("golden_gray16.png"))
            .unwrap()
            .color_type,
        PngColorType::Grayscale16
    );
    assert_eq!(
        read_png_file(media.join("golden_indexed.png"))
            .unwrap()
            .width,
        6
    );

    // Golden animated GIF (raw + composited)
    let gif_bytes = std::fs::read(media.join("golden_anim.gif")).unwrap();
    let gif_raw = decode_gif(&gif_bytes).unwrap();
    assert_eq!(gif_raw.frames.len(), 3);
    let gif_comp = decode_gif_composited(&gif_bytes).unwrap();
    assert_eq!(gif_comp.frames.len(), 3);

    // Golden WAVs (8, 16 stereo, 24, f32, f64)
    assert_eq!(
        read_wav_file(media.join("golden_pcm8.wav"))
            .unwrap()
            .samples
            .len(),
        4
    );
    let w16 = read_wav_file(media.join("golden_pcm16_stereo.wav")).unwrap();
    assert_eq!(w16.spec.channels, 2);
    assert_eq!(w16.samples.len(), 4);
    assert_eq!(
        read_wav_file(media.join("golden_pcm24.wav"))
            .unwrap()
            .spec
            .bits_per_sample,
        24
    );
    assert!(
        read_wav_file(media.join("golden_f32.wav"))
            .unwrap()
            .spec
            .is_float
    );
    assert_eq!(
        read_wav_file(media.join("golden_f64.wav"))
            .unwrap()
            .spec
            .bits_per_sample,
        64
    );

    // Golden VTK (ASCII + Binary)
    let va = read_vtk(vtk_dir.join("golden_structured_ascii.vtk")).unwrap();
    let vb = read_vtk(vtk_dir.join("golden_structured_binary.vtk")).unwrap();
    assert_eq!(va.dataset, vb.dataset);
    assert_eq!(va.point_data, vb.point_data);
}
