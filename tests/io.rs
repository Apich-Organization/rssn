//! Array I/O in `rssn::io` (CSV, JSON and, with the `npy` feature, NumPy files).

use assert_approx_eq::assert_approx_eq;
use ndarray::Array2;
use ndarray::arr2;
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::graph::Graph;
use rssn::graph::Payload;
use rssn::io::*;
use std::path::PathBuf;
use std::sync::Arc;

fn cfg() -> ProptestConfig {
    ProptestConfig {
        cases: 24,
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

/// A unique file in the system temp directory that is removed on drop.
struct Tmp(PathBuf);

impl Tmp {
    fn new(
        tag: &str,
        ext: &str,
    ) -> Self {
        use std::sync::atomic::AtomicUsize;
        use std::sync::atomic::Ordering;
        static N: AtomicUsize = AtomicUsize::new(0);
        let n = N.fetch_add(1, Ordering::SeqCst);
        Self(std::env::temp_dir().join(format!("rssn_io_{}_{n}_{tag}.{ext}", std::process::id())))
    }
}

impl Drop for Tmp {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

fn read_svg(f: &Tmp) -> String {
    std::fs::read_to_string(&f.0).unwrap_or_else(|e| panic!("{e}"))
}

fn sample() -> Array2<f64> {
    arr2(&[[1.0, -2.5, 3.25], [0.0, 1e-9, 1e12]])
}

#[test]
fn csv_round_trip_preserves_values_and_shape() {
    let f = Tmp::new("rt", "csv");
    write_csv_file(&f.0, &sample()).unwrap_or_else(|e| panic!("{e}"));
    let back = read_csv_file(&f.0).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(back, sample());
}

#[test]
fn csv_is_written_one_row_per_line_comma_separated() {
    let f = Tmp::new("fmt", "csv");
    write_csv_file(&f.0, &arr2(&[[1.0, 2.0], [3.5, 4.0]])).unwrap_or_else(|e| panic!("{e}"));
    let text = std::fs::read_to_string(&f.0).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(text, "1,2\n3.5,4\n");
}

#[test]
fn csv_reader_skips_blank_lines_and_trims_whitespace() {
    let f = Tmp::new("blank", "csv");
    std::fs::write(&f.0, "1, 2 ,3\n\n  \n4,5,6\n").unwrap_or_else(|e| panic!("{e}"));
    let a = read_csv_file(&f.0).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(a, arr2(&[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]));
}

#[test]
fn csv_reader_rejects_ragged_rows_and_garbage() {
    let f = Tmp::new("ragged", "csv");
    std::fs::write(&f.0, "1,2,3\n4,5\n").unwrap_or_else(|e| panic!("{e}"));
    let err = read_csv_file(&f.0).err().unwrap_or_default();
    assert!(err.contains("Inconsistent column count"), "{err}");
    std::fs::write(&f.0, "1,abc\n").unwrap_or_else(|e| panic!("{e}"));
    assert!(read_csv_file(&f.0).is_err());
}

#[test]
fn csv_empty_file_gives_an_empty_array() {
    let f = Tmp::new("empty", "csv");
    std::fs::write(&f.0, "").unwrap_or_else(|e| panic!("{e}"));
    let a = read_csv_file(&f.0).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(a.shape(), &[0, 0]);
}

#[test]
fn missing_files_are_reported_as_errors() {
    let missing = std::env::temp_dir().join("rssn_io_definitely_missing_file");
    assert!(read_csv_file(&missing).is_err());
    assert!(read_json_file(&missing).is_err());
    assert!(read_npy_file(&missing).is_err());
}

#[test]
fn json_round_trip_preserves_values_and_shape() {
    let f = Tmp::new("rt", "json");
    write_json_file(&f.0, &sample()).unwrap_or_else(|e| panic!("{e}"));
    let back = read_json_file(&f.0).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(back, sample());
}

#[test]
fn json_reader_rejects_malformed_input() {
    let f = Tmp::new("bad", "json");
    std::fs::write(&f.0, "{ not json").unwrap_or_else(|e| panic!("{e}"));
    assert!(read_json_file(&f.0).is_err());
}

#[test]
fn writing_into_a_missing_directory_fails() {
    let bad = std::env::temp_dir().join("rssn_io_no_such_dir").join("x");
    assert!(write_csv_file(&bad, &sample()).is_err());
    assert!(write_json_file(&bad, &sample()).is_err());
    assert!(write_npy_file(&bad, &sample()).is_err());
}

#[cfg(feature = "npy")]
#[test]
fn npy_round_trip_preserves_values_and_shape() {
    let f = Tmp::new("rt", "npy");
    write_npy_file(&f.0, &sample()).unwrap_or_else(|e| panic!("{e}"));
    let back = read_npy_file(&f.0).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(back, sample());
}

#[cfg(feature = "npy")]
#[test]
fn npy_file_has_the_numpy_magic_header() {
    let f = Tmp::new("magic", "npy");
    write_npy_file(&f.0, &arr2(&[[1.0, 2.0]])).unwrap_or_else(|e| panic!("{e}"));
    let bytes = std::fs::read(&f.0).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(&bytes[..6], b"\x93NUMPY");
}

#[cfg(feature = "npy")]
#[test]
fn npy_reader_rejects_non_npy_content() {
    let f = Tmp::new("bad", "npy");
    std::fs::write(&f.0, "this is not a numpy file").unwrap_or_else(|e| panic!("{e}"));
    assert!(read_npy_file(&f.0).is_err());
}

#[cfg(not(feature = "npy"))]
#[test]
fn npy_functions_report_the_missing_feature() {
    let f = Tmp::new("nofeat", "npy");
    let w = write_npy_file(&f.0, &sample()).err().unwrap_or_default();
    let r = read_npy_file(&f.0).err().unwrap_or_default();
    assert!(w.contains("npy") && r.contains("npy"), "{w} / {r}");
}

#[test]
fn all_formats_agree_on_the_same_data() {
    let a = arr2(&[[0.5, 1.5], [2.5, 3.5], [4.5, 5.5]]);
    let (c, j) = (Tmp::new("all", "csv"), Tmp::new("all", "json"));
    write_csv_file(&c.0, &a).unwrap_or_else(|e| panic!("{e}"));
    write_json_file(&j.0, &a).unwrap_or_else(|e| panic!("{e}"));
    let (bc, bj) = (
        read_csv_file(&c.0).unwrap_or_else(|e| panic!("{e}")),
        read_json_file(&j.0).unwrap_or_else(|e| panic!("{e}")),
    );
    for (x, y) in bc.iter().zip(bj.iter()) {
        assert_approx_eq!(x, y, 1e-15);
    }
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_csv_round_trip(rows in 1usize..6, cols in 1usize..6, seed in any::<u64>()) {
        let mut s = seed;
        let data: Vec<f64> = (0..rows * cols).map(|_| {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (s >> 11) as f64 / (1u64 << 40) as f64 - 4.0
        }).collect();
        let a = Array2::from_shape_vec((rows, cols), data).unwrap_or_else(|e| panic!("{e}"));
        let f = Tmp::new("prop", "csv");
        write_csv_file(&f.0, &a).unwrap_or_else(|e| panic!("{e}"));
        prop_assert_eq!(read_csv_file(&f.0).unwrap_or_else(|e| panic!("{e}")), a);
    }

    #[test]
    fn prop_json_round_trip(rows in 1usize..6, cols in 1usize..6, x in -1e6f64..1e6) {
        let a = Array2::from_shape_fn((rows, cols), |(i, j)| x * (i as f64 + 1.0) - j as f64);
        let f = Tmp::new("prop", "json");
        write_json_file(&f.0, &a).unwrap_or_else(|e| panic!("{e}"));
        // serde_json without `float_roundtrip` may be off by one ulp when parsing.
        let back = read_json_file(&f.0).unwrap_or_else(|e| panic!("{e}"));
        prop_assert_eq!(back.shape(), a.shape());
        for (u, v) in back.iter().zip(a.iter()) {
            prop_assert!((u - v).abs() <= 1e-14 * v.abs().max(1.0));
        }
    }
}

// ---- Contour export (IO-1) ----

#[test]
fn contour_of_a_circle_is_one_closed_polyline() {
    let lines = plot::contour_lines_2d(
        |x, y| x * x + y * y - 1.0,
        (-2.0, 2.0),
        (-2.0, 2.0),
        81,
        &[0.0],
    )
    .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(lines.len(), 1, "a circle has one contour component");
    assert_eq!(lines[0].level.to_bits(), 0.0f64.to_bits());
    let points = &lines[0].points;
    assert!(
        points.len() > 100,
        "only {} points on the circle",
        points.len()
    );
    assert_eq!(
        points.first(),
        points.last(),
        "the contour should be closed"
    );
    for &(x, y) in points {
        let r = x.hypot(y);
        assert!((r - 1.0).abs() < 0.01, "point ({x}, {y}) is at radius {r}");
    }
    // Deterministic: the same request produces the same segments.
    let again = plot::contour_lines_2d(
        |x, y| x * x + y * y - 1.0,
        (-2.0, 2.0),
        (-2.0, 2.0),
        81,
        &[0.0],
    )
    .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(lines, again);
}

#[test]
fn contour_of_a_cone_is_the_same_unit_circle() {
    let lines = plot::contour_lines_2d(|x, y| x.hypot(y), (-2.0, 2.0), (-2.0, 2.0), 81, &[1.0])
        .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(lines.len(), 1);
    assert!(lines[0].points.len() > 100);
    for &(x, y) in &lines[0].points {
        assert!(
            (x.hypot(y) - 1.0).abs() < 0.01,
            "({x}, {y}) off the cone level set"
        );
    }
}

#[test]
fn contour_saddle_cells_are_resolved_deterministically() {
    // One 2x2 cell of `x*y`: the corner mean (the bilinear saddle value) is
    // exactly the level, so the saddle decider resolves the cell
    // deterministically into two segments instead of four or none.
    let saddle = plot::contour_lines_2d(|x, y| x * y, (-1.0, 1.0), (-1.0, 1.0), 2, &[0.0])
        .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(saddle.len(), 2);
    for line in &saddle {
        assert_eq!(line.points.len(), 2);
    }
    let again = plot::contour_lines_2d(|x, y| x * y, (-1.0, 1.0), (-1.0, 1.0), 2, &[0.0])
        .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(saddle, again);

    // `x*y = 0.37` has two disconnected hyperbola branches.
    let level = 0.37;
    let branches = plot::contour_lines_2d(|x, y| x * y, (-1.0, 1.0), (-1.0, 1.0), 9, &[level])
        .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(branches.len(), 2, "each hyperbola branch is one polyline");
    let on_boundary =
        |p: &(f64, f64)| (p.0.abs() - 1.0).abs() < 1e-12 || (p.1.abs() - 1.0).abs() < 1e-12;
    for line in &branches {
        for &(x, y) in &line.points {
            assert!(
                (x * y - level).abs() < 1e-9,
                "({x}, {y}) is not on x*y = {level}"
            );
        }
        let starts_outside = on_boundary(&line.points[0]);
        let ends_outside = on_boundary(line.points.last().unwrap());
        assert!(
            starts_outside && ends_outside,
            "an open branch must end on the boundary"
        );
    }
}

#[test]
fn contour_levels_outside_the_data_range_are_empty() {
    let lines = plot::contour_lines_2d(
        |x, y| x * x + y * y,
        (-1.0, 1.0),
        (-1.0, 1.0),
        21,
        &[10.0, -10.0],
    )
    .unwrap_or_else(|e| panic!("{e}"));
    assert!(lines.is_empty());
    let f = Tmp::new("contour_empty", "svg");
    plot::plot_contour_2d(
        |x, y| x * x + y * y,
        (-1.0, 1.0),
        (-1.0, 1.0),
        21,
        &[10.0],
        &f.0,
    )
    .unwrap_or_else(|e| panic!("{e}"));
    let svg = read_svg(&f);
    assert!(svg.starts_with("<svg"));
    assert!(svg.trim_end().ends_with("</svg>"));
    assert_eq!(svg.matches("<polyline").count(), 0);
}

#[test]
fn contour_svg_has_one_polyline_per_contour_line() {
    let levels = [0.0, 0.5];
    let lines = plot::contour_lines_2d(
        |x, y| x * x + y * y - 1.0,
        (-2.0, 2.0),
        (-2.0, 2.0),
        61,
        &levels,
    )
    .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(lines.len(), 2, "one circle per level");
    let f = Tmp::new("contour", "svg");
    plot::plot_contour_2d(
        |x, y| x * x + y * y - 1.0,
        (-2.0, 2.0),
        (-2.0, 2.0),
        61,
        &levels,
        &f.0,
    )
    .unwrap_or_else(|e| panic!("{e}"));
    let svg = read_svg(&f);
    assert_eq!(svg.matches("<polyline").count(), lines.len());
}

#[test]
fn contour_rejects_invalid_requests() {
    let f = Tmp::new("contour_bad", "svg");
    assert!(plot::contour_lines_2d(|x, y| x + y, (1.0, 1.0), (0.0, 1.0), 11, &[0.0]).is_err());
    assert!(plot::contour_lines_2d(|x, y| x + y, (0.0, 1.0), (0.0, 1.0), 1, &[0.0]).is_err());
    assert!(plot::contour_lines_2d(|x, y| x + y, (0.0, 1.0), (0.0, 1.0), 11, &[]).is_err());
    assert!(plot::plot_contour_2d(|x, y| x + y, (0.0, 1.0), (0.0, 1.0), 11, &[], &f.0).is_err());
    assert!(!f.0.exists());
}

// ---- Streamline export (IO-1) ----

#[test]
fn circular_field_streamlines_stay_on_the_circle() {
    let paths = plot::streamline_paths_2d(
        |x, y| (-y, x),
        (-2.0, 2.0),
        (-2.0, 2.0),
        &[(1.0, 0.0)],
        0.01,
        2000,
    )
    .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(paths.len(), 1);
    assert_eq!(paths[0].len(), 2001, "seed plus one point per RK4 step");
    for &(x, y) in &paths[0] {
        let r = x.hypot(y);
        assert!(
            (r - 1.0).abs() < 1e-4,
            "point ({x}, {y}) drifted to radius {r}"
        );
    }
    let again = plot::streamline_paths_2d(
        |x, y| (-y, x),
        (-2.0, 2.0),
        (-2.0, 2.0),
        &[(1.0, 0.0)],
        0.01,
        2000,
    )
    .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(paths, again);
}

#[test]
#[allow(clippy::float_cmp)] // the uniform field integrates `y` exactly, so equality is the test
fn uniform_field_streamlines_are_straight_lines() {
    let seeds = [(0.0, -0.5), (0.0, 0.5)];
    let paths =
        plot::streamline_paths_2d(|_, _| (1.0, 0.0), (0.0, 1.0), (-1.0, 1.0), &seeds, 0.1, 100)
            .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(paths.len(), 2);
    for (path, &(_, y0)) in paths.iter().zip(&seeds) {
        assert!(path.len() >= 2);
        for &(_, y) in path {
            assert_eq!(y, y0, "a uniform field cannot bend the line");
        }
        for pair in path.windows(2) {
            assert!(pair[1].0 > pair[0].0, "the path must advance in +x");
        }
        assert!(
            path.last().unwrap().0 <= 1.0 + 1e-12,
            "the path must stay in range"
        );
    }
}

#[test]
fn streamline_svg_has_one_polyline_per_path() {
    let seeds = [(0.0, -0.5), (0.0, 0.0), (0.0, 0.5)];
    let paths =
        plot::streamline_paths_2d(|_, _| (1.0, 0.0), (0.0, 1.0), (-1.0, 1.0), &seeds, 0.1, 100)
            .unwrap_or_else(|e| panic!("{e}"));
    let f = Tmp::new("streamlines", "svg");
    plot::plot_streamlines_2d(
        |_, _| (1.0, 0.0),
        (0.0, 1.0),
        (-1.0, 1.0),
        &seeds,
        0.1,
        100,
        &f.0,
    )
    .unwrap_or_else(|e| panic!("{e}"));
    let svg = read_svg(&f);
    assert!(svg.starts_with("<svg"));
    assert!(svg.trim_end().ends_with("</svg>"));
    assert_eq!(svg.matches("<polyline").count(), paths.len());
}

#[test]
fn streamlines_reject_invalid_requests() {
    let field = |_, _| (1.0, 0.0);
    assert!(
        plot::streamline_paths_2d(field, (0.0, 1.0), (0.0, 1.0), &[(0.0, 0.0)], 0.0, 10).is_err()
    );
    assert!(
        plot::streamline_paths_2d(field, (0.0, 1.0), (0.0, 1.0), &[(0.0, 0.0)], f64::NAN, 10)
            .is_err()
    );
    assert!(
        plot::streamline_paths_2d(field, (0.0, 1.0), (0.0, 1.0), &[(0.0, 0.0)], 0.1, 0).is_err()
    );
    assert!(plot::streamline_paths_2d(field, (0.0, 1.0), (0.0, 1.0), &[], 0.1, 10).is_err());
    assert!(
        plot::streamline_paths_2d(field, (0.0, 1.0), (0.0, 1.0), &[(2.0, 0.0)], 0.1, 10).is_err()
    );
    assert!(
        plot::streamline_paths_2d(field, (1.0, 1.0), (0.0, 1.0), &[(0.0, 0.0)], 0.1, 10).is_err()
    );
    // A vanishing field leaves only the seed point.
    let dead = plot::streamline_paths_2d(
        |_, _| (0.0, 0.0),
        (0.0, 1.0),
        (0.0, 1.0),
        &[(0.5, 0.5)],
        0.1,
        5,
    )
    .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(dead, vec![vec![(0.5, 0.5)]]);
}

// ---- LaTeX/Typst escaping and pretty symbol coverage (IO-1) ----

#[test]
fn latex_text_escapes_tex_special_characters() {
    let mut g = Graph::new();
    let node = g.lit(Payload::Str(Arc::from(r#"a_b $5 & #1 %2 {x} \ ^ ~"#)));
    let out = rssn::io::latex::to_latex(&g, node);
    assert_eq!(
        out,
        r"\text{a\_b \$5 \& \#1 \%2 \{x\} \textbackslash{} \textasciicircum{} \textasciitilde{}}"
    );
    let plain = g.lit(Payload::Str(Arc::from("plain text")));
    assert_eq!(rssn::io::latex::to_latex(&g, plain), r"\text{plain text}");
}

#[test]
fn typst_string_literals_escape_quotes_and_backslashes() {
    let mut g = Graph::new();
    let node = g.lit(Payload::Str(Arc::from(r#"say "hi" \ there"#)));
    let out = rssn::io::typst::to_typst(&g, node);
    assert_eq!(out, r#""say \"hi\" \\ there""#);
    let plain = g.lit(Payload::Str(Arc::from("plain")));
    assert_eq!(rssn::io::typst::to_typst(&g, plain), "\"plain\"");
}

#[test]
fn pretty_prints_lowercase_omicron_and_capital_upsilon() {
    let mut g = Graph::new();
    let omicron = {
        let id = g.interner_mut().symbol("omicron");
        g.lit(Payload::Sym(id))
    };
    assert_eq!(rssn::io::pretty::to_pretty(&g, omicron), "ο");
    let upsilon = {
        let id = g.interner_mut().symbol("Upsilon");
        g.lit(Payload::Sym(id))
    };
    assert_eq!(rssn::io::pretty::to_pretty(&g, upsilon), "Υ");
}

// ---- Filled contour bands and raster export (IO-2) ----

fn read_bytes(f: &Tmp) -> Vec<u8> {
    std::fs::read(&f.0).unwrap_or_else(|e| panic!("{e}"))
}

/// Shoelace area of a closed polygon given without a repeated last vertex.
fn polygon_area(points: &[(f64, f64)]) -> f64 {
    let mut acc = 0.0;
    for i in 0..points.len() {
        let a = points[i];
        let b = points[(i + 1) % points.len()];
        acc += a.0 * b.1 - b.0 * a.1;
    }
    (acc / 2.0).abs()
}

/// Reads a quoted SVG attribute such as `y2="-3.00"` from one element line.
fn svg_attr(
    line: &str,
    name: &str,
) -> String {
    let prefix = format!("{name}=\"");
    let start = line
        .find(&prefix)
        .unwrap_or_else(|| panic!("{name} missing in {line}"))
        + prefix.len();
    let rest = &line[start..];
    let end = rest
        .find('"')
        .unwrap_or_else(|| panic!("unterminated {name} in {line}"));
    rest[..end].to_owned()
}

#[test]
fn contour_bands_partition_a_cone() {
    let bands =
        plot::contour_bands_2d(|x, y| x.hypot(y), (-2.0, 2.0), (-2.0, 2.0), 41, &[1.0, 2.0])
            .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(bands.len(), 3, "one band per gap plus the two outer bands");
    assert!(bands[0].low.is_infinite() && bands[0].low < 0.0);
    assert!(bands[2].high.is_infinite() && bands[2].high > 0.0);
    let expected = [
        std::f64::consts::PI,
        3.0 * std::f64::consts::PI,
        16.0 - 4.0 * std::f64::consts::PI,
    ];
    for (band, want) in bands.iter().zip(expected) {
        let area: f64 = band.polygons.iter().map(|p| polygon_area(p)).sum();
        assert!(
            (area - want).abs() < 0.01 * want,
            "band ({}, {}) has area {area}, expected {want}",
            band.low,
            band.high
        );
        assert!(!band.polygons.is_empty());
        for polygon in &band.polygons {
            assert!(polygon.len() >= 3, "a filled polygon needs three vertices");
            for &(x, y) in polygon {
                let value = x.hypot(y);
                assert!(
                    (band.low - 0.01..=band.high + 0.01).contains(&value),
                    "vertex ({x}, {y}) with value {value} escapes band ({}, {})",
                    band.low,
                    band.high
                );
            }
        }
    }
    // The bands are an exact partition: their areas add up to the rectangle.
    let total: f64 = bands
        .iter()
        .flat_map(|band| band.polygons.iter())
        .map(|p| polygon_area(p))
        .sum();
    assert!(
        (total - 16.0).abs() < 1e-9,
        "bands overlap or leave gaps: {total}"
    );
}

#[test]
fn contour_bands_outside_the_data_range_fill_the_rectangle() {
    let bands =
        plot::contour_bands_2d(|x, y| x * x + y * y, (-1.0, 1.0), (-1.0, 1.0), 21, &[100.0])
            .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(bands.len(), 2);
    assert!(bands[1].polygons.is_empty());
    let area: f64 = bands[0].polygons.iter().map(|p| polygon_area(p)).sum();
    assert!(
        (area - 4.0).abs() < 0.01,
        "below-level band covers {area} of 4"
    );
}

#[test]
fn contour_bands_svg_has_one_polygon_per_filled_polygon() {
    let levels = [0.5, 1.5];
    let bands = plot::contour_bands_2d(|x, y| x * x + y * y, (-2.0, 2.0), (-2.0, 2.0), 31, &levels)
        .unwrap_or_else(|e| panic!("{e}"));
    let expected: usize = bands.iter().map(|band| band.polygons.len()).sum();
    assert!(expected > 0);
    let f = Tmp::new("bands", "svg");
    plot::plot_contour_bands_2d(
        |x, y| x * x + y * y,
        (-2.0, 2.0),
        (-2.0, 2.0),
        31,
        &levels,
        &f.0,
    )
    .unwrap_or_else(|e| panic!("{e}"));
    let svg = read_svg(&f);
    assert!(svg.starts_with("<svg"));
    assert!(svg.trim_end().ends_with("</svg>"));
    assert_eq!(svg.matches("<polygon").count(), expected);
    // Deterministic: the same request produces the same bands.
    let again = plot::contour_bands_2d(|x, y| x * x + y * y, (-2.0, 2.0), (-2.0, 2.0), 31, &levels)
        .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(bands, again);
}

#[test]
fn contour_bands_reject_invalid_requests() {
    let f = Tmp::new("bands_bad", "svg");
    assert!(plot::contour_bands_2d(|x, y| x + y, (0.0, 1.0), (0.0, 1.0), 11, &[]).is_err());
    assert!(plot::contour_bands_2d(|x, y| x + y, (0.0, 1.0), (0.0, 1.0), 11, &[1.0, 1.0]).is_err());
    assert!(plot::contour_bands_2d(|x, y| x + y, (0.0, 1.0), (0.0, 1.0), 11, &[1.0, 0.0]).is_err());
    assert!(plot::contour_bands_2d(|x, y| x + y, (0.0, 1.0), (0.0, 1.0), 11, &[f64::NAN]).is_err());
    assert!(plot::contour_bands_2d(|x, y| x + y, (0.0, 1.0), (0.0, 1.0), 1, &[0.0]).is_err());
    assert!(plot::contour_bands_2d(|x, y| x + y, (1.0, 1.0), (0.0, 1.0), 11, &[0.0]).is_err());
    assert!(
        plot::plot_contour_bands_2d(|x, y| x + y, (0.0, 1.0), (0.0, 1.0), 11, &[], &f.0).is_err()
    );
    assert!(!f.0.exists());
}

#[test]
fn pgm_raster_has_the_expected_header_and_exact_samples() {
    let data = arr2(&[[0.0, 1.0], [2.0, 3.0]]);
    let f = Tmp::new("field_pgm", "pgm");
    plot::plot_scalar_field_pgm(&data, &f.0).unwrap_or_else(|e| panic!("{e}"));
    let bytes = read_bytes(&f);
    assert_eq!(&bytes[..11], b"P5\n2 2\n255\n");
    assert_eq!(&bytes[11..], &[0, 85, 170, 255]);
}

#[test]
fn ppm_raster_uses_the_heat_map_palette() {
    let data = arr2(&[[0.0, 1.0], [2.0, 3.0]]);
    let f = Tmp::new("field_ppm", "ppm");
    plot::plot_scalar_field_ppm(&data, &f.0).unwrap_or_else(|e| panic!("{e}"));
    let bytes = read_bytes(&f);
    assert_eq!(&bytes[..11], b"P6\n2 2\n255\n");
    assert_eq!(bytes.len(), 11 + 2 * 2 * 3);
    assert_eq!(
        &bytes[11..14],
        &[68, 1, 84],
        "t = 0 is the low palette stop"
    );
    assert_eq!(
        &bytes[20..23],
        &[253, 231, 37],
        "t = 1 is the high palette stop"
    );
}

#[test]
fn scalar_field_rasters_handle_non_finite_entries() {
    let data = arr2(&[[f64::NAN, 5.0], [5.0, 5.0]]);
    let g = Tmp::new("field_nan", "pgm");
    plot::plot_scalar_field_pgm(&data, &g.0).unwrap_or_else(|e| panic!("{e}"));
    let pgm = read_bytes(&g);
    // A constant finite range becomes mid-grey; NaN becomes black.
    assert_eq!(&pgm[11..], &[0, 128, 128, 128]);
    let p = Tmp::new("field_nan_rgb", "ppm");
    plot::plot_scalar_field_ppm(&data, &p.0).unwrap_or_else(|e| panic!("{e}"));
    let ppm = read_bytes(&p);
    assert_eq!(&ppm[11..14], &[0xdd, 0xdd, 0xdd]);
    assert_eq!(
        &ppm[14..17],
        &[33, 145, 140],
        "mid palette stop for the constant value"
    );
}

#[test]
fn raster_writers_encode_headers_and_validate_shapes() {
    let pgm = Tmp::new("pgm_raw", "pgm");
    raster::write_pgm_file(&pgm.0, 1, 2, &[10, 20]).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(read_bytes(&pgm), b"P5\n1 2\n255\n\x0a\x14");
    let pgm16 = Tmp::new("pgm16_raw", "pgm");
    raster::write_pgm16_file(&pgm16.0, 1, 1, &[258]).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(read_bytes(&pgm16), b"P5\n1 1\n65535\n\x01\x02");
    let ppm = Tmp::new("ppm_raw", "ppm");
    raster::write_ppm_file(&ppm.0, 1, 1, &[1, 2, 3]).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(read_bytes(&ppm), b"P6\n1 1\n255\n\x01\x02\x03");

    let bad = Tmp::new("pgm_bad", "pgm");
    assert!(raster::write_pgm_file(&bad.0, 2, 2, &[1, 2, 3]).is_err());
    assert!(raster::write_pgm_file(&bad.0, 0, 1, &[]).is_err());
    assert!(raster::write_ppm_file(&bad.0, 1, 1, &[1, 2]).is_err());
    assert!(raster::write_pgm16_file(&bad.0, 2, 1, &[1]).is_err());
    assert!(!bad.0.exists(), "a rejected raster must not create a file");
    let missing = std::env::temp_dir()
        .join("rssn_io_no_such_dir")
        .join("raster.pgm");
    assert!(raster::write_pgm_file(&missing, 1, 1, &[0]).is_err());
    let empty: Array2<f64> = Array2::zeros((0, 0));
    assert!(plot::plot_scalar_field_pgm(&empty, &bad.0).is_err());
    let no_finite = Array2::from_elem((2, 2), f64::NAN);
    assert!(plot::plot_scalar_field_ppm(&no_finite, &bad.0).is_err());
}

// ---- Quiver plot (IO-2) ----

#[test]
fn quiver_uniform_field_draws_horizontal_arrows() {
    let f = Tmp::new("quiver", "svg");
    plot::plot_quiver_2d(|_, _| (1.0, 0.0), (-1.0, 1.0), (-1.0, 1.0), 5, 1.0, &f.0)
        .unwrap_or_else(|e| panic!("{e}"));
    let svg = read_svg(&f);
    assert!(svg.starts_with("<svg"));
    assert!(svg.trim_end().ends_with("</svg>"));
    assert_eq!(svg.matches("<polygon").count(), 25, "one head per sample");
    let shafts: Vec<&str> = svg
        .lines()
        .filter(|line| line.starts_with("<line") && line.contains("stroke-width=\"1.2\""))
        .collect();
    assert_eq!(shafts.len(), 25, "one shaft per sample");
    for shaft in shafts {
        let y1: f64 = svg_attr(shaft, "y1")
            .parse()
            .unwrap_or_else(|e| panic!("{e}"));
        let y2: f64 = svg_attr(shaft, "y2")
            .parse()
            .unwrap_or_else(|e| panic!("{e}"));
        let x1: f64 = svg_attr(shaft, "x1")
            .parse()
            .unwrap_or_else(|e| panic!("{e}"));
        let x2: f64 = svg_attr(shaft, "x2")
            .parse()
            .unwrap_or_else(|e| panic!("{e}"));
        assert!(
            (y1 - y2).abs() < 1e-9,
            "a uniform field cannot tilt an arrow"
        );
        assert!(x2 > x1, "the field points along +x");
    }
}

#[test]
fn quiver_scales_lengths_and_drops_zero_vectors() {
    let f = Tmp::new("quiver_scaled", "svg");
    plot::plot_quiver_2d(|x, _| (x, 0.0), (-1.0, 1.0), (-1.0, 1.0), 3, 1.0, &f.0)
        .unwrap_or_else(|e| panic!("{e}"));
    let svg = read_svg(&f);
    assert_eq!(
        svg.matches("<polygon").count(),
        6,
        "the x = 0 column has no arrows"
    );
    let mut lengths: Vec<f64> = Vec::new();
    for shaft in svg
        .lines()
        .filter(|line| line.starts_with("<line") && line.contains("stroke-width=\"1.2\""))
    {
        let x1: f64 = svg_attr(shaft, "x1")
            .parse()
            .unwrap_or_else(|e| panic!("{e}"));
        let x2: f64 = svg_attr(shaft, "x2")
            .parse()
            .unwrap_or_else(|e| panic!("{e}"));
        lengths.push((x2 - x1).abs());
    }
    assert_eq!(lengths.len(), 6);
    for length in &lengths {
        assert!(
            (length - lengths[0]).abs() < 1e-9,
            "equal magnitudes must get equal lengths"
        );
    }
    let zero = Tmp::new("quiver_zero", "svg");
    plot::plot_quiver_2d(|_, _| (0.0, 0.0), (-1.0, 1.0), (-1.0, 1.0), 5, 1.0, &zero.0)
        .unwrap_or_else(|e| panic!("{e}"));
    let zero_svg = read_svg(&zero);
    assert_eq!(zero_svg.matches("<polygon").count(), 0);
    assert!(!zero_svg.contains("stroke-width=\"1.2\""));
}

#[test]
fn quiver_is_deterministic_and_rejects_invalid_requests() {
    let a = Tmp::new("quiver_det_a", "svg");
    let b = Tmp::new("quiver_det_b", "svg");
    plot::plot_quiver_2d(|x, y| (-y, x), (-1.0, 1.0), (-1.0, 1.0), 5, 0.8, &a.0)
        .unwrap_or_else(|e| panic!("{e}"));
    plot::plot_quiver_2d(|x, y| (-y, x), (-1.0, 1.0), (-1.0, 1.0), 5, 0.8, &b.0)
        .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(read_svg(&a), read_svg(&b));

    let f = Tmp::new("quiver_bad", "svg");
    assert!(
        plot::plot_quiver_2d(|_, _| (1.0, 0.0), (-1.0, 1.0), (-1.0, 1.0), 5, 0.0, &f.0).is_err()
    );
    assert!(
        plot::plot_quiver_2d(
            |_, _| (1.0, 0.0),
            (-1.0, 1.0),
            (-1.0, 1.0),
            5,
            f64::NAN,
            &f.0
        )
        .is_err()
    );
    assert!(
        plot::plot_quiver_2d(|_, _| (1.0, 0.0), (-1.0, 1.0), (-1.0, 1.0), 1, 1.0, &f.0).is_err()
    );
    assert!(
        plot::plot_quiver_2d(|_, _| (1.0, 0.0), (1.0, 1.0), (-1.0, 1.0), 5, 1.0, &f.0).is_err()
    );
    assert!(!f.0.exists());
}
