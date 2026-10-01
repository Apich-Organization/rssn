//! Array I/O in `rssn::io` (CSV, JSON and, with the `npy` feature, NumPy files).

use assert_approx_eq::assert_approx_eq;
use ndarray::{Array2, arr2};
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::io::*;
use std::path::PathBuf;

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
        use std::sync::atomic::{AtomicUsize, Ordering};
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
