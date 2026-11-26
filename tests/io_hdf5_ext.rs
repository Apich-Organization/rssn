//! Extended HDF5 integration tests: chunked write/read, all filters,
//! dense groups, vlen numeric datasets, soft links, and v4 chunked layouts.

use std::path::PathBuf;

use rssn::io::hdf5::*;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/hdf5")
        .join(name)
}

#[test]
fn reads_chunked_shuffle_deflate_fletcher32_fixture() {
    let f = read_hdf5(&fixture("chunked_filters.h5")).expect("read chunked_filters.h5");
    let ds = f
        .datasets
        .iter()
        .find(|d| d.path == "/deflate_shuffle_f32")
        .expect("dataset present");
    assert_eq!(ds.shape, vec![10, 12]);
    assert_eq!(ds.data.len(), 120);
    for (i, &v) in ds.data.iter().enumerate() {
        assert_eq!(v, i as f64);
    }
}

#[test]
fn reads_v4_chunked_index_types_fixture() {
    let f = read_hdf5(&fixture("v4_chunked.h5")).expect("read v4_chunked.h5");
    for name in ["/single_chunk", "/fixed_array", "/ext_array", "/btree_v2"] {
        let ds = f
            .datasets
            .iter()
            .find(|d| d.path == name)
            .unwrap_or_else(|| panic!("missing {name}"));
        assert_eq!(ds.shape, vec![4, 6]);
        for (i, &v) in ds.data.iter().enumerate() {
            assert_eq!(v, i as f64, "mismatch in {name} at {i}");
        }
    }
}

#[test]
fn reads_vlen_numeric_fixture() {
    let f = read_hdf5(&fixture("vlen_numeric.h5")).expect("read vlen_numeric.h5");
    let ds_i = f
        .datasets
        .iter()
        .find(|d| d.path == "/vlen_i32")
        .expect("vlen_i32");
    assert_eq!(ds_i.h5_type, H5Type::VlenNumeric);
    assert_eq!(
        ds_i.vlen_data,
        vec![vec![1.0, 2.0, 3.0], vec![], vec![10.0, 20.0]]
    );

    let ds_f = f
        .datasets
        .iter()
        .find(|d| d.path == "/vlen_f64")
        .expect("vlen_f64");
    assert_eq!(ds_f.vlen_data, vec![vec![1.5, 2.5], vec![3.25]]);
}

#[test]
fn reads_soft_links_fixture() {
    let f = read_hdf5(&fixture("soft_links.h5")).expect("read soft_links.h5");
    assert_eq!(f.resolve_soft_link("/alias"), Some("/grp/target"));
    assert_eq!(f.resolve_soft_link("/grp/rel_alias"), Some("/grp/target"));
    let alias_ds = f
        .datasets
        .iter()
        .find(|d| d.path == "/alias")
        .expect("alias dataset resolved");
    assert_eq!(alias_ds.data, vec![10.0, 20.0, 30.0]);
}

#[test]
fn reads_plugin_filters_fixture() {
    let f = read_hdf5(&fixture("filters_plugins.h5")).expect("read filters_plugins.h5");
    for name in ["/lz4", "/zstd", "/blosc_lz4", "/blosc_zstd"] {
        let ds = f
            .datasets
            .iter()
            .find(|d| d.path == name)
            .unwrap_or_else(|| panic!("missing {name}"));
        assert_eq!(ds.data.len(), 100);
        for (i, &v) in ds.data.iter().enumerate() {
            assert_eq!(v, i as f64, "mismatch in {name} at {i}");
        }
    }
}

#[test]
fn round_trip_chunked_filters_dense_vlen_and_soft_links() {
    let grid = Dataset::new_numeric(
        "/grid",
        vec![6, 8],
        H5Type::F64,
        (0..48).map(|i| i as f64).collect(),
    );
    let opts_chunked = H5WriteOptions {
        dense_groups: false,
        dataset_options: vec![(
            "/grid".to_owned(),
            DatasetWriteOptions {
                chunk_dims: Some(vec![3, 4]),
                filters: vec![
                    H5Filter::Shuffle,
                    H5Filter::Deflate(6),
                    H5Filter::Fletcher32,
                ],
            },
        )],
        soft_links: Vec::new(),
    };
    let bytes = write_hdf5_bytes_with_options(std::slice::from_ref(&grid), &[], &opts_chunked)
        .expect("write chunked");
    let back = read_hdf5_bytes(&bytes).expect("read chunked");
    assert_eq!(back.datasets[0].data, grid.data);

    // Dense group write + read
    let mut dense_ds = Vec::new();
    for i in 0..20 {
        dense_ds.push(Dataset::new_numeric(
            format!("/grp/d{i}"),
            vec![2],
            H5Type::F64,
            vec![i as f64, (i * 2) as f64],
        ));
    }
    let opts_dense = H5WriteOptions {
        dense_groups: true,
        ..Default::default()
    };
    let dbytes = write_hdf5_bytes_with_options(&dense_ds, &[], &opts_dense).expect("write dense");
    let dback = read_hdf5_bytes(&dbytes).expect("read dense");
    assert_eq!(dback.datasets.len(), 20);

    // Vlen numeric write + read
    let vlen_ds = Dataset::new_vlen_numeric(
        "/vlen_i32",
        vec![3],
        H5Type::I32,
        vec![vec![1.0, 2.0, 3.0], vec![], vec![10.0, 20.0]],
    );
    let vbytes = write_hdf5_bytes_with_options(&[vlen_ds], &[], &H5WriteOptions::default())
        .expect("write vlen");
    let vback = read_hdf5_bytes(&vbytes).expect("read vlen");
    assert_eq!(
        vback.datasets[0].vlen_data,
        vec![vec![1.0, 2.0, 3.0], vec![], vec![10.0, 20.0]]
    );

    // Soft link write + read
    let base_ds = Dataset::new_numeric("/data", vec![3], H5Type::F64, vec![1.0, 2.0, 3.0]);
    let sbytes = write_hdf5_bytes_with_options(
        &[base_ds],
        &[],
        &H5WriteOptions {
            soft_links: vec![SoftLink {
                path: "/alias".to_owned(),
                target: "/data".to_owned(),
            }],
            ..Default::default()
        },
    )
    .expect("write soft link");
    let sback = read_hdf5_bytes(&sbytes).expect("read soft link");
    assert_eq!(sback.resolve_soft_link("/alias"), Some("/data"));
    assert!(sback.datasets.iter().any(|d| d.path == "/alias"));

    // Optionally dump for Python cross-verification
    if let Ok(dir) = std::env::var("RSSN_VERIFY_DIR") {
        let dir = PathBuf::from(dir);
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("rssn_chunked.h5"), &bytes).unwrap();
        std::fs::write(dir.join("rssn_dense.h5"), &dbytes).unwrap();
        std::fs::write(dir.join("rssn_vlen.h5"), &vbytes).unwrap();
        std::fs::write(dir.join("rssn_soft.h5"), &sbytes).unwrap();

        let pdata = Dataset::new_numeric(
            "/lz4",
            vec![64],
            H5Type::F64,
            (0..64).map(|i| i as f64).collect(),
        );
        let p2 = Dataset::new_numeric("/zstd", vec![64], H5Type::F64, pdata.data.clone());
        let p3 = Dataset::new_numeric("/blosc_lz4", vec![64], H5Type::F64, pdata.data.clone());
        let popts = H5WriteOptions {
            dataset_options: vec![
                (
                    "/lz4".to_owned(),
                    DatasetWriteOptions {
                        chunk_dims: Some(vec![32]),
                        filters: vec![H5Filter::Lz4],
                    },
                ),
                (
                    "/zstd".to_owned(),
                    DatasetWriteOptions {
                        chunk_dims: Some(vec![32]),
                        filters: vec![H5Filter::Zstd(3)],
                    },
                ),
                (
                    "/blosc_lz4".to_owned(),
                    DatasetWriteOptions {
                        chunk_dims: Some(vec![32]),
                        filters: vec![H5Filter::Blosc {
                            clevel: 5,
                            shuffle: 1,
                            compressor: 1,
                        }],
                    },
                ),
            ],
            ..Default::default()
        };
        write_hdf5_with_options(&dir.join("rssn_plugins.h5"), &[pdata, p2, p3], &[], &popts)
            .expect("write plugins");
    }
}
