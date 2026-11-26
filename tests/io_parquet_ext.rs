//! Extended Parquet integration tests: Brotli, LZO, Bloom Filter, Page Index.

use std::path::PathBuf;

use rssn::io::frame::Column;
use rssn::io::parquet::*;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/parquet")
        .join(name)
}

#[test]
fn reads_arrow_brotli_and_page_index_fixtures() {
    let (_, cols_b) = read_parquet(&fixture("codec_brotli.parquet")).expect("read brotli");
    assert_eq!(cols_b[0].len(), 2000);

    let pbytes = std::fs::read(fixture("page_index.parquet")).expect("read page_index bytes");
    let (_, cols_p) = read_parquet_bytes(&pbytes).expect("read page_index");
    assert_eq!(cols_p[0].len(), 2000);
    let idx = read_parquet_page_index(&pbytes).expect("read page index");
    assert!(!idx.is_empty());
    assert!(idx[0][0].column_index.is_some());
    assert!(idx[0][0].offset_index.is_some());
}

#[test]
fn round_trip_brotli_lzo_bloom_and_page_index() {
    let schema = Schema {
        fields: vec![
            Field {
                name: "id".to_owned(),
                physical: Physical::Int64,
                optional: false,
            },
            Field {
                name: "tag".to_owned(),
                physical: Physical::ByteArray,
                optional: true,
            },
        ],
    };
    let ids: Vec<i64> = (0..100).collect();
    let tags: Vec<String> = (0..100).map(|i| format!("t-{}", i % 10)).collect();
    let cols = vec![
        Column::Int(ids),
        Column::new_nullable(Column::Str(tags), vec![false; 100]).unwrap(),
    ];

    for comp in [
        Compression::Brotli,
        Compression::Lzo,
        Compression::Zstd,
        Compression::Lz4,
        Compression::Snappy,
        Compression::Gzip,
    ] {
        let opts = WriteOptions {
            compression: comp,
            write_bloom_filter: true,
            write_page_index: true,
            max_page_rows: 25,
            ..Default::default()
        };
        let bytes = write_parquet_bytes_with_options(&schema, &cols, &opts)
            .unwrap_or_else(|e| panic!("write {comp:?}: {e}"));
        let (s_back, c_back) = read_parquet_bytes(&bytes).expect("read back");
        assert_eq!(s_back.fields.len(), 2);
        assert_eq!(c_back[0].len(), 100);

        let bfs = read_parquet_bloom_filters(&bytes).expect("bloom filters");
        let bf0 = bfs[0][0].as_ref().expect("id bloom filter");
        assert!(bf0.might_contain_i64(0));
        assert!(bf0.might_contain_i64(99));
        let bf1 = bfs[0][1].as_ref().expect("tag bloom filter");
        assert!(bf1.might_contain_str("t-3"));

        let pidx = read_parquet_page_index(&bytes).expect("page index");
        let off_idx = pidx[0][0].offset_index.as_ref().expect("offset index");
        assert_eq!(off_idx.page_locations.len(), 4);
        let col_idx = pidx[0][0].column_index.as_ref().expect("column index");
        assert_eq!(col_idx.min_values.len(), 4);
    }

    if let Ok(dir) = std::env::var("RSSN_VERIFY_DIR") {
        let dir = PathBuf::from(dir);
        std::fs::create_dir_all(&dir).unwrap();
        for (name, comp) in [
            ("brotli", Compression::Brotli),
            ("zstd", Compression::Zstd),
            ("lz4", Compression::Lz4),
            ("snappy", Compression::Snappy),
            ("gzip", Compression::Gzip),
        ] {
            let opts = WriteOptions {
                compression: comp,
                ..Default::default()
            };
            let bytes = write_parquet_bytes_with_options(&schema, &cols, &opts).unwrap();
            std::fs::write(dir.join(format!("rssn_{name}.parquet")), bytes).unwrap();
        }
        let small_cols = vec![
            Column::Int((0..20).collect()),
            Column::new_nullable(
                Column::Str((0..20).map(|i| format!("v{i}")).collect()),
                vec![false; 20],
            )
            .unwrap(),
        ];
        let opts_idx = WriteOptions {
            compression: Compression::Uncompressed,
            max_page_rows: 5,
            write_bloom_filter: true,
            write_page_index: true,
            ..Default::default()
        };
        let ibytes = write_parquet_bytes_with_options(&schema, &small_cols, &opts_idx).unwrap();
        std::fs::write(dir.join("rssn_indexed.parquet"), ibytes).unwrap();
    }
}
