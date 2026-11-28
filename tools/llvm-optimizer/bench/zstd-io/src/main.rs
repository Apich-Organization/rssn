#![allow(dead_code, unused_imports)]
mod io;
use io::frame::Column;
use io::hdf5;
use io::parquet;
use std::collections::hash_map::DefaultHasher;
use std::hash::Hash;
use std::hash::Hasher;
use std::hint::black_box;
use std::path::Path;
use std::time::Duration;

fn time() -> Duration {
    let mut ts = std::mem::MaybeUninit::<libc::timespec>::uninit();
    assert_eq!(
        unsafe { libc::clock_gettime(libc::CLOCK_THREAD_CPUTIME_ID, ts.as_mut_ptr()) },
        0
    );
    let ts = unsafe { ts.assume_init() };
    Duration::new(
        ts.tv_sec.try_into().unwrap(),
        ts.tv_nsec.try_into().unwrap(),
    )
}

fn generated_parquet() -> Vec<u8> {
    let rows = 48_000usize;
    let cols = vec![
        Column::Float(
            (0..rows)
                .map(|i| ((i % 141) as f64).sin() * 0.25 + (i % 12) as f64)
                .collect(),
        ),
        Column::Str((0..rows).map(|i| format!("zone-{}", i % 13)).collect()),
    ];
    let schema = parquet::Schema {
        fields: vec![
            parquet::Field {
                name: "signal".into(),
                physical: parquet::Physical::Double,
                optional: false,
            },
            parquet::Field {
                name: "station".into(),
                physical: parquet::Physical::ByteArray,
                optional: false,
            },
        ],
    };
    let options = parquet::WriteOptions {
        compression: parquet::Compression::Zstd,
        max_page_rows: 96,
        max_row_group_rows: 16_000,
        data_page_version: 2,
        ..Default::default()
    };
    let bytes = parquet::write_parquet_bytes_with_options(&schema, &cols, &options).unwrap();
    assert_eq!(parquet::read_parquet_bytes(&bytes).unwrap().1, cols);
    bytes
}

fn generated_hdf5() -> Vec<u8> {
    let rows = 12_000usize;
    let dataset = hdf5::Dataset::new_numeric(
        "/signal",
        vec![rows, 8],
        hdf5::H5Type::F64,
        (0..rows * 8)
            .map(|i| ((i % 111) as f64).sin() + (i % 15) as f64)
            .collect(),
    );
    let options = hdf5::H5WriteOptions {
        dataset_options: vec![(
            "/signal".into(),
            hdf5::DatasetWriteOptions {
                chunk_dims: Some(vec![48, 8]),
                filters: vec![hdf5::H5Filter::Shuffle, hdf5::H5Filter::Zstd(3)],
            },
        )],
        ..Default::default()
    };
    let bytes =
        hdf5::write_hdf5_bytes_with_options(std::slice::from_ref(&dataset), &[], &options).unwrap();
    assert_eq!(hdf5::read_hdf5_bytes(&bytes).unwrap().datasets, [dataset]);
    bytes
}

fn hash_column(
    c: &Column,
    s: &mut DefaultHasher,
) {
    std::mem::discriminant(c).hash(s);
    match c {
        | Column::Int(v) => v.hash(s),
        | Column::Float(v) => v.iter().map(|x| x.to_bits()).for_each(|x| x.hash(s)),
        | Column::Bool(v) => v.hash(s),
        | Column::Str(v) => v.hash(s),
        | Column::Bytes(v) => v.hash(s),
        | Column::Nullable(v, n) => {
            hash_column(v, s);
            n.hash(s);
        },
        | Column::List { values, offsets } => {
            hash_column(values, s);
            offsets.hash(s);
        },
    }
}

fn parquet_rounds(
    bytes: &[u8],
    rounds: usize,
) {
    let mut digest = DefaultHasher::new();
    for c in &parquet::read_parquet_bytes(bytes).unwrap().1 {
        hash_column(c, &mut digest);
    }
    let t = time();
    for _ in 0..rounds {
        black_box(parquet::read_parquet_bytes(black_box(bytes)).unwrap());
    }
    println!("{} {:016x}", (time() - t).as_nanos(), digest.finish());
}

fn hdf5_rounds(
    bytes: &[u8],
    rounds: usize,
) {
    let mut digest = DefaultHasher::new();
    for d in &hdf5::read_hdf5_bytes(bytes).unwrap().datasets {
        d.path.hash(&mut digest);
        d.shape.hash(&mut digest);
        d.data
            .iter()
            .map(|x| x.to_bits())
            .for_each(|x| x.hash(&mut digest));
    }
    let t = time();
    for _ in 0..rounds {
        black_box(hdf5::read_hdf5_bytes(black_box(bytes)).unwrap());
    }
    println!("{} {:016x}", (time() - t).as_nanos(), digest.finish());
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    if args[1] == "generate" {
        std::fs::write(&args[2], generated_parquet()).unwrap();
        return;
    }
    if args[1] == "generate-hdf5" {
        std::fs::write(&args[2], generated_hdf5()).unwrap();
        return;
    }
    let bytes = std::fs::read(Path::new(&args[2])).unwrap();
    let rounds = args[3].parse().unwrap();
    match args[1].as_str() {
        | "parquet" => parquet_rounds(&bytes, rounds),
        | "hdf5" => hdf5_rounds(&bytes, rounds),
        | _ => panic!("expected parquet or hdf5"),
    }
}
