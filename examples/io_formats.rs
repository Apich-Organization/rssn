//! Reading and writing the tabular and scientific file formats.
//!
//! Run with `cargo run --locked --example io_formats`.
//!
//! The example writes a small file in each format, reads it back and prints
//! what came out. The files land in the system temporary directory, so they can
//! also be opened with `h5py`, `pyarrow` or the NetCDF tools to check that
//! other implementations agree.

use std::path::Path;

use rssn::io::Colormap;
use rssn::io::GifAnimation;
use rssn::io::GifDisposal;
use rssn::io::GifFrame;
use rssn::io::PngColorType;
use rssn::io::PngImage;
use rssn::io::WavData;
use rssn::io::WavSpec;
use rssn::io::decode_gif_composited;
use rssn::io::frame::Column;
use rssn::io::frame::DType;
use rssn::io::frame::DataFrame;
use rssn::io::frame::write_csv;
use rssn::io::hdf5::AttrValue;
use rssn::io::hdf5::Attribute;
use rssn::io::hdf5::Dataset;
use rssn::io::hdf5::DatasetWriteOptions;
use rssn::io::hdf5::H5Filter;
use rssn::io::hdf5::H5Type;
use rssn::io::hdf5::H5WriteOptions;
use rssn::io::hdf5::SoftLink;
use rssn::io::hdf5::read_hdf5;
use rssn::io::hdf5::write_hdf5;
use rssn::io::hdf5::write_hdf5_with_options;
use rssn::io::parquet::Compression;
use rssn::io::parquet::Field;
use rssn::io::parquet::Physical;
use rssn::io::parquet::Schema;
use rssn::io::parquet::WriteOptions;
use rssn::io::parquet::read_parquet;
use rssn::io::parquet::read_parquet_bloom_filters;
use rssn::io::parquet::read_parquet_page_index;
use rssn::io::parquet::write_parquet;
use rssn::io::parquet::write_parquet_with_options;
use rssn::io::read_png_file;
use rssn::io::read_wav_file;
use rssn::io::vtk::VtkAttributeData;
use rssn::io::vtk::VtkDataset;
use rssn::io::vtk::VtkFile;
use rssn::io::vtk::VtkFormat;
use rssn::io::vtk::read_vtk;
use rssn::io::vtk::write_vtk;
use rssn::io::write_gif_file;
use rssn::io::write_png_file;
use rssn::io::write_wav_file;

fn main() -> Result<(), String> {
    let directory = std::env::temp_dir().join("rssn-io-formats");
    std::fs::create_dir_all(&directory)
        .map_err(|error| format!("cannot create {}: {error}", directory.display()))?;

    hdf5(&directory)?;
    parquet(&directory)?;
    frame(&directory)?;
    media_and_vtk(&directory)?;
    Ok(())
}

fn hdf5(directory: &Path) -> Result<(), String> {
    let path = directory.join("rssn_example.h5");
    let datasets = vec![
        Dataset::new_numeric(
            "/observations/temperature",
            vec![3],
            H5Type::F64,
            vec![250.0, 251.5, 253.0],
        ),
        Dataset::new_strings(
            "/observations/station",
            vec![3],
            vec!["alpha".to_owned(), "beta".to_owned(), "gamma".to_owned()],
        ),
        Dataset::new_numeric("/counts", vec![2], H5Type::I32, vec![7.0, 42.0]),
    ];
    let attributes = vec![(
        "/observations/temperature".to_owned(),
        vec![Attribute {
            name: "units".to_owned(),
            value: AttrValue::Text("K".to_owned()),
        }],
    )];
    write_hdf5(&path, &datasets, &attributes)?;
    let file = read_hdf5(&path)?;
    println!(
        "hdf5: {} groups, {} datasets",
        file.groups.len(),
        file.datasets.len()
    );
    for dataset in &file.datasets {
        if dataset.strings.is_empty() {
            println!("  {} {:?} {:?}", dataset.path, dataset.shape, dataset.data);
        } else {
            println!(
                "  {} {:?} {:?}",
                dataset.path, dataset.shape, dataset.strings
            );
        }
    }
    for (path, attributes) in &file.attributes {
        for attribute in attributes {
            println!(
                "  attribute {path}: {} = {:?}",
                attribute.name, attribute.value
            );
        }
    }
    Ok(())
}

fn parquet(directory: &Path) -> Result<(), String> {
    let path = directory.join("rssn_example.parquet");
    let schema = Schema {
        fields: vec![
            Field {
                name: "label".to_owned(),
                physical: Physical::ByteArray,
                optional: false,
            },
            Field {
                name: "value".to_owned(),
                physical: Physical::Double,
                optional: true,
            },
            Field {
                name: "xs".to_owned(),
                physical: Physical::Int32,
                optional: true,
            },
        ],
    };
    // A nested list column: [[1, 2], [3]], null, [].
    let inner = Column::new_list(
        Column::Nullable(
            Box::new(Column::Int(vec![1, 0, 3])),
            vec![false, true, false],
        ),
        vec![0, 2, 3, 3],
    )
    .map_err(|error| error.to_string())?;
    let nested = Column::new_list(
        Column::Nullable(Box::new(inner), vec![false, false, false]),
        vec![0, 2, 2, 3],
    )
    .map_err(|error| error.to_string())?;
    let columns = vec![
        Column::Str(vec!["a".to_owned(), "b".to_owned(), "c".to_owned()]),
        Column::Nullable(
            Box::new(Column::Float(vec![1.0, 0.0, 3.0])),
            vec![false, true, false],
        ),
        Column::Nullable(Box::new(nested), vec![false, true, false]),
    ];
    write_parquet(&path, &schema, &columns)?;
    // The same table with version-2 data pages.
    let v2 = directory.join("rssn_example_v2.parquet");
    let options = rssn::io::parquet::WriteOptions {
        data_page_version: 2,
        ..Default::default()
    };
    rssn::io::parquet::write_parquet_with_options(&v2, &schema, &columns, &options)?;
    let (read_schema, read_columns) = read_parquet(&path)?;
    println!("parquet: {} fields", read_schema.fields.len());
    for (field, column) in read_schema.fields.iter().zip(&read_columns) {
        println!(
            "  {} {:?} {} rows",
            field.name,
            column.dtype(),
            column.len()
        );
    }
    Ok(())
}

fn frame(directory: &Path) -> Result<(), String> {
    let frame = DataFrame::new(vec![
        (
            "label".to_owned(),
            Column::Str(vec!["a".to_owned(), "b".to_owned()]),
        ),
        (
            "xs".to_owned(),
            Column::new_list(Column::Int(vec![1, 2, 3]), vec![0, 2, 3])
                .map_err(|error| error.to_string())?,
        ),
    ])
    .map_err(|error| error.to_string())?;
    let text = write_csv(&frame, &Default::default());
    let path = directory.join("rssn_example.csv");
    std::fs::write(&path, &text).map_err(|error| error.to_string())?;
    println!(
        "frame: {} columns, csv at {}",
        frame.n_cols(),
        path.display()
    );
    assert_eq!(frame.column("xs").map(Column::dtype), Some(DType::List));
    Ok(())
}

fn media_and_vtk(directory: &Path) -> Result<(), String> {
    // 1. Chunked + Shuffle + Deflate + Fletcher32 HDF5 + SoftLink + VlenNumeric
    let h5_ext = directory.join("rssn_extended.h5");
    let ds_grid = Dataset::new_numeric(
        "/grid",
        vec![4, 6],
        H5Type::F64,
        (0..24).map(|i| i as f64).collect(),
    );
    let ds_vlen = Dataset::new_vlen_numeric(
        "/ragged",
        vec![2],
        H5Type::I32,
        vec![vec![1.0, 2.0, 3.0], vec![10.0, 20.0]],
    );
    write_hdf5_with_options(
        &h5_ext,
        &[ds_grid, ds_vlen],
        &[],
        &H5WriteOptions {
            dense_groups: true,
            dataset_options: vec![(
                "/grid".to_owned(),
                DatasetWriteOptions {
                    chunk_dims: Some(vec![2, 3]),
                    filters: vec![
                        H5Filter::Shuffle,
                        H5Filter::Deflate(6),
                        H5Filter::Fletcher32,
                    ],
                },
            )],
            soft_links: vec![SoftLink {
                path: "/alias".to_owned(),
                target: "/grid".to_owned(),
            }],
        },
    )?;
    let h5_back = read_hdf5(&h5_ext)?;
    println!(
        "hdf5 extended: {} datasets, soft_links={:?}",
        h5_back.datasets.len(),
        h5_back.soft_links
    );

    // 2. Parquet with Brotli, BloomFilter, and PageIndex
    let pq_ext = directory.join("rssn_extended.parquet");
    let schema = Schema {
        fields: vec![Field {
            name: "id".to_owned(),
            physical: Physical::Int64,
            optional: false,
        }],
    };
    let cols = vec![Column::Int((0..40).collect())];
    write_parquet_with_options(
        &pq_ext,
        &schema,
        &cols,
        &WriteOptions {
            compression: Compression::Brotli,
            max_page_rows: 10,
            write_bloom_filter: true,
            write_page_index: true,
            ..Default::default()
        },
    )?;
    let pq_bytes = std::fs::read(&pq_ext).map_err(|e| e.to_string())?;
    let bfs = read_parquet_bloom_filters(&pq_bytes)?;
    let pidx = read_parquet_page_index(&pq_bytes)?;
    println!(
        "parquet extended: bloom_present={}, pages={}",
        bfs[0][0].as_ref().is_some_and(|b| b.might_contain_i64(17)),
        pidx[0][0]
            .offset_index
            .as_ref()
            .map_or(0, |o| o.page_locations.len())
    );

    // 3. VTK Binary
    let vtk_path = directory.join("rssn_example.vtk");
    write_vtk(
        &vtk_path,
        &VtkFile {
            title: "RSSN Example".to_owned(),
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
    )?;
    let vtk_back = read_vtk(&vtk_path)?;
    println!(
        "vtk: {:?}, points={}",
        vtk_back.format,
        vtk_back.dataset.num_points()
    );

    // 4. PNG, GIF, WAV
    let png_path = directory.join("rssn_example.png");
    let mut rgb = Vec::with_capacity(4 * 4 * 3);
    for i in 0..16 {
        let (r, g, b) = Colormap::Viridis.sample_rgb((i as f64) / 15.0);
        rgb.extend_from_slice(&[r, g, b]);
    }
    write_png_file(
        &png_path,
        &PngImage {
            width: 4,
            height: 4,
            color_type: PngColorType::Rgb8,
            data: rgb,
        },
    )?;
    let png_back = read_png_file(&png_path)?;
    println!("png: {}x{}", png_back.width, png_back.height);

    let gif_path = directory.join("rssn_example.gif");
    let mut f1 = GifFrame::new(4, 4, 10, [255, 0, 0, 255].repeat(16));
    f1.dispose = GifDisposal::Background;
    let mut f2 = GifFrame::new(2, 2, 15, [0, 255, 0, 255].repeat(4));
    f2.left = 1;
    f2.top = 1;
    write_gif_file(&gif_path, &GifAnimation::new(4, 4, 0, vec![f1, f2]))?;
    let gif_bytes = std::fs::read(&gif_path).map_err(|e| e.to_string())?;
    let comp = decode_gif_composited(&gif_bytes)?;
    println!("gif: {} composited frames", comp.frames.len());

    let wav_path = directory.join("rssn_example.wav");
    write_wav_file(
        &wav_path,
        &WavData::new(
            WavSpec {
                channels: 1,
                sample_rate: 16_000,
                bits_per_sample: 24,
                is_float: false,
            },
            vec![0.0, 0.5, -0.5, 0.25],
        ),
    )?;
    let wav_back = read_wav_file(&wav_path)?;
    println!(
        "wav: {} Hz, {} samples",
        wav_back.spec.sample_rate,
        wav_back.samples.len()
    );
    Ok(())
}
