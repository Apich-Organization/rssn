# Test fixtures

[English](README.md) | [简体中文](README.zh-CN.md)

Small binary files that the `io::hdf5` and `io::parquet` tests read. They are
checked in so `cargo test` needs no extra tooling: the tests are ordinary
tests, not `#[ignore]`d ones.

Every file is produced by [`generate.py`](generate.py), which uses `h5py` and
`pyarrow` — the reference implementations of the two formats. Run it from the
repository root to rebuild them:

```sh
python tests/fixtures/generate.py
```

Keep the files small. The largest is `hdf5/dense_nested.h5` (about 1.2 MB),
which needs 40 000 links before a fractal heap grows a second level of
indirect blocks.

## hdf5/

| file | what it exercises |
| --- | --- |
| `dense_20.h5` | a dense group that fits in one fractal-heap block |
| `dense_100.h5` | several direct heap blocks |
| `dense_1000.h5` | a version-2 B-tree of depth two |
| `dense_order.h5` | `track_order=True`: creation-order flags in the link info message, creation order in each link message, and a two-byte creation index in every object header message |
| `dense_nested.h5` | 40 000 links to 100 objects: the root indirect block has more rows than the maximum direct block size, so it points at nested indirect blocks |
| `types.h5` | an `enum` dataset and a variable-length string dataset (global heap) |
| `compound.h5` | a compound dataset and a reference dataset, which the reader must skip while still reading the plain dataset beside them |
| `dense_attrs.h5` | fourteen attributes on one dataset, which forces the dense attribute layout (attribute info message, fractal heap, version-2 B-tree) |
| `chunked_filters.h5` | chunked dataset with shuffle, deflate, and fletcher32 filters |
| `v4_chunked.h5` | version-4 chunked layouts: single-chunk, fixed-array, extensible-array, and version-2 B-tree |
| `vlen_numeric.h5` | variable-length `int32` and `float64` datasets (global heap) |
| `soft_links.h5` | absolute and relative HDF5 soft links |
| `filters_plugins.h5` | plugin filters (`LZ4` 32004, `Zstd` 32015, `Blosc` 32001 with byte/bit-shuffle) |

## parquet/

| file | what it exercises |
| --- | --- |
| `lists_v1.parquet` | `list<double>` with a null list, an empty list and a null element, `list<string>`, and a `map<string,int64>`, written with version-1 data pages |
| `lists_v2.parquet` | the same data written with version-2 data pages |
| `codec_zstd.parquet`, `codec_lz4.parquet`, `codec_brotli.parquet` | 2000-row string/float table compressed with Zstandard, raw LZ4, and Brotli |
| `page_index.parquet` | 2000-row table with `ColumnIndex` and `OffsetIndex` |
| `decimal.parquet` | an `INT32` `DECIMAL(9,2)` column |
| `nested_lists.parquet`, `nested_lists_v2.parquet` | `list<list<int32>>` with a null list, an empty list, an empty inner list and a null element |
| `encodings_delta.parquet` | a `DELTA_BYTE_ARRAY` string column |
| `encodings_bss.parquet` | a `BYTE_STREAM_SPLIT` double column |

## netcdf/

| file | what it exercises |
| --- | --- |
| `cdf2.nc` | a classic-format file with 64-bit offsets (version byte 2) |
| `cdf5.nc` | a 64-bit data format file (version byte 5) |
