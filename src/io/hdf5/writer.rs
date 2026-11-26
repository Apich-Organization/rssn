//! Extended HDF5 writer: chunked datasets with filters, dense groups (v2/v3),
//! variable-length numeric datasets, and soft links.

use std::collections::BTreeMap;
use std::path::Path;

use super::Attribute;
use super::BTREE_NODE_SIZE;
use super::Dataset;
use super::H5Type;
use super::MAX_ARRAY_BYTES;
use super::SNOD_CAPACITY;
use super::SNOD_SIZE;
use super::align8;
use super::attribute_message;
use super::attributes_for;
use super::build_group_nodes;
use super::checksum::lookup3;
use super::dataspace_message;
use super::datatype_message;
use super::element_count;
use super::element_size;
use super::encode_data;
use super::encode_dataset_data;
use super::fill_value_message;
use super::filters;
use super::filters::Filter;
use super::group_btree;
use super::layout_message;
use super::local_heap;
use super::normalise_path;
use super::object_header_len;
use super::object_header_v1;
use super::path_segments;
use super::strides_of;
use super::string_datatype_message;
use super::superblock_v0;
use super::to_usize;

/// Filter specification for writing chunked HDF5 datasets.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum H5Filter {
    /// Deflate (zlib) compression at level `0..=9`.
    Deflate(u32),
    /// Byte-shuffle filter (element width is inferred from the dataset type).
    Shuffle,
    /// Fletcher-32 chunk checksum.
    Fletcher32,
    /// Szip (CCSDS 121.0) compression.
    Szip {
        /// Option mask (e.g. `32 | 8` for NN + LSB).
        options_mask: u32,
        /// Pixels per block (even, typically `16` or `32`).
        pixels_per_block: u32,
    },
    /// N-bit bit-packing filter.
    Nbit {
        /// Bit precision (`1..=64`).
        precision: u32,
        /// Bit offset (`0..`).
        offset: u32,
    },
    /// Scale-Offset filter.
    ScaleOffset {
        /// `0` = float D-scale, `1` = float E-scale, `2` = int.
        scale_type: u32,
        /// Decimal scale factor (for floats) or `0` (for ints).
        scale_factor: i32,
    },
    /// HDF5-LZ4 plugin filter (ID 32004).
    Lz4,
    /// HDF5-Zstd plugin filter (ID 32015) at level `1..=19`.
    Zstd(u32),
    /// HDF5-Blosc plugin filter (ID 32001).
    Blosc {
        /// Compression level `0..=9`.
        clevel: u32,
        /// `0` = none, `1` = byte-shuffle, `2` = bit-shuffle.
        shuffle: u32,
        /// `0` = blosclz, `1` = lz4, `2` = lz4hc, `3` = snappy, `4` = zlib, `5` = zstd.
        compressor: u32,
    },
}

/// Per-dataset write options.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct DatasetWriteOptions {
    /// Chunk shape; when `Some`, the dataset is written with chunked layout.
    pub chunk_dims: Option<Vec<usize>>,
    /// Filter pipeline applied to each chunk (requires `chunk_dims`).
    pub filters: Vec<H5Filter>,
}

/// A symbolic (soft) link inside an HDF5 file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SoftLink {
    /// Absolute path of the link inside the file, e.g. `"/alias"`.
    pub path: String,
    /// Target path string stored in the link, e.g. `"/data"` or `"sub/data"`.
    pub target: String,
}

/// Options controlling extended HDF5 writing.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct H5WriteOptions {
    /// Emit version-2 superblock, object headers v2, and dense groups
    /// (fractal heap + version-2 B-tree) instead of v0/v1 symbol tables.
    pub dense_groups: bool,
    /// Per-dataset chunking and filter options keyed by dataset path.
    pub dataset_options: Vec<(String, DatasetWriteOptions)>,
    /// Soft links to emit in the group hierarchy.
    pub soft_links: Vec<SoftLink>,
}

/// Writes an HDF5 file to disk with extended options.
///
/// # Errors
/// Fails when paths collide, values do not fit their declared types, filter/chunk
/// configuration is invalid, or the file cannot be written.
pub fn write_hdf5_with_options(
    path: &Path,
    datasets: &[Dataset],
    attributes: &[(String, Vec<Attribute>)],
    options: &H5WriteOptions,
) -> Result<(), String> {
    let bytes = write_hdf5_bytes_with_options(datasets, attributes, options)?;
    std::fs::write(path, bytes)
        .map_err(|error| format!("failed to write {}: {error}", path.display()))
}

/// Serialises an HDF5 file image with extended options.
///
/// # Errors
/// Fails when paths collide, values do not fit their declared types, or filter/chunk
/// configuration is invalid.
pub fn write_hdf5_bytes_with_options(
    datasets: &[Dataset],
    attributes: &[(String, Vec<Attribute>)],
    options: &H5WriteOptions,
) -> Result<Vec<u8>, String> {
    if options.dense_groups {
        write_dense_hdf5(datasets, attributes, options)
    } else {
        write_v0_extended_hdf5(datasets, attributes, options)
    }
}

// -------------------------------------------------------------------------
// Shared dataset preparation (contiguous, chunked, vlen-numeric)
// -------------------------------------------------------------------------

enum PreparedPayload {
    Contiguous {
        data: Vec<u8>,
        gcol: Option<Vec<u8>>,
    },
    Chunked {
        chunk_dims: Vec<usize>,
        filters: Vec<Filter>,
        chunks: Vec<(Vec<u64>, Vec<u8>)>,
    },
}

fn dataset_opts_for<'a>(
    options: &'a H5WriteOptions,
    path: &str,
) -> Option<&'a DatasetWriteOptions> {
    let norm = normalise_path(path);
    options
        .dataset_options
        .iter()
        .find(|(p, _)| normalise_path(p) == norm)
        .map(|(_, o)| o)
}

fn resolve_filters(
    h5_filters: &[H5Filter],
    elem_bytes: usize,
    is_signed: bool,
    is_float: bool,
    chunk_elems: usize,
) -> Vec<Filter> {
    h5_filters
        .iter()
        .map(|f| {
            match f {
                | H5Filter::Deflate(level) => {
                    Filter {
                        id: filters::FILTER_DEFLATE,
                        client: vec![*level],
                    }
                },
                | H5Filter::Shuffle => {
                    Filter {
                        id: filters::FILTER_SHUFFLE,
                        client: vec![elem_bytes as u32],
                    }
                },
                | H5Filter::Fletcher32 => {
                    Filter {
                        id: filters::FILTER_FLETCHER32,
                        client: Vec::new(),
                    }
                },
                | H5Filter::Szip {
                    options_mask,
                    pixels_per_block,
                } => {
                    Filter {
                        id: filters::FILTER_SZIP,
                        client: vec![
                            *options_mask,
                            *pixels_per_block,
                            (elem_bytes * 8) as u32,
                            chunk_elems as u32,
                        ],
                    }
                },
                | H5Filter::Nbit { precision, offset } => {
                    Filter {
                        id: filters::FILTER_NBIT,
                        client: vec![
                            7,
                            0,
                            u32::from(is_float),
                            elem_bytes as u32,
                            0,
                            *precision,
                            *offset,
                            u32::from(is_signed),
                        ],
                    }
                },
                | H5Filter::ScaleOffset {
                    scale_type,
                    scale_factor,
                } => {
                    Filter {
                        id: filters::FILTER_SCALEOFFSET,
                        client: vec![
                            *scale_type,
                            *scale_factor as u32,
                            chunk_elems as u32,
                            u32::from(is_float),
                            elem_bytes as u32,
                            0,
                            u32::from(is_signed),
                        ],
                    }
                },
                | H5Filter::Lz4 => {
                    Filter {
                        id: filters::FILTER_LZ4,
                        client: vec![0],
                    }
                },
                | H5Filter::Zstd(level) => {
                    Filter {
                        id: filters::FILTER_ZSTD,
                        client: vec![*level],
                    }
                },
                | H5Filter::Blosc {
                    clevel,
                    shuffle,
                    compressor,
                } => {
                    Filter {
                        id: filters::FILTER_BLOSC,
                        client: vec![
                            2,
                            2,
                            elem_bytes as u32,
                            (chunk_elems * elem_bytes) as u32,
                            *clevel,
                            *shuffle,
                            *compressor,
                        ],
                    }
                },
            }
        })
        .collect()
}

fn vlen_numeric_datatype_message(base: H5Type) -> Result<Vec<u8>, String> {
    let base_msg = datatype_message(base)?;
    let mut out = Vec::with_capacity(8 + base_msg.len());
    // Class 9 (vlen), version 1, kind 0 (sequence)
    out.push(0x19);
    out.push(0x00);
    out.extend_from_slice(&[0, 0]);
    out.extend_from_slice(&16_u32.to_le_bytes());
    out.extend_from_slice(&base_msg);
    Ok(out)
}

fn dataset_datatype_msg_ext(dataset: &Dataset) -> Result<Vec<u8>, String> {
    match dataset.h5_type {
        | H5Type::Str => Ok(string_datatype_message(super::string_width(dataset))),
        | H5Type::VlenNumeric => {
            let base = dataset.vlen_base_type.unwrap_or(H5Type::F64);
            vlen_numeric_datatype_message(base)
        },
        | other => datatype_message(other),
    }
}

/// Builds a Global Heap Collection (`GCOL`) from raw objects, returning the
/// collection bytes and the `(obj_index, elem_count)` for each entry.
fn build_gcol(
    objects: &[Vec<u8>],
    elem_counts: &[u32],
    gcol_addr: u64,
) -> (Vec<u8>, Vec<u8>) {
    let mut body = Vec::new();
    let mut descriptors = Vec::with_capacity(objects.len() * 16);
    for (idx, (obj, &seq_len)) in objects.iter().zip(elem_counts.iter()).enumerate() {
        let obj_index = (idx + 1) as u16;
        if seq_len == 0 {
            descriptors.extend_from_slice(&0_u32.to_le_bytes());
            descriptors.extend_from_slice(&0_u64.to_le_bytes());
            descriptors.extend_from_slice(&0_u32.to_le_bytes());
            continue;
        }
        body.extend_from_slice(&obj_index.to_le_bytes());
        body.extend_from_slice(&1_u16.to_le_bytes()); // refcount
        body.extend_from_slice(&[0; 4]);
        body.extend_from_slice(&(obj.len() as u64).to_le_bytes());
        body.extend_from_slice(obj);
        let pad = (8 - (obj.len() % 8)) % 8;
        body.extend(std::iter::repeat_n(0_u8, pad));

        descriptors.extend_from_slice(&seq_len.to_le_bytes());
        descriptors.extend_from_slice(&gcol_addr.to_le_bytes());
        descriptors.extend_from_slice(&u32::from(obj_index).to_le_bytes());
    }
    // Object 0 (free space) must account for the remaining space in a >=4096-byte GCOL.
    let used = 16 + body.len();
    let total_size = ((used + 16).max(4096) + 7) & !7;
    let free_space = total_size - used;
    body.extend_from_slice(&0_u16.to_le_bytes());
    body.extend_from_slice(&0_u16.to_le_bytes());
    body.extend_from_slice(&[0_u8; 4]);
    body.extend_from_slice(&(free_space as u64).to_le_bytes());
    body.extend(std::iter::repeat_n(0_u8, free_space - 16));

    let mut gcol = Vec::with_capacity(total_size);
    gcol.extend_from_slice(b"GCOL");
    gcol.push(1); // version 1
    gcol.extend_from_slice(&[0; 3]);
    gcol.extend_from_slice(&(total_size as u64).to_le_bytes());
    gcol.extend_from_slice(&body);
    (gcol, descriptors)
}

fn prepare_dataset_payload(
    dataset: &Dataset,
    opts: Option<&DatasetWriteOptions>,
    gcol_addr_placeholder: u64,
) -> Result<(PreparedPayload, Option<(Vec<Vec<u8>>, Vec<u32>)>), String> {
    let count = element_count(&dataset.shape)?;
    let count = to_usize(count, "dataset element count")?;

    let (raw, elem_bytes, is_signed, is_float) = match dataset.h5_type {
        | H5Type::Str => {
            if count != dataset.strings.len() {
                return Err(format!(
                    "dataset {} has shape {:?} but {} strings",
                    dataset.path,
                    dataset.shape,
                    dataset.strings.len()
                ));
            }
            let w = super::string_width(dataset);
            (encode_dataset_data(dataset)?, w, false, false)
        },
        | H5Type::VlenNumeric => {
            if count != dataset.vlen_data.len() {
                return Err(format!(
                    "dataset {} has shape {:?} but {} vlen sequences",
                    dataset.path,
                    dataset.shape,
                    dataset.vlen_data.len()
                ));
            }
            let base = dataset.vlen_base_type.unwrap_or(H5Type::F64);
            let mut objs = Vec::with_capacity(count);
            let mut counts = Vec::with_capacity(count);
            for seq in &dataset.vlen_data {
                objs.push(encode_data(base, seq)?);
                counts.push(seq.len() as u32);
            }
            let (gcol, descs) = build_gcol(&objs, &counts, gcol_addr_placeholder);
            return Ok((
                PreparedPayload::Contiguous {
                    data: descs,
                    gcol: Some(gcol),
                },
                Some((objs, counts)),
            ));
        },
        | other => {
            if count != dataset.data.len() {
                return Err(format!(
                    "dataset {} has shape {:?} but {} values",
                    dataset.path,
                    dataset.shape,
                    dataset.data.len()
                ));
            }
            let eb = element_size(dataset);
            let signed = matches!(other, H5Type::I8 | H5Type::I16 | H5Type::I32 | H5Type::I64);
            let float = matches!(other, H5Type::F32 | H5Type::F64);
            (encode_dataset_data(dataset)?, eb, signed, float)
        },
    };

    if raw.len() > MAX_ARRAY_BYTES {
        return Err(format!(
            "dataset {} exceeds the {MAX_ARRAY_BYTES}-byte limit",
            dataset.path
        ));
    }

    let Some(ds_opts) = opts else {
        return Ok((
            PreparedPayload::Contiguous {
                data: raw,
                gcol: None,
            },
            None,
        ));
    };
    let Some(chunk_dims) = &ds_opts.chunk_dims else {
        if !ds_opts.filters.is_empty() {
            return Err(format!(
                "dataset {} specifies filters without chunk_dims",
                dataset.path
            ));
        }
        return Ok((
            PreparedPayload::Contiguous {
                data: raw,
                gcol: None,
            },
            None,
        ));
    };

    let rank = dataset.shape.len();
    if rank == 0 || chunk_dims.len() != rank || chunk_dims.iter().any(|&d| d == 0) {
        return Err(format!(
            "dataset {} has invalid chunk_dims {:?} for shape {:?}",
            dataset.path, chunk_dims, dataset.shape
        ));
    }
    let chunk_elems: usize = chunk_dims.iter().product();
    let filters = resolve_filters(
        &ds_opts.filters,
        elem_bytes,
        is_signed,
        is_float,
        chunk_elems,
    );

    // Enumerate all n-D chunk origins
    let grid_counts: Vec<usize> = dataset
        .shape
        .iter()
        .zip(chunk_dims.iter())
        .map(|(&s, &c)| s.div_ceil(c))
        .collect();
    let total_chunks: usize = grid_counts.iter().product();
    let ds_strides = strides_of(&dataset.shape);
    let mut chunks = Vec::with_capacity(total_chunks);

    for chunk_idx in 0..total_chunks {
        let mut rem = chunk_idx;
        let mut origin = vec![0_usize; rank];
        for axis in (0..rank).rev() {
            let g = grid_counts[axis];
            origin[axis] = (rem % g) * chunk_dims[axis];
            rem /= g;
        }
        let mut chunk_buf = vec![0_u8; chunk_elems * elem_bytes];
        for flat in 0..chunk_elems {
            let mut c_rem = flat;
            let mut inside = true;
            let mut src_elem = 0_usize;
            for axis in (0..rank).rev() {
                let coord = c_rem % chunk_dims[axis];
                c_rem /= chunk_dims[axis];
                let pos = origin[axis] + coord;
                if pos >= dataset.shape[axis] {
                    inside = false;
                    break;
                }
                src_elem += pos * ds_strides[axis];
            }
            if inside {
                let s0 = src_elem * elem_bytes;
                let d0 = flat * elem_bytes;
                chunk_buf[d0..d0 + elem_bytes].copy_from_slice(&raw[s0..s0 + elem_bytes]);
            }
        }
        let encoded_chunk = filters::apply_filters_encode(chunk_buf, &filters)?;
        let mut offsets: Vec<u64> = origin.iter().map(|&o| o as u64).collect();
        offsets.push(0); // trailing dimension for element size
        chunks.push((offsets, encoded_chunk));
    }

    let mut full_chunk_dims = chunk_dims.clone();
    full_chunk_dims.push(elem_bytes);

    Ok((
        PreparedPayload::Chunked {
            chunk_dims: full_chunk_dims,
            filters,
            chunks,
        },
        None,
    ))
}

fn chunked_layout_message(
    btree_addr: u64,
    chunk_dims: &[usize],
) -> Vec<u8> {
    let mut out = Vec::with_capacity(3 + 8 + chunk_dims.len() * 4);
    out.push(3); // version 3
    out.push(2); // class 2 = chunked
    out.push(chunk_dims.len() as u8);
    out.extend_from_slice(&btree_addr.to_le_bytes());
    for &d in chunk_dims {
        out.extend_from_slice(&(d as u32).to_le_bytes());
    }
    out
}

fn filter_pipeline_message(filters: &[Filter]) -> Vec<u8> {
    // Version 2 filter pipeline (compact, no 8-byte name alignment padding)
    let mut out = Vec::new();
    out.push(2); // version 2
    out.push(filters.len() as u8);
    for f in filters {
        out.extend_from_slice(&f.id.to_le_bytes());
        if f.id >= 256 {
            out.extend_from_slice(&0_u16.to_le_bytes()); // name_length = 0
        }
        out.extend_from_slice(&0_u16.to_le_bytes()); // flags = mandatory
        out.extend_from_slice(&(f.client.len() as u16).to_le_bytes());
        for &c in &f.client {
            out.extend_from_slice(&c.to_le_bytes());
        }
    }
    out
}

/// Builds a single-level or multi-level version-1 chunk B-tree (`TREE`, node type 1).
fn build_chunk_btree(
    chunk_dims: &[usize],
    chunks: &[(Vec<u64>, Vec<u8>, u64)], // (offsets, data, file_addr)
    btree_addr: u64,
) -> Vec<u8> {
    // For simplicity and full compatibility, emit a leaf (level 0) if chunks.len() <= 2*K,
    // or a 2-level tree if larger. Using 2*K = chunks.len().max(16) in a single leaf node
    // is completely valid in HDF5 v1 B-trees (the B-tree header K can be large, or we can
    // use a standard fixed K=32 and build parent nodes if needed).
    // Let's use K = (chunks.len() / 2).max(16) so a single leaf node holds all entries
    // cleanly, or split if > 32768.
    let ndims = chunk_dims.len();
    let key_size = 8 + 8 * ndims;
    let entry_count = chunks.len();
    let k = entry_count.div_ceil(2).max(32);
    let node_size = 24 + (2 * k) * (key_size + 8) + key_size;
    let mut node = vec![0_u8; node_size];
    node[..4].copy_from_slice(b"TREE");
    node[4] = 1; // chunk B-tree
    node[5] = 0; // level 0 (leaf)
    node[6..8].copy_from_slice(&(entry_count as u16).to_le_bytes());
    node[8..16].copy_from_slice(&u64::MAX.to_le_bytes()); // left sibling
    node[16..24].copy_from_slice(&u64::MAX.to_le_bytes()); // right sibling

    let mut sorted: Vec<&(Vec<u64>, Vec<u8>, u64)> = chunks.iter().collect();
    sorted.sort_by(|a, b| a.0.cmp(&b.0));

    let mut pos = 24_usize;
    for (offsets, data, addr) in &sorted {
        node[pos..pos + 4].copy_from_slice(&(data.len() as u32).to_le_bytes());
        node[pos + 4..pos + 8].copy_from_slice(&0_u32.to_le_bytes()); // filter mask = 0
        for d in 0..ndims {
            let off = offsets.get(d).copied().unwrap_or(0);
            node[pos + 8 + d * 8..pos + 16 + d * 8].copy_from_slice(&off.to_le_bytes());
        }
        pos += key_size;
        node[pos..pos + 8].copy_from_slice(&addr.to_le_bytes());
        pos += 8;
    }
    // Upper-bound sentinel key[entry_count] so HDF5 H5B_find matches the final chunk.
    if let Some((last_offsets, _, _)) = sorted.last() {
        for d in 0..ndims.saturating_sub(1) {
            let bound = last_offsets.get(d).copied().unwrap_or(0) + chunk_dims[d] as u64;
            node[pos + 8 + d * 8..pos + 16 + d * 8].copy_from_slice(&bound.to_le_bytes());
        }
    }
    let _ = btree_addr;
    node
}

// -------------------------------------------------------------------------
// v0 Writer with Chunked, VlenNumeric, and SoftLink support
// -------------------------------------------------------------------------

#[derive(Clone, Copy)]
enum ExtChild {
    Group(usize),
    Dataset(usize),
    SoftLink(usize),
}

struct ExtGroupPlan {
    path: String,
    children: Vec<ExtChild>,
    entries: Vec<(String, u64, u64, u32)>, // name, name_off, obj_addr_or_target_off, cache_type
    heap: Vec<u8>,
    nodes: Vec<super::GroupNode>,
    root_index: usize,
    snod_count: usize,
    header_len: usize,
    header_addr: u64,
    heap_addr: u64,
    btree_addr: u64,
    snod_addr: u64,
}

/// Writes a symbol table node (`SNOD`) supporting both hard links (`cache_type = 0`)
/// and soft links (`cache_type = 2`, `target_offset` stored in `obj_addr` slot).
fn symbol_table_node_ext(entries: &[(String, u64, u64, u32)]) -> Vec<u8> {
    let mut out = vec![0_u8; SNOD_SIZE];
    out[..4].copy_from_slice(b"SNOD");
    out[4] = 1;
    out[6..8].copy_from_slice(&u16::try_from(entries.len()).unwrap_or(0).to_le_bytes());
    for (index, (_, name_offset, target, cache_type)) in entries.iter().enumerate() {
        let pos = 8 + index * 40;
        out[pos..pos + 8].copy_from_slice(&name_offset.to_le_bytes());
        if *cache_type == 2 {
            out[pos + 8..pos + 16].copy_from_slice(&u64::MAX.to_le_bytes()); // undefined obj header
            out[pos + 16..pos + 20].copy_from_slice(&2_u32.to_le_bytes()); // cache_type = 2 (soft link)
            // scratchpad: 4-byte offset into local heap for soft link value
            out[pos + 24..pos + 28].copy_from_slice(&(*target as u32).to_le_bytes());
        } else {
            out[pos + 8..pos + 16].copy_from_slice(&target.to_le_bytes());
        }
    }
    out
}

fn write_v0_extended_hdf5(
    datasets: &[Dataset],
    attributes: &[(String, Vec<Attribute>)],
    options: &H5WriteOptions,
) -> Result<Vec<u8>, String> {
    let mut prepared = Vec::with_capacity(datasets.len());
    let mut vlen_raws = Vec::with_capacity(datasets.len());
    for ds in datasets {
        let (p, v) = prepare_dataset_payload(ds, dataset_opts_for(options, &ds.path), 0)?;
        prepared.push(p);
        vlen_raws.push(v);
    }

    // Build group tree
    let mut groups: Vec<ExtGroupPlan> = vec![ExtGroupPlan {
        path: "/".to_owned(),
        children: Vec::new(),
        entries: Vec::new(),
        heap: Vec::new(),
        nodes: Vec::new(),
        root_index: 0,
        snod_count: 1,
        header_len: 0,
        header_addr: 0,
        heap_addr: 0,
        btree_addr: 0,
        snod_addr: 0,
    }];
    let mut group_index: BTreeMap<String, usize> = BTreeMap::new();
    group_index.insert("/".to_owned(), 0);
    let mut dataset_paths = Vec::with_capacity(datasets.len());

    fn ensure_v0_parent(
        groups: &mut Vec<ExtGroupPlan>,
        group_index: &mut BTreeMap<String, usize>,
        full_path: &str,
    ) -> Result<(usize, String), String> {
        let segments = path_segments(full_path)?;
        let mut parent = 0_usize;
        let mut cur = String::new();
        for segment in &segments[..segments.len() - 1] {
            cur.push('/');
            cur.push_str(segment);
            parent = match group_index.get(&cur) {
                | Some(&idx) => idx,
                | None => {
                    let idx = groups.len();
                    groups.push(ExtGroupPlan {
                        path: cur.clone(),
                        children: Vec::new(),
                        entries: Vec::new(),
                        heap: Vec::new(),
                        nodes: Vec::new(),
                        root_index: 0,
                        snod_count: 1,
                        header_len: 0,
                        header_addr: 0,
                        heap_addr: 0,
                        btree_addr: 0,
                        snod_addr: 0,
                    });
                    group_index.insert(cur.clone(), idx);
                    groups[parent].children.push(ExtChild::Group(idx));
                    idx
                },
            };
        }
        let full = format!("{cur}/{}", segments[segments.len() - 1]);
        Ok((parent, full))
    }

    for (pos, ds) in datasets.iter().enumerate() {
        let (parent, full) = ensure_v0_parent(&mut groups, &mut group_index, &ds.path)?;
        if group_index.contains_key(&full) {
            return Err(format!("dataset {} collides with a group", ds.path));
        }
        if dataset_paths.contains(&full) {
            return Err(format!("duplicate dataset name {full}"));
        }
        dataset_paths.push(full);
        groups[parent].children.push(ExtChild::Dataset(pos));
    }

    let mut soft_paths = Vec::with_capacity(options.soft_links.len());
    for (pos, link) in options.soft_links.iter().enumerate() {
        let (parent, full) = ensure_v0_parent(&mut groups, &mut group_index, &link.path)?;
        soft_paths.push(full);
        groups[parent].children.push(ExtChild::SoftLink(pos));
    }

    // Build per-group local heap + B-tree
    let child_names: Vec<Vec<String>> = groups
        .iter()
        .map(|g| {
            g.children
                .iter()
                .map(|c| {
                    match c {
                        | ExtChild::Group(i) => {
                            groups[*i].path.rsplit('/').next().unwrap_or("").to_owned()
                        },
                        | ExtChild::Dataset(i) => {
                            dataset_paths[*i]
                                .rsplit('/')
                                .next()
                                .unwrap_or("")
                                .to_owned()
                        },
                        | ExtChild::SoftLink(i) => {
                            soft_paths[*i].rsplit('/').next().unwrap_or("").to_owned()
                        },
                    }
                })
                .collect()
        })
        .collect();

    for index in 0..groups.len() {
        let mut order: Vec<usize> = (0..groups[index].children.len()).collect();
        order.sort_by(|a, b| child_names[index][*a].cmp(&child_names[index][*b]));
        let children: Vec<ExtChild> = order.iter().map(|&p| groups[index].children[p]).collect();
        let names: Vec<String> = order
            .iter()
            .map(|&p| child_names[index][p].clone())
            .collect();
        if names.windows(2).any(|w| w[0] == w[1]) {
            return Err(format!(
                "group {} has two children with the same name",
                groups[index].path
            ));
        }

        // Put both child names and soft-link target strings into the local heap
        let mut heap_strings = names.clone();
        let mut soft_target_idx = vec![None; children.len()];
        for (ci, child) in children.iter().enumerate() {
            if let ExtChild::SoftLink(si) = child {
                soft_target_idx[ci] = Some(heap_strings.len());
                heap_strings.push(options.soft_links[*si].target.clone());
            }
        }
        let (heap, offsets) = local_heap(&heap_strings);
        let entries: Vec<(String, u64, u64, u32)> = names
            .iter()
            .enumerate()
            .map(|(ci, name)| {
                if let Some(ti) = soft_target_idx[ci] {
                    (name.clone(), offsets[ci], offsets[ti], 2)
                } else {
                    (name.clone(), offsets[ci], 0, 0)
                }
            })
            .collect();
        let btree_entries: Vec<(String, u64, u64)> = entries
            .iter()
            .map(|(n, no, t, _)| (n.clone(), *no, *t))
            .collect();
        let (nodes, root_index) = build_group_nodes(&btree_entries);
        let snod_count = entries.chunks(SNOD_CAPACITY).count().max(1);
        let g = &mut groups[index];
        g.children = children;
        g.heap = heap;
        g.entries = entries;
        g.nodes = nodes;
        g.root_index = root_index;
        g.snod_count = snod_count;
    }

    // Header sizes
    for g in &mut groups {
        let attrs = attributes_for(attributes, &g.path);
        let mut messages = vec![(0x11_u16, vec![0_u8; 16])];
        for attr in &attrs {
            messages.push((0x0C, attribute_message(attr)?));
        }
        g.header_len = object_header_len(&messages);
    }

    let mut ds_header_lens = Vec::with_capacity(datasets.len());
    for (i, ds) in datasets.iter().enumerate() {
        let attrs = attributes_for(attributes, &ds.path);
        let mut msgs = vec![
            (0x01, dataspace_message(&ds.shape)),
            (0x03, dataset_datatype_msg_ext(ds)?),
            (0x05, fill_value_message()),
        ];
        match &prepared[i] {
            | PreparedPayload::Contiguous { data, .. } => {
                msgs.push((0x08, layout_message(0, data.len() as u64)));
            },
            | PreparedPayload::Chunked {
                chunk_dims, filters, ..
            } => {
                if !filters.is_empty() {
                    msgs.push((0x0B, filter_pipeline_message(filters)));
                }
                msgs.push((0x08, chunked_layout_message(0, chunk_dims)));
            },
        }
        for attr in &attrs {
            msgs.push((0x0C, attribute_message(attr)?));
        }
        ds_header_lens.push(object_header_len(&msgs));
    }

    // Assign file addresses
    let mut pos: u64 = 96;
    for g in &mut groups {
        pos = align8(pos);
        g.header_addr = pos;
        pos += g.header_len as u64;
        g.heap_addr = pos;
        if let Some(slot) = g.heap.get_mut(24..32) {
            slot.copy_from_slice(&(g.heap_addr + 32).to_le_bytes());
        }
        pos += g.heap.len() as u64;
        g.btree_addr = pos;
        pos += (g.nodes.len() * BTREE_NODE_SIZE) as u64;
        g.snod_addr = pos;
        pos += (g.snod_count * SNOD_SIZE) as u64;
    }

    enum DsLayoutAddr {
        Contiguous {
            header_addr: u64,
            data_addr: u64,
            gcol_addr: Option<u64>,
        },
        Chunked {
            header_addr: u64,
            btree_addr: u64,
            chunk_addrs: Vec<u64>,
        },
    }

    let mut ds_addrs = Vec::with_capacity(datasets.len());
    for (i, payload) in prepared.iter_mut().enumerate() {
        pos = align8(pos);
        let header_addr = pos;
        pos += ds_header_lens[i] as u64;
        match payload {
            | PreparedPayload::Contiguous { data, gcol } => {
                let gcol_addr = if let Some(g) = gcol {
                    pos = align8(pos);
                    let ga = pos;
                    pos += g.len() as u64;
                    if let Some((objs, counts)) = &vlen_raws[i] {
                        let (new_gcol, new_descs) = build_gcol(objs, counts, ga);
                        *g = new_gcol;
                        *data = new_descs;
                    }
                    Some(ga)
                } else {
                    None
                };
                pos = align8(pos);
                let data_addr = pos;
                pos += data.len() as u64;
                ds_addrs.push(DsLayoutAddr::Contiguous {
                    header_addr,
                    data_addr,
                    gcol_addr,
                });
            },
            | PreparedPayload::Chunked {
                chunk_dims, chunks, ..
            } => {
                let mut chunk_addrs = Vec::with_capacity(chunks.len());
                for (_, chunk_bytes) in chunks.iter() {
                    pos = align8(pos);
                    chunk_addrs.push(pos);
                    pos += chunk_bytes.len() as u64;
                }
                pos = align8(pos);
                let btree_addr = pos;
                let dummy_chunks: Vec<(Vec<u64>, Vec<u8>, u64)> = chunks
                    .iter()
                    .zip(chunk_addrs.iter())
                    .map(|((o, d), &a)| (o.clone(), d.clone(), a))
                    .collect();
                let btree_len = build_chunk_btree(chunk_dims, &dummy_chunks, btree_addr).len();
                pos += btree_len as u64;
                ds_addrs.push(DsLayoutAddr::Chunked {
                    header_addr,
                    btree_addr,
                    chunk_addrs,
                });
            },
        }
    }
    let eof = align8(pos);

    // Resolve hard-link addresses in group entries
    for idx in 0..groups.len() {
        let resolved: Vec<Option<u64>> = groups[idx]
            .children
            .iter()
            .map(|c| {
                match c {
                    | ExtChild::Group(t) => Some(groups[*t].header_addr),
                    | ExtChild::Dataset(t) => {
                        Some(match &ds_addrs[*t] {
                            | DsLayoutAddr::Contiguous { header_addr, .. } => *header_addr,
                            | DsLayoutAddr::Chunked { header_addr, .. } => *header_addr,
                        })
                    },
                    | ExtChild::SoftLink(_) => None,
                }
            })
            .collect();
        for (entry, addr_opt) in groups[idx].entries.iter_mut().zip(resolved) {
            if let Some(addr) = addr_opt {
                entry.2 = addr;
            }
        }
    }

    // Emit bytes
    let root = &groups[0];
    let mut out = Vec::with_capacity(eof as usize);
    out.extend_from_slice(&superblock_v0(
        root.header_addr,
        root.btree_addr,
        root.heap_addr,
        eof,
    ));

    for g in &groups {
        while out.len() < g.header_addr as usize {
            out.push(0);
        }
        let attrs = attributes_for(attributes, &g.path);
        let mut messages = vec![(
            0x11_u16,
            [g.btree_addr.to_le_bytes(), g.heap_addr.to_le_bytes()].concat(),
        )];
        for attr in &attrs {
            messages.push((0x0C, attribute_message(attr)?));
        }
        out.extend_from_slice(&object_header_v1(&messages)?);
        while out.len() < g.heap_addr as usize {
            out.push(0);
        }
        out.extend_from_slice(&g.heap);
        out.extend_from_slice(&group_btree(
            &g.nodes,
            g.root_index,
            g.btree_addr,
            g.snod_addr,
        ));
        for chunk in g.entries.chunks(SNOD_CAPACITY) {
            out.extend_from_slice(&symbol_table_node_ext(chunk));
        }
    }

    for (i, ds) in datasets.iter().enumerate() {
        let attrs = attributes_for(attributes, &ds.path);
        let mut msgs = vec![
            (0x01, dataspace_message(&ds.shape)),
            (0x03, dataset_datatype_msg_ext(ds)?),
            (0x05, fill_value_message()),
        ];
        match (&prepared[i], &ds_addrs[i]) {
            | (
                PreparedPayload::Contiguous { data, gcol },
                DsLayoutAddr::Contiguous {
                    header_addr,
                    data_addr,
                    gcol_addr,
                },
            ) => {
                msgs.push((0x08, layout_message(*data_addr, data.len() as u64)));
                for attr in &attrs {
                    msgs.push((0x0C, attribute_message(attr)?));
                }
                while out.len() < *header_addr as usize {
                    out.push(0);
                }
                out.extend_from_slice(&object_header_v1(&msgs)?);
                if let (Some(g), Some(ga)) = (gcol, gcol_addr) {
                    while out.len() < *ga as usize {
                        out.push(0);
                    }
                    out.extend_from_slice(g);
                }
                while out.len() < *data_addr as usize {
                    out.push(0);
                }
                out.extend_from_slice(data);
            },
            | (
                PreparedPayload::Chunked {
                    chunk_dims,
                    filters,
                    chunks,
                    ..
                },
                DsLayoutAddr::Chunked {
                    header_addr,
                    btree_addr,
                    chunk_addrs,
                },
            ) => {
                if !filters.is_empty() {
                    msgs.push((0x0B, filter_pipeline_message(filters)));
                }
                msgs.push((0x08, chunked_layout_message(*btree_addr, chunk_dims)));
                for attr in &attrs {
                    msgs.push((0x0C, attribute_message(attr)?));
                }
                while out.len() < *header_addr as usize {
                    out.push(0);
                }
                out.extend_from_slice(&object_header_v1(&msgs)?);
                let mut btree_items = Vec::with_capacity(chunks.len());
                for ((offsets, cbytes), &caddr) in chunks.iter().zip(chunk_addrs.iter()) {
                    while out.len() < caddr as usize {
                        out.push(0);
                    }
                    out.extend_from_slice(cbytes);
                    btree_items.push((offsets.clone(), cbytes.clone(), caddr));
                }
                while out.len() < *btree_addr as usize {
                    out.push(0);
                }
                out.extend_from_slice(&build_chunk_btree(chunk_dims, &btree_items, *btree_addr));
            },
            | _ => unreachable!(),
        }
    }

    while out.len() < eof as usize {
        out.push(0);
    }
    Ok(out)
}

// -------------------------------------------------------------------------
// Dense Group Writer (Superblock v2, OHDR v2, FRHP/FHDB, BTHD/BTLF)
// -------------------------------------------------------------------------

fn superblock_v2(
    root_addr: u64,
    eof: u64,
) -> Vec<u8> {
    let mut out = Vec::with_capacity(48);
    out.extend_from_slice(b"\x89HDF\r\n\x1a\n");
    out.push(2); // version 2
    out.push(8); // size of offsets
    out.push(8); // size of lengths
    out.push(0); // file consistency flags
    out.extend_from_slice(&0_u64.to_le_bytes()); // base address
    out.extend_from_slice(&u64::MAX.to_le_bytes()); // superblock extension address
    out.extend_from_slice(&eof.to_le_bytes()); // EOF address
    out.extend_from_slice(&root_addr.to_le_bytes()); // root group OHDR address
    let cksum = lookup3(&out);
    out.extend_from_slice(&cksum.to_le_bytes());
    out
}

fn object_header_v2(messages: &[(u8, Vec<u8>)]) -> Result<Vec<u8>, String> {
    let mut chunk_data = Vec::new();
    for (kind, data) in messages {
        let len = u16::try_from(data.len())
            .map_err(|_| "OHDR v2 message is larger than 64 KiB".to_owned())?;
        chunk_data.push(*kind);
        chunk_data.extend_from_slice(&len.to_le_bytes());
        chunk_data.push(0); // message flags
        chunk_data.extend_from_slice(data);
    }
    // Flags: 0x03 => 8-byte chunk#0 size? No: bits 0-1: 0=1B, 1=2B, 2=4B, 3=8B.
    // Let's always use 4-byte chunk#0 size (flags = 0x02) so size is predictable.
    let mut out = Vec::with_capacity(15 + chunk_data.len());
    out.extend_from_slice(b"OHDR");
    out.push(2); // version 2
    out.push(0x02); // 4-byte chunk size, no timestamps, no attr phase change
    out.extend_from_slice(&(chunk_data.len() as u32).to_le_bytes());
    out.extend_from_slice(&chunk_data);
    let cksum = lookup3(&out);
    out.extend_from_slice(&cksum.to_le_bytes());
    Ok(out)
}

fn link_message_hard(
    name: &str,
    target_addr: u64,
) -> Vec<u8> {
    let nb = name.as_bytes();
    let mut out = Vec::with_capacity(12 + nb.len());
    out.push(1); // version 1
    out.push(0x01); // flags: 2-byte name length, hard link (link_type=0 omitted)
    out.extend_from_slice(&(nb.len() as u16).to_le_bytes());
    out.extend_from_slice(nb);
    out.extend_from_slice(&target_addr.to_le_bytes());
    out
}

fn link_message_soft(
    name: &str,
    target_path: &str,
) -> Vec<u8> {
    let nb = name.as_bytes();
    let tb = target_path.as_bytes();
    let mut out = Vec::with_capacity(7 + nb.len() + tb.len());
    out.push(1); // version 1
    out.push(0x09); // flags: bit 0 (2-byte name len) | bit 3 (link_type present)
    out.push(1); // link_type = 1 (soft link)
    out.extend_from_slice(&(nb.len() as u16).to_le_bytes());
    out.extend_from_slice(nb);
    out.extend_from_slice(&(tb.len() as u16).to_le_bytes());
    out.extend_from_slice(tb);
    out
}

fn link_info_message(
    fheap_addr: u64,
    btree_addr: u64,
) -> Vec<u8> {
    let mut out = Vec::with_capacity(26);
    out.push(0); // version 0
    out.push(0); // flags: no creation order tracking
    out.extend_from_slice(&fheap_addr.to_le_bytes());
    out.extend_from_slice(&btree_addr.to_le_bytes());
    out
}

/// Builds a single-direct-block Fractal Heap (`FRHP` + `FHDB`) storing `objects`.
/// Returns `(frhp_bytes, fhdb_bytes, heap_ids)`.
fn build_fractal_heap_single_block(
    objects: &[Vec<u8>],
    frhp_addr: u64,
    fhdb_addr: u64,
) -> (Vec<u8>, Vec<u8>, Vec<[u8; 7]>) {
    // FHDB header: "FHDB"(4) + version(1) + heap_addr(8) + block_offset(4) + checksum(4) = 21 bytes
    // Wait: block_offset width is `ceil(max_heap_size_bits / 8)`.
    // Let's pick max_heap_size_bits = 32 (4 GiB) => block_offset_bytes = 4, and
    // heap_id_len = 7: 1 byte (version 0 | managed 0) + 4-byte offset + 2-byte length.
    // Wait: let's check HDF5 spec for Fractal Heap Managed Object ID:
    // byte 0: version(0) | type(0)
    // next `ceil(max_heap_size_bits / 8)` bytes: offset in heap
    // next `ceil(log2(max_direct_block_size) / 8)` or `id_len - 1 - offset_bytes` bytes: length.
    // Our reader's `parse_fractal_heap` reads:
    //   `offset_bytes = max_heap_bits.div_ceil(8)`
    //   `length_bytes = id_len - 1 - offset_bytes`
    // So with `max_heap_bits = 32` (`offset_bytes = 4`) and `id_len = 7` (`length_bytes = 2`),
    // each direct block can be up to 64 KiB (or use `max_direct_block_size = 65536` or larger if needed).
    let total_payload: usize = objects.iter().map(Vec::len).sum();
    let header_len_17 = 4 + 1 + 8 + 4; // 17 bytes
    let block_size = (header_len_17 + total_payload).next_power_of_two().max(512);
    let max_direct = block_size.max(32768).next_power_of_two();
    let mut fhdb = Vec::with_capacity(block_size);
    fhdb.extend_from_slice(b"FHDB");
    fhdb.push(0);
    fhdb.extend_from_slice(&frhp_addr.to_le_bytes());
    fhdb.extend_from_slice(&0_u32.to_le_bytes());
    let mut ids = Vec::with_capacity(objects.len());
    for obj in objects {
        let offset = fhdb.len() as u32;
        let len = obj.len() as u16;
        fhdb.extend_from_slice(obj);
        let mut id = [0_u8; 7];
        id[0] = 0;
        id[1..5].copy_from_slice(&offset.to_le_bytes());
        id[5..7].copy_from_slice(&len.to_le_bytes());
        ids.push(id);
    }
    fhdb.resize(block_size, 0);

    // Build FRHP (146 bytes with off_size=8, len_size=8)
    let mut frhp = Vec::with_capacity(146);
    frhp.extend_from_slice(b"FRHP");
    frhp.push(0); // version 0
    frhp.extend_from_slice(&7_u16.to_le_bytes()); // heap_id_len = 7
    frhp.extend_from_slice(&0_u16.to_le_bytes()); // io_filters_len = 0
    frhp.push(0x00); // flags = 0 (direct blocks not checksummed)
    frhp.extend_from_slice(&(u32::MAX).to_le_bytes()); // max_managed_obj_size (not strictly limited)
    frhp.extend_from_slice(&0_u64.to_le_bytes()); // next_huge_obj_id
    frhp.extend_from_slice(&u64::MAX.to_le_bytes()); // huge_objs_btree_addr
    frhp.extend_from_slice(&((block_size - header_len_17 - total_payload) as u64).to_le_bytes()); // free_space
    frhp.extend_from_slice(&u64::MAX.to_le_bytes()); // free_space_mgr_addr
    frhp.extend_from_slice(&(block_size as u64).to_le_bytes()); // managed_space
    frhp.extend_from_slice(&(block_size as u64).to_le_bytes()); // alloc_managed_space
    frhp.extend_from_slice(&0_u64.to_le_bytes()); // dblock_alloc_iter_offset
    frhp.extend_from_slice(&(objects.len() as u64).to_le_bytes()); // num_managed_objs
    frhp.extend_from_slice(&0_u64.to_le_bytes()); // size_huge_objs
    frhp.extend_from_slice(&0_u64.to_le_bytes()); // num_huge_objs
    frhp.extend_from_slice(&0_u64.to_le_bytes()); // size_tiny_objs
    frhp.extend_from_slice(&0_u64.to_le_bytes()); // num_tiny_objs
    frhp.extend_from_slice(&4_u16.to_le_bytes()); // table_width = 4
    frhp.extend_from_slice(&(block_size as u64).to_le_bytes()); // starting_block_size
    frhp.extend_from_slice(&(max_direct as u64).to_le_bytes()); // max_direct_block_size
    frhp.extend_from_slice(&32_u16.to_le_bytes()); // max_heap_size = 32 bits
    frhp.extend_from_slice(&1_u16.to_le_bytes()); // starting_num_rows_in_root_indirect
    frhp.extend_from_slice(&fhdb_addr.to_le_bytes()); // root_block_addr
    frhp.extend_from_slice(&0_u16.to_le_bytes()); // current_num_rows_in_root_indirect = 0 (direct!)
    let cksum = lookup3(&frhp);
    frhp.extend_from_slice(&cksum.to_le_bytes());

    (frhp, fhdb, ids)
}

/// Builds a version-2 B-tree (type 5: link name hash) with a single leaf `BTLF`.
fn build_v2_link_btree(
    entries: &[(u32, [u8; 7])], // (lookup3(name), heap_id) sorted by hash
    btlf_addr: u64,
) -> (Vec<u8>, Vec<u8>) {
    // Record size = 4 (hash) + 7 (heap_id) = 11 bytes
    let record_size = 11_usize;
    let node_size = (10 + entries.len() * record_size).max(512);

    let mut btlf = Vec::with_capacity(node_size);
    btlf.extend_from_slice(b"BTLF");
    btlf.push(0); // version 0
    btlf.push(5); // type 5 (link name index)
    for (hash, id) in entries {
        btlf.extend_from_slice(&hash.to_le_bytes());
        btlf.extend_from_slice(id);
    }
    btlf.resize(node_size - 4, 0);
    let btlf_cksum = lookup3(&btlf[..6 + entries.len() * record_size]);
    // Wait: in HDF5 BTLF, checksum is at the END of the node (`node_size - 4..node_size`),
    // computed over the entire node_size - 4 bytes? Wait: `H5B2cache.c` computes lookup3 over
    // `6 + num_records * record_size` bytes and stores it at `6 + num_records * record_size` OR at `node_size - 4`?
    // In `H5B2cache.c` (`H5B2__cache_leaf_serialize`):
    // `uint32_t metadata_chksum = H5_checksum_lookup3(image, (size_t)(p - (const uint8_t *)image));`
    // `UINT32ENCODE(p, metadata_chksum);`
    // It stores the checksum IMMEDIATELY after the used records (`6 + num_records * record_size`), NOT at `node_size - 4`!
    // Let's write it right after the records, and pad to `node_size`.
    let mut btlf_exact = Vec::with_capacity(node_size);
    btlf_exact.extend_from_slice(b"BTLF");
    btlf_exact.push(0);
    btlf_exact.push(5);
    for (hash, id) in entries {
        btlf_exact.extend_from_slice(&hash.to_le_bytes());
        btlf_exact.extend_from_slice(id);
    }
    let cksum = lookup3(&btlf_exact);
    btlf_exact.extend_from_slice(&cksum.to_le_bytes());
    btlf_exact.resize(node_size, 0);
    let _ = btlf_cksum;

    // BTHD (38 bytes)
    let mut bthd = Vec::with_capacity(38);
    bthd.extend_from_slice(b"BTHD");
    bthd.push(0); // version 0
    bthd.push(5); // type 5
    bthd.extend_from_slice(&(node_size as u32).to_le_bytes());
    bthd.extend_from_slice(&(record_size as u16).to_le_bytes());
    bthd.extend_from_slice(&0_u16.to_le_bytes()); // depth = 0
    bthd.push(100); // split_percent
    bthd.push(40); // merge_percent
    bthd.extend_from_slice(&btlf_addr.to_le_bytes());
    bthd.extend_from_slice(&(entries.len() as u16).to_le_bytes());
    bthd.extend_from_slice(&(entries.len() as u64).to_le_bytes());
    let bthd_cksum = lookup3(&bthd);
    bthd.extend_from_slice(&bthd_cksum.to_le_bytes());

    (bthd, btlf_exact)
}

fn write_dense_hdf5(
    datasets: &[Dataset],
    attributes: &[(String, Vec<Attribute>)],
    options: &H5WriteOptions,
) -> Result<Vec<u8>, String> {
    let mut prepared = Vec::with_capacity(datasets.len());
    let mut vlen_raws = Vec::with_capacity(datasets.len());
    for ds in datasets {
        let (p, v) = prepare_dataset_payload(ds, dataset_opts_for(options, &ds.path), 0)?;
        prepared.push(p);
        vlen_raws.push(v);
    }

    struct DenseGroup {
        path: String,
        children: Vec<ExtChild>,
        header_addr: u64,
        frhp_addr: u64,
        fhdb_addr: u64,
        bthd_addr: u64,
        btlf_addr: u64,
        frhp: Vec<u8>,
        fhdb: Vec<u8>,
        bthd: Vec<u8>,
        btlf: Vec<u8>,
        header_len: usize,
    }

    let mut groups: Vec<DenseGroup> = vec![DenseGroup {
        path: "/".to_owned(),
        children: Vec::new(),
        header_addr: 0,
        frhp_addr: 0,
        fhdb_addr: 0,
        bthd_addr: 0,
        btlf_addr: 0,
        frhp: Vec::new(),
        fhdb: Vec::new(),
        bthd: Vec::new(),
        btlf: Vec::new(),
        header_len: 0,
    }];
    let mut group_index: BTreeMap<String, usize> = BTreeMap::new();
    group_index.insert("/".to_owned(), 0);
    let mut dataset_paths = Vec::with_capacity(datasets.len());

    fn ensure_dense_parent(
        groups: &mut Vec<DenseGroup>,
        group_index: &mut BTreeMap<String, usize>,
        full_path: &str,
    ) -> Result<(usize, String), String> {
        let segments = path_segments(full_path)?;
        let mut parent = 0_usize;
        let mut cur = String::new();
        for segment in &segments[..segments.len() - 1] {
            cur.push('/');
            cur.push_str(segment);
            parent = match group_index.get(&cur) {
                | Some(&idx) => idx,
                | None => {
                    let idx = groups.len();
                    groups.push(DenseGroup {
                        path: cur.clone(),
                        children: Vec::new(),
                        header_addr: 0,
                        frhp_addr: 0,
                        fhdb_addr: 0,
                        bthd_addr: 0,
                        btlf_addr: 0,
                        frhp: Vec::new(),
                        fhdb: Vec::new(),
                        bthd: Vec::new(),
                        btlf: Vec::new(),
                        header_len: 0,
                    });
                    group_index.insert(cur.clone(), idx);
                    groups[parent].children.push(ExtChild::Group(idx));
                    idx
                },
            };
        }
        let full = format!("{cur}/{}", segments[segments.len() - 1]);
        Ok((parent, full))
    }

    for (pos, ds) in datasets.iter().enumerate() {
        let (parent, full) = ensure_dense_parent(&mut groups, &mut group_index, &ds.path)?;
        if group_index.contains_key(&full) {
            return Err(format!("dataset {} collides with a group", ds.path));
        }
        if dataset_paths.contains(&full) {
            return Err(format!("duplicate dataset name {full}"));
        }
        dataset_paths.push(full);
        groups[parent].children.push(ExtChild::Dataset(pos));
    }

    let mut soft_paths = Vec::with_capacity(options.soft_links.len());
    for (pos, link) in options.soft_links.iter().enumerate() {
        let (parent, full) = ensure_dense_parent(&mut groups, &mut group_index, &link.path)?;
        soft_paths.push(full);
        groups[parent].children.push(ExtChild::SoftLink(pos));
    }

    // Dataset header lengths
    let mut ds_header_lens = Vec::with_capacity(datasets.len());
    for (i, ds) in datasets.iter().enumerate() {
        let attrs = attributes_for(attributes, &ds.path);
        let mut msgs = vec![
            (0x01_u8, dataspace_message(&ds.shape)),
            (0x03_u8, dataset_datatype_msg_ext(ds)?),
            (0x05_u8, fill_value_message()),
        ];
        match &prepared[i] {
            | PreparedPayload::Contiguous { data, .. } => {
                msgs.push((0x08, layout_message(0, data.len() as u64)));
            },
            | PreparedPayload::Chunked {
                chunk_dims, filters, ..
            } => {
                if !filters.is_empty() {
                    msgs.push((0x0B, filter_pipeline_message(filters)));
                }
                msgs.push((0x08, chunked_layout_message(0, chunk_dims)));
            },
        }
        for attr in &attrs {
            msgs.push((0x0C, attribute_message(attr)?));
        }
        ds_header_lens.push(object_header_v2(&msgs)?.len());
    }

    // Group header lengths
    for g in &mut groups {
        let attrs = attributes_for(attributes, &g.path);
        let mut msgs = vec![(0x02_u8, link_info_message(0, 0))];
        for attr in &attrs {
            msgs.push((0x0C, attribute_message(attr)?));
        }
        g.header_len = object_header_v2(&msgs)?.len();
    }

    // First pass: assign dataset addresses so we know link target addresses
    let mut pos: u64 = 48; // superblock v2 is 48 bytes
    // Reserve space for groups after computing their FRHP/FHDB/BTHD/BTLF sizes
    for idx in 0..groups.len() {
        pos = align8(pos);
        groups[idx].header_addr = pos;
        pos += groups[idx].header_len as u64;

        let link_msgs: Vec<Vec<u8>> = groups[idx]
            .children
            .iter()
            .map(|c| {
                match c {
                    | ExtChild::Group(gi) => {
                        let name = groups[*gi].path.rsplit('/').next().unwrap_or("");
                        link_message_hard(name, 0)
                    },
                    | ExtChild::Dataset(di) => {
                        let name = dataset_paths[*di].rsplit('/').next().unwrap_or("");
                        link_message_hard(name, 0)
                    },
                    | ExtChild::SoftLink(si) => {
                        let name = soft_paths[*si].rsplit('/').next().unwrap_or("");
                        link_message_soft(name, &options.soft_links[*si].target)
                    },
                }
            })
            .collect();
        let names: Vec<String> = groups[idx]
            .children
            .iter()
            .map(|c| {
                match c {
                    | ExtChild::Group(gi) => {
                        groups[*gi].path.rsplit('/').next().unwrap_or("").to_owned()
                    },
                    | ExtChild::Dataset(di) => {
                        dataset_paths[*di]
                            .rsplit('/')
                            .next()
                            .unwrap_or("")
                            .to_owned()
                    },
                    | ExtChild::SoftLink(si) => {
                        soft_paths[*si].rsplit('/').next().unwrap_or("").to_owned()
                    },
                }
            })
            .collect();

        pos = align8(pos);
        let frhp_addr = pos;
        let (dummy_frhp, dummy_fhdb, dummy_ids) =
            build_fractal_heap_single_block(&link_msgs, frhp_addr, 0);
        pos += dummy_frhp.len() as u64;
        pos = align8(pos);
        let fhdb_addr = pos;
        pos += dummy_fhdb.len() as u64;

        let mut btree_entries: Vec<(u32, [u8; 7])> = names
            .iter()
            .zip(dummy_ids.iter())
            .map(|(n, id)| (lookup3(n.as_bytes()), *id))
            .collect();
        btree_entries.sort_by_key(|e| e.0);

        pos = align8(pos);
        let bthd_addr = pos;
        let (dummy_bthd, dummy_btlf) = build_v2_link_btree(&btree_entries, 0);
        pos += dummy_bthd.len() as u64;
        pos = align8(pos);
        let btlf_addr = pos;
        pos += dummy_btlf.len() as u64;

        groups[idx].frhp_addr = frhp_addr;
        groups[idx].fhdb_addr = fhdb_addr;
        groups[idx].bthd_addr = bthd_addr;
        groups[idx].btlf_addr = btlf_addr;
    }

    enum DsAddr {
        Contiguous {
            header_addr: u64,
            data_addr: u64,
            gcol_addr: Option<u64>,
        },
        Chunked {
            header_addr: u64,
            btree_addr: u64,
            chunk_addrs: Vec<u64>,
        },
    }
    let mut ds_addrs = Vec::with_capacity(datasets.len());
    for (i, payload) in prepared.iter_mut().enumerate() {
        pos = align8(pos);
        let header_addr = pos;
        pos += ds_header_lens[i] as u64;
        match payload {
            | PreparedPayload::Contiguous { data, gcol } => {
                let gcol_addr = if let Some(g) = gcol {
                    pos = align8(pos);
                    let ga = pos;
                    pos += g.len() as u64;
                    if let Some((objs, counts)) = &vlen_raws[i] {
                        let (new_gcol, new_descs) = build_gcol(objs, counts, ga);
                        *g = new_gcol;
                        *data = new_descs;
                    }
                    Some(ga)
                } else {
                    None
                };
                pos = align8(pos);
                let data_addr = pos;
                pos += data.len() as u64;
                ds_addrs.push(DsAddr::Contiguous {
                    header_addr,
                    data_addr,
                    gcol_addr,
                });
            },
            | PreparedPayload::Chunked {
                chunk_dims, chunks, ..
            } => {
                let mut chunk_addrs = Vec::with_capacity(chunks.len());
                for (_, cbytes) in chunks.iter() {
                    pos = align8(pos);
                    chunk_addrs.push(pos);
                    pos += cbytes.len() as u64;
                }
                pos = align8(pos);
                let btree_addr = pos;
                let dummy: Vec<(Vec<u64>, Vec<u8>, u64)> = chunks
                    .iter()
                    .zip(chunk_addrs.iter())
                    .map(|((o, d), &a)| (o.clone(), d.clone(), a))
                    .collect();
                pos += build_chunk_btree(chunk_dims, &dummy, btree_addr).len() as u64;
                ds_addrs.push(DsAddr::Chunked {
                    header_addr,
                    btree_addr,
                    chunk_addrs,
                });
            },
        }
    }
    let eof = align8(pos);

    // Second pass: rebuild each group's FRHP/FHDB/BTHD/BTLF with resolved addresses
    for idx in 0..groups.len() {
        let link_msgs: Vec<Vec<u8>> = groups[idx]
            .children
            .iter()
            .map(|c| {
                match c {
                    | ExtChild::Group(gi) => {
                        let name = groups[*gi].path.rsplit('/').next().unwrap_or("");
                        link_message_hard(name, groups[*gi].header_addr)
                    },
                    | ExtChild::Dataset(di) => {
                        let name = dataset_paths[*di].rsplit('/').next().unwrap_or("");
                        let addr = match &ds_addrs[*di] {
                            | DsAddr::Contiguous { header_addr, .. } => *header_addr,
                            | DsAddr::Chunked { header_addr, .. } => *header_addr,
                        };
                        link_message_hard(name, addr)
                    },
                    | ExtChild::SoftLink(si) => {
                        let name = soft_paths[*si].rsplit('/').next().unwrap_or("");
                        link_message_soft(name, &options.soft_links[*si].target)
                    },
                }
            })
            .collect();
        let names: Vec<String> = groups[idx]
            .children
            .iter()
            .map(|c| {
                match c {
                    | ExtChild::Group(gi) => {
                        groups[*gi].path.rsplit('/').next().unwrap_or("").to_owned()
                    },
                    | ExtChild::Dataset(di) => {
                        dataset_paths[*di]
                            .rsplit('/')
                            .next()
                            .unwrap_or("")
                            .to_owned()
                    },
                    | ExtChild::SoftLink(si) => {
                        soft_paths[*si].rsplit('/').next().unwrap_or("").to_owned()
                    },
                }
            })
            .collect();

        let (frhp, fhdb, ids) = build_fractal_heap_single_block(
            &link_msgs,
            groups[idx].frhp_addr,
            groups[idx].fhdb_addr,
        );
        let mut btree_entries: Vec<(u32, [u8; 7])> = names
            .iter()
            .zip(ids.iter())
            .map(|(n, id)| (lookup3(n.as_bytes()), *id))
            .collect();
        btree_entries.sort_by_key(|e| e.0);
        let (bthd, btlf) = build_v2_link_btree(&btree_entries, groups[idx].btlf_addr);

        groups[idx].frhp = frhp;
        groups[idx].fhdb = fhdb;
        groups[idx].bthd = bthd;
        groups[idx].btlf = btlf;
    }

    // Emit file bytes
    let mut out = Vec::with_capacity(eof as usize);
    out.extend_from_slice(&superblock_v2(groups[0].header_addr, eof));

    for g in &groups {
        while out.len() < g.header_addr as usize {
            out.push(0);
        }
        let attrs = attributes_for(attributes, &g.path);
        let mut msgs = vec![(0x02_u8, link_info_message(g.frhp_addr, g.bthd_addr))];
        for attr in &attrs {
            msgs.push((0x0C, attribute_message(attr)?));
        }
        out.extend_from_slice(&object_header_v2(&msgs)?);
        while out.len() < g.frhp_addr as usize {
            out.push(0);
        }
        out.extend_from_slice(&g.frhp);
        while out.len() < g.fhdb_addr as usize {
            out.push(0);
        }
        out.extend_from_slice(&g.fhdb);
        while out.len() < g.bthd_addr as usize {
            out.push(0);
        }
        out.extend_from_slice(&g.bthd);
        while out.len() < g.btlf_addr as usize {
            out.push(0);
        }
        out.extend_from_slice(&g.btlf);
    }

    for (i, ds) in datasets.iter().enumerate() {
        let attrs = attributes_for(attributes, &ds.path);
        let mut msgs = vec![
            (0x01_u8, dataspace_message(&ds.shape)),
            (0x03_u8, dataset_datatype_msg_ext(ds)?),
            (0x05_u8, fill_value_message()),
        ];
        match (&prepared[i], &ds_addrs[i]) {
            | (
                PreparedPayload::Contiguous { data, gcol },
                DsAddr::Contiguous {
                    header_addr,
                    data_addr,
                    gcol_addr,
                },
            ) => {
                msgs.push((0x08, layout_message(*data_addr, data.len() as u64)));
                for attr in &attrs {
                    msgs.push((0x0C, attribute_message(attr)?));
                }
                while out.len() < *header_addr as usize {
                    out.push(0);
                }
                out.extend_from_slice(&object_header_v2(&msgs)?);
                if let (Some(g), Some(ga)) = (gcol, gcol_addr) {
                    while out.len() < *ga as usize {
                        out.push(0);
                    }
                    out.extend_from_slice(g);
                }
                while out.len() < *data_addr as usize {
                    out.push(0);
                }
                out.extend_from_slice(data);
            },
            | (
                PreparedPayload::Chunked {
                    chunk_dims,
                    filters,
                    chunks,
                    ..
                },
                DsAddr::Chunked {
                    header_addr,
                    btree_addr,
                    chunk_addrs,
                },
            ) => {
                if !filters.is_empty() {
                    msgs.push((0x0B, filter_pipeline_message(filters)));
                }
                msgs.push((0x08, chunked_layout_message(*btree_addr, chunk_dims)));
                for attr in &attrs {
                    msgs.push((0x0C, attribute_message(attr)?));
                }
                while out.len() < *header_addr as usize {
                    out.push(0);
                }
                out.extend_from_slice(&object_header_v2(&msgs)?);
                let mut btree_items = Vec::with_capacity(chunks.len());
                for ((offsets, cbytes), &caddr) in chunks.iter().zip(chunk_addrs.iter()) {
                    while out.len() < caddr as usize {
                        out.push(0);
                    }
                    out.extend_from_slice(cbytes);
                    btree_items.push((offsets.clone(), cbytes.clone(), caddr));
                }
                while out.len() < *btree_addr as usize {
                    out.push(0);
                }
                out.extend_from_slice(&build_chunk_btree(chunk_dims, &btree_items, *btree_addr));
            },
            | _ => unreachable!(),
        }
    }

    while out.len() < eof as usize {
        out.push(0);
    }
    Ok(out)
}
