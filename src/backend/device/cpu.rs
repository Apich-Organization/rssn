//! The reference CPU devices.
//!
//! Both devices run the same code: the kernel's registers are held as
//! *lane blocks* of [`LANES`] rows, and every instruction is one tight loop
//! over a block (`dst[l] = a[l] + b[l]`), which the compiler turns into SIMD
//! for the arithmetic instructions. `cpu-parallel` splits the rows into
//! halves recursively with `rayon::join` over disjoint slices of the output
//! columns, down to [`GRAIN`] rows; there is no `unsafe` and no per-row
//! synchronisation. Every row is computed by the same instruction sequence in
//! the same order whichever device or thread runs it, so the two devices
//! produce bit-identical results, and reductions use a fixed chunking
//! ([`REDUCE_CHUNK`]) so they do too.
//!
//! Semantics are the tape's over `f64` ([`KernelIr::eval_row`] is the
//! executable specification; the tests compare against it and against
//! [`Tape::eval`](crate::backend::tape::Tape::eval)).
//!
//! Structured kernels (index arithmetic, addressed loads, a fold; see
//! [`super::kernel`]) do not fit lane blocks: a row may read anywhere and
//! loop a data-dependent number of times. They run through the reference
//! evaluator [`KernelIr::eval_at`], one row at a time, with the same
//! recursive splitting; the lane path is kept for elementwise kernels whose
//! columns all have one element per row.

use std::any::Any;
use std::sync::Arc;
use std::sync::RwLock;

use rayon::prelude::*;

use super::Buffer;
use super::CompiledKernel;
use super::Device;
use super::DeviceError;
use super::DeviceInfo;
use super::DeviceKind;
use super::DeviceOps;
use super::Dtype;
use super::KernelIr;
use super::KernelOp;
use super::ReduceOp;
use super::Storage;
use super::check_launch;
use super::check_launch_cols;
use super::kernel::Analysis;
use super::kernel::Columns;
use super::kernel::Scratch;
use super::kernel::cmp_value;
use super::kernel::column_offsets;
use super::kernel::gather_slot;
use super::kernel::select_is_true;
use super::next_uid;
use crate::backend::Intrinsic;
use crate::backend::tape::eval_powi_slice;
use crate::graph::id::FastMap;

/// Rows per register block.
pub const LANES: usize = 128;
/// Rows below which `cpu-parallel` stops splitting.
pub const GRAIN: usize = 2048;
/// Elements per reduction chunk (fixed, so results do not depend on threads).
pub const REDUCE_CHUNK: usize = 4096;

/// Limit task overhead while keeping each wide output row divisible.
const MATMUL_TILE_COLUMNS: usize = 2_048;
const MATMUL_MIN_TILE_PRODUCTS: usize = 262_144;
/// L1-resident tile dimensions for cache-blocked `matmul`.
const MATMUL_TILE_M: usize = 4;
const MATMUL_TILE_N: usize = 256;
/// Maximum compiled kernels kept in a device's hot cache; one entry is
/// evicted when it is full.
pub(super) const MAX_CACHED_KERNELS: usize = 64;

struct HostStorage(Vec<f64>);

impl Storage for HostStorage {
    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

struct CpuKernel {
    ir: KernelIr,
    analysis: Analysis,
    /// The kernel as machine code (feature `jit`), compiled on the first
    /// launch of at least [`NATIVE_ROWS`] rows; `None` inside when it cannot
    /// be lowered (see [`super::native`]).
    #[cfg(feature = "jit")]
    native: std::sync::OnceLock<Option<super::native::NativeKernel>>,
}

/// Rows from which a launch compiles and runs the kernel natively.
#[cfg(feature = "jit")]
pub const NATIVE_ROWS: usize = 256;

/// Runs rows `lo .. lo + len(out[0])` natively over disjoint chunks; `false` if any part bailed.
#[cfg(feature = "jit")]
fn run_native(
    native: &super::native::NativeKernel,
    inputs: &[&[f64]],
    params: &[f64],
    cols: &mut [&mut [f64]],
    rows: usize,
    parallel: bool,
) -> bool {
    let len = cols.first().map_or(0, |c| c.len());
    if len == 0 {
        return true;
    }
    if !parallel || !crate::execution::should_split(len.div_ceil(GRAIN), 2) {
        return native.run(inputs, params, cols, 0, rows);
    }
    if cols.len() == 1 {
        let ok = std::sync::atomic::AtomicBool::new(true);
        crate::kernels::par::for_chunks_mut(cols[0], GRAIN, 2, |block, slice| {
            if !ok.load(std::sync::atomic::Ordering::Relaxed) || crate::cancellation::is_cancelled()
            {
                return;
            }
            let lo = block * GRAIN;
            let mut single = [slice];
            if !native.run(inputs, params, &mut single, lo, rows) {
                ok.store(false, std::sync::atomic::Ordering::Relaxed);
            }
        });
        return ok.load(std::sync::atomic::Ordering::Relaxed);
    }
    let n_chunks = len.div_ceil(GRAIN);
    let mut chunk_cols: Vec<Vec<&mut [f64]>> = (0..n_chunks)
        .map(|_| Vec::with_capacity(cols.len()))
        .collect();
    for col in cols.iter_mut() {
        let mut rest = &mut **col;
        for chunk in &mut chunk_cols {
            let take = GRAIN.min(rest.len());
            let (head, tail) = rest.split_at_mut(take);
            chunk.push(head);
            rest = tail;
        }
    }
    let ok = std::sync::atomic::AtomicBool::new(true);
    let cancellation = crate::cancellation::current();
    chunk_cols.par_iter_mut().enumerate().for_each_init(
        || crate::cancellation::install_optional(cancellation.clone()),
        |_, (chunk_idx, block_cols)| {
            if !ok.load(std::sync::atomic::Ordering::Relaxed) || crate::cancellation::is_cancelled()
            {
                return;
            }
            let lo = chunk_idx * GRAIN;
            if !native.run(inputs, params, block_cols, lo, rows) {
                ok.store(false, std::sync::atomic::Ordering::Relaxed);
            }
        },
    );
    ok.load(std::sync::atomic::Ordering::Relaxed)
}

impl CompiledKernel for CpuKernel {
    fn ir(&self) -> &KernelIr {
        &self.ir
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}

/// A CPU device, sequential or over the rayon pool.
pub struct CpuDevice {
    info: DeviceInfo,
    uid: u64,
    parallel: bool,
    cache: RwLock<FastMap<u64, Arc<CpuKernel>>>,
}

impl CpuDevice {
    /// Number of kernels currently retained in the hot cache (tests only).
    #[cfg(test)]
    pub(super) fn cached_kernel_count(&self) -> usize {
        self.cache.read().map_or(0, |cache| cache.len())
    }

    /// The single-threaded device (`cpu`).
    #[must_use]
    pub fn sequential() -> Self {
        Self::make(false)
    }

    /// The rayon device (`cpu-parallel`).
    #[must_use]
    pub fn parallel() -> Self {
        Self::make(true)
    }

    fn make(parallel: bool) -> Self {
        let info = DeviceInfo {
            id: if parallel {
                "cpu-parallel"
            } else {
                "cpu"
            }
            .into(),
            name: if parallel {
                "host CPU, rayon thread pool"
            } else {
                "host CPU, one thread"
            }
            .into(),
            kind: if parallel {
                DeviceKind::CpuParallel
            } else {
                DeviceKind::Cpu
            },
            dtypes: vec![Dtype::F64],
            ops: DeviceOps {
                elementwise: true,
                reduce: true,
                matmul: true,
                structured: true,
            },
            parallelism: if parallel {
                crate::execution::cpu_threads()
            } else {
                1
            },
            max_buffer_elems: isize::MAX as usize / 8,
            unified_memory: true,
            deterministic: true,
            emulated_f64: false,
        };
        Self {
            info,
            uid: next_uid(),
            parallel,
            cache: RwLock::new(FastMap::default()),
        }
    }

    fn data<'a>(
        &self,
        buffer: &'a Buffer,
    ) -> Result<&'a [f64], DeviceError> {
        if buffer.owner() != self.uid {
            return Err(DeviceError::Mismatch(
                "buffer belongs to another device".into(),
            ));
        }
        buffer
            .storage()
            .as_any()
            .downcast_ref::<HostStorage>()
            .map(|s| s.0.as_slice())
            .ok_or_else(|| DeviceError::Mismatch("not a host buffer".into()))
    }

    fn data_mut<'a>(
        &self,
        buffer: &'a mut Buffer,
    ) -> Result<&'a mut [f64], DeviceError> {
        if buffer.owner() != self.uid {
            return Err(DeviceError::Mismatch(
                "buffer belongs to another device".into(),
            ));
        }
        buffer
            .storage_mut()
            .as_any_mut()
            .downcast_mut::<HostStorage>()
            .map(|s| s.0.as_mut_slice())
            .ok_or_else(|| DeviceError::Mismatch("not a host buffer".into()))
    }

    fn wrap(
        &self,
        v: Vec<f64>,
    ) -> Buffer {
        let len = v.len();
        Buffer::new(self.uid, Dtype::F64, len, Box::new(HostStorage(v)))
    }

    fn need_f64(dtype: Dtype) -> Result<(), DeviceError> {
        if dtype == Dtype::F64 {
            Ok(())
        } else {
            Err(DeviceError::Unsupported(
                "the CPU devices hold f64 buffers only".into(),
            ))
        }
    }
}

impl Device for CpuDevice {
    fn info(&self) -> &DeviceInfo {
        &self.info
    }

    fn uid(&self) -> u64 {
        self.uid
    }

    fn alloc(
        &self,
        len: usize,
        dtype: Dtype,
    ) -> Result<Buffer, DeviceError> {
        Self::need_f64(dtype)?;
        let mut v = Vec::new();
        v.try_reserve_exact(len)
            .map_err(|e| DeviceError::OutOfMemory(e.to_string()))?;
        v.resize(len, 0.0);
        Ok(self.wrap(v))
    }

    fn upload(
        &self,
        host: &[f64],
        dtype: Dtype,
    ) -> Result<Buffer, DeviceError> {
        Self::need_f64(dtype)?;
        Ok(self.wrap(host.to_vec()))
    }

    fn download(
        &self,
        buffer: &Buffer,
        out: &mut [f64],
    ) -> Result<(), DeviceError> {
        let data = self.data(buffer)?;
        if data.len() != out.len() {
            return Err(DeviceError::Mismatch(format!(
                "buffer has {} elements, destination {}",
                data.len(),
                out.len()
            )));
        }
        out.copy_from_slice(data);
        Ok(())
    }

    fn copy(
        &self,
        src: &Buffer,
        src_at: usize,
        dst: &mut Buffer,
        dst_at: usize,
        len: usize,
    ) -> Result<(), DeviceError> {
        super::check_copy(self.uid, src, src_at, dst, dst_at, len)?;
        let from = self.data(src)?;
        let to = self.data_mut(dst)?;
        to[dst_at..dst_at + len].copy_from_slice(&from[src_at..src_at + len]);
        Ok(())
    }

    fn compile(
        &self,
        ir: &KernelIr,
    ) -> Result<Arc<dyn CompiledKernel>, DeviceError> {
        let key = ir.structural_hash();
        {
            let cache = self
                .cache
                .read()
                .map_err(|_| DeviceError::Backend("kernel cache poisoned".into()))?;
            if let Some(kernel) = cache.get(&key).filter(|kernel| kernel.ir.same_program(ir)) {
                return Ok(kernel.clone());
            }
        }
        let analysis = ir.analyse()?;
        let mut cache = self
            .cache
            .write()
            .map_err(|_| DeviceError::Backend("kernel cache poisoned".into()))?;
        // Another caller can publish a validated kernel while we analyse.
        if let Some(k) = cache.get(&key).filter(|k| k.ir.same_program(ir)) {
            return Ok(k.clone());
        }
        let k = Arc::new(CpuKernel {
            ir: ir.clone(),
            analysis,
            #[cfg(feature = "jit")]
            native: std::sync::OnceLock::new(),
        });
        // Retain a bounded working set; returned Arc kernels remain valid after
        // eviction. One entry is evicted rather than the whole table: clearing
        // makes a program that cycles through more than `MAX_CACHED_KERNELS`
        // distinct kernels recompile on every launch, which is far worse than
        // the eviction it is meant to bound.
        if cache.len() >= MAX_CACHED_KERNELS
            && let Some(victim) = cache.keys().next().copied()
        {
            cache.remove(&victim);
        }
        cache.insert(key, k.clone());
        drop(cache);
        Ok(k)
    }

    fn launch(
        &self,
        kernel: &dyn CompiledKernel,
        rows: usize,
        input: &Buffer,
        params: &[f64],
        output: &mut Buffer,
    ) -> Result<(), DeviceError> {
        self.launch_impl(kernel, rows, input, None, params, output)
    }

    fn launch_cols(
        &self,
        kernel: &dyn CompiledKernel,
        rows: usize,
        input: &Buffer,
        lens: &[usize],
        params: &[f64],
        output: &mut Buffer,
    ) -> Result<(), DeviceError> {
        self.launch_impl(kernel, rows, input, Some(lens), params, output)
    }

    fn reduce(
        &self,
        op: ReduceOp,
        input: &Buffer,
    ) -> Result<f64, DeviceError> {
        if crate::cancellation::is_cancelled() {
            return Err(DeviceError::Backend("computation cancelled".into()));
        }
        let data = self.data(input)?;
        let fold = |c: &[f64]| c.iter().fold(op.identity(), |a, &x| op.combine(a, x));
        let n_chunks = data.len().div_ceil(REDUCE_CHUNK);
        let result = if self.parallel && crate::execution::should_split(n_chunks, 3) {
            let mut parts = crate::kernels::reuse::f64_scratch(n_chunks);
            parts.resize(n_chunks, op.identity());
            let cancellation = crate::cancellation::current();
            data.par_chunks(REDUCE_CHUNK)
                .zip(parts.par_iter_mut())
                .for_each_init(
                    || crate::cancellation::install_optional(cancellation.clone()),
                    |_, (chunk, slot)| {
                        if !crate::cancellation::is_cancelled() {
                            *slot = fold(chunk);
                        }
                    },
                );
            if crate::cancellation::is_cancelled() {
                return Err(DeviceError::Backend("computation cancelled".into()));
            }
            parts.iter().fold(op.identity(), |a, &x| op.combine(a, x))
        } else {
            let mut acc = op.identity();
            for chunk in data.chunks(REDUCE_CHUNK) {
                if crate::cancellation::is_cancelled() {
                    return Err(DeviceError::Backend("computation cancelled".into()));
                }
                acc = op.combine(acc, fold(chunk));
            }
            acc
        };
        Ok(result)
    }

    fn matmul(
        &self,
        a: &Buffer,
        b: &Buffer,
        c: &mut Buffer,
        m: usize,
        k: usize,
        n: usize,
    ) -> Result<(), DeviceError> {
        if crate::cancellation::is_cancelled() {
            return Err(DeviceError::Backend("computation cancelled".into()));
        }
        let a = self.data(a)?;
        let b = self.data(b)?;
        let c = self.data_mut(c)?;
        let sizes = [m.checked_mul(k), k.checked_mul(n), m.checked_mul(n)];
        if sizes != [Some(a.len()), Some(b.len()), Some(c.len())] {
            return Err(DeviceError::Mismatch(format!(
                "matmul {m}x{k} * {k}x{n} is incompatible with buffers of {}, {}, {} elements",
                a.len(),
                b.len(),
                c.len()
            )));
        }
        if c.is_empty() {
            return Ok(());
        }
        let threads = crate::execution::cpu_threads();
        // A wide output can have substantial work even with few rows: tile by
        // columns while preserving the exact per-(i, j) p = 0..k accumulation.
        if self.parallel
            && m < threads.saturating_mul(2)
            && n > MATMUL_TILE_COLUMNS
            && m.saturating_mul(k).saturating_mul(n) >= MATMUL_MIN_TILE_PRODUCTS
            && crate::execution::should_split(m.saturating_mul(n.div_ceil(MATMUL_TILE_COLUMNS)), 2)
        {
            if m == 1 {
                crate::kernels::par::for_chunks_mut(c, MATMUL_TILE_COLUMNS, 2, |tile, output| {
                    if !crate::cancellation::is_cancelled() {
                        let column = tile.saturating_mul(MATMUL_TILE_COLUMNS);
                        matmul_row_tile(a, b, k, n, column, output);
                    }
                });
            } else {
                let mut tiles = Vec::new();
                for (row_idx, row_slice) in c.chunks_mut(n).enumerate() {
                    for (tile_idx, col_slice) in
                        row_slice.chunks_mut(MATMUL_TILE_COLUMNS).enumerate()
                    {
                        tiles.push((
                            row_idx,
                            tile_idx.saturating_mul(MATMUL_TILE_COLUMNS),
                            col_slice,
                        ));
                    }
                }
                let cancellation = crate::cancellation::current();
                tiles.par_iter_mut().for_each_init(
                    || crate::cancellation::install_optional(cancellation.clone()),
                    |_, (row_idx, col_start, output)| {
                        if !crate::cancellation::is_cancelled() {
                            let a_row = &a[*row_idx * k..(*row_idx + 1) * k];
                            matmul_row_tile(a_row, b, k, n, *col_start, output);
                        }
                    },
                );
            }
            if crate::cancellation::is_cancelled() {
                return Err(DeviceError::Backend("computation cancelled".into()));
            }
            return Ok(());
        }
        let min_rows = 16_384_usize.div_ceil(k.saturating_mul(n).max(1));
        if self.parallel && crate::execution::should_split(m, min_rows) {
            let band_rows = if m >= threads.saturating_mul(MATMUL_TILE_M * 2) {
                MATMUL_TILE_M
            } else {
                1
            };
            let cancellation = crate::cancellation::current();
            c.par_chunks_mut(band_rows.saturating_mul(n))
                .enumerate()
                .for_each_init(
                    || crate::cancellation::install_optional(cancellation.clone()),
                    |_, (band_idx, c_band)| {
                        if !crate::cancellation::is_cancelled() {
                            let row_start = band_idx.saturating_mul(band_rows);
                            let actual_rows = c_band.len() / n;
                            matmul_band(a, b, c_band, row_start, actual_rows, k, n);
                        }
                    },
                );
            if crate::cancellation::is_cancelled() {
                return Err(DeviceError::Backend("computation cancelled".into()));
            }
        } else {
            let mut row_start = 0;
            for c_band in c.chunks_mut(MATMUL_TILE_M.saturating_mul(n)) {
                if crate::cancellation::is_cancelled() {
                    return Err(DeviceError::Backend("computation cancelled".into()));
                }
                let actual_rows = c_band.len() / n;
                matmul_band(a, b, c_band, row_start, actual_rows, k, n);
                row_start = row_start.saturating_add(actual_rows);
            }
        }
        Ok(())
    }
}

impl CpuDevice {
    fn launch_impl(
        &self,
        kernel: &dyn CompiledKernel,
        rows: usize,
        input: &Buffer,
        lens: Option<&[usize]>,
        params: &[f64],
        output: &mut Buffer,
    ) -> Result<(), DeviceError> {
        let ir = kernel.ir();
        // A kernel compiled here carries its analysis; one that came from
        // another device is analysed again, which also validates it.
        let fresh;
        let analysis = match kernel.as_any().downcast_ref::<CpuKernel>() {
            | Some(k) => &k.analysis,
            | None => {
                fresh = ir.analyse()?;
                &fresh
            },
        };
        let elementwise = if let Some(lens) = lens {
            check_launch_cols(self.uid, ir, rows, input, lens, output)?;
            ir.is_elementwise() && lens.iter().all(|&l| l == rows)
        } else {
            check_launch(self.uid, ir, rows, input, output)?;
            ir.is_elementwise()
        };
        if crate::cancellation::is_cancelled() {
            return Err(DeviceError::Backend("computation cancelled".into()));
        }
        let input = self.data(input)?;
        let out = self.data_mut(output)?;
        if rows == 0 || ir.outputs.is_empty() {
            return Ok(());
        }
        let grain = if elementwise {
            GRAIN
        } else {
            GRAIN / 16
        };
        let parallel = self.parallel && crate::execution::should_split(rows, grain + 1);
        let mut single_out: [&mut [f64]; 1];
        let mut multi_out: Vec<&mut [f64]>;
        let cols: &mut [&mut [f64]] = if ir.outputs.len() == 1 {
            single_out = [out];
            &mut single_out
        } else {
            multi_out = out.chunks_mut(rows).collect();
            &mut multi_out
        };
        let uniform_lens_stack = [rows; 16];
        let mut uniform_lens_heap = Vec::new();
        #[cfg(feature = "jit")]
        if rows >= NATIVE_ROWS
            && let Some(k) = kernel.as_any().downcast_ref::<CpuKernel>()
            && let Some(native) = k
                .native
                .get_or_init(|| super::native::NativeKernel::compile(&k.ir, &k.analysis))
                .as_ref()
        {
            let resolved_lens: &[usize] = match lens {
                | Some(lens) => lens,
                | None if ir.inputs <= uniform_lens_stack.len() => &uniform_lens_stack[..ir.inputs],
                | None => {
                    uniform_lens_heap.resize(ir.inputs, rows);
                    &uniform_lens_heap
                },
            };
            let offsets = column_offsets(resolved_lens);
            let inputs: Vec<&[f64]> = offsets
                .iter()
                .zip(resolved_lens)
                .map(|(&at, &len)| input.get(at..at.saturating_add(len)).unwrap_or(&[]))
                .collect();
            if run_native(native, &inputs, params, cols, rows, parallel) {
                return if crate::cancellation::is_cancelled() {
                    Err(DeviceError::Backend("computation cancelled".into()))
                } else {
                    Ok(())
                };
            }
            // A bail (an index overflow): interpret the whole launch.
        }
        if elementwise {
            run_range(ir, input, rows, params, cols, parallel);
        } else {
            let resolved_lens: &[usize] = match lens {
                | Some(lens) => lens,
                | None if ir.inputs <= uniform_lens_stack.len() => &uniform_lens_stack[..ir.inputs],
                | None => {
                    uniform_lens_heap.resize(ir.inputs, rows);
                    &uniform_lens_heap
                },
            };
            let offsets = column_offsets(resolved_lens);
            let view = Columns {
                data: input,
                offsets: &offsets,
                lens: resolved_lens,
            };
            run_rows(ir, analysis, view, rows, params, cols, parallel);
        }
        if crate::cancellation::is_cancelled() {
            Err(DeviceError::Backend("computation cancelled".into()))
        } else {
            Ok(())
        }
    }
}

fn matmul_row_tile(
    a_row: &[f64],
    b: &[f64],
    k: usize,
    width: usize,
    first_column: usize,
    output: &mut [f64],
) {
    let mut sub_off = 0;
    for sub_out in output.chunks_mut(MATMUL_TILE_N) {
        sub_out.fill(0.0);
        let col_start = first_column + sub_off;
        for p in 0..k {
            let a_val = a_row[p];
            let bp = &b[p * width + col_start..p * width + col_start + sub_out.len()];
            let (out_c, out_rem) = sub_out.as_chunks_mut::<8>();
            let (bp_c, bp_rem) = bp.as_chunks::<8>();
            for (d, s) in out_c.iter_mut().zip(bp_c) {
                for l in 0..8 {
                    d[l] += a_val * s[l];
                }
            }
            for (d, &s) in out_rem.iter_mut().zip(bp_rem) {
                *d += a_val * s;
            }
        }
        sub_off += sub_out.len();
    }
}

fn matmul_band(
    a: &[f64],
    b: &[f64],
    c_band: &mut [f64],
    row_start: usize,
    band_rows: usize,
    k: usize,
    n: usize,
) {
    let mut col_start = 0;
    while col_start < n {
        let col_len = MATMUL_TILE_N.min(n - col_start);
        let mut row_off = 0;
        while row_off + MATMUL_TILE_M <= band_rows {
            let i0 = row_start + row_off;
            let i1 = i0 + 1;
            let i2 = i0 + 2;
            let i3 = i0 + 3;
            let (r0, rest) = c_band[row_off * n..].split_at_mut(n);
            let (r1, rest) = rest.split_at_mut(n);
            let (r2, rest) = rest.split_at_mut(n);
            let (r3, _) = rest.split_at_mut(n);
            let c0 = &mut r0[col_start..col_start + col_len];
            let c1 = &mut r1[col_start..col_start + col_len];
            let c2 = &mut r2[col_start..col_start + col_len];
            let c3 = &mut r3[col_start..col_start + col_len];
            c0.fill(0.0);
            c1.fill(0.0);
            c2.fill(0.0);
            c3.fill(0.0);
            for p in 0..k {
                let a0 = a[i0 * k + p];
                let a1 = a[i1 * k + p];
                let a2 = a[i2 * k + p];
                let a3 = a[i3 * k + p];
                let bp = &b[p * n + col_start..p * n + col_start + col_len];
                let (bp_c, bp_rem) = bp.as_chunks::<8>();
                let (c0_c, c0_rem) = c0.as_chunks_mut::<8>();
                let (c1_c, c1_rem) = c1.as_chunks_mut::<8>();
                let (c2_c, c2_rem) = c2.as_chunks_mut::<8>();
                let (c3_c, c3_rem) = c3.as_chunks_mut::<8>();
                for (((d0, d1), (d2, d3)), s) in c0_c
                    .iter_mut()
                    .zip(c1_c.iter_mut())
                    .zip(c2_c.iter_mut().zip(c3_c.iter_mut()))
                    .zip(bp_c)
                {
                    for l in 0..8 {
                        let bv = s[l];
                        d0[l] += a0 * bv;
                        d1[l] += a1 * bv;
                        d2[l] += a2 * bv;
                        d3[l] += a3 * bv;
                    }
                }
                for j in 0..bp_rem.len() {
                    let bv = bp_rem[j];
                    c0_rem[j] += a0 * bv;
                    c1_rem[j] += a1 * bv;
                    c2_rem[j] += a2 * bv;
                    c3_rem[j] += a3 * bv;
                }
            }
            row_off += MATMUL_TILE_M;
        }
        while row_off < band_rows {
            let i = row_start + row_off;
            let ci = &mut c_band[row_off * n + col_start..row_off * n + col_start + col_len];
            ci.fill(0.0);
            for p in 0..k {
                let aip = a[i * k + p];
                let bp = &b[p * n + col_start..p * n + col_start + col_len];
                let (ci_c, ci_rem) = ci.as_chunks_mut::<8>();
                let (bp_c, bp_rem) = bp.as_chunks::<8>();
                for (d, s) in ci_c.iter_mut().zip(bp_c) {
                    for l in 0..8 {
                        d[l] += aip * s[l];
                    }
                }
                for (d, &s) in ci_rem.iter_mut().zip(bp_rem) {
                    *d += aip * s;
                }
            }
            row_off += 1;
        }
        col_start += col_len;
    }
}

fn run_rows_chunk(
    ir: &KernelIr,
    analysis: &Analysis,
    view: Columns<'_>,
    rows: usize,
    params: &[f64],
    lo: usize,
    cols: &mut [&mut [f64]],
) {
    let len = cols.first().map_or(0, |c| c.len());
    let mut scratch = Scratch::default();
    let mut out = crate::kernels::reuse::f64_scratch(ir.outputs.len());
    out.resize(ir.outputs.len(), f64::NAN);
    for at in 0..len {
        if at.is_multiple_of(256) && crate::cancellation::is_cancelled() {
            return;
        }
        ir.eval_at(
            analysis,
            lo + at,
            rows,
            view,
            params,
            &mut scratch,
            &mut out,
        );
        for (col, &value) in cols.iter_mut().zip(out.iter()) {
            col[at] = value;
        }
    }
}

/// Runs rows `0..len` of a structured kernel through the reference evaluator;
/// `cols[k]` is output column `k`.
#[allow(clippy::too_many_arguments)]
fn run_rows(
    ir: &KernelIr,
    analysis: &Analysis,
    view: Columns<'_>,
    rows: usize,
    params: &[f64],
    cols: &mut [&mut [f64]],
    parallel: bool,
) {
    let len = cols.first().map_or(0, |c| c.len());
    // A folded row can cost far more than an elementwise one, so the grain
    // is smaller than the lane path's.
    let grain = GRAIN / 16;
    if !parallel || !crate::execution::should_split(len.div_ceil(grain), 2) {
        run_rows_chunk(ir, analysis, view, rows, params, 0, cols);
        return;
    }
    if cols.len() == 1 {
        crate::kernels::par::for_chunks_mut(cols[0], grain, 2, |block, slice| {
            if !crate::cancellation::is_cancelled() {
                run_rows_chunk(
                    ir,
                    analysis,
                    view,
                    rows,
                    params,
                    block * grain,
                    &mut [slice],
                );
            }
        });
        return;
    }
    let n_chunks = len.div_ceil(grain);
    let mut chunk_cols: Vec<Vec<&mut [f64]>> = (0..n_chunks)
        .map(|_| Vec::with_capacity(cols.len()))
        .collect();
    for col in cols.iter_mut() {
        let mut rest = &mut **col;
        for chunk in &mut chunk_cols {
            let take = grain.min(rest.len());
            let (head, tail) = rest.split_at_mut(take);
            chunk.push(head);
            rest = tail;
        }
    }
    let cancellation = crate::cancellation::current();
    chunk_cols.par_iter_mut().enumerate().for_each_init(
        || crate::cancellation::install_optional(cancellation.clone()),
        |_, (chunk_idx, block_cols)| {
            if !crate::cancellation::is_cancelled() {
                run_rows_chunk(
                    ir,
                    analysis,
                    view,
                    rows,
                    params,
                    chunk_idx * grain,
                    block_cols,
                );
            }
        },
    );
}

fn run_range_chunk(
    ir: &KernelIr,
    input: &[f64],
    rows: usize,
    params: &[f64],
    lo: usize,
    cols: &mut [&mut [f64]],
) {
    let len = cols.first().map_or(0, |c| c.len());
    let register_slots = ir.ops.len() * LANES;
    let mut regs = crate::kernels::reuse::f64_scratch(register_slots);
    regs.resize(register_slots, 0.0);
    let mut at = 0;
    while at < len {
        if crate::cancellation::is_cancelled() {
            return;
        }
        let n = LANES.min(len - at);
        exec_block(ir, input, rows, params, lo + at, n, &mut regs);
        for (col, &r) in cols.iter_mut().zip(&ir.outputs) {
            let src = &regs[r as usize * LANES..r as usize * LANES + n];
            col[at..at + n].copy_from_slice(src);
        }
        at += n;
    }
}

/// Runs rows `0..len` where `cols[k]` is output column `k`.
fn run_range(
    ir: &KernelIr,
    input: &[f64],
    rows: usize,
    params: &[f64],
    cols: &mut [&mut [f64]],
    parallel: bool,
) {
    let len = cols.first().map_or(0, |c| c.len());
    if !parallel || !crate::execution::should_split(len.div_ceil(GRAIN), 2) {
        run_range_chunk(ir, input, rows, params, 0, cols);
        return;
    }
    if cols.len() == 1 {
        crate::kernels::par::for_chunks_mut(cols[0], GRAIN, 2, |block, slice| {
            if !crate::cancellation::is_cancelled() {
                run_range_chunk(ir, input, rows, params, block * GRAIN, &mut [slice]);
            }
        });
        return;
    }
    let n_chunks = len.div_ceil(GRAIN);
    let mut chunk_cols: Vec<Vec<&mut [f64]>> = (0..n_chunks)
        .map(|_| Vec::with_capacity(cols.len()))
        .collect();
    for col in cols.iter_mut() {
        let mut rest = &mut **col;
        for chunk in &mut chunk_cols {
            let take = GRAIN.min(rest.len());
            let (head, tail) = rest.split_at_mut(take);
            chunk.push(head);
            rest = tail;
        }
    }
    let cancellation = crate::cancellation::current();
    chunk_cols.par_iter_mut().enumerate().for_each_init(
        || crate::cancellation::install_optional(cancellation.clone()),
        |_, (chunk_idx, block_cols)| {
            if !crate::cancellation::is_cancelled() {
                run_range_chunk(ir, input, rows, params, chunk_idx * GRAIN, block_cols);
            }
        },
    );
}

#[inline(always)]
fn unary(
    regs: &mut [f64],
    d: usize,
    a: u32,
    n: usize,
    f: impl Fn(f64) -> f64,
) {
    let (head, tail) = regs.split_at_mut(d * LANES);
    let src = &head[a as usize * LANES..a as usize * LANES + n];
    let dst = &mut tail[..n];
    let (dst_chunks, dst_rem) = dst.as_chunks_mut::<8>();
    let (src_chunks, src_rem) = src.as_chunks::<8>();
    for (dc, sc) in dst_chunks.iter_mut().zip(src_chunks) {
        dc[0] = f(sc[0]);
        dc[1] = f(sc[1]);
        dc[2] = f(sc[2]);
        dc[3] = f(sc[3]);
        dc[4] = f(sc[4]);
        dc[5] = f(sc[5]);
        dc[6] = f(sc[6]);
        dc[7] = f(sc[7]);
    }
    for (o, &x) in dst_rem.iter_mut().zip(src_rem) {
        *o = f(x);
    }
}

#[inline(always)]
fn binary(
    regs: &mut [f64],
    d: usize,
    a: u32,
    b: u32,
    n: usize,
    f: impl Fn(f64, f64) -> f64,
) {
    let (head, tail) = regs.split_at_mut(d * LANES);
    let x = &head[a as usize * LANES..a as usize * LANES + n];
    let y = &head[b as usize * LANES..b as usize * LANES + n];
    let dst = &mut tail[..n];
    let (dst_chunks, dst_rem) = dst.as_chunks_mut::<8>();
    let (x_chunks, x_rem) = x.as_chunks::<8>();
    let (y_chunks, y_rem) = y.as_chunks::<8>();
    for ((dc, xc), yc) in dst_chunks.iter_mut().zip(x_chunks).zip(y_chunks) {
        dc[0] = f(xc[0], yc[0]);
        dc[1] = f(xc[1], yc[1]);
        dc[2] = f(xc[2], yc[2]);
        dc[3] = f(xc[3], yc[3]);
        dc[4] = f(xc[4], yc[4]);
        dc[5] = f(xc[5], yc[5]);
        dc[6] = f(xc[6], yc[6]);
        dc[7] = f(xc[7], yc[7]);
    }
    for ((o, &p), &q) in dst_rem.iter_mut().zip(x_rem).zip(y_rem) {
        *o = f(p, q);
    }
}

/// Executes the program for rows `lo..lo + n` (`n <= LANES`).
fn exec_block(
    ir: &KernelIr,
    input: &[f64],
    rows: usize,
    params: &[f64],
    lo: usize,
    n: usize,
    regs: &mut [f64],
) {
    for (d, op) in ir.ops.iter().enumerate() {
        match *op {
            | KernelOp::Const(c) => regs[d * LANES..d * LANES + n].fill(c),
            | KernelOp::Input(i) => {
                let base = i as usize * rows + lo;
                regs[d * LANES..d * LANES + n].copy_from_slice(&input[base..base + n]);
            },
            | KernelOp::Param(i) => {
                let v = params.get(i as usize).copied().unwrap_or(f64::NAN);
                regs[d * LANES..d * LANES + n].fill(v);
            },
            | KernelOp::Add(a, b) => binary(regs, d, a, b, n, |x, y| x + y),
            | KernelOp::Sub(a, b) => binary(regs, d, a, b, n, |x, y| x - y),
            | KernelOp::Mul(a, b) => binary(regs, d, a, b, n, |x, y| x * y),
            | KernelOp::Div(a, b) => binary(regs, d, a, b, n, |x, y| x / y),
            | KernelOp::Neg(a) => unary(regs, d, a, n, |x| -x),
            | KernelOp::PowI(a, e) => {
                let (head, tail) = regs.split_at_mut(d * LANES);
                let start = a as usize * LANES;
                eval_powi_slice(&mut tail[..n], &head[start..start + n], e);
            },
            | KernelOp::PowF(a, b) => binary(regs, d, a, b, n, f64::powf),
            | KernelOp::Atan2(a, b) => binary(regs, d, a, b, n, f64::atan2),
            | KernelOp::Unary(f, a) => {
                match f {
                    | Intrinsic::Sqrt => unary(regs, d, a, n, f64::sqrt),
                    | Intrinsic::Abs => unary(regs, d, a, n, f64::abs),
                    | Intrinsic::Floor => unary(regs, d, a, n, f64::floor),
                    | Intrinsic::Ceil => unary(regs, d, a, n, f64::ceil),
                    | other => unary(regs, d, a, n, |x| other.apply(x)),
                }
            },
            | KernelOp::Cmp(c, a, b) => binary(regs, d, a, b, n, |x, y| cmp_value(c, x, y)),
            | KernelOp::Select(c, a, b) => {
                let (head, tail) = regs.split_at_mut(d * LANES);
                let cc = &head[c as usize * LANES..c as usize * LANES + n];
                let x = &head[a as usize * LANES..a as usize * LANES + n];
                let y = &head[b as usize * LANES..b as usize * LANES + n];
                let dst = &mut tail[..n];
                let (dst_chunks, dst_rem) = dst.as_chunks_mut::<8>();
                let (cc_chunks, cc_rem) = cc.as_chunks::<8>();
                let (x_chunks, x_rem) = x.as_chunks::<8>();
                let (y_chunks, y_rem) = y.as_chunks::<8>();
                for (((dc, kc), xc), yc) in dst_chunks
                    .iter_mut()
                    .zip(cc_chunks)
                    .zip(x_chunks)
                    .zip(y_chunks)
                {
                    for l in 0..8 {
                        dc[l] = if select_is_true(kc[l]) {
                            xc[l]
                        } else {
                            yc[l]
                        };
                    }
                }
                for (((o, &k), &p), &q) in dst_rem.iter_mut().zip(cc_rem).zip(x_rem).zip(y_rem) {
                    *o = if select_is_true(k) { p } else { q };
                }
            },
            | KernelOp::Gather(base, idx) => {
                let (head, tail) = regs.split_at_mut(d * LANES);
                let offs = &head[idx as usize * LANES..idx as usize * LANES + n];
                for (l, (o, &off)) in tail[..n].iter_mut().zip(offs).enumerate() {
                    *o = gather_slot(base, off, ir.inputs)
                        .map_or(f64::NAN, |slot| input[slot * rows + lo + l]);
                }
            },
            // Structured instructions never reach the lane path
            // (`launch_cols` routes them to `run_rows`).
            | _ => regs[d * LANES..d * LANES + n].fill(f64::NAN),
        }
    }
}
