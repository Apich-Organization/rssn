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
use std::collections::HashMap;
use std::sync::Arc;
use std::sync::Mutex;

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
use crate::backend::tape::int_pow;

/// Rows per register block.
pub const LANES: usize = 128;
/// Rows below which `cpu-parallel` stops splitting.
pub const GRAIN: usize = 2048;
/// Elements per reduction chunk (fixed, so results do not depend on threads).
pub const REDUCE_CHUNK: usize = 4096;

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
    cache: Mutex<HashMap<u64, Arc<CpuKernel>>>,
}

impl CpuDevice {
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
                rayon::current_num_threads()
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
            cache: Mutex::new(HashMap::new()),
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
        let analysis = ir.analyse()?;
        let key = ir.structural_hash();
        let mut cache = self
            .cache
            .lock()
            .map_err(|_| DeviceError::Backend("kernel cache poisoned".into()))?;
        if let Some(k) = cache.get(&key).filter(|k| k.ir == *ir) {
            return Ok(k.clone());
        }
        let k = Arc::new(CpuKernel {
            ir: ir.clone(),
            analysis,
        });
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
        let lens = vec![rows; kernel.ir().inputs];
        self.launch_cols(kernel, rows, input, &lens, params, output)
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
        check_launch_cols(self.uid, ir, rows, input, lens, output)?;
        let input = self.data(input)?;
        let out = self.data_mut(output)?;
        if rows == 0 {
            return Ok(());
        }
        let cols: Vec<&mut [f64]> = out.chunks_mut(rows).collect();
        if ir.is_elementwise() && lens.iter().all(|&l| l == rows) {
            run_range(ir, input, rows, params, 0, cols, self.parallel);
        } else {
            let offsets = column_offsets(lens);
            let view = Columns {
                data: input,
                offsets: &offsets,
                lens,
            };
            run_rows(ir, analysis, view, rows, params, 0, cols, self.parallel);
        }
        Ok(())
    }

    fn reduce(
        &self,
        op: ReduceOp,
        input: &Buffer,
    ) -> Result<f64, DeviceError> {
        let data = self.data(input)?;
        let fold = |c: &[f64]| c.iter().fold(op.identity(), |a, &x| op.combine(a, x));
        let parts: Vec<f64> = if self.parallel {
            data.par_chunks(REDUCE_CHUNK).map(fold).collect()
        } else {
            data.chunks(REDUCE_CHUNK).map(fold).collect()
        };
        Ok(parts.iter().fold(op.identity(), |a, &x| op.combine(a, x)))
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
        let a = self.data(a)?;
        let b = self.data(b)?;
        let c = self.data_mut(c)?;
        if a.len() != m * k || b.len() != k * n || c.len() != m * n {
            return Err(DeviceError::Mismatch(format!(
                "matmul {m}x{k} * {k}x{n} needs buffers of {}, {}, {} elements",
                m * k,
                k * n,
                m * n
            )));
        }
        if m * n == 0 {
            return Ok(());
        }
        let row = |(i, ci): (usize, &mut [f64])| {
            ci.fill(0.0);
            for p in 0..k {
                let aip = a[i * k + p];
                let bp = &b[p * n..(p + 1) * n];
                for (cj, &bj) in ci.iter_mut().zip(bp) {
                    *cj += aip * bj;
                }
            }
        };
        if self.parallel {
            c.par_chunks_mut(n).enumerate().for_each(row);
        } else {
            c.chunks_mut(n).enumerate().for_each(row);
        }
        Ok(())
    }
}

/// Runs rows `lo..lo + len` of a structured kernel through the reference
/// evaluator; `cols[k]` is output column `k` restricted to that range.
#[allow(clippy::too_many_arguments)]
fn run_rows(
    ir: &KernelIr,
    analysis: &Analysis,
    view: Columns<'_>,
    rows: usize,
    params: &[f64],
    lo: usize,
    mut cols: Vec<&mut [f64]>,
    parallel: bool,
) {
    let len = cols.first().map_or(0, |c| c.len());
    // A folded row can cost far more than an elementwise one, so the grain
    // is smaller than the lane path's.
    if parallel && len > GRAIN / 16 {
        let mid = len / 2;
        let (mut left, mut right) = (Vec::new(), Vec::new());
        for c in cols {
            let (a, b) = c.split_at_mut(mid);
            left.push(a);
            right.push(b);
        }
        rayon::join(
            || run_rows(ir, analysis, view, rows, params, lo, left, true),
            || run_rows(ir, analysis, view, rows, params, lo + mid, right, true),
        );
        return;
    }
    let mut scratch = Scratch::default();
    let mut out = vec![f64::NAN; ir.outputs.len()];
    for at in 0..len {
        ir.eval_at(
            analysis,
            lo + at,
            rows,
            view,
            params,
            &mut scratch,
            &mut out,
        );
        for (col, &value) in cols.iter_mut().zip(&out) {
            col[at] = value;
        }
    }
}

/// Runs rows `lo..lo + len` where `cols[k]` is output column `k` restricted
/// to that range.
fn run_range(
    ir: &KernelIr,
    input: &[f64],
    rows: usize,
    params: &[f64],
    lo: usize,
    mut cols: Vec<&mut [f64]>,
    parallel: bool,
) {
    let len = cols.first().map_or(0, |c| c.len());
    if parallel && len > GRAIN {
        let mid = len / 2;
        let (mut left, mut right) = (Vec::new(), Vec::new());
        for c in cols {
            let (a, b) = c.split_at_mut(mid);
            left.push(a);
            right.push(b);
        }
        rayon::join(
            || run_range(ir, input, rows, params, lo, left, true),
            || run_range(ir, input, rows, params, lo + mid, right, true),
        );
        return;
    }
    let mut regs = vec![0.0f64; ir.ops.len() * LANES];
    let mut at = 0;
    while at < len {
        let n = LANES.min(len - at);
        exec_block(ir, input, rows, params, lo + at, n, &mut regs);
        for (col, &r) in cols.iter_mut().zip(&ir.outputs) {
            let src = &regs[r as usize * LANES..r as usize * LANES + n];
            col[at..at + n].copy_from_slice(src);
        }
        at += n;
    }
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
    for (o, &x) in tail[..n].iter_mut().zip(src) {
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
    for ((o, &p), &q) in tail[..n].iter_mut().zip(x).zip(y) {
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
            | KernelOp::PowI(a, e) => unary(regs, d, a, n, |x| int_pow(x, e)),
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
                for (((o, &k), &p), &q) in tail[..n].iter_mut().zip(cc).zip(x).zip(y) {
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
