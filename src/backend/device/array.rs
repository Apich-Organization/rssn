//! Array programs: everything beyond an elementwise map, as kernels.
//!
//! A device implements one entry point,
//! [`Device::launch_cols`](super::Device::launch_cols). This module shows
//! why that is enough: each classic data-parallel pattern is a function that
//! *builds a [`KernelIr`]*, and the drivers below launch those kernels through
//! the trait and nothing else. A device that runs structured kernels
//! therefore gets all of it without writing a line of shader for any of it.
//!
//! | pattern | kernel | rows | what the fold does |
//! |---|---|---|---|
//! | reduction | [`reduce_pass`] | `ceil(n / FAN_IN)` | combines one block of the column |
//! | map, then reduce (fused) | [`map_reduce_pass`] | `ceil(n / FAN_IN)` | evaluates the map inside the block; the mapped column is never stored |
//! | matrix product | [`matmul_kernel`] | `m * n` | the dot product of a row and a column |
//! | sparse matrix times vector | [`spmv_csr_kernel`] | matrix rows | the stored entries of one row |
//! | stencil, 1-D and 2-D | [`stencil1d_kernel`], [`stencil2d_kernel`] | grid points | no fold: a fixed set of neighbour loads |
//! | all-pairs interaction | [`pairwise_kernel`] | bodies | every partner of one body |
//! | adaptive ODE ensemble | [`dopri5_step`] | trajectories | no fold: one embedded Runge-Kutta attempt per row with its own step size |
//!
//! Loop-invariant work is hoisted for free: the IR assigns each register a
//! phase from its data flow, so in the all-pairs kernel the loads of "this
//! body" run once per row although the builder never says so.
//!
//! ## Determinism
//!
//! [`reduce`] always combines in the same order: blocks of [`FAN_IN`]
//! consecutive elements, level by level. Two devices with the same float
//! type therefore agree bit for bit. (A device's own
//! [`Device::reduce`](super::Device::reduce) may use another order; it only
//! promises to be deterministic for that device.)

use super::Buffer;
use super::Device;
use super::DeviceError;
use super::Fold;
use super::KernelBuilder;
use super::KernelIr;
use super::KernelOp;
use super::ReduceOp;
use super::kernel::rhs_shape;
use super::preferred_dtype;
use crate::backend::Intrinsic;
use crate::backend::tape::Cmp;
use crate::backend::tape::Reg;

/// Elements combined by one row of a reduction pass.
pub const FAN_IN: usize = 256;

fn col(i: usize) -> u32 {
    u32::try_from(i).unwrap_or(u32::MAX)
}

/// How a stencil treats neighbours outside the grid.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Boundary {
    /// Outside values are zero.
    Zero,
    /// The nearest grid value is repeated.
    Clamp,
    /// The grid wraps around.
    Periodic,
}

/// One pass of a tree reduction over column `column` of a buffer with
/// `inputs` columns: row `r` combines elements `r * FAN_IN ..` of that
/// column.
///
/// # Errors
/// [`DeviceError::InvalidKernel`] when `column >= inputs`.
pub fn reduce_pass(
    op: ReduceOp,
    inputs: usize,
    column: usize,
) -> Result<KernelIr, DeviceError> {
    let mut b = KernelBuilder::new(inputs, 0);
    let (at, ok, fan) = block_index(&mut b, column);
    let v = b.push(KernelOp::Load(col(column), at));
    let ident = b.num(op.identity());
    let x = b.push(KernelOp::Select(ok, v, ident));
    let acc = b.push(KernelOp::Acc(0));
    b.finish_fold(
        vec![acc],
        Some(Fold {
            count: fan,
            accs: vec![(op, x)],
        }),
    )
}

/// `(index, in_range, fan)` of the element a reduction row visits at the
/// inner index.
fn block_index(
    b: &mut KernelBuilder,
    column: usize,
) -> (Reg, Reg, Reg) {
    let fan = b.idx(i64::try_from(FAN_IN).unwrap_or(256));
    let row = b.push(KernelOp::Row);
    let base = b.push(KernelOp::IMul(row, fan));
    let inner = b.push(KernelOp::Inner);
    let at = b.push(KernelOp::IAdd(base, inner));
    let len = b.push(KernelOp::Len(col(column)));
    let ok = b.push(KernelOp::ICmp(Cmp::Lt, at, len));
    (at, ok, fan)
}

/// The first pass of a fused map-reduce: `map` (elementwise, any number of
/// inputs and outputs) is evaluated inside the reduction block and each of
/// its outputs is combined with `op`. The mapped columns never exist.
///
/// # Errors
/// [`DeviceError::InvalidKernel`] when `map` is not elementwise or uses
/// `gather`.
pub fn map_reduce_pass(
    map: &KernelIr,
    op: ReduceOp,
) -> Result<KernelIr, DeviceError> {
    if !map.is_elementwise() || map.ops.iter().any(|o| matches!(o, KernelOp::Gather(..))) {
        return Err(DeviceError::InvalidKernel(
            "a fused map must be elementwise and cannot use gather".into(),
        ));
    }
    if map.inputs == 0 {
        return Err(DeviceError::InvalidKernel(
            "a fused map needs at least one input column".into(),
        ));
    }
    let mut b = KernelBuilder::new(map.inputs, map.params);
    let (at, ok, fan) = block_index(&mut b, 0);
    let inputs: Vec<Reg> = (0..map.inputs)
        .map(|i| b.push(KernelOp::Load(col(i), at)))
        .collect();
    let params: Vec<Reg> = (0..map.params)
        .map(|j| b.push(KernelOp::Param(col(j))))
        .collect();
    let values = b.splice(map, &inputs, &params);
    let ident = b.num(op.identity());
    let mut accs = Vec::with_capacity(values.len());
    let mut outputs = Vec::with_capacity(values.len());
    for (k, &v) in values.iter().enumerate() {
        let x = b.push(KernelOp::Select(ok, v, ident));
        accs.push((op, x));
        outputs.push(b.push(KernelOp::Acc(col(k))));
    }
    b.finish_fold(outputs, Some(Fold { count: fan, accs }))
}

/// `C = A * B` for row-major `A` (`m x k`, column 0) and `B` (`k x n`,
/// column 1), one row of the launch per element of `C` (`m * n` rows).
/// Parameter 0 is `n`; `k` follows from the length of `B`.
///
/// # Errors
/// Never in practice; the kernel is fixed.
pub fn matmul_kernel() -> Result<KernelIr, DeviceError> {
    let mut b = KernelBuilder::new(2, 1);
    let n_num = b.push(KernelOp::Param(0));
    let n = b.push(KernelOp::ToIdx(n_num));
    let len_b = b.push(KernelOp::Len(1));
    let k = b.push(KernelOp::IDiv(len_b, n));
    let row = b.push(KernelOp::Row);
    let i = b.push(KernelOp::IDiv(row, n));
    let c = b.push(KernelOp::IMod(row, n));
    let a_base = b.push(KernelOp::IMul(i, k));
    let inner = b.push(KernelOp::Inner);
    let a_at = b.push(KernelOp::IAdd(a_base, inner));
    let b_row = b.push(KernelOp::IMul(inner, n));
    let b_at = b.push(KernelOp::IAdd(b_row, c));
    let x = b.push(KernelOp::Load(0, a_at));
    let y = b.push(KernelOp::Load(1, b_at));
    let xy = b.push(KernelOp::Mul(x, y));
    let acc = b.push(KernelOp::Acc(0));
    b.finish_fold(
        vec![acc],
        Some(Fold {
            count: k,
            accs: vec![(ReduceOp::Sum, xy)],
        }),
    )
}

/// `y = A x` for a CSR matrix: column 0 the row pointers (`rows + 1`
/// entries), column 1 the column indices, column 2 the values, column 3 the
/// vector. One row of the launch per matrix row.
///
/// Row pointers and column indices are stored as floats, so they are exact
/// below `2^24` on an `f32` device (`2^53` on an `f64` one).
///
/// # Errors
/// Never in practice; the kernel is fixed.
pub fn spmv_csr_kernel() -> Result<KernelIr, DeviceError> {
    let mut b = KernelBuilder::new(4, 0);
    let row = b.push(KernelOp::Row);
    let one = b.idx(1);
    let next = b.push(KernelOp::IAdd(row, one));
    let start_num = b.push(KernelOp::Load(0, row));
    let end_num = b.push(KernelOp::Load(0, next));
    let start = b.push(KernelOp::ToIdx(start_num));
    let end = b.push(KernelOp::ToIdx(end_num));
    let count = b.push(KernelOp::ISub(end, start));
    let inner = b.push(KernelOp::Inner);
    let at = b.push(KernelOp::IAdd(start, inner));
    let value = b.push(KernelOp::Load(2, at));
    let column_num = b.push(KernelOp::Load(1, at));
    let column = b.push(KernelOp::ToIdx(column_num));
    let x = b.push(KernelOp::Load(3, column));
    let term = b.push(KernelOp::Mul(value, x));
    let acc = b.push(KernelOp::Acc(0));
    b.finish_fold(
        vec![acc],
        Some(Fold {
            count,
            accs: vec![(ReduceOp::Sum, term)],
        }),
    )
}

/// Resolves the index `i` on an axis of extent `extent` under `boundary`:
/// the index to load and, for [`Boundary::Zero`], whether it is inside.
fn resolve(
    b: &mut KernelBuilder,
    i: Reg,
    extent: Reg,
    boundary: Boundary,
) -> (Reg, Option<Reg>) {
    let zero = b.idx(0);
    match boundary {
        | Boundary::Zero => {
            let lo = b.push(KernelOp::ICmp(Cmp::Ge, i, zero));
            let hi = b.push(KernelOp::ICmp(Cmp::Lt, i, extent));
            let inside = b.push(KernelOp::Mul(lo, hi));
            (i, Some(inside))
        },
        | Boundary::Clamp => {
            let one = b.idx(1);
            let last = b.push(KernelOp::ISub(extent, one));
            let low = b.push(KernelOp::IMax(i, zero));
            (b.push(KernelOp::IMin(low, last)), None)
        },
        | Boundary::Periodic => {
            let m = b.push(KernelOp::IMod(i, extent));
            let shifted = b.push(KernelOp::IAdd(m, extent));
            (b.push(KernelOp::IMod(shifted, extent)), None)
        },
    }
}

/// Adds `weight * value` (zero outside the grid when `inside` is given) to
/// the running sum.
fn accumulate(
    b: &mut KernelBuilder,
    sum: Option<Reg>,
    weight: f64,
    value: Reg,
    inside: Option<Reg>,
) -> Reg {
    let value = match inside {
        | Some(flag) => {
            let zero = b.num(0.0);
            b.push(KernelOp::Select(flag, value, zero))
        },
        | None => value,
    };
    let w = b.num(weight);
    let term = b.push(KernelOp::Mul(w, value));
    match sum {
        | Some(s) => b.push(KernelOp::Add(s, term)),
        | None => term,
    }
}

/// A 1-D stencil: `out[r] = sum of weight * x[r + offset]` over `points`.
///
/// # Errors
/// [`DeviceError::InvalidKernel`] when `points` is empty or an offset does
/// not fit 32 bits.
pub fn stencil1d_kernel(
    points: &[(i64, f64)],
    boundary: Boundary,
) -> Result<KernelIr, DeviceError> {
    let mut b = KernelBuilder::new(1, 0);
    let row = b.push(KernelOp::Row);
    let len = b.push(KernelOp::Len(0));
    let mut sum = None;
    for &(offset, weight) in points {
        let d = b.idx(offset);
        let i = b.push(KernelOp::IAdd(row, d));
        let (at, inside) = resolve(&mut b, i, len, boundary);
        let v = b.push(KernelOp::Load(0, at));
        sum = Some(accumulate(&mut b, sum, weight, v, inside));
    }
    let out = sum.ok_or_else(|| DeviceError::InvalidKernel("a stencil needs points".into()))?;
    b.finish(vec![out])
}

/// A 2-D stencil on a row-major `ny x nx` grid stored in one column:
/// `out[y][x] = sum of weight * u[y + dy][x + dx]` over `points`
/// (`(dx, dy, weight)`). Parameter 0 is `nx`; `ny` follows from the column
/// length. The boundary rule applies per axis.
///
/// # Errors
/// [`DeviceError::InvalidKernel`] when `points` is empty or an offset does
/// not fit 32 bits.
pub fn stencil2d_kernel(
    points: &[(i64, i64, f64)],
    boundary: Boundary,
) -> Result<KernelIr, DeviceError> {
    let mut b = KernelBuilder::new(1, 1);
    let nx_num = b.push(KernelOp::Param(0));
    let nx = b.push(KernelOp::ToIdx(nx_num));
    let len = b.push(KernelOp::Len(0));
    let ny = b.push(KernelOp::IDiv(len, nx));
    let row = b.push(KernelOp::Row);
    let x = b.push(KernelOp::IMod(row, nx));
    let y = b.push(KernelOp::IDiv(row, nx));
    let mut sum = None;
    for &(dx, dy, weight) in points {
        let ox = b.idx(dx);
        let oy = b.idx(dy);
        let xi = b.push(KernelOp::IAdd(x, ox));
        let yi = b.push(KernelOp::IAdd(y, oy));
        let (xx, in_x) = resolve(&mut b, xi, nx, boundary);
        let (yy, in_y) = resolve(&mut b, yi, ny, boundary);
        let line = b.push(KernelOp::IMul(yy, nx));
        let at = b.push(KernelOp::IAdd(line, xx));
        let v = b.push(KernelOp::Load(0, at));
        let inside = match (in_x, in_y) {
            | (Some(p), Some(q)) => Some(b.push(KernelOp::Mul(p, q))),
            | (p, q) => p.or(q),
        };
        sum = Some(accumulate(&mut b, sum, weight, v, inside));
    }
    let out = sum.ok_or_else(|| DeviceError::InvalidKernel("a stencil needs points".into()))?;
    b.finish(vec![out])
}

/// An all-pairs interaction over `columns` attribute columns.
///
/// `pair` is elementwise with `2 * columns` inputs, the attributes of body
/// `i` followed by those of body `j`, and any number of outputs; output `k`
/// of the result is the `op`-combination over every `j` of output `k` of
/// `pair`. With `skip_self` the term `j == i` is left out. The attributes of
/// body `i` are loaded once per row, not once per pair.
///
/// # Errors
/// [`DeviceError::InvalidKernel`] when `pair` does not have that shape.
pub fn pairwise_kernel(
    pair: &KernelIr,
    columns: usize,
    op: ReduceOp,
    skip_self: bool,
) -> Result<KernelIr, DeviceError> {
    if pair.inputs != columns * 2 || columns == 0 || !pair.is_elementwise() {
        return Err(DeviceError::InvalidKernel(format!(
            "a pair kernel over {columns} columns must be elementwise with {} inputs",
            columns * 2
        )));
    }
    if pair.ops.iter().any(|o| matches!(o, KernelOp::Gather(..))) {
        return Err(DeviceError::InvalidKernel(
            "a pair kernel cannot use gather".into(),
        ));
    }
    let mut b = KernelBuilder::new(columns, pair.params);
    let row = b.push(KernelOp::Row);
    let inner = b.push(KernelOp::Inner);
    let count = b.push(KernelOp::Rows);
    let mut inputs: Vec<Reg> = (0..columns)
        .map(|c| b.push(KernelOp::Load(col(c), row)))
        .collect();
    inputs.extend((0..columns).map(|c| b.push(KernelOp::Load(col(c), inner))));
    let params: Vec<Reg> = (0..pair.params)
        .map(|j| b.push(KernelOp::Param(col(j))))
        .collect();
    let values = b.splice(pair, &inputs, &params);
    let other = b.push(KernelOp::ICmp(Cmp::Ne, row, inner));
    let ident = b.num(op.identity());
    let mut accs = Vec::with_capacity(values.len());
    let mut outputs = Vec::with_capacity(values.len());
    for (k, &v) in values.iter().enumerate() {
        let x = if skip_self {
            b.push(KernelOp::Select(other, v, ident))
        } else {
            v
        };
        accs.push((op, x));
        outputs.push(b.push(KernelOp::Acc(col(k))));
    }
    b.finish_fold(outputs, Some(Fold { count, accs }))
}

// --- drivers ------------------------------------------------------------------

/// Launches `ir` over host columns of independent lengths and downloads the
/// outputs, one vector per output column.
///
/// # Errors
/// Any [`DeviceError`].
pub fn run_cols(
    device: &dyn Device,
    ir: &KernelIr,
    rows: usize,
    columns: &[&[f64]],
    params: &[f64],
) -> Result<Vec<Vec<f64>>, DeviceError> {
    let lens: Vec<usize> = columns.iter().map(|c| c.len()).collect();
    let mut flat = Vec::with_capacity(lens.iter().sum());
    for c in columns {
        flat.extend_from_slice(c);
    }
    let dtype = preferred_dtype(device.info());
    let kernel = device.compile(ir)?;
    let input = device.upload(&flat, dtype)?;
    let mut output = device.alloc(ir.outputs.len().saturating_mul(rows), dtype)?;
    device.launch_cols(&*kernel, rows, &input, &lens, params, &mut output)?;
    let mut host = vec![0.0; output.len()];
    device.download(&output, &mut host)?;
    Ok(if rows == 0 {
        vec![Vec::new(); ir.outputs.len()]
    } else {
        host.chunks(rows).map(<[f64]>::to_vec).collect()
    })
}

/// Reduces column `column` of `buffer`, whose columns have the lengths
/// `lens`, to one number, on the device.
///
/// Only the final value crosses to the host. An empty column gives the
/// identity of `op`.
///
/// # Errors
/// Any [`DeviceError`]; [`DeviceError::Unsupported`] on a device without
/// structured kernels.
pub fn reduce_column(
    device: &dyn Device,
    op: ReduceOp,
    buffer: &Buffer,
    lens: &[usize],
    column: usize,
) -> Result<f64, DeviceError> {
    let Some(&n) = lens.get(column) else {
        return Err(DeviceError::Mismatch(format!(
            "column {column} of {} columns",
            lens.len()
        )));
    };
    if n == 0 {
        return Ok(op.identity());
    }
    let first = device.compile(&reduce_pass(op, lens.len(), column)?)?;
    let rows = n.div_ceil(FAN_IN);
    let mut level = device.alloc(rows, buffer.dtype())?;
    device.launch_cols(&*first, rows, buffer, lens, &[], &mut level)?;
    finish_reduction(device, op, level, rows)
}

/// Continues a tree reduction from a one-column buffer of `n` partial
/// results down to one value and downloads it.
fn finish_reduction(
    device: &dyn Device,
    op: ReduceOp,
    mut level: Buffer,
    mut n: usize,
) -> Result<f64, DeviceError> {
    let pass = device.compile(&reduce_pass(op, 1, 0)?)?;
    while n > 1 {
        let rows = n.div_ceil(FAN_IN);
        let mut next = device.alloc(rows, level.dtype())?;
        device.launch_cols(&*pass, rows, &level, &[n], &[], &mut next)?;
        level = next;
        n = rows;
    }
    let mut host = [0.0];
    device.download(&level, &mut host)?;
    Ok(host[0])
}

/// Reduces a whole one-column buffer on the device (see [`reduce_column`]).
///
/// # Errors
/// As [`reduce_column`].
pub fn reduce(
    device: &dyn Device,
    op: ReduceOp,
    buffer: &Buffer,
) -> Result<f64, DeviceError> {
    reduce_column(device, op, buffer, &[buffer.len()], 0)
}

/// `op` over every output of `map` applied to `columns`, fused: one value
/// per output of `map`, and no mapped column is ever stored.
///
/// # Errors
/// Any [`DeviceError`]; [`DeviceError::Mismatch`] when the columns differ in
/// length or number.
pub fn map_reduce(
    device: &dyn Device,
    map: &KernelIr,
    columns: &[&[f64]],
    params: &[f64],
    op: ReduceOp,
) -> Result<Vec<f64>, DeviceError> {
    if columns.len() != map.inputs {
        return Err(DeviceError::Mismatch(format!(
            "the map has {} inputs, {} columns given",
            map.inputs,
            columns.len()
        )));
    }
    let n = columns.first().map_or(0, |c| c.len());
    if columns.iter().any(|c| c.len() != n) {
        return Err(DeviceError::Mismatch("columns differ in length".into()));
    }
    let outputs = map.outputs.len();
    if n == 0 {
        return Ok(vec![op.identity(); outputs]);
    }
    let fused = map_reduce_pass(map, op)?;
    let lens = vec![n; columns.len()];
    let mut flat = Vec::with_capacity(n.saturating_mul(columns.len()));
    for c in columns {
        flat.extend_from_slice(c);
    }
    let dtype = preferred_dtype(device.info());
    let kernel = device.compile(&fused)?;
    let input = device.upload(&flat, dtype)?;
    let rows = n.div_ceil(FAN_IN);
    let mut level = device.alloc(outputs.saturating_mul(rows), dtype)?;
    device.launch_cols(&*kernel, rows, &input, &lens, params, &mut level)?;
    let level_lens = vec![rows; outputs];
    (0..outputs)
        .map(|k| reduce_column(device, op, &level, &level_lens, k))
        .collect()
}

/// `c = a * b` for row-major `a` (`m x k`), `b` (`k x n`), `c` (`m x n`),
/// through [`matmul_kernel`]. This is what a device without a tuned
/// `matmul` of its own uses.
///
/// # Errors
/// [`DeviceError::Mismatch`] on wrong sizes; any error of the launch.
pub fn matmul(
    device: &dyn Device,
    a: &Buffer,
    b: &Buffer,
    c: &mut Buffer,
    m: usize,
    k: usize,
    n: usize,
) -> Result<(), DeviceError> {
    let sizes = (m.checked_mul(k), k.checked_mul(n), m.checked_mul(n));
    if sizes != (Some(a.len()), Some(b.len()), Some(c.len())) {
        return Err(DeviceError::Mismatch(format!(
            "matmul {m}x{k} * {k}x{n} needs buffers of {:?}, {:?}, {:?} elements",
            sizes.0, sizes.1, sizes.2
        )));
    }
    if c.is_empty() {
        return Ok(());
    }
    if k == 0 {
        // An empty sum: the product is the zero matrix.
        let zeros = device.upload(&vec![0.0; c.len()], c.dtype())?;
        *c = zeros;
        return Ok(());
    }
    // A kernel reads its columns from one buffer: join the operands there
    // with device-side copies (nothing crosses to the host on a device that
    // implements `copy`).
    let mut joined = device.alloc(a.len() + b.len(), a.dtype())?;
    device.copy(a, 0, &mut joined, 0, a.len())?;
    device.copy(b, 0, &mut joined, a.len(), b.len())?;
    let kernel = device.compile(&matmul_kernel()?)?;
    #[allow(clippy::cast_precision_loss)]
    let params = [n as f64];
    device.launch_cols(&*kernel, m * n, &joined, &[m * k, k * n], &params, c)
}

/// Step-size control of [`ensemble_dopri5`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Tolerance {
    /// Relative tolerance.
    pub rtol: f64,
    /// Absolute tolerance.
    pub atol: f64,
    /// First step; `None` chooses `(t1 - t0) / 64`.
    pub first_step: Option<f64>,
    /// Launches after which the integration gives up.
    pub max_launches: usize,
}

impl Default for Tolerance {
    fn default() -> Self {
        Self {
            rtol: 1.0e-6,
            atol: 1.0e-9,
            first_step: None,
            max_launches: 100_000,
        }
    }
}

/// What [`ensemble_dopri5`] returns.
#[derive(Clone, Debug, PartialEq)]
pub struct Ensemble {
    /// `y[i][row]`: state variable `i` of trajectory `row` at `t1`.
    pub y: Vec<Vec<f64>>,
    /// Number of step launches (each advances every unfinished trajectory by
    /// one attempted step).
    pub launches: usize,
}

/// One attempted Dormand-Prince 5(4) step per row, with per-row time and
/// step size.
///
/// Inputs: the `n_state` state columns, then `t`, then `h`. Parameters: those
/// of `rhs`, then `t_end` and `t_done` (the time from which a row counts as
/// finished), then `rtol` and `atol`. Outputs: the same `n_state + 2`
/// columns after the attempt. A row whose error estimate is within tolerance
/// advances; any other keeps its state and retries with a smaller step; a
/// finished row is left untouched. No row waits for another.
///
/// # Errors
/// [`DeviceError::InvalidKernel`] when `rhs` is not an ODE right-hand side
/// for `n_state` variables (see [`super::rk4_step`]).
#[allow(clippy::too_many_lines, clippy::many_single_char_names)]
pub fn dopri5_step(
    rhs: &KernelIr,
    n_state: usize,
) -> Result<KernelIr, DeviceError> {
    // The Dormand-Prince tableau: stage coefficients, the fifth-order
    // weights and the difference to the embedded fourth-order weights.
    const A: [&[f64]; 6] = [
        &[1.0 / 5.0],
        &[3.0 / 40.0, 9.0 / 40.0],
        &[44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0],
        &[
            19372.0 / 6561.0,
            -25360.0 / 2187.0,
            64448.0 / 6561.0,
            -212.0 / 729.0,
        ],
        &[
            9017.0 / 3168.0,
            -355.0 / 33.0,
            46732.0 / 5247.0,
            49.0 / 176.0,
            -5103.0 / 18656.0,
        ],
        &[
            35.0 / 384.0,
            0.0,
            500.0 / 1113.0,
            125.0 / 192.0,
            -2187.0 / 6784.0,
            11.0 / 84.0,
        ],
    ];
    const C: [f64; 6] = [1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0];
    const E: [f64; 7] = [
        35.0 / 384.0 - 5179.0 / 57600.0,
        0.0,
        500.0 / 1113.0 - 7571.0 / 16695.0,
        125.0 / 192.0 - 393.0 / 640.0,
        -2187.0 / 6784.0 + 92097.0 / 339_200.0,
        11.0 / 84.0 - 187.0 / 2100.0,
        -1.0 / 40.0,
    ];
    let timed = rhs_shape(rhs, n_state)?;
    let p = rhs.params;
    let mut b = KernelBuilder::new(n_state + 2, p + 4);
    let params: Vec<Reg> = (0..p).map(|j| b.push(KernelOp::Param(col(j)))).collect();
    let t_end = b.push(KernelOp::Param(col(p)));
    let t_done = b.push(KernelOp::Param(col(p + 1)));
    let rtol = b.push(KernelOp::Param(col(p + 2)));
    let atol = b.push(KernelOp::Param(col(p + 3)));
    let y: Vec<Reg> = (0..n_state)
        .map(|i| b.push(KernelOp::Input(col(i))))
        .collect();
    let t = b.push(KernelOp::Input(col(n_state)));
    let h_in = b.push(KernelOp::Input(col(n_state + 1)));
    let active = b.push(KernelOp::Cmp(Cmp::Lt, t, t_done));
    // Never step past the end.
    let left = b.push(KernelOp::Sub(t_end, t));
    let shorter = b.push(KernelOp::Cmp(Cmp::Lt, left, h_in));
    let h = b.push(KernelOp::Select(shorter, left, h_in));

    let eval = |b: &mut KernelBuilder, time: Reg, state: &[Reg]| -> Vec<Reg> {
        let mut ins = Vec::with_capacity(rhs.inputs);
        if timed {
            ins.push(time);
        }
        ins.extend_from_slice(state);
        b.splice(rhs, &ins, &params)
    };
    // y + h * sum(coefficient * k).
    let combine = |b: &mut KernelBuilder, k: &[Vec<Reg>], coefficients: &[f64]| -> Vec<Reg> {
        (0..n_state)
            .map(|i| {
                let mut sum: Option<Reg> = None;
                for (stage, &c) in k.iter().zip(coefficients) {
                    if c == 0.0 {
                        continue;
                    }
                    let cr = b.num(c);
                    let term = b.push(KernelOp::Mul(cr, stage.get(i).copied().unwrap_or(0)));
                    sum = Some(match sum {
                        | Some(s) => b.push(KernelOp::Add(s, term)),
                        | None => term,
                    });
                }
                let yi = y.get(i).copied().unwrap_or(0);
                match sum {
                    | Some(s) => {
                        let step = b.push(KernelOp::Mul(h, s));
                        b.push(KernelOp::Add(yi, step))
                    },
                    | None => yi,
                }
            })
            .collect()
    };
    let mut k: Vec<Vec<Reg>> = vec![eval(&mut b, t, &y)];
    let mut y_new = y.clone();
    for (row, &c) in A.iter().zip(&C) {
        let stage_state = combine(&mut b, &k, row);
        let cr = b.num(c);
        let ch = b.push(KernelOp::Mul(cr, h));
        let time = b.push(KernelOp::Add(t, ch));
        k.push(eval(&mut b, time, &stage_state));
        y_new = stage_state;
    }
    // `y_new` is the fifth-order solution (the last row of A is the weights).
    // Error estimate in the maximum norm, scaled by the tolerance.
    let mut worst: Option<Reg> = None;
    for i in 0..n_state {
        let mut sum: Option<Reg> = None;
        for (stage, &e) in k.iter().zip(&E) {
            if e == 0.0 {
                continue;
            }
            let er = b.num(e);
            let term = b.push(KernelOp::Mul(er, stage.get(i).copied().unwrap_or(0)));
            sum = Some(match sum {
                | Some(s) => b.push(KernelOp::Add(s, term)),
                | None => term,
            });
        }
        let Some(sum) = sum else {
            continue;
        };
        let raw = b.push(KernelOp::Mul(h, sum));
        let err = b.push(KernelOp::Unary(Intrinsic::Abs, raw));
        let a0 = b.push(KernelOp::Unary(
            Intrinsic::Abs,
            y.get(i).copied().unwrap_or(0),
        ));
        let a1 = b.push(KernelOp::Unary(
            Intrinsic::Abs,
            y_new.get(i).copied().unwrap_or(0),
        ));
        let bigger = b.push(KernelOp::Cmp(Cmp::Gt, a1, a0));
        let size = b.push(KernelOp::Select(bigger, a1, a0));
        let rel = b.push(KernelOp::Mul(rtol, size));
        let scale = b.push(KernelOp::Add(atol, rel));
        let ratio = b.push(KernelOp::Div(err, scale));
        worst = Some(match worst {
            | Some(w) => {
                // The larger of the two, and NaN if either is: a NaN `w`
                // survives because no comparison replaces it, a NaN `ratio`
                // is taken explicitly. A failed stage must reject the step.
                let larger = b.push(KernelOp::Cmp(Cmp::Gt, ratio, w));
                let m = b.push(KernelOp::Select(larger, ratio, w));
                let nan = b.push(KernelOp::Cmp(Cmp::Ne, ratio, ratio));
                b.push(KernelOp::Select(nan, ratio, m))
            },
            | None => ratio,
        });
    }
    let worst = worst.unwrap_or_else(|| b.num(0.0));
    let one = b.num(1.0);
    let accept_raw = b.push(KernelOp::Cmp(Cmp::Le, worst, one));
    let accept = b.push(KernelOp::Mul(active, accept_raw));
    // Step factor 0.9 * err^(-1/5), kept within [0.2, 5]; a NaN estimate
    // shrinks the step.
    let exponent = b.num(-0.2);
    let tiny = b.num(1.0e-30);
    let small = b.push(KernelOp::Cmp(Cmp::Lt, worst, tiny));
    let floor = b.push(KernelOp::Select(small, tiny, worst));
    let power = b.push(KernelOp::PowF(floor, exponent));
    let safety = b.num(0.9);
    let raw_factor = b.push(KernelOp::Mul(safety, power));
    let lo = b.num(0.2);
    let hi = b.num(5.0);
    let above = b.push(KernelOp::Cmp(Cmp::Ge, raw_factor, lo));
    let f1 = b.push(KernelOp::Select(above, raw_factor, lo));
    let ordered = b.push(KernelOp::Cmp(Cmp::Eq, raw_factor, raw_factor));
    let f2 = b.push(KernelOp::Select(ordered, f1, lo));
    let below = b.push(KernelOp::Cmp(Cmp::Le, f2, hi));
    let factor = b.push(KernelOp::Select(below, f2, hi));

    let mut out = Vec::with_capacity(n_state + 2);
    for i in 0..n_state {
        let old = y.get(i).copied().unwrap_or(0);
        let new = y_new.get(i).copied().unwrap_or(0);
        out.push(b.push(KernelOp::Select(accept, new, old)));
    }
    let t_new = b.push(KernelOp::Add(t, h));
    out.push(b.push(KernelOp::Select(accept, t_new, t)));
    let h_new = b.push(KernelOp::Mul(h, factor));
    out.push(b.push(KernelOp::Select(active, h_new, h_in)));
    b.finish(out)
}

/// Integrates `y' = f([t,] y; params)` from `t0` to `t1 > t0` for many
/// initial conditions at once with adaptive Dormand-Prince 5(4) steps.
///
/// Every trajectory carries its own time and step size on the device
/// ([`dopri5_step`]), so a stiff row takes small steps while the others take
/// large ones. The host only asks, every few launches, for the minimum of
/// the time column ([`reduce_column`]): one number crosses per check. The
/// result for a row does not depend on the rest of the batch.
///
/// # Errors
/// Any [`DeviceError`]; [`DeviceError::Mismatch`] for `t1 < t0`, a wrong
/// parameter count or ragged initial conditions; [`DeviceError::Backend`]
/// when some trajectory has not arrived after `tolerance.max_launches`
/// launches (a blow-up, or a tolerance the float type cannot meet).
pub fn ensemble_dopri5(
    device: &dyn Device,
    rhs: &KernelIr,
    params: &[f64],
    y0: &[&[f64]],
    t0: f64,
    t1: f64,
    tolerance: &Tolerance,
) -> Result<Ensemble, DeviceError> {
    if params.len() != rhs.params {
        return Err(DeviceError::Mismatch(format!(
            "rhs has {} parameters, {} given",
            rhs.params,
            params.len()
        )));
    }
    if !(t1 >= t0) {
        return Err(DeviceError::Mismatch(
            "ensemble_dopri5 integrates forward: t1 must not be before t0".into(),
        ));
    }
    let n = y0.len();
    let step = dopri5_step(rhs, n)?;
    let rows = y0.first().map_or(0, |c| c.len());
    if y0.iter().any(|c| c.len() != rows) {
        return Err(DeviceError::Mismatch(
            "initial conditions differ in length".into(),
        ));
    }
    if rows == 0 || t1 == t0 {
        return Ok(Ensemble {
            y: y0.iter().map(|c| c.to_vec()).collect(),
            launches: 0,
        });
    }
    let dtype = preferred_dtype(device.info());
    // A row is finished once it is within a few rounding errors of `t1`.
    let eps = match dtype {
        | super::Dtype::F32 => f64::from(f32::EPSILON),
        | super::Dtype::F64 => f64::EPSILON,
    };
    let t_done = t1 - 8.0 * eps * t1.abs().max(t0.abs()).max(1.0);
    let h0 = tolerance
        .first_step
        .filter(|h| *h > 0.0)
        .unwrap_or((t1 - t0) / 64.0);
    let mut flat = Vec::with_capacity((n + 2) * rows);
    for c in y0 {
        flat.extend_from_slice(c);
    }
    flat.extend(std::iter::repeat_n(t0, rows));
    flat.extend(std::iter::repeat_n(h0, rows));
    let mut p = params.to_vec();
    p.extend([t1, t_done, tolerance.rtol, tolerance.atol]);
    let kernel = device.compile(&step)?;
    let mut a = device.upload(&flat, dtype)?;
    let mut b = device.alloc((n + 2) * rows, dtype)?;
    let lens = vec![rows; n + 2];
    let mut launches = 0_usize;
    // Checking costs a few tiny launches and one download; do it rarely
    // while far from the end and never more than every fourth step.
    let check_every = 4_usize;
    loop {
        for _ in 0..check_every {
            device.launch(&*kernel, rows, &a, &p, &mut b)?;
            std::mem::swap(&mut a, &mut b);
            launches += 1;
        }
        let slowest = reduce_column(device, ReduceOp::Min, &a, &lens, n)?;
        if slowest >= t_done {
            break;
        }
        if launches >= tolerance.max_launches || slowest.is_nan() {
            return Err(DeviceError::Backend(format!(
                "ensemble_dopri5: after {launches} launches the slowest trajectory is at t = \
                 {slowest}, short of {t1}"
            )));
        }
    }
    let mut host = vec![0.0; (n + 2) * rows];
    device.download(&a, &mut host)?;
    Ok(Ensemble {
        y: host.chunks(rows).take(n).map(<[f64]>::to_vec).collect(),
        launches,
    })
}
