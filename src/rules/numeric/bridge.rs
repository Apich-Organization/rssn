//! Numeric kernels exposed as operators, table driven.
//!
//! Each entry decodes its arguments (numbers, rectangular nested lists,
//! symbols naming a wavelet) into [`Data`], calls a kernel of
//! [`crate::kernels`] and encodes the answer back into a term. Like the rest
//! of the `numeric` rule set the operators compute only from literal numbers:
//! a symbolic argument, a ragged list, a non-finite answer or a rejected
//! parameter leaves the call as a request.
//!
//! | group | operators |
//! |---|---|
//! | robust statistics | `mad`, `mad_sigma`, `trimmed_mean`, `winsorized_mean`, `huber_location`, `biweight_location`, `ecdf`, `kde`, `bandwidth_silverman`, `bandwidth_scott` |
//! | signal | `hann_window`, `hamming_window`, `blackman_window`, `cross_correlation`, `fft_abs`, `fft_re`, `fft_im`, `ifft_re`, `ifft_im`, `ifft_abs` |
//! | wavelets | `dwt_approx`, `dwt_detail`, `idwt` (the wavelet is a symbol: `haar`, `db2` .. `db10`, `sym4` .. `sym6`, `bior1_1` .. `bior4_4`) |
//! | finite-difference PDEs | `heat_fd`, `wave_fd`, `poisson_fd`, `advect_fd`, and the siblings `heat_fd_bc`, `wave_fd_bc`, `poisson_fd_bc` (Dirichlet or Neumann walls, per face or axis) and `advect_fd_upwind` |
//! | tensors | `tensor_norm`, `tensor_inner`, `tensor_outer`, `tensordot`, `tensor_contract` |
//!
//! The finite-difference operators are **dimension free**: the number of
//! space dimensions is the nesting depth of the initial data (a list is 1D, a
//! list of lists 2D, three levels 3D), the domain is the unit interval,
//! square or cube, and `list[i][j][k]` is the node `(x_i, y_j, z_k)`. They
//! replace what would otherwise be a `_1d`/`_2d`/`_3d` triple of each.

use rustfft::num_complex::Complex;

use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::op::core;
use crate::kernel::diag;
use crate::kernel::diag::DiagKind;
use crate::kernels::functional_analysis;
use crate::kernels::krylov;
use crate::kernels::krylov_ext;
use crate::kernels::lowrank;
use crate::kernels::pde;
use crate::kernels::robust;
use crate::kernels::signal;
use crate::kernels::tensor;
use crate::kernels::wavelets;

/// Longest signal or sample the bridge accepts.
pub(super) const MAX_LEN: usize = 100_000;

/// A decoded argument or result.
#[derive(Clone, Debug, PartialEq)]
pub(super) enum Data {
    /// A number.
    Scalar(f64),
    /// A rectangular nested list, flattened row-major.
    Array {
        /// Extent of each nesting level.
        shape: Vec<usize>,
        /// The numbers, last index fastest.
        data: Vec<f64>,
    },
    /// A symbol used as a name.
    Name(String),
    /// A list of unlike parts: a ragged list of lists (the coefficient
    /// arrays of a multilevel wavelet transform, an adjacency list) or a
    /// structured result (a decomposition's core and factors).
    Tree(Vec<Self>),
}

impl Data {
    pub(super) fn vector(data: Vec<f64>) -> Self {
        Self::Array {
            shape: vec![data.len()],
            data,
        }
    }
}

/// What an argument is, for a diagnostic message.
struct Shown<'a>(Option<&'a Data>);

impl std::fmt::Display for Shown<'_> {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        match self.0 {
            | None => f.write_str("nothing (the argument is missing)"),
            | Some(Data::Scalar(x)) => write!(f, "the number {x}"),
            | Some(Data::Array { shape, .. }) => {
                let extents: Vec<String> = shape.iter().map(ToString::to_string).collect();
                write!(
                    f,
                    "a {}-level list of shape {}",
                    shape.len(),
                    extents.join(" x ")
                )
            },
            | Some(Data::Name(n)) => write!(f, "the symbol `{n}`"),
            | Some(Data::Tree(parts)) => write!(f, "a ragged list of {} parts", parts.len()),
        }
    }
}

/// Records why argument `at` (zero-based) was rejected and returns `None`.
pub(super) fn bad<T>(
    at: usize,
    kind: DiagKind,
    message: std::fmt::Arguments<'_>,
) -> Option<T> {
    diag::fail(kind, Some(at), message)
}

/// `Some(())` when `ok`; else records the reason (blaming argument `at`, when
/// one is at fault) and returns `None`. Written `require(..)?`.
pub(super) fn require(
    ok: bool,
    at: Option<usize>,
    kind: DiagKind,
    message: std::fmt::Arguments<'_>,
) -> Option<()> {
    if ok {
        Some(())
    } else {
        diag::fail(kind, at, message)
    }
}

/// Turns the error of a kernel into a recorded reason.
pub(super) trait Reason<T> {
    /// The value, or `None` after recording the error's text as a `kind`
    /// finding (blaming argument `at` when given).
    fn reason(
        self,
        kind: DiagKind,
        at: Option<usize>,
    ) -> Option<T>;
}

impl<T, E: std::fmt::Display> Reason<T> for Result<T, E> {
    fn reason(
        self,
        kind: DiagKind,
        at: Option<usize>,
    ) -> Option<T> {
        match self {
            | Ok(value) => Some(value),
            | Err(error) => diag::fail(kind, at, format_args!("{error}")),
        }
    }
}

/// Deepest nesting of lists the decoder follows.
const MAX_DEPTH: usize = 8;

/// The decoded form of `node`, if it is a number, a list of numbers or of
/// lists (rectangular, else a [`Data::Tree`]), or a symbol.
pub(super) fn decode(
    graph: &Graph,
    node: NodeId,
) -> Option<Data> {
    decode_at(graph, node, 0)
}

fn decode_at(
    graph: &Graph,
    node: NodeId,
    depth: usize,
) -> Option<Data> {
    if depth > MAX_DEPTH {
        return None;
    }
    // A ragged list is tried only after the rectangular reading fails, and
    // only near the top, so hostile nesting cannot make the decoder
    // exponential.
    decode_uniform(graph, node, depth).or_else(|| {
        if depth > 3 {
            return None;
        }
        let list = graph
            .enodes(graph.find(node))
            .find(|&n| graph.op(n) == core::LIST)?;
        let items = graph.children(list);
        if items.is_empty() || items.len() > MAX_LEN {
            return None;
        }
        let parts: Option<Vec<Data>> = items
            .iter()
            .map(|&item| decode_at(graph, item, depth + 1))
            .collect();
        let parts = parts?;
        (!parts.iter().any(|p| matches!(p, Data::Name(_)))).then_some(Data::Tree(parts))
    })
}

fn decode_uniform(
    graph: &Graph,
    node: NodeId,
    depth: usize,
) -> Option<Data> {
    if let Some(n) = graph.number_of(node) {
        return Some(Data::Scalar(n.to_f64()));
    }
    // The argument may be a call that the engine has already proved equal to
    // a list (`ifft_re(fft_re(x), ..)`): look for the list node in its class.
    let Some(list) = graph
        .enodes(graph.find(node))
        .find(|&n| graph.op(n) == core::LIST)
    else {
        let symbol = graph.symbol_of(node)?;
        return Some(Data::Name(graph.interner().symbol_name(symbol).to_owned()));
    };
    let items = graph.children(list);
    let mut shape = vec![items.len()];
    let mut flat = Vec::new();
    let mut inner: Option<Vec<usize>> = None;
    for &item in items {
        match decode_at(graph, item, depth + 1)? {
            | Data::Scalar(x) => {
                if inner.as_ref().is_some_and(|s| !s.is_empty()) {
                    return None;
                }
                inner = Some(Vec::new());
                flat.push(x);
            },
            | Data::Array { shape: s, data } => {
                if inner.as_ref().is_some_and(|prev| *prev != s) {
                    return None;
                }
                inner = Some(s);
                flat.extend(data);
            },
            | Data::Name(_) | Data::Tree(_) => return None,
        }
        if flat.len() > MAX_LEN * 4 {
            return None;
        }
    }
    if let Some(s) = inner {
        shape.extend(s);
    }
    Some(Data::Array { shape, data: flat })
}

/// Encodes a result as a term.
fn encode(
    graph: &mut Graph,
    data: &Data,
) -> Option<NodeId> {
    match data {
        | Data::Scalar(x) => x.is_finite().then(|| graph.float(*x)),
        | Data::Name(_) => None,
        | Data::Tree(parts) => {
            let items: Option<Vec<NodeId>> = parts.iter().map(|p| encode(graph, p)).collect();
            Some(graph.node(core::LIST, &items?))
        },
        | Data::Array { shape, data } => {
            if data.iter().any(|x| !x.is_finite()) {
                return None;
            }
            let leaves: Vec<NodeId> = data.iter().map(|&x| graph.float(x)).collect();
            Some(nest(graph, shape, &leaves))
        },
    }
}

fn nest(
    graph: &mut Graph,
    shape: &[usize],
    leaves: &[NodeId],
) -> NodeId {
    let Some((&extent, rest)) = shape.split_first() else {
        return leaves[0];
    };
    let stride: usize = rest.iter().product();
    let items: Vec<NodeId> = (0..extent)
        .map(|i| {
            if rest.is_empty() {
                leaves[i]
            } else {
                nest(graph, rest, &leaves[i * stride..(i + 1) * stride])
            }
        })
        .collect();
    graph.node(core::LIST, &items)
}

// ---------------------------------------------------------------------------
// Argument access
// ---------------------------------------------------------------------------

pub(super) fn scalar(
    args: &[Data],
    at: usize,
) -> Option<f64> {
    match args.get(at) {
        | Some(Data::Scalar(x)) if x.is_finite() => Some(*x),
        | Some(Data::Scalar(x)) => {
            bad(
                at,
                DiagKind::Domain,
                format_args!("expected a finite number, got {x}"),
            )
        },
        | other => {
            bad(
                at,
                DiagKind::Shape,
                format_args!("expected a number, got {}", Shown(other)),
            )
        },
    }
}

pub(super) fn vector(
    args: &[Data],
    at: usize,
) -> Option<&[f64]> {
    match args.get(at) {
        | Some(Data::Array { shape, data }) if shape.len() == 1 => {
            require(
                data.len() <= MAX_LEN,
                Some(at),
                DiagKind::Domain,
                format_args!(
                    "a list may hold at most {MAX_LEN} numbers, got {}",
                    data.len()
                ),
            )?;
            require(
                data.iter().all(|x| x.is_finite()),
                Some(at),
                DiagKind::Domain,
                format_args!("the list holds a number that is not finite"),
            )?;
            Some(data)
        },
        | other => {
            bad(
                at,
                DiagKind::Shape,
                format_args!("expected a flat list of numbers, got {}", Shown(other)),
            )
        },
    }
}

/// A non-empty finite sample.
pub(super) fn sample(
    args: &[Data],
    at: usize,
) -> Option<&[f64]> {
    let v = vector(args, at)?;
    require(
        !v.is_empty(),
        Some(at),
        DiagKind::Domain,
        format_args!("the sample must not be empty"),
    )?;
    Some(v)
}

pub(super) fn count(
    args: &[Data],
    at: usize,
    max: usize,
) -> Option<usize> {
    let x = scalar(args, at)?;
    require(
        x >= 0.0 && x.fract() == 0.0 && x <= max as f64,
        Some(at),
        DiagKind::Domain,
        format_args!("expected an integer from 0 to {max}, got {x}"),
    )?;
    Some(x as usize)
}

pub(super) fn name(
    args: &[Data],
    at: usize,
) -> Option<&str> {
    match args.get(at) {
        | Some(Data::Name(n)) => Some(n),
        | other => {
            bad(
                at,
                DiagKind::Shape,
                format_args!("expected a symbol, got {}", Shown(other)),
            )
        },
    }
}

/// The numbers of a one-dimensional argument.
pub(super) fn numbers_of(data: &Data) -> Option<&[f64]> {
    match data {
        | Data::Array { shape, data } if shape.len() == 1 && data.iter().all(|x| x.is_finite()) => {
            Some(data)
        },
        | _ => None,
    }
}

/// The one level below an argument: the rows of a rectangular list, or the
/// parts of a [`Data::Tree`].
pub(super) fn parts_of(data: &Data) -> Option<Vec<Data>> {
    match data {
        | Data::Tree(parts) => Some(parts.clone()),
        | Data::Array { shape, data } if shape.len() >= 2 => {
            let stride: usize = shape[1..].iter().product();
            let rest = &shape[1..];
            Some(
                (0..shape[0])
                    .map(|i| {
                        let chunk = data[i * stride..(i + 1) * stride].to_vec();
                        if rest.len() == 1 {
                            Data::vector(chunk)
                        } else {
                            Data::Array {
                                shape: rest.to_vec(),
                                data: chunk,
                            }
                        }
                    })
                    .collect(),
            )
        },
        | _ => None,
    }
}

/// A list of number lists of any lengths, as vectors.
pub(super) fn rows(data: &Data) -> Option<Vec<Vec<f64>>> {
    parts_of(data)?
        .iter()
        .map(|part| numbers_of(part).map(<[f64]>::to_vec))
        .collect()
}

// reason: kept `Option` so every kernel body is a uniform `Option<Data>` expression.
#[allow(clippy::unnecessary_wraps)]
pub(super) const fn num(x: f64) -> Option<Data> {
    Some(Data::Scalar(x))
}

// reason: kept `Option` so every kernel body is a uniform `Option<Data>` expression.
#[allow(clippy::unnecessary_wraps)]
pub(super) fn vec_out(v: Vec<f64>) -> Option<Data> {
    Some(Data::vector(v))
}

// ---------------------------------------------------------------------------
// Robust statistics
// ---------------------------------------------------------------------------

fn mad(a: &[Data]) -> Option<Data> {
    num(robust::mad(sample(a, 0)?, false))
}

fn mad_sigma(a: &[Data]) -> Option<Data> {
    num(robust::mad(sample(a, 0)?, true))
}

fn trimmed_mean(a: &[Data]) -> Option<Data> {
    num(robust::trimmed_mean(sample(a, 0)?, scalar(a, 1)?))
}

fn winsorized_mean(a: &[Data]) -> Option<Data> {
    num(robust::winsorized_mean(sample(a, 0)?, scalar(a, 1)?))
}

fn huber_location(a: &[Data]) -> Option<Data> {
    let k = scalar(a, 1)?;
    require(
        k > 0.0,
        Some(1),
        DiagKind::Domain,
        format_args!("the tuning constant must be positive, got {k}"),
    )?;
    num(robust::huber_location(sample(a, 0)?, k))
}

fn biweight_location(a: &[Data]) -> Option<Data> {
    let c = scalar(a, 1)?;
    require(
        c > 0.0,
        Some(1),
        DiagKind::Domain,
        format_args!("the tuning constant must be positive, got {c}"),
    )?;
    num(robust::biweight_location(sample(a, 0)?, c))
}

fn ecdf(a: &[Data]) -> Option<Data> {
    num(robust::ecdf(sample(a, 0)?, scalar(a, 1)?))
}

fn kde(a: &[Data]) -> Option<Data> {
    let h = scalar(a, 2)?;
    require(
        h > 0.0,
        Some(2),
        DiagKind::Domain,
        format_args!("the bandwidth must be positive, got {h}"),
    )?;
    num(robust::kde(
        sample(a, 0)?,
        scalar(a, 1)?,
        h,
        robust::Kernel::Gaussian,
    ))
}

fn bandwidth_silverman(a: &[Data]) -> Option<Data> {
    num(robust::bandwidth_silverman(sample(a, 0)?))
}

fn bandwidth_scott(a: &[Data]) -> Option<Data> {
    num(robust::bandwidth_scott(sample(a, 0)?))
}

// ---------------------------------------------------------------------------
// Signal
// ---------------------------------------------------------------------------

fn hann_window(a: &[Data]) -> Option<Data> {
    vec_out(signal::hann_window(count(a, 0, MAX_LEN)?))
}

fn hamming_window(a: &[Data]) -> Option<Data> {
    vec_out(signal::hamming_window(count(a, 0, MAX_LEN)?))
}

fn blackman_window(a: &[Data]) -> Option<Data> {
    vec_out(signal::blackman_window(count(a, 0, MAX_LEN)?))
}

fn cross_correlation(a: &[Data]) -> Option<Data> {
    let (x, y) = (vector(a, 0)?, vector(a, 1)?);
    vec_out(signal::cross_correlation(x, y))
}

fn spectrum(a: &[Data]) -> Option<Vec<Complex<f64>>> {
    let x = vector(a, 0)?;
    let mut input: Vec<Complex<f64>> = x.iter().map(|&v| Complex::new(v, 0.0)).collect();
    Some(signal::fft(&mut input))
}

fn fft_abs(a: &[Data]) -> Option<Data> {
    vec_out(spectrum(a)?.iter().map(|c| c.norm()).collect())
}

fn fft_re(a: &[Data]) -> Option<Data> {
    vec_out(spectrum(a)?.iter().map(|c| c.re).collect())
}

fn fft_im(a: &[Data]) -> Option<Data> {
    vec_out(spectrum(a)?.iter().map(|c| c.im).collect())
}

fn inverse(a: &[Data]) -> Option<Vec<Complex<f64>>> {
    let (re, im) = (vector(a, 0)?, vector(a, 1)?);
    require(
        re.len() == im.len(),
        Some(1),
        DiagKind::Shape,
        format_args!(
            "the imaginary part has {} entries but the real part has {}",
            im.len(),
            re.len()
        ),
    )?;
    let mut input: Vec<Complex<f64>> = re
        .iter()
        .zip(im)
        .map(|(&r, &i)| Complex::new(r, i))
        .collect();
    Some(signal::ifft(&mut input))
}

fn ifft_re(a: &[Data]) -> Option<Data> {
    vec_out(inverse(a)?.iter().map(|c| c.re).collect())
}

fn ifft_im(a: &[Data]) -> Option<Data> {
    vec_out(inverse(a)?.iter().map(|c| c.im).collect())
}

fn ifft_abs(a: &[Data]) -> Option<Data> {
    vec_out(inverse(a)?.iter().map(|c| c.norm()).collect())
}

// ---------------------------------------------------------------------------
// Wavelets
// ---------------------------------------------------------------------------

/// The wavelet a symbol names: the lower-cased variant name (`haar`, `db4`,
/// `sym5`, `bior1_3`).
pub(super) fn wavelet(
    a: &[Data],
    at: usize,
) -> Option<wavelets::Wavelet> {
    let wanted = name(a, at)?.to_ascii_lowercase();
    let found = wavelets::Wavelet::ALL
        .into_iter()
        .find(|w| format!("{w:?}").to_ascii_lowercase() == wanted);
    if found.is_none() {
        return bad(
            at,
            DiagKind::Domain,
            format_args!(
                "unknown wavelet `{wanted}` (use haar, db2 .. db10, sym4 .. sym6 or bior1_1 .. bior4_4)"
            ),
        );
    }
    found
}

/// A wavelet transform's failure as a recorded reason.
pub(super) fn wavelet_reason<T>(
    result: Result<T, wavelets::WaveletError>,
    what: &str,
) -> Option<T> {
    match result {
        | Ok(value) => Some(value),
        | Err(wavelets::WaveletError::OddLength) => {
            bad_op(
                DiagKind::Shape,
                format_args!("{what}: the signal length must be even and nonzero"),
            )
        },
        | Err(wavelets::WaveletError::TooManyLevels) => {
            bad_op(
                DiagKind::Domain,
                format_args!("{what}: more levels than the signal length allows"),
            )
        },
        | Err(wavelets::WaveletError::Shape) => {
            bad_op(
                DiagKind::Shape,
                format_args!(
                    "{what}: the coefficient lists must be nonempty, of equal length and fit together"
                ),
            )
        },
    }
}

/// Records a finding that blames no single argument.
pub(super) fn bad_op<T>(
    kind: DiagKind,
    message: std::fmt::Arguments<'_>,
) -> Option<T> {
    diag::fail(kind, None, message)
}

fn dwt_approx(a: &[Data]) -> Option<Data> {
    let (x, w) = (vector(a, 0)?, wavelet(a, 1)?);
    vec_out(wavelet_reason(wavelets::dwt(x, w), "dwt_approx")?.0)
}

fn dwt_detail(a: &[Data]) -> Option<Data> {
    let (x, w) = (vector(a, 0)?, wavelet(a, 1)?);
    vec_out(wavelet_reason(wavelets::dwt(x, w), "dwt_detail")?.1)
}

fn idwt(a: &[Data]) -> Option<Data> {
    let (approx, detail, w) = (vector(a, 0)?, vector(a, 1)?, wavelet(a, 2)?);
    vec_out(wavelet_reason(wavelets::idwt(approx, detail, w), "idwt")?)
}

// ---------------------------------------------------------------------------
// Finite-difference PDEs, any number of space dimensions
// ---------------------------------------------------------------------------

/// Most nodes a solve may use per space dimension, and most node-steps.
const MAX_NODES: [usize; 3] = [20_000, 16_384, 4_096];
const MAX_WORK: usize = 5_000_000;

fn field(
    a: &[Data],
    at: usize,
) -> Option<(&[usize], &[f64])> {
    match a.get(at) {
        | Some(Data::Array { shape, data }) if (1..=3).contains(&shape.len()) => {
            require(
                shape.iter().all(|&n| n >= 3),
                Some(at),
                DiagKind::Domain,
                format_args!(
                    "every extent of a grid needs at least 3 nodes, got {}",
                    Shown(a.get(at))
                ),
            )?;
            require(
                data.len() <= MAX_NODES[shape.len() - 1],
                Some(at),
                DiagKind::Domain,
                format_args!(
                    "a {}D grid is limited to {} nodes, got {}",
                    shape.len(),
                    MAX_NODES[shape.len() - 1],
                    data.len()
                ),
            )?;
            require(
                data.iter().all(|x| x.is_finite()),
                Some(at),
                DiagKind::Domain,
                format_args!("the grid holds a number that is not finite"),
            )?;
            Some((shape, data))
        },
        | Some(Data::Array { shape, .. }) => {
            bad(
                at,
                DiagKind::Shape,
                format_args!(
                    "a field is a 1D, 2D or 3D nested list, got {} nesting levels",
                    shape.len()
                ),
            )
        },
        | other => {
            bad(
                at,
                DiagKind::Shape,
                format_args!("expected a nested list (the field), got {}", Shown(other)),
            )
        },
    }
}

/// Nearest node of a coordinate in `[0, 1]` on `n` nodes (`periodic`: `n`
/// cells, node `i` at `i / n`).
fn node(
    x: f64,
    n: usize,
    periodic: bool,
) -> usize {
    let cells = if periodic { n } else { n - 1 } as f64;
    ((x * cells).round().max(0.0) as usize).min(n - 1)
}

/// A sampler of the nested list at a point of the unit box.
fn sampler<'a>(
    shape: &'a [usize],
    data: &'a [f64],
    periodic: bool,
) -> impl Fn(&[f64]) -> f64 + 'a {
    move |p: &[f64]| {
        let mut at = 0;
        for (d, &n) in shape.iter().enumerate() {
            at = at * n + node(p[d], n, periodic);
        }
        data[at]
    }
}

fn unit_grids(shape: &[usize]) -> (pde::Grid1D, pde::Grid2D, pde::Grid3D) {
    let n = |d: usize| shape.get(d).copied().unwrap_or(3);
    (
        pde::Grid1D {
            x0: 0.0,
            x1: 1.0,
            n: n(0),
        },
        pde::Grid2D {
            x0: 0.0,
            x1: 1.0,
            y0: 0.0,
            y1: 1.0,
            nx: n(0),
            ny: n(1),
        },
        pde::Grid3D {
            x0: 0.0,
            x1: 1.0,
            y0: 0.0,
            y1: 1.0,
            z0: 0.0,
            z1: 1.0,
            nx: n(0),
            ny: n(1),
            nz: n(2),
        },
    )
}

/// The solution as a nested list `out[i][j][k]`.
fn solution(sol: pde::PdeSolution) -> Data {
    match sol {
        | pde::PdeSolution::D1(f) => Data::vector(f.u),
        | pde::PdeSolution::D2(f) => {
            let mut data = Vec::with_capacity(f.u.len());
            for i in 0..f.nx {
                for j in 0..f.ny {
                    data.push(f.at(i, j));
                }
            }
            Data::Array {
                shape: vec![f.nx, f.ny],
                data,
            }
        },
        | pde::PdeSolution::D3(f) => {
            let mut data = Vec::with_capacity(f.u.len());
            for i in 0..f.nx {
                for j in 0..f.ny {
                    for k in 0..f.nz {
                        data.push(f.at(i, j, k));
                    }
                }
            }
            Data::Array {
                shape: vec![f.nx, f.ny, f.nz],
                data,
            }
        },
    }
}

/// How the boundary value of a finite-difference solve is imposed on every
/// face of the unit box.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Condition {
    /// `u = value`.
    Dirichlet,
    /// `du/dn = value`, `n` the outward normal (`0` is an insulated wall).
    Neumann,
}

/// The condition and value of every face of the unit box, in the order of the
/// kernel: `[x-, x+, y-, y+, z-, z+]` truncated to `2 * rank` faces. `x` is the
/// first nesting level of the field, so the faces `x-` and `x+` are
/// `list[0][..]` and `list[n - 1][..]`, and `y`, `z` follow.
type Faces = Vec<(Condition, f64)>;

/// The conditions a boundary symbol names for a field of `rank` space
/// dimensions, one per face.
///
/// * `dirichlet` or `neumann`: every face;
/// * a pattern of `rank` letters from `d`/`n` (any case): one per axis, both
///   faces of the axis, e.g. `dn` in 2D is Dirichlet on the `x` faces and
///   Neumann on the `y` faces;
/// * a pattern of `2 * rank` letters: one per face in the order `x-, x+, y-,
///   y+, z-, z+`, e.g. `ddnn` in 2D.
fn conditions(
    a: &[Data],
    at: usize,
    rank: usize,
) -> Option<Vec<Condition>> {
    let symbol = name(a, at)?.to_ascii_lowercase();
    let letter = |c: char| {
        match c {
            | 'd' => Some(Condition::Dirichlet),
            | 'n' => Some(Condition::Neumann),
            | _ => None,
        }
    };
    let letters: Option<Vec<Condition>> = symbol.chars().map(letter).collect();
    match (symbol.as_str(), letters) {
        | ("dirichlet", _) => Some(vec![Condition::Dirichlet; 2 * rank]),
        | ("neumann", _) => Some(vec![Condition::Neumann; 2 * rank]),
        | (_, Some(each)) if each.len() == rank => {
            Some(each.iter().flat_map(|&c| [c, c]).collect())
        },
        | (_, Some(each)) if each.len() == 2 * rank => Some(each),
        | _ => {
            bad(
                at,
                DiagKind::Domain,
                format_args!(
                    "unknown boundary condition `{symbol}`: use dirichlet, neumann, or a pattern \
                     of d/n letters, {rank} (one per axis) or {} (one per face x-, x+, y-, y+, \
                     z-, z+) of them for a {rank}D field",
                    2 * rank
                ),
            )
        },
    }
}

/// The boundary value of every face: a number for all of them, a list of
/// `rank` numbers (one per axis) or of `2 * rank` (one per face).
fn face_values(
    a: &[Data],
    at: usize,
    rank: usize,
) -> Option<Vec<f64>> {
    match a.get(at) {
        | Some(Data::Scalar(_)) => Some(vec![scalar(a, at)?; 2 * rank]),
        | Some(Data::Array { shape, .. }) if shape.len() == 1 => {
            let each = vector(a, at)?;
            if each.len() == rank {
                Some(each.iter().flat_map(|&v| [v, v]).collect())
            } else if each.len() == 2 * rank {
                Some(each.to_vec())
            } else {
                bad(
                    at,
                    DiagKind::Shape,
                    format_args!(
                        "a {rank}D field takes one boundary value, {rank} (per axis) or {} (per \
                         face), got {}",
                        2 * rank,
                        each.len()
                    ),
                )
            }
        },
        | other => {
            bad(
                at,
                DiagKind::Shape,
                format_args!(
                    "expected a boundary value (a number or a list of numbers), got {}",
                    Shown(other)
                ),
            )
        },
    }
}

/// The faces of a `_bc` operator: the condition symbol at `kind_at`, the value
/// at `value_at`.
fn faces(
    a: &[Data],
    kind_at: usize,
    value_at: usize,
    rank: usize,
) -> Option<Faces> {
    let kinds = conditions(a, kind_at, rank)?;
    let values = face_values(a, value_at, rank)?;
    Some(kinds.into_iter().zip(values).collect())
}

/// Dirichlet `value` on every face.
fn all_dirichlet(
    value: f64,
    rank: usize,
) -> Faces {
    vec![(Condition::Dirichlet, value); 2 * rank]
}

fn bcs<const N: usize>(faces: &[(Condition, f64)]) -> [pde::Boundary; N] {
    std::array::from_fn(|i| {
        match faces[i] {
            | (Condition::Dirichlet, value) => pde::Boundary::dirichlet_const(value),
            | (Condition::Neumann, value) => pde::Boundary::neumann_const(value),
        }
    })
}

fn steps_within(
    a: &[Data],
    at: usize,
    nodes: usize,
) -> Option<usize> {
    let steps = scalar(a, at)?;
    require(
        steps >= 1.0 && steps.fract() == 0.0 && steps <= 1e6,
        Some(at),
        DiagKind::Domain,
        format_args!("the step count must be an integer from 1 to 1000000, got {steps}"),
    )?;
    let steps = steps as usize;
    require(
        steps
            .checked_mul(nodes)
            .is_some_and(|work| work <= MAX_WORK),
        Some(at),
        DiagKind::Domain,
        format_args!(
            "{steps} steps on {nodes} nodes exceeds the work limit of {MAX_WORK} node-steps"
        ),
    )?;
    Some(steps)
}

/// A solver's error as a recorded reason, in the vocabulary of [`DiagKind`].
fn pde_reason<T>(
    result: Result<T, pde::PdeError>,
    solver: &str,
) -> Option<T> {
    let error = match result {
        | Ok(value) => return Some(value),
        | Err(error) => error,
    };
    let kind = match error {
        | pde::PdeError::InvalidInput(_)
        | pde::PdeError::Unstable { .. }
        | pde::PdeError::NonFinite => DiagKind::Domain,
        | pde::PdeError::Unsupported(_) => DiagKind::Unsupported,
        | pde::PdeError::SolveFailed => DiagKind::Singular,
        | pde::PdeError::Incompatible { .. } => DiagKind::Incompatible,
        | pde::PdeError::NotConverged { .. } | pde::PdeError::Ode(_) => DiagKind::NotConverged,
    };
    bad_op(kind, format_args!("{solver}: {error}"))
}

/// `heat_fd(alpha, u0, boundary_value, t_end, steps)`: Dirichlet walls.
fn heat_fd(a: &[Data]) -> Option<Data> {
    let (shape, _) = field(a, 1)?;
    let faces = all_dirichlet(scalar(a, 2)?, shape.len());
    heat(a, &faces, 3)
}

/// `heat_fd_bc(alpha, u0, kind, boundary_value, t_end, steps)`: `kind` is
/// `dirichlet`, `neumann` or a per-axis or per-face pattern of `d`/`n`
/// letters (see [`conditions`]); `boundary_value` is a number, a list per axis
/// or a list per face.
fn heat_fd_bc(a: &[Data]) -> Option<Data> {
    let (shape, _) = field(a, 1)?;
    let faces = faces(a, 2, 3, shape.len())?;
    heat(a, &faces, 4)
}

/// The heat solve; `t_at` is the position of `t_end`, `steps` follows it.
fn heat(
    a: &[Data],
    faces: &[(Condition, f64)],
    t_at: usize,
) -> Option<Data> {
    let (alpha, t_end) = (scalar(a, 0)?, scalar(a, t_at)?);
    let (shape, data) = field(a, 1)?;
    let steps = steps_within(a, t_at + 1, data.len())?;
    let (g1, g2, g3) = unit_grids(shape);
    let at = sampler(shape, data, false);
    let problem = match shape.len() {
        | 1 => {
            let u0 = |x: f64| at(&[x]);
            return Some(solution(pde_reason(
                pde::pde_solver(&pde::PdeProblem::Heat1D {
                    alpha,
                    grid: g1,
                    bc: bcs(faces),
                    u0: &u0,
                    t_end,
                    steps,
                }),
                "heat",
            )?));
        },
        | 2 => {
            let u0 = |x: f64, y: f64| at(&[x, y]);
            return Some(solution(pde_reason(
                pde::pde_solver(&pde::PdeProblem::Heat2D {
                    alpha,
                    grid: g2,
                    bc: bcs(faces),
                    u0: &u0,
                    t_end,
                    steps,
                }),
                "heat",
            )?));
        },
        | _ => {
            let u0 = |x: f64, y: f64, z: f64| at(&[x, y, z]);
            pde::pde_solver(&pde::PdeProblem::Heat3D {
                alpha,
                grid: g3,
                bc: bcs(faces),
                u0: &u0,
                t_end,
                steps,
            })
        },
    };
    Some(solution(pde_reason(problem, "heat")?))
}

/// `wave_fd(c, u0, v0, boundary_value, t_end, steps)`: Dirichlet walls.
fn wave_fd(a: &[Data]) -> Option<Data> {
    let (shape, _) = field(a, 1)?;
    let faces = all_dirichlet(scalar(a, 3)?, shape.len());
    wave(a, &faces, 4)
}

/// `wave_fd_bc(c, u0, v0, kind, boundary_value, t_end, steps)`: `kind` and
/// `boundary_value` as for [`heat_fd_bc`]; `dirichlet` is a fixed end,
/// `neumann` a free, reflecting one.
fn wave_fd_bc(a: &[Data]) -> Option<Data> {
    let (shape, _) = field(a, 1)?;
    let faces = faces(a, 3, 4, shape.len())?;
    wave(a, &faces, 5)
}

/// The wave solve; `t_at` is the position of `t_end`, `steps` follows it.
fn wave(
    a: &[Data],
    faces: &[(Condition, f64)],
    t_at: usize,
) -> Option<Data> {
    let (c, t_end) = (scalar(a, 0)?, scalar(a, t_at)?);
    let (shape, data) = field(a, 1)?;
    let (vshape, vdata) = field(a, 2)?;
    require(
        vshape == shape,
        Some(2),
        DiagKind::Shape,
        format_args!("the initial velocity has shape {vshape:?} but the displacement {shape:?}"),
    )?;
    let steps = steps_within(a, t_at + 1, data.len())?;
    let (g1, g2, g3) = unit_grids(shape);
    let (at, vt) = (sampler(shape, data, false), sampler(shape, vdata, false));
    let solved = match shape.len() {
        | 1 => {
            let (u0, v0) = (|x: f64| at(&[x]), |x: f64| vt(&[x]));
            pde::pde_solver(&pde::PdeProblem::Wave1D {
                c,
                grid: g1,
                bc: bcs(faces),
                u0: &u0,
                v0: &v0,
                t_end,
                steps,
            })
        },
        | 2 => {
            let (u0, v0) = (|x: f64, y: f64| at(&[x, y]), |x: f64, y: f64| vt(&[x, y]));
            pde::pde_solver(&pde::PdeProblem::Wave2D {
                c,
                grid: g2,
                bc: bcs(faces),
                u0: &u0,
                v0: &v0,
                t_end,
                steps,
            })
        },
        | _ => {
            let (u0, v0) = (
                |x: f64, y: f64, z: f64| at(&[x, y, z]),
                |x: f64, y: f64, z: f64| vt(&[x, y, z]),
            );
            pde::pde_solver(&pde::PdeProblem::Wave3D {
                c,
                grid: g3,
                bc: bcs(faces),
                u0: &u0,
                v0: &v0,
                t_end,
                steps,
            })
        },
    };
    Some(solution(pde_reason(solved, "wave")?))
}

/// `poisson_fd(f, boundary_value)`: `-Laplace(u) = f` with `u = boundary_value`
/// on the whole boundary of the unit box.
fn poisson_fd(a: &[Data]) -> Option<Data> {
    let (shape, _) = field(a, 0)?;
    let faces = all_dirichlet(scalar(a, 1)?, shape.len());
    poisson(a, &faces)
}

/// `poisson_fd_bc(f, kind, boundary_value)`: `-Laplace(u) = f` on the unit box
/// with `kind` and `boundary_value` as for [`heat_fd_bc`]: Dirichlet
/// (`u = value`) or Neumann (`du/dn = value`, `n` the outward normal) per face
/// or per axis.
///
/// With at least one Dirichlet face the problem is regular. With Neumann on
/// **every** face it is singular (`u + c` solves it with `u`) and solvable only
/// when the source and the boundary fluxes balance: `integral f + integral of
/// the outward flux = 0`, here its discrete form with trapezoid weights (the
/// weights that make the ghost-node operator symmetric), to a relative `1e-8`.
/// Data that violate it are an `Incompatible` finding; otherwise the system is
/// solved in the space orthogonal to the constants and the **zero-mean**
/// solution is returned: the trapezoid-weighted mean of the answer is zero (the
/// discrete `integral u = 0`). Add the constant of your choice.
fn poisson_fd_bc(a: &[Data]) -> Option<Data> {
    let (shape, _) = field(a, 0)?;
    let faces = faces(a, 1, 2, shape.len())?;
    poisson(a, &faces)
}

/// Relative tolerance of the conjugate-gradient solve of the pure-Neumann
/// problem.
const NEUMANN_CG_TOL: f64 = 1e-11;

fn poisson(
    a: &[Data],
    faces: &[(Condition, f64)],
) -> Option<Data> {
    let (shape, data) = field(a, 0)?;
    let (g1, g2, g3) = unit_grids(shape);
    let at = sampler(shape, data, false);
    let solver = if faces.iter().all(|&(c, _)| c == Condition::Neumann) {
        pde::PoissonSolver::MeanFree { tol: NEUMANN_CG_TOL }
    } else {
        pde::PoissonSolver::SparseLu
    };
    let solved = match shape.len() {
        | 1 => {
            let f = |x: f64| at(&[x]);
            pde::pde_solver(&pde::PdeProblem::Poisson1D {
                grid: g1,
                bc: bcs(faces),
                f: &f,
                solver,
            })
        },
        | 2 => {
            let f = |x: f64, y: f64| at(&[x, y]);
            pde::pde_solver(&pde::PdeProblem::Poisson2D {
                grid: g2,
                bc: bcs(faces),
                f: &f,
                solver,
            })
        },
        | _ => {
            let f = |x: f64, y: f64, z: f64| at(&[x, y, z]);
            pde::pde_solver(&pde::PdeProblem::Poisson3D {
                grid: g3,
                bc: bcs(faces),
                f: &f,
                solver,
            })
        },
    };
    Some(solution(pde_reason(solved, "poisson")?))
}

/// `advect_fd(velocity, u0, t_end, steps)`: periodic advection (Lax-Wendroff)
/// of `u0` on the unit box; `velocity` has one component per space dimension.
fn advect_fd(a: &[Data]) -> Option<Data> {
    advect(a, pde::AdvectionScheme::LaxWendroff)
}

/// `advect_fd_upwind(velocity, u0, t_end, steps)`: the same problem by the
/// first-order upwind scheme, which is monotone (no over- or undershoot at a
/// front) where Lax-Wendroff is second order and oscillates there.
fn advect_fd_upwind(a: &[Data]) -> Option<Data> {
    advect(a, pde::AdvectionScheme::Upwind)
}

fn advect(
    a: &[Data],
    scheme: pde::AdvectionScheme,
) -> Option<Data> {
    let velocity = vector(a, 0)?;
    let (shape, data) = field(a, 1)?;
    require(
        velocity.len() == shape.len(),
        Some(0),
        DiagKind::Shape,
        format_args!(
            "the velocity needs one component per space dimension ({}), got {}",
            shape.len(),
            velocity.len()
        ),
    )?;
    let t_end = scalar(a, 2)?;
    let steps = steps_within(a, 3, data.len())?;
    let (g1, g2, g3) = unit_grids(shape);
    let at = sampler(shape, data, true);
    let solved = match shape.len() {
        | 1 => {
            let u0 = |x: f64| at(&[x]);
            pde::pde_solver(&pde::PdeProblem::Advection1D {
                a: velocity[0],
                grid: g1,
                scheme,
                u0: &u0,
                t_end,
                steps,
            })
        },
        | 2 => {
            let u0 = |x: f64, y: f64| at(&[x, y]);
            pde::pde_solver(&pde::PdeProblem::Advection2D {
                a: [velocity[0], velocity[1]],
                grid: g2,
                scheme,
                u0: &u0,
                t_end,
                steps,
            })
        },
        | _ => {
            let u0 = |x: f64, y: f64, z: f64| at(&[x, y, z]);
            pde::pde_solver(&pde::PdeProblem::Advection3D {
                a: [velocity[0], velocity[1], velocity[2]],
                grid: g3,
                scheme,
                u0: &u0,
                t_end,
                steps,
            })
        },
    };
    Some(solution(pde_reason(solved, "advection")?))
}

// ---------------------------------------------------------------------------
// Krylov solvers on dense matrices
// ---------------------------------------------------------------------------

/// Largest system the dense Krylov operators accept.
pub(super) const MAX_SYSTEM: usize = 1_500;

/// A square matrix and a right-hand side of the same order.
pub(super) fn system(a: &[Data]) -> Option<(usize, &[f64], &[f64])> {
    let Some(Data::Array { shape, data }) = a.first() else {
        return bad(
            0,
            DiagKind::Shape,
            format_args!("expected a square matrix, got {}", Shown(a.first())),
        );
    };
    let b = vector(a, 1)?;
    let n = b.len();
    require(
        n > 0 && n <= MAX_SYSTEM,
        Some(1),
        DiagKind::Domain,
        format_args!("the system order must lie in 1..={MAX_SYSTEM}, got {n}"),
    )?;
    require(
        shape.as_slice() == [n, n],
        Some(0),
        DiagKind::Shape,
        format_args!(
            "the matrix must be {n} x {n} to match the right-hand side, got {}",
            Shown(a.first())
        ),
    )?;
    require(
        data.iter().all(|x| x.is_finite()),
        Some(0),
        DiagKind::Domain,
        format_args!("the matrix holds a number that is not finite"),
    )?;
    Some((n, data.as_slice(), b))
}

pub(super) fn dense_apply(
    n: usize,
    m: &[f64],
) -> impl Fn(&[f64], &mut [f64]) + '_ {
    move |x: &[f64], y: &mut [f64]| {
        for (i, out) in y.iter_mut().enumerate().take(n) {
            *out = m[i * n..(i + 1) * n]
                .iter()
                .zip(x)
                .map(|(a, b)| a * b)
                .sum();
        }
    }
}

pub(super) const fn identity_precondition(
    r: &[f64],
    z: &mut [f64],
) {
    z.copy_from_slice(r);
}

/// The solution of an iteration that reached `tol`; anything else is not an answer.
pub(super) fn converged(
    result: krylov::IterResult,
    tol: f64,
) -> Option<Data> {
    if !result.converged {
        return bad_op(
            DiagKind::NotConverged,
            format_args!(
                "the iteration reached a relative residual of {:e} after {} steps, short of the \
                 tolerance {tol:e} (the matrix may be singular, indefinite for this method, or \
                 badly conditioned)",
                result.residual, result.iterations
            ),
        );
    }
    require(
        result.x.iter().all(|v| v.is_finite()),
        None,
        DiagKind::Domain,
        format_args!("the solution is not finite"),
    )?;
    Some(Data::vector(result.x))
}

pub(super) fn tolerance(a: &[Data]) -> Option<f64> {
    let tol = scalar(a, 2)?;
    require(
        tol > 0.0 && tol < 1.0,
        Some(2),
        DiagKind::Domain,
        format_args!("the tolerance must lie strictly between 0 and 1, got {tol}"),
    )?;
    Some(tol)
}

/// `minres_solve(A, b, tol)`: symmetric, possibly indefinite `A`.
fn minres_solve(a: &[Data]) -> Option<Data> {
    let (n, m, b) = system(a)?;
    let tol = tolerance(a)?;
    converged(
        krylov_ext::minres(
            dense_apply(n, m),
            identity_precondition,
            b,
            None,
            tol,
            10 * n + 100,
        ),
        tol,
    )
}

/// `fgmres_solve(A, b, tol)`: restarted flexible GMRES for a general `A`.
fn fgmres_solve(a: &[Data]) -> Option<Data> {
    let (n, m, b) = system(a)?;
    let tol = tolerance(a)?;
    converged(
        krylov_ext::fgmres(
            dense_apply(n, m),
            identity_precondition,
            b,
            None,
            n.min(50),
            tol,
            10 * n + 100,
        ),
        tol,
    )
}

/// `bicgstab_l_solve(A, b, tol)`: BiCGSTAB(2) for a general `A`.
fn bicgstab_l_solve(a: &[Data]) -> Option<Data> {
    let (n, m, b) = system(a)?;
    let tol = tolerance(a)?;
    converged(
        krylov_ext::bicgstab_l(
            dense_apply(n, m),
            identity_precondition,
            2,
            b,
            None,
            tol,
            10 * n + 100,
        ),
        tol,
    )
}

// ---------------------------------------------------------------------------
// Tensors
// ---------------------------------------------------------------------------

/// Most elements a tensor operator may produce.
pub(super) const MAX_TENSOR: usize = 1_000_000;

pub(super) fn tensor_arg(
    a: &[Data],
    at: usize,
) -> Option<ndarray::ArrayD<f64>> {
    match a.get(at) {
        | Some(Data::Array { shape, data }) => {
            require(
                data.iter().all(|x| x.is_finite()),
                Some(at),
                DiagKind::Domain,
                format_args!("the tensor holds a number that is not finite"),
            )?;
            ndarray::ArrayD::from_shape_vec(ndarray::IxDyn(shape), data.clone())
                .reason(DiagKind::Shape, Some(at))
        },
        | other => {
            bad(
                at,
                DiagKind::Shape,
                format_args!("expected a tensor (a nested list), got {}", Shown(other)),
            )
        },
    }
}

pub(super) fn from_tensor(t: &ndarray::ArrayD<f64>) -> Option<Data> {
    if t.ndim() == 0 {
        return t.iter().next().map(|&x| Data::Scalar(x));
    }
    Some(Data::Array {
        shape: t.shape().to_vec(),
        data: t.iter().copied().collect(),
    })
}

/// An axis list: a vector of small non-negative integers.
fn axes(
    a: &[Data],
    at: usize,
) -> Option<Vec<usize>> {
    let list = vector(a, at)?;
    require(
        list.iter()
            .all(|&x| x >= 0.0 && x.fract() == 0.0 && x < 16.0),
        Some(at),
        DiagKind::Domain,
        format_args!("axes are integers from 0 to 15, got {list:?}"),
    )?;
    Some(list.iter().map(|&x| x as usize).collect())
}

/// Both operands of a product are small enough that the result is bounded.
fn bounded_product(
    x: &ndarray::ArrayD<f64>,
    y: &ndarray::ArrayD<f64>,
) -> Option<()> {
    require(
        x.len()
            .checked_mul(y.len())
            .is_some_and(|n| n <= MAX_TENSOR),
        None,
        DiagKind::Domain,
        format_args!(
            "the operands ({} and {} elements) could give a result of more than {MAX_TENSOR} \
             elements",
            x.len(),
            y.len()
        ),
    )
}

fn tensor_norm(a: &[Data]) -> Option<Data> {
    num(tensor::norm(&tensor_arg(a, 0)?))
}

fn tensor_inner(a: &[Data]) -> Option<Data> {
    let (x, y) = (tensor_arg(a, 0)?, tensor_arg(a, 1)?);
    num(tensor::inner_product(&x, &y).reason(DiagKind::Shape, Some(1))?)
}

fn tensor_outer(a: &[Data]) -> Option<Data> {
    let (x, y) = (tensor_arg(a, 0)?, tensor_arg(a, 1)?);
    bounded_product(&x, &y)?;
    from_tensor(&tensor::outer_product(&x, &y).reason(DiagKind::Domain, None)?)
}

/// `tensordot(x, y, axes_x, axes_y)`: the sum over the paired axes, whose
/// extents must agree; an axis past the rank or named twice is rejected.
fn tensordot(a: &[Data]) -> Option<Data> {
    let (x, y) = (tensor_arg(a, 0)?, tensor_arg(a, 1)?);
    let (ax, ay) = (axes(a, 2)?, axes(a, 3)?);
    require(
        ax.len() == ay.len(),
        Some(3),
        DiagKind::Shape,
        format_args!(
            "{} axes of the first tensor are paired with {} of the second",
            ax.len(),
            ay.len()
        ),
    )?;
    for (list, t, at, which) in [(&ax, &x, 2, "first"), (&ay, &y, 3, "second")] {
        let rank = t.ndim();
        require(
            list.iter().all(|&axis| axis < rank),
            Some(at),
            DiagKind::Domain,
            format_args!("an axis is past the rank of the {which} tensor (rank {rank}): {list:?}"),
        )?;
        let mut sorted = list.clone();
        sorted.sort_unstable();
        sorted.dedup();
        require(
            sorted.len() == list.len(),
            Some(at),
            DiagKind::Domain,
            format_args!("an axis of the {which} tensor is named twice: {list:?}"),
        )?;
    }
    for (&p, &q) in ax.iter().zip(&ay) {
        require(
            x.shape()[p] == y.shape()[q],
            Some(2),
            DiagKind::Shape,
            format_args!(
                "the paired axes {p} (extent {}) and {q} (extent {}) must have equal extents",
                x.shape()[p],
                y.shape()[q]
            ),
        )?;
    }
    // The result has at most `len(x) * len(y)` elements.
    bounded_product(&x, &y)?;
    from_tensor(&tensor::tensordot(&x, &y, &ax, &ay).reason(DiagKind::Shape, None)?)
}

/// `tensor_contract(T, i, j)`: the generalised trace over axes `i` and `j`,
/// a tensor of two modes fewer (a scalar for a matrix).
fn tensor_contract(a: &[Data]) -> Option<Data> {
    let t = tensor_arg(a, 0)?;
    let (i, j) = (count(a, 1, 15)?, count(a, 2, 15)?);
    require(
        i < t.ndim() && j < t.ndim(),
        Some(1),
        DiagKind::Domain,
        format_args!(
            "the axes {i} and {j} must both be below the rank {} of the tensor",
            t.ndim()
        ),
    )?;
    require(
        i != j,
        Some(2),
        DiagKind::Domain,
        format_args!("a trace needs two different axes, got {i} twice"),
    )?;
    require(
        t.shape()[i] == t.shape()[j],
        Some(1),
        DiagKind::Shape,
        format_args!(
            "the traced axes {i} and {j} must have equal extents, got {} and {}",
            t.shape()[i],
            t.shape()[j]
        ),
    )?;
    from_tensor(&tensor::contract(&t, i, j).reason(DiagKind::Shape, None)?)
}

// ---------------------------------------------------------------------------
// Functions as sampled points
// ---------------------------------------------------------------------------

/// A function given as `list(list(x, y), ...)`, ascending in `x`.
fn points(
    a: &[Data],
    at: usize,
) -> Option<Vec<(f64, f64)>> {
    match a.get(at) {
        | Some(Data::Array { shape, data }) if shape.len() == 2 && shape[1] == 2 => {
            require(
                shape[0] >= 2 && shape[0] <= MAX_LEN,
                Some(at),
                DiagKind::Domain,
                format_args!(
                    "a sampled function has 2 to {MAX_LEN} points, got {}",
                    shape[0]
                ),
            )?;
            require(
                data.iter().all(|x| x.is_finite()),
                Some(at),
                DiagKind::Domain,
                format_args!("the points hold a number that is not finite"),
            )?;
            Some(data.chunks(2).map(|p| (p[0], p[1])).collect())
        },
        | other => {
            bad(
                at,
                DiagKind::Shape,
                format_args!("expected a list of (x, y) pairs, got {}", Shown(other)),
            )
        },
    }
}

fn from_points(p: &[(f64, f64)]) -> Data {
    Data::Array {
        shape: vec![p.len(), 2],
        data: p.iter().flat_map(|&(x, y)| [x, y]).collect(),
    }
}

fn function_l1_norm(a: &[Data]) -> Option<Data> {
    num(functional_analysis::l1_norm(&points(a, 0)?))
}

fn function_l2_norm(a: &[Data]) -> Option<Data> {
    num(functional_analysis::l2_norm(&points(a, 0)?))
}

fn function_inf_norm(a: &[Data]) -> Option<Data> {
    num(functional_analysis::infinity_norm(&points(a, 0)?))
}

fn function_inner(a: &[Data]) -> Option<Data> {
    let (f, g) = (points(a, 0)?, points(a, 1)?);
    num(functional_analysis::inner_product(&f, &g).reason(DiagKind::Domain, Some(1))?)
}

fn function_project(a: &[Data]) -> Option<Data> {
    let (f, g) = (points(a, 0)?, points(a, 1)?);
    Some(from_points(
        &functional_analysis::project(&f, &g).reason(DiagKind::Domain, Some(1))?,
    ))
}

fn function_normalize(a: &[Data]) -> Option<Data> {
    Some(from_points(&functional_analysis::normalize(&points(a, 0)?)))
}

// ---------------------------------------------------------------------------
// Low-rank decompositions of dense matrices
// ---------------------------------------------------------------------------

/// A rectangular matrix of at most `MAX_SYSTEM` rows and columns.
pub(super) fn matrix(
    a: &[Data],
    at: usize,
) -> Option<(usize, usize, &[f64])> {
    matrix_of(a.get(at), at)
}

/// A matrix argument already picked out (`at` is the position blamed).
pub(super) fn matrix_of(
    item: Option<&Data>,
    at: usize,
) -> Option<(usize, usize, &[f64])> {
    match item {
        | Some(Data::Array { shape, data }) if shape.len() == 2 => {
            require(
                shape[0] > 0 && shape[1] > 0 && shape[0] <= MAX_SYSTEM && shape[1] <= MAX_SYSTEM,
                Some(at),
                DiagKind::Domain,
                format_args!(
                    "a matrix has 1 to {MAX_SYSTEM} rows and columns, got {} x {}",
                    shape[0], shape[1]
                ),
            )?;
            require(
                data.iter().all(|x| x.is_finite()),
                Some(at),
                DiagKind::Domain,
                format_args!("the matrix holds a number that is not finite"),
            )?;
            Some((shape[0], shape[1], data))
        },
        | other => {
            bad(
                at,
                DiagKind::Shape,
                format_args!("expected a matrix (a list of rows), got {}", Shown(other)),
            )
        },
    }
}

/// `randomized_singular_values(A, k)`: the `k` largest singular values by
/// randomised range finding (fixed seed, so the answer is reproducible).
fn randomized_singular_values(a: &[Data]) -> Option<Data> {
    let (m, n, data) = matrix(a, 0)?;
    let k = count(a, 1, m.min(n))?;
    let forward = |x: &[f64], y: &mut [f64]| {
        for (i, out) in y.iter_mut().enumerate().take(m) {
            *out = data[i * n..(i + 1) * n]
                .iter()
                .zip(x)
                .map(|(p, q)| p * q)
                .sum();
        }
    };
    let backward = |x: &[f64], y: &mut [f64]| {
        for v in y.iter_mut() {
            *v = 0.0;
        }
        for (i, xi) in x.iter().enumerate().take(m) {
            for (out, p) in y.iter_mut().zip(&data[i * n..(i + 1) * n]) {
                *out += p * xi;
            }
        }
    };
    let svd = lowrank::randomized_svd(
        &forward,
        &backward,
        m,
        n,
        k,
        &lowrank::RandomizedOptions::default(),
    )
    .reason(DiagKind::Domain, None)?;
    vec_out(svd.s.into_iter().take(k).collect())
}

/// `lobpcg_eigenvalues(A, k, tol)`: the `k` smallest eigenvalues of a
/// symmetric matrix (`3 k <= n`).
fn lobpcg_eigenvalues(a: &[Data]) -> Option<Data> {
    let (m, n, data) = matrix(a, 0)?;
    require(
        m == n,
        Some(0),
        DiagKind::Shape,
        format_args!("the matrix must be square, got {m} x {n}"),
    )?;
    let k = count(a, 1, n / 3)?;
    let tol = tolerance(a)?;
    let symmetric = (0..n).all(|i| {
        (0..i).all(|j| {
            (data[i * n + j] - data[j * n + i]).abs() <= 1e-12 * (1.0 + data[i * n + j].abs())
        })
    });
    require(
        symmetric,
        Some(0),
        DiagKind::Domain,
        format_args!("the matrix must be symmetric"),
    )?;
    let apply = dense_apply(n, data);
    let pairs = lowrank::lobpcg(&apply, None, n, k, tol, 20 * n + 200, 0x5EED)
        .reason(DiagKind::NotConverged, None)?;
    vec_out(pairs.values)
}

// ---------------------------------------------------------------------------
// The table
// ---------------------------------------------------------------------------

/// How an operator reads its arguments, for the signature overlay and the
/// arity: `s` scalar, `v` vector, `n` nested list (rank 1 to 3), `w` symbol.
pub(super) struct Entry {
    pub(super) name: &'static str,
    pub(super) params: &'static str,
    pub(super) eval: fn(&[Data]) -> Option<Data>,
}

/// Every operator of the bridge.
pub(super) const ENTRIES: &[Entry] = &[
    Entry {
        name: "mad",
        params: "v",
        eval: mad,
    },
    Entry {
        name: "mad_sigma",
        params: "v",
        eval: mad_sigma,
    },
    Entry {
        name: "trimmed_mean",
        params: "vs",
        eval: trimmed_mean,
    },
    Entry {
        name: "winsorized_mean",
        params: "vs",
        eval: winsorized_mean,
    },
    Entry {
        name: "huber_location",
        params: "vs",
        eval: huber_location,
    },
    Entry {
        name: "biweight_location",
        params: "vs",
        eval: biweight_location,
    },
    Entry {
        name: "ecdf",
        params: "vs",
        eval: ecdf,
    },
    Entry {
        name: "kde",
        params: "vss",
        eval: kde,
    },
    Entry {
        name: "bandwidth_silverman",
        params: "v",
        eval: bandwidth_silverman,
    },
    Entry {
        name: "bandwidth_scott",
        params: "v",
        eval: bandwidth_scott,
    },
    Entry {
        name: "hann_window",
        params: "s",
        eval: hann_window,
    },
    Entry {
        name: "hamming_window",
        params: "s",
        eval: hamming_window,
    },
    Entry {
        name: "blackman_window",
        params: "s",
        eval: blackman_window,
    },
    Entry {
        name: "cross_correlation",
        params: "vv",
        eval: cross_correlation,
    },
    Entry {
        name: "fft_abs",
        params: "v",
        eval: fft_abs,
    },
    Entry {
        name: "fft_re",
        params: "v",
        eval: fft_re,
    },
    Entry {
        name: "fft_im",
        params: "v",
        eval: fft_im,
    },
    Entry {
        name: "ifft_re",
        params: "vv",
        eval: ifft_re,
    },
    Entry {
        name: "ifft_im",
        params: "vv",
        eval: ifft_im,
    },
    Entry {
        name: "ifft_abs",
        params: "vv",
        eval: ifft_abs,
    },
    Entry {
        name: "dwt_approx",
        params: "vw",
        eval: dwt_approx,
    },
    Entry {
        name: "dwt_detail",
        params: "vw",
        eval: dwt_detail,
    },
    Entry {
        name: "idwt",
        params: "vvw",
        eval: idwt,
    },
    Entry {
        name: "heat_fd",
        params: "snsss",
        eval: heat_fd,
    },
    Entry {
        name: "heat_fd_bc",
        params: "snwsss",
        eval: heat_fd_bc,
    },
    Entry {
        name: "wave_fd",
        params: "snnsss",
        eval: wave_fd,
    },
    Entry {
        name: "wave_fd_bc",
        params: "snnwsss",
        eval: wave_fd_bc,
    },
    Entry {
        name: "poisson_fd",
        params: "ns",
        eval: poisson_fd,
    },
    Entry {
        name: "poisson_fd_bc",
        params: "nws",
        eval: poisson_fd_bc,
    },
    Entry {
        name: "advect_fd",
        params: "vnss",
        eval: advect_fd,
    },
    Entry {
        name: "advect_fd_upwind",
        params: "vnss",
        eval: advect_fd_upwind,
    },
    Entry {
        name: "minres_solve",
        params: "vvs",
        eval: minres_solve,
    },
    Entry {
        name: "fgmres_solve",
        params: "vvs",
        eval: fgmres_solve,
    },
    Entry {
        name: "bicgstab_l_solve",
        params: "vvs",
        eval: bicgstab_l_solve,
    },
    Entry {
        name: "randomized_singular_values",
        params: "vs",
        eval: randomized_singular_values,
    },
    Entry {
        name: "lobpcg_eigenvalues",
        params: "vss",
        eval: lobpcg_eigenvalues,
    },
    Entry {
        name: "tensor_norm",
        params: "n",
        eval: tensor_norm,
    },
    Entry {
        name: "tensor_inner",
        params: "nn",
        eval: tensor_inner,
    },
    Entry {
        name: "tensor_outer",
        params: "nn",
        eval: tensor_outer,
    },
    Entry {
        name: "tensordot",
        params: "nnvv",
        eval: tensordot,
    },
    Entry {
        name: "tensor_contract",
        params: "nss",
        eval: tensor_contract,
    },
    Entry {
        name: "function_l1_norm",
        params: "n",
        eval: function_l1_norm,
    },
    Entry {
        name: "function_l2_norm",
        params: "n",
        eval: function_l2_norm,
    },
    Entry {
        name: "function_inf_norm",
        params: "n",
        eval: function_inf_norm,
    },
    Entry {
        name: "function_inner",
        params: "nn",
        eval: function_inner,
    },
    Entry {
        name: "function_project",
        params: "nn",
        eval: function_project,
    },
    Entry {
        name: "function_normalize",
        params: "n",
        eval: function_normalize,
    },
];

/// Whether every number of a result is finite, the condition for encoding it.
fn finite(data: &Data) -> bool {
    match data {
        | Data::Scalar(x) => x.is_finite(),
        | Data::Array { data, .. } => data.iter().all(|x| x.is_finite()),
        | Data::Name(_) => false,
        | Data::Tree(parts) => parts.iter().all(finite),
    }
}

/// Why the call `node` of the operator `entry` does not reduce, or `None` when
/// it does (or when the operator is not a bridged one).
///
/// The arguments are decoded as the kernel does and the evaluator is run again
/// under a [`diag::capture`], so the answer is the reason the evaluator itself
/// recorded; nothing is added to the graph.
pub(super) fn explain_entry(
    entry: &Entry,
    graph: &Graph,
    node: NodeId,
) -> Option<diag::Diag> {
    let name = entry.name;
    let children = graph.children(node);
    if children.len() != entry.params.len() {
        return Some(
            diag::Diag::new(
                DiagKind::Arity,
                format!(
                    "takes {} arguments, got {}",
                    entry.params.len(),
                    children.len()
                ),
            )
            .in_op(name),
        );
    }
    let mut args = Vec::with_capacity(children.len());
    for (index, &child) in children.iter().enumerate() {
        let Some(value) = decode(graph, child) else {
            return Some(
                diag::Diag::new(
                    DiagKind::Symbolic,
                    "the argument is not a literal number, rectangular list or symbol (it is an \
                     expression, or a ragged or too deeply nested list); only literals are \
                     computed",
                )
                .at(index)
                .in_op(name),
            );
        };
        args.push(value);
    }
    let (value, found) = diag::capture(|| (entry.eval)(&args));
    match value {
        | Some(out) if finite(&out) => None,
        | Some(_) => {
            Some(diag::Diag::new(DiagKind::Domain, "the result is not finite").in_op(name))
        },
        | None => {
            Some(
                found
                    .unwrap_or_else(|| {
                        diag::Diag::new(
                            DiagKind::Unspecified,
                            "the arguments were declined and the operator gave no reason",
                        )
                    })
                    .in_op(name),
            )
        },
    }
}

/// The reduction function of one bridged operator.
pub(super) type BridgeFn = fn(&[Data]) -> Option<Data>;

/// The kernel reducing every operator of [`ENTRIES`].
pub(super) struct Bridge {
    pub(super) ops: Vec<(OpId, BridgeFn)>,
}

impl Kernel for Bridge {
    fn ops(&self) -> Vec<OpId> {
        self.ops.iter().map(|&(op, _)| op).collect()
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let op = cx.graph.op(node);
        let Some(&(_, eval)) = self.ops.iter().find(|&&(o, _)| o == op) else {
            return Outcome::Pass;
        };
        let args: Option<Vec<Data>> = cx
            .graph
            .children(node)
            .to_vec()
            .into_iter()
            .map(|child| decode(cx.graph, child))
            .collect();
        let Some(args) = args else {
            return Outcome::Pass;
        };
        match eval(&args).and_then(|out| encode(cx.graph, &out)) {
            | Some(term) => Outcome::Equal(term),
            | None => Outcome::Pass,
        }
    }
}
