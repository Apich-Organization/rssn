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
//! | signal | `hann_window`, `hamming_window`, `blackman_window`, `cross_correlation`, `fft_abs`, `fft_re`, `fft_im`, `ifft_re`, `ifft_im` |
//! | wavelets | `dwt_approx`, `dwt_detail`, `idwt` (the wavelet is a symbol: `haar`, `db2` .. `db10`, `sym4` .. `sym6`, `bior1_1` .. `bior4_4`) |
//! | finite-difference PDEs | `heat_fd`, `wave_fd`, `poisson_fd`, `advect_fd` |
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
    match args.get(at)? {
        | Data::Scalar(x) if x.is_finite() => Some(*x),
        | _ => None,
    }
}

pub(super) fn vector(
    args: &[Data],
    at: usize,
) -> Option<&[f64]> {
    match args.get(at)? {
        | Data::Array { shape, data }
            if shape.len() == 1 && data.len() <= MAX_LEN && data.iter().all(|x| x.is_finite()) =>
        {
            Some(data)
        },
        | _ => None,
    }
}

/// A non-empty finite sample.
pub(super) fn sample(
    args: &[Data],
    at: usize,
) -> Option<&[f64]> {
    vector(args, at).filter(|v| !v.is_empty())
}

pub(super) fn count(
    args: &[Data],
    at: usize,
    max: usize,
) -> Option<usize> {
    let x = scalar(args, at)?;
    (x >= 0.0 && x.fract() == 0.0 && x <= max as f64).then_some(x as usize)
}

pub(super) fn name(
    args: &[Data],
    at: usize,
) -> Option<&str> {
    match args.get(at)? {
        | Data::Name(n) => Some(n),
        | _ => None,
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
    (k > 0.0).then_some(())?;
    num(robust::huber_location(sample(a, 0)?, k))
}

fn biweight_location(a: &[Data]) -> Option<Data> {
    let c = scalar(a, 1)?;
    (c > 0.0).then_some(())?;
    num(robust::biweight_location(sample(a, 0)?, c))
}

fn ecdf(a: &[Data]) -> Option<Data> {
    num(robust::ecdf(sample(a, 0)?, scalar(a, 1)?))
}

fn kde(a: &[Data]) -> Option<Data> {
    let h = scalar(a, 2)?;
    (h > 0.0).then_some(())?;
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
    if re.len() != im.len() {
        return None;
    }
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
    wavelets::Wavelet::ALL
        .into_iter()
        .find(|w| format!("{w:?}").to_ascii_lowercase() == wanted)
}

fn dwt_approx(a: &[Data]) -> Option<Data> {
    let (x, w) = (vector(a, 0)?, wavelet(a, 1)?);
    vec_out(wavelets::dwt(x, w).ok()?.0)
}

fn dwt_detail(a: &[Data]) -> Option<Data> {
    let (x, w) = (vector(a, 0)?, wavelet(a, 1)?);
    vec_out(wavelets::dwt(x, w).ok()?.1)
}

fn idwt(a: &[Data]) -> Option<Data> {
    let (approx, detail, w) = (vector(a, 0)?, vector(a, 1)?, wavelet(a, 2)?);
    vec_out(wavelets::idwt(approx, detail, w).ok()?)
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
    match a.get(at)? {
        | Data::Array { shape, data }
            if (1..=3).contains(&shape.len())
                && shape.iter().all(|&n| n >= 3)
                && data.len() <= MAX_NODES[shape.len() - 1]
                && data.iter().all(|x| x.is_finite()) =>
        {
            Some((shape, data))
        },
        | _ => None,
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

fn bcs<const N: usize>(value: f64) -> [pde::Boundary; N] {
    std::array::from_fn(|_| pde::Boundary::dirichlet_const(value))
}

fn steps_within(
    steps: f64,
    nodes: usize,
) -> Option<usize> {
    let steps = if steps >= 1.0 && steps.fract() == 0.0 && steps <= 1e6 {
        steps as usize
    } else {
        return None;
    };
    (steps.checked_mul(nodes)? <= MAX_WORK).then_some(steps)
}

/// `heat_fd(alpha, u0, boundary_value, t_end, steps)`.
fn heat_fd(a: &[Data]) -> Option<Data> {
    let (alpha, bv, t_end) = (scalar(a, 0)?, scalar(a, 2)?, scalar(a, 3)?);
    let (shape, data) = field(a, 1)?;
    let steps = steps_within(scalar(a, 4)?, data.len())?;
    let (g1, g2, g3) = unit_grids(shape);
    let at = sampler(shape, data, false);
    let problem = match shape.len() {
        | 1 => {
            let u0 = |x: f64| at(&[x]);
            return Some(solution(
                pde::pde_solver(&pde::PdeProblem::Heat1D {
                    alpha,
                    grid: g1,
                    bc: bcs(bv),
                    u0: &u0,
                    t_end,
                    steps,
                })
                .ok()?,
            ));
        },
        | 2 => {
            let u0 = |x: f64, y: f64| at(&[x, y]);
            return Some(solution(
                pde::pde_solver(&pde::PdeProblem::Heat2D {
                    alpha,
                    grid: g2,
                    bc: bcs(bv),
                    u0: &u0,
                    t_end,
                    steps,
                })
                .ok()?,
            ));
        },
        | _ => {
            let u0 = |x: f64, y: f64, z: f64| at(&[x, y, z]);
            pde::pde_solver(&pde::PdeProblem::Heat3D {
                alpha,
                grid: g3,
                bc: bcs(bv),
                u0: &u0,
                t_end,
                steps,
            })
        },
    };
    Some(solution(problem.ok()?))
}

/// `wave_fd(c, u0, v0, boundary_value, t_end, steps)`.
fn wave_fd(a: &[Data]) -> Option<Data> {
    let (c, bv, t_end) = (scalar(a, 0)?, scalar(a, 3)?, scalar(a, 4)?);
    let (shape, data) = field(a, 1)?;
    let (vshape, vdata) = field(a, 2)?;
    (vshape == shape).then_some(())?;
    let steps = steps_within(scalar(a, 5)?, data.len())?;
    let (g1, g2, g3) = unit_grids(shape);
    let (at, vt) = (sampler(shape, data, false), sampler(shape, vdata, false));
    let solved = match shape.len() {
        | 1 => {
            let (u0, v0) = (|x: f64| at(&[x]), |x: f64| vt(&[x]));
            pde::pde_solver(&pde::PdeProblem::Wave1D {
                c,
                grid: g1,
                bc: bcs(bv),
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
                bc: bcs(bv),
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
                bc: bcs(bv),
                u0: &u0,
                v0: &v0,
                t_end,
                steps,
            })
        },
    };
    Some(solution(solved.ok()?))
}

/// `poisson_fd(f, boundary_value)`: `-Laplace(u) = f` with `u = boundary_value`
/// on the whole boundary of the unit box.
fn poisson_fd(a: &[Data]) -> Option<Data> {
    let bv = scalar(a, 1)?;
    let (shape, data) = field(a, 0)?;
    let (g1, g2, g3) = unit_grids(shape);
    let at = sampler(shape, data, false);
    let solver = pde::PoissonSolver::SparseLu;
    let solved = match shape.len() {
        | 1 => {
            let f = |x: f64| at(&[x]);
            pde::pde_solver(&pde::PdeProblem::Poisson1D {
                grid: g1,
                bc: bcs(bv),
                f: &f,
                solver,
            })
        },
        | 2 => {
            let f = |x: f64, y: f64| at(&[x, y]);
            pde::pde_solver(&pde::PdeProblem::Poisson2D {
                grid: g2,
                bc: bcs(bv),
                f: &f,
                solver,
            })
        },
        | _ => {
            let f = |x: f64, y: f64, z: f64| at(&[x, y, z]);
            pde::pde_solver(&pde::PdeProblem::Poisson3D {
                grid: g3,
                bc: bcs(bv),
                f: &f,
                solver,
            })
        },
    };
    Some(solution(solved.ok()?))
}

/// `advect_fd(velocity, u0, t_end, steps)`: periodic advection (Lax-Wendroff)
/// of `u0` on the unit box; `velocity` has one component per space dimension.
fn advect_fd(a: &[Data]) -> Option<Data> {
    let velocity = vector(a, 0)?;
    let (shape, data) = field(a, 1)?;
    (velocity.len() == shape.len()).then_some(())?;
    let t_end = scalar(a, 2)?;
    let steps = steps_within(scalar(a, 3)?, data.len())?;
    let (g1, g2, g3) = unit_grids(shape);
    let at = sampler(shape, data, true);
    let scheme = pde::AdvectionScheme::LaxWendroff;
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
    Some(solution(solved.ok()?))
}

// ---------------------------------------------------------------------------
// Krylov solvers on dense matrices
// ---------------------------------------------------------------------------

/// Largest system the dense Krylov operators accept.
pub(super) const MAX_SYSTEM: usize = 1_500;

/// A square matrix and a right-hand side of the same order.
pub(super) fn system(a: &[Data]) -> Option<(usize, &[f64], &[f64])> {
    let Data::Array { shape, data } = a.first()? else {
        return None;
    };
    let b = vector(a, 1)?;
    let n = b.len();
    (n > 0 && n <= MAX_SYSTEM && shape.as_slice() == [n, n] && data.iter().all(|x| x.is_finite()))
        .then_some((n, data.as_slice(), b))
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
pub(super) fn converged(result: krylov::IterResult) -> Option<Data> {
    (result.converged && result.x.iter().all(|v| v.is_finite())).then(|| Data::vector(result.x))
}

pub(super) fn tolerance(a: &[Data]) -> Option<f64> {
    let tol = scalar(a, 2)?;
    (tol > 0.0 && tol < 1.0).then_some(tol)
}

/// `minres_solve(A, b, tol)`: symmetric, possibly indefinite `A`.
fn minres_solve(a: &[Data]) -> Option<Data> {
    let (n, m, b) = system(a)?;
    converged(krylov_ext::minres(
        dense_apply(n, m),
        identity_precondition,
        b,
        None,
        tolerance(a)?,
        10 * n + 100,
    ))
}

/// `fgmres_solve(A, b, tol)`: restarted flexible GMRES for a general `A`.
fn fgmres_solve(a: &[Data]) -> Option<Data> {
    let (n, m, b) = system(a)?;
    converged(krylov_ext::fgmres(
        dense_apply(n, m),
        identity_precondition,
        b,
        None,
        n.min(50),
        tolerance(a)?,
        10 * n + 100,
    ))
}

/// `bicgstab_l_solve(A, b, tol)`: BiCGSTAB(2) for a general `A`.
fn bicgstab_l_solve(a: &[Data]) -> Option<Data> {
    let (n, m, b) = system(a)?;
    converged(krylov_ext::bicgstab_l(
        dense_apply(n, m),
        identity_precondition,
        2,
        b,
        None,
        tolerance(a)?,
        10 * n + 100,
    ))
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
    match a.get(at)? {
        | Data::Array { shape, data } if data.iter().all(|x| x.is_finite()) => {
            ndarray::ArrayD::from_shape_vec(ndarray::IxDyn(shape), data.clone()).ok()
        },
        | _ => None,
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
    vector(a, at)?
        .iter()
        .map(|&x| (x >= 0.0 && x.fract() == 0.0 && x < 16.0).then_some(x as usize))
        .collect()
}

fn tensor_norm(a: &[Data]) -> Option<Data> {
    num(tensor::norm(&tensor_arg(a, 0)?))
}

fn tensor_inner(a: &[Data]) -> Option<Data> {
    num(tensor::inner_product(&tensor_arg(a, 0)?, &tensor_arg(a, 1)?).ok()?)
}

fn tensor_outer(a: &[Data]) -> Option<Data> {
    let (x, y) = (tensor_arg(a, 0)?, tensor_arg(a, 1)?);
    (x.len().checked_mul(y.len())? <= MAX_TENSOR).then_some(())?;
    from_tensor(&tensor::outer_product(&x, &y).ok()?)
}

fn tensordot(a: &[Data]) -> Option<Data> {
    let (x, y) = (tensor_arg(a, 0)?, tensor_arg(a, 1)?);
    let (ax, ay) = (axes(a, 2)?, axes(a, 3)?);
    // The result has at most `len(x) * len(y)` elements.
    (x.len().checked_mul(y.len())? <= MAX_TENSOR).then_some(())?;
    from_tensor(&tensor::tensordot(&x, &y, &ax, &ay).ok()?)
}

// ---------------------------------------------------------------------------
// Functions as sampled points
// ---------------------------------------------------------------------------

/// A function given as `list(list(x, y), ...)`, ascending in `x`.
fn points(
    a: &[Data],
    at: usize,
) -> Option<Vec<(f64, f64)>> {
    match a.get(at)? {
        | Data::Array { shape, data }
            if shape.len() == 2
                && shape[1] == 2
                && shape[0] >= 2
                && shape[0] <= MAX_LEN
                && data.iter().all(|x| x.is_finite()) =>
        {
            Some(data.chunks(2).map(|p| (p[0], p[1])).collect())
        },
        | _ => None,
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
    num(functional_analysis::inner_product(&points(a, 0)?, &points(a, 1)?).ok()?)
}

fn function_project(a: &[Data]) -> Option<Data> {
    Some(from_points(
        &functional_analysis::project(&points(a, 0)?, &points(a, 1)?).ok()?,
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
    match a.get(at)? {
        | Data::Array { shape, data }
            if shape.len() == 2
                && shape[0] > 0
                && shape[1] > 0
                && shape[0] <= MAX_SYSTEM
                && shape[1] <= MAX_SYSTEM
                && data.iter().all(|x| x.is_finite()) =>
        {
            Some((shape[0], shape[1], data))
        },
        | _ => None,
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
    .ok()?;
    vec_out(svd.s.into_iter().take(k).collect())
}

/// `lobpcg_eigenvalues(A, k, tol)`: the `k` smallest eigenvalues of a
/// symmetric matrix (`3 k <= n`).
fn lobpcg_eigenvalues(a: &[Data]) -> Option<Data> {
    let (m, n, data) = matrix(a, 0)?;
    (m == n).then_some(())?;
    let k = count(a, 1, n / 3)?;
    let tol = tolerance(&[a[0].clone(), a[1].clone(), a.get(2)?.clone()])?;
    let symmetric = (0..n).all(|i| {
        (0..i).all(|j| {
            (data[i * n + j] - data[j * n + i]).abs() <= 1e-12 * (1.0 + data[i * n + j].abs())
        })
    });
    symmetric.then_some(())?;
    let apply = dense_apply(n, data);
    let pairs = lowrank::lobpcg(&apply, None, n, k, tol, 20 * n + 200, 0x5EED).ok()?;
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
        name: "wave_fd",
        params: "snnsss",
        eval: wave_fd,
    },
    Entry {
        name: "poisson_fd",
        params: "ns",
        eval: poisson_fd,
    },
    Entry {
        name: "advect_fd",
        params: "vnss",
        eval: advect_fd,
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
