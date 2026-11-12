//! Sparse matrices with exact and symbolic entries.
//!
//! A sparse matrix is the term `sparse(m, n, list(list(i, j, v), ...))`: its
//! shape and its non-zero entries, with one-based indices like every list in
//! the kernel. The constructor normalises itself (entries sorted by position,
//! repeated positions added, zeros dropped), so two equal sparse matrices are
//! the same term.
//!
//! The algorithms are those of [`crate::kernels::sparse_exact`], chosen by
//! what the entries are:
//!
//! * every entry a rational number: modular methods. `sp_det` combines
//!   determinants modulo word-sized primes up to Hadamard's bound and
//!   `sp_solve` uses Dixon's `p`-adic lifting, so no fraction grows during an
//!   elimination;
//! * anything else: the same sparsity-preserving elimination over the
//!   *domain of terms*, in which arithmetic builds expressions, brings them to
//!   a normal form with the running engine and tests for zero through it.
//!
//! | operator | value |
//! |---|---|
//! | `sparse(m, n, entries)` | the matrix; `entries` is a list of `list(i, j, v)` |
//! | `sp_from_dense(A)`, `sp_to_dense(S)` | conversion from and to a list of rows |
//! | `sp_identity(n)`, `sp_diag(list(d1, ...))` | constructors |
//! | `sp_dims(S)`, `sp_nnz(S)` | `list(m, n)`; the number of stored entries |
//! | `sp_transpose(S)`, `sp_add(S, T)`, `sp_scale(c, S)`, `sp_matmul(S, T)` | arithmetic, sparse results |
//! | `sp_matvec(S, v)` | the dense vector `S v` |
//! | `sp_rank(S)`, `sp_det(S)` | rank; determinant |
//! | `sp_solve(S, b)` | the unique solution of `S x = b`; `list()` when there is none; unreduced when there are several |
//! | `sp_solve_general(S, b)` | `list(particular, list(basis...))`, or `list()` when inconsistent |
//! | `sp_nullspace(S)` | a basis of `{x : S x = 0}` |
//! | `sp_rank_mod(S, p)`, `sp_det_mod(S, p)`, `sp_solve_mod(S, b, p)` | the same over the field of integers modulo a prime `p < 2^62` |
//! | `sp_minpoly_mod(S, p)` | the minimal polynomial over that field by Wiedemann's algorithm, as coefficients from the constant term up |
//!
//! An operator whose arguments are not of this form (a shape mismatch, a
//! composite modulus, an entry outside the declared shape) is left
//! unreduced.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::ToPrimitive;

use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::Tier;
use crate::graph::dim::DimRule;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::kernels::sparse_exact as sx;
use crate::kernels::sparse_exact::Domain;
use crate::kernels::sparse_exact::Gfp;
use crate::kernels::sparse_exact::Rationals;
use crate::kernels::sparse_exact::Solution;
use crate::kernels::sparse_exact::SparseMat;
use crate::rules::linalg;
use crate::rules::poly::best;

use super::linalg::linalg as linalg_rules;

/// The sparse-matrix rule set.
#[must_use]
pub fn sparse_symbolic() -> RuleSet {
    RuleSet::new("sparse_symbolic", install).needs(linalg_rules())
}

/// The field of simplified terms: elements are nodes, arithmetic builds
/// nodes and normalises them with the running engine.
struct Terms<'c, 'a> {
    cx: &'c mut Cx<'a>,
}

impl Domain for Terms<'_, '_> {
    type E = NodeId;

    fn zero(&mut self) -> NodeId {
        self.cx.graph.int(0)
    }

    fn one(&mut self) -> NodeId {
        self.cx.graph.int(1)
    }

    fn is_zero(
        &mut self,
        a: &NodeId,
    ) -> bool {
        if let Some(n) = self.cx.graph.number_of(*a) {
            return n.is_zero();
        }
        let value = linalg::normal(self.cx, *a);
        self.cx.graph.number_of(value).is_some_and(Number::is_zero)
    }

    fn add(
        &mut self,
        a: &NodeId,
        b: &NodeId,
    ) -> NodeId {
        let sum = linalg::add(self.cx.graph, &[*a, *b]);
        linalg::normal(self.cx, sum)
    }

    fn sub(
        &mut self,
        a: &NodeId,
        b: &NodeId,
    ) -> NodeId {
        let negated = linalg::neg(self.cx.graph, *b);
        let sum = linalg::add(self.cx.graph, &[*a, negated]);
        linalg::normal(self.cx, sum)
    }

    fn mul(
        &mut self,
        a: &NodeId,
        b: &NodeId,
    ) -> NodeId {
        let product = linalg::mul(self.cx.graph, &[*a, *b]);
        linalg::normal(self.cx, product)
    }

    fn div(
        &mut self,
        a: &NodeId,
        b: &NodeId,
    ) -> Option<NodeId> {
        if self.is_zero(b) {
            return None;
        }
        let inverse = linalg::reciprocal(self.cx.graph, *b);
        let quotient = linalg::mul(self.cx.graph, &[*a, inverse]);
        Some(linalg::normal(self.cx, quotient))
    }
}

/// What a sparse operator computes.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Kind {
    Normalize,
    FromDense,
    ToDense,
    Identity,
    Diag,
    Dims,
    Nnz,
    Transpose,
    Add,
    Scale,
    Matmul,
    Matvec,
    Rank,
    Det,
    Solve,
    SolveGeneral,
    Nullspace,
    RankMod,
    DetMod,
    SolveMod,
    MinpolyMod,
}

struct Sparse {
    /// The operator this kernel answers.
    op: OpId,
    /// The `sparse` constructor.
    ctor: OpId,
    kind: Kind,
}

/// A matrix read from a term: over the rationals when every entry is a
/// number, over terms otherwise.
enum Matrix {
    Exact(SparseMat<BigRational>),
    Symbolic(SparseMat<NodeId>),
}

fn index(
    graph: &Graph,
    node: NodeId,
) -> Option<usize> {
    let r = graph.number_of(node)?.to_rational()?;
    if !r.is_integer() {
        return None;
    }
    r.to_integer().to_usize()
}

fn rational(
    graph: &Graph,
    node: NodeId,
) -> Option<BigRational> {
    graph.number_of(node).and_then(Number::to_rational)
}

fn int_node(
    graph: &mut Graph,
    value: usize,
) -> NodeId {
    graph.num(Number::Int(BigInt::from(value)))
}

fn residue_node(
    graph: &mut Graph,
    value: u64,
) -> NodeId {
    graph.num(Number::Int(BigInt::from(value)))
}

impl Sparse {
    /// The shape and raw `(row, column, value)` entries (zero-based) of a
    /// `sparse(...)` term.
    fn raw(
        &self,
        graph: &mut Graph,
        node: NodeId,
    ) -> Option<(usize, usize, Vec<(usize, usize, NodeId)>)> {
        let class = graph.find(node);
        let term = graph
            .enodes(class)
            .filter(|&n| graph.op(n) == self.ctor)
            .min()?;
        let children = graph.children(term).to_vec();
        let [m, n, entries] = children.as_slice() else {
            return None;
        };
        let rows = index(graph, *m)?;
        let cols = index(graph, *n)?;
        let mut out = Vec::new();
        for item in linalg::vector(graph, *entries)? {
            let triple = linalg::vector(graph, item)?;
            let [i, j, v] = triple.as_slice() else {
                return None;
            };
            let (i, j) = (index(graph, *i)?, index(graph, *j)?);
            if i == 0 || j == 0 || i > rows || j > cols {
                return None;
            }
            out.push((i - 1, j - 1, best(graph, *v).unwrap_or(*v)));
        }
        Some((rows, cols, out))
    }

    fn read(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Option<Matrix> {
        let (rows, cols, entries) = self.raw(cx.graph, node)?;
        let exact: Option<Vec<(usize, usize, BigRational)>> = entries
            .iter()
            .map(|&(i, j, v)| rational(cx.graph, v).map(|q| (i, j, q)))
            .collect();
        if let Some(triplets) = exact {
            return SparseMat::from_triplets(&mut Rationals, rows, cols, triplets)
                .ok()
                .map(Matrix::Exact);
        }
        let mut terms = Terms { cx };
        SparseMat::from_triplets(&mut terms, rows, cols, entries)
            .ok()
            .map(Matrix::Symbolic)
    }

    /// A matrix with residue entries modulo `p`, or `None` when an entry is
    /// not an integer or a fraction whose denominator is invertible.
    fn read_mod(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
        field: Gfp,
    ) -> Option<SparseMat<u64>> {
        let (rows, cols, entries) = self.raw(cx.graph, node)?;
        let mut dom = field;
        let mut triplets = Vec::with_capacity(entries.len());
        for (i, j, v) in entries {
            triplets.push((i, j, residue(&rational(cx.graph, v)?, field)?));
        }
        SparseMat::from_triplets(&mut dom, rows, cols, triplets).ok()
    }

    fn build(
        &self,
        graph: &mut Graph,
        rows: usize,
        cols: usize,
        entries: Vec<(usize, usize, NodeId)>,
    ) -> NodeId {
        let items: Vec<NodeId> = entries
            .into_iter()
            .map(|(i, j, v)| {
                let i = int_node(graph, i + 1);
                let j = int_node(graph, j + 1);
                linalg::list(graph, &[i, j, v])
            })
            .collect();
        let m = int_node(graph, rows);
        let n = int_node(graph, cols);
        let list = linalg::list(graph, &items);
        graph.node(self.ctor, &[m, n, list])
    }

    fn write(
        &self,
        graph: &mut Graph,
        matrix: &Matrix,
    ) -> NodeId {
        match matrix {
            | Matrix::Exact(m) => {
                let entries = m
                    .triplets()
                    .map(|(i, j, v)| (i, j, v.clone()))
                    .collect::<Vec<_>>()
                    .into_iter()
                    .map(|(i, j, v)| (i, j, linalg::rational_node(graph, &v)))
                    .collect();
                self.build(graph, m.rows, m.cols, entries)
            },
            | Matrix::Symbolic(m) => {
                let entries = m.triplets().map(|(i, j, v)| (i, j, *v)).collect();
                self.build(graph, m.rows, m.cols, entries)
            },
        }
    }

    #[allow(clippy::too_many_lines)]
    fn compute(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
        args: &[NodeId],
    ) -> Option<NodeId> {
        match (self.kind, args) {
            | (Kind::Normalize, _) => {
                let matrix = self.read(cx, node)?;
                let canonical = self.write(cx.graph, &matrix);
                (canonical != node).then_some(canonical)
            },
            | (Kind::FromDense, [dense]) => {
                let rows = linalg::vector(cx.graph, *dense)?;
                let mut entries = Vec::new();
                let mut cols = 0;
                for (i, row) in rows.iter().enumerate() {
                    let items = linalg::vector(cx.graph, *row)?;
                    cols = cols.max(items.len());
                    for (j, v) in items.into_iter().enumerate() {
                        entries.push((i, j, v));
                    }
                }
                Some(self.build(cx.graph, rows.len(), cols, entries))
            },
            | (Kind::ToDense, [s]) => {
                let rows = match self.read(cx, *s)? {
                    | Matrix::Exact(m) => {
                        m.to_dense(&mut Rationals)
                            .iter()
                            .map(|row| {
                                row.iter()
                                    .map(|v| linalg::rational_node(cx.graph, v))
                                    .collect::<Vec<_>>()
                            })
                            .collect::<Vec<_>>()
                    },
                    | Matrix::Symbolic(m) => m.to_dense(&mut Terms { cx }),
                };
                Some(linalg::matrix_term(cx.graph, &rows))
            },
            | (Kind::Identity, [n]) => {
                let n = index(cx.graph, *n)?;
                let one = cx.graph.int(1);
                Some(self.build(cx.graph, n, n, (0..n).map(|i| (i, i, one)).collect()))
            },
            | (Kind::Diag, [d]) => {
                let items = linalg::vector(cx.graph, *d)?;
                let n = items.len();
                let entries = items
                    .into_iter()
                    .enumerate()
                    .map(|(i, v)| (i, i, v))
                    .collect();
                Some(self.build(cx.graph, n, n, entries))
            },
            | (Kind::Dims, [s]) => {
                let (rows, cols, _) = self.raw(cx.graph, *s)?;
                let m = int_node(cx.graph, rows);
                let n = int_node(cx.graph, cols);
                Some(linalg::list(cx.graph, &[m, n]))
            },
            | (Kind::Nnz, [s]) => {
                let count = match self.read(cx, *s)? {
                    | Matrix::Exact(m) => m.nnz(),
                    | Matrix::Symbolic(m) => m.nnz(),
                };
                Some(int_node(cx.graph, count))
            },
            | (Kind::Transpose, [s]) => {
                let t = match self.read(cx, *s)? {
                    | Matrix::Exact(m) => Matrix::Exact(m.transpose()),
                    | Matrix::Symbolic(m) => Matrix::Symbolic(m.transpose()),
                };
                Some(self.write(cx.graph, &t))
            },
            | (Kind::Add | Kind::Matmul, [a, b]) => {
                let left = self.read(cx, *a)?;
                let right = self.read(cx, *b)?;
                let product = self.kind == Kind::Matmul;
                let result = match (left, right) {
                    | (Matrix::Exact(x), Matrix::Exact(y)) => {
                        Matrix::Exact(if product {
                            x.matmul(&mut Rationals, &y)?
                        } else {
                            x.add(&mut Rationals, &y)?
                        })
                    },
                    | (x, y) => {
                        let x = self.symbolic(cx.graph, x);
                        let y = self.symbolic(cx.graph, y);
                        let mut terms = Terms { cx };
                        Matrix::Symbolic(if product {
                            x.matmul(&mut terms, &y)?
                        } else {
                            x.add(&mut terms, &y)?
                        })
                    },
                };
                Some(self.write(cx.graph, &result))
            },
            | (Kind::Scale, [c, s]) => {
                let c = best(cx.graph, *c)?;
                if cx.graph.op(c) == core::LIST {
                    return None;
                }
                let result = match (self.read(cx, *s)?, rational(cx.graph, c)) {
                    | (Matrix::Exact(m), Some(q)) => Matrix::Exact(m.scale(&mut Rationals, &q)),
                    | (m, _) => {
                        let m = self.symbolic(cx.graph, m);
                        Matrix::Symbolic(m.scale(&mut Terms { cx }, &c))
                    },
                };
                Some(self.write(cx.graph, &result))
            },
            | (Kind::Matvec, [s, v]) => {
                let v = linalg::vector(cx.graph, *v)?;
                let exact: Option<Vec<BigRational>> =
                    v.iter().map(|&e| rational(cx.graph, e)).collect();
                let items = match (self.read(cx, *s)?, exact) {
                    | (Matrix::Exact(m), Some(x)) => {
                        m.matvec(&mut Rationals, &x)?
                            .iter()
                            .map(|q| linalg::rational_node(cx.graph, q))
                            .collect::<Vec<_>>()
                    },
                    | (m, _) => {
                        let m = self.symbolic(cx.graph, m);
                        m.matvec(&mut Terms { cx }, &v)?
                    },
                };
                Some(linalg::list(cx.graph, &items))
            },
            | (Kind::Rank, [s]) => {
                let rank = match self.read(cx, *s)? {
                    | Matrix::Exact(m) => sx::rank(&mut Rationals, &m),
                    | Matrix::Symbolic(m) => sx::rank(&mut Terms { cx }, &m),
                };
                Some(int_node(cx.graph, rank))
            },
            | (Kind::Det, [s]) => {
                match self.read(cx, *s)? {
                    | Matrix::Exact(m) => {
                        let d = sx::det_rational(&m)?;
                        Some(linalg::rational_node(cx.graph, &d))
                    },
                    | Matrix::Symbolic(m) => sx::det(&mut Terms { cx }, &m),
                }
            },
            | (Kind::Solve | Kind::SolveGeneral, [s, b]) => self.solve(cx, *s, *b),
            | (Kind::Nullspace, [s]) => {
                let basis: Vec<Vec<NodeId>> = match self.read(cx, *s)? {
                    | Matrix::Exact(m) => {
                        sx::nullspace(&mut Rationals, &m)
                            .iter()
                            .map(|v| {
                                v.iter()
                                    .map(|q| linalg::rational_node(cx.graph, q))
                                    .collect()
                            })
                            .collect()
                    },
                    | Matrix::Symbolic(m) => sx::nullspace(&mut Terms { cx }, &m),
                };
                Some(linalg::matrix_term(cx.graph, &basis))
            },
            | (Kind::RankMod, [s, p]) => {
                let mut field = prime(cx.graph, *p)?;
                let m = self.read_mod(cx, *s, field)?;
                Some(int_node(cx.graph, sx::rank(&mut field, &m)))
            },
            | (Kind::DetMod, [s, p]) => {
                let mut field = prime(cx.graph, *p)?;
                let m = self.read_mod(cx, *s, field)?;
                let d = sx::det(&mut field, &m)?;
                Some(residue_node(cx.graph, d))
            },
            | (Kind::MinpolyMod, [s, p]) => {
                let field = prime(cx.graph, *p)?;
                let m = self.read_mod(cx, *s, field)?;
                let poly = sx::minpoly_wiedemann(&m, field, 4, 0x9E37_79B9_7F4A_7C15)?;
                let items: Vec<NodeId> = poly.iter().map(|&c| residue_node(cx.graph, c)).collect();
                Some(linalg::list(cx.graph, &items))
            },
            | (Kind::SolveMod, [s, b, p]) => {
                let mut field = prime(cx.graph, *p)?;
                let m = self.read_mod(cx, *s, field)?;
                let b: Vec<u64> = linalg::vector(cx.graph, *b)?
                    .iter()
                    .map(|&e| rational(cx.graph, e).and_then(|q| residue(&q, field)))
                    .collect::<Option<_>>()?;
                let items: Vec<NodeId> = match sx::solve(&mut field, &m, &b)? {
                    | Solution::Inconsistent => Vec::new(),
                    | Solution::Family {
                        particular,
                        nullspace,
                    } if nullspace.is_empty() => {
                        particular
                            .iter()
                            .map(|&c| residue_node(cx.graph, c))
                            .collect()
                    },
                    | Solution::Family { .. } => return None,
                };
                Some(linalg::list(cx.graph, &items))
            },
            | _ => None,
        }
    }

    /// The matrix over terms, whichever way it was read.
    fn symbolic(
        &self,
        graph: &mut Graph,
        matrix: Matrix,
    ) -> SparseMat<NodeId> {
        let _ = self;
        match matrix {
            | Matrix::Symbolic(m) => m,
            | Matrix::Exact(m) => {
                SparseMat {
                    rows: m.rows,
                    cols: m.cols,
                    data: m
                        .data
                        .iter()
                        .map(|row| {
                            row.iter()
                                .map(|(j, v)| (*j, linalg::rational_node(graph, v)))
                                .collect()
                        })
                        .collect(),
                }
            },
        }
    }

    fn solve(
        &self,
        cx: &mut Cx<'_>,
        s: NodeId,
        b: NodeId,
    ) -> Option<NodeId> {
        let b = linalg::vector(cx.graph, b)?;
        let exact: Option<Vec<BigRational>> = b.iter().map(|&e| rational(cx.graph, e)).collect();
        let solution: Solution<NodeId> = match (self.read(cx, s)?, exact) {
            | (Matrix::Exact(m), Some(rhs)) => {
                // Lifting handles the common case (square, non-singular)
                // without rational elimination; anything else is analysed by
                // the general eliminator.
                let found = match sx::solve_dixon(&m, &rhs) {
                    | Some(x) => {
                        Solution::Family {
                            particular: x,
                            nullspace: Vec::new(),
                        }
                    },
                    | None => sx::solve(&mut Rationals, &m, &rhs)?,
                };
                let node = |graph: &mut Graph, v: &[BigRational]| -> Vec<NodeId> {
                    v.iter().map(|q| linalg::rational_node(graph, q)).collect()
                };
                match found {
                    | Solution::Inconsistent => Solution::Inconsistent,
                    | Solution::Family {
                        particular,
                        nullspace,
                    } => {
                        Solution::Family {
                            particular: node(cx.graph, &particular),
                            nullspace: nullspace.iter().map(|v| node(cx.graph, v)).collect(),
                        }
                    },
                }
            },
            | (m, _) => {
                let m = self.symbolic(cx.graph, m);
                sx::solve(&mut Terms { cx }, &m, &b)?
            },
        };
        match solution {
            | Solution::Inconsistent => Some(linalg::list(cx.graph, &[])),
            | Solution::Family {
                particular,
                nullspace,
            } => {
                if self.kind == Kind::SolveGeneral {
                    let x = linalg::list(cx.graph, &particular);
                    let basis = linalg::matrix_term(cx.graph, &nullspace);
                    Some(linalg::list(cx.graph, &[x, basis]))
                } else if nullspace.is_empty() {
                    Some(linalg::list(cx.graph, &particular))
                } else {
                    // Several solutions: `sp_solve` promises the unique one.
                    None
                }
            },
        }
    }
}

/// The prime field named by a term, when the term is a prime below `2^62`.
fn prime(
    graph: &Graph,
    node: NodeId,
) -> Option<Gfp> {
    let p = graph
        .number_of(node)?
        .to_rational()
        .filter(BigRational::is_integer)?
        .to_integer()
        .to_u64()?;
    (p < (1 << 62) && sx::is_prime_u64(p)).then_some(Gfp(p))
}

/// The residue of a rational modulo `p`, or `None` when the denominator is
/// divisible by `p`.
fn residue(
    q: &BigRational,
    field: Gfp,
) -> Option<u64> {
    let den = field.inv(field.reduce(q.denom()))?;
    Some(field.mul_mod(field.reduce(q.numer()), den))
}

impl Kernel for Sparse {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        self.compute(cx, node, &args)
            .map_or(Outcome::Pass, Outcome::Equal)
    }

    fn revisit(&self) -> bool {
        // An argument may become readable only after another request reduces.
        self.kind != Kind::Normalize
    }
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    // The constructor is data, not a request: it is never "unreduced".
    let ctor = i.op(OpDescriptor::new("sparse", Arity::Fixed(3)).cost(4))?;
    i.dim(ctor, DimRule::Container);
    i.kernel(
        "sparse_symbolic/normalize",
        Tier::Normalize,
        Sparse {
            op: ctor,
            ctor,
            kind: Kind::Normalize,
        },
    );
    let table: [(&str, u8, Kind, DimRule); 20] = [
        ("sp_from_dense", 1, Kind::FromDense, DimRule::Container),
        ("sp_to_dense", 1, Kind::ToDense, DimRule::Container),
        ("sp_identity", 1, Kind::Identity, DimRule::Container),
        ("sp_diag", 1, Kind::Diag, DimRule::Container),
        ("sp_dims", 1, Kind::Dims, DimRule::NonMetric),
        ("sp_nnz", 1, Kind::Nnz, DimRule::NonMetric),
        ("sp_transpose", 1, Kind::Transpose, DimRule::Container),
        ("sp_add", 2, Kind::Add, DimRule::Container),
        ("sp_scale", 2, Kind::Scale, DimRule::Container),
        ("sp_matmul", 2, Kind::Matmul, DimRule::Container),
        ("sp_matvec", 2, Kind::Matvec, DimRule::Container),
        ("sp_rank", 1, Kind::Rank, DimRule::NonMetric),
        ("sp_det", 1, Kind::Det, DimRule::Container),
        ("sp_solve", 2, Kind::Solve, DimRule::Container),
        (
            "sp_solve_general",
            2,
            Kind::SolveGeneral,
            DimRule::Container,
        ),
        ("sp_nullspace", 1, Kind::Nullspace, DimRule::Container),
        ("sp_rank_mod", 2, Kind::RankMod, DimRule::NonMetric),
        ("sp_det_mod", 2, Kind::DetMod, DimRule::NonMetric),
        ("sp_solve_mod", 3, Kind::SolveMod, DimRule::NonMetric),
        ("sp_minpoly_mod", 2, Kind::MinpolyMod, DimRule::NonMetric),
    ];
    for (name, arity, kind, rule) in table {
        let op = i.op(OpDescriptor::new(name, Arity::Fixed(arity))
            .flags(OpFlags::HEAVY)
            .cost(100))?;
        i.dim(op, rule);
        i.kernel(
            &format!("sparse_symbolic/{name}"),
            Tier::Reduce,
            Sparse { op, ctor, kind },
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::sparse_symbolic;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[sparse_symbolic()], src)
    }

    const A: &str =
        "sparse(3, 3, list(list(1, 1, 2), list(2, 2, 3), list(3, 1, 1), list(3, 3, 4)))";

    #[test]
    fn the_constructor_is_canonical() {
        assert_eq!(
            run("sparse(2, 2, list(list(2, 1, 5), list(1, 1, 1), list(2, 1, -5), list(1, 2, 0)))"),
            "sparse(2, 2, list(list(1, 1, 1)))"
        );
        assert_eq!(
            run("sp_from_dense(list(list(0, a), list(0, 0)))"),
            "sparse(2, 2, list(list(1, 2, a)))"
        );
        assert_eq!(
            run(&format!("sp_to_dense({A})")),
            "list(list(2, 0, 0), list(0, 3, 0), list(1, 0, 4))"
        );
        assert_eq!(run(&format!("sp_dims({A})")), "list(3, 3)");
        assert_eq!(run(&format!("sp_nnz({A})")), "4");
        assert_eq!(
            run("sp_identity(2)"),
            "sparse(2, 2, list(list(1, 1, 1), list(2, 2, 1)))"
        );
        assert_eq!(
            run("sp_diag(list(x, 0, 3))"),
            "sparse(3, 3, list(list(1, 1, x), list(3, 3, 3)))"
        );
    }

    #[test]
    fn arithmetic_keeps_results_sparse() {
        assert_eq!(
            run(&format!("sp_transpose({A})")),
            "sparse(3, 3, list(list(1, 1, 2), list(1, 3, 1), list(2, 2, 3), list(3, 3, 4)))"
        );
        assert_eq!(
            run(&format!("sp_to_dense(sp_matmul({A}, sp_identity(3)))")),
            "list(list(2, 0, 0), list(0, 3, 0), list(1, 0, 4))"
        );
        assert_eq!(run(&format!("sp_nnz(sp_add({A}, sp_scale(-1, {A})))")), "0");
        assert_eq!(
            run(&format!("sp_matvec({A}, list(1, 1, 1))")),
            "list(2, 3, 5)"
        );
        assert_eq!(
            run("sp_matvec(sparse(2, 2, list(list(1, 2, a), list(2, 1, b))), list(x, y))"),
            "list(a*y, b*x)"
        );
    }

    #[test]
    fn exact_matrices_use_the_modular_algorithms() {
        assert_eq!(run(&format!("sp_det({A})")), "24");
        assert_eq!(run(&format!("sp_rank({A})")), "3");
        assert_eq!(
            run(&format!("sp_solve({A}, list(2, 3, 5))")),
            "list(1, 1, 1)"
        );
        assert_eq!(
            run("sp_solve(sparse(2, 2, list(list(1, 1, 1/2), list(2, 2, 3))), list(1, 1))"),
            "list(2, 1/3)"
        );
        // Singular: no solution, or a family.
        let singular = "sparse(2, 2, list(list(1, 1, 1), list(2, 1, 1)))";
        assert_eq!(run(&format!("sp_det({singular})")), "0");
        assert_eq!(run(&format!("sp_solve({singular}, list(1, 2))")), "list()");
        assert_eq!(
            run(&format!("sp_solve_general({singular}, list(3, 3))")),
            "list(list(3, 0), list(list(0, 1)))"
        );
        assert_eq!(
            run(&format!("sp_nullspace({singular})")),
            "list(list(0, 1))"
        );
        // Several solutions: `sp_solve` does not pick one.
        let (_, reduced) = reduce_with(
            &[sparse_symbolic()],
            &format!("sp_solve({singular}, list(3, 3))"),
            &[],
        );
        assert!(!reduced);
    }

    #[test]
    fn symbolic_entries_use_the_same_elimination() {
        assert_eq!(
            run("sp_det(sparse(2, 2, list(list(1, 1, a), list(2, 2, b))))"),
            "a*b"
        );
        assert_eq!(
            run("sp_det(sparse(2, 2, list(list(1, 2, a), list(2, 1, b))))"),
            "-a*b"
        );
        assert_eq!(
            run(
                "sp_rank(sparse(2, 2, list(list(1, 1, x), list(1, 2, x^2), list(2, 1, 1), list(2, 2, x))))"
            ),
            "1"
        );
        assert_eq!(
            run("sp_solve(sparse(2, 2, list(list(1, 1, a), list(2, 2, 2))), list(a^2, 4))"),
            "list(a, 2)"
        );
    }

    #[test]
    fn prime_field_operators() {
        assert_eq!(run(&format!("sp_det_mod({A}, 7)")), "3");
        assert_eq!(run(&format!("sp_rank_mod({A}, 2)")), "2");
        assert_eq!(
            run(&format!("sp_solve_mod({A}, list(2, 3, 5), 101)")),
            "list(1, 1, 1)"
        );
        // The identity has minimal polynomial x - 1.
        assert_eq!(run("sp_minpoly_mod(sp_identity(4), 101)"), "list(100, 1)");
        // A composite modulus is not a field.
        let (_, reduced) = reduce_with(&[sparse_symbolic()], &format!("sp_det_mod({A}, 6)"), &[]);
        assert!(!reduced);
    }

    #[test]
    fn malformed_matrices_are_left_alone() {
        for src in [
            "sp_nnz(sparse(2, 2, list(list(3, 1, 1))))",
            "sp_nnz(sparse(2, 2, list(list(0, 1, 1))))",
            "sp_det(sparse(2, 3, list(list(1, 1, 1))))",
            "sp_matmul(sp_identity(2), sp_identity(3))",
            "sp_nnz(7)",
        ] {
            let (_, reduced) = reduce_with(&[sparse_symbolic()], src, &[]);
            assert!(!reduced, "{src}");
        }
    }
}
