//! Session-owned exact chain complexes, simplicial complexes and presentations.
//! Kernels own every mathematical construction; adapters validate finite work,
//! preserve radius dimensions, and never deserialize native object internals.
use super::codecs::Arg;
use super::codecs::Output;
use super::numeric;
use super::runtime_dimensions;
use crate::graph::Arity;
use crate::graph::OpDescriptor;
use crate::graph::RuleError;
use crate::graph::op::EvalDyn;
use crate::graph::rule::Installer;
use crate::kernel::call::CallCtx;
use crate::kernel::call::OpError;
use crate::kernel::native_objects::NativeObjects;
use crate::kernel::signature::ArityInfo;
use crate::kernel::signature::Determinism;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::Effect;
use crate::kernel::signature::EffectRow;
use crate::kernel::signature::NatExpr;
use crate::kernel::signature::Purity;
use crate::kernel::signature::Signature;
use crate::kernel::signature::Stability;
use crate::kernel::signature::ThreadSafety;
use crate::kernel::signature::TypeExpr;
use crate::kernel::signature::TypeParam;
use crate::kernel::signature::TypeParamKind;
use crate::kernel::value::DynValue;
use crate::kernel::value::RecordValue;
use crate::kernel::value::TensorValue;
use crate::kernel::value::Value;
use crate::kernel::value::VariantValue;
use crate::kernels::homology::ChainComplex;
use crate::kernels::homology::Presentation;
use crate::kernels::homology::{
    self,
};
use crate::kernels::topology::Simplex;
use num_rational::BigRational;
use num_traits::ToPrimitive;
use std::sync::Arc;

const COMPLEX: &str = "rssn/native/homology/complex/v1";
const CHAINS: &str = "rssn/native/homology/chains/v1";
const PRESENTATION: &str = "rssn/native/homology/presentation/v1";
const FILTRATION: &str = "rssn/native/homology/filtration/v1";
const LIMIT: usize = 4096;
fn error(text: &str) -> OpError {
    OpError::Domain(text.into())
}
fn bank<'a>(ctx: &'a CallCtx<'_>) -> Result<&'a NativeObjects, OpError> {
    ctx.native_objects
        .as_ref()
        .ok_or_else(|| error("homology handles require a session call context"))
}
fn record(fields: Vec<(&str, Value)>) -> Value {
    Value::Record(RecordValue {
        fields: fields
            .into_iter()
            .map(|(name, value)| (name.into(), value))
            .collect(),
    })
}
fn field<'a>(
    value: &'a Value,
    name: &str,
    names: &[&str],
) -> Result<&'a Value, OpError> {
    let Value::Record(r) = value else {
        return Err(error("expected a homology record"));
    };
    if r.fields.len() != names.len()
        || names
            .iter()
            .any(|name| r.fields.iter().filter(|(key, _)| key == name).count() != 1)
    {
        return Err(error(
            "homology record has missing, duplicate or extra fields",
        ));
    }
    r.get(name).ok_or_else(|| error("missing homology field"))
}
fn integer(value: &Value) -> Result<usize, OpError> {
    let Value::Int(n) = value else {
        return Err(error("homology count must be an exact nonnegative integer"));
    };
    n.to_usize()
        .ok_or_else(|| error("homology count does not fit usize"))
}
fn count(
    value: &Value,
    max: usize,
) -> Result<usize, OpError> {
    integer(value).and_then(|n| {
        if n <= max {
            Ok(n)
        } else {
            Err(error("homology work parameter exceeds its budget"))
        }
    })
}
fn raw_tensor(
    value: &Value,
    rank: usize,
) -> Result<&TensorValue, OpError> {
    let Value::Tensor(t) = value else {
        return Err(error("homology expected a tensor"));
    };
    if t.shape.len() != rank
        || t.shape.iter().any(|n| *n > 65_536)
        || t.shape.iter().try_fold(1usize, |a, b| a.checked_mul(*b)) != Some(t.data.len())
        || t.data.len() > 65_536
    {
        return Err(error(
            "homology tensor shape or allocation budget is invalid",
        ));
    }
    Ok(t)
}
fn canonical(dimension: &DimExpr) -> Result<DimExpr, OpError> {
    runtime_dimensions::expression(&runtime_dimensions::canonical(dimension)?)
}
fn tensor_dimension(value: &Value) -> Result<DimExpr, OpError> {
    let tensor = raw_tensor(value, 1)?;
    match &tensor.elem {
        | TypeExpr::Dim(_, dim) => canonical(dim),
        | _ => Ok(DimExpr::Dimless),
    }
}
fn measured(
    value: f64,
    dimension: &DimExpr,
) -> Value {
    if *dimension == DimExpr::Dimless {
        Value::Float(value)
    } else {
        Value::Dyn(DynValue {
            value: Box::new(Value::Float(value)),
            dim: dimension.clone(),
        })
    }
}
fn tensor_type(
    rank: u64,
    elem: TypeExpr,
) -> TypeExpr {
    TypeExpr::Tensor {
        rank: NatExpr::Lit(rank),
        shape: None,
        elem: Box::new(elem),
    }
}
fn option_type(ty: TypeExpr) -> TypeExpr {
    TypeExpr::Variant(vec![("None".into(), vec![]), ("Some".into(), vec![ty])])
}
fn option(value: Option<Value>) -> Value {
    Value::Variant(match value {
        | Some(value) => {
            VariantValue {
                tag: "Some".into(),
                payload: vec![value],
            }
        },
        | None => {
            VariantValue {
                tag: "None".into(),
                payload: vec![],
            }
        },
    })
}
fn measured_vector(
    values: Vec<f64>,
    dimension: &DimExpr,
) -> Result<Value, OpError> {
    let mut result = values.encode()?;
    if let Value::Tensor(tensor) = &mut result {
        tensor.elem = TypeExpr::Float.with_dim(dimension.clone());
        tensor.data = tensor
            .data
            .drain(..)
            .map(|value| {
                let Value::Float(v) = value else {
                    unreachable!()
                };
                measured(v, dimension)
            })
            .collect();
    }
    Ok(result)
}
fn group_type() -> TypeExpr {
    TypeExpr::Record(vec![
        ("rank".into(), TypeExpr::Int),
        ("torsion".into(), tensor_type(1, TypeExpr::Int)),
    ])
}
fn group(value: homology::Group) -> Result<Value, OpError> {
    Ok(record(vec![
        ("rank", Value::int(value.rank)),
        ("torsion", value.torsion.encode()?),
    ]))
}
#[allow(clippy::indexing_slicing, clippy::arithmetic_side_effects)] // Validated finite simplex/count bounds make these storage calculations safe.
fn validate_complex(complex: &[Simplex]) -> Result<usize, OpError> {
    if complex.len() > LIMIT {
        return Err(error("simplicial complex exceeds 4096 simplex budget"));
    }
    let mut faces = 0usize;
    let mut words = 0usize;
    for simplex in complex {
        if simplex.len() > 12 || simplex.iter().any(|v| *v > 1_000_000) {
            return Err(error("simplex exceeds dimension or vertex-label budget"));
        }
        faces = faces
            .checked_add(1usize << simplex.len())
            .filter(|n| *n <= 65_536)
            .ok_or_else(|| error("simplicial closure exceeds finite face budget"))?;
        words = words
            .checked_add(simplex.len())
            .ok_or_else(|| error("simplex storage overflows"))?;
    }
    Ok(words * std::mem::size_of::<usize>() + complex.len() * 24)
}
fn store_complex(
    ctx: &CallCtx<'_>,
    complex: Vec<Simplex>,
) -> Result<Value, OpError> {
    let bytes = validate_complex(&complex)?;
    bank(ctx)?.insert(COMPLEX, vec![], complex, bytes)
}
fn read_complex(
    ctx: &CallCtx<'_>,
    value: &Value,
) -> Result<Vec<Simplex>, OpError> {
    bank(ctx)?.with::<Vec<Simplex>, _>(value, COMPLEX, |complex| Ok(complex.clone()))
}
#[allow(clippy::indexing_slicing, clippy::arithmetic_side_effects)] // Validated finite simplex/count bounds make these storage calculations safe.
fn chain_budget(complex: &[Simplex]) -> Result<(), OpError> {
    let mut sizes = [0usize; 13];
    for simplex in complex {
        if !simplex.is_empty() {
            sizes[simplex.len() - 1] += 1;
        }
    }
    let mut work = 0usize;
    for size in sizes {
        if size > 256 {
            return Err(error(
                "chain-group rank exceeds 256 exact elimination budget",
            ));
        }
        work = work
            .checked_add(size * size * size)
            .filter(|n| *n <= 50_000_000)
            .ok_or_else(|| error("chain exact-elimination work exceeds budget"))?;
    }
    Ok(())
}
#[allow(clippy::indexing_slicing, clippy::arithmetic_side_effects)] // Validated finite simplex/count bounds make these storage calculations safe.
fn store_chains(
    ctx: &CallCtx<'_>,
    chains: ChainComplex,
) -> Result<Value, OpError> {
    let entries = chains.d.iter().flatten().map(Vec::len).sum::<usize>();
    if entries > 65_536 {
        return Err(error("chain boundary matrices exceed storage budget"));
    }
    bank(ctx)?.insert(CHAINS, vec![], chains, entries * 8 + 4096)
}
#[allow(clippy::indexing_slicing, clippy::arithmetic_side_effects)] // Validated finite simplex/count bounds make these storage calculations safe.
fn store_presentation(
    ctx: &CallCtx<'_>,
    value: Presentation,
) -> Result<Value, OpError> {
    let words = value.relators.iter().map(Vec::len).sum::<usize>();
    if value.count > 256
        || words > 4000
        || value
            .relators
            .iter()
            .flatten()
            .any(|v| *v == 0 || v.unsigned_abs() > value.count as u64)
    {
        return Err(error(
            "presentation exceeds generator/word budget or contains invalid generator indices",
        ));
    }
    bank(ctx)?.insert(PRESENTATION, vec![], value, words * 8 + 8192)
}
fn filtration_type(dimension: DimExpr) -> TypeExpr {
    TypeExpr::Record(vec![
        ("payload".into(), TypeExpr::Opaque(FILTRATION.into())),
        ("value_unit".into(), TypeExpr::Float.with_dim(dimension)),
    ])
}
fn bars_type(dimension: DimExpr) -> TypeExpr {
    TypeExpr::Record(vec![
        ("dimensions".into(), tensor_type(1, TypeExpr::Int)),
        (
            "births".into(),
            tensor_type(1, TypeExpr::Float.with_dim(dimension.clone())),
        ),
        (
            "deaths".into(),
            tensor_type(1, TypeExpr::Float.with_dim(dimension)),
        ),
        ("essential".into(), tensor_type(1, TypeExpr::Bool)),
    ])
}
fn interval_type(dimension: DimExpr) -> TypeExpr {
    TypeExpr::Record(vec![
        (
            "births".into(),
            tensor_type(1, TypeExpr::Float.with_dim(dimension.clone())),
        ),
        (
            "deaths".into(),
            tensor_type(1, TypeExpr::Float.with_dim(dimension)),
        ),
        ("essential".into(), tensor_type(1, TypeExpr::Bool)),
    ])
}
#[allow(clippy::indexing_slicing, clippy::arithmetic_side_effects)] // Validated finite simplex/count bounds make these storage calculations safe.
fn store_filtration(
    ctx: &CallCtx<'_>,
    value: Vec<(f64, Simplex)>,
    dim: &DimExpr,
) -> Result<Value, OpError> {
    validate_complex(&value.iter().map(|(_, s)| s.clone()).collect::<Vec<_>>())?;
    if value.len() > LIMIT || value.iter().any(|(v, _)| !v.is_finite()) {
        return Err(error("filtration exceeds finite radius/simplex budget"));
    }
    let charge = value.iter().map(|(_, s)| s.len() * 8 + 40).sum();
    let handle = bank(ctx)?.insert(FILTRATION, vec![dim.clone()], value, charge)?;
    Ok(record(vec![
        ("payload", handle),
        ("value_unit", measured(1.0, dim)),
    ]))
}
fn read_filtration(
    ctx: &CallCtx<'_>,
    value: &Value,
) -> Result<(Vec<(f64, Simplex)>, DimExpr), OpError> {
    let names = ["payload", "value_unit"];
    let handle = field(value, "payload", &names)?;
    let witness = field(value, "value_unit", &names)?;
    if !crate::kernels::exact_eq(numeric(witness, 0)?, 1.0) {
        return Err(error("filtration unit witness must be one coherent unit"));
    }
    let dimensions = bank(ctx)?.dimensions(handle)?;
    let dim = dimensions
        .first()
        .ok_or_else(|| error("filtration lost its radius dimension"))?
        .clone();
    let witness_dim = match witness {
        | Value::Dyn(v) => canonical(&v.dim)?,
        | Value::Quantity(q) => canonical(&q.unit.dimension)?,
        | _ => DimExpr::Dimless,
    };
    if canonical(&dim)? != witness_dim {
        return Err(error(
            "filtration unit witness conflicts with its native object",
        ));
    }
    Ok((
        bank(ctx)?.with::<Vec<(f64, Simplex)>, _>(handle, FILTRATION, |v| Ok(v.clone()))?,
        dim,
    ))
}
fn encode_bars(
    bars: Vec<homology::Bar>,
    dim: &DimExpr,
) -> Result<Value, OpError> {
    let mut dimensions = Vec::new();
    let mut births = Vec::new();
    let mut deaths = Vec::new();
    let mut essential = Vec::new();
    for bar in bars {
        dimensions.push(bar.dim);
        births.push(bar.birth);
        let infinite = crate::kernels::exact_eq(bar.death, f64::INFINITY);
        deaths.push(if infinite { 0.0 } else { bar.death });
        essential.push(infinite);
    }
    Ok(record(vec![
        ("dimensions", dimensions.encode()?),
        ("births", measured_vector(births, dim)?),
        ("deaths", measured_vector(deaths, dim)?),
        ("essential", essential.encode()?),
    ]))
}

#[derive(Clone, Copy)]
enum Method {
    Build,
    Export,
    Chain,
    ComplexChains,
    Relative,
    Reduced,
    Homology,
    Cohomology,
    ModP,
    BettiModP,
    TopDegree,
    Full,
    Sphere,
    Cone,
    Suspension,
    Join,
    Wedge,
    Product,
    Klein,
    Dunce,
    Projective,
    Union,
    Intersection,
    Subcomplex,
    Betti,
    Mayer,
    Pi1,
    Simplify,
    Abelianization,
    PresentationExport,
    CubicalBetti,
    CubicalChains,
    Smith,
    RankModP,
    IsPrime,
    MatrixRank,
    Distance,
    Rips,
    Filtration,
    Bars,
    Bottleneck,
}
struct Adapter(Method);
use Method::Abelianization;
use Method::Bars;
use Method::Betti;
use Method::BettiModP;
use Method::Bottleneck;
use Method::Build;
use Method::Chain;
use Method::Cohomology;
use Method::ComplexChains;
use Method::Cone;
use Method::CubicalBetti;
use Method::CubicalChains;
use Method::Distance;
use Method::Dunce;
use Method::Export;
use Method::Filtration;
use Method::Full;
use Method::Homology;
use Method::Intersection;
use Method::IsPrime;
use Method::Join;
use Method::Klein;
use Method::MatrixRank;
use Method::Mayer;
use Method::ModP;
use Method::Pi1;
use Method::PresentationExport;
use Method::Product;
use Method::Projective;
use Method::RankModP;
use Method::Reduced;
use Method::Relative;
use Method::Rips;
use Method::Simplify;
use Method::Smith;
use Method::Sphere;
use Method::Subcomplex;
use Method::Suspension;
use Method::TopDegree;
use Method::Union;
use Method::Wedge;
impl EvalDyn for Adapter {
    #[allow(clippy::too_many_lines)] // One signature table keeps the native homology object API consistent.
    fn signature(&self) -> Signature {
        let complex = || TypeExpr::Opaque(COMPLEX.into());
        let chains = || TypeExpr::Opaque(CHAINS.into());
        let presentation = || TypeExpr::Opaque(PRESENTATION.into());
        let d = DimExpr::Var("D".into());
        let int = || TypeExpr::Int;
        let matrix = || tensor_type(2, TypeExpr::Int);
        let csr = || <Vec<Vec<usize>> as Arg>::ty();
        let (params, ret) = match self.0 {
            | Build => (vec![csr()], complex()),
            | Export => (vec![complex()], csr()),
            | Chain => (vec![tensor_type(1, complex())], chains()),
            | ComplexChains => (vec![complex()], chains()),
            | Relative => (vec![complex(), complex()], chains()),
            | Reduced | Homology | Cohomology => (vec![chains(), int()], group_type()),
            | ModP => (vec![chains(), int(), int()], int()),
            | BettiModP => (vec![chains(), int()], tensor_type(1, int())),
            | TopDegree => (vec![chains()], option_type(int())),
            | Full | Sphere | Dunce | Projective => (vec![int()], complex()),
            | Klein => (vec![int(), int()], complex()),
            | Cone | Suspension => (vec![complex()], complex()),
            | Join | Union | Intersection => (vec![complex(), complex()], complex()),
            | Wedge => (vec![complex(), complex(), int(), int()], complex()),
            | Product => (vec![complex(), complex(), int()], option_type(complex())),
            | Subcomplex => (vec![complex(), complex()], TypeExpr::Bool),
            | Betti => (vec![complex()], tensor_type(1, int())),
            | Mayer => {
                (
                    vec![complex(), complex()],
                    option_type(tensor_type(1, int())),
                )
            },
            | Pi1 => (vec![complex()], option_type(presentation())),
            | Simplify => (vec![presentation()], presentation()),
            | Abelianization => (vec![presentation()], group_type()),
            | PresentationExport => {
                (
                    vec![presentation()],
                    TypeExpr::Record(vec![
                        ("count".into(), int()),
                        (
                            "relators".into(),
                            TypeExpr::Record(vec![
                                ("offsets".into(), tensor_type(1, int())),
                                ("values".into(), tensor_type(1, int())),
                            ]),
                        ),
                        ("edges".into(), tensor_type(2, int())),
                    ]),
                )
            },
            | CubicalBetti => {
                (
                    vec![tensor_type(2, TypeExpr::Bool), TypeExpr::Bool],
                    TypeExpr::Record(vec![("item0".into(), int()), ("item1".into(), int())]),
                )
            },
            | CubicalChains => {
                (
                    vec![tensor_type(2, TypeExpr::Bool), TypeExpr::Bool],
                    chains(),
                )
            },
            | Smith => (vec![matrix()], tensor_type(1, int())),
            | RankModP => (vec![matrix(), int()], int()),
            | IsPrime => (vec![int()], TypeExpr::Bool),
            | MatrixRank => (vec![tensor_type(2, TypeExpr::Real.with_dim(d))], int()),
            | Distance => {
                (
                    vec![tensor_type(2, TypeExpr::Float.with_dim(d.clone()))],
                    tensor_type(2, TypeExpr::Float.with_dim(d)),
                )
            },
            | Rips => {
                (
                    vec![
                        tensor_type(2, TypeExpr::Float.with_dim(d.clone())),
                        TypeExpr::Float.with_dim(d.clone()),
                        int(),
                        int(),
                    ],
                    option_type(filtration_type(d)),
                )
            },
            | Filtration => {
                (
                    vec![TypeExpr::Record(vec![
                        (
                            "radii".into(),
                            tensor_type(1, TypeExpr::Float.with_dim(d.clone())),
                        ),
                        ("simplices".into(), csr()),
                    ])],
                    filtration_type(d),
                )
            },
            | Bars => {
                (
                    vec![filtration_type(d.clone()), int()],
                    option_type(bars_type(d)),
                )
            },
            | Bottleneck => {
                (
                    vec![interval_type(d.clone()), interval_type(d.clone())],
                    TypeExpr::Float.with_dim(d),
                )
            },
        };
        let mut sig = Signature::new(params, ret);
        if matches!(
            self.0,
            Distance | Rips | Filtration | Bars | Bottleneck | MatrixRank
        ) {
            sig.type_params.push(TypeParam {
                name: "D".into(),
                kind: TypeParamKind::Dim,
                bounds: vec![],
            });
        }
        sig.effects = EffectRow {
            effects: vec![Effect::Stateful, Effect::Diverges],
        };
        sig
    }

    #[allow(
        clippy::too_many_lines,
        clippy::indexing_slicing,
        clippy::arithmetic_side_effects
    )] // Exact arity, object kinds and finite shape/work budgets precede native calls.
    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let arity = self.signature().params.len();
        if args.len() != arity {
            return Err(OpError::Arity {
                expected: ArityInfo::fixed(arity as u8),
                got: args.len(),
            });
        }
        match self.0 {
            | Build => {
                let raw = Vec::<Vec<usize>>::decode(&args[0], 0)?;
                validate_complex(&raw)?;
                store_complex(ctx, crate::kernels::topology::close_complex(&raw))
            },
            | Export => read_complex(ctx, &args[0])?.encode(),
            | Full | Sphere | Dunce | Projective | Klein => {
                let n = count(
                    &args[0],
                    if matches!(self.0, Full | Sphere) {
                        7
                    } else {
                        64
                    },
                )?;
                let result = match self.0 {
                    | Full => homology::full_simplex(n),
                    | Sphere => homology::sphere(n),
                    | Dunce => {
                        if n < 3 {
                            return Err(error(
                                "dunce-cap boundary subdivision must be at least three",
                            ));
                        }
                        homology::dunce_cap(n)
                    },
                    | Projective => {
                        if n != 0 && n < 3 {
                            return Err(error(
                                "projective-plane subdivision must be zero or at least three",
                            ));
                        }
                        homology::projective_plane(n)
                    },
                    | _ => {
                        let m = count(&args[1], 64)?;
                        if n < 3 || m < 3 || n * m > 128 {
                            return Err(error(
                                "Klein triangulation needs subdivisions at least three within 128-cell budget",
                            ));
                        }
                        homology::klein_bottle(n, m)
                    },
                };
                store_complex(ctx, result)
            },
            | Cone | Suspension => {
                let complex = read_complex(ctx, &args[0])?;
                if complex.len() > 1024 {
                    return Err(error("cone/suspension input exceeds construction budget"));
                }
                store_complex(
                    ctx,
                    if matches!(self.0, Cone) {
                        homology::cone(&complex)
                    } else {
                        homology::suspension(&complex)
                    },
                )
            },
            | Union | Intersection | Join | Wedge | Product | Subcomplex | Mayer => {
                let a = read_complex(ctx, &args[0])?;
                let b = read_complex(ctx, &args[1])?;
                if matches!(self.0, Subcomplex) {
                    return Ok(Value::Bool(homology::is_subcomplex(&a, &b)));
                }
                if matches!(self.0, Join | Product)
                    && a.len().checked_mul(b.len()).is_none_or(|n| n > 4096)
                {
                    return Err(error("complex combination exceeds finite product budget"));
                }
                match self.0 {
                    | Union => store_complex(ctx, homology::union(&a, &b)),
                    | Intersection => store_complex(ctx, homology::intersection(&a, &b)),
                    | Join => store_complex(ctx, homology::join(&a, &b)),
                    | Wedge => {
                        let va = integer(&args[2])?;
                        let vb = integer(&args[3])?;
                        if !a.contains(&vec![va]) || !b.contains(&vec![vb]) {
                            return Err(error(
                                "wedge points must be vertices of their respective complexes",
                            ));
                        }
                        store_complex(ctx, homology::wedge(&a, &b, va, vb))
                    },
                    | Product => {
                        let limit = count(&args[2], LIMIT)?;
                        let vertices_a = a.iter().filter(|simplex| simplex.len() == 1).count();
                        let vertices_b = b.iter().filter(|simplex| simplex.len() == 1).count();
                        let width = a.iter().map(Vec::len).max().unwrap_or(0)
                            + b.iter().map(Vec::len).max().unwrap_or(0);
                        if vertices_a
                            .checked_mul(vertices_b)
                            .and_then(|n| n.checked_mul(limit))
                            .and_then(|n| n.checked_mul(width))
                            .is_none_or(|n| n > 1_048_576)
                        {
                            return Err(error(
                                "complex product frontier exceeds its checked work-array budget",
                            ));
                        }
                        let out = homology::product(&a, &b, limit);
                        Ok(option(out.map(|v| store_complex(ctx, v)).transpose()?))
                    },
                    | _ => {
                        let united = homology::union(&a, &b);
                        chain_budget(&united)?;
                        Ok(option(
                            homology::mayer_vietoris(&a, &b)
                                .map(Output::encode)
                                .transpose()?,
                        ))
                    },
                }
            },
            | Betti | Pi1 | ComplexChains => {
                let complex = read_complex(ctx, &args[0])?;
                chain_budget(&complex)?;
                match self.0 {
                    | Betti => homology::betti_list(&complex).encode(),
                    | Pi1 => {
                        Ok(option(
                            homology::pi1_presentation(&complex)
                                .map(|p| store_presentation(ctx, p))
                                .transpose()?,
                        ))
                    },
                    | _ => store_chains(ctx, homology::complex_chains(&complex)),
                }
            },
            | Chain => {
                let handles = raw_tensor(&args[0], 1)?;
                if handles.data.len() > 13 {
                    return Err(error("chain basis exceeds degree budget"));
                }
                let mut basis = Vec::new();
                for (degree, handle) in handles.data.iter().enumerate() {
                    let simplices = read_complex(ctx, handle)?
                        .into_iter()
                        .filter(|s| s.len() == degree + 1)
                        .collect::<Vec<_>>();
                    if simplices.len() > 256 {
                        return Err(error("chain basis exceeds exact elimination budget"));
                    }
                    basis.push(simplices);
                }
                store_chains(ctx, homology::chain_complex(&basis))
            },
            | Relative => {
                let complex = read_complex(ctx, &args[0])?;
                let sub = read_complex(ctx, &args[1])?;
                if !homology::is_subcomplex(&sub, &complex) {
                    return Err(error("relative chains require an actual subcomplex"));
                }
                chain_budget(&complex)?;
                store_chains(ctx, homology::relative_chains(&complex, &sub))
            },
            | Reduced | Homology | Cohomology | ModP | BettiModP | TopDegree => {
                bank(ctx)?.with::<ChainComplex, _>(&args[0], CHAINS, |chains| {
                    match self.0 {
                        | TopDegree => Ok(option(chains.top_degree().map(Value::int))),
                        | BettiModP => chains.betti_mod_p(field_prime(&args[1], 1)?).encode(),
                        | _ => {
                            let k = count(&args[1], 64)?;
                            match self.0 {
                                | Reduced => group(homology::reduced_homology(chains, k)),
                                | Homology => group(chains.homology(k)),
                                | Cohomology => group(chains.cohomology(k)),
                                | _ => {
                                    Ok(Value::int(
                                        chains.homology_mod_p(k, field_prime(&args[2], 2)?),
                                    ))
                                },
                            }
                        },
                    }
                })
            },
            | Simplify | Abelianization | PresentationExport => {
                let presentation = bank(ctx)?.with::<Presentation, _>(
                    &args[0],
                    PRESENTATION,
                    |p| Ok(p.clone()),
                )?;
                match self.0 {
                    | Simplify => {
                        if presentation.count > 16
                            || presentation.relators.iter().map(Vec::len).sum::<usize>() > 256
                        {
                            return Err(error(
                                "presentation simplification is bounded to sixteen generators and 256 input word coefficients",
                            ));
                        }
                        store_presentation(ctx, homology::simplify_presentation(&presentation))
                    },
                    | Abelianization => {
                        matrix_budget(presentation.relators.len(), presentation.count)?;
                        group(homology::abelianization(&presentation))
                    },
                    | _ => {
                        let mut offsets = vec![0usize];
                        let mut values = Vec::new();
                        for relator in presentation.relators {
                            values.extend(relator);
                            offsets.push(values.len());
                        }
                        let edges = presentation
                            .edges
                            .into_iter()
                            .map(|(a, b)| [a, b])
                            .collect::<Vec<_>>();
                        Ok(record(vec![
                            ("count", Value::int(presentation.count)),
                            (
                                "relators",
                                record(vec![
                                    ("offsets", offsets.encode()?),
                                    ("values", values.encode()?),
                                ]),
                            ),
                            ("edges", edges.encode()?),
                        ]))
                    },
                }
            },
            | CubicalBetti | CubicalChains => {
                let image = decode_image(
                    &args[0],
                    if matches!(self.0, CubicalChains) {
                        64
                    } else {
                        4096
                    },
                )?;
                let eight = bool::decode(&args[1], 1)?;
                if matches!(self.0, CubicalBetti) {
                    homology::cubical_betti(&image, eight).encode()
                } else {
                    let chains = homology::cubical_chains(&image, eight);
                    store_chains(ctx, chains)
                }
            },
            | Smith | RankModP => {
                let data = Vec::<Vec<i64>>::decode(&args[0], 0)?;
                matrix_budget(data.len(), data.first().map_or(0, Vec::len))?;
                if matches!(self.0, Smith) {
                    homology::smith_invariants(&data).encode()
                } else {
                    let prime = field_prime(&args[1], 1)?;
                    Ok(Value::int(homology::rank_mod_p(&data, prime)))
                }
            },
            | IsPrime => Ok(Value::Bool(homology::is_prime(u64::decode(&args[0], 0)?))),
            | MatrixRank => {
                let tensor = raw_tensor(&args[0], 2)?;
                matrix_budget(tensor.shape[0], tensor.shape[1])?;
                let data = (0..tensor.shape[0])
                    .map(|i| {
                        (0..tensor.shape[1])
                            .map(|j| BigRational::decode(&tensor.data[i * tensor.shape[1] + j], 0))
                            .collect::<Result<Vec<_>, _>>()
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                Ok(Value::int(homology::matrix_rank(&data)))
            },
            | Distance => {
                let t = raw_tensor(&args[0], 2)?;
                let points = Vec::<Vec<f64>>::decode(&args[0], 0)?;
                if points.len() > 256 {
                    return Err(error("point distance matrix exceeds pair budget"));
                }
                let dim = element_dimension(&t.elem)?;
                measured_matrix(homology::distance_matrix(&points), &dim)
            },
            | Rips => {
                let t = raw_tensor(&args[0], 2)?;
                let dist = Vec::<Vec<f64>>::decode(&args[0], 0)?;
                if dist.len() > 128 || t.shape[0] != t.shape[1] {
                    return Err(error(
                        "Rips distances must be a square matrix of order at most 128",
                    ));
                }
                let dim = element_dimension(&t.elem)?;
                let epsilon = numeric(&args[1], 1)?;
                if epsilon < 0.0 || scalar_dimension(&args[1])? != dim {
                    return Err(error(
                        "Rips radius must be nonnegative and share distance units",
                    ));
                }
                for (i, row) in dist.iter().enumerate() {
                    for (j, v) in row.iter().enumerate() {
                        if *v < 0.0
                            || !crate::kernels::exact_eq(*v, dist[j][i])
                            || (i == j && !crate::kernels::exact_eq(*v, 0.0))
                        {
                            return Err(error(
                                "distance matrix must be symmetric, nonnegative, and zero on its diagonal",
                            ));
                        }
                    }
                }
                let degree = count(&args[2], 8)?;
                let limit = count(&args[3], LIMIT)?;
                if limit < dist.len() {
                    return Err(error("Rips simplex limit must cover the initial vertices"));
                }
                Ok(option(
                    homology::rips_filtration(&dist, epsilon, degree, limit)
                        .map(|f| store_filtration(ctx, f, &dim))
                        .transpose()?,
                ))
            },
            | Filtration => {
                let radii = field(&args[0], "radii", &["radii", "simplices"])?;
                let simplices = field(&args[0], "simplices", &["radii", "simplices"])?;
                let dim = tensor_dimension(radii)?;
                let radii = Vec::<f64>::decode(radii, 0)?;
                let simplices = Vec::<Vec<usize>>::decode(simplices, 0)?;
                if radii.len() != simplices.len() {
                    return Err(error(
                        "filtration radii and simplices must have equal counts",
                    ));
                }
                store_filtration(ctx, radii.into_iter().zip(simplices).collect(), &dim)
            },
            | Bars => {
                let (filtration, dim) = read_filtration(ctx, &args[0])?;
                let degree = count(&args[1], 8)?;
                Ok(option(
                    homology::persistence_bars(&filtration, degree)
                        .map(|bars| encode_bars(bars, &dim))
                        .transpose()?,
                ))
            },
            | Bottleneck => {
                let (a, dim) = decode_intervals(&args[0])?;
                let (b, other) = decode_intervals(&args[1])?;
                if dim != other {
                    return Err(error("bottleneck diagrams must share endpoint dimensions"));
                }
                if a.len() + b.len() > 128 {
                    return Err(error(
                        "bottleneck matching exceeds 128-interval work budget",
                    ));
                }
                let distance = homology::bottleneck_distance(&a, &b);
                if distance.is_nan() {
                    return Err(error("native bottleneck distance returned NaN"));
                }
                Ok(measured(distance, &dim))
            },
        }
    }
}
fn element_dimension(ty: &TypeExpr) -> Result<DimExpr, OpError> {
    match ty {
        | TypeExpr::Dim(_, dim) => canonical(dim),
        | TypeExpr::Float | TypeExpr::Real | TypeExpr::Int => Ok(DimExpr::Dimless),
        | _ => Err(error("homology values must be numeric coordinates")),
    }
}
fn scalar_dimension(value: &Value) -> Result<DimExpr, OpError> {
    match value {
        | Value::Dyn(d) => canonical(&d.dim),
        | Value::Quantity(q) => canonical(&q.unit.dimension),
        | Value::Float(_) | Value::Int(_) | Value::Real(_) => Ok(DimExpr::Dimless),
        | _ => Err(error("expected a scalar coordinate")),
    }
}
fn measured_matrix(
    values: Vec<Vec<f64>>,
    dim: &DimExpr,
) -> Result<Value, OpError> {
    let mut value = values.encode()?;
    if let Value::Tensor(t) = &mut value {
        t.elem = TypeExpr::Float.with_dim(dim.clone());
        t.data = t
            .data
            .drain(..)
            .map(|v| {
                let Value::Float(v) = v else { unreachable!() };
                measured(v, dim)
            })
            .collect();
    }
    Ok(value)
}
fn field_prime(
    value: &Value,
    argument: usize,
) -> Result<u64, OpError> {
    let prime = u64::decode(value, argument)?;
    if !homology::is_prime(prime) {
        return Err(error(
            "the coefficient modulus must be a prime unsigned 64-bit integer",
        ));
    }
    Ok(prime)
}

fn matrix_budget(
    rows: usize,
    cols: usize,
) -> Result<(), OpError> {
    if rows > 128
        || cols > 128
        || rows
            .checked_mul(cols)
            .and_then(|n| n.checked_mul(rows.min(cols)))
            .is_none_or(|n| n > 2_097_152)
    {
        return Err(error("exact matrix elimination is bounded to 128 x 128"));
    }
    Ok(())
}
#[allow(clippy::indexing_slicing, clippy::arithmetic_side_effects)] // Exact checked rectangular shape bounds every row index.
fn decode_image(
    value: &Value,
    max: usize,
) -> Result<Vec<Vec<bool>>, OpError> {
    let t = raw_tensor(value, 2)?;
    if t.data.len() > max {
        return Err(error("cubical image exceeds finite topology work budget"));
    }
    let cols = t.shape[1];
    (0..t.shape[0])
        .map(|r| {
            (0..cols)
                .map(|c| bool::decode(&t.data[r * cols + c], 0))
                .collect()
        })
        .collect()
}
fn decode_intervals(value: &Value) -> Result<(Vec<(f64, f64)>, DimExpr), OpError> {
    let names = ["births", "deaths", "essential"];
    let births = field(value, "births", &names)?;
    let deaths = field(value, "deaths", &names)?;
    let flags = field(value, "essential", &names)?;
    let dim = tensor_dimension(births)?;
    if tensor_dimension(deaths)? != dim {
        return Err(error("persistence birth/death dimensions differ"));
    }
    let births = Vec::<f64>::decode(births, 0)?;
    let deaths = Vec::<f64>::decode(deaths, 0)?;
    let essential = Vec::<bool>::decode(flags, 0)?;
    if births.len() != deaths.len() || births.len() != essential.len() || births.len() > 128 {
        return Err(error(
            "persistence endpoint/mask counts disagree or exceed matching budget",
        ));
    }
    let intervals = births
        .into_iter()
        .zip(deaths)
        .zip(essential)
        .map(|((b, d), e)| {
            if !e && d < b {
                Err(error("finite persistence death precedes birth"))
            } else {
                Ok((b, if e { f64::INFINITY } else { d }))
            }
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok((intervals, dim))
}

#[allow(clippy::too_many_lines)] // Exact source provenance is declared alongside each adapter method.
pub(super) fn install(installer: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, path, method) in [
        (
            "num_homology_complex",
            "crate::kernels::topology::close_complex",
            Build,
        ),
        ("num_homology_complex_export", "", Export),
        (
            "num_homology_chain_complex",
            "crate::kernels::homology::chain_complex",
            Chain,
        ),
        (
            "num_homology_complex_chains",
            "crate::kernels::homology::complex_chains",
            ComplexChains,
        ),
        (
            "num_homology_relative_chains",
            "crate::kernels::homology::relative_chains",
            Relative,
        ),
        (
            "num_homology_reduced_homology",
            "crate::kernels::homology::reduced_homology",
            Reduced,
        ),
        (
            "num_homology_chain_homology",
            "crate::kernels::homology::ChainComplex::homology",
            Homology,
        ),
        (
            "num_homology_chain_cohomology",
            "crate::kernels::homology::ChainComplex::cohomology",
            Cohomology,
        ),
        (
            "num_homology_chain_homology_mod_p",
            "crate::kernels::homology::ChainComplex::homology_mod_p",
            ModP,
        ),
        (
            "num_homology_chain_betti_mod_p",
            "crate::kernels::homology::ChainComplex::betti_mod_p",
            BettiModP,
        ),
        (
            "num_homology_chain_top_degree",
            "crate::kernels::homology::ChainComplex::top_degree",
            TopDegree,
        ),
        (
            "num_homology_full_simplex",
            "crate::kernels::homology::full_simplex",
            Full,
        ),
        (
            "num_homology_sphere",
            "crate::kernels::homology::sphere",
            Sphere,
        ),
        ("num_homology_cone", "crate::kernels::homology::cone", Cone),
        (
            "num_homology_suspension",
            "crate::kernels::homology::suspension",
            Suspension,
        ),
        ("num_homology_join", "crate::kernels::homology::join", Join),
        (
            "num_homology_wedge",
            "crate::kernels::homology::wedge",
            Wedge,
        ),
        (
            "num_homology_product",
            "crate::kernels::homology::product",
            Product,
        ),
        (
            "num_homology_klein_bottle",
            "crate::kernels::homology::klein_bottle",
            Klein,
        ),
        (
            "num_homology_dunce_cap",
            "crate::kernels::homology::dunce_cap",
            Dunce,
        ),
        (
            "num_homology_projective_plane",
            "crate::kernels::homology::projective_plane",
            Projective,
        ),
        (
            "num_homology_union",
            "crate::kernels::homology::union",
            Union,
        ),
        (
            "num_homology_intersection",
            "crate::kernels::homology::intersection",
            Intersection,
        ),
        (
            "num_homology_is_subcomplex",
            "crate::kernels::homology::is_subcomplex",
            Subcomplex,
        ),
        (
            "num_homology_betti_list",
            "crate::kernels::homology::betti_list",
            Betti,
        ),
        (
            "num_homology_mayer_vietoris",
            "crate::kernels::homology::mayer_vietoris",
            Mayer,
        ),
        (
            "num_homology_pi1_presentation",
            "crate::kernels::homology::pi1_presentation",
            Pi1,
        ),
        (
            "num_homology_simplify_presentation",
            "crate::kernels::homology::simplify_presentation",
            Simplify,
        ),
        (
            "num_homology_abelianization",
            "crate::kernels::homology::abelianization",
            Abelianization,
        ),
        ("num_homology_presentation_export", "", PresentationExport),
        (
            "num_homology_cubical_betti",
            "crate::kernels::homology::cubical_betti",
            CubicalBetti,
        ),
        (
            "num_homology_cubical_chains",
            "crate::kernels::homology::cubical_chains",
            CubicalChains,
        ),
        (
            "num_homology_smith_invariants",
            "crate::kernels::homology::smith_invariants",
            Smith,
        ),
        (
            "num_homology_rank_mod_p",
            "crate::kernels::homology::rank_mod_p",
            RankModP,
        ),
        (
            "num_homology_is_prime",
            "crate::kernels::homology::is_prime",
            IsPrime,
        ),
        (
            "num_homology_matrix_rank",
            "crate::kernels::homology::matrix_rank",
            MatrixRank,
        ),
        (
            "num_homology_distance_matrix",
            "crate::kernels::homology::distance_matrix",
            Distance,
        ),
        (
            "num_homology_rips_filtration",
            "crate::kernels::homology::rips_filtration",
            Rips,
        ),
        ("num_homology_filtration", "", Filtration),
        (
            "num_homology_persistence_bars",
            "crate::kernels::homology::persistence_bars",
            Bars,
        ),
        (
            "num_homology_bottleneck_distance",
            "crate::kernels::homology::bottleneck_distance",
            Bottleneck,
        ),
    ] {
        let adapter = Adapter(method);
        let signature = adapter.signature();
        let mut descriptor=OpDescriptor::new(name,Arity::Fixed(signature.params.len()as u8)).category("numerical_vector")
 .doc("Native exact algebraic topology with session-owned, kind/version checked complex, chain and presentation handles. Complex inputs/exports use CSR offsets/indices. Handles are invalid in other sessions and after release; all objects own their data. Finite construction bounds: 4096 simplices, 12 vertices per simplex, 65536 face/storage words, 256 generators per chain group and 128 x 128 exact matrix elimination. Presentation simplification accepts sixteen generators and 256 input word coefficients with the native 4000-word intermediate limit; product frontiers are bounded to 1048576 coefficient slots. Cubical chains accept 64 pixels; direct cubical Betti accepts 4096. Filtrations carry a canonical value_unit witness verified against stored radius dimensions. Bars use explicit essential flags and zero death placeholders; bottleneck inputs are births/deaths/essential records and unequal essential counts give positive infinity. Products return None on their native simplex limit; malformed inputs and adapter budget exhaustion are errors. Pure source algorithms are invoked in RSSN kernels; opaque payloads contain only session handles, never serialized pointers.")
 .stability(Stability::Experimental).purity(if matches!(method,Build|Chain|ComplexChains|Relative|Full|Sphere|Cone|Suspension|Join|Wedge|Product|Klein|Dunce|Projective|Union|Intersection|Pi1|Simplify|CubicalChains|Rips|Filtration){Purity::Effectful}else{Purity::ReadOnly}).determinism(Determinism::EnvDependent).thread_safety(ThreadSafety::SessionBound)
 .backend(crate::kernel::signature::BackendSupport{interpreter:true,jit:false,ffi:false,gpu:false})
 .signature(signature).eval_dyn(Arc::new(adapter));
        let documentation = match method {
            | BettiModP => {
                Some(
                    "Return the Betti numbers of a session-owned chain complex in ascending degree order over a prime field. The second argument is an exact prime in the unsigned 64-bit range. Each boundary matrix is reduced once, with packed binary elimination for characteristic two.",
                )
            },
            | ModP => {
                Some(
                    "Return the dimension of one homology group over a prime field. Pass a session-owned chain complex, a nonnegative degree and a prime in the unsigned 64-bit range.",
                )
            },
            | RankModP => {
                Some(
                    "Return the rank of a rectangular integer matrix of at most 128 by 128 over a prime field. The modulus is an exact prime in the unsigned 64-bit range. Negative entries and characteristic two are supported.",
                )
            },
            | IsPrime => {
                Some(
                    "Determine whether an exact unsigned 64-bit integer is prime, using deterministic Miller–Rabin witnesses.",
                )
            },
            | _ => None,
        };
        if let Some(documentation) = documentation {
            descriptor = descriptor.doc(documentation);
        }
        if !path.is_empty() {
            descriptor = descriptor.rust_path(path);
        }
        let op = installer.op(descriptor)?;
        installer.dim(op, crate::graph::dim::DimRule::Container);
    }
    Ok(())
}
