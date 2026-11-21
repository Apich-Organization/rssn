// Explicit portable adapters for opaque Rust models. Descriptors store defining
// data, never pointers or unchecked internal tree nodes. Queries reconstruct the
// native model; this is intentionally a portable format, not a persistent cache.
const MODEL_LIMIT: usize = 1_000_000;
fn model_error(error: impl std::fmt::Display) -> OpError {
    OpError::Domain(error.to_string())
}
fn model_entries(rows: &[Vec<f64>]) -> Result<usize, OpError> {
    let count = rows
        .iter()
        .try_fold(1usize, |n, row| n.checked_add(row.len())?.checked_add(1))
        .ok_or_else(|| model_error("model size overflows"))?;
    if count > MODEL_LIMIT {
        return Err(model_error("model exceeds portable allocation limit"));
    }
    Ok(count)
}
fn pack_model(rows: &[Vec<f64>]) -> Result<Vec<u8>, OpError> {
    let count = model_entries(rows)?;
    let mut bytes = Vec::with_capacity(count * 8);
    bytes.extend_from_slice(&(rows.len() as u64).to_le_bytes());
    for row in rows {
        bytes.extend_from_slice(&(row.len() as u64).to_le_bytes());
        for &v in row {
            if !v.is_finite() {
                return Err(model_error("model values must be finite"));
            }
            bytes.extend_from_slice(&v.to_le_bytes());
        }
    }
    Ok(bytes)
}
fn unpack_model(bytes: &[u8]) -> Result<Vec<Vec<f64>>, OpError> {
    if bytes.len() > MODEL_LIMIT * 8 || !bytes.len().is_multiple_of(8) {
        return Err(model_error("invalid portable model size"));
    }
    let mut chunks = bytes.as_chunks::<8>().0.iter();
    let read = |chunk: Option<&[u8; 8]>| -> Result<[u8; 8], OpError> {
        chunk
            .copied()
            .ok_or_else(|| model_error("truncated portable model"))
    };
    let n = usize::try_from(u64::from_le_bytes(read(chunks.next())?)).map_err(model_error)?;
    if n > chunks.len() {
        return Err(model_error("invalid portable row count"));
    }
    let mut rows = Vec::with_capacity(n);
    for _ in 0..n {
        let len = usize::try_from(u64::from_le_bytes(read(chunks.next())?)).map_err(model_error)?;
        if len > chunks.len() {
            return Err(model_error("invalid portable row length"));
        }
        let mut row = Vec::with_capacity(len);
        for _ in 0..len {
            let v = f64::from_le_bytes(read(chunks.next())?);
            if !v.is_finite() {
                return Err(model_error("non-finite portable model"));
            }
            row.push(v);
        }
        rows.push(row);
    }
    if chunks.next().is_some() {
        return Err(model_error("trailing portable model data"));
    }
    Ok(rows)
}
fn model_field<'a>(
    value: &'a Value,
    field: &str,
    fields: &[&str],
) -> Result<&'a Value, OpError> {
    let Value::Record(record) = value else {
        return Err(model_error("expected portable model record"));
    };
    if record.fields.len() != fields.len()
        || fields
            .iter()
            .any(|name| record.fields.iter().filter(|(key, _)| key == name).count() != 1)
    {
        return Err(model_error("invalid portable model fields"));
    }
    Ok(&record
        .fields
        .iter()
        .find(|(key, _)| key == field)
        .ok_or_else(|| model_error("missing model field"))?
        .1)
}
macro_rules! portable_payload {
    ($name:ident,$kind:literal) => {
        struct $name(Vec<Vec<f64>>);
        impl Arg for $name {
            fn ty() -> TypeExpr {
                TypeExpr::Opaque($kind.into())
            }

            fn decode(
                value: &Value,
                _index: usize,
            ) -> Result<Self, OpError> {
                let Value::Opaque(payload) = value else {
                    return Err(model_error("expected opaque model payload"));
                };
                if payload.kind != $kind {
                    return Err(model_error("opaque model kind/version mismatch"));
                }
                Ok(Self(unpack_model(&payload.bytes)?))
            }
        }
        impl Output for $name {
            fn ty() -> TypeExpr {
                <Self as Arg>::ty()
            }

            fn encode(self) -> Result<Value, OpError> {
                Ok(Value::Opaque(crate::kernel::value::OpaqueValue {
                    kind: $kind.into(),
                    bytes: pack_model(&self.0)?,
                }))
            }
        }
    };
}
portable_payload!(KdPayload, "rssn/kdtree-defining-data/v1");
portable_payload!(SplinePayload, "rssn/bspline-defining-data/v1");
struct KdModel {
    payload: KdPayload,
    coordinate_unit: f64,
}
record!(KdModel {
    payload: KdPayload,
    coordinate_unit: f64
});
impl Arg for KdModel {
    fn ty() -> TypeExpr {
        <Self as Output>::ty()
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let fields = ["payload", "coordinate_unit"];
        let unit = numeric(model_field(value, "coordinate_unit", &fields)?, index)?;
        if !crate::kernels::exact_eq(unit, 1.0) {
            return Err(model_error("coordinate_unit must be one canonical unit"));
        }
        Ok(Self {
            payload: KdPayload::decode(model_field(value, "payload", &fields)?, index)?,
            coordinate_unit: unit,
        })
    }
}
struct SplineModel {
    payload: SplinePayload,
    coordinate_unit: f64,
    value_unit: f64,
}
record!(SplineModel {
    payload: SplinePayload,
    coordinate_unit: f64,
    value_unit: f64
});
impl Arg for SplineModel {
    fn ty() -> TypeExpr {
        <Self as Output>::ty()
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let fields = ["payload", "coordinate_unit", "value_unit"];
        let x = numeric(model_field(value, "coordinate_unit", &fields)?, index)?;
        let y = numeric(model_field(value, "value_unit", &fields)?, index)?;
        if !crate::kernels::exact_eq(x, 1.0) || !crate::kernels::exact_eq(y, 1.0) {
            return Err(model_error("model witnesses must be one canonical unit"));
        }
        Ok(Self {
            payload: SplinePayload::decode(model_field(value, "payload", &fields)?, index)?,
            coordinate_unit: x,
            value_unit: y,
        })
    }
}
struct Neighbors {
    indices: Vec<usize>,
    distances: Vec<f64>,
}
record!(Neighbors {indices:Vec<usize>,distances:Vec<f64>});
fn neighbors(items: impl IntoIterator<Item = (usize, f64)>) -> Neighbors {
    let (indices, distances) = items.into_iter().unzip();
    Neighbors { indices, distances }
}
binding!(kd_build,"num_kdtree_build","numerical_vector","Validate a rectangular point tensor and store it as a portable k-d-tree descriptor. Coordinates retain units. Query operators construct the native tree from the defining points.","crate::kernels::geometry::KdTree::build",DimensionPolicy::SpatialBuild,
(points:Vec<Vec<f64>>)->Result<crate::rules::numerical_abi::codecs::KdModel,OpError>,|points:Vec<Vec<f64>>|{super::model_entries(&points)?;crate::kernels::geometry::KdTree::validate_defining_points(&points).map_err(super::model_error)?;Ok(super::KdModel{payload:super::KdPayload(points),coordinate_unit:1.0})});
binding!(kd_nearest,"num_kdtree_nearest","numerical_vector","Exact k-nearest queries on a portable k-d-tree descriptor. Returns index and distance tensors; query coordinates share descriptor units.","crate::kernels::geometry::KdTree::k_nearest",DimensionPolicy::SpatialQuery,
(model:crate::rules::numerical_abi::codecs::KdModel,query:Vec<f64>,k:usize)->Result<crate::rules::numerical_abi::codecs::Neighbors,OpError>,|model:super::KdModel,query:Vec<f64>,k:usize|{if k>super::MODEL_LIMIT{return Err(super::model_error("neighbor count exceeds allocation limit"));} let tree=crate::kernels::geometry::KdTree::build(model.payload.0).map_err(super::model_error)?;Ok(super::neighbors(tree.k_nearest(&query,k).map_err(super::model_error)?))});
binding!(kd_single,"num_kdtree_single_nearest","numerical_vector","Nearest query on the portable native k-d tree, preserving coordinate units. A portable tree requires at least one point.","crate::kernels::geometry::KdTree::nearest",DimensionPolicy::SpatialQuery,
(model:crate::rules::numerical_abi::codecs::KdModel,query:Vec<f64>)->Result<crate::rules::numerical_abi::codecs::Neighbors,OpError>,|model:super::KdModel,query:Vec<f64>|{let tree=crate::kernels::geometry::KdTree::build(model.payload.0).map_err(super::model_error)?;Ok(super::neighbors(tree.nearest(&query).map_err(super::model_error)?))});
binding!(
    kd_nearest_batch,
    "num_kdtree_nearest_batch",
    "numerical_vector",
    "Finds one nearest point per query row in a portable k-d tree, reconstructing the tree once. Indices follow query order and distances retain the descriptor's coordinate dimension. An empty query matrix returns empty index and distance vectors.",
    "crate::kernels::geometry::KdTree::nearest_batch",
    DimensionPolicy::SpatialQuery,
    (model: crate::rules::numerical_abi::codecs::KdModel, queries: Vec<Vec<f64>>) -> Result<crate::rules::numerical_abi::codecs::Neighbors, OpError>,
    |model: super::KdModel, queries: Vec<Vec<f64>>| {
        super::model_entries(&queries)?;
        let tree = crate::kernels::geometry::KdTree::build(model.payload.0)
            .map_err(super::model_error)?;
        let matches = tree.nearest_batch(&queries).map_err(super::model_error)?;
        Ok(super::neighbors(matches.into_iter().flatten()))
    }
);
binding!(kd_radius,"num_kdtree_radius","numerical_vector","Radius queries with coordinate and radius units checked against the descriptor. Returns indices and measured distances.","crate::kernels::geometry::KdTree::radius_search",DimensionPolicy::SpatialRadius,
(model:crate::rules::numerical_abi::codecs::KdModel,query:Vec<f64>,radius:f64)->Result<crate::rules::numerical_abi::codecs::Neighbors,OpError>,|model:super::KdModel,query:Vec<f64>,radius:f64|{if radius<0.0{return Err(super::model_error("radius must be nonnegative"));} let tree=crate::kernels::geometry::KdTree::build(model.payload.0).map_err(super::model_error)?;Ok(super::neighbors(tree.radius_search(&query,radius).map_err(super::model_error)?))});
binding!(bspline_build,"num_bspline_build","numerical_vector","Fit a portable B-spline model. Axis and ordinate units are carried independently by canonical witnesses.","crate::kernels::interp_ext::BSpline::interpolate",DimensionPolicy::SplineBuild,
(xs:Vec<f64>,ys:Vec<f64>,degree:usize)->Result<crate::rules::numerical_abi::codecs::SplineModel,OpError>,|xs:Vec<f64>,ys:Vec<f64>,degree:usize|{if degree>64{return Err(super::model_error("portable spline degree must be at most 64"));}
if xs.len().checked_mul(xs.len()).is_none_or(|n|n>super::MODEL_LIMIT){return Err(super::model_error("portable spline fit exceeds its dense allocation budget"));}let model=crate::kernels::interp_ext::BSpline::interpolate(&xs,&ys,degree).map_err(super::model_error)?;Ok(super::SplineModel{payload:super::SplinePayload(vec![model.knots,model.coeffs,vec![degree as f64]]),coordinate_unit:1.0,value_unit:1.0})});
binding!(bspline_query,"num_bspline_query","numerical_vector","Evaluate a fitted portable B-spline. Validates opaque kind/version, payload bounds and knot structure; query and returned values retain their independent physical units.","crate::kernels::interp_ext::BSpline::eval",DimensionPolicy::SplineQuery,
(model:crate::rules::numerical_abi::codecs::SplineModel,x:f64)->Result<f64,OpError>,|model: super::SplineModel, x: f64| {
    let [knots, coeffs, degree]: [Vec<f64>; 3] = model.payload.0.try_into()
        .map_err(|_| super::model_error("invalid portable spline degree"))?;
    let [degree] = degree.as_slice() else {
        return Err(super::model_error("invalid portable spline degree"));
    };
    if !(0.0..=64.0).contains(degree) || degree.fract() != 0.0 {
        return Err(super::model_error("invalid portable spline degree"));
    }
    crate::kernels::interp_ext::BSpline { knots, coeffs, degree: *degree as usize }
        .eval(x).map_err(super::model_error)
});
binding!(spline_clamped,"num_spline_clamped","numerical_vector","Clamped cubic spline with explicit endpoint derivatives, each in ordinate per abscissa units.","crate::kernels::interp_ext::cubic_spline",DimensionPolicy::ClampedSpline,
(xs:Vec<f64>,ys:Vec<f64>,x:f64,left:f64,right:f64)->Result<f64,OpError>,|xs:Vec<f64>,ys:Vec<f64>,x:f64,left:f64,right:f64|crate::kernels::interp_ext::cubic_spline(&xs,&ys,x,crate::kernels::interp_ext::SplineBc::Clamped(left,right)).map_err(super::model_error));
binding!(spline_notaknot,"num_spline_not_a_knot","numerical_vector","Not-a-knot cubic spline with independent coordinate and ordinate dimensions.","crate::kernels::interp_ext::cubic_spline",DimensionPolicy::Interpolation,
(xs:Vec<f64>,ys:Vec<f64>,x:f64)->Result<f64,OpError>,|xs:Vec<f64>,ys:Vec<f64>,x:f64|crate::kernels::interp_ext::cubic_spline(&xs,&ys,x,crate::kernels::interp_ext::SplineBc::NotAKnot).map_err(super::model_error));
binding!(rbf_scaled,"num_rbf_scaled","numerical_vector","RBF interpolation with independent coordinate and value units. Positive coordinate scale nondimensionalizes distances. Kernel names: gaussian, multiquadric, inverse_multiquadric, thin_plate, cubic, linear.","crate::kernels::interp_ext::rbf_interp",DimensionPolicy::RbfScaled,
(points:Vec<Vec<f64>>,values:Vec<f64>,query:Vec<f64>,kernel:String,scale:f64)->Result<f64,OpError>,|points:Vec<Vec<f64>>,values:Vec<f64>,query:Vec<f64>,kernel:String,scale:f64|{use crate::kernels::interp_ext::RbfKernel as K;if scale<=0.0{return Err(super::model_error("RBF coordinate scale must be positive"));}let k=match kernel.as_str(){"gaussian"=>K::Gaussian(1.0),"multiquadric"=>K::Multiquadric(1.0),"inverse_multiquadric"=>K::InverseMultiquadric(1.0),"thin_plate"=>K::ThinPlateSpline,"cubic"=>K::Cubic,"linear"=>K::Linear,_=>return Err(super::model_error("unknown RBF kernel"))};let p:Vec<Vec<f64>>=points.iter().map(|r|r.iter().map(|x|x/scale).collect()).collect();let refs:Vec<&[f64]>=p.iter().map(Vec::as_slice).collect();let q:Vec<f64>=query.iter().map(|x|x/scale).collect();crate::kernels::interp_ext::rbf_interp(&refs,&values,&q,k).map_err(super::model_error)});
/// Explicit adapters for models whose native Rust APIs are opaque or generic.
pub const OBJECT_EXPORTS: &[&CodecExport] = &[
    &kd_build::EXPORT,
    &kd_nearest::EXPORT,
    &kd_nearest_batch::EXPORT,
    &kd_single::EXPORT,
    &kd_radius::EXPORT,
    &bspline_build::EXPORT,
    &bspline_query::EXPORT,
    &spline_clamped::EXPORT,
    &spline_notaknot::EXPORT,
    &rbf_scaled::EXPORT,
    &ad_gradient_scaled::EXPORT,
    &ad_hessian_scaled::EXPORT,
    &polyfit_normalized::EXPORT,
    &polyval_normalized::EXPORT,
    &roots_normalized::EXPORT,
    &det_units::EXPORT,
    &log_density_scaled::EXPORT,
    &percentile_copy::EXPORT,
];
binding!(ad_gradient_scaled,"num_ad_gradient_scaled","numerical_vector","Exact forward/reverse AD of a dimensionless expression in normalized coordinates x/coordinate_scale. value_scale assigns physical output units. Returns physical gradient in output/coordinate units. Homogeneous physical coordinates; mode is forward or reverse.","crate::backend::autodiff::reverse_gradient",DimensionPolicy::ScaledGradient,
(expression:String,names:Vec<String>,point:Vec<f64>,coordinate_scale:f64,value_scale:f64,mode:String)->Result<Vec<f64>,OpError>,|expression:String,names:Vec<String>,point:Vec<f64>,coordinate_scale:f64,value_scale:f64,mode:String|{if coordinate_scale<=0.0||value_scale<=0.0{return Err(super::model_error("AD reference scales must be positive"));}let point:Vec<_>=point.iter().map(|v|v/coordinate_scale).collect();let g=match mode.as_str(){"forward"=>crate::backend::autodiff::forward_gradient(&expression,&names,&point),"reverse"=>crate::backend::autodiff::reverse_gradient(&expression,&names,&point),_=>return Err(super::model_error("AD mode must be forward or reverse"))}.map_err(super::model_error)?;Ok(g.into_iter().map(|v|v*value_scale/coordinate_scale).collect())});
binding!(ad_hessian_scaled,"num_ad_hessian_scaled","numerical_vector","Exact second-order AD on normalized coordinates, restored to physical output/coordinate-squared units. Expression output is multiplied by value_scale.","crate::backend::autodiff::hessian",DimensionPolicy::ScaledHessian,
(expression:String,names:Vec<String>,point:Vec<f64>,coordinate_scale:f64,value_scale:f64)->Result<crate::kernels::dense::Mat,OpError>,|expression:String,names:Vec<String>,point:Vec<f64>,coordinate_scale:f64,value_scale:f64|{if coordinate_scale<=0.0||value_scale<=0.0{return Err(super::model_error("AD reference scales must be positive"));}let point:Vec<_>=point.iter().map(|v|v/coordinate_scale).collect();let mut h=crate::backend::autodiff::hessian(&expression,&names,&point).map_err(super::model_error)?;for v in &mut h.data {*v *= (value_scale/coordinate_scale)/coordinate_scale;} Ok(h)});
binding!(polyfit_normalized,"num_polyfit_normalized","numerical_vector","Fit a polynomial in x/coordinate_scale. Positive scale carries coordinate units; all returned coefficients retain ordinate units, so powers never require a heterogeneous coefficient tensor.","crate::kernels::interp_ext::polyfit",DimensionPolicy::NormalizedFit,
(xs:Vec<f64>,ys:Vec<f64>,degree:usize,coordinate_scale:f64)->Result<Vec<f64>,OpError>,|xs:Vec<f64>,ys:Vec<f64>,degree:usize,coordinate_scale:f64|{if coordinate_scale<=0.0{return Err(super::model_error("polynomial coordinate scale must be positive"));}let x:Vec<_>=xs.iter().map(|v|v/coordinate_scale).collect();crate::kernels::interp_ext::polyfit(&x,&ys,degree).map_err(super::model_error)});
binding!(polyval_normalized,"num_polyval_normalized","numerical_vector","Evaluate coefficients with ordinate units at normalized x/coordinate_scale; returns ordinate units.","crate::kernels::interp_ext::polyval",DimensionPolicy::NormalizedPolynomial,
(coefficients:Vec<f64>,x:f64,coordinate_scale:f64)->Result<f64,OpError>,|coefficients:Vec<f64>,x:f64,coordinate_scale:f64|{if coordinate_scale<=0.0{return Err(super::model_error("polynomial coordinate scale must be positive"));}Ok(crate::kernels::interp_ext::polyval(&coefficients,x/coordinate_scale))});
binding!(roots_normalized,"num_polynomial_roots_scaled","numerical_vector","Aberth-Ehrlich complex roots of a polynomial expressed in powers of x/coordinate_scale. Coefficients may share arbitrary ordinate units; roots carry coordinate units.","crate::kernels::rootfind::polynomial_roots_aberth",DimensionPolicy::NormalizedRoots,
(coefficients:Vec<f64>,coordinate_scale:f64)->Result<Vec<num_complex::Complex64>,OpError>,|coefficients:Vec<f64>,coordinate_scale:f64|{if coordinate_scale<=0.0{return Err(super::model_error("root coordinate scale must be positive"));}Ok(crate::kernels::rootfind::polynomial_roots_aberth(&coefficients).map_err(|e|super::model_error(format!("polynomial root solve failed: {e:?}")))?.into_iter().map(|(re,im)|num_complex::Complex64::new(re*coordinate_scale,im*coordinate_scale)).collect())});

binding!(det_units,"num_dense_det_units","numerical_matrix","Determinant with runtime matrix-order dimension exponent. The Value result retains entry-units raised to matrix order; only square matrices are accepted.","crate::kernels::dense::Lu::det",DimensionPolicy::DeterminantUnits,
(a:dense::Mat)->Result<f64,OpError>,|a:dense::Mat|{if a.rows!=a.cols{return Err(super::model_error("determinant requires a square matrix"));}match dense::lu_factor(&a) {Ok(lu)=>Ok(lu.det()),Err(dense::DenseError::Singular)=>Ok(0.0),Err(error)=>Err(dense_error(error))}});
binding!(log_density_scaled,"num_multivariate_normal_log_pdf_scaled","numerical_stats","Multivariate log density relative to coordinate_scale^n, where n is the vector length. An explicit positive scale makes the logarithm dimensionless and defines the reference measure.","crate::kernels::distributions::multivariate_normal_log_pdf",DimensionPolicy::NormalizedLogDensity,
(mean:Vec<f64>,covariance:dense::Mat,x:Vec<f64>,coordinate_scale:f64)->Result<f64,OpError>,|mean:Vec<f64>,covariance:dense::Mat,x:Vec<f64>,coordinate_scale:f64|{if coordinate_scale<=0.0{return Err(super::model_error("density reference scale must be positive"));}let m:Vec<_>=mean.iter().map(|v|v/coordinate_scale).collect();let mut c=covariance;for v in &mut c.data{*v=(*v/coordinate_scale)/coordinate_scale;}let x:Vec<_>=x.iter().map(|v|v/coordinate_scale).collect();crate::kernels::distributions::multivariate_normal_log_pdf(&m,&c,&x).map_err(super::model_error)});

binding!(percentile_copy,"num_stats_percentile","numerical_stats","Measured percentile with a private mutable copy of the input. The native sorting/selection routine cannot mutate the caller's tensor; percentile p is dimensionless in 0 through 100.","crate::kernels::stats::percentile",DimensionPolicy::Preserve,
(data:Vec<f64>,p:f64)->f64,|mut data:Vec<f64>,p:f64| crate::kernels::stats::percentile(&mut data,p));
