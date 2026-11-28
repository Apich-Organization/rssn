//! True parallel native quadrature using owned immutable compiled programs.
use crate::graph::RuleError;
use crate::graph::dim::DimRule;
use crate::graph::rule::Installer;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::Exp;
use crate::kernel::signature::Rational;
use crate::kernel::signature::Signature;
use crate::kernel::signature::TypeExpr;
use crate::kernels::quadrature as q;
use crate::rules::numerical_abi::codecs::Arg;
use crate::rules::numerical_abi::codecs::Output;
use crate::rules::numerical_abi::native_backend::CompiledSnapshot;
use crate::rules::numerical_abi::native_backend::PROGRAM;
use crate::rules::numerical_abi::native_backend::compiled_snapshot;
use crate::rules::numerical_abi::native_object_api::Binding;
use crate::rules::numerical_abi::native_object_api::Dynamic;
use crate::rules::numerical_abi::native_object_api::Input;
use crate::rules::numerical_abi::native_object_api::object_binding;
use crate::rules::numerical_abi::native_object_api::{
    self as object,
};
pub(super) fn limit(n: usize) -> Result<(), OpError> {
    if n == 0 || n > 1_000_000 {
        Err(object::error(
            "parallel quadrature work limit must be1..1000000",
        ))
    } else {
        Ok(())
    }
}
pub(super) fn bounds(
    ctx: &CallCtx<'_>,
    lo: &Input,
    hi: &Input,
    function: &Input,
    budget: usize,
) -> Result<(Vec<f64>, Vec<f64>, CompiledSnapshot, DimExpr), OpError> {
    limit(budget)?;
    let d = object::dimension(&lo.0)?;
    object::same_dimension(&d, &object::dimension(&hi.0)?)?;
    let lo = <Vec<f64> as Arg>::decode(&lo.0, 0)?;
    let hi = <Vec<f64> as Arg>::decode(&hi.0, 1)?;
    if lo.is_empty()
        || lo.len() > 16
        || lo.len() != hi.len()
        || lo.iter().zip(&hi).any(|(a, b)| a >= b)
    {
        return Err(object::error(
            "quadrature bounds need1..16 finite ordered same-size coordinates",
        ));
    }
    let snapshot = compiled_snapshot(ctx, &function.0, budget)?;
    if snapshot.input_dimensions.len() != lo.len() {
        return Err(object::error(
            "compiled quadrature function arity differs from bounds",
        ));
    }
    for input in &snapshot.input_dimensions {
        object::same_dimension(&d, input)?;
    }
    let dimension = DimExpr::Mul(vec![
        snapshot.output_dimension.clone(),
        DimExpr::Pow(Box::new(d), Exp::Rat(Rational::new(lo.len() as i64, 1))),
    ]);
    Ok((lo, hi, snapshot, dimension))
}
pub(super) fn output(
    v: q::Integral,
    d: DimExpr,
) -> Result<Dynamic, OpError> {
    let error = if v.error.is_finite() {
        crate::kernel::Value::Variant(crate::kernel::value::VariantValue {
            tag: "Some".into(),
            payload: vec![object::quantity(v.error, d.clone())?],
        })
    } else {
        crate::kernel::Value::Variant(crate::kernel::value::VariantValue {
            tag: "None".into(),
            payload: vec![],
        })
    };
    Ok(Dynamic(object::record(vec![
        ("value", object::quantity(v.value, d)?),
        ("error", error),
        ("evaluations", v.evaluations.encode()?),
        ("converged", v.converged.encode()?),
    ])))
}
object_binding!(mc,"num_analysis_monte_carlo_parallel","crate::kernels::quadrature::monte_carlo_par","True native parallel Monte Carlo over physical bounds using an owned immutable compiled program; same-unit coordinates, arbitrary output units, actual native error/status and explicit sample/evaluation limits.",false,(lo:Input,hi:Input,samples:usize,seed:u64,budget:usize,function:Input)->Dynamic,|ctx:&CallCtx<'_>,lo:Input,hi:Input,samples:usize,seed:u64,budget:usize,function:Input|{
 limit(samples)?;let(lo,hi,f,d)=bounds(ctx,&lo,&hi,&function,budget)?;let out=q::monte_carlo_par(|x|f.sample(x),&lo,&hi,samples,seed);f.finish()?;output(out,d)
});
object_binding!(qmc,"num_analysis_quasi_monte_carlo_parallel","crate::kernels::quadrature::quasi_monte_carlo_par","True native parallel Halton/Sobol cubature with checked physical compiled-function references; sequence0=Halton,1=Sobol.",false,(lo:Input,hi:Input,samples:usize,sequence:u64,budget:usize,function:Input)->Dynamic,|ctx:&CallCtx<'_>,lo:Input,hi:Input,samples:usize,sequence:u64,budget:usize,function:Input|{
 limit(samples)?;let(lo,hi,f,d)=bounds(ctx,&lo,&hi,&function,budget)?;let sequence=match sequence{0=>q::LowDiscrepancy::Halton,1 if lo.len()<=8=>q::LowDiscrepancy::Sobol,_=>return Err(object::error("sequence must be0(Halton)or1(Sobol,max8dimensions)"))};let out=q::quasi_monte_carlo_par(|x|f.sample(x),&lo,&hi,samples,sequence);f.finish()?;output(out,d)
});
object_binding!(nested,"num_analysis_cubature_nested_parallel","crate::kernels::quadrature::cubature_nested_par","True native parallel nested cubature; absolute tolerance has output times coordinate-volume dimension.",false,(lo:Input,hi:Input,tolerance:Input,budget:usize,function:Input)->Dynamic,|ctx:&CallCtx<'_>,lo:Input,hi:Input,tolerance:Input,budget:usize,function:Input|{
 let(lo,hi,f,d)=bounds(ctx,&lo,&hi,&function,budget)?;if lo.len()>4{return Err(object::error("nested cubature is bounded to4coordinates"));}object::same_dimension(&d,&object::dimension(&tolerance.0)?)?;let tolerance=super::super::numeric(&tolerance.0,2)?;if tolerance<=0.0{return Err(object::error("absolute cubature tolerance must be positive"));}let out=q::cubature_nested_par(|x|f.sample(x),&lo,&hi,tolerance);f.finish()?;output(out,d)
});
object_binding!(gk,"num_analysis_integrate_gk_parallel","crate::kernels::quadrature::integrate_gk_par","True native parallel Gauss-Kronrod quadrature using physical coherent bounds and compiled reference scales; rule0=G7K15,1=G10K21.",false,(lo:Input,hi:Input,rule:u64,absolute:Input,relative:f64,panels:usize,budget:usize,function:Input)->Dynamic,|ctx:&CallCtx<'_>,lo:Input,hi:Input,rule:u64,absolute:Input,relative:f64,panels:usize,budget:usize,function:Input|{
 limit(panels)?;let coordinate=object::dimension(&lo.0)?;object::same_dimension(&coordinate,&object::dimension(&hi.0)?)?;
 let a=super::super::numeric(&lo.0,0)?;let b=super::super::numeric(&hi.0,1)?;if a>=b{return Err(object::error("parallel quadrature requires ordered finite bounds"));}
 let f=compiled_snapshot(ctx,&function.0,budget)?;if f.input_dimensions.len()!=1{return Err(object::error("compiled quadrature integrand needs one argument"));}object::same_dimension(&coordinate,&f.input_dimensions[0])?;
 let d=DimExpr::Mul(vec![coordinate,f.output_dimension.clone()]);object::same_dimension(&d,&object::dimension(&absolute.0)?)?;
 let absolute=super::super::numeric(&absolute.0,3)?;if absolute<0.0||relative<0.0||absolute==0.0&&relative==0.0{return Err(object::error("quadrature tolerances must be nonnegative and one positive"));}
 let rule=match rule{0=>q::GkRule::G7K15,1=>q::GkRule::G10K21,_=>return Err(object::error("rule must be0or1"))};let out=q::integrate_gk_par(|x|f.sample(&[x]),a,b,rule,absolute,relative,panels);f.finish()?;output(out,d)
});
object_binding!(mc_serial,"num_analysis_monte_carlo_compiled","crate::kernels::quadrature::monte_carlo","Native serial Monte Carlo over physical bounds using an owned immutable compiled program; same-unit coordinates, arbitrary output units, actual native error/status and explicit sample/evaluation limits.",false,(lo:Input,hi:Input,samples:usize,seed:u64,budget:usize,function:Input)->Dynamic,|ctx:&CallCtx<'_>,lo:Input,hi:Input,samples:usize,seed:u64,budget:usize,function:Input|{
 limit(samples)?;let(lo,hi,f,d)=bounds(ctx,&lo,&hi,&function,budget)?;let out=q::monte_carlo(|x|f.sample(x),&lo,&hi,samples,seed);f.finish()?;output(out,d)
});

object_binding!(qmc_serial,"num_analysis_quasi_monte_carlo_compiled","crate::kernels::quadrature::quasi_monte_carlo","Native serial Halton/Sobol cubature with checked physical compiled-function references; sequence0=Halton,1=Sobol.",false,(lo:Input,hi:Input,samples:usize,sequence:u64,budget:usize,function:Input)->Dynamic,|ctx:&CallCtx<'_>,lo:Input,hi:Input,samples:usize,sequence:u64,budget:usize,function:Input|{
 limit(samples)?;let(lo,hi,f,d)=bounds(ctx,&lo,&hi,&function,budget)?;let sequence=match sequence{0=>q::LowDiscrepancy::Halton,1 if lo.len()<=8=>q::LowDiscrepancy::Sobol,_=>return Err(object::error("sequence must be0(Halton)or1(Sobol,max8dimensions)"))};let out=q::quasi_monte_carlo(|x|f.sample(x),&lo,&hi,samples,sequence);f.finish()?;output(out,d)
});

object_binding!(nested_serial,"num_analysis_cubature_nested_compiled","crate::kernels::quadrature::cubature_nested","Native serial nested cubature; absolute tolerance has output times coordinate-volume dimension.",false,(lo:Input,hi:Input,tolerance:Input,budget:usize,function:Input)->Dynamic,|ctx:&CallCtx<'_>,lo:Input,hi:Input,tolerance:Input,budget:usize,function:Input|{
 let(lo,hi,f,d)=bounds(ctx,&lo,&hi,&function,budget)?;if lo.len()>4{return Err(object::error("nested cubature is bounded to4coordinates"));}object::same_dimension(&d,&object::dimension(&tolerance.0)?)?;let tolerance=super::super::numeric(&tolerance.0,2)?;if tolerance<=0.0{return Err(object::error("absolute cubature tolerance must be positive"));}let out=q::cubature_nested(|x|f.sample(x),&lo,&hi,tolerance);f.finish()?;output(out,d)
});

object_binding!(antithetic,"num_analysis_monte_carlo_antithetic_compiled","crate::kernels::quadrature::monte_carlo_antithetic","Native antithetic Monte Carlo with physical compiled-function references and bounded pair count.",false,(lo:Input,hi:Input,pairs:usize,seed:u64,budget:usize,function:Input)->Dynamic,|ctx:&CallCtx<'_>,lo:Input,hi:Input,pairs:usize,seed:u64,budget:usize,function:Input|{
 limit(pairs)?;if pairs>500_000{return Err(object::error("antithetic pair count exceeds one million samples"));}let(lo,hi,f,d)=bounds(ctx,&lo,&hi,&function,budget)?;let out=q::monte_carlo_antithetic(|x|f.sample(x),&lo,&hi,pairs,seed);f.finish()?;output(out,d)
});
object_binding!(stratified,"num_analysis_monte_carlo_stratified_compiled","crate::kernels::quadrature::monte_carlo_stratified","Native stratified Monte Carlo with checked exponential grid/sample cost and independent physical output units.",false,(lo:Input,hi:Input,strata:usize,per_stratum:usize,seed:u64,budget:usize,function:Input)->Dynamic,|ctx:&CallCtx<'_>,lo:Input,hi:Input,strata:usize,per_stratum:usize,seed:u64,budget:usize,function:Input|{
 limit(strata)?;limit(per_stratum)?;let(lo,hi,f,d)=bounds(ctx,&lo,&hi,&function,budget)?;
 if strata.checked_pow(lo.len()as u32).and_then(|n|n.checked_mul(per_stratum)).is_none_or(|n|n>1_000_000){return Err(object::error("stratification grid exceeds sample budget"));}
 let out=q::monte_carlo_stratified(|x|f.sample(x),&lo,&hi,strata,per_stratum,seed);f.finish()?;output(out,d)
});
object_binding!(genz,"num_analysis_cubature_genz_malik_compiled","crate::kernels::quadrature::cubature_genz_malik","Native Genz-Malik cubature: physical same-unit coordinates, independently measured output, explicit evaluation limit and volume-dimensioned absolute tolerance.",false,(lo:Input,hi:Input,absolute:Input,relative:f64,budget:usize,function:Input)->Dynamic,|ctx:&CallCtx<'_>,lo:Input,hi:Input,absolute:Input,relative:f64,budget:usize,function:Input|{
 let(lo,hi,f,d)=bounds(ctx,&lo,&hi,&function,budget)?;object::same_dimension(&d,&object::dimension(&absolute.0)?)?;let absolute=super::super::numeric(&absolute.0,2)?;
 if absolute<0.0||relative<0.0||absolute==0.0&&relative==0.0{return Err(object::error("cubature tolerances must be nonnegative and one positive"));}
 let out=q::cubature_genz_malik(|x|f.sample(x),&lo,&hi,absolute,relative,budget);f.finish()?;output(out,d)
});
object_binding!(gk_panel,"num_analysis_gk_panel_parallel","crate::kernels::quadrature::gk_panel_par","Actual parallel native Gauss-Kronrod panel; returns physical integral/error and evaluation count without inventing a convergence verdict.",false,(lo:Input,hi:Input,rule:u64,budget:usize,function:Input)->Dynamic,|ctx:&CallCtx<'_>,lo:Input,hi:Input,rule:u64,budget:usize,function:Input|{
 let d=object::dimension(&lo.0)?;object::same_dimension(&d,&object::dimension(&hi.0)?)?;let a=super::super::numeric(&lo.0,0)?;let b=super::super::numeric(&hi.0,1)?;
 let f=compiled_snapshot(ctx,&function.0,budget)?;if f.input_dimensions.len()!=1{return Err(object::error("Gauss-Kronrod panel requires scalar compiled input"));}object::same_dimension(&d,&f.input_dimensions[0])?;
 let rule=match rule{0=>q::GkRule::G7K15,1=>q::GkRule::G10K21,_=>return Err(object::error("rule must be0or1"))};let d=DimExpr::Mul(vec![d,f.output_dimension.clone()]);
 let(v,e,n)=q::gk_panel_par(&|x|f.sample(&[x]),a,b,rule);f.finish()?;
 Ok(Dynamic(object::record(vec![("value",object::quantity(v,d.clone())?),("error",object::quantity(e,d)?),("evaluations",Output::encode(n)?)])))
});

#[derive(Clone, Copy)]
enum FixedRule {
    GaussLegendre,
    ClenshawCurtis,
}

impl FixedRule {
    fn integrate(
        self,
        ctx: &CallCtx<'_>,
        lo: &Input,
        hi: &Input,
        order: usize,
        output_reference: &Input,
        function: &Input,
    ) -> Result<Dynamic, OpError> {
        if !(1..=4096).contains(&order) {
            return Err(object::error("fixed quadrature order must be in 1..=4096"));
        }
        let coordinate = object::dimension(&lo.0)?;
        object::same_dimension(&coordinate, &object::dimension(&hi.0)?)?;
        let a = crate::rules::numerical_abi::numeric(&lo.0, 0)?;
        let b = crate::rules::numerical_abi::numeric(&hi.0, 1)?;
        let evaluations = match self {
            | Self::GaussLegendre => order,
            | Self::ClenshawCurtis => order.max(2).next_multiple_of(2) + 1,
        };
        let program = compiled_snapshot(ctx, &function.0, evaluations)?;
        let [input] = program.input_dimensions.as_slice() else {
            return Err(object::error(
                "fixed quadrature requires a scalar compiled integrand",
            ));
        };
        object::same_dimension(&coordinate, input)?;
        crate::rules::numerical_abi::numeric(&output_reference.0, 3)?;
        let integrand = object::dimension(&output_reference.0)?;
        object::same_dimension(&integrand, &program.output_dimension)?;
        let units = DimExpr::Mul(vec![coordinate, integrand]);
        if crate::kernels::exact_eq(a, b) {
            return object::quantity(0.0, units).map(Dynamic);
        }
        let f = |x| program.sample(&[x]);
        let value = match self {
            | Self::GaussLegendre => q::gauss_legendre_par(f, a, b, order),
            | Self::ClenshawCurtis => q::clenshaw_curtis_par(f, a, b, order),
        };
        program.finish()?;
        object::quantity(value, units).map(Dynamic)
    }
}

object_binding!(
    gauss_fixed,
    "num_analysis_gauss_legendre_parallel",
    "crate::kernels::quadrature::gauss_legendre_par",
    "Fixed Gauss–Legendre quadrature of a scalar compiled integrand. Bounds share the input units; the supplied output reference witnesses the integrand units. Order is 1..=4096.",
    false,
    (lo: Input, hi: Input, order: usize, output_reference: Input, function: Input) -> Dynamic,
    |ctx: &CallCtx<'_>, lo, hi, order, output_reference, function| {
        super::FixedRule::GaussLegendre.integrate(ctx, &lo, &hi, order, &output_reference, &function)
    }
);

object_binding!(
    clenshaw_fixed,
    "num_analysis_clenshaw_curtis_parallel",
    "crate::kernels::quadrature::clenshaw_curtis_par",
    "Fixed Clenshaw–Curtis quadrature of a scalar compiled integrand. Order is rounded up to an even number of at least two; the return unit is coordinate times integrand.",
    false,
    (lo: Input, hi: Input, order: usize, output_reference: Input, function: Input) -> Dynamic,
    |ctx: &CallCtx<'_>, lo, hi, order, output_reference, function| {
        super::FixedRule::ClenshawCurtis.integrate(ctx, &lo, &hi, order, &output_reference, &function)
    }
);

object_binding!(
    sum_compensated,
    "num_analysis_sum_compensated",
    "crate::kernels::par::sum_kahan_fixed",
    "Deterministic Neumaier-compensated sum of a finite measured vector. Positive block and grain parameters control accumulation and parallel scheduling. The return unit matches the input.",
    false,
    (values: Input, block: usize, grain: usize) -> Dynamic,
    |_ctx: &CallCtx<'_>, values: Input, block: usize, grain: usize| {
        if block == 0 || grain == 0 {
            return Err(object::error("compensated sum block and grain must be positive"));
        }
        let units = object::dimension(&values.0)?;
        let values = <Vec<f64> as Arg>::decode(&values.0, 0)?;
        let sum =
            crate::kernels::par::sum_kahan_fixed(values.len(), block, grain, |index| values[index]);
        object::quantity(sum, units).map(Dynamic)
    }
);

fn fixed_signature() -> Signature {
    let coordinate = DimExpr::Var("Coordinate".into());
    let integrand = DimExpr::Var("Integrand".into());
    let bounds = TypeExpr::Float.with_dim(coordinate.clone());
    Signature::new(
        vec![
            bounds.clone(),
            bounds,
            TypeExpr::Int,
            TypeExpr::Float.with_dim(integrand.clone()),
            TypeExpr::Opaque(PROGRAM.into()),
        ],
        TypeExpr::Float.with_dim(DimExpr::Mul(vec![coordinate, integrand])),
    )
}

fn sum_signature() -> Signature {
    let sample = TypeExpr::Float.with_dim(DimExpr::Var("Sample".into()));
    Signature::new(
        vec![
            object::tensor_type(1, sample.clone()),
            TypeExpr::Int,
            TypeExpr::Int,
        ],
        sample,
    )
}

static GAUSS_FIXED: Binding = Binding {
    signature: fixed_signature,
    ..gauss_fixed::BINDING
};
static CLENSHAW_FIXED: Binding = Binding {
    signature: fixed_signature,
    ..clenshaw_fixed::BINDING
};
static SUM_COMPENSATED: Binding = Binding {
    signature: sum_signature,
    ..sum_compensated::BINDING
};

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    object::install(
        i,
        &[
            &mc::BINDING,
            &qmc::BINDING,
            &nested::BINDING,
            &gk::BINDING,
            &mc_serial::BINDING,
            &qmc_serial::BINDING,
            &nested_serial::BINDING,
            &antithetic::BINDING,
            &stratified::BINDING,
            &genz::BINDING,
            &gk_panel::BINDING,
        ],
        "numerical_vector",
    )?;
    object::install_with_dimension(
        i,
        &[&GAUSS_FIXED, &CLENSHAW_FIXED],
        "numerical_vector",
        DimRule::Integral { f: 3, x: 0 },
    )?;
    object::install_with_dimension(i, &[&SUM_COMPENSATED], "numerical_vector", DimRule::Same(0))
}
