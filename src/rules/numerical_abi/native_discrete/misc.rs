#![allow(clippy::wildcard_imports)] // Checked exact/discrete adapters.
use super::*;
use crate::kernels::number_exact;
use crate::kernels::ordering;
use crate::kernels::ordering_ml;
use crate::kernels::polynomial::Polynomial;
use crate::kernels::quadratic_sieve;
use crate::kernels::quadratic_sieve::QsConfig;
use crate::kernels::quadratic_sieve::QsLimits;
use crate::kernels::quadratic_sieve::QsStats;
use crate::kernels::real_roots;
use num_bigint::BigInt;
fn adjacency(data: Vec<Vec<usize>>) -> Result<Vec<Vec<usize>>, OpError> {
    limit(data.len(), 512)?;
    let n = data.len();
    let mut edges = std::collections::BTreeSet::new();
    for (u, row) in data.iter().enumerate() {
        for &v in row {
            if v >= n {
                return Err(error("adjacency neighbour outside graph"));
            }
            edges.insert((u, v));
        }
    }
    limit(edges.len(), 65536)?;
    if edges.iter().any(|(u, v)| !edges.contains(&(*v, *u))) {
        return Err(error("ordering requires an undirected adjacency relation"));
    }
    Ok(data)
}
fn nd_options(value: &Input) -> Result<ordering_ml::NdOptions, OpError> {
    let leaf = <usize as Arg>::decode(api::field(&value.0, "leaf")?, 0)?;
    let imbalance = numeric(api::field(&value.0, "imbalance")?, 0)?;
    let passes = <usize as Arg>::decode(api::field(&value.0, "fm_passes")?, 0)?;
    if leaf == 0 || leaf > 4096 || !(0.0..1.0).contains(&imbalance) || passes > 32 {
        return Err(error("nested dissection options exceed domain/work bounds"));
    }
    Ok(ordering_ml::NdOptions {
        leaf,
        imbalance,
        fm_passes: passes,
    })
}
fn positive_tolerance(
    value: &Input,
    scale: &Input,
) -> Result<(f64, f64, DimExpr), OpError> {
    let (r, d) = api::reference(&scale.0)?;
    let (t, e) = scalar(value)?;
    api::same_dimension(&d, &e)?;
    let normalized = t / r;
    if normalized <= 0.0 || !normalized.is_finite() {
        return Err(error(
            "root tolerance is not positive/representable at its coordinate reference",
        ));
    }
    Ok((normalized, r, d))
}
fn root_polynomial(p: &Polynomial) -> Result<(), OpError> {
    if p.coeffs.is_empty() || p.coeffs.iter().all(|v| crate::kernels::exact_eq(*v, 0.0)) {
        return Err(error("zero polynomial has no finite isolated root set"));
    }
    if crate::kernels::exact_eq(p.coeffs[0], 0.0) {
        return Err(error(
            "root polynomial must have a nonzero leading coefficient",
        ));
    }
    Ok(())
}
fn sieve_limits(limits: &QsLimits) -> Result<(), OpError> {
    if limits.bound > 20000
        || limits.interval > 200000
        || limits.polynomials > 4096
        || limits.max_sieve_positions > 10_000_000
    {
        return Err(error(
            "quadratic sieve request exceeds factor-base/allocation/work budget",
        ));
    }
    Ok(())
}
impl NativeType for QsLimits {
    const KIND: &'static str = "rssn/native/discrete/qs-limits/v1";
    const MAX_BYTES: usize = 128;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }
}
object_binding!(number_exact_factor,"num_native_number_exact_factor","crate::kernels::number_exact::factor","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(n:i64)->Vec<[i64;2]>,|_ctx:&CallCtx<'_>,n:i64|{Ok(number_exact::factor(n).into_iter().map(|(p,e)|[p,i64::from(e)]).collect())});
object_binding!(number_exact_is_prime,"num_native_number_exact_is_prime","crate::kernels::number_exact::is_prime","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(n:u64)->bool,|_ctx:&CallCtx<'_>,n:u64|{Ok(number_exact::is_prime(n))});
object_binding!(number_exact_gcd,"num_native_number_exact_gcd","crate::kernels::number_exact::gcd","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:i64,b:i64)->i64,|_ctx:&CallCtx<'_>,a:i64,b:i64|{Ok(number_exact::gcd(a,b))});
object_binding!(number_exact_lcm,"num_native_number_exact_lcm","crate::kernels::number_exact::lcm","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:i64,b:i64)->i64,|_ctx:&CallCtx<'_>,a:i64,b:i64|{Ok(number_exact::lcm(a,b))});
object_binding!(number_exact_is_negative,"num_native_number_exact_is_negative","crate::kernels::number_exact::is_negative","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(n:BigInt)->bool,|_ctx:&CallCtx<'_>,n:BigInt|{Ok(number_exact::is_negative(&n))});
object_binding!(number_exact_primes,"num_native_number_exact_primes","crate::kernels::number_exact::primes","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(limit_n:i64)->Vec<i64>,|_ctx:&CallCtx<'_>,limit_n:i64|{if limit_n>1_000_000{return Err(error("prime sieve limit exceeds allocation budget"));}Ok(number_exact::primes(limit_n))});
object_binding!(number_exact_binomial,"num_native_number_exact_binomial","crate::kernels::number_exact::binomial","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(n:i64,k:i64)->BigInt,|_ctx:&CallCtx<'_>,n:i64,k:i64|{if n>=0&&k>=0&&k<=n&&k.min(n-k)>4096{return Err(error("exact binomial multiplication budget exceeded"));}Ok(number_exact::binomial(n,k))});
object_binding!(number_exact_permutations,"num_native_number_exact_permutations","crate::kernels::number_exact::permutations","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(n:i64)->BigInt,|_ctx:&CallCtx<'_>,n:i64|{if n>10000{return Err(error("factorial multiplication budget exceeded"));}Ok(number_exact::permutations(n))});
object_binding!(ordering_amd_order,"num_native_ordering_amd_order","crate::kernels::ordering::amd_order","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(adj:Vec<Vec<usize>>)->Vec<usize>,|_ctx:&CallCtx<'_>,adj:Vec<Vec<usize>>|{let adj=adjacency(adj)?;Ok(ordering::amd_order(&adj))});
object_binding!(ordering_nested_dissection_order,"num_native_ordering_nested_dissection_order","crate::kernels::ordering::nested_dissection_order","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(adj:Vec<Vec<usize>>)->Vec<usize>,|_ctx:&CallCtx<'_>,adj:Vec<Vec<usize>>|{let adj=adjacency(adj)?;Ok(ordering::nested_dissection_order(&adj))});
object_binding!(ordering_ml_hybrid_nested_dissection_order,"num_native_ordering_ml_hybrid_nested_dissection_order","crate::kernels::ordering_ml::hybrid_nested_dissection_order","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(adj:Vec<Vec<usize>>)->Vec<usize>,|_ctx:&CallCtx<'_>,adj:Vec<Vec<usize>>|{let adj=adjacency(adj)?;Ok(ordering_ml::hybrid_nested_dissection_order(&adj))});
object_binding!(ordering_ml_multilevel_nested_dissection_order,"num_native_ordering_ml_multilevel_nested_dissection_order","crate::kernels::ordering_ml::multilevel_nested_dissection_order","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(adj:Vec<Vec<usize>>)->Vec<usize>,|_ctx:&CallCtx<'_>,adj:Vec<Vec<usize>>|{let adj=adjacency(adj)?;Ok(ordering_ml::multilevel_nested_dissection_order(&adj))});
object_binding!(ordering_cholesky_fill,"num_native_ordering_cholesky_fill","crate::kernels::ordering::cholesky_fill","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(adj:Vec<Vec<usize>>,permutation:Vec<usize>)->usize,|_ctx:&CallCtx<'_>,adj:Vec<Vec<usize>>,permutation:Vec<usize>|{let adj=adjacency(adj)?;if permutation.len()!=adj.len(){return Err(error("ordering permutation dimension mismatch"));}let mut check=permutation.clone();check.sort_unstable();if check!=(0..adj.len()).collect::<Vec<_>>(){return Err(error("ordering requires a bijective permutation"));}Ok(ordering::cholesky_fill(&adj,&permutation))});
object_binding!(ordering_ml_multilevel_nested_dissection_with,"num_native_ordering_ml_multilevel_nested_dissection_with","crate::kernels::ordering_ml::multilevel_nested_dissection_with","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(adj:Vec<Vec<usize>>,options:Input)->Vec<usize>,|_ctx:&CallCtx<'_>,adj:Vec<Vec<usize>>,options:Input|{let adj=adjacency(adj)?;let options=nd_options(&options)?;Ok(ordering_ml::multilevel_nested_dissection_with(&adj,&options))});
object_binding!(ordering_ml_vertex_separator_bisection,"num_native_ordering_ml_vertex_separator_bisection","crate::kernels::ordering_ml::vertex_separator_bisection","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(adj:Vec<Vec<usize>>,options:Input)->(Vec<usize>,Vec<usize>,Vec<usize>),|_ctx:&CallCtx<'_>,adj:Vec<Vec<usize>>,options:Input|{let adj=adjacency(adj)?;let options=nd_options(&options)?;Ok(ordering_ml::vertex_separator_bisection(&adj,&options))});
object_binding!(polynomial_normalized_new,"num_native_polynomial_normalized_new","crate::kernels::polynomial::Polynomial::new","Native polynomial in explicitly normalized dimensionless coordinate u, highest degree first. Real-root methods require an explicit positive coordinate reference and restore physical roots/tolerances.",true,(coefficients:Vec<f64>)->NewObject<Polynomial>,|_ctx:&CallCtx<'_>,coefficients:Vec<f64>|{limit(coefficients.len(),129)?;Ok(NewObject::normalized(Polynomial::new(coefficients)))});
object_binding!(real_roots_sturm_sequence,"num_native_real_roots_sturm_sequence","crate::kernels::real_roots::sturm_sequence","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(poly:Handle<Polynomial>)->Many<NewObject<Polynomial>>,|ctx:&CallCtx<'_>,poly:Handle<Polynomial>|{poly.with(ctx,|p|{if p.coeffs.first().is_some_and(|v|crate::kernels::exact_eq(*v,0.0))&&p.coeffs.len()>1{return Err(error("Sturm polynomial leading coefficient must be nonzero"));}Ok(Many(real_roots::sturm_sequence(p).into_iter().map(NewObject::normalized).collect()))})});
object_binding!(real_roots_find_roots,"num_native_real_roots_find_roots","crate::kernels::real_roots::find_roots","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(poly:Handle<Polynomial>,tolerance:Input,coordinate_scale:Input)->Dynamic,|ctx:&CallCtx<'_>,poly:Handle<Polynomial>,tolerance:Input,coordinate_scale:Input|{let(t,r,d)=positive_tolerance(&tolerance,&coordinate_scale)?;poly.with(ctx,|p|{root_polynomial(p)?;quantities(real_roots::find_roots(p,t).map_err(error)?.into_iter().map(|x|x*r).collect(),&(d))})});
object_binding!(real_roots_isolate_real_roots,"num_native_real_roots_isolate_real_roots","crate::kernels::real_roots::isolate_real_roots","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(poly:Handle<Polynomial>,precision:Input,coordinate_scale:Input)->Dynamic,|ctx:&CallCtx<'_>,poly:Handle<Polynomial>,precision:Input,coordinate_scale:Input|{let(t,r,d)=positive_tolerance(&precision,&coordinate_scale)?;poly.with(ctx,|p|{root_polynomial(p)?;let values=real_roots::isolate_real_roots(p,t).map_err(error)?.into_iter().map(|(a,b)|vec![a*r,b*r]).collect();quantities_matrix(values,&(d))})});
object_binding!(real_roots_refine_root_bisection,"num_native_real_roots_refine_root_bisection","crate::kernels::real_roots::refine_root_bisection","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(poly:Handle<Polynomial>,interval:Input,tolerance:Input,coordinate_scale:Input)->Dynamic,|ctx:&CallCtx<'_>,poly:Handle<Polynomial>,interval:Input,tolerance:Input,coordinate_scale:Input|{let(t,r,d)=positive_tolerance(&tolerance,&coordinate_scale)?;api::same_dimension(&d,&api::dimension(&interval.0)?)?;let interval=<[f64;2] as Arg>::decode(&interval.0,1)?;if interval[0]>interval[1]{return Err(error("root interval endpoints must increase"));}poly.with(ctx,|p|{root_polynomial(p)?;extended_quantity(real_roots::refine_root_bisection(p,(interval[0]/r,interval[1]/r),t)*r,d)})});
object_binding!(qs_limits_for_digits,"num_native_qs_limits_for_digits","crate::kernels::quadratic_sieve::QsLimits::for_digits","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(digits:u32)->NewObject<QsLimits>,|_ctx:&CallCtx<'_>,digits:u32|{Ok(NewObject::normalized(QsLimits::for_digits(digits)))});
object_binding!(qs_is_probable_prime,"num_native_qs_is_probable_prime","crate::kernels::quadratic_sieve::is_probable_prime","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(n:BigInt)->bool,|_ctx:&CallCtx<'_>,n:BigInt|{if n.bits()>4096{return Err(error("primality input exceeds bit-work budget"));}Ok(quadratic_sieve::is_probable_prime(&n))});
object_binding!(qs_factor,"num_native_qs_factor","crate::kernels::quadratic_sieve::qs_factor","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(n:BigInt,limits:Handle<QsLimits>)->Option<BigInt>,|ctx:&CallCtx<'_>,n:BigInt,limits:Handle<QsLimits>|{if n.bits()>4096{return Err(error("factor input exceeds bit budget"));}limits.with(ctx,|limits|{sieve_limits(limits)?;Ok(quadratic_sieve::qs_factor(&n,*limits))})});
object_binding!(qs_factor_with,"num_native_qs_factor_with","crate::kernels::quadratic_sieve::qs_factor_with","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(n:BigInt,limits:Handle<QsLimits>,large_primes:bool)->Dynamic,|ctx:&CallCtx<'_>,n:BigInt,limits:Handle<QsLimits>,large_primes:bool|{if n.bits()>4096{return Err(error("factor input exceeds bit budget"));}limits.with(ctx,|limits|{sieve_limits(limits)?;let mut stats=QsStats::default();let factor=quadratic_sieve::qs_factor_with(&n,QsConfig{limits:*limits,large_primes},&mut stats);Ok(Dynamic(api::record(vec![("factor",Output::encode(factor)?),("stats",api::record(vec![("polynomials",Value::int(stats.polynomials)),("candidates",Value::int(stats.candidates)),("full_relations",Value::int(stats.full_relations)),("partial_relations",Value::int(stats.partial_relations)),("combined_relations",Value::int(stats.combined_relations)),("sieve_positions",Value::int(stats.sieve_positions)),("matrix_rows",Value::int(stats.matrix_rows)),("matrix_columns",Value::int(stats.matrix_columns)),("dependencies",Value::int(stats.dependencies))]))])))})});
object_binding!(qs_limits_new,"num_native_qs_limits_new","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(bound:u64,interval:u64,polynomials:usize,large_prime_factor:u64,extra_relations:usize,max_sieve_positions:u64)->NewObject<QsLimits>,|_ctx:&CallCtx<'_>,bound:u64,interval:u64,polynomials:usize,large_prime_factor:u64,extra_relations:usize,max_sieve_positions:u64|{let limits=QsLimits{bound,interval,polynomials,large_prime_factor,extra_relations,max_sieve_positions};sieve_limits(&limits)?;Ok(NewObject::normalized(limits))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &number_exact_factor::BINDING,
            &number_exact_is_prime::BINDING,
            &number_exact_gcd::BINDING,
            &number_exact_lcm::BINDING,
            &number_exact_is_negative::BINDING,
            &number_exact_primes::BINDING,
            &number_exact_binomial::BINDING,
            &number_exact_permutations::BINDING,
            &ordering_amd_order::BINDING,
            &ordering_nested_dissection_order::BINDING,
            &ordering_ml_hybrid_nested_dissection_order::BINDING,
            &ordering_ml_multilevel_nested_dissection_order::BINDING,
            &ordering_cholesky_fill::BINDING,
            &ordering_ml_multilevel_nested_dissection_with::BINDING,
            &ordering_ml_vertex_separator_bisection::BINDING,
            &polynomial_normalized_new::BINDING,
            &real_roots_sturm_sequence::BINDING,
            &real_roots_find_roots::BINDING,
            &real_roots_isolate_real_roots::BINDING,
            &real_roots_refine_root_bisection::BINDING,
            &qs_limits_for_digits::BINDING,
            &qs_is_probable_prime::BINDING,
            &qs_factor::BINDING,
            &qs_factor_with::BINDING,
            &qs_limits_new::BINDING,
        ],
        "numerical_vector",
    )
}
