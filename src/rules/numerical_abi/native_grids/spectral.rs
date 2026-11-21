//! Native Fourier/Chebyshev discretizations and normalized spectral coefficients.
use super::DimExpr;
use super::Dynamic;
use super::Exp;
use super::Faces;
use super::Frame;
use super::Handle;
use super::Input;
use super::Installer;
use super::LIMIT;
use super::NativeType;
use super::NewObject;
use super::OpError;
use super::Output;
use super::Rational;
use super::RuleError;
use super::TensorValue;
use super::TypeExpr;
use super::Value;
use super::api;
use super::canonical;
use super::err;
use super::grid_size;
use super::object_binding;
use super::payload;
use super::structured_ops::checked_snapshot;
use super::structured_ops::dims_forcing;
use super::structured_ops::space;
use super::structured_ops::st;
use crate::kernel::CallCtx;
use crate::rules::numerical_abi::numeric;
use crate::sim::spectral as sp;
use num_complex::Complex64 as C;
#[derive(Clone)]
struct Fourier {
    data: sp::SpectralGrid,
    frame: Frame,
}
impl NativeType for Fourier {
    const KIND: &'static str = "rssn/native/fourier-grid/v1";
    const MAX_BYTES: usize = 4096;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        if self.data.size() > LIMIT {
            return Err(err("Fourier grid allocation budget exceeded"));
        }
        Ok(())
    }
}
#[derive(Clone)]
struct Cheb {
    data: sp::ChebGrid,
    bcs: Faces,
}
impl NativeType for Cheb {
    const KIND: &'static str = "rssn/native/chebyshev-grid/v1";
    const MAX_BYTES: usize = 131072;

    fn validate(&self) -> Result<(), OpError> {
        self.bcs.validate()?;
        if self.data.size() > 512 {
            return Err(err("Chebyshev collocation work budget exceeded"));
        }
        Ok(())
    }
}
#[derive(Clone)]
struct Spectrum {
    data: Vec<C>,
    extents: Vec<usize>,
    frame: Frame,
}
impl NativeType for Spectrum {
    const KIND: &'static str = "rssn/native/spectrum/v1";
    const MAX_BYTES: usize = LIMIT * 16 + 4096;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        if self.data.len() > LIMIT
            || self
                .extents
                .iter()
                .try_fold(1usize, |n, m| n.checked_mul(*m))
                != Some(self.data.len())
            || self
                .data
                .iter()
                .any(|z| !z.re.is_finite() || !z.im.is_finite())
        {
            return Err(err("invalid finite spectrum shape or allocation budget"));
        }
        Ok(())
    }
}
fn spectral_match(
    g: &Fourier,
    s: &Spectrum,
) -> Result<(), OpError> {
    api::same_dimension(&g.frame.dims[0], &s.frame.dims[0])?;
    if s.extents != g.data.n[..g.data.dim] || g.frame.length.to_bits() != s.frame.length.to_bits() {
        return Err(err(
            "spectrum geometry reference does not match Fourier grid",
        ));
    }
    Ok(())
}
fn spectral_payload(s: &Spectrum) -> Result<Dynamic, OpError> {
    Ok(Dynamic(api::record(vec![
        ("shape", s.extents.clone().encode()?),
        ("data", s.data.clone().encode()?),
        (
            "reference",
            api::quantity(s.frame.value, s.frame.dims[2].clone())?,
        ),
    ])))
}
fn reciprocal_length(
    f: &Frame,
    power: i64,
) -> Result<DimExpr, OpError> {
    canonical(DimExpr::Pow(
        Box::new(f.dims[0].clone()),
        Exp::Rat(Rational::new(-power, 1)),
    ))
}

struct FourierForcing(Input);

impl api::NativeArg for FourierForcing {
    fn ty() -> TypeExpr {
        api::tensor_type(2, TypeExpr::Float.with_dim(DimExpr::Var("Forcing".into())))
    }

    fn decode(
        value: &Value,
        index: usize,
        ctx: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        <Input as api::NativeArg>::decode(value, index, ctx).map(Self)
    }
}

struct FourierShift(Input);

impl api::NativeArg for FourierShift {
    fn ty() -> TypeExpr {
        TypeExpr::Float.with_dim(DimExpr::Base("Length".into()).pow(Exp::Rat(Rational::int(-2))))
    }

    fn decode(
        value: &Value,
        index: usize,
        ctx: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        <Input as api::NativeArg>::decode(value, index, ctx).map(Self)
    }
}

object_binding!(
    spectral_fourier_helmholtz_batch,
    "num_grid_fourier_helmholtz_batch",
    "crate::sim::spectral::PeriodicHelmholtzPlan::solve_batch",
    "Solves (-Laplace + sigma) u = f on a periodic Fourier grid for one or more forcing rows in input order. Give measured forcing values in field/length^2 and sigma in inverse length squared. Each row has exactly grid.size points with the x axis fastest. The result reports the batch size followed by the Fourier grid extents, flattened solutions divided by value_reference, and the grid length, time, and field references. Singular compatible Fourier modes use zero amplitude.",
    false,
    (grid: Handle<Fourier>, sources: FourierForcing, sigma: FourierShift) -> Dynamic,
    |ctx, grid, sources, sigma| grid.with(ctx, |g| {
        let frame = &g.frame;
        let FourierForcing(Input(sources)) = sources;
        let FourierShift(Input(sigma)) = sigma;
        let shift_dimension = reciprocal_length(frame, 2)?;
        let forcing_dimension = canonical(DimExpr::Mul(vec![
            frame.dims[2].clone(),
            shift_dimension.clone(),
        ]))?;
        api::same_dimension(&api::dimension(&sources)?, &forcing_dimension)?;
        api::same_dimension(&api::dimension(&sigma)?, &shift_dimension)?;
        let Value::Tensor(input) = &sources else {
            return Err(err("expected a rank-two spectral forcing tensor"));
        };
        let [batch, points] = input.shape.as_slice() else {
            return Err(err("expected a rank-two spectral forcing tensor"));
        };
        if *batch == 0 || *points != g.data.size()
            || batch.checked_mul(*points).is_none_or(|n| n > LIMIT || n != input.data.len())
        {
            return Err(err("periodic Helmholtz forcing shape or batch budget exceeded"));
        }
        let sources = input.data.chunks_exact(*points).map(|row| {
            row.iter().map(|value| {
                api::scale_without_underflow(
                    numeric(value, 1)?, &[frame.length, frame.length], &[frame.value],
                )
            }).collect::<Result<Vec<_>, _>>()
        }).collect::<Result<Vec<_>, _>>()?;
        let sigma = api::scale_without_underflow(
            numeric(&sigma, 2)?, &[frame.length, frame.length], &[],
        )?;
        let solution = g.data.helmholtz_plan(sigma).map_err(err)?.solve_batch(&sources).map_err(err)?;
        let mut shape = vec![*batch];
        shape.extend_from_slice(&g.data.n[..g.data.dim]);
        payload(frame, solution.into_iter().flatten().collect(), shape)
    })
);
object_binding!(spectral_fourier_new,"num_grid_fourier_new","crate::sim::spectral::SpectralGrid::new","Native spectral discretization with explicit coordinate/time/field references. Transforms retain the field domain; wave numbers and derivatives retain inverse-length powers. Native sizes, callback domains and solve work are checked.",true,(frame:Handle<Frame>,extents:Vec<usize>,lengths:Vec<f64>)->NewObject<Fourier>,|ctx,frame,extents,lengths|{let _=ctx;grid_size(&extents)?;if lengths.len()!=extents.len()||lengths.iter().any(|x|*x<=0.0){return Err(err("Fourier lengths must match positive active axes"));}frame.with(ctx,|f|Ok(NewObject::new(Fourier{data:sp::SpectralGrid::new(&extents,&lengths).map_err(err)?,frame:f.clone()},f.dims.to_vec())))});
object_binding!(spectral_fourier_size,"num_grid_fourier_size","crate::sim::spectral::SpectralGrid::size","Native spectral discretization with explicit coordinate/time/field references. Transforms retain the field domain; wave numbers and derivatives retain inverse-length powers. Native sizes, callback domains and solve work are checked.",true,(grid:Handle<Fourier>)->usize,|ctx,grid|{let _=ctx;grid.with(ctx,|g|Ok(g.data.size()))});
object_binding!(spectral_fourier_h,"num_grid_fourier_h","crate::sim::spectral::SpectralGrid::h","Native spectral discretization with explicit coordinate/time/field references. Transforms retain the field domain; wave numbers and derivatives retain inverse-length powers. Native sizes, callback domains and solve work are checked.",true,(grid:Handle<Fourier>,axis:usize)->Dynamic,|ctx,grid,axis|{let _=ctx;grid.with(ctx,|g|{if axis>=g.data.dim{return Err(err("Fourier axis is inactive"));}Ok(Dynamic(api::quantity(g.data.h(axis)*g.frame.length,g.frame.dims[0].clone())?))})});
object_binding!(spectral_fourier_coord,"num_grid_fourier_coord","crate::sim::spectral::SpectralGrid::coord","Native spectral discretization with explicit coordinate/time/field references. Transforms retain the field domain; wave numbers and derivatives retain inverse-length powers. Native sizes, callback domains and solve work are checked.",true,(grid:Handle<Fourier>,index:usize)->Dynamic,|ctx,grid,index|{let _=ctx;grid.with(ctx,|g|{if index>=g.data.size(){return Err(err("Fourier node index out of range"));}let data=g.data.coord(index).into_iter().map(|v|crate::rules::numerical_abi::output_scalar(v*g.frame.length)).collect::<Result<_,_>>()?;Ok(Dynamic(Value::Tensor(TensorValue{shape:vec![3],elem:TypeExpr::Float.with_dim(g.frame.dims[0].clone()),data})))})});
object_binding!(spectral_fourier_k2,"num_grid_fourier_k2","crate::sim::spectral::SpectralGrid::k2","Native spectral discretization with explicit coordinate/time/field references. Transforms retain the field domain; wave numbers and derivatives retain inverse-length powers. Native sizes, callback domains and solve work are checked.",true,(grid:Handle<Fourier>,index:usize)->Dynamic,|ctx,grid,index|{let _=ctx;grid.with(ctx,|g|{if index>=g.data.size(){return Err(err("Fourier index out of range"));}Ok(Dynamic(api::quantity(g.data.k2(index)/(g.frame.length*g.frame.length),reciprocal_length(&g.frame,2)?)?))})});
object_binding!(spectral_fourier_kvec,"num_grid_fourier_kvec","crate::sim::spectral::SpectralGrid::kvec","Native spectral discretization with explicit coordinate/time/field references. Transforms retain the field domain; wave numbers and derivatives retain inverse-length powers. Native sizes, callback domains and solve work are checked.",true,(grid:Handle<Fourier>,index:usize)->Dynamic,|ctx,grid,index|{let _=ctx;grid.with(ctx,|g|{if index>=g.data.size(){return Err(err("Fourier index out of range"));}let data=g.data.kvec(index).into_iter().map(|v|crate::rules::numerical_abi::output_scalar(v/g.frame.length)).collect::<Result<_,_>>()?;Ok(Dynamic(Value::Tensor(TensorValue{shape:vec![3],elem:TypeExpr::Float.with_dim(reciprocal_length(&g.frame,1)?),data})))})});
object_binding!(spectral_fourier_wavenumber,"num_grid_fourier_wavenumber","crate::sim::spectral::SpectralGrid::wavenumber","Native spectral discretization with explicit coordinate/time/field references. Transforms retain the field domain; wave numbers and derivatives retain inverse-length powers. Native sizes, callback domains and solve work are checked.",true,(grid:Handle<Fourier>,axis:usize,index:usize)->Dynamic,|ctx,grid,axis,index|{let _=ctx;grid.with(ctx,|g|{if axis>=g.data.dim||index>=g.data.n[axis]{return Err(err("Fourier frequency index out of range"));}Ok(Dynamic(api::quantity(g.data.wavenumber(axis,index)/g.frame.length,reciprocal_length(&g.frame,1)?)?))})});
object_binding!(spectral_fourier_dealias_mask,"num_grid_fourier_dealias_mask","crate::sim::spectral::SpectralGrid::dealias_mask","Native spectral discretization with explicit coordinate/time/field references. Transforms retain the field domain; wave numbers and derivatives retain inverse-length powers. Native sizes, callback domains and solve work are checked.",true,(grid:Handle<Fourier>)->Vec<bool>,|ctx,grid|{let _=ctx;grid.with(ctx,|g|Ok(g.data.dealias_mask()))});
object_binding!(spectral_fourier_to_spec,"num_grid_fourier_to_spec","crate::sim::spectral::SpectralGrid::to_spec","Native spectral discretization with explicit coordinate/time/field references. Transforms retain the field domain; wave numbers and derivatives retain inverse-length powers. Native sizes, callback domains and solve work are checked.",true,(grid:Handle<Fourier>,values:Vec<f64>)->NewObject<Spectrum>,|ctx,grid,values|{let _=ctx;grid.with(ctx,|g|{if values.len()!=g.data.size(){return Err(err("Fourier input length mismatch"));}Ok(NewObject::new(Spectrum{data:g.data.to_spec(&values),extents:g.data.n[..g.data.dim].to_vec(),frame:g.frame.clone()},g.frame.dims.to_vec()))})});
object_binding!(spectral_fourier_to_spec_c,"num_grid_fourier_to_spec_c","crate::sim::spectral::SpectralGrid::to_spec_c","Native spectral discretization with explicit coordinate/time/field references. Transforms retain the field domain; wave numbers and derivatives retain inverse-length powers. Native sizes, callback domains and solve work are checked.",true,(grid:Handle<Fourier>,values:Vec<C>)->NewObject<Spectrum>,|ctx,grid,values|{let _=ctx;grid.with(ctx,|g|{if values.len()!=g.data.size(){return Err(err("Fourier input length mismatch"));}Ok(NewObject::new(Spectrum{data:g.data.to_spec_c(&values),extents:g.data.n[..g.data.dim].to_vec(),frame:g.frame.clone()},g.frame.dims.to_vec()))})});
object_binding!(spectral_fourier_to_phys,"num_grid_fourier_to_phys","crate::sim::spectral::SpectralGrid::to_phys","Native spectral discretization with explicit coordinate/time/field references. Transforms retain the field domain; wave numbers and derivatives retain inverse-length powers. Native sizes, callback domains and solve work are checked.",true,(grid:Handle<Fourier>,spectrum:Handle<Spectrum>)->Dynamic,|ctx,grid,spectrum|{let _=ctx;grid.with(ctx,|g|spectrum.with(ctx,|s|{spectral_match(g,s)?;let data=g.data.to_phys(&s.data);payload(&s.frame,data,s.extents.clone())}))});
object_binding!(spectral_fourier_to_phys_c,"num_grid_fourier_to_phys_c","crate::sim::spectral::SpectralGrid::to_phys_c","Native spectral discretization with explicit coordinate/time/field references. Transforms retain the field domain; wave numbers and derivatives retain inverse-length powers. Native sizes, callback domains and solve work are checked.",true,(grid:Handle<Fourier>,spectrum:Handle<Spectrum>)->Dynamic,|ctx,grid,spectrum|{let _=ctx;grid.with(ctx,|g|spectrum.with(ctx,|s|{spectral_match(g,s)?;let data=g.data.to_phys_c(&s.data);Ok(Dynamic(api::record(vec![("shape",s.extents.clone().encode()?),("data",data.encode()?),("reference",api::quantity(s.frame.value,s.frame.dims[2].clone())?)])))}))});
object_binding!(spectral_fourier_deriv,"num_grid_fourier_deriv","crate::sim::spectral::SpectralGrid::deriv","Native spectral discretization with explicit coordinate/time/field references. Transforms retain the field domain; wave numbers and derivatives retain inverse-length powers. Native sizes, callback domains and solve work are checked.",true,(grid:Handle<Fourier>,spectrum:Handle<Spectrum>,axis:usize)->NewObject<Spectrum>,|ctx,grid,spectrum,axis|{let _=ctx;grid.with(ctx,|g|spectrum.with(ctx,|s|{spectral_match(g,s)?;if axis>=g.data.dim{return Err(err("Fourier derivative axis is inactive"));}let mut frame=s.frame.clone();frame.value/=frame.length;frame.dims[2]=canonical(DimExpr::Mul(vec![frame.dims[2].clone(),DimExpr::Inverse(Box::new(frame.dims[0].clone()))]))?;Ok(NewObject::new(Spectrum{data:g.data.deriv(&s.data,axis),extents:s.extents.clone(),frame:frame.clone()},frame.dims.to_vec()))}))});
object_binding!(spectral_fourier_helmholtz,"num_grid_fourier_helmholtz","crate::sim::spectral::SpectralGrid::helmholtz","Native spectral discretization with explicit coordinate/time/field references. Transforms retain the field domain; wave numbers and derivatives retain inverse-length powers. Native sizes, callback domains and solve work are checked.",true,(grid:Handle<Fourier>,source:Vec<f64>,sigma:f64)->Dynamic,|ctx,grid,source,sigma|{let _=ctx;grid.with(ctx,|g|{if source.len()!=g.data.size(){return Err(err("spectral forcing length mismatch"));}payload(&g.frame,g.data.helmholtz(sigma,&source).map_err(err)?,g.data.n[..g.data.dim].to_vec())})});
object_binding!(spectral_spectrum_data,"num_grid_spectrum_data","","Native spectral discretization with explicit coordinate/time/field references. Transforms retain the field domain; wave numbers and derivatives retain inverse-length powers. Native sizes, callback domains and solve work are checked.",true,(spectrum:Handle<Spectrum>)->Dynamic,|ctx,spectrum|{let _=ctx;spectrum.with(ctx,spectral_payload)});
object_binding!(spectral_fft_nd,"num_grid_fft_nd","crate::sim::spectral::fft_nd","Native spectral discretization with explicit coordinate/time/field references. Transforms retain the field domain; wave numbers and derivatives retain inverse-length powers. Native sizes, callback domains and solve work are checked.",true,(spectrum:Handle<Spectrum>,inverse:bool)->NewObject<Spectrum>,|ctx,spectrum,inverse|{let _=ctx;spectrum.with(ctx,|s|{let mut out=s.clone();sp::fft_nd(&mut out.data,&out.extents,inverse).map_err(err)?;Ok(NewObject::new(out,s.frame.dims.to_vec()))})});
object_binding!(spectral_chebyshev_new,"num_grid_chebyshev_new","crate::sim::spectral::ChebGrid::new","Native spectral discretization with explicit coordinate/time/field references. Transforms retain the field domain; wave numbers and derivatives retain inverse-length powers. Native sizes, callback domains and solve work are checked.",true,(bcs:Handle<Faces>,degrees:Vec<usize>,lo:Vec<f64>,hi:Vec<f64>)->NewObject<Cheb>,|ctx,bcs,degrees,lo,hi|{let _=ctx;if degrees.is_empty()||degrees.len()>3||degrees.iter().any(|d|*d<2||*d>128)||degrees.iter().try_fold(1usize,|n,d|n.checked_mul(d+1)).is_none_or(|n|n>512)||lo.len()!=degrees.len()||hi.len()!=degrees.len()||lo.iter().zip(&hi).any(|(l,h)|l>=h){return Err(err("Chebyshev degree/bounds/work budget invalid"));}bcs.with(ctx,|b|Ok(NewObject::dependent(Cheb{data:sp::ChebGrid::new(&degrees,&lo,&hi).map_err(err)?,bcs:b.clone()},b.frame.dims.to_vec(),vec![bcs.value.clone()])))});
object_binding!(spectral_chebyshev_size,"num_grid_chebyshev_size","crate::sim::spectral::ChebGrid::size","Native spectral discretization with explicit coordinate/time/field references. Transforms retain the field domain; wave numbers and derivatives retain inverse-length powers. Native sizes, callback domains and solve work are checked.",true,(grid:Handle<Cheb>)->usize,|ctx,grid|{let _=ctx;grid.with(ctx,|g|Ok(g.data.size()))});
object_binding!(spectral_chebyshev_node,"num_grid_chebyshev_node","crate::sim::spectral::ChebGrid::node","Native spectral discretization with explicit coordinate/time/field references. Transforms retain the field domain; wave numbers and derivatives retain inverse-length powers. Native sizes, callback domains and solve work are checked.",true,(grid:Handle<Cheb>,axis:usize,index:usize)->Dynamic,|ctx,grid,axis,index|{let _=ctx;grid.with(ctx,|g|{if axis>=g.data.dim||index>g.data.deg[axis]{return Err(err("Chebyshev node index out of range"));}Ok(Dynamic(api::quantity(g.data.node(axis,index)*g.bcs.frame.length,g.bcs.frame.dims[0].clone())?))})});
object_binding!(spectral_chebyshev_coord,"num_grid_chebyshev_coord","crate::sim::spectral::ChebGrid::coord","Native spectral discretization with explicit coordinate/time/field references. Transforms retain the field domain; wave numbers and derivatives retain inverse-length powers. Native sizes, callback domains and solve work are checked.",true,(grid:Handle<Cheb>,index:usize)->Dynamic,|ctx,grid,index|{let _=ctx;grid.with(ctx,|g|{if index>=g.data.size(){return Err(err("Chebyshev node out of range"));}let data=g.data.coord(index).into_iter().map(|v|crate::rules::numerical_abi::output_scalar(v*g.bcs.frame.length)).collect::<Result<_,_>>()?;Ok(Dynamic(Value::Tensor(TensorValue{shape:vec![3],elem:TypeExpr::Float.with_dim(g.bcs.frame.dims[0].clone()),data})))})});
object_binding!(spectral_chebyshev_sample,"num_grid_chebyshev_sample","crate::sim::spectral::ChebGrid::sample","Native spectral discretization with explicit coordinate/time/field references. Transforms retain the field domain; wave numbers and derivatives retain inverse-length powers. Native sizes, callback domains and solve work are checked.",true,(grid:Handle<Cheb>,program:Input,budget:usize)->Dynamic,|ctx,grid,program,budget|{let _=ctx;grid.with(ctx,|g|{let frame=&g.bcs.frame;let domain=frame.dims[2].clone();let s=checked_snapshot(ctx,&program,budget,&vec![frame.dims[0].clone();3],&domain)?;let source=|x|space(&s,frame,x,frame.value);let out=g.data.sample(source);s.finish()?;payload(frame,out,vec![g.data.size()])})});
object_binding!(spectral_chebyshev_max_error,"num_grid_chebyshev_max_error","crate::sim::spectral::ChebGrid::max_error","Native spectral discretization with explicit coordinate/time/field references. Transforms retain the field domain; wave numbers and derivatives retain inverse-length powers. Native sizes, callback domains and solve work are checked.",true,(grid:Handle<Cheb>,program:Input,budget:usize,values:Vec<f64>)->Dynamic,|ctx,grid,program,budget,values|{let _=ctx;grid.with(ctx,|g|{let frame=&g.bcs.frame;let domain=frame.dims[2].clone();let s=checked_snapshot(ctx,&program,budget,&vec![frame.dims[0].clone();3],&domain)?;let source=|x|space(&s,frame,x,frame.value);if values.len()!=g.data.size(){return Err(err("Chebyshev field length mismatch"));}let out=g.data.max_error(&values,source);s.finish()?;Ok(Dynamic(api::quantity(out*frame.value,frame.dims[2].clone())?))})});
object_binding!(spectral_chebyshev_solve_helmholtz,"num_grid_chebyshev_solve_helmholtz","crate::sim::spectral::ChebGrid::solve_helmholtz","Native spectral discretization with explicit coordinate/time/field references. Transforms retain the field domain; wave numbers and derivatives retain inverse-length powers. Native sizes, callback domains and solve work are checked.",true,(grid:Handle<Cheb>,program:Input,budget:usize,sigma:f64)->Dynamic,|ctx,grid,program,budget,sigma|{let _=ctx;grid.with(ctx,|g|{let frame=&g.bcs.frame;let frame_reference_m2_0_1=frame.reference([-2, 0, 1])?;let domain=canonical(DimExpr::Mul(vec![frame.dims[2].clone(),reciprocal_length(frame,2)?]))?;let s=checked_snapshot(ctx,&program,budget,&vec![frame.dims[0].clone();3],&domain)?;let source=|x|space(&s,frame,x,frame_reference_m2_0_1);let out=g.data.solve_helmholtz(&g.bcs.data,sigma,source).map_err(err)?;s.finish()?;g.bcs.finish_for(g.data.dim)?;payload(frame,out,vec![g.data.size()])})});
object_binding!(spectral_chebyshev_heat,"num_grid_chebyshev_heat","crate::sim::spectral::ChebGrid::heat","Native spectral discretization with explicit coordinate/time/field references. Transforms retain the field domain; wave numbers and derivatives retain inverse-length powers. Native sizes, callback domains and solve work are checked.",true,(grid:Handle<Cheb>,initial:Input,source:Input,budget:usize,alpha:f64,t_end:f64,steps:usize,theta:f64)->Dynamic,|ctx,grid,initial,source,budget,alpha,t_end,steps,theta|{let _=ctx;grid.with(ctx,|g|{let f=&g.bcs.frame;let f_reference_0_m1_1=if matches!(&source.0,Value::Variant(v)if v.tag=="None"&&v.payload.is_empty()){1.0}else{f.reference([0, -1, 1])?};
if steps==0||steps>8192||alpha<0.0||t_end<=0.0||!(0.0..=1.0).contains(&theta)||g.data.size().checked_mul(steps).is_none_or(|n|n>1_048_576){return Err(err("Chebyshev heat control/work budget exceeded"));}let init=checked_snapshot(ctx,&initial,budget,&vec![f.dims[0].clone();3],&f.dims[2])?;let src=if matches!(&source.0,Value::Variant(v)if v.tag=="None"&&v.payload.is_empty()){None}else{Some(checked_snapshot(ctx,&source,budget,&[f.dims[1].clone(),f.dims[0].clone(),f.dims[0].clone(),f.dims[0].clone()],&dims_forcing(f,1)?)?)};let initfn=|x|space(&init,f,x,f.value);let sourcefn=|t,x|src.as_ref().map_or(0.0,|s|st(s,f,t,x,f_reference_0_m1_1));let out=g.data.heat(&g.bcs.data,alpha,initfn,Some(&sourcefn),t_end,steps,theta).map_err(err)?;init.finish()?;if let Some(s)=src{s.finish()?;}g.bcs.finish_for(g.data.dim)?;payload(f,out,vec![g.data.size()])})});
object_binding!(spectral_cheb_diff_matrix,"num_grid_cheb_diff_matrix","crate::sim::spectral::cheb_diff_matrix","Native spectral discretization with explicit coordinate/time/field references. Transforms retain the field domain; wave numbers and derivatives retain inverse-length powers. Native sizes, callback domains and solve work are checked.",true,(degree:usize,length:f64)->Vec<Vec<f64>>,|ctx,degree,length|{let _=ctx;if degree==0||degree>128||length<=0.0{return Err(err("Chebyshev derivative degree/length budget invalid"));}Ok(sp::cheb_diff_matrix(degree,length))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(i, BINDINGS, "numerical_vector")
}
static BINDINGS: &[&api::Binding] = &[
    &spectral_etdrk4::BINDING,
    &spectral_fourier_new::BINDING,
    &spectral_fourier_size::BINDING,
    &spectral_fourier_h::BINDING,
    &spectral_fourier_coord::BINDING,
    &spectral_fourier_k2::BINDING,
    &spectral_fourier_kvec::BINDING,
    &spectral_fourier_wavenumber::BINDING,
    &spectral_fourier_dealias_mask::BINDING,
    &spectral_fourier_to_spec::BINDING,
    &spectral_fourier_to_spec_c::BINDING,
    &spectral_fourier_to_phys::BINDING,
    &spectral_fourier_to_phys_c::BINDING,
    &spectral_fourier_deriv::BINDING,
    &spectral_fourier_helmholtz::BINDING,
    &spectral_fourier_helmholtz_batch::BINDING,
    &spectral_spectrum_data::BINDING,
    &spectral_fft_nd::BINDING,
    &spectral_chebyshev_new::BINDING,
    &spectral_chebyshev_size::BINDING,
    &spectral_chebyshev_node::BINDING,
    &spectral_chebyshev_coord::BINDING,
    &spectral_chebyshev_sample::BINDING,
    &spectral_chebyshev_max_error::BINDING,
    &spectral_chebyshev_solve_helmholtz::BINDING,
    &spectral_chebyshev_heat::BINDING,
    &spectral_cheb_diff_matrix::BINDING,
];
object_binding!(spectral_etdrk4,"num_grid_etdrk4","crate::sim::spectral::etdrk4","Native ETDRK4 with a genuinely coupled compiled nonlinear map: two scalar programs per complex mode (real then imaginary), each taking coherent [time,re0,im0,...]. At most31 modes fit the compiled ABI's63 inputs. Linear symbols are normalized inverse-time coefficients; nonlinear outputs must have field/time units. Native x-fastest spectral extent ordering and exact ETDRK4 stepping are preserved.",true,(initial:Handle<Spectrum>,linear:Vec<C>,programs:Input,budget:usize,t_end:f64,steps:usize)->NewObject<Spectrum>,|ctx,initial,linear,programs,budget,t_end,steps|{
 initial.with(ctx,|s|{let s_frame_reference_0_m1_1=s.frame.reference([0, -1, 1])?;let n=s.data.len();if n>31||linear.len()!=n||steps==0||steps>4096||t_end<=0.0||n.checked_mul(steps).is_none_or(|w|w>32768){return Err(err("coupled spectral ETDRK4 shape/work budget exceeded"));}let Value::Tensor(t)=&programs.0 else{return Err(err("expected real/imaginary compiled map programs"));};
if t.shape!=[2*n]{return Err(err("one real and one imaginary program required per mode"));}let mut input_dimensions=vec![s.frame.dims[1].clone()];input_dimensions.extend(vec![s.frame.dims[2].clone();2*n]);let output=dims_forcing(&s.frame,1)?;let mut callbacks=Vec::new();for program in &t.data{callbacks.push(checked_snapshot(ctx,&Input(program.clone()),budget,&input_dimensions,&output)?);}let nonlinear=|time:f64,u:&[C]|{let mut coherent=Vec::with_capacity(1+2*n);coherent.push(time*s.frame.time);for z in u{coherent.push(z.re*s.frame.value);coherent.push(z.im*s.frame.value);}callbacks.as_chunks::<2>().0.iter().map(|pair|C::new(pair[0].sample(&coherent)/(s_frame_reference_0_m1_1),pair[1].sample(&coherent)/(s_frame_reference_0_m1_1))).collect::<Vec<_>>()};let data=sp::etdrk4(s.data.clone(),&linear,&nonlinear,t_end,steps);for callback in callbacks{callback.finish()?;}Ok(NewObject::new(Spectrum{data,extents:s.extents.clone(),frame:s.frame.clone()},s.frame.dims.to_vec()))})
});
