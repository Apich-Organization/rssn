//! Native electrostatic PIC with explicit SI geometry, charge and transverse permittivity.
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::Exp;
use crate::kernel::signature::Rational;
use crate::kernel::signature::TypeExpr;
use crate::rules::numerical_abi::codecs::Arg;
use crate::rules::numerical_abi::native_object_api::Dynamic;
use crate::rules::numerical_abi::native_object_api::Handle;
use crate::rules::numerical_abi::native_object_api::Input;
use crate::rules::numerical_abi::native_object_api::NativeType;
use crate::rules::numerical_abi::native_object_api::NewObject;
use crate::rules::numerical_abi::native_object_api::object_binding;
use crate::rules::numerical_abi::native_object_api::{
    self as api,
};
use crate::sim::pic as p;
#[derive(Clone)]
struct Frame {
    length: f64,
    time: f64,
    charge: f64,
    mass: f64,
}
fn length() -> DimExpr {
    DimExpr::Base("Length".into())
}
fn time() -> DimExpr {
    DimExpr::Base("Time".into())
}
fn mass() -> DimExpr {
    DimExpr::Base("Mass".into())
}
fn charge() -> DimExpr {
    DimExpr::Mul(vec![DimExpr::Base("Current".into()), time()])
}
fn inv(d: DimExpr) -> DimExpr {
    DimExpr::Inverse(Box::new(d))
}
fn power(
    d: DimExpr,
    n: i64,
) -> DimExpr {
    DimExpr::Pow(Box::new(d), Exp::Rat(Rational::int(n)))
}
fn velocity() -> DimExpr {
    DimExpr::Mul(vec![length(), inv(time())])
}
fn electric() -> DimExpr {
    DimExpr::Mul(vec![mass(), length(), inv(charge()), power(time(), -2)])
}
impl Frame {
    fn new(
        length_ref: &Input,
        time_ref: &Input,
        charge_ref: &Input,
        epsilon_area: &Input,
    ) -> Result<Self, OpError> {
        let (l, ld) = api::reference(&length_ref.0)?;
        let (t, td) = api::reference(&time_ref.0)?;
        let (q, qd) = api::reference(&charge_ref.0)?;
        let (e, ed) = api::reference(&epsilon_area.0)?;
        api::same_dimension(&ld, &length())?;
        api::same_dimension(&td, &time())?;
        api::same_dimension(&qd, &charge())?;
        api::same_dimension(
            &ed,
            &DimExpr::Mul(vec![
                power(charge(), 2),
                power(time(), 2),
                inv(mass()),
                inv(length()),
            ]),
        )?;
        let m = api::scale_product(1.0, &[q, q, t, t], &[e, l])?;
        Ok(Self {
            length: l,
            time: t,
            charge: q,
            mass: m,
        })
    }

    fn e_scale(&self) -> Result<f64, OpError> {
        api::scale_product(
            1.0,
            &[self.mass, self.length],
            &[self.charge, self.time, self.time],
        )
    }
}
#[allow(
    clippy::needless_pass_by_value,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn scalar(
    v: &Input,
    scale: f64,
    dim: DimExpr,
) -> Result<f64, OpError> {
    api::same_dimension(&api::dimension(&v.0)?, &dim)?;
    api::scale_product(
        crate::rules::numerical_abi::numeric(&v.0, 0)?,
        &[],
        &[scale],
    )
}
#[allow(
    clippy::needless_pass_by_value,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn vector(
    v: &Input,
    scale: f64,
    dim: DimExpr,
) -> Result<Vec<f64>, OpError> {
    api::same_dimension(&api::dimension(&v.0)?, &dim)?;
    let values = <Vec<f64> as Arg>::decode(&v.0, 0)?;
    if values.len() > 65536 {
        return Err(api::error("PIC vector budget exceeds65536"));
    }
    values
        .into_iter()
        .map(|v| api::scale_product(v, &[], &[scale]))
        .collect()
}
fn output(
    v: &[f64],
    scale: f64,
    dim: DimExpr,
) -> Result<Dynamic, OpError> {
    Ok(Dynamic(Value::Tensor(crate::kernel::value::TensorValue {
        shape: vec![v.len()],
        elem: TypeExpr::Dim(Box::new(TypeExpr::Float), dim),
        data: v
            .iter()
            .map(|v| api::scale_product(*v, &[scale], &[]).map(Value::Float))
            .collect::<Result<_, _>>()?,
    })))
}
fn shape(name: &str) -> Result<p::Shape, OpError> {
    match name {
        | "Linear" | "CIC" => Ok(p::Shape::Linear),
        | "Cubic" => Ok(p::Shape::Cubic),
        | _ => Err(api::error("PIC shape must beLinear/Cubic")),
    }
}
#[derive(Clone)]
struct Grid {
    data: p::Grid,
    frame: Frame,
}
impl NativeType for Grid {
    const KIND: &'static str = "rssn/native/analysis/pic/grid/v1";
    const MAX_BYTES: usize = 1024;

    fn validate(&self) -> Result<(), OpError> {
        if self.data.cells() > 65536 {
            return Err(api::error("PIC grid cell budget exceeds65536"));
        }
        Ok(())
    }
}
#[derive(Clone)]
struct Config {
    data: p::PicConfig,
    frame: Frame,
}
impl NativeType for Config {
    const KIND: &'static str = "rssn/native/analysis/pic/config/v1";
    const MAX_BYTES: usize = 1024;

    fn validate(&self) -> Result<(), OpError> {
        if self.data.cells > 65536 {
            return Err(api::error("PIC config cell budget exceeds65536"));
        }
        Ok(())
    }
}
struct Pic {
    data: p::Pic,
    frame: Frame,
}
impl NativeType for Pic {
    const KIND: &'static str = "rssn/native/analysis/pic/system/v1";
    const MAX_BYTES: usize = 8 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        if self.data.positions().len() > 65536 || self.data.grid().cells() > 65536 {
            return Err(api::error("PIC particle/grid budget exceeds65536"));
        }
        Ok(())
    }
}
object_binding!(shape_cells,"num_analysis_pic_shape_cells","crate::sim::pic::Shape::cells","Actual native electrostatic PIC operation. SI length/time/charge domains and an explicit positive permittivity-times-transverse-area reference define the native unit-coefficient1D Poisson model. Mass and electric references follow Q²T²/(epsilon_area L) and Q/epsilon_area; no SI constant is guessed.",false,(name:String)->usize,|_ctx,name|Ok(shape(&name)?.cells()));
object_binding!(shape_order,"num_analysis_pic_shape_order","crate::sim::pic::Shape::order","Actual native electrostatic PIC operation. SI length/time/charge domains and an explicit positive permittivity-times-transverse-area reference define the native unit-coefficient1D Poisson model. Mass and electric references follow Q²T²/(epsilon_area L) and Q/epsilon_area; no SI constant is guessed.",false,(name:String)->usize,|_ctx,name|Ok(shape(&name)?.order()));
object_binding!(grid_new,"num_analysis_pic_grid_new","crate::sim::pic::Grid::new","Actual native electrostatic PIC operation. SI length/time/charge domains and an explicit positive permittivity-times-transverse-area reference define the native unit-coefficient1D Poisson model. Mass and electric references follow Q²T²/(epsilon_area L) and Q/epsilon_area; no SI constant is guessed.",false,(cells:usize,box_length:Input,length_ref:Input,time_ref:Input,charge_ref:Input,epsilon_area:Input)->NewObject<Grid>,|_ctx,cells,box_length,length_ref,time_ref,charge_ref,epsilon_area|{if cells>65536{return Err(api::error("PIC cell budget exceeds65536"));}let frame=Frame::new(&length_ref,&time_ref,&charge_ref,&epsilon_area)?;let data=p::Grid::new(cells,scalar(&box_length,frame.length,length())?).map_err(api::error)?;Ok(NewObject::new(Grid{data,frame},vec![length(),time(),charge(),mass()]))});
object_binding!(grid_cells,"num_analysis_pic_grid_cells","crate::sim::pic::Grid::cells","Actual native electrostatic PIC operation. SI length/time/charge domains and an explicit positive permittivity-times-transverse-area reference define the native unit-coefficient1D Poisson model. Mass and electric references follow Q²T²/(epsilon_area L) and Q/epsilon_area; no SI constant is guessed.",false,(grid:Handle<Grid>)->usize,|ctx,grid|grid.with(ctx,|g|Ok(g.data.cells())));
object_binding!(grid_length,"num_analysis_pic_grid_length","crate::sim::pic::Grid::length","Actual native electrostatic PIC operation. SI length/time/charge domains and an explicit positive permittivity-times-transverse-area reference define the native unit-coefficient1D Poisson model. Mass and electric references follow Q²T²/(epsilon_area L) and Q/epsilon_area; no SI constant is guessed.",false,(grid:Handle<Grid>)->Dynamic,|ctx,grid|grid.with(ctx,|g|Ok(Dynamic(api::quantity(api::scale_product(g.data.length(),&[g.frame.length],&[])?,length())?))));
object_binding!(grid_spacing,"num_analysis_pic_grid_spacing","crate::sim::pic::Grid::spacing","Actual native electrostatic PIC operation. SI length/time/charge domains and an explicit positive permittivity-times-transverse-area reference define the native unit-coefficient1D Poisson model. Mass and electric references follow Q²T²/(epsilon_area L) and Q/epsilon_area; no SI constant is guessed.",false,(grid:Handle<Grid>)->Dynamic,|ctx,grid|grid.with(ctx,|g|Ok(Dynamic(api::quantity(api::scale_product(g.data.spacing(),&[g.frame.length],&[])?,length())?))));
object_binding!(deposit,"num_analysis_pic_deposit","crate::sim::pic::Grid::deposit","Actual native electrostatic PIC operation. SI length/time/charge domains and an explicit positive permittivity-times-transverse-area reference define the native unit-coefficient1D Poisson model. Mass and electric references follow Q²T²/(epsilon_area L) and Q/epsilon_area; no SI constant is guessed.",false,(grid:Handle<Grid>,positions:Input,charges:Input,name:String)->Dynamic,|ctx,grid,positions,charges,name|grid.with(ctx,|g|{let positions=vector(&positions,g.frame.length,length())?;let charges=vector(&charges,g.frame.charge,charge())?;if positions.len()!=charges.len(){return Err(api::error("PIC deposit particle shape mismatch"));}let rho=g.data.deposit(&positions,&charges,shape(&name)?);output(&rho,api::scale_product(g.frame.charge,&[],&[g.frame.length])?,DimExpr::Mul(vec![charge(),inv(length())]))}));
object_binding!(gather,"num_analysis_pic_gather","crate::sim::pic::Grid::gather","Actual native electrostatic PIC operation. SI length/time/charge domains and an explicit positive permittivity-times-transverse-area reference define the native unit-coefficient1D Poisson model. Mass and electric references follow Q²T²/(epsilon_area L) and Q/epsilon_area; no SI constant is guessed.",false,(grid:Handle<Grid>,field:Input,positions:Input,name:String,field_reference:Input)->Dynamic,|ctx,grid,field,positions,name,field_reference|grid.with(ctx,|g|{let(scale,dim)=api::reference(&field_reference.0)?;let field=vector(&field,scale,dim.clone())?;let positions=vector(&positions,g.frame.length,length())?;output(&g.data.gather(&field,&positions,shape(&name)?).map_err(api::error)?,scale,dim)}));
object_binding!(solve,"num_analysis_pic_solve","crate::sim::pic::Grid::solve","Actual native electrostatic PIC operation. SI length/time/charge domains and an explicit positive permittivity-times-transverse-area reference define the native unit-coefficient1D Poisson model. Mass and electric references follow Q²T²/(epsilon_area L) and Q/epsilon_area; no SI constant is guessed.",false,(grid:Handle<Grid>,rho:Input,name:String)->Dynamic,|ctx,grid,rho,name|grid.with(ctx,|g|{let rho=vector(&rho,api::scale_product(g.frame.charge,&[],&[g.frame.length])?,DimExpr::Mul(vec![charge(),inv(length())]))?;let r=g.data.solve(&rho,shape(&name)?).map_err(api::error)?;let e=g.frame.e_scale()?;Ok(Dynamic(api::record(vec![("potential",output(&r.phi,api::scale_product(e,&[g.frame.length],&[])?,DimExpr::Mul(vec![electric(),length()]))?.0),("electric_field",output(&r.e,e,electric())?.0)])))}));
object_binding!(config_new,"num_analysis_pic_config_new","crate::sim::pic::PicConfig::new","Actual native electrostatic PIC operation. SI length/time/charge domains and an explicit positive permittivity-times-transverse-area reference define the native unit-coefficient1D Poisson model. Mass and electric references follow Q²T²/(epsilon_area L) and Q/epsilon_area; no SI constant is guessed.",false,(grid:Handle<Grid>,name:String)->NewObject<Config>,|ctx,grid,name|grid.with(ctx,|g|{let data=p::PicConfig::new(g.data.cells(),g.data.length(),shape(&name)?).map_err(api::error)?;Ok(NewObject::new(Config{data,frame:g.frame.clone()},vec![length(),time(),charge(),mass()]))}));
object_binding!(new,"num_analysis_pic_new","crate::sim::pic::Pic::new","Actual native electrostatic PIC operation. SI length/time/charge domains and an explicit positive permittivity-times-transverse-area reference define the native unit-coefficient1D Poisson model. Mass and electric references follow Q²T²/(epsilon_area L) and Q/epsilon_area; no SI constant is guessed.",false,(config:Handle<Config>,positions:Input,velocities:Input,charges:Input,masses:Input)->NewObject<Pic>,|ctx,config,positions,velocities,charges,masses|config.with(ctx,|c|{let frame=&c.frame;let positions=vector(&positions,frame.length,length())?;let velocities=vector(&velocities,api::scale_product(frame.length,&[],&[frame.time])?,velocity())?;let charges=vector(&charges,frame.charge,charge())?;let masses=vector(&masses,frame.mass,mass())?;let data=p::Pic::new(&c.data,positions,velocities,charges,masses).map_err(api::error)?;Ok(NewObject::new(Pic{data,frame:frame.clone()},vec![length(),time(),charge(),mass()]))}));
object_binding!(step,"num_analysis_pic_step","crate::sim::pic::Pic::step","Actual native electrostatic PIC operation. SI length/time/charge domains and an explicit positive permittivity-times-transverse-area reference define the native unit-coefficient1D Poisson model. Mass and electric references follow Q²T²/(epsilon_area L) and Q/epsilon_area; no SI constant is guessed.",false,(pic:Handle<Pic>,h:Input)->(),|ctx,pic,h|pic.with_mut(ctx,|p|{let h=scalar(&h,p.frame.time,time())?;p.data.step(h).map_err(api::error)}));
object_binding!(positions,"num_analysis_pic_positions","crate::sim::pic::Pic::positions","Actual native electrostatic PIC operation. SI length/time/charge domains and an explicit positive permittivity-times-transverse-area reference define the native unit-coefficient1D Poisson model. Mass and electric references follow Q²T²/(epsilon_area L) and Q/epsilon_area; no SI constant is guessed.",false,(pic:Handle<Pic>)->Dynamic,|ctx,pic|pic.with(ctx,|p|output(p.data.positions(),p.frame.length,length())));
object_binding!(velocities,"num_analysis_pic_velocities","crate::sim::pic::Pic::velocities","Actual native electrostatic PIC operation. SI length/time/charge domains and an explicit positive permittivity-times-transverse-area reference define the native unit-coefficient1D Poisson model. Mass and electric references follow Q²T²/(epsilon_area L) and Q/epsilon_area; no SI constant is guessed.",false,(pic:Handle<Pic>)->Dynamic,|ctx,pic|pic.with(ctx,|p|output(p.data.velocities(),api::scale_product(p.frame.length,&[],&[p.frame.time])?,velocity())));
object_binding!(field,"num_analysis_pic_field","crate::sim::pic::Pic::field","Actual native electrostatic PIC operation. SI length/time/charge domains and an explicit positive permittivity-times-transverse-area reference define the native unit-coefficient1D Poisson model. Mass and electric references follow Q²T²/(epsilon_area L) and Q/epsilon_area; no SI constant is guessed.",false,(pic:Handle<Pic>)->Dynamic,|ctx,pic|pic.with(ctx,|p|output(p.data.field(),p.frame.e_scale()?,electric())));
object_binding!(grid,"num_analysis_pic_grid","crate::sim::pic::Pic::grid","Actual native electrostatic PIC operation. SI length/time/charge domains and an explicit positive permittivity-times-transverse-area reference define the native unit-coefficient1D Poisson model. Mass and electric references follow Q²T²/(epsilon_area L) and Q/epsilon_area; no SI constant is guessed.",false,(pic:Handle<Pic>)->NewObject<Grid>,|ctx,pic|pic.with(ctx,|p|Ok(NewObject::new(Grid{data:*p.data.grid(),frame:p.frame.clone()},vec![length(),time(),charge(),mass()]))));
object_binding!(shape,"num_analysis_pic_shape","crate::sim::pic::Pic::shape","Actual native electrostatic PIC operation. SI length/time/charge domains and an explicit positive permittivity-times-transverse-area reference define the native unit-coefficient1D Poisson model. Mass and electric references follow Q²T²/(epsilon_area L) and Q/epsilon_area; no SI constant is guessed.",false,(pic:Handle<Pic>)->String,|ctx,pic|pic.with(ctx,|p|Ok(format!("{:?}",p.data.shape()))));
object_binding!(total_charge,"num_analysis_pic_total_charge","crate::sim::pic::Pic::total_charge","Actual native electrostatic PIC operation. SI length/time/charge domains and an explicit positive permittivity-times-transverse-area reference define the native unit-coefficient1D Poisson model. Mass and electric references follow Q²T²/(epsilon_area L) and Q/epsilon_area; no SI constant is guessed.",false,(pic:Handle<Pic>)->Dynamic,|ctx,pic|pic.with(ctx,|p|Ok(Dynamic(api::quantity(api::scale_product(p.data.total_charge(),&[p.frame.charge],&[])?,charge())?))));
object_binding!(total_momentum,"num_analysis_pic_total_momentum","crate::sim::pic::Pic::total_momentum","Actual native electrostatic PIC operation. SI length/time/charge domains and an explicit positive permittivity-times-transverse-area reference define the native unit-coefficient1D Poisson model. Mass and electric references follow Q²T²/(epsilon_area L) and Q/epsilon_area; no SI constant is guessed.",false,(pic:Handle<Pic>)->Dynamic,|ctx,pic|pic.with(ctx,|p|Ok(Dynamic(api::quantity(api::scale_product(p.data.total_momentum(),&[api::scale_product(1.0,&[p.frame.mass,p.frame.length],&[p.frame.time])?],&[])?,DimExpr::Mul(vec![mass(),velocity()]))?))));
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &shape_cells::BINDING,
            &shape_order::BINDING,
            &grid_new::BINDING,
            &grid_cells::BINDING,
            &grid_length::BINDING,
            &grid_spacing::BINDING,
            &deposit::BINDING,
            &gather::BINDING,
            &solve::BINDING,
            &config_new::BINDING,
            &new::BINDING,
            &step::BINDING,
            &positions::BINDING,
            &velocities::BINDING,
            &field::BINDING,
            &grid::BINDING,
            &shape::BINDING,
            &total_charge::BINDING,
            &total_momentum::BINDING,
        ],
        "numerical_vector",
    )
}
