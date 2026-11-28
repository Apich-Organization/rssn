//! Session-native continuum geometry and boundary models, in an explicit reference frame.
use super::codecs::Output;
use super::native_backend::CompiledSnapshot;
use super::native_backend::compiled_snapshot;
use super::native_object_api::Dynamic;
use super::native_object_api::Handle;
use super::native_object_api::Input;
use super::native_object_api::NativeArg;
use super::native_object_api::NativeType;
use super::native_object_api::NewObject;
use super::native_object_api::object_binding;
use super::native_object_api::{
    self as api,
};
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::call::CallCtx;
use crate::kernel::call::OpError;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::Exp;
use crate::kernel::signature::Rational;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::TensorValue;
use crate::kernel::value::Value;
use crate::sim::boundary;
use crate::sim::structured;
const LIMIT: usize = 65536;
fn err(e: impl std::fmt::Display) -> OpError {
    api::error(e)
}
#[allow(
    clippy::needless_pass_by_value,
    reason = "accepts owned temporary dimension expressions built by adapters"
)]
fn canonical(d: DimExpr) -> Result<DimExpr, OpError> {
    crate::rules::numerical_abi::runtime_dimensions::expression(
        &crate::rules::numerical_abi::runtime_dimensions::canonical(&d)?,
    )
}
#[derive(Clone)]
struct Frame {
    length: f64,
    time: f64,
    value: f64,
    dims: [DimExpr; 3],
}
impl NativeType for Frame {
    const KIND: &'static str = "rssn/native/pde-frame/v1";
    const MAX_BYTES: usize = 2048;

    fn validate(&self) -> Result<(), OpError> {
        if [self.length, self.time, self.value]
            .iter()
            .any(|s| !s.is_finite() || *s <= 0.0)
        {
            return Err(err("reference scales must be finite positive"));
        }
        api::same_dimension(&self.dims[0], &DimExpr::Base("Length".into()))?;
        api::same_dimension(&self.dims[1], &DimExpr::Base("Time".into()))?;
        Ok(())
    }
}
impl Frame {
    // Validate only references used by a given operation. Keep factors separate
    // so a representable quotient does not overflow in an intermediate product.
    fn reference(
        &self,
        powers: [i32; 3],
    ) -> Result<f64, OpError> {
        let mut numerator = Vec::new();
        let mut denominator = Vec::new();
        for (scale, exponent) in [self.length, self.time, self.value].into_iter().zip(powers) {
            if exponent.unsigned_abs() > 16 {
                return Err(err("derived reference exponent budget exceeded"));
            }
            let factors = if exponent >= 0 {
                &mut numerator
            } else {
                &mut denominator
            };
            factors.extend(std::iter::repeat_n(scale, exponent.unsigned_abs() as usize));
        }
        let scale = api::scale_product(1.0, &numerator, &denominator)?;
        if scale <= 0.0 {
            return Err(err("derived physical reference underflows Float"));
        }
        Ok(scale)
    }

    fn gradient(&self) -> Result<DimExpr, OpError> {
        canonical(DimExpr::Mul(vec![
            self.dims[2].clone(),
            DimExpr::Inverse(Box::new(self.dims[0].clone())),
        ]))
    }

    fn same(
        &self,
        b: &Self,
    ) -> Result<(), OpError> {
        for (a, b) in self.dims.iter().zip(&b.dims) {
            api::same_dimension(a, b)?;
        }
        if [self.length, self.time, self.value]
            .iter()
            .zip([b.length, b.time, b.value])
            .any(|(a, b)| a.to_bits() != b.to_bits())
        {
            return Err(err("native PDE model reference frames differ"));
        }
        Ok(())
    }
}
fn payload(
    frame: &Frame,
    values: Vec<f64>,
    shape: Vec<usize>,
) -> Result<Dynamic, OpError> {
    if shape.iter().try_fold(1usize, |a, b| a.checked_mul(*b)) != Some(values.len())
        || values.len() > LIMIT
    {
        return Err(err("PDE field shape/allocation budget exceeded"));
    }
    Ok(Dynamic(api::record(vec![
        ("shape", shape.encode()?),
        ("data", values.encode()?),
        (
            "length_reference",
            api::quantity(frame.length, frame.dims[0].clone())?,
        ),
        (
            "time_reference",
            api::quantity(frame.time, frame.dims[1].clone())?,
        ),
        (
            "value_reference",
            api::quantity(frame.value, frame.dims[2].clone())?,
        ),
    ])))
}
#[derive(Clone)]
struct Datum {
    data: boundary::BcValue,
    frame: Frame,
    value_scale: f64,
    value_dimension: DimExpr,
    snapshot: Option<CompiledSnapshot>,
}
impl NativeType for Datum {
    const KIND: &'static str = "rssn/native/bc-value/v1";
    const MAX_BYTES: usize = 4096;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        if !self.value_scale.is_finite() || self.value_scale <= 0.0 {
            return Err(err("invalid boundary output reference"));
        }
        Ok(())
    }
}
impl Datum {
    fn callback_value(
        snapshot: CompiledSnapshot,
        frame: &Frame,
        output_scale: f64,
    ) -> Result<boundary::BcValue, OpError> {
        let length = frame.length;
        let time = frame.time;
        match snapshot.input_dimensions.len() {
            | 1 => {
                Ok(boundary::BcValue::time(move |t| {
                    snapshot.sample(&[t * time]) / output_scale
                }))
            },
            | 3 => {
                Ok(boundary::BcValue::space(move |x| {
                    snapshot.sample(&x.map(|value| value * length)) / output_scale
                }))
            },
            | 4 => {
                Ok(boundary::BcValue::spacetime(move |t, x| {
                    snapshot.sample(&[t * time, x[0] * length, x[1] * length, x[2] * length])
                        / output_scale
                }))
            },
            | _ => Err(err("boundary callback arity mismatch")),
        }
    }

    fn finish(&self) -> Result<(), OpError> {
        self.snapshot
            .as_ref()
            .map_or(Ok(()), CompiledSnapshot::finish)
    }

    fn require(
        &self,
        gradient: bool,
    ) -> Result<(), OpError> {
        let d = if gradient {
            self.frame.gradient()?
        } else {
            self.frame.dims[2].clone()
        };
        api::same_dimension(&d, &self.value_dimension)?;
        let s = if gradient {
            self.frame.reference([-1, 0, 1])?
        } else {
            self.frame.value
        };
        if s.to_bits() != self.value_scale.to_bits() {
            return Err(err(
                "boundary output reference must match field reference (or field/length for Neumann)",
            ));
        }
        Ok(())
    }
}
#[derive(Clone)]
struct Boundary {
    data: boundary::BoundaryCondition,
    frame: Frame,
    datums: Vec<Datum>,
}
impl NativeType for Boundary {
    const KIND: &'static str = "rssn/native/scalar-boundary/v1";
    const MAX_BYTES: usize = 8192;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()
    }
}
#[derive(Clone)]
struct Faces {
    data: boundary::FaceBcs,
    frame: Frame,
    datums: [Vec<Datum>; 6],
}
impl NativeType for Faces {
    const KIND: &'static str = "rssn/native/scalar-faces/v1";
    const MAX_BYTES: usize = 65536;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        if self.datums.iter().map(Vec::len).sum::<usize>() > 64 {
            return Err(err("boundary callback retention budget exceeded"));
        }
        Ok(())
    }
}
impl Faces {
    fn finish_for(
        &self,
        active_axes: usize,
    ) -> Result<(), OpError> {
        for d in self.datums.iter().take(active_axes.min(3) * 2).flatten() {
            d.finish()?;
        }
        Ok(())
    }
}
#[derive(Clone)]
struct Grid {
    data: structured::Grid,
    frame: Frame,
    bcs: Faces,
}
impl NativeType for Grid {
    const KIND: &'static str = "rssn/native/structured-grid/v1";
    const MAX_BYTES: usize = 131072;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        if self.data.len() > LIMIT {
            return Err(err("structured grid allocation budget exceeded"));
        }
        Ok(())
    }
}
fn face(i: usize) -> Result<boundary::Face, OpError> {
    boundary::Face::ALL
        .get(i)
        .copied()
        .ok_or_else(|| err("face index must be zero through five"))
}
fn grid_size(n: &[usize]) -> Result<(), OpError> {
    if n.is_empty()
        || n.len() > 3
        || n.iter().any(|n| *n < 2 || *n > 4096)
        || n.iter()
            .try_fold(1usize, |a, b| a.checked_mul(*b))
            .is_none_or(|n| n > LIMIT)
    {
        return Err(err(
            "grid requires one through three axes within allocation budget",
        ));
    }
    Ok(())
}
object_binding!(grid_frame_new,"num_grid_frame_new","","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",false,(length:Input,time:Input,value:Input)->NewObject<Frame>,|ctx,length,time,value|{let _=ctx;let(l,ld)=api::reference(&length.0)?;let(t,td)=api::reference(&time.0)?;let(v,vd)=api::reference(&value.0)?;let f=Frame{length:l,time:t,value:v,dims:[ld,td,vd]};f.validate()?;Ok(NewObject::new(f.clone(),f.dims.to_vec()))});
object_binding!(grid_frame_references,"num_grid_frame_references","","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(frame:Handle<Frame>)->Dynamic,|ctx,frame|{let _=ctx;frame.with(ctx,|f|Ok(Dynamic(api::record(vec![("length",api::quantity(f.length,f.dims[0].clone())?),("time",api::quantity(f.time,f.dims[1].clone())?),("value",api::quantity(f.value,f.dims[2].clone())?)]))))});
object_binding!(grid_datum_constant,"num_grid_datum_constant","crate::sim::boundary::BcValue::constant","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",false,(frame:Handle<Frame>,value:NormalizedFloat,reference:Input)->NewObject<Datum>,|ctx,frame,value,reference|{let _=ctx;let(s,d)=api::reference(&reference.0)?;frame.with(ctx,|f|Ok(NewObject::new(Datum{data:boundary::BcValue::constant(value.0),frame:f.clone(),value_scale:s,value_dimension:d.clone(),snapshot:None},vec![d])))});
object_binding!(grid_datum_time,"num_grid_datum_time","crate::sim::boundary::BcValue::time","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",false,(frame:Handle<Frame>,program:Input,budget:usize,reference:Input)->NewObject<Datum>,|ctx,frame,program,budget,reference|{let _=ctx;let snapshot=compiled_snapshot(ctx,&program.0,budget)?;let(s,d)=api::reference(&reference.0)?;frame.with(ctx,|f|{let expected=vec![f.dims[1].clone()];if snapshot.input_dimensions.len()!=1{return Err(err("boundary callback arity mismatch"));}for(a,b)in snapshot.input_dimensions.iter().zip(expected){api::same_dimension(a,&b)?;}api::same_dimension(&snapshot.output_dimension,&d)?;let data=Datum::callback_value(snapshot.clone(),f,s)?;Ok(NewObject::dependent(Datum{data,frame:f.clone(),value_scale:s,value_dimension:d.clone(),snapshot:Some(snapshot)},vec![d],vec![program.0.clone()]))})});
object_binding!(grid_datum_space,"num_grid_datum_space","crate::sim::boundary::BcValue::space","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",false,(frame:Handle<Frame>,program:Input,budget:usize,reference:Input)->NewObject<Datum>,|ctx,frame,program,budget,reference|{let _=ctx;let snapshot=compiled_snapshot(ctx,&program.0,budget)?;let(s,d)=api::reference(&reference.0)?;frame.with(ctx,|f|{let expected=vec![f.dims[0].clone();3];if snapshot.input_dimensions.len()!=3{return Err(err("boundary callback arity mismatch"));}for(a,b)in snapshot.input_dimensions.iter().zip(expected){api::same_dimension(a,&b)?;}api::same_dimension(&snapshot.output_dimension,&d)?;let data=Datum::callback_value(snapshot.clone(),f,s)?;Ok(NewObject::dependent(Datum{data,frame:f.clone(),value_scale:s,value_dimension:d.clone(),snapshot:Some(snapshot)},vec![d],vec![program.0.clone()]))})});
object_binding!(grid_datum_spacetime,"num_grid_datum_spacetime","crate::sim::boundary::BcValue::spacetime","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",false,(frame:Handle<Frame>,program:Input,budget:usize,reference:Input)->NewObject<Datum>,|ctx,frame,program,budget,reference|{let _=ctx;let snapshot=compiled_snapshot(ctx,&program.0,budget)?;let(s,d)=api::reference(&reference.0)?;frame.with(ctx,|f|{let expected=vec![f.dims[1].clone(),f.dims[0].clone(),f.dims[0].clone(),f.dims[0].clone()];if snapshot.input_dimensions.len()!=4{return Err(err("boundary callback arity mismatch"));}for(a,b)in snapshot.input_dimensions.iter().zip(expected){api::same_dimension(a,&b)?;}api::same_dimension(&snapshot.output_dimension,&d)?;let data=Datum::callback_value(snapshot.clone(),f,s)?;Ok(NewObject::dependent(Datum{data,frame:f.clone(),value_scale:s,value_dimension:d.clone(),snapshot:Some(snapshot)},vec![d],vec![program.0.clone()]))})});
object_binding!(grid_datum_eval,"num_grid_datum_eval","crate::sim::boundary::BcValue::eval","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(datum:Handle<Datum>,time:f64,point:Vec<f64>)->Dynamic,|ctx,datum,time,point|{let _=ctx;if point.len()!=3{return Err(err("boundary point must have three coordinates"));}datum.with(ctx,|d|{let v=d.data.eval(time,[point[0],point[1],point[2]]);d.finish()?;Ok(Dynamic(api::quantity(v*d.value_scale,d.value_dimension.clone())?))})});
object_binding!(grid_dirichlet,"num_grid_dirichlet","crate::sim::boundary::BoundaryCondition::dirichlet","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(frame:Handle<Frame>,value:f64)->NewObject<Boundary>,|ctx,frame,value|{let _=ctx;frame.with(ctx,|f|Ok(NewObject::new(Boundary{data:boundary::BoundaryCondition::dirichlet(value),frame:f.clone(),datums:vec![]},f.dims.to_vec())))});
object_binding!(grid_neumann,"num_grid_neumann","crate::sim::boundary::BoundaryCondition::neumann","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(frame:Handle<Frame>,value:f64)->NewObject<Boundary>,|ctx,frame,value|{let _=ctx;frame.with(ctx,|f|Ok(NewObject::new(Boundary{data:boundary::BoundaryCondition::neumann(value),frame:f.clone(),datums:vec![]},f.dims.to_vec())))});
object_binding!(grid_robin,"num_grid_robin","crate::sim::boundary::BoundaryCondition::robin","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(frame:Handle<Frame>,alpha:f64,beta:f64,value:f64)->NewObject<Boundary>,|ctx,frame,alpha,beta,value|{let _=ctx;if alpha==0.0&&beta==0.0{return Err(err("Robin alpha and beta cannot both vanish"));}frame.with(ctx,|f|Ok(NewObject::new(Boundary{data:boundary::BoundaryCondition::robin(alpha,beta,value),frame:f.clone(),datums:vec![]},f.dims.to_vec())))});
object_binding!(grid_boundary_from_datum,"num_grid_boundary_from_datum","","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(datum:Handle<Datum>,kind:String,alpha:f64,beta:f64)->NewObject<Boundary>,|ctx,datum,kind,alpha,beta|{let _=ctx;datum.with(ctx,|d|{let data=match kind.as_str(){"Dirichlet"=>{d.require(false)?;boundary::BoundaryCondition::Dirichlet(d.data.clone())},"Neumann"=>{d.require(true)?;boundary::BoundaryCondition::Neumann(d.data.clone())},"Robin"=>{d.require(false)?;if alpha==0.0&&beta==0.0{return Err(err("degenerate Robin condition"));}boundary::BoundaryCondition::Robin{alpha,beta,g:d.data.clone()}},_=>return Err(err("boundary kind must be Dirichlet, Neumann or Robin"))};Ok(NewObject::dependent(Boundary{data,frame:d.frame.clone(),datums:vec![d.clone()]},d.frame.dims.to_vec(),vec![datum.value.clone()]))})});
object_binding!(grid_boundary_kind,"num_grid_boundary_kind","","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(frame:Handle<Frame>,kind:String)->NewObject<Boundary>,|ctx,frame,kind|{let _=ctx;frame.with(ctx,|f|{let data=match kind.as_str(){"Periodic"=>boundary::BoundaryCondition::Periodic,"Symmetry"=>boundary::BoundaryCondition::Symmetry,"Absorbing"=>boundary::BoundaryCondition::Absorbing,_=>return Err(err("kind must be Periodic, Symmetry or Absorbing"))};Ok(NewObject::new(Boundary{data,frame:f.clone(),datums:vec![]},f.dims.to_vec()))})});
object_binding!(grid_boundary_is_dirichlet,"num_grid_boundary_is_dirichlet","crate::sim::boundary::BoundaryCondition::is_dirichlet","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(boundary:Handle<Boundary>)->bool,|ctx,boundary|{let _=ctx;boundary.with(ctx,|b|Ok(b.data.is_dirichlet()))});
object_binding!(grid_boundary_is_periodic,"num_grid_boundary_is_periodic","crate::sim::boundary::BoundaryCondition::is_periodic","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(boundary:Handle<Boundary>)->bool,|ctx,boundary|{let _=ctx;boundary.with(ctx,|b|Ok(b.data.is_periodic()))});
object_binding!(grid_face_axis,"num_grid_face_axis","crate::sim::boundary::Face::axis","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(index:usize)->usize,|ctx,index|{let _=ctx;Ok(face(index)?.axis())});
object_binding!(grid_face_is_max,"num_grid_face_is_max","crate::sim::boundary::Face::is_max","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(index:usize)->bool,|ctx,index|{let _=ctx;Ok(face(index)?.is_max())});
object_binding!(grid_face_index,"num_grid_face_index","crate::sim::boundary::Face::index","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(index:usize)->usize,|ctx,index|{let _=ctx;Ok(face(index)?.index())});
object_binding!(grid_face_of_axis,"num_grid_face_of_axis","crate::sim::boundary::Face::of_axis","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(axis:usize,max:bool)->usize,|ctx,axis,max|{let _=ctx;if axis>2{return Err(err("axis must be zero through two"));}Ok(boundary::Face::of_axis(axis,max).index())});
object_binding!(grid_faces_uniform,"num_grid_faces_uniform","crate::sim::boundary::FaceBcs::uniform","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(boundary:Handle<Boundary>)->NewObject<Faces>,|ctx,boundary|{let _=ctx;boundary.with(ctx,|b|Ok(NewObject::dependent(Faces{data:boundary::FaceBcs::uniform(&b.data),frame:b.frame.clone(),datums:std::array::from_fn(|_|b.datums.clone())},b.frame.dims.to_vec(),vec![boundary.value.clone()])))});
object_binding!(grid_faces_dirichlet_zero,"num_grid_faces_dirichlet_zero","crate::sim::boundary::FaceBcs::dirichlet_zero","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(frame:Handle<Frame>)->NewObject<Faces>,|ctx,frame|{let _=ctx;frame.with(ctx,|f|Ok(NewObject::new(Faces{data:boundary::FaceBcs::dirichlet_zero(),frame:f.clone(),datums:std::array::from_fn(|_|vec![])},f.dims.to_vec())))});
object_binding!(grid_faces_with,"num_grid_faces_with","crate::sim::boundary::FaceBcs::with","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(bcs:Handle<Faces>,boundary:Handle<Boundary>,index:usize)->NewObject<Faces>,|ctx,bcs,boundary,index|{let _=ctx;bcs.with(ctx,|b|boundary.with(ctx,|c|{b.frame.same(&c.frame)?;let mut datums=b.datums.clone();let selected=face(index)?.index();datums[selected].clone_from(&c.datums);Ok(NewObject::dependent(Faces{data:b.data.clone().with(face(index)?,c.data.clone()),frame:b.frame.clone(),datums},b.frame.dims.to_vec(),vec![bcs.value.clone(),boundary.value.clone()]))}))});
object_binding!(grid_faces_with_axis,"num_grid_faces_with_axis","crate::sim::boundary::FaceBcs::with_axis","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(bcs:Handle<Faces>,boundary:Handle<Boundary>,axis:usize)->NewObject<Faces>,|ctx,bcs,boundary,axis|{let _=ctx;if axis>2{return Err(err("axis must be zero through two"));}bcs.with(ctx,|b|boundary.with(ctx,|c|{b.frame.same(&c.frame)?;let mut datums=b.datums.clone();datums[2*axis].clone_from(&c.datums);datums[2*axis+1].clone_from(&c.datums);Ok(NewObject::dependent(Faces{data:b.data.clone().with_axis(axis,&c.data),frame:b.frame.clone(),datums},b.frame.dims.to_vec(),vec![bcs.value.clone(),boundary.value.clone()]))}))});
object_binding!(grid_faces_homogeneous,"num_grid_faces_homogeneous","crate::sim::boundary::FaceBcs::homogeneous","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(bcs:Handle<Faces>)->NewObject<Faces>,|ctx,bcs|{let _=ctx;bcs.with(ctx,|b|Ok(NewObject::new(Faces{data:b.data.homogeneous(),frame:b.frame.clone(),datums:std::array::from_fn(|_|vec![])},b.frame.dims.to_vec())))});
object_binding!(grid_faces_get,"num_grid_faces_get","crate::sim::boundary::FaceBcs::get","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(bcs:Handle<Faces>,index:usize)->NewObject<Boundary>,|ctx,bcs,index|{let _=ctx;bcs.with(ctx,|b|Ok(NewObject::dependent(Boundary{data:b.data.get(face(index)?).clone(),frame:b.frame.clone(),datums:b.datums[face(index)?.index()].clone()},b.frame.dims.to_vec(),vec![bcs.value.clone()])) )});
object_binding!(grid_faces_axis_periodic,"num_grid_faces_axis_periodic","crate::sim::boundary::FaceBcs::axis_periodic","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(bcs:Handle<Faces>,axis:usize)->bool,|ctx,bcs,axis|{let _=ctx;if axis>2{return Err(err("axis must be zero through two"));}bcs.with(ctx,|b|b.data.axis_periodic(axis).map_err(err))});
object_binding!(grid_structured_new,"num_grid_structured_new","crate::sim::structured::Grid::new","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(bcs:Handle<Faces>,n:Vec<usize>,lo:Vec<f64>,hi:Vec<f64>)->NewObject<Grid>,|ctx,bcs,n,lo,hi|{let _=ctx;grid_size(&n)?;bcs.with(ctx,|b|{let data=structured::Grid::new(&n,&lo,&hi,&b.data).map_err(err)?;Ok(NewObject::dependent(Grid{data,frame:b.frame.clone(),bcs:b.clone()},b.frame.dims.to_vec(),vec![bcs.value.clone()]))})});
object_binding!(grid_structured_unit,"num_grid_structured_unit","crate::sim::structured::Grid::unit","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(bcs:Handle<Faces>,dimension:usize,n:usize)->NewObject<Grid>,|ctx,bcs,dimension,n|{let _=ctx;if dimension==0||dimension>3{return Err(err("grid requires one through three dimensions"));}grid_size(&vec![n;dimension])?;bcs.with(ctx,|b|{let data=structured::Grid::unit(dimension,n,&b.data).map_err(err)?;Ok(NewObject::dependent(Grid{data,frame:b.frame.clone(),bcs:b.clone()},b.frame.dims.to_vec(),vec![bcs.value.clone()]))})});
object_binding!(grid_structured_len,"num_grid_structured_len","crate::sim::structured::Grid::len","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(grid:Handle<Grid>)->usize,|ctx,grid|{let _=ctx;grid.with(ctx,|g|Ok(g.data.len()))});
object_binding!(grid_structured_is_empty,"num_grid_structured_is_empty","crate::sim::structured::Grid::is_empty","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(grid:Handle<Grid>)->bool,|ctx,grid|{let _=ctx;grid.with(ctx,|g|Ok(g.data.is_empty()))});
object_binding!(grid_structured_h,"num_grid_structured_h","crate::sim::structured::Grid::h","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(grid:Handle<Grid>,axis:usize)->Dynamic,|ctx,grid,axis|{let _=ctx;grid.with(ctx,|g|{if axis>=g.data.dim{return Err(err("axis is inactive"));}Ok(Dynamic(api::quantity(g.data.h(axis)*g.frame.length,g.frame.dims[0].clone())?))})});
object_binding!(grid_structured_h_min,"num_grid_structured_h_min","crate::sim::structured::Grid::h_min","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(grid:Handle<Grid>)->Dynamic,|ctx,grid|{let _=ctx;grid.with(ctx,|g|{Ok(Dynamic(api::quantity(g.data.h_min()*g.frame.length,g.frame.dims[0].clone())?))})});
object_binding!(grid_structured_index,"num_grid_structured_index","crate::sim::structured::Grid::index","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(grid:Handle<Grid>,ijk:Vec<usize>)->usize,|ctx,grid,ijk|{let _=ctx;grid.with(ctx,|g|{let i:[usize;3]=ijk.try_into().map_err(|_|err("grid index needs three coordinates"))?;if i.iter().zip(g.data.n).any(|(i,n)|*i>=n){return Err(err("grid index out of range"));}Ok(g.data.index(i))})});
object_binding!(grid_structured_ijk,"num_grid_structured_ijk","crate::sim::structured::Grid::ijk","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(grid:Handle<Grid>,index:usize)->Vec<usize>,|ctx,grid,index|{let _=ctx;grid.with(ctx,|g|{if index>=g.data.len(){return Err(err("grid index out of range"));}Ok(g.data.ijk(index).to_vec())})});
object_binding!(grid_structured_coord,"num_grid_structured_coord","crate::sim::structured::Grid::coord","Native continuum model in an explicit positive Length/Time/field reference frame. Numeric coordinates, boundary coefficients and field arrays are normalized in that frame. Session capabilities, model frames, physical callback domains and work/allocation budgets are checked.",true,(grid:Handle<Grid>,index:usize)->Dynamic,|ctx,grid,index|{let _=ctx;grid.with(ctx,|g|{if index>=g.data.len(){return Err(err("grid index out of range"));}let data=g.data.coord(index).into_iter().map(|x|super::super::output_scalar(x*g.frame.length)).collect::<Result<Vec<_>,_>>()?;Ok(Dynamic(Value::Tensor(TensorValue{shape:vec![3],elem:TypeExpr::Float.with_dim(g.frame.dims[0].clone()),data})))})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(i, BINDINGS, "numerical_vector")?;
    structured_ops::install(i)?;
    transport::install(i)?;
    boundary_models::install(i)?;
    spectral::install(i)?;
    mechanics::install(i)?;
    fdm::install(i)?;
    fem::install(i)?;
    multigrid::install(i)?;
    navier::install(i)?;
    transforms::install(i)?;
    dg::install(i)?;
    triangles::install(i)?;
    euler::install(i)?;
    cutcell::install(i)?;
    flow_models::install(i)?;
    pde_facade::install(i)?;
    bem::install(i)?;
    gas::install(i)?;
    acoustics::install(i)?;
    elastic::install(i)?;
    fem_legacy::install(i)?;
    elastic_quad::install(i)?;
    legacy_physics::install(i)
}
static BINDINGS: &[&api::Binding] = &[
    &grid_frame_new::BINDING,
    &grid_frame_references::BINDING,
    &grid_datum_constant::BINDING,
    &grid_datum_time::BINDING,
    &grid_datum_space::BINDING,
    &grid_datum_spacetime::BINDING,
    &grid_datum_eval::BINDING,
    &grid_dirichlet::BINDING,
    &grid_neumann::BINDING,
    &grid_robin::BINDING,
    &grid_boundary_from_datum::BINDING,
    &grid_boundary_kind::BINDING,
    &grid_boundary_is_dirichlet::BINDING,
    &grid_boundary_is_periodic::BINDING,
    &grid_face_axis::BINDING,
    &grid_face_is_max::BINDING,
    &grid_face_index::BINDING,
    &grid_face_of_axis::BINDING,
    &grid_faces_uniform::BINDING,
    &grid_faces_dirichlet_zero::BINDING,
    &grid_faces_with::BINDING,
    &grid_faces_with_axis::BINDING,
    &grid_faces_homogeneous::BINDING,
    &grid_faces_get::BINDING,
    &grid_faces_axis_periodic::BINDING,
    &grid_structured_new::BINDING,
    &grid_structured_unit::BINDING,
    &grid_structured_len::BINDING,
    &grid_structured_is_empty::BINDING,
    &grid_structured_h::BINDING,
    &grid_structured_h_min::BINDING,
    &grid_structured_index::BINDING,
    &grid_structured_ijk::BINDING,
    &grid_structured_coord::BINDING,
];

#[path = "structured.rs"]
mod structured_ops;

struct NormalizedFloat(f64);
impl NativeArg for NormalizedFloat {
    fn ty() -> TypeExpr {
        TypeExpr::Float
    }

    fn decode(
        v: &Value,
        i: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        api::normalized(v)?;
        Ok(Self(super::numeric(v, i)?))
    }
}

mod transport;

mod boundary_models;

mod spectral;

mod mechanics;

mod fdm;

mod fem;

mod multigrid;

mod navier;

mod transforms;

mod legacy_physics;

mod dg;

mod bem;

mod triangles;

mod gas;

mod acoustics;

mod euler;

mod cutcell;

mod elastic;

mod flow_models;

mod pde_facade;

mod fem_legacy;

mod elastic_quad;
