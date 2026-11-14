//! Native marker, fluid and elasticity boundary value models.
use super::super::native_object_api::Maybe;
use super::Boundary;
use super::Datum;
use super::Frame;
use super::Handle;
use super::Installer;
use super::NativeType;
use super::NewObject;
use super::OpError;
use super::RuleError;
use super::api;
use super::boundary;
use super::err;
use super::face;
use super::object_binding;
#[derive(Clone)]
pub(super) struct Markers {
    pub(super) data: boundary::MarkerBcs<boundary::BoundaryCondition>,
    pub(super) frame: Frame,
    pub(super) callbacks: std::collections::BTreeMap<i32, Vec<Datum>>,
    pub(super) default: Vec<Datum>,
}
impl NativeType for Markers {
    const KIND: &'static str = "rssn/native/marker-bcs/v1";
    const MAX_BYTES: usize = 262144;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        if self.callbacks.len() > 256 {
            return Err(err("boundary marker budget exceeded"));
        }
        Ok(())
    }
}
#[derive(Clone)]
pub(super) struct Fluid {
    pub data: boundary::FluidBoundary,
    pub frame: Frame,
}
impl NativeType for Fluid {
    const KIND: &'static str = "rssn/native/fluid-boundary/v1";
    const MAX_BYTES: usize = 4096;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()
    }
}
#[derive(Clone)]
pub(super) struct FluidFaces {
    pub data: boundary::FluidFaceBcs,
    pub frame: Frame,
}
impl NativeType for FluidFaces {
    const KIND: &'static str = "rssn/native/fluid-face-bcs/v1";
    const MAX_BYTES: usize = 32768;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()
    }
}
impl NativeType for boundary::ElasticBoundary {
    const KIND: &'static str = "rssn/native/elastic-boundary/v1";
    const MAX_BYTES: usize = 4096;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }
}
object_binding!(boundary_marker_new,"num_grid_marker_new","crate::sim::boundary::MarkerBcs<T>::new","Native boundary model. Scalar/fluid quantities use the explicit grid reference frame; standalone elasticity boundary coefficients are normalized and require the solver material/reference frame. Marker and capability bounds are checked.",true,(frame:Handle<Frame>)->NewObject<Markers>,|ctx,frame|{let _=ctx;frame.with(ctx,|f|Ok(NewObject::new(Markers{data:boundary::MarkerBcs::new(),frame:f.clone(),callbacks:std::collections::BTreeMap::new(),default:vec![]},f.dims.to_vec())))});
object_binding!(boundary_marker_with_default,"num_grid_marker_with_default","crate::sim::boundary::MarkerBcs<T>::with_default","Native boundary model. Scalar/fluid quantities use the explicit grid reference frame; standalone elasticity boundary coefficients are normalized and require the solver material/reference frame. Marker and capability bounds are checked.",true,(boundary:Handle<Boundary>)->NewObject<Markers>,|ctx,boundary|{let _=ctx;boundary.with(ctx,|b|Ok(NewObject::dependent(Markers{data:boundary::MarkerBcs::with_default(b.data.clone()),frame:b.frame.clone(),callbacks:std::collections::BTreeMap::new(),default:b.datums.clone()},b.frame.dims.to_vec(),vec![boundary.value.clone()])))});
object_binding!(boundary_marker_set,"num_grid_marker_set","crate::sim::boundary::MarkerBcs<T>::set","Native boundary model. Scalar/fluid quantities use the explicit grid reference frame; standalone elasticity boundary coefficients are normalized and require the solver material/reference frame. Marker and capability bounds are checked.",true,(markers:Handle<Markers>,marker:i32,boundary:Handle<Boundary>)->NewObject<Markers>,|ctx,markers,marker,boundary|{let _=ctx;markers.with(ctx,|m|boundary.with(ctx,|b|{m.frame.same(&b.frame)?;if m.callbacks.len()>=256&&!m.callbacks.contains_key(&marker){return Err(err("boundary marker allocation budget exceeded"));}let mut callbacks=m.callbacks.clone();callbacks.insert(marker,b.datums.clone());Ok(NewObject::dependent(Markers{data:m.data.clone().set(marker,b.data.clone()),frame:m.frame.clone(),callbacks,default:m.default.clone()},m.frame.dims.to_vec(),vec![markers.value.clone(),boundary.value.clone()]))}))});
object_binding!(boundary_marker_get,"num_grid_marker_get","crate::sim::boundary::MarkerBcs<T>::get","Native boundary model. Scalar/fluid quantities use the explicit grid reference frame; standalone elasticity boundary coefficients are normalized and require the solver material/reference frame. Marker and capability bounds are checked.",true,(markers:Handle<Markers>,marker:i32)->Maybe<NewObject<Boundary>>,|ctx,markers,marker|{let _=ctx;markers.with(ctx,|m|Ok(Maybe(m.data.get(marker).map(|data|NewObject::dependent(Boundary{data:data.clone(),frame:m.frame.clone(),datums:m.callbacks.get(&marker).unwrap_or(&m.default).clone()},m.frame.dims.to_vec(),vec![markers.value.clone()])))))});
object_binding!(boundary_fluid_wall,"num_grid_fluid_wall","crate::sim::boundary::FluidBoundary::wall","Native boundary model. Scalar/fluid quantities use the explicit grid reference frame; standalone elasticity boundary coefficients are normalized and require the solver material/reference frame. Marker and capability bounds are checked.",true,(frame:Handle<Frame>)->NewObject<Fluid>,|ctx,frame|{let _=ctx;frame.with(ctx,|f|{Ok(NewObject::new(Fluid{data:boundary::FluidBoundary::wall(),frame:f.clone()},f.dims.to_vec()))})});
object_binding!(boundary_fluid_outflow,"num_grid_fluid_outflow","crate::sim::boundary::FluidBoundary::outflow","Native boundary model. Scalar/fluid quantities use the explicit grid reference frame; standalone elasticity boundary coefficients are normalized and require the solver material/reference frame. Marker and capability bounds are checked.",true,(frame:Handle<Frame>)->NewObject<Fluid>,|ctx,frame|{let _=ctx;frame.with(ctx,|f|{Ok(NewObject::new(Fluid{data:boundary::FluidBoundary::outflow(),frame:f.clone()},f.dims.to_vec()))})});
object_binding!(boundary_fluid_moving_wall,"num_grid_fluid_moving_wall","crate::sim::boundary::FluidBoundary::moving_wall","Native boundary model. Scalar/fluid quantities use the explicit grid reference frame; standalone elasticity boundary coefficients are normalized and require the solver material/reference frame. Marker and capability bounds are checked.",true,(frame:Handle<Frame>,velocity:Vec<f64>)->NewObject<Fluid>,|ctx,frame,velocity|{let _=ctx;frame.with(ctx,|f|{let v:[f64;3]=velocity.try_into().map_err(|_|err("fluid velocity needs three components"))?;Ok(NewObject::new(Fluid{data:boundary::FluidBoundary::moving_wall(v),frame:f.clone()},f.dims.to_vec()))})});
object_binding!(boundary_fluid_inflow,"num_grid_fluid_inflow","crate::sim::boundary::FluidBoundary::inflow","Native boundary model. Scalar/fluid quantities use the explicit grid reference frame; standalone elasticity boundary coefficients are normalized and require the solver material/reference frame. Marker and capability bounds are checked.",true,(frame:Handle<Frame>,velocity:Vec<f64>)->NewObject<Fluid>,|ctx,frame,velocity|{let _=ctx;frame.with(ctx,|f|{let v:[f64;3]=velocity.try_into().map_err(|_|err("fluid velocity needs three components"))?;Ok(NewObject::new(Fluid{data:boundary::FluidBoundary::inflow(v),frame:f.clone()},f.dims.to_vec()))})});
object_binding!(boundary_fluid_kind,"num_grid_fluid_kind","","Native boundary model. Scalar/fluid quantities use the explicit grid reference frame; standalone elasticity boundary coefficients are normalized and require the solver material/reference frame. Marker and capability bounds are checked.",true,(frame:Handle<Frame>,kind:String)->NewObject<Fluid>,|ctx,frame,kind|{let _=ctx;frame.with(ctx,|f|{let data=match kind.as_str(){"Slip"=>boundary::FluidBoundary::Slip,"Periodic"=>boundary::FluidBoundary::Periodic,_=>return Err(err("fluid kind must be Slip or Periodic"))};Ok(NewObject::new(Fluid{data,frame:f.clone()},f.dims.to_vec()))})});
object_binding!(boundary_fluid_faces_uniform,"num_grid_fluid_faces_uniform","crate::sim::boundary::FluidFaceBcs::uniform","Native boundary model. Scalar/fluid quantities use the explicit grid reference frame; standalone elasticity boundary coefficients are normalized and require the solver material/reference frame. Marker and capability bounds are checked.",true,(boundary:Handle<Fluid>)->NewObject<FluidFaces>,|ctx,boundary|{let _=ctx;boundary.with(ctx,|b|Ok(NewObject::new(FluidFaces{data:boundary::FluidFaceBcs::uniform(&b.data),frame:b.frame.clone()},b.frame.dims.to_vec())))});
object_binding!(boundary_fluid_faces_with,"num_grid_fluid_faces_with","crate::sim::boundary::FluidFaceBcs::with","Native boundary model. Scalar/fluid quantities use the explicit grid reference frame; standalone elasticity boundary coefficients are normalized and require the solver material/reference frame. Marker and capability bounds are checked.",true,(bcs:Handle<FluidFaces>,boundary:Handle<Fluid>,index:usize)->NewObject<FluidFaces>,|ctx,bcs,boundary,index|{let _=ctx;bcs.with(ctx,|b|boundary.with(ctx,|c|{b.frame.same(&c.frame)?;Ok(NewObject::new(FluidFaces{data:b.data.clone().with(face(index)?,c.data.clone()),frame:b.frame.clone()},b.frame.dims.to_vec()))}))});
object_binding!(boundary_fluid_faces_with_axis,"num_grid_fluid_faces_with_axis","crate::sim::boundary::FluidFaceBcs::with_axis","Native boundary model. Scalar/fluid quantities use the explicit grid reference frame; standalone elasticity boundary coefficients are normalized and require the solver material/reference frame. Marker and capability bounds are checked.",true,(bcs:Handle<FluidFaces>,boundary:Handle<Fluid>,axis:usize)->NewObject<FluidFaces>,|ctx,bcs,boundary,axis|{let _=ctx;if axis>2{return Err(err("axis must be zero through two"));}bcs.with(ctx,|b|boundary.with(ctx,|c|{b.frame.same(&c.frame)?;Ok(NewObject::new(FluidFaces{data:b.data.clone().with_axis(axis,&c.data),frame:b.frame.clone()},b.frame.dims.to_vec()))}))});
object_binding!(boundary_fluid_faces_get,"num_grid_fluid_faces_get","crate::sim::boundary::FluidFaceBcs::get","Native boundary model. Scalar/fluid quantities use the explicit grid reference frame; standalone elasticity boundary coefficients are normalized and require the solver material/reference frame. Marker and capability bounds are checked.",true,(bcs:Handle<FluidFaces>,index:usize)->NewObject<Fluid>,|ctx,bcs,index|{let _=ctx;bcs.with(ctx,|b|Ok(NewObject::new(Fluid{data:b.data.get(face(index)?).clone(),frame:b.frame.clone()},b.frame.dims.to_vec())))});
object_binding!(boundary_elastic_clamped,"num_grid_elastic_clamped","crate::sim::boundary::ElasticBoundary::clamped","Native boundary model. Scalar/fluid quantities use the explicit grid reference frame; standalone elasticity boundary coefficients are normalized and require the solver material/reference frame. Marker and capability bounds are checked.",true,()->NewObject<boundary::ElasticBoundary>,|ctx,|{let _=ctx;Ok(NewObject::normalized(boundary::ElasticBoundary::clamped()))});
object_binding!(boundary_elastic_traction,"num_grid_elastic_traction","crate::sim::boundary::ElasticBoundary::traction","Native boundary model. Scalar/fluid quantities use the explicit grid reference frame; standalone elasticity boundary coefficients are normalized and require the solver material/reference frame. Marker and capability bounds are checked.",true,(traction:Vec<f64>)->NewObject<boundary::ElasticBoundary>,|ctx,traction|{let _=ctx;let t:[f64;3]=traction.try_into().map_err(|_|err("traction needs three components"))?;Ok(NewObject::normalized(boundary::ElasticBoundary::traction(t)))});
object_binding!(boundary_elastic_roller,"num_grid_elastic_roller","crate::sim::boundary::ElasticBoundary::roller","Native boundary model. Scalar/fluid quantities use the explicit grid reference frame; standalone elasticity boundary coefficients are normalized and require the solver material/reference frame. Marker and capability bounds are checked.",true,(component:usize)->NewObject<boundary::ElasticBoundary>,|ctx,component|{let _=ctx;if component>2{return Err(err("elastic displacement component must be zero through two"));}Ok(NewObject::normalized(boundary::ElasticBoundary::roller(component)))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(i, BINDINGS, "numerical_vector")
}
static BINDINGS: &[&api::Binding] = &[
    &boundary_marker_new::BINDING,
    &boundary_marker_with_default::BINDING,
    &boundary_marker_set::BINDING,
    &boundary_marker_get::BINDING,
    &boundary_fluid_wall::BINDING,
    &boundary_fluid_outflow::BINDING,
    &boundary_fluid_moving_wall::BINDING,
    &boundary_fluid_inflow::BINDING,
    &boundary_fluid_kind::BINDING,
    &boundary_fluid_faces_uniform::BINDING,
    &boundary_fluid_faces_with::BINDING,
    &boundary_fluid_faces_with_axis::BINDING,
    &boundary_fluid_faces_get::BINDING,
    &boundary_elastic_clamped::BINDING,
    &boundary_elastic_traction::BINDING,
    &boundary_elastic_roller::BINDING,
];
