//! Native finite-element parts in an explicit Length/Time/Mass/Temperature reference frame.
use super::DimExpr;
use super::Dynamic;
use super::Exp;
use super::Handle;
use super::Input;
use super::Installer;
use super::NativeType;
use super::NewObject;
use super::OpError;
use super::Output;
use super::Rational;
use super::RuleError;
use super::Value;
use super::api;
use super::canonical;
use super::err;
use super::object_binding;
use crate::kernels::matrix::Matrix;
use crate::sim::physics_fea as fea;

pub(super) const REFERENCE_DOMAINS: [&str; 4] = ["Length", "Time", "Mass", "Temperature"];

#[derive(Clone)]
pub(super) struct References {
    pub(super) scales: [f64; 4],
    pub(super) dims: [DimExpr; 4],
}
impl NativeType for References {
    const KIND: &'static str = "rssn/native/mechanical-references/v1";
    const MAX_BYTES: usize = 4096;

    fn validate(&self) -> Result<(), OpError> {
        for (i, name) in REFERENCE_DOMAINS.into_iter().enumerate() {
            api::same_dimension(&self.dims[i], &DimExpr::base(name))?;
            if !self.scales[i].is_finite() || self.scales[i] <= 0.0 {
                return Err(err("mechanical references must be finite positive"));
            }
        }
        Ok(())
    }
}
impl References {
    pub(super) fn derived(
        &self,
        powers: [i64; 4],
    ) -> Result<(f64, DimExpr), OpError> {
        let mut numerator = Vec::new();
        let mut denominator = Vec::new();
        let mut dimensions = Vec::new();
        for (i, p) in powers.into_iter().enumerate() {
            if p.unsigned_abs() > 16 {
                return Err(err("mechanical reference exponent budget exceeded"));
            }
            let factors = if p >= 0 {
                &mut numerator
            } else {
                &mut denominator
            };
            factors.extend(std::iter::repeat_n(
                self.scales[i],
                p.unsigned_abs() as usize,
            ));
            if p != 0 {
                dimensions.push(DimExpr::Pow(
                    Box::new(self.dims[i].clone()),
                    Exp::Rat(Rational::new(p, 1)),
                ));
            }
        }
        let value = api::scale_product(1.0, &numerator, &denominator)?;
        if value <= 0.0 {
            return Err(err("mechanical derived reference overflow/underflow"));
        }
        Ok((value, canonical(DimExpr::Mul(dimensions))?))
    }

    fn same(
        &self,
        b: &Self,
    ) -> Result<(), OpError> {
        for (a, b) in self.scales.iter().zip(b.scales) {
            if a.to_bits() != b.to_bits() {
                return Err(err("mechanical reference frames differ"));
            }
        }
        Ok(())
    }
}
pub(super) trait Part: 'static {
    const KIND: &'static str;
    fn valid(&self) -> Result<(), OpError>;
}
#[derive(Clone)]
pub(super) struct Model<T> {
    pub(super) data: T,
    pub(super) frame: References,
}
impl<T: Part> NativeType for Model<T> {
    const KIND: &'static str = T::KIND;
    const MAX_BYTES: usize = 16384;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        self.data.valid()
    }
}
fn positive(values: &[f64]) -> Result<(), OpError> {
    if values.iter().any(|v| !v.is_finite() || *v <= 0.0) {
        return Err(err(
            "finite positive element geometry/material coefficients required",
        ));
    }
    Ok(())
}
impl Part for fea::Material {
    const KIND: &'static str = "rssn/native/fea-material/v1";

    fn valid(&self) -> Result<(), OpError> {
        positive(&[self.youngs_modulus])?;
        if [
            self.youngs_modulus,
            self.poissons_ratio,
            self.density,
            self.thermal_conductivity,
            self.thermal_expansion,
            self.yield_strength,
        ]
        .iter()
        .any(|v| !v.is_finite())
            || self.poissons_ratio <= -1.0
            || self.poissons_ratio >= 0.5
            || self.density < 0.0
            || self.thermal_conductivity < 0.0
            || self.yield_strength < 0.0
        {
            return Err(err("invalid isotropic material coefficients"));
        }
        Ok(())
    }
}
impl Part for fea::Node2D {
    const KIND: &'static str = "rssn/native/fea-node2/v1";

    fn valid(&self) -> Result<(), OpError> {
        if self.id > 1_000_000 || !self.x.is_finite() || !self.y.is_finite() {
            return Err(err("invalid finite node or identifier budget"));
        }
        Ok(())
    }
}
impl Part for fea::Node3D {
    const KIND: &'static str = "rssn/native/fea-node3/v1";

    fn valid(&self) -> Result<(), OpError> {
        if self.id > 1_000_000 || [self.x, self.y, self.z].iter().any(|v| !v.is_finite()) {
            return Err(err("invalid finite node or identifier budget"));
        }
        Ok(())
    }
}
impl Part for fea::LinearElement1D {
    const KIND: &'static str = "rssn/native/fea-line/v1";

    fn valid(&self) -> Result<(), OpError> {
        positive(&[self.length, self.youngs_modulus, self.area])
    }
}
impl Part for fea::BeamElement2D {
    const KIND: &'static str = "rssn/native/fea-beam/v1";

    fn valid(&self) -> Result<(), OpError> {
        positive(&[
            self.length,
            self.youngs_modulus,
            self.area,
            self.moment_of_inertia,
        ])
    }
}
impl Part for fea::ThermalElement1D {
    const KIND: &'static str = "rssn/native/fea-thermal-line/v1";

    fn valid(&self) -> Result<(), OpError> {
        positive(&[self.length, self.conductivity, self.area])
    }
}
impl Part for fea::ThermalTriangle2D {
    const KIND: &'static str = "rssn/native/fea-thermal-triangle/v1";

    fn valid(&self) -> Result<(), OpError> {
        positive(&[self.area(), self.thickness, self.conductivity])
    }
}
impl Part for fea::TriangleElement2D {
    const KIND: &'static str = "rssn/native/fea-triangle/v1";

    fn valid(&self) -> Result<(), OpError> {
        if self.nodes.iter().any(|n| *n > 1_000_000) {
            return Err(err("element node identifier budget exceeded"));
        }
        positive(&[self.area(), self.thickness])?;
        self.material.valid()
    }
}
impl Part for fea::TetrahedronElement3D {
    const KIND: &'static str = "rssn/native/fea-tetrahedron/v1";

    fn valid(&self) -> Result<(), OpError> {
        if self.nodes.iter().any(|n| *n > 1_000_000) {
            return Err(err("element node identifier budget exceeded"));
        }
        positive(&[self.volume()])?;
        self.material.valid()
    }
}
fn matrix(
    m: Matrix<f64>,
    frame: &References,
    powers: [i64; 4],
) -> Result<Dynamic, OpError> {
    let (s, d) = frame.derived(powers)?;
    Ok(Dynamic(api::record(vec![
        ("rows", m.rows().encode()?),
        ("cols", m.cols().encode()?),
        ("data", m.into_data().encode()?),
        ("entry_reference", api::quantity(s, d)?),
    ])))
}
fn scalar(
    value: f64,
    frame: &References,
    powers: [i64; 4],
) -> Result<Dynamic, OpError> {
    let (s, d) = frame.derived(powers)?;
    Ok(Dynamic(api::quantity(value * s, d)?))
}
fn mixed_beam(
    m: Matrix<f64>,
    frame: &References,
    mass: bool,
) -> Result<Dynamic, OpError> {
    let (s, d) = frame.derived(if mass {
        [2, 0, 1, 0]
    } else {
        [2, -2, 1, 0]
    })?;
    let mut references = Vec::new();
    for i in 0..6 {
        references.push((
            format!("dof{i}"),
            if i % 3 == 2 {
                Value::Float(1.0)
            } else {
                api::quantity(frame.scales[0], frame.dims[0].clone())?
            },
        ));
    }
    Ok(Dynamic(api::record(vec![
        ("rows", m.rows().encode()?),
        ("cols", m.cols().encode()?),
        ("data", m.into_data().encode()?),
        ("quadratic_form_reference", api::quantity(s, d)?),
        (
            "dof_references",
            Value::Record(crate::kernel::value::RecordValue { fields: references }),
        ),
    ])))
}
fn material_from_si(
    mut m: fea::Material,
    f: &References,
) -> Result<fea::Material, OpError> {
    let stress = f.derived([-1, -2, 1, 0])?.0;
    m.youngs_modulus /= stress;
    m.yield_strength /= stress;
    m.density /= f.derived([-3, 0, 1, 0])?.0;
    m.thermal_conductivity /= f.derived([1, -3, 1, -1])?.0;
    m.thermal_expansion /= f.derived([0, 0, 0, -1])?.0;
    Ok(m)
}
object_binding!(mechanical_mechanical_references,"num_grid_mechanical_references","","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",false,(length:Input,time:Input,mass:Input,temperature:Input)->NewObject<References>,|ctx,length,time,mass,temperature|{let _=ctx;let(l,dl)=api::reference(&length.0)?;let(t,dt)=api::reference(&time.0)?;let(m,dm)=api::reference(&mass.0)?;let(k,dk)=api::reference(&temperature.0)?;let f=References{scales:[l,t,m,k],dims:[dl,dt,dm,dk]};f.validate()?;Ok(NewObject::new(f.clone(),f.dims.to_vec()))});
object_binding!(mechanical_material_new,"num_grid_material_new","crate::sim::physics_fea::Material::new","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(frame:Handle<References>,young:f64,poisson:f64,density:f64,conductivity:f64,expansion:f64,yield_strength:f64)->NewObject<Model<fea::Material>>,|ctx,frame,young,poisson,density,conductivity,expansion,yield_strength|{let _=ctx;frame.with(ctx,|f|Ok(NewObject::new(Model{data:fea::Material::new(young,poisson,density,conductivity,expansion,yield_strength),frame:f.clone()},f.dims.to_vec())))});
object_binding!(mechanical_material_steel,"num_grid_material_steel","crate::sim::physics_fea::Material::steel","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(frame:Handle<References>)->NewObject<Model<fea::Material>>,|ctx,frame|{let _=ctx;frame.with(ctx,|f|Ok(NewObject::new(Model{data:material_from_si(fea::Material::steel(),f)?,frame:f.clone()},f.dims.to_vec())))});
object_binding!(mechanical_material_aluminum,"num_grid_material_aluminum","crate::sim::physics_fea::Material::aluminum","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(frame:Handle<References>)->NewObject<Model<fea::Material>>,|ctx,frame|{let _=ctx;frame.with(ctx,|f|Ok(NewObject::new(Model{data:material_from_si(fea::Material::aluminum(),f)?,frame:f.clone()},f.dims.to_vec())))});
object_binding!(mechanical_material_copper,"num_grid_material_copper","crate::sim::physics_fea::Material::copper","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(frame:Handle<References>)->NewObject<Model<fea::Material>>,|ctx,frame|{let _=ctx;frame.with(ctx,|f|Ok(NewObject::new(Model{data:material_from_si(fea::Material::copper(),f)?,frame:f.clone()},f.dims.to_vec())))});
object_binding!(mechanical_material_shear_modulus,"num_grid_material_shear_modulus","crate::sim::physics_fea::Material::shear_modulus","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(material:Handle<Model<fea::Material>>)->Dynamic,|ctx,material|{let _=ctx;material.with(ctx,|m|scalar(m.data.shear_modulus(),&m.frame,[-1,-2,1,0]))});
object_binding!(mechanical_material_bulk_modulus,"num_grid_material_bulk_modulus","crate::sim::physics_fea::Material::bulk_modulus","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(material:Handle<Model<fea::Material>>)->Dynamic,|ctx,material|{let _=ctx;material.with(ctx,|m|scalar(m.data.bulk_modulus(),&m.frame,[-1,-2,1,0]))});
object_binding!(mechanical_node2_new,"num_grid_node2_new","crate::sim::physics_fea::Node2D::new","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(frame:Handle<References>,id:usize,x:f64,y:f64)->NewObject<Model<fea::Node2D>>,|ctx,frame,id,x,y|{let _=ctx;frame.with(ctx,|f|Ok(NewObject::new(Model{data:fea::Node2D::new(id,x,y),frame:f.clone()},f.dims.to_vec())))});
object_binding!(mechanical_node3_new,"num_grid_node3_new","crate::sim::physics_fea::Node3D::new","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(frame:Handle<References>,id:usize,x:f64,y:f64,z:f64)->NewObject<Model<fea::Node3D>>,|ctx,frame,id,x,y,z|{let _=ctx;frame.with(ctx,|f|Ok(NewObject::new(Model{data:fea::Node3D::new(id,x,y,z),frame:f.clone()},f.dims.to_vec())))});
object_binding!(mechanical_node2_distance,"num_grid_node2_distance","crate::sim::physics_fea::Node2D::distance_to","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(a:Handle<Model<fea::Node2D>>,b:Handle<Model<fea::Node2D>>)->Dynamic,|ctx,a,b|{let _=ctx;a.with(ctx,|a|b.with(ctx,|b|{a.frame.same(&b.frame)?;scalar(a.data.distance_to(&b.data),&a.frame,[1,0,0,0])}))});
object_binding!(mechanical_line_element_new,"num_grid_line_element_new","","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(frame:Handle<References>,length:f64,young:f64,area:f64)->NewObject<Model<fea::LinearElement1D>>,|ctx,frame,length,young,area|{let _=ctx;frame.with(ctx,|f|Ok(NewObject::new(Model{data:fea::LinearElement1D{length,youngs_modulus:young,area},frame:f.clone()},f.dims.to_vec())))});
object_binding!(mechanical_beam_new,"num_grid_beam_new","crate::sim::physics_fea::BeamElement2D::new","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(frame:Handle<References>,length:f64,young:f64,area:f64,inertia:f64,angle:f64)->NewObject<Model<fea::BeamElement2D>>,|ctx,frame,length,young,area,inertia,angle|{let _=ctx;frame.with(ctx,|f|Ok(NewObject::new(Model{data:fea::BeamElement2D::new(length,young,area,inertia,angle),frame:f.clone()},f.dims.to_vec())))});
object_binding!(mechanical_thermal_line_new,"num_grid_thermal_line_new","crate::sim::physics_fea::ThermalElement1D::new","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(frame:Handle<References>,length:f64,conductivity:f64,area:f64)->NewObject<Model<fea::ThermalElement1D>>,|ctx,frame,length,conductivity,area|{let _=ctx;frame.with(ctx,|f|Ok(NewObject::new(Model{data:fea::ThermalElement1D::new(length,conductivity,area),frame:f.clone()},f.dims.to_vec())))});
object_binding!(mechanical_thermal_triangle_new,"num_grid_thermal_triangle_new","crate::sim::physics_fea::ThermalTriangle2D::new","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(frame:Handle<References>,coords:Vec<[f64;2]>,thickness:f64,conductivity:f64)->NewObject<Model<fea::ThermalTriangle2D>>,|ctx,frame,coords,thickness,conductivity|{let _=ctx;let c:[[f64;2];3]=coords.try_into().map_err(|_|err("thermal triangle needs three coordinates"))?;frame.with(ctx,|f|Ok(NewObject::new(Model{data:fea::ThermalTriangle2D::new(c.map(|p|(p[0],p[1])),thickness,conductivity),frame:f.clone()},f.dims.to_vec())))});
object_binding!(mechanical_triangle_new,"num_grid_triangle_new","crate::sim::physics_fea::TriangleElement2D::new","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(material:Handle<Model<fea::Material>>,nodes:Vec<usize>,coords:Vec<[f64;2]>,thickness:f64,plane_stress:bool)->NewObject<Model<fea::TriangleElement2D>>,|ctx,material,nodes,coords,thickness,plane_stress|{let _=ctx;let n:[usize;3]=nodes.try_into().map_err(|_|err("triangle needs three node IDs"))?;let c:[[f64;2];3]=coords.try_into().map_err(|_|err("triangle needs three coordinates"))?;material.with(ctx,|m|Ok(NewObject::new(Model{data:fea::TriangleElement2D::new(n,c.map(|p|(p[0],p[1])),thickness,m.data,plane_stress),frame:m.frame.clone()},m.frame.dims.to_vec())))});
object_binding!(mechanical_tetrahedron_new,"num_grid_tetrahedron_new","crate::sim::physics_fea::TetrahedronElement3D::new","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(material:Handle<Model<fea::Material>>,nodes:Vec<usize>,coords:Vec<[f64;3]>)->NewObject<Model<fea::TetrahedronElement3D>>,|ctx,material,nodes,coords|{let _=ctx;let n:[usize;4]=nodes.try_into().map_err(|_|err("tetrahedron needs four node IDs"))?;let c:[[f64;3];4]=coords.try_into().map_err(|_|err("tetrahedron needs four coordinates"))?;material.with(ctx,|m|Ok(NewObject::new(Model{data:fea::TetrahedronElement3D::new(n,c.map(|p|(p[0],p[1],p[2])),m.data),frame:m.frame.clone()},m.frame.dims.to_vec())))});
object_binding!(mechanical_line_local_stiffness_matrix,"num_grid_line_local_stiffness_matrix","crate::sim::physics_fea::LinearElement1D::local_stiffness_matrix","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(element:Handle<Model<fea::LinearElement1D>>)->Dynamic,|ctx,element|{let _=ctx;element.with(ctx,|m|matrix(m.data.local_stiffness_matrix(),&m.frame,[0, -2, 1, 0]))});
object_binding!(mechanical_triangle_local_stiffness_matrix,"num_grid_triangle_local_stiffness_matrix","crate::sim::physics_fea::TriangleElement2D::local_stiffness_matrix","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(element:Handle<Model<fea::TriangleElement2D>>)->Dynamic,|ctx,element|{let _=ctx;element.with(ctx,|m|matrix(m.data.local_stiffness_matrix(),&m.frame,[0, -2, 1, 0]))});
object_binding!(mechanical_triangle_b_matrix,"num_grid_triangle_b_matrix","crate::sim::physics_fea::TriangleElement2D::b_matrix","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(element:Handle<Model<fea::TriangleElement2D>>)->Dynamic,|ctx,element|{let _=ctx;element.with(ctx,|m|matrix(m.data.b_matrix(),&m.frame,[-1, 0, 0, 0]))});
object_binding!(mechanical_triangle_constitutive_matrix,"num_grid_triangle_constitutive_matrix","crate::sim::physics_fea::TriangleElement2D::constitutive_matrix","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(element:Handle<Model<fea::TriangleElement2D>>)->Dynamic,|ctx,element|{let _=ctx;element.with(ctx,|m|matrix(m.data.constitutive_matrix(),&m.frame,[-1, -2, 1, 0]))});
object_binding!(mechanical_tetrahedron_local_stiffness_matrix,"num_grid_tetrahedron_local_stiffness_matrix","crate::sim::physics_fea::TetrahedronElement3D::local_stiffness_matrix","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(element:Handle<Model<fea::TetrahedronElement3D>>)->Dynamic,|ctx,element|{let _=ctx;element.with(ctx,|m|matrix(m.data.local_stiffness_matrix(),&m.frame,[0, -2, 1, 0]))});
object_binding!(mechanical_tetrahedron_b_matrix,"num_grid_tetrahedron_b_matrix","crate::sim::physics_fea::TetrahedronElement3D::b_matrix","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(element:Handle<Model<fea::TetrahedronElement3D>>)->Dynamic,|ctx,element|{let _=ctx;element.with(ctx,|m|matrix(m.data.b_matrix(),&m.frame,[-1, 0, 0, 0]))});
object_binding!(mechanical_thermal_line_conductivity_matrix,"num_grid_thermal_line_conductivity_matrix","crate::sim::physics_fea::ThermalElement1D::conductivity_matrix","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(element:Handle<Model<fea::ThermalElement1D>>)->Dynamic,|ctx,element|{let _=ctx;element.with(ctx,|m|matrix(m.data.conductivity_matrix(),&m.frame,[2, -3, 1, -1]))});
object_binding!(mechanical_thermal_triangle_conductivity_matrix,"num_grid_thermal_triangle_conductivity_matrix","crate::sim::physics_fea::ThermalTriangle2D::conductivity_matrix","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(element:Handle<Model<fea::ThermalTriangle2D>>)->Dynamic,|ctx,element|{let _=ctx;element.with(ctx,|m|matrix(m.data.conductivity_matrix(),&m.frame,[2, -3, 1, -1]))});
object_binding!(mechanical_beam_local_stiffness_matrix,"num_grid_beam_local_stiffness_matrix","crate::sim::physics_fea::BeamElement2D::local_stiffness_matrix","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(element:Handle<Model<fea::BeamElement2D>>)->Dynamic,|ctx,element|{let _=ctx;element.with(ctx,|m|mixed_beam(m.data.local_stiffness_matrix(),&m.frame,false))});
object_binding!(mechanical_beam_global_stiffness_matrix,"num_grid_beam_global_stiffness_matrix","crate::sim::physics_fea::BeamElement2D::global_stiffness_matrix","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(element:Handle<Model<fea::BeamElement2D>>)->Dynamic,|ctx,element|{let _=ctx;element.with(ctx,|m|mixed_beam(m.data.global_stiffness_matrix(),&m.frame,false))});
object_binding!(mechanical_beam_transformation,"num_grid_beam_transformation","crate::sim::physics_fea::BeamElement2D::transformation_matrix","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(element:Handle<Model<fea::BeamElement2D>>)->Dynamic,|ctx,element|{let _=ctx;element.with(ctx,|m|matrix(m.data.transformation_matrix(),&m.frame,[0,0,0,0]))});
object_binding!(mechanical_beam_mass_matrix,"num_grid_beam_mass_matrix","crate::sim::physics_fea::BeamElement2D::mass_matrix","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(element:Handle<Model<fea::BeamElement2D>>,density:f64)->Dynamic,|ctx,element,density|{let _=ctx;if density<0.0{return Err(err("density must be nonnegative"));}element.with(ctx,|m|mixed_beam(m.data.mass_matrix(density),&m.frame,true))});
object_binding!(mechanical_triangle_area,"num_grid_triangle_area","crate::sim::physics_fea::TriangleElement2D::area","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(element:Handle<Model<fea::TriangleElement2D>>)->Dynamic,|ctx,element|{let _=ctx;element.with(ctx,|m|scalar(m.data.area(),&m.frame,[2, 0, 0, 0]))});
object_binding!(mechanical_thermal_triangle_area,"num_grid_thermal_triangle_area","crate::sim::physics_fea::ThermalTriangle2D::area","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(element:Handle<Model<fea::ThermalTriangle2D>>)->Dynamic,|ctx,element|{let _=ctx;element.with(ctx,|m|scalar(m.data.area(),&m.frame,[2, 0, 0, 0]))});
object_binding!(mechanical_tetrahedron_volume,"num_grid_tetrahedron_volume","crate::sim::physics_fea::TetrahedronElement3D::volume","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(element:Handle<Model<fea::TetrahedronElement3D>>)->Dynamic,|ctx,element|{let _=ctx;element.with(ctx,|m|scalar(m.data.volume(),&m.frame,[3, 0, 0, 0]))});
object_binding!(mechanical_tetrahedron_signed_volume,"num_grid_tetrahedron_signed_volume","crate::sim::physics_fea::TetrahedronElement3D::signed_volume","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(element:Handle<Model<fea::TetrahedronElement3D>>)->Dynamic,|ctx,element|{let _=ctx;element.with(ctx,|m|scalar(m.data.signed_volume(),&m.frame,[3, 0, 0, 0]))});
object_binding!(mechanical_tetrahedron_dof_map,"num_grid_tetrahedron_dof_map","crate::sim::physics_fea::TetrahedronElement3D::dof_map","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(element:Handle<Model<fea::TetrahedronElement3D>>)->Vec<usize>,|ctx,element|{let _=ctx;element.with(ctx,|m|Ok(m.data.dof_map().to_vec()))});
object_binding!(mechanical_tetrahedron_shape_gradients,"num_grid_tetrahedron_shape_gradients","crate::sim::physics_fea::TetrahedronElement3D::shape_gradients","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(element:Handle<Model<fea::TetrahedronElement3D>>)->Dynamic,|ctx,element|{let _=ctx;element.with(ctx,|m|{let a=m.data.shape_gradients();let rows=a.len();let cols=a[0].len();matrix(Matrix::new(rows,cols,a.into_iter().flatten().collect()),&m.frame,[-1, 0, 0, 0])})});
object_binding!(mechanical_tetrahedron_elasticity_matrix,"num_grid_tetrahedron_elasticity_matrix","crate::sim::physics_fea::TetrahedronElement3D::elasticity_matrix","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(element:Handle<Model<fea::TetrahedronElement3D>>)->Dynamic,|ctx,element|{let _=ctx;element.with(ctx,|m|{let a=m.data.elasticity_matrix();let rows=a.len();let cols=a[0].len();matrix(Matrix::new(rows,cols,a.into_iter().flatten().collect()),&m.frame,[-1, -2, 1, 0])})});
object_binding!(mechanical_triangle_compute_stress,"num_grid_triangle_compute_stress","crate::sim::physics_fea::TriangleElement2D::compute_stress","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(element:Handle<Model<fea::TriangleElement2D>>,displacements:Vec<f64>)->Dynamic,|ctx,element,displacements|{let _=ctx;if displacements.len()!=6{return Err(err("triangle displacement vector needs six entries"));}element.with(ctx,|m|{let out=m.data.compute_stress(&displacements);let(s,d)=m.frame.derived([-1,-2,1,0])?;Ok(Dynamic(api::record(vec![("data",out.encode()?),("reference",api::quantity(s,d)?)])))})});
object_binding!(mechanical_triangle_von_mises,"num_grid_triangle_von_mises","crate::sim::physics_fea::TriangleElement2D::von_mises_stress","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(frame:Handle<References>,stress:Vec<f64>)->Dynamic,|ctx,frame,stress|{let _=ctx;if stress.len()!=3{return Err(err("plane stress vector needs three entries"));}frame.with(ctx,|f|scalar(fea::TriangleElement2D::von_mises_stress(&stress),f,[-1,-2,1,0]))});
object_binding!(mechanical_max_shear_stress,"num_grid_max_shear_stress","crate::sim::physics_fea::max_shear_stress","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(frame:Handle<References>,sigma1:f64,sigma2:f64)->Dynamic,|ctx,frame,sigma1,sigma2|{let _=ctx;frame.with(ctx,|f|scalar(fea::max_shear_stress(sigma1,sigma2),f,[-1,-2,1,0]))});
object_binding!(mechanical_principal_stresses,"num_grid_principal_stresses","crate::sim::physics_fea::principal_stresses","Native finite-element model in an explicit Length/Time/Mass/Temperature frame. All numeric coefficients/coordinates/displacements are normalized by the corresponding derived references. Results retain physical reference metadata; beam matrices expose heterogeneous translation/rotation DOF references and their quadratic-form reference.",true,(frame:Handle<References>,stress:Vec<f64>)->Dynamic,|ctx,frame,stress|{let _=ctx;if stress.len()!=3{return Err(err("plane stress vector needs three entries"));}frame.with(ctx,|f|{let(a,b,c)=fea::principal_stresses(&stress);let(s,d)=f.derived([-1,-2,1,0])?;Ok(Dynamic(api::record(vec![("sigma1",api::quantity(a*s,d.clone())?),("sigma2",api::quantity(b*s,d)?),("angle",c.encode()?)])))})});
object_binding!(mechanical_safety_factor_von_mises,"num_grid_safety_factor_von_mises","crate::sim::physics_fea::safety_factor_von_mises","Native dimensionless safety factor from three normalized plane stresses and yield strength in the same stress reference; zero equivalent stress preserves the native positive-infinite result.",true,(stress:Vec<f64>,yield_strength:f64)->Dynamic,|ctx,stress,yield_strength|{let _=ctx;if stress.len()!=3||yield_strength<=0.0{return Err(err("stress size or yield strength invalid"));}Ok(Dynamic(Value::Float(fea::safety_factor_von_mises(&stress,yield_strength))))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(i, BINDINGS, "numerical_vector")
}
static BINDINGS: &[&api::Binding] = &[
    &mechanical_assemble_stiffness_matrix::BINDING,
    &mechanical_assemble_1d_stiffness_matrix::BINDING,
    &mechanical_assemble_global_stiffness_matrix::BINDING,
    &mechanical_assemble_2d_stiffness_matrix::BINDING,
    &mechanical_assemble_3d_stiffness_matrix::BINDING,
    &mechanical_mechanical_system_new::BINDING,
    &mechanical_mechanical_system_data::BINDING,
    &mechanical_solve_static_structural::BINDING,
    &mechanical_apply_boundary_conditions_penalty::BINDING,
    &mechanical_compute_element_strain::BINDING,
    &mechanical_create_rectangular_mesh::BINDING,
    &mechanical_refine_mesh::BINDING,
    &mechanical_fea_mesh_data::BINDING,
    &mechanical_mechanical_references::BINDING,
    &mechanical_material_new::BINDING,
    &mechanical_material_steel::BINDING,
    &mechanical_material_aluminum::BINDING,
    &mechanical_material_copper::BINDING,
    &mechanical_material_shear_modulus::BINDING,
    &mechanical_material_bulk_modulus::BINDING,
    &mechanical_node2_new::BINDING,
    &mechanical_node3_new::BINDING,
    &mechanical_node2_distance::BINDING,
    &mechanical_line_element_new::BINDING,
    &mechanical_beam_new::BINDING,
    &mechanical_thermal_line_new::BINDING,
    &mechanical_thermal_triangle_new::BINDING,
    &mechanical_triangle_new::BINDING,
    &mechanical_tetrahedron_new::BINDING,
    &mechanical_line_local_stiffness_matrix::BINDING,
    &mechanical_triangle_local_stiffness_matrix::BINDING,
    &mechanical_triangle_b_matrix::BINDING,
    &mechanical_triangle_constitutive_matrix::BINDING,
    &mechanical_tetrahedron_local_stiffness_matrix::BINDING,
    &mechanical_tetrahedron_b_matrix::BINDING,
    &mechanical_thermal_line_conductivity_matrix::BINDING,
    &mechanical_thermal_triangle_conductivity_matrix::BINDING,
    &mechanical_beam_local_stiffness_matrix::BINDING,
    &mechanical_beam_global_stiffness_matrix::BINDING,
    &mechanical_beam_transformation::BINDING,
    &mechanical_beam_mass_matrix::BINDING,
    &mechanical_triangle_area::BINDING,
    &mechanical_thermal_triangle_area::BINDING,
    &mechanical_tetrahedron_volume::BINDING,
    &mechanical_tetrahedron_signed_volume::BINDING,
    &mechanical_tetrahedron_dof_map::BINDING,
    &mechanical_tetrahedron_shape_gradients::BINDING,
    &mechanical_tetrahedron_elasticity_matrix::BINDING,
    &mechanical_triangle_compute_stress::BINDING,
    &mechanical_triangle_von_mises::BINDING,
    &mechanical_max_shear_stress::BINDING,
    &mechanical_principal_stresses::BINDING,
    &mechanical_safety_factor_von_mises::BINDING,
];

#[derive(Clone)]
struct MechanicalSystem {
    data: Matrix<f64>,
    frame: References,
    dofs: Vec<(f64, DimExpr)>,
}
impl NativeType for MechanicalSystem {
    const KIND: &'static str = "rssn/native/mechanical-system/v1";
    const MAX_BYTES: usize = 4 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        let degree_count = self.data.rows();
        if degree_count > 64
            || degree_count != self.data.cols()
            || self.dofs.len() != degree_count
            || self.data.data().iter().any(|x| !x.is_finite())
            || self.dofs.iter().any(|(s, _)| !s.is_finite() || *s <= 0.0)
        {
            return Err(err(
                "invalid finite mechanical system shape/reference budget",
            ));
        }
        Ok(())
    }
}
fn dof_references(
    v: &Input,
    n: usize,
) -> Result<Vec<(f64, DimExpr)>, OpError> {
    let Value::Record(r) = &v.0 else {
        return Err(err("DOF references must be a record named dof0..dofN"));
    };
    if n > 64 || r.fields.len() != n {
        return Err(err("DOF reference count/allocation budget mismatch"));
    }
    (0..n)
        .map(|i| api::reference(api::field(&v.0, &format!("dof{i}"))?))
        .collect()
}
fn reference_record(dofs: &[(f64, DimExpr)]) -> Result<Value, OpError> {
    Ok(Value::Record(crate::kernel::value::RecordValue {
        fields: dofs
            .iter()
            .enumerate()
            .map(|(i, (s, d))| Ok((format!("dof{i}"), api::quantity(*s, d.clone())?)))
            .collect::<Result<Vec<_>, OpError>>()?,
    }))
}
fn system_output(s: &MechanicalSystem) -> Result<Dynamic, OpError> {
    let (e, d) = s.frame.derived([2, -2, 1, 0])?;
    Ok(Dynamic(api::record(vec![
        ("rows", s.data.rows().encode()?),
        ("data", s.data.data().clone().encode()?),
        ("dof_references", reference_record(&s.dofs)?),
        ("energy_reference", api::quantity(e, d)?),
    ])))
}
type NativeElement = (Matrix<f64>, Vec<usize>);
fn elements(
    input: &Input,
    n: usize,
    fixed: Option<usize>,
) -> Result<Vec<NativeElement>, OpError> {
    api::normalized(&input.0)?;
    let Value::Tensor(t) = &input.0 else {
        return Err(err("elements must be a record vector"));
    };
    if t.shape != [t.data.len()] || t.data.len() > 128 || n > 64 {
        return Err(err("assembly element/DOF allocation budget exceeded"));
    }
    let mut out = Vec::new();
    for value in &t.data {
        let Value::Record(record) = value else {
            return Err(err("assembly element must contain matrix/dofs"));
        };
        if record.fields.len() != 2 {
            return Err(err("assembly element requires exactly matrix and dofs"));
        }
        let m = crate::rules::numerical_abi::matrix(api::field(value, "matrix")?, 1)?;
        let map = <Vec<usize> as super::super::codecs::Arg>::decode(api::field(value, "dofs")?, 1)?;
        if m.rows > 64
            || m.cols > 64
            || m.rows < map.len()
            || m.cols < map.len()
            || fixed.is_some_and(|k| map.len() != k)
            || map.iter().any(|i| *i >= n)
        {
            return Err(err("local matrix/DOF map shape or index invalid"));
        }
        out.push((Matrix::new(m.rows, m.cols, m.data), map));
    }
    Ok(out)
}
#[derive(Clone)]
struct Mesh {
    nodes: Vec<fea::Node2D>,
    triangles: Vec<[usize; 3]>,
    frame: References,
}
impl NativeType for Mesh {
    const KIND: &'static str = "rssn/native/fea-triangle-mesh/v1";
    const MAX_BYTES: usize = 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        if self.nodes.len() > 4096
            || self.triangles.len() > 8192
            || self
                .triangles
                .iter()
                .flatten()
                .any(|i| *i >= self.nodes.len())
        {
            return Err(err("finite-element mesh allocation/index budget exceeded"));
        }
        for n in &self.nodes {
            n.valid()?;
        }
        Ok(())
    }
}
object_binding!(mechanical_assemble_stiffness_matrix,"num_grid_assemble_stiffness_matrix","crate::sim::physics_fea::assemble_stiffness_matrix","Actual native finite-element assembly/structural method in an explicit mechanical frame. Matrices and loads are normalized in the declared DOF references; the energy reference and heterogeneous DOF references determine every physical coefficient. Assembly, geometry, bounds and mutation results are checked before retention.",false,(frame:Handle<References>,num_dofs:usize,local_elements:Input,dof_references_value:Input)->NewObject<MechanicalSystem>,|ctx,frame,num_dofs,local_elements,dof_references_value|{let _=ctx;let dofs=dof_references(&dof_references_value,num_dofs)?;let elements=elements(&local_elements,num_dofs,None)?;frame.with(ctx,|f|{Ok(NewObject::new(MechanicalSystem{data:fea::assemble_stiffness_matrix(num_dofs,&elements),frame:f.clone(),dofs:dofs.clone()},dofs.into_iter().map(|(_,d)|d).collect()))})});
object_binding!(mechanical_assemble_1d_stiffness_matrix,"num_grid_assemble_1d_stiffness_matrix","crate::sim::physics_fea::assemble_1d_stiffness_matrix","Actual native finite-element assembly/structural method in an explicit mechanical frame. Matrices and loads are normalized in the declared DOF references; the energy reference and heterogeneous DOF references determine every physical coefficient. Assembly, geometry, bounds and mutation results are checked before retention.",false,(frame:Handle<References>,num_dofs:usize,local_elements:Input,dof_references_value:Input)->NewObject<MechanicalSystem>,|ctx,frame,num_dofs,local_elements,dof_references_value|{let _=ctx;let dofs=dof_references(&dof_references_value,num_dofs)?;let elements=elements(&local_elements,num_dofs,Some(2))?;frame.with(ctx,|f|{let elements=elements.into_iter().map(|(m,d)|(m,d[0],d[1])).collect::<Vec<_>>();Ok(NewObject::new(MechanicalSystem{data:fea::assemble_1d_stiffness_matrix(num_dofs,&elements),frame:f.clone(),dofs:dofs.clone()},dofs.into_iter().map(|(_,d)|d).collect()))})});
object_binding!(mechanical_assemble_global_stiffness_matrix,"num_grid_assemble_global_stiffness_matrix","crate::sim::physics_fea::assemble_global_stiffness_matrix","Actual native finite-element assembly/structural method in an explicit mechanical frame. Matrices and loads are normalized in the declared DOF references; the energy reference and heterogeneous DOF references determine every physical coefficient. Assembly, geometry, bounds and mutation results are checked before retention.",false,(frame:Handle<References>,num_dofs:usize,local_elements:Input,dof_references_value:Input)->NewObject<MechanicalSystem>,|ctx,frame,num_dofs,local_elements,dof_references_value|{let _=ctx;let dofs=dof_references(&dof_references_value,num_dofs)?;let elements=elements(&local_elements,num_dofs,Some(2))?;frame.with(ctx,|f|{let elements=elements.into_iter().map(|(m,d)|(m,d[0],d[1])).collect::<Vec<_>>();Ok(NewObject::new(MechanicalSystem{data:fea::assemble_global_stiffness_matrix(num_dofs,&elements),frame:f.clone(),dofs:dofs.clone()},dofs.into_iter().map(|(_,d)|d).collect()))})});
object_binding!(mechanical_assemble_2d_stiffness_matrix,"num_grid_assemble_2d_stiffness_matrix","crate::sim::physics_fea::assemble_2d_stiffness_matrix","Actual native finite-element assembly/structural method in an explicit mechanical frame. Matrices and loads are normalized in the declared DOF references; the energy reference and heterogeneous DOF references determine every physical coefficient. Assembly, geometry, bounds and mutation results are checked before retention.",false,(frame:Handle<References>,num_dofs:usize,local_elements:Input,dof_references_value:Input)->NewObject<MechanicalSystem>,|ctx,frame,num_dofs,local_elements,dof_references_value|{let _=ctx;let dofs=dof_references(&dof_references_value,num_dofs)?;let elements=elements(&local_elements,num_dofs,Some(6))?;frame.with(ctx,|f|{let elements=elements.into_iter().map(|(m,d)|Ok((m,<Vec<usize> as TryInto<[usize;6]>>::try_into(d).map_err(|_|err("DOF shape mismatch"))?))).collect::<Result<Vec<_>,OpError>>()?;Ok(NewObject::new(MechanicalSystem{data:fea::assemble_2d_stiffness_matrix(num_dofs,&elements),frame:f.clone(),dofs:dofs.clone()},dofs.into_iter().map(|(_,d)|d).collect()))})});
object_binding!(mechanical_assemble_3d_stiffness_matrix,"num_grid_assemble_3d_stiffness_matrix","crate::sim::physics_fea::assemble_3d_stiffness_matrix","Actual native finite-element assembly/structural method in an explicit mechanical frame. Matrices and loads are normalized in the declared DOF references; the energy reference and heterogeneous DOF references determine every physical coefficient. Assembly, geometry, bounds and mutation results are checked before retention.",false,(frame:Handle<References>,num_dofs:usize,local_elements:Input,dof_references_value:Input)->NewObject<MechanicalSystem>,|ctx,frame,num_dofs,local_elements,dof_references_value|{let _=ctx;let dofs=dof_references(&dof_references_value,num_dofs)?;let elements=elements(&local_elements,num_dofs,Some(12))?;frame.with(ctx,|f|{let elements=elements.into_iter().map(|(m,d)|Ok((m,<Vec<usize> as TryInto<[usize;12]>>::try_into(d).map_err(|_|err("DOF shape mismatch"))?))).collect::<Result<Vec<_>,OpError>>()?;Ok(NewObject::new(MechanicalSystem{data:fea::assemble_3d_stiffness_matrix(num_dofs,&elements),frame:f.clone(),dofs:dofs.clone()},dofs.into_iter().map(|(_,d)|d).collect()))})});
object_binding!(mechanical_mechanical_system_new,"num_grid_mechanical_system_new","","Actual native finite-element assembly/structural method in an explicit mechanical frame. Matrices and loads are normalized in the declared DOF references; the energy reference and heterogeneous DOF references determine every physical coefficient. Assembly, geometry, bounds and mutation results are checked before retention.",false,(frame:Handle<References>,matrix_value:Input,dof_references_value:Input)->NewObject<MechanicalSystem>,|ctx,frame,matrix_value,dof_references_value|{let _=ctx;api::normalized(&matrix_value.0)?;let m=crate::rules::numerical_abi::matrix(&matrix_value.0,1)?;if m.rows!=m.cols||m.rows>64{return Err(err("mechanical system must be bounded square matrix"));}let dofs=dof_references(&dof_references_value,m.rows)?;frame.with(ctx,|f|Ok(NewObject::new(MechanicalSystem{data:Matrix::new(m.rows,m.cols,m.data),frame:f.clone(),dofs:dofs.clone()},dofs.into_iter().map(|(_,d)|d).collect())))});
object_binding!(mechanical_mechanical_system_data,"num_grid_mechanical_system_data","","Actual native finite-element assembly/structural method in an explicit mechanical frame. Matrices and loads are normalized in the declared DOF references; the energy reference and heterogeneous DOF references determine every physical coefficient. Assembly, geometry, bounds and mutation results are checked before retention.",true,(system:Handle<MechanicalSystem>)->Dynamic,|ctx,system|{let _=ctx;system.with(ctx,system_output)});
object_binding!(mechanical_solve_static_structural,"num_grid_solve_static_structural","crate::sim::physics_fea::solve_static_structural","Actual native finite-element assembly/structural method in an explicit mechanical frame. Matrices and loads are normalized in the declared DOF references; the energy reference and heterogeneous DOF references determine every physical coefficient. Assembly, geometry, bounds and mutation results are checked before retention.",true,(system:Handle<MechanicalSystem>,forces:Vec<f64>,fixed_indices:Vec<usize>,fixed_values:Vec<f64>)->Dynamic,|ctx,system,forces,fixed_indices,fixed_values|{let _=ctx;system.with(ctx,|s|{if forces.len()!=s.dofs.len()||fixed_indices.len()!=fixed_values.len()||fixed_indices.iter().any(|i|*i>=s.dofs.len()){return Err(err("structural load/constraint shape or index invalid"));}let fixed=fixed_indices.into_iter().zip(fixed_values).collect::<Vec<_>>();let out=fea::solve_static_structural(s.data.clone(),forces,&fixed).map_err(err)?;Ok(Dynamic(api::record(vec![("displacements",out.encode()?),("references",reference_record(&s.dofs)?)])))})});
object_binding!(mechanical_apply_boundary_conditions_penalty,"num_grid_apply_boundary_conditions_penalty","crate::sim::physics_fea::apply_boundary_conditions_penalty","Actual native finite-element assembly/structural method in an explicit mechanical frame. Matrices and loads are normalized in the declared DOF references; the energy reference and heterogeneous DOF references determine every physical coefficient. Assembly, geometry, bounds and mutation results are checked before retention.",true,(system:Handle<MechanicalSystem>,forces:Vec<f64>,fixed_indices:Vec<usize>,fixed_values:Vec<f64>,penalty:f64)->Dynamic,|ctx,system,forces,fixed_indices,fixed_values,penalty|{let _=ctx;if penalty<=0.0{return Err(err("penalty must be positive"));}let mut candidate=system.with(ctx,|s|Ok(s.clone()))?;if forces.len()!=candidate.dofs.len()||fixed_indices.len()!=fixed_values.len()||fixed_indices.iter().any(|i|*i>=candidate.dofs.len()){return Err(err("structural load/constraint shape or index invalid"));}let fixed=fixed_indices.into_iter().zip(fixed_values).collect::<Vec<_>>();let mut forces=forces;fea::apply_boundary_conditions_penalty(&mut candidate.data,&mut forces,&fixed,penalty);candidate.validate()?;if forces.iter().any(|v|!v.is_finite()){return Err(err("penalty force arithmetic overflow"));}let result=Dynamic(api::record(vec![("forces",forces.encode()?),("system",system_output(&candidate)?.0)]));system.with_mut(ctx,|s|{*s=candidate;Ok(())})?;Ok(result)});
object_binding!(mechanical_compute_element_strain,"num_grid_compute_element_strain","crate::sim::physics_fea::compute_element_strain","Native normalized strain-coordinate contraction; B and displacement columns use one explicitly normalized coordinate basis, returning dimensionless strain coordinates. Shapes and work are bounded.",true,(b_matrix:Vec<Vec<f64>>,displacements:Vec<f64>)->Vec<f64>,|ctx,b_matrix,displacements|{let _=ctx;if b_matrix.len()>16||b_matrix.iter().any(|r|r.len()!=displacements.len())||displacements.len()>64{return Err(err("strain matrix/displacement shape budget invalid"));}let rows=b_matrix.len();let cols=displacements.len();let m=Matrix::new(rows,cols,b_matrix.into_iter().flatten().collect());Ok(fea::compute_element_strain(&m,&displacements))});
object_binding!(mechanical_create_rectangular_mesh,"num_grid_create_rectangular_mesh","crate::sim::physics_fea::create_rectangular_mesh","Actual native finite-element assembly/structural method in an explicit mechanical frame. Matrices and loads are normalized in the declared DOF references; the energy reference and heterogeneous DOF references determine every physical coefficient. Assembly, geometry, bounds and mutation results are checked before retention.",true,(frame:Handle<References>,width:f64,height:f64,nx:usize,ny:usize)->NewObject<Mesh>,|ctx,frame,width,height,nx,ny|{let _=ctx;if width<=0.0||height<=0.0||nx==0||ny==0||nx>64||ny>64||(nx+1).checked_mul(ny+1).is_none_or(|n|n>4096){return Err(err("rectangular mesh geometry/allocation budget invalid"));}frame.with(ctx,|f|{let(nodes,triangles)=fea::create_rectangular_mesh(width,height,nx,ny);Ok(NewObject::new(Mesh{nodes,triangles,frame:f.clone()},f.dims.to_vec()))})});
object_binding!(mechanical_refine_mesh,"num_grid_refine_mesh","crate::sim::physics_fea::refine_mesh","Actual native finite-element assembly/structural method in an explicit mechanical frame. Matrices and loads are normalized in the declared DOF references; the energy reference and heterogeneous DOF references determine every physical coefficient. Assembly, geometry, bounds and mutation results are checked before retention.",true,(mesh:Handle<Mesh>)->NewObject<Mesh>,|ctx,mesh|{let _=ctx;mesh.with(ctx,|m|{if m.triangles.len()>2048||m.nodes.len()+3*m.triangles.len()>4096{return Err(err("refined mesh allocation budget exceeded"));}let(nodes,triangles)=fea::refine_mesh(&m.nodes,&m.triangles);Ok(NewObject::new(Mesh{nodes,triangles,frame:m.frame.clone()},m.frame.dims.to_vec()))})});
object_binding!(mechanical_fea_mesh_data,"num_grid_fea_mesh_data","","Actual native finite-element assembly/structural method in an explicit mechanical frame. Matrices and loads are normalized in the declared DOF references; the energy reference and heterogeneous DOF references determine every physical coefficient. Assembly, geometry, bounds and mutation results are checked before retention.",true,(mesh:Handle<Mesh>)->Dynamic,|ctx,mesh|{let _=ctx;mesh.with(ctx,|m|Ok(Dynamic(api::record(vec![("nodes",m.nodes.iter().map(|n|[n.x,n.y]).collect::<Vec<_>>().encode()?),("triangles",m.triangles.clone().encode()?),("length_reference",api::quantity(m.frame.scales[0],m.frame.dims[0].clone())?)]))))});
