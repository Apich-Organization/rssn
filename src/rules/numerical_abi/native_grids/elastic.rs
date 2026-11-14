//! Native small/large deformation models with separate mechanical dimensions.
use super::super::native_object_api::Many;
use super::super::native_object_api::Maybe;
use super::CallCtx;
use super::CompiledSnapshot;
use super::DimExpr;
use super::Dynamic;
use super::Exp;
use super::Frame;
use super::Handle;
use super::Input;
use super::Installer;
use super::NativeArg;
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
use super::boundary;
use super::err;
use super::fem::Mesh;
use super::fem::Space;
use super::mechanics::References;
use super::object_binding;
use super::structured;
use super::structured_ops::checked_snapshot;
use crate::sim::elasticity_fem as el;
use crate::sim::hyperelastic as hyp;
#[derive(Clone)]
struct Bcs {
    data: boundary::MarkerBcs<boundary::ElasticBoundary>,
    refs: References,
}
impl NativeType for Bcs {
    const KIND: &'static str = "rssn/native/elastic-marker-bcs/v1";
    const MAX_BYTES: usize = 1048576;

    fn validate(&self) -> Result<(), OpError> {
        self.refs.validate()
    }
}
#[derive(Clone)]
struct Material {
    data: hyp::HyperMaterial,
    refs: References,
}
impl NativeType for Material {
    const KIND: &'static str = "rssn/native/hyper-material/v1";
    const MAX_BYTES: usize = 4096;

    fn validate(&self) -> Result<(), OpError> {
        self.refs.validate()?;
        material(self.data.young, self.data.poisson)?;
        if self.data.density <= 0.0 || !self.data.density.is_finite() {
            return Err(err("positive density required"));
        }
        Ok(())
    }
}
struct Solution {
    data: el::OwnedElasticSolution,
    refs: References,
}
impl NativeType for Solution {
    const KIND: &'static str = "rssn/native/elastic-solution/v1";
    const MAX_BYTES: usize = 16 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.refs.validate()?;
        self.data.with(|s| {
            if s.u.len() > 1024 || s.u.iter().any(|u| !u.is_finite()) {
                return Err(err("elastic solution invariant failed"));
            }
            Ok(())
        })
    }
}
#[derive(Clone)]
struct Contact {
    marker: i32,
    gap: CompiledSnapshot,
    normal: [CompiledSnapshot; 3],
    penalty: f64,
    refs: References,
}
impl NativeType for Contact {
    const KIND: &'static str = "rssn/native/penalty-contact/v1";
    const MAX_BYTES: usize = 262144;

    fn validate(&self) -> Result<(), OpError> {
        self.refs.validate()?;
        if std::iter::once(&self.gap)
            .chain(self.normal.iter())
            .try_fold(0usize, |n, s| n.checked_add(s.retained_bytes))
            .is_none_or(|n| n > 131072)
        {
            return Err(err("penalty contact compiled storage cap exceeded"));
        }
        if self.penalty <= 0.0 || !self.penalty.is_finite() {
            return Err(err("positive contact stiffness required"));
        }
        Ok(())
    }
}
#[derive(Clone)]
struct Hyper {
    mesh: Mesh,
    material: Material,
    bcs: Bcs,
    force: [CompiledSnapshot; 3],
    contacts: Vec<Contact>,
}
impl NativeType for Hyper {
    const KIND: &'static str = "rssn/native/hyper-problem/v1";
    const MAX_BYTES: usize = 16 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.mesh.validate()?;
        self.material.validate()?;
        check_refs(&self.material.refs, &self.bcs.refs)?;
        check_frame(&self.mesh.frame, &self.material.refs)?;
        let retained = self
            .force
            .iter()
            .chain(
                self.contacts
                    .iter()
                    .flat_map(|c| std::iter::once(&c.gap).chain(c.normal.iter())),
            )
            .try_fold(0usize, |n, s| n.checked_add(s.retained_bytes));
        if retained.is_none_or(|n| n > 4194304) {
            return Err(err("hyperelastic compiled storage cap exceeded"));
        }
        if self.mesh.data.dim < 2
            || self.mesh.data.nodes.len() * self.mesh.data.dim > 384
            || self.mesh.data.cells.len() > 512
            || self.contacts.len() > 16
        {
            return Err(err("hyperelastic DOF/contact budget exceeded"));
        }
        Ok(())
    }
}
struct State {
    data: hyp::DynamicState,
    problem: Hyper,
}
impl NativeType for State {
    const KIND: &'static str = "rssn/native/hyper-state/v1";
    const MAX_BYTES: usize = 20 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.problem.validate()?;
        let n = self.problem.mesh.data.nodes.len() * self.problem.mesh.data.dim;
        if self.data.u.len() != n
            || self.data.v.len() != n
            || self.data.a.len() != n
            || self
                .data
                .u
                .iter()
                .chain(&self.data.v)
                .chain(&self.data.a)
                .any(|x| !x.is_finite())
            || !self.data.t.is_finite()
        {
            return Err(err("hyperelastic dynamic state invalid"));
        }
        Ok(())
    }
}
fn material(
    e: f64,
    nu: f64,
) -> Result<(), OpError> {
    if e <= 0.0 || !e.is_finite() || !nu.is_finite() || nu <= -1.0 || nu >= 0.5 {
        return Err(err("material requires E>0 and -1<nu<1/2"));
    }
    Ok(())
}
fn model(s: &str) -> Result<el::Model, OpError> {
    match s {
        | "PlaneStress" => Ok(el::Model::PlaneStress),
        | "PlaneStrain" => Ok(el::Model::PlaneStrain),
        | "Solid3d" => Ok(el::Model::Solid3d),
        | _ => {
            Err(err(
                "elastic model must be PlaneStress, PlaneStrain or Solid3d",
            ))
        },
    }
}
fn check_refs(
    a: &References,
    b: &References,
) -> Result<(), OpError> {
    for (x, y) in a.scales.iter().zip(b.scales) {
        if x.to_bits() != y.to_bits() {
            return Err(err("mechanical references differ"));
        }
    }
    for (x, y) in a.dims.iter().zip(&b.dims) {
        api::same_dimension(x, y)?;
    }
    Ok(())
}
fn check_frame(
    a: &Frame,
    b: &References,
) -> Result<(), OpError> {
    api::same_dimension(&a.dims[0], &b.dims[0])?;
    api::same_dimension(&a.dims[2], &b.dims[0])?;
    if a.length.to_bits() != b.scales[0].to_bits() || a.value.to_bits() != b.scales[0].to_bits() {
        return Err(err(
            "elastic mesh displacement/coordinate references must equal mechanical Length reference",
        ));
    }
    Ok(())
}
fn fields(
    f: &References,
    data: Vec<f64>,
    powers: [i64; 4],
) -> Result<Dynamic, OpError> {
    let (s, d) = f.derived(powers)?;
    Ok(Dynamic(api::record(vec![
        ("normalized", data.encode()?),
        ("reference", api::quantity(s, d)?),
    ])))
}
fn mat3(v: Vec<Vec<f64>>) -> Result<[[f64; 3]; 3], OpError> {
    let rows = v
        .into_iter()
        .map(|r| {
            r.try_into()
                .map_err(|_| err("3x3 deformation gradient required"))
        })
        .collect::<Result<Vec<[f64; 3]>, OpError>>()?;
    rows.try_into()
        .map_err(|_| err("3x3 deformation gradient required"))
}
fn force(
    ctx: &CallCtx<'_>,
    programs: &Many<Input>,
    refs: &References,
    budget: usize,
) -> Result<[CompiledSnapshot; 3], OpError> {
    if programs.0.len() != 3 {
        return Err(err("three compiled body-force components required"));
    }
    let (_, d) = refs.derived([-2, -2, 1, 0])?;
    programs
        .0
        .iter()
        .map(|p| checked_snapshot(ctx, p, budget, &vec![refs.dims[0].clone(); 3], &d))
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .map_err(|_| err("three body-force components required"))
}
fn f_at(
    s: &[CompiledSnapshot; 3],
    refs: &References,
    x: [f64; 3],
) -> [f64; 3] {
    let scale = refs.derived([-2, -2, 1, 0]).map_or(f64::NAN, |x| x.0);
    let x = x.map(|x| x * refs.scales[0]);
    std::array::from_fn(|i| s[i].sample(&x) / scale)
}
fn finish(s: &[CompiledSnapshot; 3]) -> Result<(), OpError> {
    for f in s {
        f.finish()?;
    }
    Ok(())
}
impl Hyper {
    fn with<R>(
        &self,
        run: impl FnOnce(&hyp::HyperProblem<'_>) -> Result<R, OpError>,
    ) -> Result<R, OpError> {
        let refs = &self.material.refs;
        let force = |x| f_at(&self.force, refs, x);
        let contacts = self
            .contacts
            .iter()
            .map(|c| {
                let g = c.gap.clone();
                let n = c.normal.clone();
                let l = refs.scales[0];
                hyp::PenaltyContact {
                    marker: c.marker,
                    penalty: c.penalty,
                    gap: Box::new(move |x| g.sample(&x.map(|x| x * l)) / l),
                    normal: Box::new(move |x| {
                        let x = x.map(|x| x * l);
                        std::array::from_fn(|i| n[i].sample(&x))
                    }),
                }
            })
            .collect();
        let problem = hyp::HyperProblem {
            mesh: &self.mesh.data,
            material: self.material.data,
            bcs: &self.bcs.data,
            force: &force,
            contacts,
        };
        let out = run(&problem);
        finish(&self.force)?;
        for c in &self.contacts {
            c.gap.finish()?;
            for n in &c.normal {
                n.finish()?;
            }
        }
        out
    }
}
object_binding!(elastic_bcs,"num_grid_elastic_bcs","","Native elasticity call in an explicit mechanical reference frame. Numeric displacements/coordinates use Length, body forces use Mass/(Length^2 Time^2), stress/Young modulus use Mass/(Length Time^2), and density uses Mass/Length^3. Reduced 2D hyperelastic energies are per unit out-of-plane thickness. Outputs retain distinct dimensional reference witnesses. Native DOF, nonlinear work, callbacks and allocation budgets are checked.",true,(refs:Handle<References>,markers:Vec<i64>,conditions:Many<Handle<boundary::ElasticBoundary>>)->NewObject<Bcs>,|ctx:&CallCtx<'_>,refs:Handle<References>,markers:Vec<i64>,conditions:Many<Handle<boundary::ElasticBoundary>>|{refs.with(ctx,|r|{if markers.len()!=conditions.0.len()||markers.len()>256{return Err(err("elastic marker size invalid"));}let mut data=boundary::MarkerBcs::new();for(marker,c)in markers.into_iter().zip(conditions.0){let marker=i32::try_from(marker).map_err(err)?;let c=c.with(ctx,|b|Ok(b.clone()))?;data=data.set(marker,c);}Ok(NewObject::new(Bcs{data,refs:r.clone()},r.dims.to_vec()))})});
object_binding!(hyper_material,"num_grid_hyper_material","","Native elasticity call in an explicit mechanical reference frame. Numeric displacements/coordinates use Length, body forces use Mass/(Length^2 Time^2), stress/Young modulus use Mass/(Length Time^2), and density uses Mass/Length^3. Reduced 2D hyperelastic energies are per unit out-of-plane thickness. Outputs retain distinct dimensional reference witnesses. Native DOF, nonlinear work, callbacks and allocation budgets are checked.",true,(refs:Handle<References>,kind:String,young:f64,poisson:f64,density:f64)->NewObject<Material>,|ctx:&CallCtx<'_>,refs:Handle<References>,kind:String,young:f64,poisson:f64,density:f64|{refs.with(ctx,|r|{let model=match kind.as_str(){"StVenantKirchhoff"=>hyp::HyperModel::StVenantKirchhoff,"NeoHookean"=>hyp::HyperModel::NeoHookean,_=>return Err(err("hyperelastic model invalid"))};Ok(NewObject::new(Material{data:hyp::HyperMaterial{model,young,poisson,density},refs:r.clone()},r.dims.to_vec()))})});
object_binding!(elastic_constitutive,"num_grid_elastic_constitutive","crate::sim::elasticity_fem::constitutive_matrix","Native elasticity call in an explicit mechanical reference frame. Numeric displacements/coordinates use Length, body forces use Mass/(Length^2 Time^2), stress/Young modulus use Mass/(Length Time^2), and density uses Mass/Length^3. Reduced 2D hyperelastic energies are per unit out-of-plane thickness. Outputs retain distinct dimensional reference witnesses. Native DOF, nonlinear work, callbacks and allocation budgets are checked.",true,(refs:Handle<References>,kind:String,young:f64,poisson:f64,thickness:f64)->Dynamic,|ctx:&CallCtx<'_>,refs:Handle<References>,kind:String,young:f64,poisson:f64,thickness:f64|{material(young,poisson)?;if thickness<=0.0{return Err(err("positive thickness required"));}refs.with(ctx,|r|{let matrix=el::constitutive_matrix(model(&kind)?,&el::Material{young,poisson,thickness});let n=matrix.len();Ok(Dynamic(api::record(vec![("shape",vec![n,n].encode()?),("entries",fields(r,matrix.into_iter().flatten().collect(),[-1,-2,1,0])?.0)])))})});
object_binding!(elastic_solve,"num_grid_elastic_solve","crate::sim::elasticity_fem::solve","Native elasticity call in an explicit mechanical reference frame. Numeric displacements/coordinates use Length, body forces use Mass/(Length^2 Time^2), stress/Young modulus use Mass/(Length Time^2), and density uses Mass/Length^3. Reduced 2D hyperelastic energies are per unit out-of-plane thickness. Outputs retain distinct dimensional reference witnesses. Native DOF, nonlinear work, callbacks and allocation budgets are checked.",true,(space:Handle<Space>,bcs:Handle<Bcs>,kind:String,young:f64,poisson:f64,thickness:f64,body:Many<Input>,budget:usize)->NewObject<Solution>,|ctx:&CallCtx<'_>,space:Handle<Space>,bcs:Handle<Bcs>,kind:String,young:f64,poisson:f64,thickness:f64,body:Many<Input>,budget:usize|{space.with(ctx,|s|bcs.with(ctx,|bc|{check_frame(&s.frame,&bc.refs)?;if s.data.mesh.dim<2||s.data.size()*s.data.mesh.dim>512||thickness<=0.0{return Err(err("elastic dimension/thickness/DOF budget invalid"));}let model=model(&kind)?;let f=force(ctx,&body,&bc.refs,budget)?;let force=|x|f_at(&f,&bc.refs,x);material(young,poisson)?;let p=el::ElasticityProblem{model,material:el::Material{young,poisson,thickness},force:&force,bcs:&bc.data};let result=el::solve(&s.data,&p,structured::LinearSolver::Auto);finish(&f)?;Ok(NewObject::new(Solution{data:result.map_err(err)?.into_owned(),refs:bc.refs.clone()},bc.refs.dims.to_vec()))}))});
object_binding!(elastic_solve_heterogeneous,"num_grid_elastic_solve_heterogeneous","crate::sim::elasticity_fem::solve_heterogeneous","Native elasticity call in an explicit mechanical reference frame. Numeric displacements/coordinates use Length, body forces use Mass/(Length^2 Time^2), stress/Young modulus use Mass/(Length Time^2), and density uses Mass/Length^3. Reduced 2D hyperelastic energies are per unit out-of-plane thickness. Outputs retain distinct dimensional reference witnesses. Native DOF, nonlinear work, callbacks and allocation budgets are checked.",true,(space:Handle<Space>,bcs:Handle<Bcs>,kind:String,young:Input,poisson:Input,thickness:f64,body:Many<Input>,budget:usize)->NewObject<Solution>,|ctx:&CallCtx<'_>,space:Handle<Space>,bcs:Handle<Bcs>,kind:String,young:Input,poisson:Input,thickness:f64,body:Many<Input>,budget:usize|{space.with(ctx,|s|bcs.with(ctx,|bc|{check_frame(&s.frame,&bc.refs)?;if s.data.mesh.dim<2||s.data.size()*s.data.mesh.dim>512||thickness<=0.0{return Err(err("elastic dimension/thickness/DOF budget invalid"));}let model=model(&kind)?;let f=force(ctx,&body,&bc.refs,budget)?;let force=|x|f_at(&f,&bc.refs,x);let(es,ed)=bc.refs.derived([-1,-2,1,0])?;let input=vec![bc.refs.dims[0].clone();3];let e=checked_snapshot(ctx,&young,budget,&input,&ed)?;let nu=checked_snapshot(ctx,&poisson,budget,&input,&DimExpr::Dimless)?;let at_e=|x:[f64;3]|e.sample(&x.map(|x|x*bc.refs.scales[0]))/es;let at_nu=|x:[f64;3]|nu.sample(&x.map(|x|x*bc.refs.scales[0]));let p=el::ElasticityProblemVar{model,young:&at_e,poisson:&at_nu,thickness,force:&force,bcs:&bc.data};let result=el::solve_heterogeneous(&s.data,&p,structured::LinearSolver::Auto);e.finish()?;nu.finish()?;finish(&f)?;Ok(NewObject::new(Solution{data:result.map_err(err)?.into_owned(),refs:bc.refs.clone()},bc.refs.dims.to_vec()))}))});
object_binding!(elastic_displacement,"num_grid_elastic_displacement","crate::sim::elasticity_fem::ElasticSolution<'_>::displacement","Native elasticity call in an explicit mechanical reference frame. Numeric displacements/coordinates use Length, body forces use Mass/(Length^2 Time^2), stress/Young modulus use Mass/(Length Time^2), and density uses Mass/Length^3. Reduced 2D hyperelastic energies are per unit out-of-plane thickness. Outputs retain distinct dimensional reference witnesses. Native DOF, nonlinear work, callbacks and allocation budgets are checked.",true,(solution:Handle<Solution>,dof:usize)->Dynamic,|ctx:&CallCtx<'_>,solution:Handle<Solution>,dof:usize|{solution.with(ctx,|s|s.data.with(|v|{if dof>=v.space.size(){return Err(err("elastic DOF out of range"));}fields(&s.refs,v.displacement(dof).to_vec(),[1,0,0,0])}))});
object_binding!(elastic_eval,"num_grid_elastic_eval","crate::sim::elasticity_fem::ElasticSolution<'_>::eval","Native elasticity call in an explicit mechanical reference frame. Numeric displacements/coordinates use Length, body forces use Mass/(Length^2 Time^2), stress/Young modulus use Mass/(Length Time^2), and density uses Mass/Length^3. Reduced 2D hyperelastic energies are per unit out-of-plane thickness. Outputs retain distinct dimensional reference witnesses. Native DOF, nonlinear work, callbacks and allocation budgets are checked.",true,(solution:Handle<Solution>,point:Vec<f64>)->Maybe<Dynamic>,|ctx:&CallCtx<'_>,solution:Handle<Solution>,point:Vec<f64>|{let point:[f64;3]=point.try_into().map_err(|_|err("three coordinates required"))?;solution.with(ctx,|s|s.data.with(|v|Ok(Maybe(v.eval(point).map(|u|fields(&s.refs,u.to_vec(),[1,0,0,0])).transpose()?))))});
object_binding!(elastic_stress,"num_grid_elastic_stress","crate::sim::elasticity_fem::ElasticSolution<'_>::stress_at_cell","Native elasticity call in an explicit mechanical reference frame. Numeric displacements/coordinates use Length, body forces use Mass/(Length^2 Time^2), stress/Young modulus use Mass/(Length Time^2), and density uses Mass/Length^3. Reduced 2D hyperelastic energies are per unit out-of-plane thickness. Outputs retain distinct dimensional reference witnesses. Native DOF, nonlinear work, callbacks and allocation budgets are checked.",true,(solution:Handle<Solution>,cell:usize)->Dynamic,|ctx:&CallCtx<'_>,solution:Handle<Solution>,cell:usize|{solution.with(ctx,|s|s.data.with(|v|{if cell>=v.space.mesh.cells.len(){return Err(err("elastic cell out of range"));}fields(&s.refs,v.stress_at_cell(cell),[-1,-2,1,0])}))});
object_binding!(elastic_von_mises,"num_grid_elastic_von_mises","crate::sim::elasticity_fem::ElasticSolution<'_>::von_mises","Native elasticity call in an explicit mechanical reference frame. Numeric displacements/coordinates use Length, body forces use Mass/(Length^2 Time^2), stress/Young modulus use Mass/(Length Time^2), and density uses Mass/Length^3. Reduced 2D hyperelastic energies are per unit out-of-plane thickness. Outputs retain distinct dimensional reference witnesses. Native DOF, nonlinear work, callbacks and allocation budgets are checked.",true,(refs:Handle<References>,stress:Vec<f64>)->Dynamic,|ctx:&CallCtx<'_>,refs:Handle<References>,stress:Vec<f64>|{if stress.len()!=3&&stress.len()!=6{return Err(err("stress Voigt tensor needs 3 or 6 entries"));}refs.with(ctx,|r|{let(s,d)=r.derived([-1,-2,1,0])?;Ok(Dynamic(api::quantity(el::ElasticSolution::von_mises(&stress)*s,d)?))})});
object_binding!(hyper_lame,"num_grid_hyper_lame","crate::sim::hyperelastic::HyperMaterial::lame","Native elasticity call in an explicit mechanical reference frame. Numeric displacements/coordinates use Length, body forces use Mass/(Length^2 Time^2), stress/Young modulus use Mass/(Length Time^2), and density uses Mass/Length^3. Reduced 2D hyperelastic energies are per unit out-of-plane thickness. Outputs retain distinct dimensional reference witnesses. Native DOF, nonlinear work, callbacks and allocation budgets are checked.",true,(material:Handle<Material>)->Dynamic,|ctx:&CallCtx<'_>,material:Handle<Material>|{material.with(ctx,|m|{let(a,b)=m.data.lame();fields(&m.refs,vec![a,b],[-1,-2,1,0])})});
object_binding!(hyper_piola,"num_grid_hyper_piola","crate::sim::hyperelastic::first_piola","Native elasticity call in an explicit mechanical reference frame. Numeric displacements/coordinates use Length, body forces use Mass/(Length^2 Time^2), stress/Young modulus use Mass/(Length Time^2), and density uses Mass/Length^3. Reduced 2D hyperelastic energies are per unit out-of-plane thickness. Outputs retain distinct dimensional reference witnesses. Native DOF, nonlinear work, callbacks and allocation budgets are checked.",true,(material:Handle<Material>,gradient:Vec<Vec<f64>>)->Dynamic,|ctx:&CallCtx<'_>,material:Handle<Material>,gradient:Vec<Vec<f64>>|{let gradient=mat3(gradient)?;material.with(ctx,|m|{let p=hyp::first_piola(&m.data,&gradient).map_err(err)?;Ok(Dynamic(api::record(vec![("shape",vec![3usize,3].encode()?),("entries",fields(&m.refs,p.into_iter().flatten().collect(),[-1,-2,1,0])?.0)])))})});
object_binding!(hyper_tangent,"num_grid_hyper_tangent","crate::sim::hyperelastic::tangent_moduli","Native elasticity call in an explicit mechanical reference frame. Numeric displacements/coordinates use Length, body forces use Mass/(Length^2 Time^2), stress/Young modulus use Mass/(Length Time^2), and density uses Mass/Length^3. Reduced 2D hyperelastic energies are per unit out-of-plane thickness. Outputs retain distinct dimensional reference witnesses. Native DOF, nonlinear work, callbacks and allocation budgets are checked.",true,(material:Handle<Material>,gradient:Vec<Vec<f64>>)->Dynamic,|ctx:&CallCtx<'_>,material:Handle<Material>,gradient:Vec<Vec<f64>>|{let gradient=mat3(gradient)?;material.with(ctx,|m|{let p=hyp::tangent_moduli(&m.data,&gradient).map_err(err)?;Ok(Dynamic(api::record(vec![("shape",vec![3usize,3,3,3].encode()?),("entries",fields(&m.refs,p.into_iter().flatten().flatten().flatten().collect(),[-1,-2,1,0])?.0)])))})});
object_binding!(hyper_contact,"num_grid_hyper_contact","","Native elasticity call in an explicit mechanical reference frame. Numeric displacements/coordinates use Length, body forces use Mass/(Length^2 Time^2), stress/Young modulus use Mass/(Length Time^2), and density uses Mass/Length^3. Reduced 2D hyperelastic energies are per unit out-of-plane thickness. Outputs retain distinct dimensional reference witnesses. Native DOF, nonlinear work, callbacks and allocation budgets are checked.",true,(refs:Handle<References>,marker:i32,gap:Input,normal:Many<Input>,penalty:f64,budget:usize)->NewObject<Contact>,|ctx:&CallCtx<'_>,refs:Handle<References>,marker:i32,gap:Input,normal:Many<Input>,penalty:f64,budget:usize|{refs.with(ctx,|r|{if normal.0.len()!=3{return Err(err("three normal components required"));}let inputs=vec![r.dims[0].clone();3];let gap=checked_snapshot(ctx,&gap,budget,&inputs,&r.dims[0])?;let normal=normal.0.iter().map(|p|checked_snapshot(ctx,p,budget,&inputs,&DimExpr::Dimless)).collect::<Result<Vec<_>,_>>()?.try_into().map_err(|_|err("three normal components required"))?;Ok(NewObject::new(Contact{marker,gap,normal,penalty,refs:r.clone()},r.dims.to_vec()))})});
object_binding!(hyper_problem,"num_grid_hyper_problem","","Native elasticity call in an explicit mechanical reference frame. Numeric displacements/coordinates use Length, body forces use Mass/(Length^2 Time^2), stress/Young modulus use Mass/(Length Time^2), and density uses Mass/Length^3. Reduced 2D hyperelastic energies are per unit out-of-plane thickness. Outputs retain distinct dimensional reference witnesses. Native DOF, nonlinear work, callbacks and allocation budgets are checked.",true,(mesh:Handle<Mesh>,material:Handle<Material>,bcs:Handle<Bcs>,body:Many<Input>,contacts:Many<Handle<Contact>>,budget:usize)->NewObject<Hyper>,|ctx:&CallCtx<'_>,mesh:Handle<Mesh>,material:Handle<Material>,bcs:Handle<Bcs>,body:Many<Input>,contacts:Many<Handle<Contact>>,budget:usize|{mesh.with(ctx,|mesh|material.with(ctx,|material|bcs.with(ctx,|bcs|{if contacts.0.len()>16{return Err(err("contact budget exceeded"));}check_refs(&material.refs,&bcs.refs)?;let force=force(ctx,&body,&material.refs,budget)?;let contacts=contacts.0.into_iter().map(|c|c.with(ctx,|c|{check_refs(&material.refs,&c.refs)?;Ok(c.clone())})).collect::<Result<Vec<_>,_>>()?;Ok(NewObject::new(Hyper{mesh:mesh.clone(),material:material.clone(),bcs:bcs.clone(),force,contacts},material.refs.dims.to_vec()))})))});
object_binding!(hyper_static,"num_grid_hyper_static","crate::sim::hyperelastic::solve_static","Native elasticity call in an explicit mechanical reference frame. Numeric displacements/coordinates use Length, body forces use Mass/(Length^2 Time^2), stress/Young modulus use Mass/(Length Time^2), and density uses Mass/Length^3. Reduced 2D hyperelastic energies are per unit out-of-plane thickness. Outputs retain distinct dimensional reference witnesses. Native DOF, nonlinear work, callbacks and allocation budgets are checked.",true,(problem:Handle<Hyper>,steps:usize,tolerance:f64,iterations:usize)->Dynamic,|ctx:&CallCtx<'_>,problem:Handle<Hyper>,steps:usize,tolerance:f64,iterations:usize|{if steps==0||steps>64||iterations==0||iterations>128||steps*iterations>4096||tolerance<=0.0||tolerance>1.0{return Err(err("hyperelastic Newton/load budget invalid"));}problem.with(ctx,|p|{if steps.checked_mul(iterations).and_then(|n|n.checked_mul(p.mesh.data.cells.len())).and_then(|n|n.checked_mul(p.mesh.data.dim.pow(4))).is_none_or(|n|n>50000000){return Err(err("hyperelastic static assembly work cap exceeded"));}p.with(|native|{let s=hyp::solve_static(native,steps,tolerance,iterations,structured::LinearSolver::Auto).map_err(err)?;let dim=p.mesh.data.dim as i64;Ok(Dynamic(api::record(vec![("displacement",fields(&p.material.refs,s.u,[1,0,0,0])?.0),("residuals",fields(&p.material.refs,s.residuals,[dim-2,-2,1,0])?.0),("iterations",s.iterations.encode()?)])))})})});

object_binding!(hyper_rest,"num_grid_hyper_rest","crate::sim::hyperelastic::rest_state","Native elasticity call in an explicit mechanical reference frame. Numeric displacements/coordinates use Length, body forces use Mass/(Length^2 Time^2), stress/Young modulus use Mass/(Length Time^2), and density uses Mass/Length^3. Reduced 2D hyperelastic energies are per unit out-of-plane thickness. Outputs retain distinct dimensional reference witnesses. Native DOF, nonlinear work, callbacks and allocation budgets are checked.",true,(problem:Handle<Hyper>,displacement:Vec<f64>)->NewObject<State>,|ctx:&CallCtx<'_>,problem:Handle<Hyper>,displacement:Vec<f64>|{problem.with(ctx,|p|{if displacement.len()!=p.mesh.data.nodes.len()*p.mesh.data.dim{return Err(err("initial displacement shape mismatch"));}p.with(|native|Ok(NewObject::new(State{data:hyp::rest_state(native,&displacement,structured::LinearSolver::Auto).map_err(err)?,problem:p.clone()},p.material.refs.dims.to_vec())))})});
object_binding!(hyper_kinetic,"num_grid_hyper_kinetic","crate::sim::hyperelastic::kinetic_energy","Native elasticity call in an explicit mechanical reference frame. Numeric displacements/coordinates use Length, body forces use Mass/(Length^2 Time^2), stress/Young modulus use Mass/(Length Time^2), and density uses Mass/Length^3. Reduced 2D hyperelastic energies are per unit out-of-plane thickness. Outputs retain distinct dimensional reference witnesses. Native DOF, nonlinear work, callbacks and allocation budgets are checked.",true,(problem:Handle<Hyper>,velocity:Vec<f64>)->Dynamic,|ctx:&CallCtx<'_>,problem:Handle<Hyper>,velocity:Vec<f64>|{problem.with(ctx,|p|{if velocity.len()!=p.mesh.data.nodes.len()*p.mesh.data.dim{return Err(err("hyperelastic field length mismatch"));}p.with(|native|{let value=hyp::kinetic_energy(native,&velocity);let(s,d)=p.material.refs.derived([p.mesh.data.dim as i64-1,-2,1,0])?;Ok(Dynamic(api::quantity(value*s,d)?))})})});
object_binding!(hyper_strain,"num_grid_hyper_strain","crate::sim::hyperelastic::strain_energy","Native elasticity call in an explicit mechanical reference frame. Numeric displacements/coordinates use Length, body forces use Mass/(Length^2 Time^2), stress/Young modulus use Mass/(Length Time^2), and density uses Mass/Length^3. Reduced 2D hyperelastic energies are per unit out-of-plane thickness. Outputs retain distinct dimensional reference witnesses. Native DOF, nonlinear work, callbacks and allocation budgets are checked.",true,(problem:Handle<Hyper>,displacement:Vec<f64>)->Dynamic,|ctx:&CallCtx<'_>,problem:Handle<Hyper>,displacement:Vec<f64>|{problem.with(ctx,|p|{if displacement.len()!=p.mesh.data.nodes.len()*p.mesh.data.dim{return Err(err("hyperelastic field length mismatch"));}p.with(|native|{let value=hyp::strain_energy(native,&displacement).map_err(err)?;let(s,d)=p.material.refs.derived([p.mesh.data.dim as i64-1,-2,1,0])?;Ok(Dynamic(api::quantity(value*s,d)?))})})});
object_binding!(hyper_dynamic,"num_grid_hyper_dynamic","crate::sim::hyperelastic::solve_dynamic","Native elasticity call in an explicit mechanical reference frame. Numeric displacements/coordinates use Length, body forces use Mass/(Length^2 Time^2), stress/Young modulus use Mass/(Length Time^2), and density uses Mass/Length^3. Reduced 2D hyperelastic energies are per unit out-of-plane thickness. Outputs retain distinct dimensional reference witnesses. Native DOF, nonlinear work, callbacks and allocation budgets are checked.",true,(state:Handle<State>,scheme:String,beta:f64,gamma:f64,rho_infinity:f64,dt:f64,steps:usize,tolerance:f64,observer:Observer,observer_budget:usize)->NewObject<State>,|ctx:&CallCtx<'_>,state:Handle<State>,scheme:String,beta:f64,gamma:f64,rho_infinity:f64,dt:f64,steps:usize,tolerance:f64,observer:Observer,observer_budget:usize|{if observer_budget<steps||observer_budget>100000||dt<=0.0||steps==0||steps>1024||tolerance<=0.0||tolerance>1.0{return Err(err("hyperelastic dynamic time/work controls invalid"));}let scheme=match scheme.as_str(){"Newmark" if beta>0.0&&gamma>=0.5=>hyp::TimeScheme::Newmark{beta,gamma},"GeneralizedAlpha" if (0.0..=1.0).contains(&rho_infinity)=>hyp::TimeScheme::GeneralizedAlpha{rho_inf:rho_infinity},_=>return Err(err("invalid hyperelastic integration scheme/parameters"))};state.with(ctx,|s|{let n=s.data.u.len();if n.checked_mul(n).and_then(|n|n.checked_mul(steps)).and_then(|n|n.checked_mul(30)).is_none_or(|n|n>50000000){return Err(err("dynamic work budget exceeded"));}s.problem.with(|native|{let failure=std::cell::RefCell::new(None);let mut remaining=observer_budget;let mut notify=|state:&hyp::DynamicState|{let result=(||{if remaining==0{return Err(err("dynamic observer budget exhausted"));}remaining-=1;let value=state_value(&s.problem.material.refs,state)?;match ctx.invoke_callback(observer.0,&[value])?{Value::Bool(true)=>Ok(()),Value::Bool(false)=>Err(err("dynamic observer requested cancellation")),_=>Err(err("dynamic observer must return Bool acknowledgement"))}})();if let Err(error)=result{*failure.borrow_mut()=Some(error);std::panic::resume_unwind(Box::new(ObserverStop));}};let result=std::panic::catch_unwind(std::panic::AssertUnwindSafe(||hyp::solve_dynamic(native,scheme,&s.data,dt,steps,tolerance,structured::LinearSolver::Auto,&mut notify)));if let Some(error)=failure.into_inner(){return Err(error);}let data=match result{Ok(result)=>result.map_err(err)?,Err(payload)=>std::panic::resume_unwind(payload)};Ok(NewObject::new(State{data,problem:s.problem.clone()},s.problem.material.refs.dims.to_vec()))})})});
object_binding!(hyper_state_fields,"num_grid_hyper_state_fields","","Native elasticity call in an explicit mechanical reference frame. Numeric displacements/coordinates use Length, body forces use Mass/(Length^2 Time^2), stress/Young modulus use Mass/(Length Time^2), and density uses Mass/Length^3. Reduced 2D hyperelastic energies are per unit out-of-plane thickness. Outputs retain distinct dimensional reference witnesses. Native DOF, nonlinear work, callbacks and allocation budgets are checked.",true,(state:Handle<State>)->Dynamic,|ctx:&CallCtx<'_>,state:Handle<State>|{state.with(ctx,|s|{let f=&s.problem.material.refs;Ok(Dynamic(api::record(vec![("displacement",fields(f,s.data.u.clone(),[1,0,0,0])?.0),("velocity",fields(f,s.data.v.clone(),[1,-1,0,0])?.0),("acceleration",fields(f,s.data.a.clone(),[1,-2,0,0])?.0),("time",api::quantity(s.data.t*f.scales[1],f.dims[1].clone())?)])))})});
struct Observer(crate::kernel::value::CallbackRef);
struct ObserverStop;
impl NativeArg for Observer {
    fn ty() -> TypeExpr {
        let l = DimExpr::Base("Length".into());
        let t = DimExpr::Base("Time".into());
        let derivative = |p: i64| {
            DimExpr::Mul(vec![
                l.clone(),
                DimExpr::Pow(Box::new(t.clone()), Exp::Rat(Rational::new(p, 1))),
            ])
        };
        TypeExpr::Function(
            vec![TypeExpr::Record(vec![
                (
                    "displacement".into(),
                    api::tensor_type(1, TypeExpr::Float.with_dim(l.clone())),
                ),
                (
                    "velocity".into(),
                    api::tensor_type(1, TypeExpr::Float.with_dim(derivative(-1))),
                ),
                (
                    "acceleration".into(),
                    api::tensor_type(1, TypeExpr::Float.with_dim(derivative(-2))),
                ),
                ("time".into(), TypeExpr::Float.with_dim(t)),
            ])],
            Box::new(TypeExpr::Bool),
        )
    }

    fn decode(
        value: &Value,
        _: usize,
        ctx: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let Value::Function(h) = value else {
            return Err(err("typed dynamic observer callback required"));
        };
        let actual = ctx.callback_signature(*h)?;
        if actual.params.len() != 1 || actual.ret != TypeExpr::Bool {
            return Err(err("dynamic observer signature mismatch"));
        }
        let TypeExpr::Function(parameters, _) = Self::ty() else {
            unreachable!()
        };
        if !crate::kernel::overload::matches_type(&parameters[0], &actual.params[0]).accepted() {
            return Err(err("dynamic observer requires measured mechanical state"));
        }
        Ok(Self(*h))
    }
}
fn measured_vector(
    f: &References,
    v: &[f64],
    p: [i64; 4],
) -> Result<Value, OpError> {
    let (s, d) = f.derived(p)?;
    Ok(Value::Tensor(TensorValue {
        shape: vec![v.len()],
        elem: TypeExpr::Float.with_dim(d),
        data: v
            .iter()
            .map(|v| crate::rules::numerical_abi::output_scalar(v * s))
            .collect::<Result<_, _>>()?,
    }))
}
fn state_value(
    f: &References,
    s: &hyp::DynamicState,
) -> Result<Value, OpError> {
    Ok(api::record(vec![
        ("displacement", measured_vector(f, &s.u, [1, 0, 0, 0])?),
        ("velocity", measured_vector(f, &s.v, [1, -1, 0, 0])?),
        ("acceleration", measured_vector(f, &s.a, [1, -2, 0, 0])?),
        ("time", api::quantity(s.t * f.scales[1], f.dims[1].clone())?),
    ]))
}
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &elastic_bcs::BINDING,
            &hyper_material::BINDING,
            &elastic_constitutive::BINDING,
            &elastic_solve::BINDING,
            &elastic_solve_heterogeneous::BINDING,
            &elastic_displacement::BINDING,
            &elastic_eval::BINDING,
            &elastic_stress::BINDING,
            &elastic_von_mises::BINDING,
            &hyper_lame::BINDING,
            &hyper_piola::BINDING,
            &hyper_tangent::BINDING,
            &hyper_contact::BINDING,
            &hyper_problem::BINDING,
            &hyper_static::BINDING,
            &hyper_rest::BINDING,
            &hyper_kinetic::BINDING,
            &hyper_strain::BINDING,
            &hyper_dynamic::BINDING,
            &hyper_state_fields::BINDING,
        ],
        "numerical_vector",
    )
}
