//! Native FVM meshes with explicit sampled-field/boundary-mask generic adapters.
use super::Arg;
use super::DimExpr;
use super::Dynamic;
use super::Frame;
use super::Handle;
use super::Input;
use super::Installer;
use super::LIMIT;
use super::NativeType;
use super::NewObject;
use super::OpError;
use super::RuleError;
use super::Value;
use super::api;
use super::array;
use super::canonical;
use super::compiled_snapshot;
use super::dim;
use super::err;
use super::object_binding;
use super::quantity;
use super::vector;
use crate::sim::physics_fvm as fv;
enum Body {
    D1(fv::Mesh),
    D2(fv::Mesh2D),
    D3(fv::Mesh3D),
}
struct Mesh {
    native: Body,
    frame: Frame,
}
impl NativeType for Mesh {
    const KIND: &'static str = "rssn/native/legacy-fvm-mesh/v1";
    const MAX_BYTES: usize = 1_048_576;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        let (s, h) = self.shape();
        let n = extent(&s)?;
        let values = self.cells();
        if values.len() != n
            || values.iter().any(|v| !v.value.is_finite())
            || h.iter().any(|h| !h.is_finite() || *h <= 0.0)
        {
            return Err(err("FVM grid data/shape/spacing invalid"));
        }
        Ok(())
    }
}
impl Mesh {
    fn shape(&self) -> (Vec<usize>, Vec<f64>) {
        match &self.native {
            | Body::D1(m) => (vec![m.num_cells()], vec![m.dx]),
            | Body::D2(m) => (vec![m.width, m.height], vec![m.dx, m.dy]),
            | Body::D3(m) => (vec![m.width, m.height, m.depth], vec![m.dx, m.dy, m.dz]),
        }
    }

    fn cells(&self) -> &[fv::Cell] {
        match &self.native {
            | Body::D1(m) => &m.cells,
            | Body::D2(m) => &m.cells,
            | Body::D3(m) => &m.cells,
        }
    }

    fn cells_mut(&mut self) -> &mut [fv::Cell] {
        match &mut self.native {
            | Body::D1(m) => &mut m.cells,
            | Body::D2(m) => &mut m.cells,
            | Body::D3(m) => &mut m.cells,
        }
    }

    fn restore(
        &mut self,
        data: Vec<f64>,
    ) -> Result<Dynamic, OpError> {
        if data.len() != self.cells().len() || data.iter().any(|x| !x.is_finite()) {
            return Err(err("native FVM returned malformed/nonfinite field"));
        }
        let values = data
            .iter()
            .copied()
            .map(|v| api::scale_product(v, &[self.frame.value], &[]))
            .collect::<Result<Vec<_>, _>>()?;
        let output = array(values, vec![data.len()], self.frame.dims[2].clone())?;
        // Encoding can fail even after the normalized solver succeeds.
        for (cell, value) in self.cells_mut().iter_mut().zip(data) {
            cell.value = value;
        }
        Ok(Dynamic(output))
    }

    fn controls(
        &self,
        velocity: &Value,
        dt: &Value,
        steps: usize,
    ) -> Result<(Vec<f64>, f64), OpError> {
        let (s, h) = self.shape();
        let v = vector(velocity, &dim([1, -1, 0, 0])?)?
            .into_iter()
            .map(|v| api::scale_product(v, &[self.frame.time], &[self.frame.length]))
            .collect::<Result<Vec<_>, _>>()?;
        let dt = quantity(dt, &self.frame.dims[1])? / self.frame.time;
        if v.len() != s.len()
            || !dt.is_finite()
            || dt <= 0.0
            || steps > 65536
            || steps
                .checked_mul(self.cells().len())
                .is_none_or(|n| n > 20_000_000)
            || dt * v.iter().zip(&h).map(|(v, h)| v.abs() / h).sum::<f64>() > 1.0
        {
            return Err(err("FVM advection dimensions/CFL/work budget invalid"));
        }
        Ok((v, dt))
    }
}
fn extent(s: &[usize]) -> Result<usize, OpError> {
    if !(1..=3).contains(&s.len()) || s.iter().any(|n| *n < 3) {
        return Err(err("FVM sampled grids need1..3 axes with at least3cells"));
    }
    s.iter()
        .try_fold(1usize, |n, d| n.checked_mul(*d))
        .filter(|n| *n <= LIMIT)
        .ok_or_else(|| err("FVM mesh allocation budget exceeded"))
}
fn init(
    f: &Frame,
    shape: &[usize],
    length: &Value,
    values: &Value,
) -> Result<(Vec<f64>, Vec<f64>), OpError> {
    let n = extent(shape)?;
    let lengths = vector(length, &f.dims[0])?;
    if lengths.len() != shape.len() {
        return Err(err("one physical domain length per mesh axis is required"));
    }
    let h = lengths
        .into_iter()
        .zip(shape)
        .map(|(v, n)| api::scale_product(v, &[], &[f.length, *n as f64]))
        .collect::<Result<Vec<_>, _>>()?;
    let data = vector(values, &f.dims[2])?
        .into_iter()
        .map(|v| api::scale_product(v, &[], &[f.value]))
        .collect::<Result<Vec<_>, _>>()?;
    if h.len() != shape.len() || h.iter().any(|h| *h <= 0.0) || data.len() != n {
        return Err(err("sampled mesh initializer length/shape mismatch"));
    }
    Ok((h, data))
}
fn mask(
    v: &Value,
    n: &[usize],
) -> Result<Vec<bool>, OpError> {
    let Value::Tensor(t) = v else {
        return Err(err("boundary mask must be a Boolean tensor"));
    };
    let cells = extent(n)?;
    if t.data.len() != cells {
        return Err(err("boundary mask length mismatch"));
    }
    let is_3d = if let [width, height, depth] = n {
        if t.shape != [cells] && t.shape != [*depth, *height, *width] {
            return Err(err(
                "3D boundary mask shape must be flat or [depth, height, width]",
            ));
        }
        true
    } else {
        false
    };
    api::normalized(v)?;
    let mask = t
        .data
        .iter()
        .map(|v| <bool as Arg>::decode(v, 0))
        .collect::<Result<Vec<_>, _>>()?;
    if is_3d {
        return Ok(mask);
    }
    for (i, selected) in mask.iter().enumerate() {
        let mut index = i;
        for len in n {
            let c = index % len;
            index /= len;
            if (c == 0 || c + 1 == *len) && !*selected {
                return Err(err(
                    "native boundary predicate must cover every exterior cell",
                ));
            }
        }
    }
    Ok(mask)
}
fn water(
    h: &[f64],
    momenta: &[Vec<f64>],
    shape: &[usize],
    spacing: &[f64],
    dt: f64,
    steps: usize,
    g: f64,
) -> Result<(), OpError> {
    if extent(shape)? != h.len()
        || momenta.len() != shape.len()
        || momenta.iter().any(|p| p.len() != h.len())
        || spacing.len() != shape.len()
        || spacing.iter().any(|h| *h <= 0.0)
        || h.iter().any(|h| *h < 0.0)
        || g <= 0.0
        || !dt.is_finite()
        || dt <= 0.0
        || steps > 65536
        || steps.checked_mul(h.len()).is_none_or(|n| n > 20_000_000)
    {
        return Err(err(
            "shallow-water dimensions/physical controls/work budget invalid",
        ));
    }
    for (i, depth) in h.iter().enumerate() {
        if *depth == 0.0 && momenta.iter().any(|p| p[i] != 0.0) {
            return Err(err("dry cells require zero momentum"));
        }
        let c = (g * depth).sqrt();
        let rate = momenta
            .iter()
            .zip(spacing)
            .map(|(p, dx)| {
                (if *depth > 0.0 {
                    (p[i] / depth).abs()
                } else {
                    0.0
                } + c)
                    / dx
            })
            .sum::<f64>();
        if dt * rate > 0.8 {
            return Err(err("shallow-water explicit CFL bound exceeded"));
        }
    }
    Ok(())
}
object_binding!(
    mesh1_new, "num_grid_legacy_fvm_mesh1_new", "crate::sim::physics_fvm::Mesh::new",
    "Creates a retained one-dimensional finite-volume mesh directly from sampled cell averages. Supply a one-entry cell-count vector with at least three cells and a one-entry measured Length vector for the total domain extent.\n\nSupply one value in the frame's field dimension per cell in increasing coordinate order. The adapter uses those samples directly and stores a uniform cell spacing from the total extent and cell count.\n\nInvalid shapes, dimensions, samples and nonpositive or unrepresentable spacing produce errors. The mesh handle can be reused for advection and data retrieval.",
    false, (frame:Handle<Frame>,shape:Vec<usize>,lengths:Input,values:Input) -> NewObject<Mesh>,
    |ctx,frame,shape,lengths,values| frame.with(ctx, |f| {
        if shape.len() != 1 {
            return Err(err("mesh rank mismatch"));
        }
        let (spacing, data) = init(f, &shape, &lengths.0, &values.0)?;
        // The values are already sampled by cell. Reconstructing a normalized
        // whole-domain extent can overflow even when each spacing is finite.
        let native = Body::D1(fv::Mesh {
            cells: data.into_iter().map(|value| fv::Cell { value }).collect(),
            dx: spacing[0],
        });
        Ok(NewObject::new(Mesh { native, frame: f.clone() }, f.dims.to_vec()))
    })
);
object_binding!(mesh2_new,"num_grid_legacy_fvm_mesh2_new","crate::sim::physics_fvm::Mesh2D::new","Actual native finite-volume specialization. Measured field, velocity, length and time retain explicit reference units. Sampled initializer arrays and finite boundary masks instantiate generic native functions; masks must cover exterior cells. Returned values update the retained mesh. Native CFL, work and allocation bounds are checked.",false,(frame:Handle<Frame>,shape:Vec<usize>,lengths:Input,values:Input)->NewObject<Mesh>,|ctx,frame,shape,lengths,values|{let _=ctx;frame.with(ctx,|f|{if shape.len()!=2{return Err(err("mesh rank mismatch"));}let(h,data)=init(f,&shape,&lengths.0,&values.0)?;let native=Body::D2(fv::Mesh2D::new(shape[0],shape[1],(h[0]*shape[0] as f64,h[1]*shape[1] as f64),|x,y|data[((x/h[0]).floor() as usize).min(shape[0]-1)+shape[0]*((y/h[1]).floor() as usize).min(shape[1]-1)]));Ok(NewObject::new(Mesh{native,frame:f.clone()},f.dims.to_vec()))})});
object_binding!(mesh3_new,"num_grid_legacy_fvm_mesh3_new","crate::sim::physics_fvm::Mesh3D::new","Actual native finite-volume specialization. Measured field, velocity, length and time retain explicit reference units. Sampled initializer arrays and finite boundary masks instantiate generic native functions; masks must cover exterior cells. Returned values update the retained mesh. Native CFL, work and allocation bounds are checked.",false,(frame:Handle<Frame>,shape:Vec<usize>,lengths:Input,values:Input)->NewObject<Mesh>,|ctx,frame,shape,lengths,values|{let _=ctx;frame.with(ctx,|f|{if shape.len()!=3{return Err(err("mesh rank mismatch"));}let(h,data)=init(f,&shape,&lengths.0,&values.0)?;let native=Body::D3(fv::Mesh3D::new(shape[0],shape[1],shape[2],(h[0]*shape[0] as f64,h[1]*shape[1] as f64,h[2]*shape[2] as f64),|x,y,z|data[((x/h[0]).floor() as usize).min(shape[0]-1)+shape[0]*(((y/h[1]).floor() as usize).min(shape[1]-1)+shape[1]*((z/h[2]).floor() as usize).min(shape[2]-1))]));Ok(NewObject::new(Mesh{native,frame:f.clone()},f.dims.to_vec()))})});
object_binding!(num_cells,"num_grid_legacy_fvm_num_cells","crate::sim::physics_fvm::Mesh::num_cells","Actual native finite-volume specialization. Measured field, velocity, length and time retain explicit reference units. Sampled initializer arrays and finite boundary masks instantiate generic native functions; masks must cover exterior cells. Returned values update the retained mesh. Native CFL, work and allocation bounds are checked.",false,(a:Handle<Mesh>)->usize,|ctx,a|{let _=ctx;a.with(ctx,|a|{let Body::D1(m)=&a.native else{return Err(err("num_cells source method is specialized to the native1D mesh"));};Ok(m.num_cells())})});
object_binding!(mesh_data,"num_grid_legacy_fvm_mesh_data","","Actual native finite-volume specialization. Measured field, velocity, length and time retain explicit reference units. Sampled initializer arrays and finite boundary masks instantiate generic native functions; masks must cover exterior cells. Returned values update the retained mesh. Native CFL, work and allocation bounds are checked.",false,(a:Handle<Mesh>)->Dynamic,|ctx,a|{let _=ctx;a.with(ctx,|a|{let values=a.cells().iter().map(|v|api::scale_product(v.value,&[a.frame.value],&[])).collect::<Result<Vec<_>,_>>()?;Ok(Dynamic(array(values,vec![a.cells().len()],a.frame.dims[2].clone())?))})});
object_binding!(
    advection1, "num_grid_legacy_fvm_advection1", "crate::sim::physics_fvm::solve_advection_1d",
    "Advances a retained 1D finite-volume mesh with an upwind flux. Supply a one-entry measured Length/Time velocity, a positive measured Time step, an integer step count, and the left and right ghost-cell values in the mesh field dimension. Both ghost values are validated.\n\nCFL must be at most one. The upwind face uses the ghost value on the side from which the flow enters. Step count and cell-step work are bounded.\n\nThe measured field is returned in increasing cell order and replaces the retained field on success. An invalid input or unrepresentable output leaves the previous field intact.",
    false, (a:Handle<Mesh>,velocity:Input,dt:Input,steps:usize,left:Input,right:Input) -> Dynamic,
    |ctx,a,velocity,dt,steps,left,right| a.with_mut(ctx, |a| {
        let (velocity, dt) = a.controls(&velocity.0, &dt.0, steps)?;
        let boundary = |value: &Value| {
            api::scale_without_underflow(quantity(value, &a.frame.dims[2])?, &[], &[a.frame.value])
        };
        let left = boundary(&left.0)?;
        let right = boundary(&right.0)?;
        let Body::D1(mesh) = &mut a.native else {
            return Err(err("expected 1D mesh"));
        };
        let output = fv::solve_advection_1d(mesh, velocity[0], dt, steps, || (left, right));
        a.restore(output)
    })
);
object_binding!(advection2,"num_grid_legacy_fvm_advection2","crate::sim::physics_fvm::solve_advection_2d","Actual native finite-volume specialization. Measured field, velocity, length and time retain explicit reference units. Sampled initializer arrays and finite boundary masks instantiate generic native functions; masks must cover exterior cells. Returned values update the retained mesh. Native CFL, work and allocation bounds are checked.",false,(a:Handle<Mesh>,velocity:Input,dt:Input,steps:usize,boundary_mask:Input)->Dynamic,|ctx,a,velocity,dt,steps,boundary_mask|{let _=ctx;a.with_mut(ctx,|a|{let(v,dt)=a.controls(&velocity.0,&dt.0,steps)?;let(n,_)=a.shape();let mask=mask(&boundary_mask.0,&n)?;let Body::D2(m)=&mut a.native else{return Err(err("mesh rank mismatch"));};let out=fv::solve_advection_2d(m,(v[0],v[1]),dt,steps,|x,y,nx,_ny|mask[x+nx*y]);a.restore(out)})});
object_binding!(advection3,"num_grid_legacy_fvm_advection3","crate::sim::physics_fvm::solve_advection_3d","Advances the retained 3D finite-volume mesh using measured velocity and timestep. The Boolean mask is flat in x-fastest order or shaped [depth, height, width]. True cells are set to zero on each step. Unmarked exterior cells use zero inflow and free outflow. Returns a flat measured field in x-fastest order. CFL and work bounds are checked.",false,(a:Handle<Mesh>,velocity:Input,dt:Input,steps:usize,boundary_mask:Input)->Dynamic,|ctx,a,velocity,dt,steps,boundary_mask|{let _=ctx;a.with_mut(ctx,|a|{let(v,dt)=a.controls(&velocity.0,&dt.0,steps)?;let(n,_)=a.shape();let mask=mask(&boundary_mask.0,&n)?;let Body::D3(m)=&mut a.native else{return Err(err("mesh rank mismatch"));};let out=fv::solve_advection_3d(m,(v[0],v[1],v[2]),dt,steps,|x,y,z,nx,ny,_nz|mask[x+nx*(y+ny*z)]);a.restore(out)})});
object_binding!(burgers1,"num_grid_legacy_fvm_burgers1","crate::sim::physics_fvm::solve_burgers_1d","Actual native finite-volume specialization. Measured field, velocity, length and time retain explicit reference units. Sampled initializer arrays and finite boundary masks instantiate generic native functions; masks must cover exterior cells. Returned values update the retained mesh. Native CFL, work and allocation bounds are checked.",false,(a:Handle<Mesh>,dt:Input,steps:usize)->Dynamic,|ctx,a,dt,steps|{let _=ctx;a.with_mut(ctx,|a|{api::same_dimension(&a.frame.dims[2],&dim([1,-1,0,0])?)?;let dt=quantity(&dt.0,&dim([0,1,0,0])?)?;let Body::D1(m)=&a.native else{return Err(err("Burgers1 requires1D mesh"));};let mut physical=fv::Mesh{cells:m.cells.iter().map(|c|Ok(fv::Cell{value:api::scale_product(c.value,&[a.frame.value],&[])?})).collect::<Result<_,OpError>>()?,dx:api::scale_product(m.dx,&[a.frame.length],&[])?};let max=physical.cells.iter().map(|v|v.value.abs()).fold(0.0,f64::max);if dt<=0.0||steps>65536||steps.checked_mul(m.cells.len()).is_none_or(|n|n>20_000_000)||dt*max/physical.dx>0.8{return Err(err("Burgers timestep/CFL/work budget invalid"));}let out=fv::solve_burgers_1d(&mut physical,dt,steps).into_iter().map(|v|api::scale_product(v,&[],&[a.frame.value])).collect::<Result<Vec<_>,_>>()?;a.restore(out)})});
object_binding!(lax_friedrichs_flux,"num_grid_legacy_fvm_lax_friedrichs_flux","crate::sim::physics_fvm::lax_friedrichs_flux","Actual native finite-volume specialization. Measured field, velocity, length and time retain explicit reference units. Sampled initializer arrays and finite boundary masks instantiate generic native functions; masks must cover exterior cells. Returned values update the retained mesh. Native CFL, work and allocation bounds are checked.",false,(left:Input,right:Input,dt:Input,dx:Input,flux_program:Input)->Dynamic,|ctx,left,right,dt,dx,flux_program|{let _=ctx;{let d=api::dimension(&left.0)?;let l=quantity(&left.0,&d)?;let r=quantity(&right.0,&d)?;let dt=quantity(&dt.0,&dim([0,1,0,0])?)?;let dx=quantity(&dx.0,&dim([1,0,0,0])?)?;if dt<=0.0||dx<=0.0{return Err(err("flux spacing/step must be positive"));}let f=compiled_snapshot(ctx,&flux_program.0,2)?;if f.input_dimensions.len()!=1{return Err(err("flux program must have1state input"));}api::same_dimension(&f.input_dimensions[0],&d)?;let outdim=canonical(DimExpr::Mul(vec![d,dim([1,-1,0,0])?]))?;api::same_dimension(&f.output_dimension,&outdim)?;let result=fv::lax_friedrichs_flux(l,r,dt,dx,|v|f.sample(&[v]));f.finish()?;Ok(Dynamic(api::quantity(result,outdim)?))}});
object_binding!(shallow_water1,"num_grid_legacy_fvm_shallow_water1","crate::sim::physics_fvm::solve_shallow_water_1d","Actual native finite-volume specialization. Measured field, velocity, length and time retain explicit reference units. Sampled initializer arrays and finite boundary masks instantiate generic native functions; masks must cover exterior cells. Returned values update the retained mesh. Native CFL, work and allocation bounds are checked.",false,(depth:Input,momentum_x:Input,shape:Vec<usize>,spacing:Input,dt:Input,steps:usize,gravity:Input)->Dynamic,|ctx,depth,momentum_x,shape,spacing,dt,steps,gravity|{let _=ctx;{let h=vector(&depth.0,&dim([1,0,0,0])?)?;let momentum_x=vector(&momentum_x.0,&dim([2,-1,0,0])?)?;let spacing=vector(&spacing.0,&dim([1,0,0,0])?)?;let dt=quantity(&dt.0,&dim([0,1,0,0])?)?;let g=quantity(&gravity.0,&dim([1,-2,0,0])?)?;if shape.len()!=1{return Err(err("shallow-water rank mismatch"));}water(&h,std::slice::from_ref(&momentum_x),&shape,&spacing,dt,steps,g)?;let out=fv::solve_shallow_water_1d(h,momentum_x,spacing[0],dt,steps,g);let n=out.len();Ok(Dynamic(api::record(vec![("depth",array(out.iter().map(|v|v.h).collect(),vec![n],dim([1,0,0,0])?)?),("momentum_x",array(out.iter().map(|v|v.hu).collect(),vec![n],dim([2,-1,0,0])?)?)]))) }});
object_binding!(shallow_water2,"num_grid_legacy_fvm_shallow_water2","crate::sim::physics_fvm::solve_shallow_water_2d","Actual native finite-volume specialization. Measured field, velocity, length and time retain explicit reference units. Sampled initializer arrays and finite boundary masks instantiate generic native functions; masks must cover exterior cells. Returned values update the retained mesh. Native CFL, work and allocation bounds are checked.",false,(depth:Input,momentum_x:Input,momentum_y:Input,shape:Vec<usize>,spacing:Input,dt:Input,steps:usize,gravity:Input)->Dynamic,|ctx,depth,momentum_x,momentum_y,shape,spacing,dt,steps,gravity|{let _=ctx;{let h=vector(&depth.0,&dim([1,0,0,0])?)?;let momentum_x=vector(&momentum_x.0,&dim([2,-1,0,0])?)?;let momentum_y=vector(&momentum_y.0,&dim([2,-1,0,0])?)?;let spacing=vector(&spacing.0,&dim([1,0,0,0])?)?;let dt=quantity(&dt.0,&dim([0,1,0,0])?)?;let g=quantity(&gravity.0,&dim([1,-2,0,0])?)?;if shape.len()!=2{return Err(err("shallow-water rank mismatch"));}water(&h,&[momentum_x.clone(),momentum_y.clone()],&shape,&spacing,dt,steps,g)?;let out=fv::solve_shallow_water_2d(&h,&momentum_x,&momentum_y,shape.try_into().map_err(|_|err("rank mismatch"))?,spacing.try_into().map_err(|_|err("spacing mismatch"))?,dt,steps,g).map_err(err)?;let n=out.len();Ok(Dynamic(api::record(vec![("depth",array(out.iter().map(|v|v.h).collect(),vec![n],dim([1,0,0,0])?)?),("momentum_x",array(out.iter().map(|v|v.hu).collect(),vec![n],dim([2,-1,0,0])?)?),("momentum_y",array(out.iter().map(|v|v.hv).collect(),vec![n],dim([2,-1,0,0])?)?)]))) }});
object_binding!(shallow_water3,"num_grid_legacy_fvm_shallow_water3","crate::sim::physics_fvm::solve_shallow_water_3d","Actual native finite-volume specialization. Measured field, velocity, length and time retain explicit reference units. Sampled initializer arrays and finite boundary masks instantiate generic native functions; masks must cover exterior cells. Returned values update the retained mesh. Native CFL, work and allocation bounds are checked.",false,(depth:Input,momentum_x:Input,momentum_y:Input,momentum_z:Input,shape:Vec<usize>,spacing:Input,dt:Input,steps:usize,gravity:Input)->Dynamic,|ctx,depth,momentum_x,momentum_y,momentum_z,shape,spacing,dt,steps,gravity|{let _=ctx;{let h=vector(&depth.0,&dim([1,0,0,0])?)?;let momentum_x=vector(&momentum_x.0,&dim([2,-1,0,0])?)?;let momentum_y=vector(&momentum_y.0,&dim([2,-1,0,0])?)?;let momentum_z=vector(&momentum_z.0,&dim([2,-1,0,0])?)?;let spacing=vector(&spacing.0,&dim([1,0,0,0])?)?;let dt=quantity(&dt.0,&dim([0,1,0,0])?)?;let g=quantity(&gravity.0,&dim([1,-2,0,0])?)?;if shape.len()!=3{return Err(err("shallow-water rank mismatch"));}water(&h,&[momentum_x.clone(),momentum_y.clone(),momentum_z.clone()],&shape,&spacing,dt,steps,g)?;let out=fv::solve_shallow_water_3d(&h,&momentum_x,&momentum_y,&momentum_z,shape.try_into().map_err(|_|err("rank mismatch"))?,spacing.try_into().map_err(|_|err("spacing mismatch"))?,dt,steps,g).map_err(err)?;let n=out.len();Ok(Dynamic(api::record(vec![("depth",array(out.iter().map(|v|v.h).collect(),vec![n],dim([1,0,0,0])?)?),("momentum_x",array(out.iter().map(|v|v.hu).collect(),vec![n],dim([2,-1,0,0])?)?),("momentum_y",array(out.iter().map(|v|v.hv).collect(),vec![n],dim([2,-1,0,0])?)?),("momentum_z",array(out.iter().map(|v|v.hw).collect(),vec![n],dim([2,-1,0,0])?)?)]))) }});
object_binding!(scenario_advection1,"num_grid_legacy_fvm_scenario_advection1","crate::sim::physics_fvm::simulate_1d_advection_scenario","Actual native finite-volume specialization. Measured field, velocity, length and time retain explicit reference units. Sampled initializer arrays and finite boundary masks instantiate generic native functions; masks must cover exterior cells. Returned values update the retained mesh. Native CFL, work and allocation bounds are checked.",false,()->Vec<f64>,|ctx,|{let _=ctx;Ok(fv::simulate_1d_advection_scenario())});
object_binding!(scenario_advection2,"num_grid_legacy_fvm_scenario_advection2","crate::sim::physics_fvm::simulate_2d_advection_scenario","Actual native finite-volume specialization. Measured field, velocity, length and time retain explicit reference units. Sampled initializer arrays and finite boundary masks instantiate generic native functions; masks must cover exterior cells. Returned values update the retained mesh. Native CFL, work and allocation bounds are checked.",false,()->Vec<f64>,|ctx,|{let _=ctx;Ok(fv::simulate_2d_advection_scenario())});
object_binding!(scenario_advection3,"num_grid_legacy_fvm_scenario_advection3","crate::sim::physics_fvm::simulate_3d_advection_scenario","Runs the built-in 3D Gaussian advection example on a 30-by-30-by-30 unit cube. The outer layer is set to zero each step. The returned dimensionless field is flattened in x-fastest order.",false,()->Vec<f64>,|ctx,|{let _=ctx;Ok(fv::simulate_3d_advection_scenario())});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &mesh1_new::BINDING,
            &mesh2_new::BINDING,
            &mesh3_new::BINDING,
            &num_cells::BINDING,
            &mesh_data::BINDING,
            &advection1::BINDING,
            &advection2::BINDING,
            &advection3::BINDING,
            &burgers1::BINDING,
            &lax_friedrichs_flux::BINDING,
            &shallow_water1::BINDING,
            &shallow_water2::BINDING,
            &shallow_water3::BINDING,
            &scenario_advection1::BINDING,
            &scenario_advection2::BINDING,
            &scenario_advection3::BINDING,
        ],
        "numerical_vector",
    )
}
