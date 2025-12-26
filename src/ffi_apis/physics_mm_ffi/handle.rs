//! Handle-based FFI API for physics MM (Meshless Methods) functions.

use crate::numerical::matrix::Matrix;
use crate::physics::physics_mm;

/// Simulates the dam break scenario and returns the final particle positions as a Matrix handle (Nx2).
#[no_mangle]
pub extern "C" fn rssn_physics_mm_simulate_dam_break() -> *mut Matrix<f64> {
    let results = physics_mm::simulate_dam_break_2d_scenario();
    let rows = results.len();
    if rows == 0 {
        return std::ptr::null_mut();
    }
    let mut flattened = Vec::with_capacity(rows * 2);
    for (x, y) in results {
        flattened.push(x);
        flattened.push(y);
    }
    Box::into_raw(Box::new(Matrix::new(rows, 2, flattened)))
}

/// Creates a new SPH system.
#[no_mangle]
pub extern "C" fn rssn_physics_mm_sph_create(
    h: f64,
    bounds_x: f64,
    bounds_y: f64,
) -> *mut physics_mm::SPHSystem {
    let system = physics_mm::SPHSystem {
        particles: Vec::new(),
        poly6: physics_mm::Poly6Kernel::new(h),
        spiky: physics_mm::SpikyKernel::new(h),
        gravity: physics_mm::Vector2D { x: 0.0, y: -9.8 },
        viscosity: 0.01,
        gas_const: 2000.0,
        rest_density: 1000.0,
        bounds: physics_mm::Vector2D {
            x: bounds_x,
            y: bounds_y,
        },
    };
    Box::into_raw(Box::new(system))
}

/// Frees an SPH system.
#[no_mangle]
pub unsafe extern "C" fn rssn_physics_mm_sph_free(system: *mut physics_mm::SPHSystem) {
    if !system.is_null() {
        let _ = Box::from_raw(system);
    }
}

/// Adds a particle to the SPH system.
#[no_mangle]
pub unsafe extern "C" fn rssn_physics_mm_sph_add_particle(
    system: *mut physics_mm::SPHSystem,
    x: f64,
    y: f64,
    vx: f64,
    vy: f64,
    mass: f64,
) {
    if let Some(sys) = system.as_mut() {
        sys.particles.push(physics_mm::Particle {
            pos: physics_mm::Vector2D { x, y },
            vel: physics_mm::Vector2D { x: vx, y: vy },
            force: physics_mm::Vector2D::default(),
            density: 0.0,
            pressure: 0.0,
            mass,
        });
    }
}

/// Updates the SPH system by one time step.
#[no_mangle]
pub unsafe extern "C" fn rssn_physics_mm_sph_update(system: *mut physics_mm::SPHSystem, dt: f64) {
    if let Some(sys) = system.as_mut() {
        sys.update(dt);
    }
}

/// Returns the number of particles in the SPH system.
#[no_mangle]
pub unsafe extern "C" fn rssn_physics_mm_sph_get_particle_count(
    system: *mut physics_mm::SPHSystem,
) -> usize {
    if let Some(sys) = system.as_ref() {
        sys.particles.len()
    } else {
        0
    }
}

/// Gets particle positions as a Matrix (Nx2).
#[no_mangle]
pub unsafe extern "C" fn rssn_physics_mm_sph_get_positions(
    system: *mut physics_mm::SPHSystem,
) -> *mut Matrix<f64> {
    if let Some(sys) = system.as_ref() {
        let rows = sys.particles.len();
        let mut flattened = Vec::with_capacity(rows * 2);
        for p in &sys.particles {
            flattened.push(p.pos.x);
            flattened.push(p.pos.y);
        }
        Box::into_raw(Box::new(Matrix::new(rows, 2, flattened)))
    } else {
        std::ptr::null_mut()
    }
}
