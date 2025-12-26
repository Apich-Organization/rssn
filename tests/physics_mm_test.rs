//! Unit and property-based tests for the physics MM (Meshless Methods) module.

use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use rssn::physics::physics_mm::*;

#[test]

fn test_sph_density_pressure_basic() {

    let h = 0.5;

    let mut system = SPHSystem {
        particles: Vec::new(),
        poly6: Poly6Kernel::new(h),
        spiky: SpikyKernel::new(h),
        gravity: Vector2D::new(
            0.0, 0.0,
        ),
        viscosity: 0.1,
        gas_const: 1000.0,
        rest_density: 1000.0,
        bounds: Vector2D::new(
            10.0, 10.0,
        ),
    };

    // Put two particles very close
    system
        .particles
        .push(Particle {
            pos: Vector2D::new(
                0.0, 0.0,
            ),
            vel: Vector2D::default(),
            force: Vector2D::default(),
            density: 0.0,
            pressure: 0.0,
            mass: 1.0,
        });

    system
        .particles
        .push(Particle {
            pos: Vector2D::new(
                0.05, 0.0,
            ),
            vel: Vector2D::default(),
            force: Vector2D::default(),
            density: 0.0,
            pressure: 0.0,
            mass: 1.0,
        });

    system.compute_density_pressure();

    assert!(
        system.particles[0].density
            > 0.0
    );

    assert_eq!(
        system.particles[0].density,
        system.particles[1].density
    );
}

#[test]

fn test_simulate_dam_break_smoke() {

    let results =
        simulate_dam_break_2d_scenario(
        );

    // 20*10 = 200 particles
    assert_eq!(results.len(), 200);
}

proptest! {
    #[test]
    fn prop_sph_update_stability(dt in 0.001..0.01f64) {
        let h = 0.1;
        let mut system = SPHSystem {
            particles: vec![Particle {
                pos: Vector2D::new(0.5, 0.5),
                vel: Vector2D::new(0.1, 0.1),
                force: Vector2D::default(),
                density: 1000.0, // avoid division by zero if not computed
                pressure: 0.0,
                mass: 1.0,
            }],
            poly6: Poly6Kernel::new(h),
            spiky: SpikyKernel::new(h),
            gravity: Vector2D::new(0.0, -9.8),
            viscosity: 0.01,
            gas_const: 2000.0,
            rest_density: 1000.0,
            bounds: Vector2D::new(1.0, 1.0),
        };

        system.update(dt);
        for p in &system.particles {
            prop_assert!(p.pos.x.is_finite());
            prop_assert!(p.pos.y.is_finite());
        }
    }
}
