//! Handle-based FFI API for physics sim linear elasticity functions.

use crate::numerical::matrix::Matrix;
use crate::physics::physics_sim::linear_elasticity;

/// Runs the 2D cantilever beam scenario and returns the displacement results as a Matrix handle (Nx2).
#[no_mangle]

pub extern "C" fn rssn_physics_sim_linear_elasticity_simulate_cantilever(
) -> *mut Matrix<f64> {

    // This scenario currently saves to CSV, I'll modify it slightly or use the core function
    // For now, let's just run a basic setup and return the displacements
    let beam_length = 10.0;

    let beam_height = 2.0;

    let nx = 10;

    let ny = 2;

    let mut nodes = Vec::new();

    for j in 0 ..= ny {

        for i in 0 ..= nx {

            nodes.push((
                i as f64 * beam_length
                    / nx as f64,
                j as f64 * beam_height
                    / ny as f64,
            ));
        }
    }

    let mut elements = Vec::new();

    for j in 0 .. ny {

        for i in 0 .. nx {

            let n1 = j * (nx + 1) + i;

            let n2 =
                j * (nx + 1) + i + 1;

            let n3 = (j + 1) * (nx + 1)
                + i
                + 1;

            let n4 =
                (j + 1) * (nx + 1) + i;

            elements
                .push([n1, n2, n3, n4]);
        }
    }

    let fixed_nodes: Vec<usize> = (0
        ..= ny)
        .map(|j| j * (nx + 1))
        .collect();

    let loads = vec![(
        (ny / 2) * (nx + 1) + nx,
        0.0,
        -1e3,
    )];

    let params = linear_elasticity::ElasticityParameters {
        nodes,
        elements,
        youngs_modulus : 1e7,
        poissons_ratio : 0.3,
        fixed_nodes,
        loads,
    };

    match linear_elasticity::run_elasticity_simulation(&params) {
        | Ok(d) => {

            let n = d.len() / 2;

            Box::into_raw(Box::new(
                Matrix::new(n, 2, d),
            ))
        },
        | Err(_) => std::ptr::null_mut(),
    }
}
