use rssn::physics::physics_sim::navier_stokes_fluid::{run_channel_flow, NavierStokesParameters};
use ndarray::Array2;

#[cfg(feature = "output")]
use rssn::output::plotting::{plot_heatmap_2d, PlotConfig, plot_vector_field_2d};

fn main() {

    println!(
        "Starting Fluid Dynamics: \
         Karman Vortex Street..."
    );

    // 1. Setup
    // Grid size must be compatible with MG solver (2^k + 1)
    // We'll use 129x129 (k=7)
    const N: usize = 129;

    const RE: f64 = 500.0; // Moderate Reynolds number for instability
    const DT: f64 = 0.005;

    const TOTAL_TIME: f64 = 2.0;

    let iterations =
        (TOTAL_TIME / DT) as usize;

    // 2. Create Obstacle (Cylinder in the middle-left)
    let mut obstacle =
        Array2::from_elem(
            (N, N),
            false,
        );

    let cx = N / 4;

    let cy = N / 2;

    let rad = N / 10;

    println!(
        "Creating cylinder obstacle \
         at ({}, {}) with radius {}",
        cx, cy, rad
    );

    for j in 0 .. N {

        for i in 0 .. N {

            let dx =
                i as i32 - cx as i32;

            let dy =
                j as i32 - cy as i32;

            if dx * dx + dy * dy
                <= (rad * rad) as i32
            {

                obstacle[[j, i]] = true;
            }
        }
    }

    // 3. Run Simulation
    println!(
        "Running simulation for {} \
         iterations (Re={})...",
        iterations, RE
    );

    // Using our newly created solver
    let result = run_channel_flow(
        N,
        N,
        RE,
        DT,
        iterations,
        &obstacle,
    );

    match result {
        | Ok((u, v, p)) => {

            println!(
                "Simulation \
                 successful."
            );

            // 4. Visualization
            #[cfg(feature = "output")]
            {

                println!(
                    "Generating \
                     plots..."
                );

                let mut config =
                    PlotConfig::default(
                    );

                config.width = 1024;

                config.height = 1024;

                // Velocity Magnitude
                let mut speed =
                    Array2::zeros((
                        N, N,
                    ));

                for j in 0 .. N {

                    for i in 0 .. N {

                        let s = (u
                            [[j, i]]
                        .powi(2)
                            + v[[
                                j, i,
                            ]]
                            .powi(2))
                        .sqrt();

                        speed[[
                            N - 1 - j,
                            i,
                        ]] = s; // Flip Y for plotting
                    }
                }

                config.caption = format!("Velocity Magnitude (Re={})", RE);

                if let Err(e) = plot_heatmap_2d(&speed, "karman_velocity_mag.png", Some(config.clone())) {
                     eprintln!("Error plotting velocity: {}", e);
                } else {
                     println!("Saved karman_velocity_mag.png");
                }

                // Vorticity (Curl) = dv/dx - du/dy
                let mut vorticity =
                    Array2::zeros((
                        N, N,
                    ));

                // Compute curl
                for j in 1 .. N - 1 {

                    for i in 1 .. N - 1
                    {

                        let dv_dx =
                            (v[[
                                j,
                                i + 1,
                            ]] - v[[
                                j,
                                i - 1,
                            ]]) / 2.0;

                        let du_dy =
                            (u[[
                                j + 1,
                                i,
                            ]] - u[[
                                j - 1,
                                i,
                            ]]) / 2.0;

                        vorticity[[
                            N - 1 - j,
                            i,
                        ]] = dv_dx
                            - du_dy;
                    }
                }

                config.caption = format!(
                    "Vorticity (Re={})",
                    RE
                );

                if let Err(e) = plot_heatmap_2d(&vorticity, "karman_vorticity.png", Some(config.clone())) {
                     eprintln!("Error plotting vorticity: {}", e);
                } else {
                     println!("Saved karman_vorticity.png");
                }
            }
        },
        | Err(e) => {

            eprintln!(
                "Simulation failed: {}",
                e
            );
        },
    }
}
