#[cfg(feature = "output")]
use ndarray::Array2;
#[cfg(feature = "output")]
use rssn::output::plotting::PlotConfig;
#[cfg(feature = "output")]
use rssn::output::plotting::plot_heatmap_2d;
use rssn::physics::physics_fdm::Dimensions;
use rssn::physics::physics_fdm::FdmGrid;
use rssn::physics::physics_fdm::FdmSolverConfig3D;

fn main() {

    println!(
        "Starting 3D Wave Equation \
         Ripple Example with Progress \
         Frames..."
    );

    // 1. Configuration for the simulation
    const WIDTH: usize = 64;

    const HEIGHT: usize = 64;

    const DEPTH: usize = 64;

    const C: f64 = 1.0; // Wave speed

    let config = FdmSolverConfig3D {
        width: WIDTH,
        height: HEIGHT,
        depth: DEPTH,
        dx: 1.0,
        dy: 1.0,
        dz: 1.0,
        dt: 0.15, /* Stability condition: dt < dx / (c * sqrt(3)) */
        steps: 100,
    };

    // 2. Initial state setup
    let dims = Dimensions::D3(
        WIDTH,
        HEIGHT,
        DEPTH,
    );

    let mut u_prev =
        FdmGrid::new(dims.clone());

    let mut u_curr =
        FdmGrid::new(dims.clone());

    let mut u_next =
        FdmGrid::new(dims.clone());

    // Initial condition: A Gaussian ripple in the center of the volume
    let initial_u =
        |x: usize,
         y: usize,
         z: usize| {

            let dx_cen = x as f64
                - (WIDTH / 2) as f64;

            let dy_cen = y as f64
                - (HEIGHT / 2) as f64;

            let dz_cen = z as f64
                - (DEPTH / 2) as f64;

            let dist2 = dx_cen * dx_cen
                + dy_cen * dy_cen
                + dz_cen * dz_cen;

            (-dist2 / 50.0).exp()
        };

    // Initialize u_curr and u_prev
    for z in 0 .. DEPTH {

        for y in 0 .. HEIGHT {

            for x in 0 .. WIDTH {

                let val =
                    initial_u(x, y, z);

                u_curr[(x, y, z)] = val;

                u_prev[(x, y, z)] = val;
            }
        }
    }

    let s_x = (C * config.dt
        / config.dx)
        .powi(2);

    let s_y = (C * config.dt
        / config.dy)
        .powi(2);

    let s_z = (C * config.dt
        / config.dz)
        .powi(2);

    println!(
        "Simulating and saving \
         middle-slice frames..."
    );

    for step in 0 .. config.steps {

        // Compute u_next (manual loop for frame saving)
        // Note: For real performance, use the solve_wave_equation_3d in library which uses rayon
        for z in 1 .. DEPTH - 1 {

            for y in 1 .. HEIGHT - 1 {

                for x in 1 .. WIDTH - 1
                {

                    let lap_x = u_curr
                        [(x + 1, y, z)]
                        - 2.0
                            * u_curr[(
                                x, y, z,
                            )]
                        + u_curr[(
                            x - 1,
                            y,
                            z,
                        )];

                    let lap_y = u_curr
                        [(x, y + 1, z)]
                        - 2.0
                            * u_curr[(
                                x, y, z,
                            )]
                        + u_curr[(
                            x,
                            y - 1,
                            z,
                        )];

                    let lap_z = u_curr
                        [(x, y, z + 1)]
                        - 2.0
                            * u_curr[(
                                x, y, z,
                            )]
                        + u_curr[(
                            x,
                            y,
                            z - 1,
                        )];

                    u_next[(x, y, z)] =
                        2.0 * u_curr
                            [(x, y, z)]
                            - u_prev[(
                                x, y, z,
                            )]
                            + s_x
                                * lap_x
                            + s_y
                                * lap_y
                            + s_z
                                * lap_z;
                }
            }
        }

        // Swap grids
        std::mem::swap(
            &mut u_prev,
            &mut u_curr,
        );

        std::mem::swap(
            &mut u_curr,
            &mut u_next,
        );

        // Save a frame (middle XY slice) every 10 steps
        if step % 10 == 0 {

            #[cfg(feature = "output")]
            {

                let frame_num =
                    step / 10;

                let plot_path = format!("ripple_3d_slice_frame_{:02}.png", frame_num);

                let z_slice = DEPTH / 2;

                let array = fdm_grid_3d_slice_to_array2(&u_curr, z_slice);

                let mut plot_config =
                    PlotConfig::default(
                    );

                plot_config.caption = format!(
                    "3D Wave Ripple - \
                     XY Slice (z={}) \
                     - Step {}",
                    z_slice, step
                );

                if let Err(e) =
                    plot_heatmap_2d(
                        &array,
                        &plot_path,
                        Some(
                            plot_config,
                        ),
                    )
                {

                    eprintln!("Failed to save frame {}: {}", frame_num, e);
                } else {

                    println!(
                        "Saved frame \
                         {}",
                        plot_path
                    );
                }
            }
        }
    }

    println!("Simulation finished.");
}

/// Helper function to extract a 2D XY slice from a 3D FdmGrid as ndarray::Array2
#[cfg(feature = "output")]

fn fdm_grid_3d_slice_to_array2(
    grid: &FdmGrid<f64>,
    z: usize,
) -> Array2<f64> {

    let (width, height, _) = match grid
        .dimensions()
    {
        | Dimensions::D3(w, h, d) => {
            (*w, *h, *d)
        },
        | _ => {
            panic!("Expected a 3D grid")
        },
    };

    let mut array =
        Array2::zeros((height, width));

    for y in 0 .. height {

        for x in 0 .. width {

            array[[y, x]] =
                grid[(x, y, z)];
        }
    }

    array
}
