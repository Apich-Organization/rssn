#[cfg(feature = "output")]
use ndarray::Array2;
#[cfg(feature = "output")]
use rssn::output::plotting::PlotConfig;
#[cfg(feature = "output")]
use rssn::output::plotting::plot_surface_2d;
use rssn::physics::physics_fdm::Dimensions;
use rssn::physics::physics_fdm::FdmGrid;
use rssn::physics::physics_fdm::FdmSolverConfig2D;
use rssn::physics::physics_fdm::solve_wave_equation_2d;

fn main() {

    println!(
        "Starting 2D Wave Equation \
         Surface (3D Vis) Example..."
    );

    // 1. Configuration for the simulation
    const WIDTH: usize = 80;

    const HEIGHT: usize = 80;

    const C: f64 = 1.0;

    let config = FdmSolverConfig2D {
        width: WIDTH,
        height: HEIGHT,
        dx: 1.0,
        dy: 1.0,
        dt: 0.2,
        steps: 200,
    };

    // 2. Initial state setup: Two overlapping Gaussian ripples
    let initial_u =
        |x: usize, y: usize| {

            let dx1 = x as f64
                - (WIDTH as f64 * 0.4);

            let dy1 = y as f64
                - (HEIGHT as f64 * 0.4);

            let r1 = (-(dx1 * dx1
                + dy1 * dy1)
                / 40.0)
                .exp();

            let dx2 = x as f64
                - (WIDTH as f64 * 0.7);

            let dy2 = y as f64
                - (HEIGHT as f64 * 0.7);

            let r2 = (-(dx2 * dx2
                + dy2 * dy2)
                / 30.0)
                .exp();

            r1 + r2
        };

    let dims =
        Dimensions::D2(WIDTH, HEIGHT);

    let mut u_prev =
        FdmGrid::new(dims.clone());

    let mut u_curr =
        FdmGrid::new(dims.clone());

    let mut u_next =
        FdmGrid::new(dims.clone());

    for y in 0 .. HEIGHT {

        for x in 0 .. WIDTH {

            let val = initial_u(x, y);

            u_curr[(x, y)] = val;

            u_prev[(x, y)] = val;
        }
    }

    let s_x = (C * config.dt
        / config.dx)
        .powi(2);

    let s_y = (C * config.dt
        / config.dy)
        .powi(2);

    println!(
        "Simulating and saving \
         surface frames..."
    );

    for step in 0 .. config.steps {

        for y in 1 .. HEIGHT - 1 {

            for x in 1 .. WIDTH - 1 {

                let lap_x = u_curr
                    [(x + 1, y)]
                    - 2.0
                        * u_curr
                            [(x, y)]
                    + u_curr
                        [(x - 1, y)];

                let lap_y = u_curr
                    [(x, y + 1)]
                    - 2.0
                        * u_curr
                            [(x, y)]
                    + u_curr
                        [(x, y - 1)];

                u_next[(x, y)] = 2.0
                    * u_curr[(x, y)]
                    - u_prev[(x, y)]
                    + s_x * lap_x
                    + s_y * lap_y;
            }
        }

        std::mem::swap(
            &mut u_prev,
            &mut u_curr,
        );

        std::mem::swap(
            &mut u_curr,
            &mut u_next,
        );

        // Save a surface plot frame every 20 steps
        if step % 20 == 0 {

            #[cfg(feature = "output")]
            {

                let frame_num =
                    step / 20;

                let plot_path = format!("wave_surface_frame_{:02}.png", frame_num);

                let array =
                    fdm_grid_to_ndarray(
                        &u_curr,
                    );

                let mut plot_config =
                    PlotConfig::default(
                    );

                plot_config.caption = format!(
                    "Wave Surface - \
                     Step {}",
                    step
                );

                if let Err(e) =
                    plot_surface_2d(
                        &array,
                        &plot_path,
                        Some(
                            plot_config,
                        ),
                    )
                {

                    eprintln!("Failed to save surface {}: {}", frame_num, e);
                } else {

                    println!("Saved surface plot {}", plot_path);
                }
            }
        }
    }

    println!("Simulation finished.");
}

/// Helper function to convert FdmGrid to ndarray::Array2
#[cfg(feature = "output")]

fn fdm_grid_to_ndarray(
    grid: &FdmGrid<f64>
) -> Array2<f64> {

    let (width, height) = match grid
        .dimensions()
    {
        | Dimensions::D2(w, h) => {
            (*w, *h)
        },
        | _ => {

            panic!("Expected a 2D grid")
        },
    };

    let mut array =
        Array2::zeros((height, width));

    for y in 0 .. height {

        for x in 0 .. width {

            array[[y, x]] =
                grid[(x, y)];
        }
    }

    array
}
