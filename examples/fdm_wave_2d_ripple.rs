#[cfg(feature = "output")]
use ndarray::Array2;
#[cfg(feature = "output")]
use rssn::output::plotting::PlotConfig;
#[cfg(feature = "output")]
use rssn::output::plotting::plot_heatmap_2d;
use rssn::physics::physics_fdm::Dimensions;
use rssn::physics::physics_fdm::FdmGrid;
use rssn::physics::physics_fdm::FdmSolverConfig2D;
use rssn::physics::physics_fdm::solve_wave_equation_2d;

fn main() {

    println!(
        "Starting 2D Wave Equation \
         Ripple Example..."
    );

    // 1. Configuration for the simulation
    const WIDTH: usize = 100;

    const HEIGHT: usize = 100;

    const C: f64 = 1.0; // Wave speed

    let config = FdmSolverConfig2D {
        width: WIDTH,
        height: HEIGHT,
        dx: 1.0,
        dy: 1.0,
        dt: 0.25, /* Stability condition: dt < dx / (c * sqrt(2)) */
        steps: 150,
    };

    // 2. Initial state setup
    let dims =
        Dimensions::D2(WIDTH, HEIGHT);

    let mut u_prev =
        FdmGrid::new(dims.clone());

    let mut u_curr =
        FdmGrid::new(dims.clone());

    let mut u_next =
        FdmGrid::new(dims.clone());

    // Initial condition: A Gaussian ripple in the center
    let initial_u =
        |x: usize, y: usize| {

            let dx_cen = x as f64
                - (WIDTH / 2) as f64;

            let dy_cen = y as f64
                - (HEIGHT / 2) as f64;

            let dist2 = dx_cen * dx_cen
                + dy_cen * dy_cen;

            (-dist2 / 50.0).exp()
        };

    // Initialize u_curr and u_prev
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
         frames..."
    );

    for step in 0 .. config.steps {

        // Compute u_next
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

        // Boundary conditions (fixed at 0)
        // (Already initialized to 0 by FdmGrid::new)

        // Swap grids
        std::mem::swap(
            &mut u_prev,
            &mut u_curr,
        );

        std::mem::swap(
            &mut u_curr,
            &mut u_next,
        );

        // Save a frame every 10 steps
        if step % 10 == 0 {

            #[cfg(feature = "output")]
            {

                let frame_num =
                    step / 10;

                let plot_path = format!(
                    "ripple_frame_{:\
                     02}.png",
                    frame_num
                );

                let array =
                    fdm_grid_to_ndarray(
                        &u_curr,
                    );

                let mut plot_config =
                    PlotConfig::default(
                    );

                plot_config.caption = format!(
                    "Wave Ripple - \
                     Step {}",
                    step
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
