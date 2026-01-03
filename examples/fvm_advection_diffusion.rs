#[cfg(feature = "output")]
use ndarray::Array2;
#[cfg(feature = "output")]
use rssn::output::plotting::PlotConfig;
#[cfg(feature = "output")]
use rssn::output::plotting::plot_heatmap_2d;
use rssn::physics::physics_fvm::Mesh2D;
use rssn::physics::physics_fvm::solve_advection_2d;

fn main() {

    println!(
        "Starting FVM Advection \
         Transport Example..."
    );

    // 1. Configuration
    const WIDTH: usize = 120;

    const HEIGHT: usize = 120;

    const DOMAIN_SIZE: (f64, f64) =
        (1.0, 1.0);

    const VELOCITY: (f64, f64) =
        (0.5, 0.3); // Moving top-right
    const TOTAL_TIME: f64 = 1.0;

    const DT: f64 = 0.005;

    // 2. Initialize Mesh with Pollutant (Gaussian Blob)
    println!(
        "Initializing {}x{} mesh...",
        WIDTH, HEIGHT
    );

    let mut mesh = Mesh2D::new(
        WIDTH,
        HEIGHT,
        DOMAIN_SIZE,
        |x, y| {

            // Gaussian blob centered at (0.2, 0.2)
            let cx = 0.2;

            let cy = 0.2;

            let sigma = 0.05;

            let dist_sq = (x - cx)
                .powi(2)
                + (y - cy).powi(2);

            (-dist_sq
                / (2.0 * sigma * sigma))
                .exp()
        },
    );

    // 3. Simulation Loop
    let num_frames = 20;

    let steps_per_frame = (TOTAL_TIME
        / DT
        / num_frames as f64)
        as usize;

    println!(
        "Simulating {} frames, {} \
         steps per frame...",
        num_frames, steps_per_frame
    );

    for frame in 0 ..= num_frames {

        // Save current state
        #[cfg(feature = "output")]
        {

            let disp_grid =
                mesh_to_ndarray(&mesh);

            let path = format!("fvm_advection_frame_{:02}.png", frame);

            let mut conf =
                PlotConfig::default();

            conf.caption = format!(
                "Pollutant Transport \
                 - T={:.2}",
                frame as f64
                    * DT
                    * steps_per_frame
                        as f64
            );

            // Use high res as requested in previous turns
            conf.width = 1024;

            conf.height = 1024;

            if let Err(e) =
                plot_heatmap_2d(
                    &disp_grid,
                    &path,
                    Some(conf),
                )
            {

                eprintln!(
                    "Failed to plot \
                     frame {}: {}",
                    frame, e
                );
            } else {

                println!(
                    "Saved {}",
                    path
                );
            }
        }

        if frame == num_frames {

            break;
        }

        // Run Solver Step
        // Note: solve_advection_2d returns the new state but doesn't update mesh automatically
        let bc = |i, j, w, h| {
            i == 0
                || i == w - 1
                || j == 0
                || j == h - 1
        };

        let new_values =
            solve_advection_2d(
                &mut mesh,
                VELOCITY,
                DT,
                steps_per_frame,
                bc,
            );

        // Update mesh for next iteration
        for (i, cell) in mesh
            .cells
            .iter_mut()
            .enumerate()
        {

            cell.value = new_values[i];
        }
    }

    println!("Simulation complete.");
}

#[cfg(feature = "output")]

fn mesh_to_ndarray(
    mesh: &Mesh2D
) -> Array2<f64> {

    let mut grid = Array2::zeros((
        mesh.height,
        mesh.width,
    ));

    for j in 0 .. mesh.height {

        for i in 0 .. mesh.width {

            // mesh.cells is likely row-major?
            // "idx / width" in physics_fvm implies row-major (j is row index)
            let idx =
                j * mesh.width + i;

            // For plotting, we often want (0,0) at bottom-left, but heatmap usually does matrix indexing (0,0 at top-left)
            // Let's just map usually.
            // If we want cartesian visualization, we might need to flip Y.
            // Let's verify physics_fvm ordering:
            // "let j = idx / width;" -> j is y-index.
            // "center_y = (j + 0.5) * dy" -> y increases with j.
            // So j=0 is bottom (y=0), j=height-1 is top (y=1).
            // Image/Matrix usually has (0,0) at top-left.
            // To visualize nicely as "Map", we probably want to flip vertical or just map it directly.
            // We'll flip Y for plotting so it looks like a Cartesian map.
            grid[[
                mesh.height - 1 - j,
                i,
            ]] = mesh.cells[idx].value;
        }
    }

    grid
}
