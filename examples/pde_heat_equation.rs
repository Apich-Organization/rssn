//! 2D Heat Equation PDE Solver Example
//!
//! This example demonstrates how to solve a partial differential equation (PDE),
//! specifically the 2D heat equation, on a grid.

use ndarray::Array2;
use rssn::output::plotting::PlotConfig;
use rssn::output::plotting::plot_heatmap_2d;

// NOTE: The rssn library currently contains a placeholder for a general PDE solver.
// This example demonstrates a direct implementation of a finite difference method
// for the 2D heat equation to showcase how such a simulation can be performed.
// The tools used here, like ndarray, are part of the broader ecosystem that rssn
// aims to integrate with.

/// Updates the grid for one time step using the finite difference method.

fn update_grid(
    grid: &Array2<f64>,
    alpha: f64,
    dt: f64,
    dx: f64,
    dy: f64,
) -> Array2<f64> {

    let mut new_grid = grid.clone();

    let (height, width) = grid.dim();

    // Stability constant. For the explicit method, we need c <= 0.5 for stability.
    let c_x = alpha * dt / (dx * dx);

    let c_y = alpha * dt / (dy * dy);

    for i in 1 .. height - 1 {

        for j in 1 .. width - 1 {

            let laplacian = c_x
                * (grid[[i + 1, j]]
                    - 2.0
                        * grid[[i, j]]
                    + grid[[i - 1, j]])
                + c_y
                    * (grid
                        [[i, j + 1]]
                        - 2.0
                            * grid[[
                                i, j,
                            ]]
                        + grid[[
                            i,
                            j - 1,
                        ]]);

            new_grid[[i, j]] = grid
                [[i, j]]
                + laplacian;
        }
    }

    new_grid
}

fn main() {

    println!(
        "=== 2D Heat Equation PDE \
         Example ===\n"
    );

    println!(
        "NOTE: The rssn library \
         currently contains a \
         placeholder for a PDE solver."
    );

    println!(
        "This example demonstrates a \
         direct implementation of a \
         finite difference method for \
         the 2D heat equation.\n"
    );

    // Grid and simulation parameters
    let (width, height) = (50, 50); // Grid dimensions
    let alpha = 0.1; // Thermal diffusivity
    let dx = 0.1; // Spatial step size in x
    let dy = 0.1; // Spatial step size in y
    // Time step, chosen for stability (dt <= dx^2 / (2*alpha))
    let dt = dx * dx / (4.0 * alpha);

    let time_steps = 2000;

    // Create the grid and initialize it
    let mut grid = Array2::<f64>::zeros(
        (height, width),
    );

    // Initial condition: a hot spot in the middle
    let center_h = height / 2;

    let center_w = width / 2;

    grid[[center_h, center_w]] = 100.0;

    grid[[
        center_h + 1,
        center_w,
    ]] = 100.0;

    grid[[
        center_h,
        center_w + 1,
    ]] = 100.0;

    grid[[
        center_h + 1,
        center_w + 1,
    ]] = 100.0;

    // Plot initial state
    println!(
        "Generating plot for initial \
         state (t=0)..."
    );

    let plot_config = PlotConfig {
        width: 800,
        height: 800,
        caption: "2D Heat Equation \
                  (t=0)"
            .to_string(),
        ..Default::default()
    };

    if let Err(e) = plot_heatmap_2d(
        &grid,
        "heat_equation_step_0.png",
        Some(plot_config),
    ) {

        eprintln!(
            "Error generating plot: {}",
            e
        );
    } else {

        println!(
            "Plot saved to \
             heat_equation_step_0.png"
        );
    }

    // Time-stepping loop to solve the PDE
    for step in 0 .. time_steps {

        grid = update_grid(
            &grid, alpha, dt, dx, dy,
        );

        // Generate a plot at certain intervals
        if (step + 1) % 500 == 0 {

            let time =
                (step + 1) as f64 * dt;

            println!(
                "\nGenerating plot \
                 for state at t = \
                 {:.2}...",
                time
            );

            let plot_path = format!("heat_equation_step_{}.png", step + 1);

            let plot_config = PlotConfig {
                width: 800,
                height: 800,
                caption: format!("2D Heat Equation (t={:.2})", time),
                ..Default::default()
            };

            if let Err(e) =
                plot_heatmap_2d(
                    &grid,
                    &plot_path,
                    Some(plot_config),
                )
            {

                eprintln!(
                    "Error generating \
                     plot: {}",
                    e
                );
            } else {

                println!(
                    "Plot saved to {}",
                    plot_path
                );
            }
        }
    }

    println!(
        "\n=== Example Complete ==="
    );
}
