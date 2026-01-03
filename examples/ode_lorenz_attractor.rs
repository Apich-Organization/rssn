//! Lorenz Attractor ODE Solver Example
//!
//! This example demonstrates how to solve a system of ordinary differential equations (ODEs)
//! for the Lorenz attractor using the numerical ODE solver functionality in `rssn`.

use rssn::input::parser::parse_expr;
use rssn::numerical::ode::OdeSolverMethod;
use rssn::numerical::ode::solve_ode_system;
use rssn::output::plotting::PlotConfig;
use rssn::output::plotting::plot_3d_path_from_points;

fn main() {

    println!(
        "=== Lorenz Attractor ODE \
         Solver Example ===\n"
    );

    // Lorenz attractor parameters
    let sigma = 10.0;

    let rho = 28.0;

    let beta = 4.0 / 3.0;

    // Define the Lorenz system of ODEs as strings.
    // The solver expects variables y0, y1, y2, ... for the state vector
    // and 'x' for the independent variable (time in this case).
    // Let: y0 = x, y1 = y, y2 = z
    let ode_system_strs = [
        format!(
            "{} * (y1 - y0)",
            sigma
        ), // dx/dt = sigma * (y - x)
        format!(
            "y0 * ({} - y2) - y1",
            rho
        ), // dy/dt = x * (rho - z) - y
        format!(
            "y0 * y1 - {} * y2",
            beta
        ), // dz/dt = x * y - beta * z
    ];

    println!("Lorenz System of ODEs:");

    println!(
        "  dx/dt = {}",
        ode_system_strs[0]
    );

    println!(
        "  dy/dt = {}",
        ode_system_strs[1]
    );

    println!(
        "  dz/dt = {}",
        ode_system_strs[2]
    );

    // Parse the strings into expressions
    let funcs = ode_system_strs
        .iter()
        .map(|s| match parse_expr(s) {
            Ok(("", expr)) => expr,
            Ok((rem, _)) => panic!("Unparsed input for expression: '{}'", rem),
            Err(e) => panic!("Failed to parse expression '{}': {:?}", s, e),
        })
        .collect::<Vec<_>>();

    // Initial conditions [x(0), y(0), z(0)]
    let y0 = &[0.0, 1.0, 1.05];

    // Time range
    let x_range = (0.0, 50000.0);

    // Number of steps
    let num_steps = 5000000;

    println!(
        "\nInitial Conditions: [x, y, \
         z] = {:?}",
        y0
    );

    println!(
        "Time Range: {:?}",
        x_range
    );

    println!(
        "Number of Steps: {}",
        num_steps
    );

    println!("Solver: RungeKutta4\n");

    // Solve the ODE system
    match solve_ode_system(
        &funcs,
        y0,
        x_range,
        num_steps,
        OdeSolverMethod::RungeKutta4,
    ) {
        | Ok(results) => {

            println!(
                "Successfully solved \
                 the ODE system."
            );

            println!(
                "Showing the first 5 \
                 and last 5 results \
                 (of {} total):",
                results.len()
            );

            // Print header
            println!(
                "{:<5} | {:<20} | \
                 {:<20} | {:<20}",
                "Step", "x", "y", "z"
            );

            println!("{:-<71}", "");

            // Print first 5 results
            for (i, res) in results
                .iter()
                .take(5)
                .enumerate()
            {

                println!(
                    "{:<5} | {:<20.8} \
                     | {:<20.8} | \
                     {:<20.8}",
                    i,
                    res[0],
                    res[1],
                    res[2]
                );
            }

            println!("...");

            // Print last 5 results
            for (i, res) in results
                .iter()
                .rev()
                .take(5)
                .rev()
                .enumerate()
            {

                let step =
                    results.len() - 5
                        + i;

                println!(
                    "{:<5} | {:<20.8} \
                     | {:<20.8} | \
                     {:<20.8}",
                    step,
                    res[0],
                    res[1],
                    res[2]
                );
            }

            // Plot the results
            println!(
                "\nGenerating plot..."
            );

            let plot_path =
                "lorenz_attractor.png";

            let plot_config = PlotConfig {
                caption: "Lorenz Attractor".to_string(),
                ..Default::default()
            };

            if let Err(e) =
                plot_3d_path_from_points(
                    &results,
                    plot_path,
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
        },
        | Err(e) => {

            eprintln!(
                "Error solving ODE \
                 system: {}",
                e
            );
        },
    }

    println!(
        "\n=== Example Complete ==="
    );
}
