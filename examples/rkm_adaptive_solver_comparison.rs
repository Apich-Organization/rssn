use std::time::Instant;

#[cfg(feature = "output")]
use rssn::output::plotting::PlotConfig;
#[cfg(feature = "output")]
use rssn::output::plotting::plot_series_2d;
use rssn::physics::physics_rkm::CashKarp45;
use rssn::physics::physics_rkm::DormandPrince54;
use rssn::physics::physics_rkm::OdeSystem;
use rssn::physics::physics_rkm::VanDerPolSystem;
use rssn::physics::physics_rkm::solve_rk4;

fn main() {

    println!(
        "Starting Runge-Kutta \
         Adaptive Solver Comparison..."
    );

    // Van der Pol oscillator with mu=5.0 (stiff enough to benefit from adaptation)
    let system = VanDerPolSystem {
        mu: 5.0,
    };

    let y0 = &[2.0, 0.0];

    let t_span = (0.0, 30.0);

    let tolerance = (1e-6, 1e-6); // (Absolute, Relative)
    let dt_initial = 0.01;

    // 1. Classical RK4 (Fixed Step)
    // We use a small enough step to be somewhat accurate
    println!(
        "Solving with static RK4 \
         (dt={})...",
        dt_initial
    );

    let start_rk4 = Instant::now();

    let rk4_res = solve_rk4(
        &system,
        y0,
        t_span,
        dt_initial,
    );

    let rk4_duration =
        start_rk4.elapsed();

    println!(
        "RK4 finished in {:?}, steps: \
         {}",
        rk4_duration,
        rk4_res.len()
    );

    // 2. Dormand-Prince 5(4) (Adaptive)
    let solver_dp =
        DormandPrince54::new();

    println!(
        "Solving with Dormand-Prince \
         5(4) (Adaptive, tol={:?})...",
        tolerance
    );

    let start_dp = Instant::now();

    let dp_res = solver_dp.solve(
        &system,
        y0,
        t_span,
        dt_initial,
        tolerance,
    );

    let dp_duration =
        start_dp.elapsed();

    println!(
        "DP 5(4) finished in {:?}, \
         steps: {}",
        dp_duration,
        dp_res.len()
    );

    // 3. Cash-Karp 4(5) (Adaptive)
    let solver_ck =
        CashKarp45::default();

    println!(
        "Solving with Cash-Karp 4(5) \
         (Adaptive, tol={:?})...",
        tolerance
    );

    let start_ck = Instant::now();

    let ck_res = solver_ck.solve(
        &system,
        y0,
        t_span,
        dt_initial,
        tolerance,
    );

    let ck_duration =
        start_ck.elapsed();

    println!(
        "CK 4(5) finished in {:?}, \
         steps: {}",
        ck_duration,
        ck_res.len()
    );

    // Results Analysis
    println!("\nPerformance Summary:");

    println!(
        "RK4: {:?} (Fixed {} steps)",
        rk4_duration,
        rk4_res.len()
    );

    println!(
        "DP 5(4): {:?} ({} steps)",
        dp_duration,
        dp_res.len()
    );

    println!(
        "CK 4(5): {:?} ({} steps)",
        ck_duration,
        ck_res.len()
    );

    // Visualization
    #[cfg(feature = "output")]
    {

        println!(
            "\nGenerating comparison \
             plot..."
        );

        // Helper to convert solver output to (t, y) points
        let to_points = |res: &[(
            f64,
            Vec<f64>,
        )]| {

            res.iter().map(|(t, y)| (*t, y[0])).collect::<Vec<(f64, f64)>>()
        };

        let series = vec![
            (
                "RK4 (Fixed)"
                    .to_string(),
                to_points(&rk4_res),
            ),
            (
                "Dormand-Prince 5(4)"
                    .to_string(),
                to_points(&dp_res),
            ),
            (
                "Cash-Karp 4(5)"
                    .to_string(),
                to_points(&ck_res),
            ),
        ];

        let mut plot_config =
            PlotConfig::default();

        plot_config.caption =
            "Van der Pol Oscillator \
             (mu=5.0) Solver \
             Comparison"
                .to_string();

        plot_config.width = 3840;

        plot_config.height = 2160;

        plot_config.caption_font_size =
            80;

        plot_config.label_font_size =
            40;

        let plot_path =
            "rkm_solver_comparison.png";

        match plot_series_2d(
            &series,
            plot_path,
            Some(plot_config),
        ) {
            | Ok(_) => {

                println!(
                    "Successfully \
                     saved comparison \
                     plot to {}",
                    plot_path
                )
            },
            | Err(e) => {

                eprintln!(
                    "Failed to save \
                     plot: {}",
                    e
                )
            },
        }
    }

    #[cfg(not(feature = "output"))]

    println!(
        "\nOutput feature not \
         enabled. No plot generated."
    );
}
