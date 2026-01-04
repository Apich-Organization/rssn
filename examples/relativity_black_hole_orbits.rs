use plotters::prelude::*;
use rssn::physics::physics_sim::geodesic_relativity::{GeodesicParameters, run_geodesic_simulation};

fn main() -> Result<
    (),
    Box<dyn std::error::Error>,
> {

    let black_hole_mass = 1.0;

    // 1. Define simulation scenarios
    let stable_orbit_params =
        GeodesicParameters {
            black_hole_mass,
            initial_state: [
                10.0, 0.0, 0.0, 0.035,
            ],
            proper_time_end: 1500.0,
            initial_dt: 0.1,
        };

    let plunging_orbit_params =
        GeodesicParameters {
            black_hole_mass,
            initial_state: [
                10.0, 0.0, 0.0, 0.02,
            ],
            proper_time_end: 500.0,
            initial_dt: 0.1,
        };

    let photon_orbit_params =
        GeodesicParameters {
            black_hole_mass,
            initial_state: [
                10.0, -1.0, 0.0, 0.03,
            ],
            proper_time_end: 50.0,
            initial_dt: 0.01,
        };

    let scenarios = vec![
        (
            "Stable Orbit",
            stable_orbit_params,
            &CYAN,
        ),
        (
            "Plunging Orbit",
            plunging_orbit_params,
            &RED,
        ),
        (
            "Photon Orbit",
            photon_orbit_params,
            &YELLOW,
        ),
    ];

    // 2. Setup plotting area
    let output_file =
        "black_hole_orbits.png";

    let root = BitMapBackend::new(
        output_file,
        (1024, 1024),
    )
    .into_drawing_area();

    root.fill(&BLACK)?;

    let mut chart =
        ChartBuilder::on(&root)
            .caption(
                "Schwarzschild Black \
                 Hole Orbits",
                ("sans-serif", 50)
                    .into_font()
                    .color(&WHITE),
            )
            .margin(20)
            .x_label_area_size(40)
            .y_label_area_size(40)
            .build_cartesian_2d(
                -15.0 .. 15.0,
                -15.0 .. 15.0,
            )?;

    chart
        .configure_mesh()
        .disable_mesh()
        .axis_style(&WHITE)
        .label_style(
            ("sans-serif", 15)
                .into_font()
                .color(&WHITE),
        )
        .draw()?;

    // 3. Draw Event Horizon (2M radius)
    let horizon_radius =
        2.0 * black_hole_mass;

    let horizon_points: Vec<(
        f64,
        f64,
    )> = (0 ..= 100)
        .map(|i| {

            let angle = i as f64
                * 2.0
                * std::f64::consts::PI
                / 100.0;

            (
                horizon_radius
                    * angle.cos(),
                horizon_radius
                    * angle.sin(),
            )
        })
        .collect();

    chart
        .draw_series(std::iter::once(
            Polygon::new(
                horizon_points.clone(),
                WHITE
                    .mix(0.2)
                    .filled(),
            ),
        ))?
        .label("Event Horizon (2M)")
        .legend(|(x, y)| {

            Circle::new(
                (x, y),
                5,
                WHITE
                    .mix(0.2)
                    .filled(),
            )
        });

    chart.draw_series(
        std::iter::once(
            PathElement::new(
                horizon_points,
                WHITE.stroke_width(1),
            ),
        ),
    )?;

    // 4. Run simulations and plot results
    for (name, params, color) in
        scenarios
    {

        println!(
            "Simulating {}...",
            name
        );

        let path =
            run_geodesic_simulation(
                &params,
            );

        chart
            .draw_series(
                LineSeries::new(
                    path,
                    color.stroke_width(
                        2,
                    ),
                ),
            )?
            .label(name)
            .legend(move |(x, y)| {

                PathElement::new(
                    vec![
                        (x, y),
                        (x + 20, y),
                    ],
                    color.stroke_width(
                        2,
                    ),
                )
            });
    }

    // 5. Draw Legend
    chart
        .configure_series_labels()
        .background_style(
            &BLACK.mix(0.8),
        )
        .border_style(&WHITE)
        .label_font(
            ("sans-serif", 20)
                .into_font()
                .color(&WHITE),
        )
        .draw()?;

    root.present()?;

    println!(
        "Simulation complete. Plot \
         saved to {}",
        output_file
    );

    Ok(())
}
