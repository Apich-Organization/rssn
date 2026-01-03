use std::fs::File;
use std::io::Read;

use plotters::style::Color;
use plotters::style::RED;
use rssn::output::plotting::PlotConfig;
use rssn::output::plotting::plot_series_2d;
use rssn::physics::physics_sim::ising_statistical;

fn main() -> Result<
    (),
    Box<dyn std::error::Error>,
> {

    // 1. Run the simulation (this generates the CSV and NPY files)
    ising_statistical::simulate_ising_phase_transition_scenario()
        .map_err(|e| e.to_string())?;

    // 2. Read the generated CSV file
    let mut file = File::open(
        "ising_magnetization_vs_temp.\
         csv",
    )?;

    let mut contents = String::new();

    file.read_to_string(&mut contents)?;

    let mut data = Vec::new();

    for line in contents
        .lines()
        .skip(1)
    {

        let parts: Vec<&str> = line
            .split(',')
            .collect();

        if parts.len() == 2 {

            let temp: f64 =
                parts[0].parse()?;

            let mag: f64 =
                parts[1].parse()?;

            data.push((temp, mag));
        }
    }

    // 3. Plot Magnetization vs Temperature using rssn::output::plotting
    let output_file =
        "ising_phase_transition.png";

    let critical_temp_line = vec![
        (2.269, 0.0),
        (2.269, 1.0),
    ];

    let series = vec![
        (
            "Magnetization".to_string(),
            data,
        ),
        (
            "Critical Temp (Tc ≈ 2.27)"
                .to_string(),
            critical_temp_line,
        ),
    ];

    let mut config =
        PlotConfig::default();

    config.caption = "Ising Model 2D: \
                      Magnetization \
                      vs Temperature"
        .to_string();

    config.width = 800;

    config.height = 600;

    config.line_color = RED.to_rgba(); // Base color, though plot_series_2d cycles colors

    plot_series_2d(
        &series,
        output_file,
        Some(config),
    )
    .map_err(|e| {
        format!(
            "Plotting error: {}",
            e
        )
    })?;

    println!(
        "Plot saved to {}",
        output_file
    );

    println!(
        "Simulation and plotting \
         complete."
    );

    Ok(())
}
