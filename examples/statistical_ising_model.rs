use rssn::physics::physics_sim::ising_statistical;
use plotters::prelude::*;
use std::fs::File;
use std::io::Read;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // 1. Run the simulation (this generates the CSV and NPY files)
    ising_statistical::simulate_ising_phase_transition_scenario()?;

    // 2. Read the generated CSV file
    let mut file = File::open("ising_magnetization_vs_temp.csv")?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;

    let mut data = Vec::new();
    for line in contents.lines().skip(1) {
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() == 2 {
            let temp: f64 = parts[0].parse()?;
            let mag: f64 = parts[1].parse()?;
            data.push((temp, mag));
        }
    }

    // 3. Plot Magnetization vs Temperature
    let output_file = "ising_phase_transition.png";
    let root = BitMapBackend::new(output_file, (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Ising Model 2D: Magnetization vs Temperature", ("sans-serif", 30).into_font())
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(0.0..4.5, 0.0..1.1)?;

    chart.configure_mesh()
        .x_desc("Temperature (T)")
        .y_desc("Magnetization |M|")
        .draw()?;

    chart.draw_series(LineSeries::new(
        data.clone(),
        &RED,
    ))?
    .label("Magnetization")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
    
    // Critical temperature for 2D Ising model is approx 2.269
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(2.269, 0.0), (2.269, 1.0)],
        BLACK.stroke_width(1).style_dash(),
    )))?
    .label("Critical Temp (Tc ≈ 2.27)")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLACK.stroke_width(1).style_dash()));

    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;
    println!("Plot saved to {}", output_file);

    println!("Simulation and plotting complete.");

    Ok(())
}
