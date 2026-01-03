use rssn::physics::physics_sim::schrodinger_quantum::{
    run_schrodinger_simulation, SchrodingerParameters,
};
use num_complex::Complex;
use std::f64::consts::PI;

#[cfg(feature = "output")]
use rssn::output::plotting::{plot_heatmap_2d, PlotConfig, plot_surface_2d};
#[cfg(feature = "output")]
use ndarray::Array2;

fn main() {
    println!("Starting Quantum Tunneling Demo...");

    // 1. Simulation Constants
    const NX: usize = 256;
    const NY: usize = 128;
    const LX: f64 = 20.0;
    const LY: f64 = 10.0;
    const TIME_STEPS: usize = 400;
    const DT: f64 = 0.05;
    const MASS: f64 = 1.0;
    const HBAR: f64 = 1.0;

    // 2. Setup Potential Barrier
    // A vertical barrier in the middle of the domain
    let mut potential = vec![0.0; NX * NY];
    let barrier_x_center = NX / 2;
    let barrier_width = 10; // in grid points
    let barrier_height = 20.0; // V0

    println!("Setting up potential barrier (Height: {}, Width: {} grid points)...", barrier_height, barrier_width);
    
    for j in 0..NY {
        for i in 0..NX {
            if i >= barrier_x_center - barrier_width/2 && i <= barrier_x_center + barrier_width/2 {
                potential[j * NX + i] = barrier_height;
            }
        }
    }

    // 3. Setup Initial Wave Packet
    // Moving to the right, towards the barrier
    let mut initial_psi = vec![Complex::default(); NX * NY];
    let packet_x0 = LX / 4.0;
    let packet_y0 = LY / 2.0;
    let kx = 3.0; // Momentum in X (Energy approx hbar^2 * k^2 / 2m = 1*9/2 = 4.5)
                  // Note: Barrier height 20.0 > 4.5, so purely classical reflection expected.
                  // Quantum tunneling depends on barrier width. This might be "thick" barrier, 
                  // effectively reflecting most. Let's try to see if we see tunneling or just reflection.
                  // For visible tunneling, maybe lower barrier or thinner.
    let ky = 0.0;
    let sigma = 1.0;

    println!("Initializing Wave Packet (k_x={:.2})...", kx);

    let dx = LX / NX as f64;
    let dy = LY / NY as f64;

    for j in 0..NY {
        for i in 0..NX {
            let x = i as f64 * dx;
            let y = j as f64 * dy;
            
            let dist_sq = (x - packet_x0).powi(2) + (y - packet_y0).powi(2);
            let envelope = (-dist_sq / (2.0 * sigma * sigma)).exp();
            let phase = kx * x + ky * y;
            
            initial_psi[j * NX + i] = Complex::from_polar(envelope, phase);
        }
    }

    // Normalize
    let norm_sq: f64 = initial_psi.iter().map(|c| c.norm_sqr()).sum();
    let norm = (norm_sq * dx * dy).sqrt();
    for c in initial_psi.iter_mut() {
        *c /= norm;
    }

    // 4. Run Simulation
    let params = SchrodingerParameters {
        nx: NX,
        ny: NY,
        lx: LX,
        ly: LY,
        dt: DT,
        time_steps: TIME_STEPS,
        hbar: HBAR,
        mass: MASS,
        potential: potential.clone(),
    };

    println!("Running Simulation ({} steps)...", TIME_STEPS);
    let snapshots = run_schrodinger_simulation(&params, &mut initial_psi).unwrap();
    println!("Simulation finished. Generated {} snapshots.", snapshots.len());

    // 5. Visualization
    #[cfg(feature = "output")]
    {
        println!("Generating frames...");
        let mut plot_config = PlotConfig::default();
        plot_config.width = 1920;
        plot_config.height = 1080;
        plot_config.label_font_size = 30;
        plot_config.caption_font_size = 60;
        
        // Save every snapshot
        for (i, density) in snapshots.iter().enumerate() {
            let frame_idx = i * 10; // snapshot taken every 10 steps in lib
            let path = format!("quantum_tunneling_frame_{:03}.png", i);
            plot_config.caption = format!("Quantum Tunneling - Step {}", frame_idx);

            // We can plot heatmap
            if let Err(e) = plot_heatmap_2d(density, &path, Some(plot_config.clone())) {
                eprintln!("Error saving frame {}: {}", i, e);
            } else {
                println!("Saved {}", path);
            }
        }
        
        // Also save a 3D surface of the moment of impact (approx middle)
        if !snapshots.is_empty() {
             let mid_idx = snapshots.len() / 2;
             let path_3d = "quantum_tunneling_3d.png";
             plot_config.caption = "Wave Packet Hitting Barrier (3D)".to_string();
             plot_config.scale = 0.5;
             if let Err(e) = plot_surface_2d(&snapshots[mid_idx], path_3d, Some(plot_config)) {
                  eprintln!("Error saving 3d plot: {}", e);
             } else {
                  println!("Saved 3D visualization to {}", path_3d);
             }
        }
    }
    
    #[cfg(not(feature = "output"))]
    println!("Output disabled.");
}
