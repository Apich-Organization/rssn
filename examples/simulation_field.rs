//! Two field solvers from `sim`: a 1-D acoustic pulse on a rigid box, and a
//! 1-D Yee FDTD pulse, each checked against what physics says must happen.
//!
//! Run with `cargo run --release --example simulation_field`.

use rssn::sim::acoustics::Acoustics;
use rssn::sim::acoustics::AcousticsConfig;
use rssn::sim::yee_fdtd::Comp;
use rssn::sim::yee_fdtd::CurrentSource;
use rssn::sim::yee_fdtd::Yee;
use rssn::sim::yee_fdtd::YeeConfig;

/// Index of the largest magnitude in `field[lo..hi]`.
fn peak(
    field: &[f64],
    lo: usize,
    hi: usize,
) -> usize {
    (lo..hi).fold(lo, |best, i| {
        if field[i].abs() > field[best].abs() {
            i
        } else {
            best
        }
    })
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // --- linear acoustics: a pressure bump in a rigid 2 m box -------------
    let (cells, h, c) = (200_usize, 0.01_f64, 1.0_f64);
    let mut wave = Acoustics::new(AcousticsConfig::new(&[cells], &[h], c)?)?;
    // Rigid (Neumann) walls everywhere, so the energy is an invariant.
    wave.set_pressure(|p| (-((p[0] - 1.0) / 0.08).powi(2)).exp());

    let energy_0 = wave.energy();
    wave.run(100); // 100 * 0.9 h / c = 0.9 s of travel
    let energy_1 = wave.energy();

    let (left, right) = (peak(&wave.p, 0, cells / 2), peak(&wave.p, cells / 2, cells));
    let x = |i: usize| (i as f64 + 0.5) * h;
    println!("acoustics, 1-D, c = 1, rigid walls");
    println!(
        "  energy {energy_0:.9} -> {energy_1:.9}  (drift {:.1e})",
        (energy_1 - energy_0).abs() / energy_0
    );
    println!(
        "  the bump split into two pulses at x = {:.2} and {:.2}  (want 0.10 and 1.90)",
        x(left),
        x(right)
    );
    println!(
        "  dt = {:.4}, peak pressure {:.4}",
        wave.dt,
        wave.p[right].abs()
    );

    // --- Yee FDTD: a current pulse radiates at c = 1 ----------------------
    let mut yee = Yee::new(YeeConfig::new(&[cells], &[h])?)?;
    yee.add_source(CurrentSource {
        cell: [cells / 2, 0, 0],
        comp: Comp::Ey,
        waveform: Box::new(|t| (-((t - 0.15) / 0.03).powi(2)).exp()),
    });
    // PEC walls reflect, so measure the outgoing pulse before it arrives.
    yee.run(20);
    let first = peak(yee.field(Comp::Ey), cells / 2 + 1, cells - 10);
    yee.run(40);
    let second = peak(yee.field(Comp::Ey), cells / 2 + 1, cells - 10);

    let steps = 40.0;
    let moved = (second - first) as f64;
    println!("\nYee FDTD, 1-D, vacuum");
    println!("  pulse cell {first} -> {second} over {steps:.0} steps");
    println!(
        "  {:.3} cells/step, i.e. speed {:.3} c  (cfl = 0.9, c = 1)",
        moved / steps,
        moved * h / (steps * yee.dt)
    );
    Ok(())
}
