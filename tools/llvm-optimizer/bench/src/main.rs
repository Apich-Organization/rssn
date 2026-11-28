use rustfft::FftPlanner;
use rustfft::num_complex::Complex64;
use std::hint::black_box;
use std::time::Instant;

fn main() {
    let mut args = std::env::args().skip(1);
    let size: usize = args.next().expect("FFT size").parse().unwrap();
    let rounds: usize = args.next().expect("rounds").parse().unwrap();
    let fft = FftPlanner::new().plan_fft_forward(size);
    let batch = (4096 / size).clamp(1, 32);
    let input: Vec<_> = (0..size * batch)
        .map(|i| Complex64::new((i % 127) as f64 / 127.0, (i % 59) as f64 / 59.0))
        .collect();
    let mut output = vec![Complex64::default(); input.len()];
    let mut scratch = vec![Complex64::default(); fft.get_immutable_scratch_len()];

    for _ in 0..64 {
        fft.process_immutable_with_scratch(&input, &mut output, &mut scratch);
    }
    let start = Instant::now();
    for _ in 0..rounds {
        fft.process_immutable_with_scratch(
            black_box(&input),
            black_box(&mut output),
            black_box(&mut scratch),
        );
    }
    let nanos = start.elapsed().as_nanos();

    // Include the entire transformed batch in the bitwise comparison.
    let checksum = output.iter().fold(0_u64, |hash, z| {
        (hash.rotate_left(7) ^ z.re.to_bits()).wrapping_mul(0x9e3779b97f4a7c15) ^ z.im.to_bits()
    });
    println!("{nanos} {checksum:016x}");
}
