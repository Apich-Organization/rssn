#![allow(dead_code)]
include!(env!("LANE_SOURCE"));

use std::hint::black_box;
mod clock;

use self::clock::now;

fn run(
    ir: &KernelIr,
    rows: usize,
    input: &[f64],
    regs: &mut [f64],
    out: &mut [f64],
) {
    for first in (0..rows).step_by(LANES) {
        let n = (rows - first).min(LANES);
        exec_block(ir, input, rows, &[1.25], first, n, regs);
        for (c, &r) in ir.outputs.iter().enumerate() {
            let start = c * rows + first;
            let source = r as usize * LANES;
            out[start..start + n].copy_from_slice(&regs[source..source + n]);
        }
    }
}
fn verify_pow() {
    let special = [
        0,
        0x8000_0000_0000_0000,
        1,
        0x8000_0000_0000_0001,
        0x0010_0000_0000_0000,
        0x7fef_ffff_ffff_ffff,
        0x7ff0_0000_0000_0000,
        0xfff0_0000_0000_0000,
        0x7ff0_0000_0000_1234,
        0xfff8_1234_0000_0055,
        0x3ff0_0000_0000_0000,
        0xbff0_0000_0000_0000,
    ];
    let mut bits = 129458565026211803u64;
    let inputs: Vec<_> = (0..260)
        .map(|i| {
            bits ^= bits << 13;
            bits ^= bits >> 7;
            bits ^= bits << 17;
            f64::from_bits(if i < special.len() {
                special[i]
            } else {
                bits
            })
        })
        .collect();
    let exponents = [
        i32::MIN,
        i32::MIN + 1,
        -1023,
        -32,
        -17,
        -3,
        -2,
        -1,
        0,
        1,
        2,
        3,
        7,
        16,
        63,
        128,
        1023,
        i32::MAX,
    ];
    for e in exponents {
        for rows in 0..=260 {
            let ir = KernelIr {
                ops: vec![KernelOp::Input(0), KernelOp::PowI(0, e)],
                inputs: 1,
                outputs: vec![1],
            };
            let mut regs = vec![0.0; 2 * LANES];
            let mut out = vec![0.0; rows];
            run(black_box(&ir), rows, &inputs[..rows], &mut regs, &mut out);
            for (x, y) in inputs.iter().zip(out) {
                let v = int_pow(*x, e);
                assert_eq!(
                    v.to_bits(),
                    y.to_bits(),
                    "bits={:x} exponent={e} rows={rows}",
                    x.to_bits()
                );
            }
        }
    }
    println!("all cases exact");
}
fn main() {
    let a = std::env::args().collect::<Vec<_>>();
    if a[1] == "check" {
        verify_pow();
        return;
    }
    let ir = match a[1].as_str() {
        | "everything" => everything(),
        | "smooth" => smooth(),
        | "powi" => KernelIr {
            ops: vec![
                KernelOp::Input(0),
                KernelOp::Input(1),
                KernelOp::PowI(0, -3),
                KernelOp::PowI(1, 11),
                KernelOp::Add(2, 3),
            ],
            inputs: 2,
            outputs: vec![4],
        },
        | "control" => KernelIr {
            ops: vec![
                KernelOp::Input(0),
                KernelOp::Input(1),
                KernelOp::Mul(0, 0),
                KernelOp::Mul(1, 1),
                KernelOp::Add(2, 3),
            ],
            inputs: 2,
            outputs: vec![4],
        },
        | _ => unreachable!(),
    };
    let rows: usize = a[2].parse().unwrap();
    let rounds: usize = a[3].parse().unwrap();
    let mut bits = 0x1357_2468_feeb_cdefu64;
    let input: Vec<f64> = (0..rows * ir.inputs)
        .map(|_| {
            bits ^= bits << 13;
            bits ^= bits >> 7;
            bits ^= bits << 17;
            0.375 + ((bits >> 11) as f64) * (1.5 / ((1u64 << 53) as f64))
        })
        .collect();
    let mut regs = vec![0.0; ir.ops.len() * LANES];
    let mut out = vec![0.0; rows * ir.outputs.len()];
    run(black_box(&ir), rows, &input, &mut regs, &mut out);
    let start = now();
    for _ in 0..rounds {
        run(black_box(&ir), rows, black_box(&input), &mut regs, &mut out);
        black_box(&out);
    }
    let elapsed = now() - start;
    let hash = out.iter().fold(0xcbf2_9ce4_8422_2325u64, |h, x| {
        (h ^ x.to_bits()).wrapping_mul(1099511628211)
    });
    println!("{elapsed} {hash:016x}");
}
