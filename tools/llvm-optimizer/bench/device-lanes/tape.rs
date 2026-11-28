#![allow(dead_code)]
use std::hint::black_box;
include!(env!("TAPE_SOURCE"));

mod clock;

use self::clock::now;

// This runner executes one uncancelled worker and reuses only its scratch.
mod cancellation {
    pub fn is_cancelled() -> bool {
        false
    }
}
mod kernels {
    pub mod reuse {
        use std::{
            cell::RefCell,
            ops::{Deref, DerefMut},
        };
        thread_local! {static V: RefCell<Vec<f64>>=const {RefCell::new(Vec::new())};}
        pub(crate) struct Hold(Vec<f64>);
        impl Deref for Hold {
            type Target = Vec<f64>;
            fn deref(&self) -> &Self::Target {
                &self.0
            }
        }
        impl DerefMut for Hold {
            fn deref_mut(&mut self) -> &mut Self::Target {
                &mut self.0
            }
        }
        impl Drop for Hold {
            fn drop(&mut self) {
                self.0.clear();
                V.with(|v| *v.borrow_mut() = std::mem::take(&mut self.0))
            }
        }
        pub(crate) fn f64_scratch(n: usize) -> Hold {
            let mut v = V.with(|v| std::mem::take(&mut *v.borrow_mut()));
            if v.capacity() < n {
                v.reserve(n);
            }
            Hold(v)
        }
    }
}

fn power_program(exponents: &[i32]) -> Tape {
    let mut insts = vec![Inst::Input(0)];
    insts.extend(exponents.iter().map(|&e| Inst::PowI(0, e)));
    Tape {
        insts,
        args: vec![],
        inputs: 1,
        params: 0,
        outputs: (1..=exponents.len() as u32).collect(),
    }
}

fn verify() {
    let powers = [
        -2147483648,
        -31,
        -3,
        -2,
        -1,
        0,
        1,
        2,
        3,
        16,
        127,
        2147483647,
    ];
    let tape = power_program(&powers);
    let edge = [
        0.0,
        -0.0,
        1.25,
        -1.25,
        f64::from_bits(1),
        -f64::from_bits(1),
        f64::MIN_POSITIVE,
        f64::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::from_bits(0x7ff8_0000_0000_0042),
        f64::from_bits(0x7ff0_0000_0000_2020),
    ];
    let xs: Vec<_> = edge.iter().copied().cycle().take(520).collect();
    for rows in 0..=520 {
        let mut out = vec![vec![0.0; rows]; powers.len()];
        let mut columns: Vec<_> = out.iter_mut().map(Vec::as_mut_slice).collect();
        let present = rows.saturating_sub(3);
        tape.eval_columns_range(&[&xs[..present]], &[], &mut columns, 0, rows);
        for (column, e) in out.iter().zip(powers) {
            for (r, &v) in column.iter().enumerate() {
                let x = xs
                    .get(r)
                    .copied()
                    .filter(|_| r < present)
                    .unwrap_or(f64::NAN);
                assert_eq!(
                    v.to_bits(),
                    int_pow(x, e).to_bits(),
                    "row {r}/{rows}, exponent {e}"
                );
            }
        }
    }
    println!("all cases exact");
}

fn main() {
    let a: Vec<_> = std::env::args().collect();
    if a[1] == "check" {
        verify();
        return;
    }
    let tape = match a[1].as_str() {
        | "everything" => everything(),
        | "smooth" => smooth(),
        | "powi" => power_program(&[-3, 11]),
        | "control" => Tape {
            insts: vec![
                Inst::Input(0),
                Inst::Input(1),
                Inst::Mul(0, 0),
                Inst::Mul(1, 1),
                Inst::Add(2, 3),
            ],
            args: vec![],
            inputs: 2,
            params: 0,
            outputs: vec![4],
        },
        | _ => panic!("unknown case"),
    };
    let rows: usize = a[2].parse().unwrap();
    let rounds: usize = a[3].parse().unwrap();
    let mut state = 0x1357_2468_feeb_cdefu64;
    let x: Vec<f64> = (0..rows * 2)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            0.375 + ((state >> 11) as f64) * (1.5 / ((1u64 << 53) as f64))
        })
        .collect();
    let cols: Vec<_> = x.chunks(rows).take(tape.inputs).collect();
    let mut out = vec![vec![0.0; rows]; tape.outputs.len()];
    let mut view: Vec<_> = out.iter_mut().map(Vec::as_mut_slice).collect();
    tape.eval_columns_range(&cols, &[1.25], &mut view, 0, rows);
    let start = now();
    for _ in 0..rounds {
        black_box(&tape).eval_columns_range(black_box(&cols), &[1.25], &mut view, 0, rows);
        black_box(&view);
    }
    let time = now() - start;
    let digest = out.iter().flatten().fold(0xcbf2_9ce4_8422_2325u64, |h, x| {
        (h ^ x.to_bits()).wrapping_mul(1099511628211)
    });
    println!("{time} {digest:016x}");
}
