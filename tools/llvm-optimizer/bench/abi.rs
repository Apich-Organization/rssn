use std::{hint::black_box, time::Instant};

// Three fat slices and two scalars are the same shape as the observed
// stack-parameter case, with a small leaf body to isolate the call boundary.
#[inline(never)]
unsafe fn rssn_abi_mix(
    a: &[u64],
    b: &[u64],
    c: &[u64],
    reverse: bool,
    key: u64,
) -> u64 {
    let i = (key as usize) & 7;
    // The benchmark passes only slices with at least eight elements.
    let (x, y, z) = unsafe {
        (
            *a.get_unchecked(i),
            *b.get_unchecked(i),
            *c.get_unchecked(i),
        )
    };
    let n = ((a.len() as u64) ^ (b.len() as u64).wrapping_mul(11))
        .wrapping_add((c.len() as u64).wrapping_mul(17));
    let xy = if reverse {
        y.wrapping_sub(x)
    } else {
        x.wrapping_add(y)
    };
    (xy ^ z).rotate_left((key & 63) as u32).wrapping_add(n)
}

#[inline(never)]
fn sample(
    data: &[u64; 72],
    skew: usize,
    rounds: u64,
) -> u64 {
    let a = &data[..17 + (skew & 3)];
    let b = &data[22..41 + ((skew >> 2) & 3)];
    let c = &data[46..69 + ((skew >> 4) & 3)];
    let mut result = 0_u64;
    for i in 0..rounds {
        result = result.wrapping_add(unsafe { rssn_abi_mix(a, b, c, i & 1 != 0, i) });
    }
    result
}

fn main() {
    let rounds: u64 = std::env::args().nth(1).expect("rounds").parse().unwrap();
    let data: [u64; 72] = std::array::from_fn(|v| (v as u64).wrapping_mul(0x9e3779b97f4a7c15));
    let data = black_box(&data);
    let skew = black_box(29);
    black_box(sample(data, skew, 10_000));
    let start = Instant::now();
    let value = sample(data, skew, rounds);
    println!("{} {value:016x}", start.elapsed().as_nanos());
}
