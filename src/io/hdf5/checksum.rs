//! HDF5 Jenkins `lookup3` checksum (`H5_checksum_lookup3`).
//!
//! Used by HDF5 superblock v2/v3, object header v2, fractal heap headers and
//! direct/indirect blocks, version-2 B-trees, fixed arrays, and extensible arrays.

#[inline]
const fn rot(
    x: u32,
    k: u32,
) -> u32 {
    x.rotate_left(k)
}

#[inline]
const fn mix(
    a: &mut u32,
    b: &mut u32,
    c: &mut u32,
) {
    *a = a.wrapping_sub(*c);
    *a ^= rot(*c, 4);
    *c = c.wrapping_add(*b);

    *b = b.wrapping_sub(*a);
    *b ^= rot(*a, 6);
    *a = a.wrapping_add(*c);

    *c = c.wrapping_sub(*b);
    *c ^= rot(*b, 8);
    *b = b.wrapping_add(*a);

    *a = a.wrapping_sub(*c);
    *a ^= rot(*c, 16);
    *c = c.wrapping_add(*b);

    *b = b.wrapping_sub(*a);
    *b ^= rot(*a, 19);
    *a = a.wrapping_add(*c);

    *c = c.wrapping_sub(*b);
    *c ^= rot(*b, 4);
    *b = b.wrapping_add(*a);
}

#[inline]
const fn final_mix(
    a: &mut u32,
    b: &mut u32,
    c: &mut u32,
) {
    *c ^= *b;
    *c = c.wrapping_sub(rot(*b, 14));

    *a ^= *c;
    *a = a.wrapping_sub(rot(*c, 11));

    *b ^= *a;
    *b = b.wrapping_sub(rot(*a, 25));

    *c ^= *b;
    *c = c.wrapping_sub(rot(*b, 16));

    *a ^= *c;
    *a = a.wrapping_sub(rot(*c, 4));

    *b ^= *a;
    *b = b.wrapping_sub(rot(*a, 14));

    *c ^= *b;
    *c = c.wrapping_sub(rot(*b, 24));
}

/// Computes the HDF5 Jenkins `lookup3` checksum (`hashlittle` with `initval = 0`).
pub(crate) fn lookup3(data: &[u8]) -> u32 {
    let mut length = data.len();
    let initval = 0xdead_beef_u32.wrapping_add(length as u32);
    let mut a = initval;
    let mut b = initval;
    let mut c = initval;

    let mut pos = 0_usize;
    while length > 12 {
        a = a.wrapping_add(u32::from_le_bytes([
            data[pos],
            data[pos + 1],
            data[pos + 2],
            data[pos + 3],
        ]));
        b = b.wrapping_add(u32::from_le_bytes([
            data[pos + 4],
            data[pos + 5],
            data[pos + 6],
            data[pos + 7],
        ]));
        c = c.wrapping_add(u32::from_le_bytes([
            data[pos + 8],
            data[pos + 9],
            data[pos + 10],
            data[pos + 11],
        ]));
        mix(&mut a, &mut b, &mut c);
        pos += 12;
        length -= 12;
    }

    let rem = &data[pos..];
    match length {
        | 12 => {
            c = c.wrapping_add(u32::from_le_bytes([rem[8], rem[9], rem[10], rem[11]]));
            b = b.wrapping_add(u32::from_le_bytes([rem[4], rem[5], rem[6], rem[7]]));
            a = a.wrapping_add(u32::from_le_bytes([rem[0], rem[1], rem[2], rem[3]]));
        },
        | 11 => {
            c = c.wrapping_add(u32::from(rem[10]) << 16);
            c = c.wrapping_add(u32::from(rem[9]) << 8);
            c = c.wrapping_add(u32::from(rem[8]));
            b = b.wrapping_add(u32::from_le_bytes([rem[4], rem[5], rem[6], rem[7]]));
            a = a.wrapping_add(u32::from_le_bytes([rem[0], rem[1], rem[2], rem[3]]));
        },
        | 10 => {
            c = c.wrapping_add(u32::from(rem[9]) << 8);
            c = c.wrapping_add(u32::from(rem[8]));
            b = b.wrapping_add(u32::from_le_bytes([rem[4], rem[5], rem[6], rem[7]]));
            a = a.wrapping_add(u32::from_le_bytes([rem[0], rem[1], rem[2], rem[3]]));
        },
        | 9 => {
            c = c.wrapping_add(u32::from(rem[8]));
            b = b.wrapping_add(u32::from_le_bytes([rem[4], rem[5], rem[6], rem[7]]));
            a = a.wrapping_add(u32::from_le_bytes([rem[0], rem[1], rem[2], rem[3]]));
        },
        | 8 => {
            b = b.wrapping_add(u32::from_le_bytes([rem[4], rem[5], rem[6], rem[7]]));
            a = a.wrapping_add(u32::from_le_bytes([rem[0], rem[1], rem[2], rem[3]]));
        },
        | 7 => {
            b = b.wrapping_add(u32::from(rem[6]) << 16);
            b = b.wrapping_add(u32::from(rem[5]) << 8);
            b = b.wrapping_add(u32::from(rem[4]));
            a = a.wrapping_add(u32::from_le_bytes([rem[0], rem[1], rem[2], rem[3]]));
        },
        | 6 => {
            b = b.wrapping_add(u32::from(rem[5]) << 8);
            b = b.wrapping_add(u32::from(rem[4]));
            a = a.wrapping_add(u32::from_le_bytes([rem[0], rem[1], rem[2], rem[3]]));
        },
        | 5 => {
            b = b.wrapping_add(u32::from(rem[4]));
            a = a.wrapping_add(u32::from_le_bytes([rem[0], rem[1], rem[2], rem[3]]));
        },
        | 4 => {
            a = a.wrapping_add(u32::from_le_bytes([rem[0], rem[1], rem[2], rem[3]]));
        },
        | 3 => {
            a = a.wrapping_add(u32::from(rem[2]) << 16);
            a = a.wrapping_add(u32::from(rem[1]) << 8);
            a = a.wrapping_add(u32::from(rem[0]));
        },
        | 2 => {
            a = a.wrapping_add(u32::from(rem[1]) << 8);
            a = a.wrapping_add(u32::from(rem[0]));
        },
        | 1 => {
            a = a.wrapping_add(u32::from(rem[0]));
        },
        | 0 => return c,
        | _ => unreachable!(),
    }

    final_mix(&mut a, &mut b, &mut c);
    c
}

#[cfg(test)]
mod tests {
    use super::lookup3;

    #[test]
    fn lookup3_known_vectors() {
        assert_eq!(lookup3(b""), 0xdead_beef);
        assert_eq!(lookup3(b"Four score and seven years ago"), 0x1777_0551);
        assert_eq!(lookup3(b"a"), 1490454280);
    }
}
