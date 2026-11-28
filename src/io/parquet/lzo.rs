//! LZO1X-1 block compression/decompression and Hadoop-LZO framing for Parquet (`CODEC_LZO = 3`).

use std::collections::HashMap;

/// Compresses `input` into a raw LZO1X-1 block, then wraps it in Hadoop-LZO
/// framing (`[u32 BE uncompressed_len][u32 BE compressed_len][lzo1x_block]`).
pub(super) fn compress_hadoop(input: &[u8]) -> Vec<u8> {
    if input.is_empty() {
        return Vec::new();
    }
    let block = compress_lzo1x(input);
    let mut out = Vec::with_capacity(8 + block.len());
    out.extend_from_slice(&(input.len() as u32).to_be_bytes());
    out.extend_from_slice(&(block.len() as u32).to_be_bytes());
    out.extend_from_slice(&block);
    out
}

/// Decompresses a Parquet LZO page (supports both Hadoop-LZO framing and raw LZO1X).
pub(super) fn decompress_hadoop(
    input: &[u8],
    expected: usize,
) -> Result<Vec<u8>, String> {
    if expected == 0 {
        return Ok(Vec::new());
    }
    // Check if input uses Hadoop-LZO framing:
    // [4-byte BE total uncompressed chunk size][4-byte BE compressed sub-block size][payload]...
    if input.len() >= 8 {
        let mut pos = 0_usize;
        let mut out = Vec::with_capacity(expected);
        let mut framed_ok = true;
        while out.len() < expected && pos + 8 <= input.len() {
            let u_len =
                u32::from_be_bytes([input[pos], input[pos + 1], input[pos + 2], input[pos + 3]])
                    as usize;
            let c_len = u32::from_be_bytes([
                input[pos + 4],
                input[pos + 5],
                input[pos + 6],
                input[pos + 7],
            ]) as usize;
            pos += 8;
            if u_len == 0 || u_len > expected - out.len() || pos + c_len > input.len() {
                framed_ok = false;
                break;
            }
            if c_len == u_len {
                out.extend_from_slice(&input[pos..pos + c_len]);
            } else {
                let dec = decompress_lzo1x(&input[pos..pos + c_len], u_len)?;
                out.extend_from_slice(&dec);
            }
            pos += c_len;
        }
        if framed_ok && out.len() == expected {
            return Ok(out);
        }
    }
    decompress_lzo1x(input, expected)
}

/// Compresses `input` as a single LZO1X-1 block.
pub(super) fn compress_lzo1x(input: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(input.len() + input.len() / 16 + 64);
    let len = input.len();
    if len <= 13 {
        emit_initial_literals(&mut out, input);
        out.extend_from_slice(&[17, 0, 0]); // M4 EOF marker
        return out;
    }

    let mut dict: HashMap<u32, usize> = HashMap::new();
    let mut ip = 0_usize;
    let mut anchor = 0_usize;
    let limit = len.saturating_sub(5);

    while ip < limit {
        let key = u32::from_le_bytes([input[ip], input[ip + 1], input[ip + 2], input[ip + 3]]);
        let prev = dict.insert(key, ip);
        let Some(m_pos) = prev else {
            ip += 1;
            continue;
        };
        let dist = ip - m_pos;
        if dist == 0 || dist > 49_151 {
            ip += 1;
            continue;
        }
        // Check match
        if input[m_pos..m_pos + 3] != input[ip..ip + 3] {
            ip += 1;
            continue;
        }
        let mut m_len = 3_usize;
        while ip + m_len < len && input[m_pos + m_len] == input[ip + m_len] {
            m_len += 1;
        }
        if dist > 16_384 && m_len < 4 {
            ip += 1;
            continue;
        }

        let lit_len = ip - anchor;
        if anchor == 0 && lit_len > 0 {
            emit_initial_literals(&mut out, &input[anchor..ip]);
        } else if lit_len > 0 {
            emit_literals(&mut out, &input[anchor..ip]);
        }

        emit_match(&mut out, dist, m_len);
        ip += m_len;
        anchor = ip;
    }

    if anchor < len {
        if anchor == 0 {
            emit_initial_literals(&mut out, input);
        } else {
            emit_literals(&mut out, &input[anchor..]);
        }
    }
    // LZO1X EOF marker: M4 with distance = 16384 (encoded as [17, 0, 0])
    out.extend_from_slice(&[17, 0, 0]);
    out
}

fn emit_initial_literals(
    out: &mut Vec<u8>,
    lits: &[u8],
) {
    let t = lits.len();
    if t <= 238 {
        out.push((t + 17) as u8);
    } else {
        out.push(0);
        let mut rem = t - 238;
        while rem > 255 {
            out.push(0);
            rem -= 255;
        }
        out.push(rem as u8);
    }
    out.extend_from_slice(lits);
}

fn emit_literals(
    out: &mut Vec<u8>,
    lits: &[u8],
) {
    let t = lits.len();
    if t <= 3 {
        // Attach to previous token's low 2 bits if possible, else emit standalone
        // Standalone literal run: 0x00..=0x0f -> t - 3
        // Note: when t <= 3, we can encode as M1/literal? Actually in LZO1X,
        // token `0..=15` means `t = 3 + token` (or `0` + extra bytes for > 18).
        // If `t <= 3`, wait: we didn't pack it into the previous match's 2-bit tail,
        // so we can patch the previous match's last byte (low 2 bits)!
        // Let's check: every match ends with 2 bytes `[d_low_6_or_8, ...]`, where the
        // first byte of the 2-byte LE distance has `t` in bits 0..2!
        // Since we always emit `0` in bits 0..2 of the match trailer, we can OR `t` into
        // `out[out.len() - 2]`!
        let idx = out.len() - 2;
        out[idx] |= (t as u8) & 0x03;
    } else if t <= 18 {
        out.push((t - 3) as u8);
    } else {
        out.push(0);
        let mut rem = t - 18;
        while rem > 255 {
            out.push(0);
            rem -= 255;
        }
        out.push(rem as u8);
    }
    out.extend_from_slice(lits);
}

fn emit_match(
    out: &mut Vec<u8>,
    dist: usize,
    m_len: usize,
) {
    if dist <= 16_384 {
        // M3 (2..33)
        let d = dist - 1;
        if m_len <= 33 {
            out.push((32 | (m_len - 2)) as u8);
        } else {
            out.push(32);
            let mut rem = m_len - 33;
            while rem > 255 {
                out.push(0);
                rem -= 255;
            }
            out.push(rem as u8);
        }
        out.push(((d & 0x3f) << 2) as u8);
        out.push((d >> 6) as u8);
    } else {
        // M4 (16385..=49151)
        let d = dist - 16_384;
        let bit14 = ((d >> 11) & 0x08) as u8;
        if m_len <= 9 {
            out.push(16 | bit14 | ((m_len - 2) as u8));
        } else {
            out.push(16 | bit14);
            let mut rem = m_len - 9;
            while rem > 255 {
                out.push(0);
                rem -= 255;
            }
            out.push(rem as u8);
        }
        out.push(((d & 0x3f) << 2) as u8);
        out.push((d >> 6) as u8);
    }
}

/// Decompresses a raw LZO1X block into `expected` bytes.
pub(super) fn decompress_lzo1x(
    input: &[u8],
    expected: usize,
) -> Result<Vec<u8>, String> {
    if input.is_empty() {
        return if expected == 0 {
            Ok(Vec::new())
        } else {
            Err("empty LZO1X input".to_owned())
        };
    }
    let mut out = Vec::with_capacity(expected);
    let mut ip = 0_usize;

    let read_ext = |ip: &mut usize, base: usize| -> Result<usize, String> {
        let mut len = base;
        loop {
            let b = *input.get(*ip).ok_or("LZO1X stream truncated")?;
            *ip += 1;
            if b == 0 {
                len = len.checked_add(255).ok_or("LZO1X length overflow")?;
            } else {
                len = len
                    .checked_add(usize::from(b))
                    .ok_or("LZO1X length overflow")?;
                return Ok(len);
            }
        }
    };

    let mut t = usize::from(input[ip]);
    ip += 1;
    if t > 17 {
        t -= 17;
        let end = ip + t;
        let lit = input
            .get(ip..end)
            .ok_or("LZO1X initial literals truncated")?;
        out.extend_from_slice(lit);
        ip = end;
        t = usize::from(
            *input
                .get(ip)
                .ok_or("LZO1X truncated after initial literals")?,
        );
        ip += 1;
    }

    'outer: loop {
        if t < 16 {
            let lit_len = if t == 0 {
                read_ext(&mut ip, 15)? + 3
            } else {
                t + 3
            };
            let end = ip + lit_len;
            let lit = input.get(ip..end).ok_or("LZO1X literal run truncated")?;
            out.extend_from_slice(lit);
            ip = end;
            t = usize::from(*input.get(ip).ok_or("LZO1X truncated after literals")?);
            ip += 1;
            if t < 16 {
                // M1 match after literal run
                let m_pos = out
                    .len()
                    .checked_sub(
                        1 + 0x0800
                            + (t >> 2)
                            + (usize::from(*input.get(ip).ok_or("LZO1X M1 truncated")?) << 2),
                    )
                    .ok_or("LZO1X M1 back-reference out of range")?;
                ip += 1;
                for i in 0..3 {
                    let b = out[m_pos + i];
                    out.push(b);
                }
                let tail = t & 3;
                if tail > 0 {
                    let end = ip + tail;
                    let lit = input.get(ip..end).ok_or("LZO1X M1 tail truncated")?;
                    out.extend_from_slice(lit);
                    ip = end;
                }
                t = usize::from(*input.get(ip).ok_or("LZO1X truncated")?);
                ip += 1;
                continue;
            }
        }

        loop {
            let (m_len, dist, tail) = if t >= 64 {
                // M2
                let b1 = usize::from(*input.get(ip).ok_or("LZO1X M2 truncated")?);
                ip += 1;
                let dist = 1 + ((t >> 2) & 0x07) + (b1 << 3);
                let m_len = (t >> 5) - 1 + 2;
                (m_len, dist, t & 3)
            } else if t >= 32 {
                // M3
                let mut m_len = t & 0x1f;
                if m_len == 0 {
                    m_len = read_ext(&mut ip, 31)?;
                }
                m_len += 2;
                let b0 = usize::from(*input.get(ip).ok_or("LZO1X M3 truncated")?);
                let b1 = usize::from(*input.get(ip + 1).ok_or("LZO1X M3 truncated")?);
                ip += 2;
                let dist = 1 + (b0 >> 2) + (b1 << 6);
                (m_len, dist, b0 & 3)
            } else if t >= 16 {
                // M4
                let mut m_len = t & 0x07;
                if m_len == 0 {
                    m_len = read_ext(&mut ip, 7)?;
                }
                m_len += 2;
                let b0 = usize::from(*input.get(ip).ok_or("LZO1X M4 truncated")?);
                let b1 = usize::from(*input.get(ip + 1).ok_or("LZO1X M4 truncated")?);
                ip += 2;
                let d_low = (b0 >> 2) + (b1 << 6);
                let dist = 16_384 + ((t & 0x08) << 11) + d_low;
                if dist == 16_384 {
                    break 'outer; // EOF
                }
                (m_len, dist, b0 & 3)
            } else {
                // M1
                let b1 = usize::from(*input.get(ip).ok_or("LZO1X M1 truncated")?);
                ip += 1;
                let dist = 1 + (t >> 2) + (b1 << 2);
                (2, dist, t & 3)
            };

            if dist > out.len() {
                return Err("LZO1X match references before start of output".to_owned());
            }
            for src in (out.len() - dist..).take(m_len) {
                let b = out[src];
                out.push(b);
            }
            if tail > 0 {
                let end = ip + tail;
                let lit = input.get(ip..end).ok_or("LZO1X match tail truncated")?;
                out.extend_from_slice(lit);
                ip = end;
                t = usize::from(*input.get(ip).ok_or("LZO1X truncated")?);
                ip += 1;
                if t >= 16 {
                    continue;
                }
                break;
            }
            t = usize::from(*input.get(ip).ok_or("LZO1X truncated")?);
            ip += 1;
            if t < 16 {
                break;
            }
        }
    }

    if out.len() != expected {
        return Err(format!(
            "LZO1X decompressed {} bytes, expected {expected}",
            out.len()
        ));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lzo_round_trip_repetitive_and_small() {
        for data in [
            Vec::new(),
            b"abc".to_vec(),
            b"the quick brown fox jumps over the lazy dog".to_vec(),
            (0..2000)
                .flat_map(|i| format!("row-{}-val-{}\n", i % 13, i).into_bytes())
                .collect::<Vec<u8>>(),
        ] {
            let enc = compress_hadoop(&data);
            let dec = decompress_hadoop(&enc, data.len()).unwrap();
            assert_eq!(dec, data);
        }
    }
}
