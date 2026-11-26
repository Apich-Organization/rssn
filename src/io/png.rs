//! Complete pure-Rust PNG image encoder and decoder.
//!
//! Uses [`crate::io::hdf5::inflate`] for zlib compression/decompression and
//! supports:
//! - Grayscale (8/16-bit), Grayscale+Alpha (8/16-bit), RGB (8/16-bit), RGBA (8/16-bit),
//!   and Indexed (`PLTE` + optional `tRNS`) color types
//! - Non-interlaced and Adam7 interlaced (`interlace = 1`) decoding
//! - All 5 PNG scanline filter types (`None`, `Sub`, `Up`, `Average`, `Paeth`)
//! - `tEXt` metadata chunks and `pHYs` pixel-density chunks
//! - Strict chunk ordering (`IHDR` first, `PLTE` before `IDAT` for indexed, `IEND` last) and CRC-32 validation

use std::path::Path;

use super::hdf5::inflate;

const PNG_MAGIC: &[u8; 8] = b"\x89PNG\r\n\x1a\n";

/// Supported PNG pixel formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PngColorType {
    /// 8-bit grayscale (1 byte/pixel).
    Grayscale8,
    /// 16-bit big-endian grayscale (2 bytes/pixel).
    Grayscale16,
    /// 8-bit grayscale + alpha (2 bytes/pixel).
    GrayscaleAlpha8,
    /// 16-bit big-endian grayscale + alpha (4 bytes/pixel).
    GrayscaleAlpha16,
    /// 24-bit RGB (3 bytes/pixel).
    Rgb8,
    /// 48-bit big-endian RGB (6 bytes/pixel).
    Rgb16,
    /// 32-bit RGBA (4 bytes/pixel).
    Rgba8,
    /// 64-bit big-endian RGBA (8 bytes/pixel).
    Rgba16,
    /// 8-bit palette-indexed (1 byte/pixel, expanded to RGBA8 on decode when `palette` is present).
    Indexed8,
}

impl PngColorType {
    const fn bytes_per_pixel(self) -> usize {
        match self {
            | Self::Grayscale8 | Self::Indexed8 => 1,
            | Self::Grayscale16 | Self::GrayscaleAlpha8 => 2,
            | Self::Rgb8 => 3,
            | Self::GrayscaleAlpha16 | Self::Rgba8 => 4,
            | Self::Rgb16 => 6,
            | Self::Rgba16 => 8,
        }
    }

    const fn bit_depth_and_type(self) -> (u8, u8) {
        match self {
            | Self::Grayscale8 => (8, 0),
            | Self::Grayscale16 => (16, 0),
            | Self::GrayscaleAlpha8 => (8, 4),
            | Self::GrayscaleAlpha16 => (16, 4),
            | Self::Rgb8 => (8, 2),
            | Self::Rgb16 => (16, 2),
            | Self::Rgba8 => (8, 6),
            | Self::Rgba16 => (16, 6),
            | Self::Indexed8 => (8, 3),
        }
    }
}

/// A decoded or ready-to-encode PNG image.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PngImage {
    /// Image width in pixels.
    pub width: u32,
    /// Image height in pixels.
    pub height: u32,
    /// Pixel format.
    pub color_type: PngColorType,
    /// Row-major pixel buffer.
    pub data: Vec<u8>,
}

fn crc32(data: &[u8]) -> u32 {
    let mut crc = 0xffff_ffff_u32;
    for &b in data {
        crc ^= u32::from(b);
        for _ in 0..8 {
            let odd = crc & 1;
            crc >>= 1;
            if odd != 0 {
                crc ^= 0xedb8_8320;
            }
        }
    }
    !crc
}

fn push_chunk(
    out: &mut Vec<u8>,
    tag: &[u8; 4],
    payload: &[u8],
) {
    out.extend_from_slice(&(payload.len() as u32).to_be_bytes());
    let crc_start = out.len();
    out.extend_from_slice(tag);
    out.extend_from_slice(payload);
    let cksum = crc32(&out[crc_start..]);
    out.extend_from_slice(&cksum.to_be_bytes());
}

/// Encodes a [`PngImage`] into PNG bytes (with optional `tEXt` metadata).
///
/// # Errors
/// Fails when dimensions are zero, `color_type` is `Indexed8` (use `Rgb8`/`Rgba8` or
/// `encode_png_with_metadata`), or `data.len()` does not match `width * height * bpp`.
pub fn encode_png(img: &PngImage) -> Result<Vec<u8>, String> {
    encode_png_with_metadata(img, &[])
}

/// Encodes a [`PngImage`] with `tEXt` key/value metadata chunks.
///
/// # Errors
/// Fails when dimensions are zero or `data.len()` mismatches.
pub fn encode_png_with_metadata(
    img: &PngImage,
    metadata: &[(&str, &str)],
) -> Result<Vec<u8>, String> {
    if img.width == 0 || img.height == 0 {
        return Err("PNG dimensions must be positive".to_owned());
    }
    if img.color_type == PngColorType::Indexed8 {
        return Err("encode_png requires non-indexed color_type (use Rgb8/Rgba8)".to_owned());
    }
    let bpp = img.color_type.bytes_per_pixel();
    let row_bytes = (img.width as usize)
        .checked_mul(bpp)
        .ok_or_else(|| "PNG row size overflows".to_owned())?;
    let expected = row_bytes
        .checked_mul(img.height as usize)
        .ok_or_else(|| "PNG buffer size overflows".to_owned())?;
    if img.data.len() != expected {
        return Err(format!(
            "PNG expected {expected} bytes for {}x{}, got {}",
            img.width,
            img.height,
            img.data.len()
        ));
    }

    let mut filtered = Vec::with_capacity((row_bytes + 1) * (img.height as usize));
    for row in img.data.chunks_exact(row_bytes) {
        filtered.push(1); // Sub filter
        for i in 0..row_bytes {
            let left = if i >= bpp { row[i - bpp] } else { 0 };
            filtered.push(row[i].wrapping_sub(left));
        }
    }
    let idat = inflate::zlib_compress(&filtered);

    let (bit_depth, color_code) = img.color_type.bit_depth_and_type();
    let mut ihdr = Vec::with_capacity(13);
    ihdr.extend_from_slice(&img.width.to_be_bytes());
    ihdr.extend_from_slice(&img.height.to_be_bytes());
    ihdr.push(bit_depth);
    ihdr.push(color_code);
    ihdr.push(0);
    ihdr.push(0);
    ihdr.push(0);

    let mut out = Vec::with_capacity(64 + idat.len());
    out.extend_from_slice(PNG_MAGIC);
    push_chunk(&mut out, b"IHDR", &ihdr);
    for (k, v) in metadata {
        let mut t = Vec::with_capacity(k.len() + 1 + v.len());
        t.extend_from_slice(k.as_bytes());
        t.push(0);
        t.extend_from_slice(v.as_bytes());
        push_chunk(&mut out, b"tEXt", &t);
    }
    push_chunk(&mut out, b"IDAT", &idat);
    push_chunk(&mut out, b"IEND", &[]);
    Ok(out)
}

/// Encodes `img` and writes it to `path`.
///
/// # Errors
/// As [`encode_png`], plus file I/O errors.
pub fn write_png_file<P: AsRef<Path>>(
    path: P,
    img: &PngImage,
) -> Result<(), String> {
    let bytes = encode_png(img)?;
    std::fs::write(path.as_ref(), &bytes)
        .map_err(|e| format!("cannot write {}: {e}", path.as_ref().display()))
}

fn paeth(
    a: u8,
    b: u8,
    c: u8,
) -> u8 {
    let a = i32::from(a);
    let b = i32::from(b);
    let c = i32::from(c);
    let p = a + b - c;
    let pa = (p - a).abs();
    let pb = (p - b).abs();
    let pc = (p - c).abs();
    if pa <= pb && pa <= pc {
        a as u8
    } else if pb <= pc {
        b as u8
    } else {
        c as u8
    }
}

fn unfilter_pass(
    raw: &[u8],
    width: usize,
    height: usize,
    bit_depth: u8,
    channels: usize,
) -> Result<(Vec<u8>, usize), String> {
    if width == 0 || height == 0 {
        return Ok((Vec::new(), 0));
    }
    let bits_per_pixel = usize::from(bit_depth) * channels;
    let bpp = bits_per_pixel.div_ceil(8);
    let row_bytes = (width * bits_per_pixel).div_ceil(8);
    let needed = (row_bytes + 1) * height;
    if raw.len() < needed {
        return Err(format!(
            "PNG scanline data truncated: got {} bytes, need {needed}",
            raw.len()
        ));
    }
    let mut out = vec![0_u8; row_bytes * height];
    for y in 0..height {
        let src_row = &raw[y * (row_bytes + 1)..(y + 1) * (row_bytes + 1)];
        let filter = src_row[0];
        let encoded = &src_row[1..];
        let dst_start = y * row_bytes;
        for x in 0..row_bytes {
            let a = if x >= bpp {
                out[dst_start + x - bpp]
            } else {
                0
            };
            let b = if y > 0 {
                out[(y - 1) * row_bytes + x]
            } else {
                0
            };
            let c = if y > 0 && x >= bpp {
                out[(y - 1) * row_bytes + x - bpp]
            } else {
                0
            };
            let recon = match filter {
                | 0 => encoded[x],
                | 1 => encoded[x].wrapping_add(a),
                | 2 => encoded[x].wrapping_add(b),
                | 3 => encoded[x].wrapping_add(((u16::from(a) + u16::from(b)) / 2) as u8),
                | 4 => encoded[x].wrapping_add(paeth(a, b, c)),
                | other => return Err(format!("unsupported PNG filter {other}")),
            };
            out[dst_start + x] = recon;
        }
    }
    Ok((out, needed))
}

/// Decodes PNG bytes into a [`PngImage`] (expanding `Indexed` palettes to `Rgba8`/`Rgb8`
/// and de-interlacing Adam7 images).
///
/// # Errors
/// Fails when the stream is malformed, chunk order is invalid, or CRC/zlib checks fail.
pub fn decode_png(bytes: &[u8]) -> Result<PngImage, String> {
    if bytes.get(..8) != Some(PNG_MAGIC.as_slice()) {
        return Err("not a PNG file (bad signature)".to_owned());
    }
    let mut pos = 8_usize;
    let mut width = 0_u32;
    let mut height = 0_u32;
    let mut bit_depth = 0_u8;
    let mut ct_code = 0_u8;
    let mut interlace = 0_u8;
    let mut palette: Option<Vec<[u8; 3]>> = None;
    let mut trns: Option<Vec<u8>> = None;
    let mut idat = Vec::new();
    let mut saw_ihdr = false;
    let mut saw_idat = false;
    let mut saw_iend = false;

    while pos + 12 <= bytes.len() {
        let len = u32::from_be_bytes([bytes[pos], bytes[pos + 1], bytes[pos + 2], bytes[pos + 3]])
            as usize;
        let tag = &bytes[pos + 4..pos + 8];
        let payload_start = pos + 8;
        let payload_end = payload_start
            .checked_add(len)
            .ok_or_else(|| "PNG chunk length overflow".to_owned())?;
        let payload = bytes
            .get(payload_start..payload_end)
            .ok_or_else(|| "PNG chunk truncated".to_owned())?;
        let crc_bytes = bytes
            .get(payload_end..payload_end + 4)
            .ok_or_else(|| "PNG chunk CRC truncated".to_owned())?;
        let expected_crc =
            u32::from_be_bytes([crc_bytes[0], crc_bytes[1], crc_bytes[2], crc_bytes[3]]);
        if crc32(&bytes[pos + 4..payload_end]) != expected_crc {
            return Err("PNG chunk CRC mismatch".to_owned());
        }
        pos = payload_end + 4;

        if !saw_ihdr && tag != b"IHDR" {
            return Err("PNG must start with IHDR chunk".to_owned());
        }
        match tag {
            | b"IHDR" => {
                if saw_ihdr || payload.len() < 13 {
                    return Err("invalid or duplicate PNG IHDR".to_owned());
                }
                width = u32::from_be_bytes([payload[0], payload[1], payload[2], payload[3]]);
                height = u32::from_be_bytes([payload[4], payload[5], payload[6], payload[7]]);
                bit_depth = payload[8];
                ct_code = payload[9];
                interlace = payload[12];
                if width == 0 || height == 0 || interlace > 1 {
                    return Err("invalid PNG dimensions or interlace method".to_owned());
                }
                saw_ihdr = true;
            },
            | b"PLTE" => {
                if saw_idat || payload.len() % 3 != 0 || payload.is_empty() {
                    return Err("invalid PNG PLTE chunk".to_owned());
                }
                palette = Some(
                    payload
                        .chunks_exact(3)
                        .map(|c| [c[0], c[1], c[2]])
                        .collect(),
                );
            },
            | b"tRNS" => {
                if saw_idat {
                    return Err("PNG tRNS after IDAT".to_owned());
                }
                trns = Some(payload.to_vec());
            },
            | b"IDAT" => {
                saw_idat = true;
                idat.extend_from_slice(payload);
            },
            | b"IEND" => {
                saw_iend = true;
                break;
            },
            | _ => {},
        }
    }

    if !saw_ihdr || !saw_idat || !saw_iend {
        return Err("PNG missing required IHDR, IDAT, or IEND chunk".to_owned());
    }
    let channels: usize = match ct_code {
        | 0 | 3 => 1,
        | 2 => 3,
        | 4 => 2,
        | 6 => 4,
        | _ => return Err(format!("unsupported PNG color_type {ct_code}")),
    };
    let raw = inflate::zlib_decompress(&idat)?;

    let w = width as usize;
    let h = height as usize;
    let pixel_bytes = if interlace == 0 {
        let (buf, consumed) = unfilter_pass(&raw, w, h, bit_depth, channels)?;
        if consumed != raw.len() {
            return Err("PNG decompressed size mismatch".to_owned());
        }
        buf
    } else {
        // Adam7 de-interlacing
        const X0: [usize; 7] = [0, 4, 0, 2, 0, 1, 0];
        const Y0: [usize; 7] = [0, 0, 4, 0, 2, 0, 1];
        const DX: [usize; 7] = [8, 8, 4, 4, 2, 2, 1];
        const DY: [usize; 7] = [8, 8, 8, 4, 4, 2, 2];
        let bpp = (usize::from(bit_depth) * channels).div_ceil(8);
        let mut full = vec![0_u8; w * h * bpp];
        let mut offset = 0_usize;
        for pass in 0..7 {
            let pw = if w > X0[pass] {
                (w - X0[pass]).div_ceil(DX[pass])
            } else {
                0
            };
            let ph = if h > Y0[pass] {
                (h - Y0[pass]).div_ceil(DY[pass])
            } else {
                0
            };
            if pw == 0 || ph == 0 {
                continue;
            }
            let (pass_buf, consumed) = unfilter_pass(&raw[offset..], pw, ph, bit_depth, channels)?;
            offset += consumed;
            for py in 0..ph {
                for px in 0..pw {
                    let gx = X0[pass] + px * DX[pass];
                    let gy = Y0[pass] + py * DY[pass];
                    let src = (py * pw + px) * bpp;
                    let dst = (gy * w + gx) * bpp;
                    full[dst..dst + bpp].copy_from_slice(&pass_buf[src..src + bpp]);
                }
            }
        }
        full
    };

    let (color_type, final_data) = match (bit_depth, ct_code) {
        | (8, 0) => (PngColorType::Grayscale8, pixel_bytes),
        | (16, 0) => (PngColorType::Grayscale16, pixel_bytes),
        | (8, 4) => (PngColorType::GrayscaleAlpha8, pixel_bytes),
        | (16, 4) => (PngColorType::GrayscaleAlpha16, pixel_bytes),
        | (8, 2) => (PngColorType::Rgb8, pixel_bytes),
        | (16, 2) => (PngColorType::Rgb16, pixel_bytes),
        | (8, 6) => (PngColorType::Rgba8, pixel_bytes),
        | (16, 6) => (PngColorType::Rgba16, pixel_bytes),
        | (8, 3) => {
            let pal = palette.ok_or("indexed PNG missing PLTE chunk")?;
            if let Some(alpha) = trns {
                let mut rgba = Vec::with_capacity(w * h * 4);
                for &idx in &pixel_bytes {
                    let rgb = pal.get(usize::from(idx)).copied().unwrap_or([0, 0, 0]);
                    let a = alpha.get(usize::from(idx)).copied().unwrap_or(255);
                    rgba.extend_from_slice(&[rgb[0], rgb[1], rgb[2], a]);
                }
                (PngColorType::Rgba8, rgba)
            } else {
                let mut rgb_out = Vec::with_capacity(w * h * 3);
                for &idx in &pixel_bytes {
                    let rgb = pal.get(usize::from(idx)).copied().unwrap_or([0, 0, 0]);
                    rgb_out.extend_from_slice(&rgb);
                }
                (PngColorType::Rgb8, rgb_out)
            }
        },
        | _ => {
            return Err(format!(
                "unsupported PNG format: bit_depth={bit_depth}, color_type={ct_code}"
            ));
        },
    };

    Ok(PngImage {
        width,
        height,
        color_type,
        data: final_data,
    })
}

/// Reads and decodes a PNG file from `path`.
///
/// # Errors
/// Fails when the file cannot be read or decoded.
pub fn read_png_file<P: AsRef<Path>>(path: P) -> Result<PngImage, String> {
    let bytes = std::fs::read(path.as_ref())
        .map_err(|e| format!("cannot read {}: {e}", path.as_ref().display()))?;
    decode_png(&bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn png_round_trip_all_supported_color_types() {
        for (ct, bpp) in [
            (PngColorType::Grayscale8, 1),
            (PngColorType::Grayscale16, 2),
            (PngColorType::GrayscaleAlpha8, 2),
            (PngColorType::GrayscaleAlpha16, 4),
            (PngColorType::Rgb8, 3),
            (PngColorType::Rgb16, 6),
            (PngColorType::Rgba8, 4),
            (PngColorType::Rgba16, 8),
        ] {
            let img = PngImage {
                width: 4,
                height: 3,
                color_type: ct,
                data: (0..(12 * bpp)).map(|i| (i % 251) as u8).collect(),
            };
            let bytes = encode_png_with_metadata(&img, &[("Software", "RSSN")]).unwrap();
            assert_eq!(&bytes[..8], PNG_MAGIC);
            let back = decode_png(&bytes).unwrap();
            assert_eq!(back, img);
        }
    }
}
