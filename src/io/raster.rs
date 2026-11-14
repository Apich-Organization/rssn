//! Minimal pure-Rust raster export: binary PGM (P5) and PPM (P6).
//!
//! The netpbm formats are written directly, byte for byte, so no image
//! encoder dependency is needed. Rows are stored top to bottom.
//!
//! # Limitations
//!
//! - Only binary PGM (P5, 8- or 16-bit) and PPM (P6, 8-bit): no ASCII P2/P3
//!   variants, no 16-bit PPM, no alpha channel, no comments or metadata, no
//!   colour management.
//! - No PNG, JPEG, GIF or any compressed format: those need an encoder, i.e. a
//!   new dependency, so they stay a documented gap.
//! - The API is raw samples with explicit dimensions: there is no resampling,
//!   scaling, text or annotation. Use [`super::plot`] for annotated vector
//!   output and [`super::plot::plot_scalar_field_pgm`] /
//!   [`super::plot::plot_scalar_field_ppm`] to rasterise a sampled scalar
//!   field.
//! - Invalid dimensions or sample counts are rejected before anything is
//!   written, so a failed call leaves no partial file behind.

use std::path::Path;

/// Rejects empty or inconsistent raster shapes before a file is touched.
fn check_image(
    width: usize,
    height: usize,
    samples: usize,
    channels: usize,
) -> Result<(), String> {
    if width == 0 || height == 0 {
        return Err(format!(
            "raster dimensions must be positive, got {width} x {height}"
        ));
    }
    let expected = width
        .checked_mul(height)
        .and_then(|pixels| pixels.checked_mul(channels))
        .ok_or_else(|| "raster dimensions overflow".to_owned())?;
    if samples != expected {
        return Err(format!(
            "expected {expected} samples for a {width} x {height} raster with {channels} \
             channel(s), got {samples}"
        ));
    }
    Ok(())
}

/// Writes `bytes` to `path`, mapping I/O errors to a message naming the file.
fn write_bytes(
    path: &Path,
    bytes: &[u8],
) -> Result<(), String> {
    std::fs::write(path, bytes).map_err(|e| format!("cannot write {}: {e}", path.display()))
}

/// Writes an 8-bit grayscale PGM (P5) file.
///
/// `samples` is row-major, one byte per pixel, top row first. The header is
/// `P5`, the dimensions and the maximum value `255`.
///
/// # Errors
/// Fails for a zero dimension, a sample count other than `width * height`, or
/// when the file cannot be written.
pub fn write_pgm_file<P: AsRef<Path>>(
    path: P,
    width: usize,
    height: usize,
    samples: &[u8],
) -> Result<(), String> {
    check_image(width, height, samples.len(), 1)?;
    let mut out = Vec::with_capacity(samples.len() + 32);
    out.extend_from_slice(format!("P5\n{width} {height}\n255\n").as_bytes());
    out.extend_from_slice(samples);
    write_bytes(path.as_ref(), &out)
}

/// Writes a 16-bit grayscale PGM (P5) file.
///
/// `samples` is row-major, one `u16` per pixel, top row first; samples are
/// stored big-endian as the format requires. The header maximum is `65535`.
///
/// # Errors
/// As [`write_pgm_file`].
pub fn write_pgm16_file<P: AsRef<Path>>(
    path: P,
    width: usize,
    height: usize,
    samples: &[u16],
) -> Result<(), String> {
    check_image(width, height, samples.len(), 1)?;
    let mut out = Vec::with_capacity(samples.len() * 2 + 32);
    out.extend_from_slice(format!("P5\n{width} {height}\n65535\n").as_bytes());
    for &sample in samples {
        out.extend_from_slice(&sample.to_be_bytes());
    }
    write_bytes(path.as_ref(), &out)
}

/// Writes an 8-bit RGB PPM (P6) file.
///
/// `samples` is row-major, three bytes (red, green, blue) per pixel, top row
/// first. The header is `P6`, the dimensions and the maximum value `255`.
///
/// # Errors
/// Fails for a zero dimension, a sample count other than
/// `width * height * 3`, or when the file cannot be written.
pub fn write_ppm_file<P: AsRef<Path>>(
    path: P,
    width: usize,
    height: usize,
    samples: &[u8],
) -> Result<(), String> {
    check_image(width, height, samples.len(), 3)?;
    let mut out = Vec::with_capacity(samples.len() + 32);
    out.extend_from_slice(format!("P6\n{width} {height}\n255\n").as_bytes());
    out.extend_from_slice(samples);
    write_bytes(path.as_ref(), &out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp(
        tag: &str,
        ext: &str,
    ) -> std::path::PathBuf {
        std::env::temp_dir().join(format!("rssn_raster_{}_{tag}.{ext}", std::process::id()))
    }

    #[test]
    fn pgm_header_and_samples_round_trip() {
        let p = tmp("pgm", "pgm");
        write_pgm_file(&p, 2, 2, &[0, 85, 170, 255]).unwrap();
        let bytes = std::fs::read(&p).unwrap();
        assert!(bytes.starts_with(b"P5\n2 2\n255\n"));
        assert_eq!(&bytes[11..], &[0, 85, 170, 255]);
        let _r = std::fs::remove_file(&p);
    }

    #[test]
    fn invalid_shapes_are_rejected_without_a_file() {
        let p = tmp("bad", "pgm");
        assert!(write_pgm_file(&p, 2, 2, &[1, 2, 3]).is_err());
        assert!(write_pgm_file(&p, 0, 2, &[]).is_err());
        assert!(write_ppm_file(&p, 1, 1, &[1, 2]).is_err());
        assert!(write_pgm16_file(&p, 1, 1, &[1, 2]).is_err());
        assert!(!p.exists());
    }
}
