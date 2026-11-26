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
//! - No PNG, JPEG or other compressed formats beyond GIF: GIF animation I/O
//!   is provided in [`super::gif`] and re-exported here.
//! - The API is raw samples with explicit dimensions: there is no resampling,
//!   scaling, text or annotation. Use [`super::plot`] for annotated vector
//!   output and [`super::plot::plot_scalar_field_pgm`] /
//!   [`super::plot::plot_scalar_field_ppm`] to rasterise a sampled scalar
//!   field.
//! - Invalid dimensions or sample counts are rejected before opening the
//!   destination. Errors from opening, writing and flushing are reported.

use std::fmt;
use std::fs::File;
use std::io;
use std::io::BufWriter;
use std::io::Write;
use std::path::Path;

pub use super::gif::GifAnimation;
pub use super::gif::GifFrame;
pub use super::gif::decode_gif;
pub use super::gif::encode_gif;
pub use super::gif::read_gif_file;
pub use super::gif::write_gif_file;
pub use super::png::PngColorType;
pub use super::png::PngImage;
pub use super::png::decode_png;
pub use super::png::encode_png;
pub use super::png::read_png_file;
pub use super::png::write_png_file;

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

fn write_image(
    path: &Path,
    header: fmt::Arguments<'_>,
    sample_bytes: usize,
    write_samples: impl FnOnce(&mut BufWriter<File>) -> io::Result<()>,
) -> Result<(), String> {
    File::create(path)
        .and_then(|file| {
            // Keep small images in one write while limiting the extra memory
            // for larger inputs. A long byte slice bypasses the buffer.
            let mut output =
                BufWriter::with_capacity(sample_bytes.min(8192).saturating_add(64), file);
            output.write_fmt(header)?;
            write_samples(&mut output)?;
            output.flush()
        })
        .map_err(|e| format!("cannot write {}: {e}", path.display()))
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
    write_image(
        path.as_ref(),
        format_args!("P5\n{width} {height}\n255\n"),
        samples.len(),
        |output| output.write_all(samples),
    )
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
    write_image(
        path.as_ref(),
        format_args!("P5\n{width} {height}\n65535\n"),
        samples.len().saturating_mul(2),
        |output| {
            use crate::kernels::reuse::ScratchPool;

            const CHUNK_SAMPLES: usize = 524_288;
            static WORK: ScratchPool<[u8; 2]> = ScratchPool::new(2, CHUNK_SAMPLES * 2);
            let count = samples.len().min(CHUNK_SAMPLES);
            let mut buffer = WORK.checkout_for(count);
            buffer.resize(count, [0; 2]);
            for chunk in samples.chunks(buffer.len()) {
                #[allow(clippy::indexing_slicing, reason = "each chunk fits the buffer")]
                let converted = &mut buffer[..chunk.len()];
                for (bytes, sample) in converted.iter_mut().zip(chunk) {
                    *bytes = sample.to_be_bytes();
                }
                output.write_all(converted.as_flattened())?;
            }
            Ok(())
        },
    )
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
    write_image(
        path.as_ref(),
        format_args!("P6\n{width} {height}\n255\n"),
        samples.len(),
        |output| output.write_all(samples),
    )
}

/// Writes a 16-bit RGB PPM (P6) file (big-endian `u16` per channel).
///
/// # Errors
/// Fails for a zero dimension, a sample count other than `width * height * 3`,
/// or when the file cannot be written.
pub fn write_ppm16_file<P: AsRef<Path>>(
    path: P,
    width: usize,
    height: usize,
    samples: &[u16],
) -> Result<(), String> {
    check_image(width, height, samples.len(), 3)?;
    write_image(
        path.as_ref(),
        format_args!("P6\n{width} {height}\n65535\n"),
        samples.len().saturating_mul(2),
        |output| {
            let mut buf = Vec::with_capacity(samples.len().min(65_536) * 2);
            for chunk in samples.chunks(65_536) {
                buf.clear();
                for &s in chunk {
                    buf.extend_from_slice(&s.to_be_bytes());
                }
                output.write_all(&buf)?;
            }
            Ok(())
        },
    )
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
