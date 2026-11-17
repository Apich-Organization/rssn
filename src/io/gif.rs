//! GIF animation encoding and decoding.
//!
//! Provides [`GifFrame`], [`GifAnimation`], [`encode_gif`], [`write_gif_file`],
//! [`decode_gif`], and [`read_gif_file`] backed by the pure-Rust `gif` crate.

use std::io::Cursor;
use std::path::Path;

/// A single RGBA frame in a GIF animation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GifFrame {
    /// Frame width in pixels.
    pub width: u16,
    /// Frame height in pixels.
    pub height: u16,
    /// Frame display delay in centiseconds (`1/100` of a second).
    pub delay_cs: u16,
    /// Row-major RGBA pixel buffer (`width * height * 4` bytes).
    pub rgba: Vec<u8>,
}

/// An animated GIF image composed of one or more [`GifFrame`] values.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GifAnimation {
    /// Logical screen width in pixels.
    pub width: u16,
    /// Logical screen height in pixels.
    pub height: u16,
    /// Number of animation loops (`0` means infinite looping).
    pub loop_count: u16,
    /// Ordered sequence of frames.
    pub frames: Vec<GifFrame>,
}

/// Encodes a [`GifAnimation`] into GIF bytes.
///
/// # Errors
/// Fails when dimensions are zero, `frames` is empty, any frame has invalid
/// dimensions or an `rgba` buffer length other than `width * height * 4`, or
/// the GIF encoder returns an error.
pub fn encode_gif(anim: &GifAnimation) -> Result<Vec<u8>, String> {
    if anim.width == 0 || anim.height == 0 {
        return Err(format!(
            "GIF dimensions must be positive, got {} x {}",
            anim.width, anim.height
        ));
    }
    if anim.frames.is_empty() {
        return Err("GIF animation must contain at least one frame".to_owned());
    }

    for (idx, frame) in anim.frames.iter().enumerate() {
        if frame.width == 0 || frame.height == 0 {
            return Err(format!(
                "GIF frame {idx} dimensions must be positive, got {} x {}",
                frame.width, frame.height
            ));
        }
        let expected = usize::from(frame.width)
            .checked_mul(usize::from(frame.height))
            .and_then(|pixels| pixels.checked_mul(4))
            .ok_or_else(|| "GIF frame dimensions overflow".to_owned())?;
        if frame.rgba.len() != expected {
            return Err(format!(
                "GIF frame {idx} expected {expected} RGBA bytes for {} x {}, got {}",
                frame.width,
                frame.height,
                frame.rgba.len()
            ));
        }
    }

    let mut out = Vec::new();
    {
        let mut encoder = gif::Encoder::new(&mut out, anim.width, anim.height, &[])
            .map_err(|e| format!("failed to initialize GIF encoder: {e}"))?;
        let repeat = if anim.loop_count == 0 {
            gif::Repeat::Infinite
        } else {
            gif::Repeat::Finite(anim.loop_count)
        };
        encoder
            .set_repeat(repeat)
            .map_err(|e| format!("failed to set GIF repeat count: {e}"))?;

        for frame in &anim.frames {
            let mut pixels = frame.rgba.clone();
            let mut gif_frame =
                gif::Frame::from_rgba_speed(frame.width, frame.height, &mut pixels, 10);
            gif_frame.delay = frame.delay_cs;
            encoder
                .write_frame(&gif_frame)
                .map_err(|e| format!("failed to encode GIF frame: {e}"))?;
        }
    }

    Ok(out)
}

/// Encodes `anim` and writes the resulting GIF bytes to `path`.
///
/// # Errors
/// As [`encode_gif`], plus file-writing errors.
pub fn write_gif_file(
    path: &Path,
    anim: &GifAnimation,
) -> Result<(), String> {
    let bytes = encode_gif(anim)?;
    std::fs::write(path, &bytes).map_err(|e| format!("cannot write {}: {e}", path.display()))
}

/// Decodes GIF bytes into a [`GifAnimation`] with RGBA frames.
///
/// # Errors
/// Fails when `bytes` is not a valid GIF stream or contains no frames.
pub fn decode_gif(bytes: &[u8]) -> Result<GifAnimation, String> {
    let mut options = gif::DecodeOptions::new();
    options.set_color_output(gif::ColorOutput::RGBA);
    let mut decoder = options
        .read_info(Cursor::new(bytes))
        .map_err(|e| format!("failed to decode GIF header: {e}"))?;

    let width = decoder.width();
    let height = decoder.height();
    if width == 0 || height == 0 {
        return Err(format!(
            "decoded GIF dimensions must be positive, got {width} x {height}"
        ));
    }

    let loop_count = match decoder.repeat() {
        | gif::Repeat::Infinite => 0,
        | gif::Repeat::Finite(n) => n,
    };

    let mut frames = Vec::new();
    while let Some(frame) = decoder
        .read_next_frame()
        .map_err(|e| format!("failed to decode GIF frame: {e}"))?
    {
        frames.push(GifFrame {
            width: frame.width,
            height: frame.height,
            delay_cs: frame.delay,
            rgba: frame.buffer.to_vec(),
        });
    }

    if frames.is_empty() {
        return Err("GIF stream contains no frames".to_owned());
    }

    Ok(GifAnimation {
        width,
        height,
        loop_count,
        frames,
    })
}

/// Reads a GIF file from `path` and decodes it into a [`GifAnimation`].
///
/// # Errors
/// Fails when the file cannot be read or decoded as a GIF.
pub fn read_gif_file(path: &Path) -> Result<GifAnimation, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    decode_gif(&bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encode_and_decode_gif_round_trip() {
        let frame1 = GifFrame {
            width: 2,
            height: 2,
            delay_cs: 10,
            rgba: vec![
                255, 0, 0, 255, // red
                0, 255, 0, 255, // green
                0, 0, 255, 255, // blue
                255, 255, 255, 255, // white
            ],
        };
        let frame2 = GifFrame {
            width: 2,
            height: 2,
            delay_cs: 20,
            rgba: vec![
                0, 0, 0, 255, // black
                255, 255, 0, 255, // yellow
                0, 255, 255, 255, // cyan
                255, 0, 255, 255, // magenta
            ],
        };
        let anim = GifAnimation {
            width: 2,
            height: 2,
            loop_count: 3,
            frames: vec![frame1, frame2],
        };
        let bytes = encode_gif(&anim).unwrap();
        assert!(bytes.starts_with(b"GIF89a") || bytes.starts_with(b"GIF87a"));

        let decoded = decode_gif(&bytes).unwrap();
        assert_eq!(decoded.width, 2);
        assert_eq!(decoded.height, 2);
        assert_eq!(decoded.frames.len(), 2);
        assert_eq!(decoded.frames.first().map(|f| f.delay_cs), Some(10));
        assert_eq!(decoded.frames.get(1).map(|f| f.delay_cs), Some(20));
    }

    #[test]
    fn write_and_read_gif_file_round_trip() {
        let path = std::env::temp_dir().join(format!("rssn_gif_test_{}.gif", std::process::id()));
        let anim = GifAnimation {
            width: 1,
            height: 1,
            loop_count: 0,
            frames: vec![GifFrame {
                width: 1,
                height: 1,
                delay_cs: 5,
                rgba: vec![12, 34, 56, 255],
            }],
        };
        write_gif_file(&path, &anim).unwrap();
        let decoded = read_gif_file(&path).unwrap();
        let _r = std::fs::remove_file(&path);
        assert_eq!(decoded.width, 1);
        assert_eq!(decoded.height, 1);
        assert_eq!(decoded.frames.len(), 1);
    }

    #[test]
    fn invalid_gif_inputs_are_rejected() {
        let empty = GifAnimation {
            width: 2,
            height: 2,
            loop_count: 0,
            frames: vec![],
        };
        assert!(encode_gif(&empty).is_err());

        let bad_buf = GifAnimation {
            width: 2,
            height: 2,
            loop_count: 0,
            frames: vec![GifFrame {
                width: 2,
                height: 2,
                delay_cs: 10,
                rgba: vec![0; 4],
            }],
        };
        assert!(encode_gif(&bad_buf).is_err());
    }
}
