//! GIF animation encoding and decoding.
//!
//! Provides [`GifFrame`], [`GifAnimation`], [`encode_gif`], [`write_gif_file`],
//! [`decode_gif`], and [`read_gif_file`] backed by the pure-Rust `gif` crate.

use std::io::Cursor;
use std::path::Path;

/// GIF89a frame disposal method.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum GifDisposal {
    /// No disposal specified.
    #[default]
    Any,
    /// Leave frame in place.
    Keep,
    /// Restore to background (transparent) before next frame.
    Background,
    /// Restore to previous canvas before this frame.
    Previous,
}

impl From<GifDisposal> for gif::DisposalMethod {
    fn from(d: GifDisposal) -> Self {
        match d {
            | GifDisposal::Any => Self::Any,
            | GifDisposal::Keep => Self::Keep,
            | GifDisposal::Background => Self::Background,
            | GifDisposal::Previous => Self::Previous,
        }
    }
}

impl From<gif::DisposalMethod> for GifDisposal {
    fn from(d: gif::DisposalMethod) -> Self {
        match d {
            | gif::DisposalMethod::Any => Self::Any,
            | gif::DisposalMethod::Keep => Self::Keep,
            | gif::DisposalMethod::Background => Self::Background,
            | gif::DisposalMethod::Previous => Self::Previous,
        }
    }
}

/// A single RGBA frame in a GIF animation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GifFrame {
    /// Frame width in pixels.
    pub width: u16,
    /// Frame height in pixels.
    pub height: u16,
    /// X offset on the logical screen.
    pub left: u16,
    /// Y offset on the logical screen.
    pub top: u16,
    /// Frame display delay in centiseconds (`1/100` of a second).
    pub delay_cs: u16,
    /// Disposal method after displaying this frame.
    pub dispose: GifDisposal,
    /// Row-major RGBA pixel buffer (`width * height * 4` bytes).
    pub rgba: Vec<u8>,
}

impl GifFrame {
    /// Creates a full-origin (`left = 0, top = 0, dispose = Any`) RGBA frame.
    #[must_use]
    pub fn new(
        width: u16,
        height: u16,
        delay_cs: u16,
        rgba: Vec<u8>,
    ) -> Self {
        Self {
            width,
            height,
            left: 0,
            top: 0,
            delay_cs,
            dispose: GifDisposal::Any,
            rgba,
        }
    }
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
    /// Optional GIF89a comment extension text.
    pub comment: Option<String>,
    /// Ordered sequence of frames.
    pub frames: Vec<GifFrame>,
}

impl GifAnimation {
    /// Creates a [`GifAnimation`] without a comment extension.
    #[must_use]
    pub fn new(
        width: u16,
        height: u16,
        loop_count: u16,
        frames: Vec<GifFrame>,
    ) -> Self {
        Self {
            width,
            height,
            loop_count,
            comment: None,
            frames,
        }
    }
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
        if u32::from(frame.left) + u32::from(frame.width) > u32::from(anim.width)
            || u32::from(frame.top) + u32::from(frame.height) > u32::from(anim.height)
        {
            return Err(format!(
                "GIF frame {idx} ({}+{}, {}+{}) exceeds logical screen {}x{}",
                frame.left, frame.width, frame.top, frame.height, anim.width, anim.height
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

        if let Some(comment) = &anim.comment {
            let _ = comment; // Written before trailer below as GIF89a Comment Extension
        }

        for frame in &anim.frames {
            let mut pixels = frame.rgba.clone();
            let mut gif_frame =
                gif::Frame::from_rgba_speed(frame.width, frame.height, &mut pixels, 10);
            gif_frame.left = frame.left;
            gif_frame.top = frame.top;
            gif_frame.delay = frame.delay_cs;
            gif_frame.dispose = frame.dispose.into();
            encoder
                .write_frame(&gif_frame)
                .map_err(|e| format!("failed to encode GIF frame: {e}"))?;
        }
    }

    if let Some(comment) = &anim.comment {
        if !comment.is_empty() && out.last().copied() == Some(0x3B) {
            out.pop();
            out.push(0x21);
            out.push(0xFE);
            for chunk in comment.as_bytes().chunks(255) {
                out.push(chunk.len() as u8);
                out.extend_from_slice(chunk);
            }
            out.push(0x00);
            out.push(0x3B);
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
    if bytes.len() < 14 || (!bytes.starts_with(b"GIF89a") && !bytes.starts_with(b"GIF87a")) {
        return Err("invalid GIF signature (expected GIF89a or GIF87a)".to_owned());
    }
    if bytes.last().copied() != Some(0x3B) {
        return Err("GIF stream missing 0x3B trailer".to_owned());
    }
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
            left: frame.left,
            top: frame.top,
            delay_cs: frame.delay,
            dispose: frame.dispose.into(),
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
        comment: None,
        frames,
    })
}

/// Decodes GIF bytes and composites every sub-frame (`left`, `top`, `dispose`)
/// onto the full `width x height` logical screen canvas.
///
/// # Errors
/// As [`decode_gif`].
pub fn decode_gif_composited(bytes: &[u8]) -> Result<GifAnimation, String> {
    let raw = decode_gif(bytes)?;
    let w = usize::from(raw.width);
    let h = usize::from(raw.height);
    let mut canvas = vec![0_u8; w * h * 4];
    let mut out_frames = Vec::with_capacity(raw.frames.len());

    for f in &raw.frames {
        let prev = if f.dispose == GifDisposal::Previous {
            Some(canvas.clone())
        } else {
            None
        };
        let fw = usize::from(f.width);
        let fh = usize::from(f.height);
        let fx = usize::from(f.left);
        let fy = usize::from(f.top);
        for row in 0..fh {
            for col in 0..fw {
                let gx = fx + col;
                let gy = fy + row;
                if gx < w && gy < h {
                    let s = (row * fw + col) * 4;
                    let d = (gy * w + gx) * 4;
                    if f.rgba[s + 3] > 0 {
                        canvas[d..d + 4].copy_from_slice(&f.rgba[s..s + 4]);
                    }
                }
            }
        }
        out_frames.push(GifFrame::new(
            raw.width,
            raw.height,
            f.delay_cs,
            canvas.clone(),
        ));
        match f.dispose {
            | GifDisposal::Background => {
                for row in 0..fh {
                    for col in 0..fw {
                        let gx = fx + col;
                        let gy = fy + row;
                        if gx < w && gy < h {
                            let d = (gy * w + gx) * 4;
                            canvas[d..d + 4].fill(0);
                        }
                    }
                }
            },
            | GifDisposal::Previous => {
                if let Some(p) = prev {
                    canvas = p;
                }
            },
            | GifDisposal::Any | GifDisposal::Keep => {},
        }
    }

    Ok(GifAnimation {
        width: raw.width,
        height: raw.height,
        loop_count: raw.loop_count,
        comment: raw.comment,
        frames: out_frames,
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

/// Options for rendering scientific 2D fields into animated GIFs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GifPlotOptions {
    /// Colormap used to map scalar values to RGB.
    pub colormap: super::colormap::Colormap,
    /// Integer upscaling factor per grid cell (minimum `1`).
    pub scale: u16,
    /// Delay per frame in centiseconds (`1/100` s).
    pub delay_cs: u16,
    /// Loop count (`0` = infinite).
    pub loop_count: u16,
    /// Optional fixed `(min, max)` color scale across all frames.
    pub fixed_range: Option<(f64, f64)>,
}

impl Default for GifPlotOptions {
    fn default() -> Self {
        Self {
            colormap: super::colormap::Colormap::Viridis,
            scale: 4,
            delay_cs: 8,
            loop_count: 0,
            fixed_range: None,
        }
    }
}

/// Encodes a sequence of 2D scalar fields (`Array2<f64>`) as an animated GIF.
///
/// # Errors
/// Fails when `frames` is empty, any frame is empty, dimensions mismatch, or encoding fails.
pub fn encode_scalar_field_gif(
    frames: &[ndarray::Array2<f64>],
    opts: &GifPlotOptions,
) -> Result<Vec<u8>, String> {
    let first = frames
        .first()
        .ok_or_else(|| "scalar field GIF requires at least one frame".to_owned())?;
    let (rows, cols) = first.dim();
    if rows == 0 || cols == 0 {
        return Err("scalar field GIF frames cannot be empty".to_owned());
    }
    let scale = usize::from(opts.scale.max(1));
    let width = u16::try_from(cols * scale).map_err(|_| "GIF width overflows u16".to_owned())?;
    let height = u16::try_from(rows * scale).map_err(|_| "GIF height overflows u16".to_owned())?;

    let (vmin, vmax) = if let Some((lo, hi)) = opts.fixed_range {
        (lo, hi)
    } else {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for f in frames {
            for &v in f.iter() {
                if v.is_finite() {
                    lo = lo.min(v);
                    hi = hi.max(v);
                }
            }
        }
        if !lo.is_finite() || !hi.is_finite() {
            (0.0, 1.0)
        } else {
            (lo, hi)
        }
    };
    let span = (vmax - vmin).abs();

    let mut gif_frames = Vec::with_capacity(frames.len());
    for (idx, field) in frames.iter().enumerate() {
        if field.dim() != (rows, cols) {
            return Err(format!(
                "GIF frame {idx} shape {:?} does not match first frame {:?}",
                field.dim(),
                (rows, cols)
            ));
        }
        let mut rgba = vec![0_u8; usize::from(width) * usize::from(height) * 4];
        for r in 0..rows {
            for c in 0..cols {
                let v = field[(r, c)];
                let (cr, cg, cb) = if v.is_finite() {
                    let t = if span > 0.0 {
                        (v - vmin) / span
                    } else {
                        0.5
                    };
                    opts.colormap.sample_rgb(t)
                } else {
                    (221, 221, 221)
                };
                for sy in 0..scale {
                    for sx in 0..scale {
                        let px = c * scale + sx;
                        let py = r * scale + sy;
                        let off = (py * usize::from(width) + px) * 4;
                        rgba[off] = cr;
                        rgba[off + 1] = cg;
                        rgba[off + 2] = cb;
                        rgba[off + 3] = 255;
                    }
                }
            }
        }
        gif_frames.push(GifFrame::new(width, height, opts.delay_cs, rgba));
    }

    encode_gif(&GifAnimation::new(
        width,
        height,
        opts.loop_count,
        gif_frames,
    ))
}

/// Encodes a sequence of 2D scalar fields and writes the animated GIF to `path`.
///
/// # Errors
/// As [`encode_scalar_field_gif`], plus file I/O errors.
pub fn write_scalar_field_gif<P: AsRef<Path>>(
    path: P,
    frames: &[ndarray::Array2<f64>],
    opts: &GifPlotOptions,
) -> Result<(), String> {
    let bytes = encode_scalar_field_gif(frames, opts)?;
    std::fs::write(path.as_ref(), &bytes)
        .map_err(|e| format!("cannot write {}: {e}", path.as_ref().display()))
}

fn set_pixel(
    rgba: &mut [u8],
    width: usize,
    height: usize,
    x: i32,
    y: i32,
    color: [u8; 4],
) {
    if x >= 0 && y >= 0 && (x as usize) < width && (y as usize) < height {
        let off = ((y as usize) * width + (x as usize)) * 4;
        rgba[off..off + 4].copy_from_slice(&color);
    }
}

fn draw_line(
    rgba: &mut [u8],
    width: usize,
    height: usize,
    mut x0: i32,
    mut y0: i32,
    x1: i32,
    y1: i32,
    color: [u8; 4],
) {
    let dx = (x1 - x0).abs();
    let dy = -(y1 - y0).abs();
    let sx = if x0 < x1 { 1 } else { -1 };
    let sy = if y0 < y1 { 1 } else { -1 };
    let mut err = dx + dy;
    loop {
        set_pixel(rgba, width, height, x0, y0, color);
        if x0 == x1 && y0 == y1 {
            break;
        }
        let e2 = 2 * err;
        if e2 >= dy {
            err += dy;
            x0 += sx;
        }
        if e2 <= dx {
            err += dx;
            y0 += sy;
        }
    }
}

fn draw_disc(
    rgba: &mut [u8],
    width: usize,
    height: usize,
    cx: i32,
    cy: i32,
    radius: i32,
    color: [u8; 4],
) {
    for dy in -radius..=radius {
        for dx in -radius..=radius {
            if dx * dx + dy * dy <= radius * radius {
                set_pixel(rgba, width, height, cx + dx, cy + dy, color);
            }
        }
    }
}

/// Encodes a time-series of 1D `(x, y)` curves into an animated GIF.
///
/// # Errors
/// Fails when `frames` is empty or dimensions are too small.
pub fn encode_series_animation_gif(
    frames: &[Vec<(f64, f64)>],
    width: u16,
    height: u16,
    delay_cs: u16,
) -> Result<Vec<u8>, String> {
    if frames.is_empty() {
        return Err("series animation GIF requires at least one frame".to_owned());
    }
    if width < 32 || height < 32 {
        return Err("series animation GIF dimensions must be at least 32x32".to_owned());
    }
    let mut xmin = f64::INFINITY;
    let mut xmax = f64::NEG_INFINITY;
    let mut ymin = f64::INFINITY;
    let mut ymax = f64::NEG_INFINITY;
    for f in frames {
        for &(x, y) in f {
            if x.is_finite() && y.is_finite() {
                xmin = xmin.min(x);
                xmax = xmax.max(x);
                ymin = ymin.min(y);
                ymax = ymax.max(y);
            }
        }
    }
    if !xmin.is_finite() {
        xmin = 0.0;
        xmax = 1.0;
        ymin = 0.0;
        ymax = 1.0;
    }
    let xspan = (xmax - xmin).max(1e-12);
    let yspan = (ymax - ymin).max(1e-12);
    let w = usize::from(width);
    let h = usize::from(height);
    let pad = 16_i32;
    let pw = (w as i32) - 2 * pad;
    let ph = (h as i32) - 2 * pad;

    let mut gif_frames = Vec::with_capacity(frames.len());
    for f in frames {
        let mut rgba = vec![255_u8; w * h * 4];
        // Axes box
        let border = [180, 180, 180, 255];
        draw_line(&mut rgba, w, h, pad, pad, pad + pw, pad, border);
        draw_line(&mut rgba, w, h, pad, pad + ph, pad + pw, pad + ph, border);
        draw_line(&mut rgba, w, h, pad, pad, pad, pad + ph, border);
        draw_line(&mut rgba, w, h, pad + pw, pad, pad + pw, pad + ph, border);

        let pts: Vec<(i32, i32)> = f
            .iter()
            .filter(|(x, y)| x.is_finite() && y.is_finite())
            .map(|&(x, y)| {
                let px = pad + (((x - xmin) / xspan) * f64::from(pw)).round() as i32;
                let py = pad + ph - (((y - ymin) / yspan) * f64::from(ph)).round() as i32;
                (px, py)
            })
            .collect();
        for seg in pts.windows(2) {
            draw_line(
                &mut rgba,
                w,
                h,
                seg[0].0,
                seg[0].1,
                seg[1].0,
                seg[1].1,
                [31, 119, 180, 255],
            );
        }
        gif_frames.push(GifFrame::new(width, height, delay_cs, rgba));
    }

    encode_gif(&GifAnimation::new(width, height, 0, gif_frames))
}

/// Writes a 1D curve animation GIF to `path`.
///
/// # Errors
/// As [`encode_series_animation_gif`], plus file I/O errors.
pub fn write_series_animation_gif<P: AsRef<Path>>(
    path: P,
    frames: &[Vec<(f64, f64)>],
    width: u16,
    height: u16,
    delay_cs: u16,
) -> Result<(), String> {
    let bytes = encode_series_animation_gif(frames, width, height, delay_cs)?;
    std::fs::write(path.as_ref(), &bytes)
        .map_err(|e| format!("cannot write {}: {e}", path.as_ref().display()))
}

/// Encodes a 2D particle trajectory animation (MD / N-body / SPH) into an animated GIF.
///
/// # Errors
/// Fails when `frames` is empty or dimensions are too small.
pub fn encode_particles_gif(
    frames: &[Vec<[f64; 2]>],
    width: u16,
    height: u16,
    delay_cs: u16,
    trail_len: usize,
) -> Result<Vec<u8>, String> {
    if frames.is_empty() {
        return Err("particle animation GIF requires at least one frame".to_owned());
    }
    if width < 32 || height < 32 {
        return Err("particle animation GIF dimensions must be at least 32x32".to_owned());
    }
    let mut xmin = f64::INFINITY;
    let mut xmax = f64::NEG_INFINITY;
    let mut ymin = f64::INFINITY;
    let mut ymax = f64::NEG_INFINITY;
    for f in frames {
        for &[x, y] in f {
            if x.is_finite() && y.is_finite() {
                xmin = xmin.min(x);
                xmax = xmax.max(x);
                ymin = ymin.min(y);
                ymax = ymax.max(y);
            }
        }
    }
    if !xmin.is_finite() {
        xmin = -1.0;
        xmax = 1.0;
        ymin = -1.0;
        ymax = 1.0;
    }
    let xspan = (xmax - xmin).max(1e-12);
    let yspan = (ymax - ymin).max(1e-12);
    let w = usize::from(width);
    let h = usize::from(height);
    let pad = 16_i32;
    let pw = (w as i32) - 2 * pad;
    let ph = (h as i32) - 2 * pad;
    let map_pt = |p: [f64; 2]| -> (i32, i32) {
        let px = pad + (((p[0] - xmin) / xspan) * f64::from(pw)).round() as i32;
        let py = pad + ph - (((p[1] - ymin) / yspan) * f64::from(ph)).round() as i32;
        (px, py)
    };

    let mut gif_frames = Vec::with_capacity(frames.len());
    for (fi, f) in frames.iter().enumerate() {
        let mut rgba = vec![255_u8; w * h * 4];
        // Trails
        if trail_len > 0 && fi > 0 {
            let start = fi.saturating_sub(trail_len);
            for t in start..fi {
                let prev = &frames[t];
                let next = &frames[t + 1];
                for (p0, p1) in prev.iter().zip(next.iter()) {
                    if p0[0].is_finite()
                        && p0[1].is_finite()
                        && p1[0].is_finite()
                        && p1[1].is_finite()
                    {
                        let (x0, y0) = map_pt(*p0);
                        let (x1, y1) = map_pt(*p1);
                        draw_line(&mut rgba, w, h, x0, y0, x1, y1, [160, 195, 230, 255]);
                    }
                }
            }
        }
        // Current particles
        for &p in f {
            if p[0].is_finite() && p[1].is_finite() {
                let (px, py) = map_pt(p);
                draw_disc(&mut rgba, w, h, px, py, 3, [214, 39, 40, 255]);
            }
        }
        gif_frames.push(GifFrame::new(width, height, delay_cs, rgba));
    }

    encode_gif(&GifAnimation::new(width, height, 0, gif_frames))
}

/// Writes a 2D particle trajectory animation GIF to `path`.
///
/// # Errors
/// As [`encode_particles_gif`], plus file I/O errors.
pub fn write_particles_gif<P: AsRef<Path>>(
    path: P,
    frames: &[Vec<[f64; 2]>],
    width: u16,
    height: u16,
    delay_cs: u16,
    trail_len: usize,
) -> Result<(), String> {
    let bytes = encode_particles_gif(frames, width, height, delay_cs, trail_len)?;
    std::fs::write(path.as_ref(), &bytes)
        .map_err(|e| format!("cannot write {}: {e}", path.as_ref().display()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encode_and_decode_gif_round_trip() {
        let frame1 = GifFrame::new(
            2,
            2,
            10,
            vec![
                255, 0, 0, 255, // red
                0, 255, 0, 255, // green
                0, 0, 255, 255, // blue
                255, 255, 255, 255, // white
            ],
        );
        let frame2 = GifFrame::new(
            2,
            2,
            20,
            vec![
                0, 0, 0, 255, // black
                255, 255, 0, 255, // yellow
                0, 255, 255, 255, // cyan
                255, 0, 255, 255, // magenta
            ],
        );
        let anim = GifAnimation::new(2, 2, 3, vec![frame1, frame2]);
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
        let anim = GifAnimation::new(1, 1, 0, vec![GifFrame::new(1, 1, 5, vec![12, 34, 56, 255])]);
        write_gif_file(&path, &anim).unwrap();
        let decoded = read_gif_file(&path).unwrap();
        let _r = std::fs::remove_file(&path);
        assert_eq!(decoded.width, 1);
        assert_eq!(decoded.height, 1);
        assert_eq!(decoded.frames.len(), 1);
    }

    #[test]
    fn invalid_gif_inputs_are_rejected() {
        let empty = GifAnimation::new(2, 2, 0, vec![]);
        assert!(encode_gif(&empty).is_err());

        let bad_buf = GifAnimation::new(2, 2, 0, vec![GifFrame::new(2, 2, 10, vec![0; 4])]);
        assert!(encode_gif(&bad_buf).is_err());

        let mut oob_frame = GifFrame::new(2, 2, 10, vec![255; 16]);
        oob_frame.left = 2;
        assert!(encode_gif(&GifAnimation::new(3, 3, 0, vec![oob_frame])).is_err());
    }
}
