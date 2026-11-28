//! Complete RIFF/WAVE audio reader and writer.
//!
//! Supports:
//! - PCM integer (`8`-bit unsigned, `16`-bit signed, `24`-bit signed 3-byte LE, `32`-bit signed)
//! - IEEE float (`32`-bit `f32`, `64`-bit `f64`) with mandatory `fact` chunk
//! - `WAVE_FORMAT_EXTENSIBLE` (`0xFFFE`) with `channel_mask` and `SubFormat` GUID
//! - `LIST`/`INFO` metadata chunks (`INAM`, `IART`, `ICMT`, `ISFT`, etc.)
//! - Strict RIFF size, `block_align`, `byte_rate`, and `data` alignment validation

use std::collections::BTreeMap;
use std::path::Path;

const KSDATAFORMAT_SUBTYPE_PCM: [u8; 16] = [
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10, 0x00, 0x80, 0x00, 0x00, 0xAA, 0x00, 0x38, 0x9B, 0x71,
];

const KSDATAFORMAT_SUBTYPE_IEEE_FLOAT: [u8; 16] = [
    0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10, 0x00, 0x80, 0x00, 0x00, 0xAA, 0x00, 0x38, 0x9B, 0x71,
];

/// Audio stream format specification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WavSpec {
    /// Number of channels (`1` = mono, `2` = stereo, etc.).
    pub channels: u16,
    /// Sample rate in Hz (e.g. `44100`, `48000`).
    pub sample_rate: u32,
    /// Bits per sample (`8`, `16`, `24`, `32` for PCM; `32`, `64` for IEEE float).
    pub bits_per_sample: u16,
    /// `true` for IEEE float, `false` for integer PCM.
    pub is_float: bool,
}

impl Default for WavSpec {
    fn default() -> Self {
        Self {
            channels: 1,
            sample_rate: 44_100,
            bits_per_sample: 16,
            is_float: false,
        }
    }
}

/// Decoded or ready-to-encode WAV audio with normalized `[-1.0, 1.0]` samples.
#[derive(Debug, Clone, PartialEq)]
pub struct WavData {
    /// Stream format.
    pub spec: WavSpec,
    /// Interleaved samples normalized to `[-1.0, 1.0]`.
    pub samples: Vec<f64>,
    /// Optional speaker channel mask (triggers `WAVE_FORMAT_EXTENSIBLE` when set).
    pub channel_mask: Option<u32>,
    /// Optional `LIST`/`INFO` metadata tags (e.g. `"INAM"`, `"IART"`).
    pub info: BTreeMap<String, String>,
}

impl WavData {
    /// Creates a [`WavData`] without extensible channel mask or metadata.
    #[must_use]
    pub const fn new(
        spec: WavSpec,
        samples: Vec<f64>,
    ) -> Self {
        Self {
            spec,
            samples,
            channel_mask: None,
            info: BTreeMap::new(),
        }
    }
}

/// Encodes [`WavData`] into RIFF/WAVE bytes.
///
/// # Errors
/// Fails when `channels == 0`, `sample_rate == 0`, sample count is not a multiple
/// of `channels`, or `(bits_per_sample, is_float)` is unsupported.
pub fn encode_wav(wav: &WavData) -> Result<Vec<u8>, String> {
    if wav.spec.channels == 0 || wav.spec.sample_rate == 0 {
        return Err("WAV channels and sample_rate must be positive".to_owned());
    }
    if !wav
        .samples
        .len()
        .is_multiple_of(usize::from(wav.spec.channels))
    {
        return Err("WAV sample count must be a multiple of channels".to_owned());
    }
    let bps: usize = match (wav.spec.bits_per_sample, wav.spec.is_float) {
        | (8, false) => 1,
        | (16, false) => 2,
        | (24, false) => 3,
        | (32, _) => 4,
        | (64, true) => 8,
        | _ => {
            return Err(format!(
                "unsupported WAV format: {} bits, float={}",
                wav.spec.bits_per_sample, wav.spec.is_float
            ));
        },
    };
    let use_extensible = wav.channel_mask.is_some();
    let format_tag: u16 = if use_extensible {
        0xFFFE
    } else if wav.spec.is_float {
        3
    } else {
        1
    };
    let block_align = wav.spec.channels * (bps as u16);
    let byte_rate = wav.spec.sample_rate * u32::from(block_align);
    let data_bytes = wav.samples.len() * bps;
    let frame_count = (wav.samples.len() / usize::from(wav.spec.channels)) as u32;

    let mut body = Vec::with_capacity(64 + data_bytes);
    body.extend_from_slice(b"WAVE");

    // `fmt ` chunk
    body.extend_from_slice(b"fmt ");
    let fmt_size: u32 = if use_extensible { 40 } else { 16 };
    body.extend_from_slice(&fmt_size.to_le_bytes());
    body.extend_from_slice(&format_tag.to_le_bytes());
    body.extend_from_slice(&wav.spec.channels.to_le_bytes());
    body.extend_from_slice(&wav.spec.sample_rate.to_le_bytes());
    body.extend_from_slice(&byte_rate.to_le_bytes());
    body.extend_from_slice(&block_align.to_le_bytes());
    body.extend_from_slice(&wav.spec.bits_per_sample.to_le_bytes());
    if use_extensible {
        body.extend_from_slice(&22_u16.to_le_bytes()); // cbSize
        body.extend_from_slice(&wav.spec.bits_per_sample.to_le_bytes()); // validBitsPerSample
        body.extend_from_slice(&wav.channel_mask.unwrap_or(0).to_le_bytes());
        let guid = if wav.spec.is_float {
            KSDATAFORMAT_SUBTYPE_IEEE_FLOAT
        } else {
            KSDATAFORMAT_SUBTYPE_PCM
        };
        body.extend_from_slice(&guid);
    }

    // `fact` chunk (required by RIFF/WAVE spec for non-PCM formats)
    if wav.spec.is_float || use_extensible {
        body.extend_from_slice(b"fact");
        body.extend_from_slice(&4_u32.to_le_bytes());
        body.extend_from_slice(&frame_count.to_le_bytes());
    }

    // Optional `LIST`/`INFO` metadata chunk
    if !wav.info.is_empty() {
        let mut info_payload = Vec::new();
        info_payload.extend_from_slice(b"INFO");
        for (tag, val) in &wav.info {
            if tag.len() != 4 {
                continue;
            }
            let mut text = val.as_bytes().to_vec();
            text.push(0); // NUL-terminated
            info_payload.extend_from_slice(tag.as_bytes());
            info_payload.extend_from_slice(&(text.len() as u32).to_le_bytes());
            info_payload.extend_from_slice(&text);
            if text.len() % 2 != 0 {
                info_payload.push(0);
            }
        }
        body.extend_from_slice(b"LIST");
        body.extend_from_slice(&(info_payload.len() as u32).to_le_bytes());
        body.extend_from_slice(&info_payload);
        if info_payload.len() % 2 != 0 {
            body.push(0);
        }
    }

    // `data` chunk
    body.extend_from_slice(b"data");
    body.extend_from_slice(&(data_bytes as u32).to_le_bytes());
    match (wav.spec.bits_per_sample, wav.spec.is_float) {
        | (8, false) => {
            for &s in &wav.samples {
                let c = if s.is_finite() {
                    s.clamp(-1.0, 1.0)
                } else {
                    0.0
                };
                let q = ((c * 127.0) + 128.0).round().clamp(0.0, 255.0) as u8;
                body.push(q);
            }
        },
        | (16, false) => {
            for &s in &wav.samples {
                let c = if s.is_finite() {
                    s.clamp(-1.0, 1.0)
                } else {
                    0.0
                };
                let q = (c * 32767.0).round().clamp(-32768.0, 32767.0) as i16;
                body.extend_from_slice(&q.to_le_bytes());
            }
        },
        | (24, false) => {
            for &s in &wav.samples {
                let c = if s.is_finite() {
                    s.clamp(-1.0, 1.0)
                } else {
                    0.0
                };
                let q = (c * 8_388_607.0).round().clamp(-8_388_608.0, 8_388_607.0) as i32;
                let le = q.to_le_bytes();
                body.extend_from_slice(&le[..3]);
            }
        },
        | (32, false) => {
            for &s in &wav.samples {
                let c = if s.is_finite() {
                    s.clamp(-1.0, 1.0)
                } else {
                    0.0
                };
                let q = (c * 2_147_483_647.0)
                    .round()
                    .clamp(-2_147_483_648.0, 2_147_483_647.0) as i32;
                body.extend_from_slice(&q.to_le_bytes());
            }
        },
        | (32, true) => {
            for &s in &wav.samples {
                body.extend_from_slice(&(s as f32).to_le_bytes());
            }
        },
        | (64, true) => {
            for &s in &wav.samples {
                body.extend_from_slice(&s.to_le_bytes());
            }
        },
        | _ => unreachable!(),
    }
    if !data_bytes.is_multiple_of(2) {
        body.push(0);
    }

    let riff_size =
        u32::try_from(body.len()).map_err(|_| "WAV data exceeds 4 GiB RIFF limit".to_owned())?;
    let mut out = Vec::with_capacity(8 + body.len());
    out.extend_from_slice(b"RIFF");
    out.extend_from_slice(&riff_size.to_le_bytes());
    out.extend_from_slice(&body);
    Ok(out)
}

/// Decodes RIFF/WAVE bytes into [`WavData`] with strict header & chunk validation.
///
/// # Errors
/// Fails when the header is invalid, sizes are inconsistent, or the audio format is unsupported.
pub fn decode_wav(bytes: &[u8]) -> Result<WavData, String> {
    if bytes.len() < 44 || bytes.get(..4) != Some(b"RIFF") || bytes.get(8..12) != Some(b"WAVE") {
        return Err("invalid RIFF/WAVE header".to_owned());
    }
    let riff_size = u32::from_le_bytes([bytes[4], bytes[5], bytes[6], bytes[7]]) as usize;
    if riff_size + 8 > bytes.len() {
        return Err(format!(
            "RIFF chunk declares {} bytes, but buffer is only {} bytes",
            riff_size + 8,
            bytes.len()
        ));
    }

    let mut pos = 12_usize;
    let end_limit = (riff_size + 8).min(bytes.len());
    let mut spec: Option<WavSpec> = None;
    let mut channel_mask: Option<u32> = None;
    let mut data_slice: Option<&[u8]> = None;
    let mut info = BTreeMap::new();

    while pos + 8 <= end_limit {
        let id = &bytes[pos..pos + 4];
        let size = u32::from_le_bytes([
            bytes[pos + 4],
            bytes[pos + 5],
            bytes[pos + 6],
            bytes[pos + 7],
        ]) as usize;
        pos += 8;
        let payload = bytes
            .get(pos..pos.saturating_add(size))
            .ok_or_else(|| "WAV chunk truncated".to_owned())?;
        if id == b"fmt " {
            if payload.len() < 16 {
                return Err("WAV fmt chunk truncated".to_owned());
            }
            let fmt = u16::from_le_bytes([payload[0], payload[1]]);
            let channels = u16::from_le_bytes([payload[2], payload[3]]);
            let sample_rate = u32::from_le_bytes([payload[4], payload[5], payload[6], payload[7]]);
            let byte_rate = u32::from_le_bytes([payload[8], payload[9], payload[10], payload[11]]);
            let block_align = u16::from_le_bytes([payload[12], payload[13]]);
            let bits_per_sample = u16::from_le_bytes([payload[14], payload[15]]);
            if channels == 0 || sample_rate == 0 {
                return Err("WAV channels and sample_rate must be positive".to_owned());
            }
            let expected_align = channels * (bits_per_sample / 8);
            if expected_align == 0 || block_align != expected_align {
                return Err(format!(
                    "WAV block_align {block_align} != expected {expected_align}"
                ));
            }
            let expected_rate = sample_rate * u32::from(block_align);
            if byte_rate != expected_rate {
                return Err(format!(
                    "WAV byte_rate {byte_rate} != expected {expected_rate}"
                ));
            }
            let is_float = match fmt {
                | 1 => false,
                | 3 => true,
                | 0xFFFE => {
                    if payload.len() < 40 {
                        return Err("WAVE_FORMAT_EXTENSIBLE fmt chunk < 40 bytes".to_owned());
                    }
                    channel_mask = Some(u32::from_le_bytes([
                        payload[20],
                        payload[21],
                        payload[22],
                        payload[23],
                    ]));
                    let guid = &payload[24..40];
                    if guid == KSDATAFORMAT_SUBTYPE_PCM {
                        false
                    } else if guid == KSDATAFORMAT_SUBTYPE_IEEE_FLOAT {
                        true
                    } else {
                        return Err("unsupported WAVE_FORMAT_EXTENSIBLE SubFormat GUID".to_owned());
                    }
                },
                | _ => return Err(format!("unsupported WAV format tag {fmt}")),
            };
            spec = Some(WavSpec {
                channels,
                sample_rate,
                bits_per_sample,
                is_float,
            });
        } else if id == b"data" {
            data_slice = Some(payload);
        } else if id == b"LIST" && payload.get(..4) == Some(b"INFO") {
            let mut ip = 4_usize;
            while ip + 8 <= payload.len() {
                let tag = String::from_utf8_lossy(&payload[ip..ip + 4]).into_owned();
                let tlen = u32::from_le_bytes([
                    payload[ip + 4],
                    payload[ip + 5],
                    payload[ip + 6],
                    payload[ip + 7],
                ]) as usize;
                ip += 8;
                if ip + tlen > payload.len() {
                    break;
                }
                let raw_str = &payload[ip..ip + tlen];
                let trimmed = raw_str.split(|&b| b == 0).next().unwrap_or(raw_str);
                info.insert(tag, String::from_utf8_lossy(trimmed).into_owned());
                ip += tlen + (tlen % 2);
            }
        }
        pos += size + (size % 2);
    }

    let spec = spec.ok_or_else(|| "WAV missing fmt chunk".to_owned())?;
    let data = data_slice.ok_or_else(|| "WAV missing data chunk".to_owned())?;
    let frame_bytes = usize::from(spec.channels) * usize::from(spec.bits_per_sample / 8);
    if data.len() % frame_bytes != 0 {
        return Err("WAV data chunk length is not aligned to frame size".to_owned());
    }

    let samples = match (spec.bits_per_sample, spec.is_float) {
        | (8, false) => {
            data.iter()
                .map(|&b| (f64::from(b) - 128.0) / 128.0)
                .collect()
        },
        | (16, false) => {
            data.as_chunks::<2>()
                .0
                .iter()
                .map(|c| f64::from(i16::from_le_bytes(*c)) / 32768.0)
                .collect()
        },
        | (24, false) => {
            data.as_chunks::<3>()
                .0
                .iter()
                .map(|c| {
                    let sign = if (c[2] & 0x80) != 0 {
                        0xFF
                    } else {
                        0x00
                    };
                    let v = i32::from_le_bytes([c[0], c[1], c[2], sign]);
                    f64::from(v) / 8_388_608.0
                })
                .collect()
        },
        | (32, false) => {
            data.as_chunks::<4>()
                .0
                .iter()
                .map(|c| f64::from(i32::from_le_bytes(*c)) / 2_147_483_648.0)
                .collect()
        },
        | (32, true) => {
            data.as_chunks::<4>()
                .0
                .iter()
                .map(|c| f64::from(f32::from_le_bytes(*c)))
                .collect()
        },
        | (64, true) => {
            data.as_chunks::<8>()
                .0
                .iter()
                .map(|c| f64::from_le_bytes(*c))
                .collect()
        },
        | _ => {
            return Err(format!(
                "unsupported WAV sample width: {} bits (float={})",
                spec.bits_per_sample, spec.is_float
            ));
        },
    };
    Ok(WavData {
        spec,
        samples,
        channel_mask,
        info,
    })
}

/// Writes [`WavData`] to `path`.
///
/// # Errors
/// As [`encode_wav`], plus file I/O errors.
pub fn write_wav_file<P: AsRef<Path>>(
    path: P,
    wav: &WavData,
) -> Result<(), String> {
    let bytes = encode_wav(wav)?;
    std::fs::write(path.as_ref(), &bytes)
        .map_err(|e| format!("cannot write {}: {e}", path.as_ref().display()))
}

/// Reads and decodes a WAV file from `path`.
///
/// # Errors
/// Fails when the file cannot be read or decoded.
pub fn read_wav_file<P: AsRef<Path>>(path: P) -> Result<WavData, String> {
    let bytes = std::fs::read(path.as_ref())
        .map_err(|e| format!("cannot read {}: {e}", path.as_ref().display()))?;
    decode_wav(&bytes)
}

/// Auto-normalizes a mono `f64` signal and writes it as a 16-bit PCM WAV file.
///
/// # Errors
/// Fails when `signal` is empty or the file cannot be written.
pub fn write_signal_wav<P: AsRef<Path>>(
    path: P,
    signal: &[f64],
    sample_rate: u32,
) -> Result<(), String> {
    if signal.is_empty() {
        return Err("signal cannot be empty".to_owned());
    }
    let peak = signal
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .map(f64::abs)
        .fold(0.0_f64, f64::max);
    let scale = if peak > 0.0 {
        0.95 / peak
    } else {
        1.0
    };
    let samples = signal
        .iter()
        .map(|&v| {
            if v.is_finite() {
                v * scale
            } else {
                0.0
            }
        })
        .collect();
    write_wav_file(
        path,
        &WavData::new(
            WavSpec {
                channels: 1,
                sample_rate,
                bits_per_sample: 16,
                is_float: false,
            },
            samples,
        ),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wav_all_bit_depths_extensible_and_info_round_trip() {
        for (bits, is_float) in [
            (8, false),
            (16, false),
            (24, false),
            (32, false),
            (32, true),
            (64, true),
        ] {
            let mut wav = WavData::new(
                WavSpec {
                    channels: 2,
                    sample_rate: 48_000,
                    bits_per_sample: bits,
                    is_float,
                },
                vec![0.0, 0.5, -0.5, 0.25],
            );
            wav.info.insert("INAM".to_owned(), "RSSN Tone".to_owned());
            if bits == 24 {
                wav.channel_mask = Some(0x03);
            }
            let bytes = encode_wav(&wav).unwrap();
            let back = decode_wav(&bytes).unwrap();
            assert_eq!(back.spec, wav.spec);
            assert_eq!(back.info.get("INAM").map(String::as_str), Some("RSSN Tone"));
            for (a, b) in back.samples.iter().zip(&wav.samples) {
                assert!(
                    (a - b).abs() < 0.02,
                    "bits={bits} float={is_float}: {a} vs {b}"
                );
            }
        }
    }
}
