//! Harvests the `///` documentation of public functions at build time.
//!
//! The Kernel Manifest describes every operator, and an operator is usually a
//! thin registration around a public Rust function whose rustdoc already says
//! what it computes, what its parameters mean and what it returns. This scanner
//! reads `src/**/*.rs` once per build, with plain line scanning (no `syn`, so no
//! build-dependency and no parse of the 370k lines), and writes a table that
//! `kernel::rustdoc` embeds with `include_str!`:
//!
//! ```text
//! <path> US <name> US <doc> RS <path> US ...
//! ```
//!
//! `path` is the module path of the function (`crate::kernels::special::erf`,
//! or `crate::x::Type::method` for an inherent method), the key an operator's
//! `rust_path` uses. Only `pub fn` items with at least one `///` line are
//! recorded. Hidden doctest lines (`# ...` inside a code fence) are dropped
//! because they are not part of what a reader sees.
//!
//! The scan is deterministic (directories are sorted) and the output file is
//! rewritten only when its content changes.

use std::fs;
use std::io;
use std::path::Path;

const UNIT: char = '\u{1f}';
const RECORD: char = '\u{1e}';
/// A single function's doc is cut here; a longer one is an essay, not a manifest entry.
const MAX_DOC: usize = 12 * 1024;

struct Entry {
    path: String,
    name: String,
    doc: String,
}

/// Scans `src` and writes the table to `target` when it changed.
pub fn generate(
    src: &Path,
    target: &Path,
) -> io::Result<()> {
    let mut entries = Vec::new();
    // Without a source tree (the build script run on its own, from an
    // archive of the data files) there is nothing to harvest: the manifest
    // keeps the declared docs only.
    if src.is_dir() {
        walk(src, src, &mut entries)?;
    }
    let mut out = String::new();
    for entry in &entries {
        out.push_str(&entry.path);
        out.push(UNIT);
        out.push_str(&entry.name);
        out.push(UNIT);
        out.push_str(&entry.doc);
        out.push(RECORD);
    }
    if fs::read_to_string(target).ok().as_deref() != Some(out.as_str()) {
        fs::write(target, out)?;
    }
    Ok(())
}

fn walk(
    root: &Path,
    dir: &Path,
    entries: &mut Vec<Entry>,
) -> io::Result<()> {
    let mut children: Vec<_> = fs::read_dir(dir)?
        .filter_map(Result::ok)
        .map(|entry| entry.path())
        .collect();
    children.sort();
    for path in children {
        if path.is_dir() {
            // Binaries are not library API.
            if path.file_name().is_some_and(|name| name == "bin") {
                continue;
            }
            walk(root, &path, entries)?;
        } else if path.extension().is_some_and(|ext| ext == "rs") {
            let stem = path
                .file_stem()
                .map(|stem| stem.to_string_lossy().into_owned())
                .unwrap_or_default();
            if stem == "main" || stem == "tests" || stem.ends_with("_tests") {
                continue;
            }
            let Some(module) = module_path(root, &path) else {
                continue;
            };
            let Ok(source) = fs::read_to_string(&path) else {
                continue;
            };
            scan(&module, &source, entries);
        }
    }
    Ok(())
}

/// `src/a/b.rs` and `src/a/b/mod.rs` are `crate::a::b`; `src/lib.rs` is `crate`.
fn module_path(
    root: &Path,
    file: &Path,
) -> Option<String> {
    let relative = file.strip_prefix(root).ok()?;
    let mut parts: Vec<String> = relative
        .iter()
        .map(|part| part.to_string_lossy().into_owned())
        .collect();
    let last = parts.pop()?;
    let stem = last.strip_suffix(".rs")?;
    if stem != "mod" && stem != "lib" {
        parts.push(stem.to_owned());
    }
    let mut path = String::from("crate");
    for part in parts {
        path.push_str("::");
        path.push_str(&part);
    }
    Some(path)
}

fn scan(
    module: &str,
    source: &str,
    entries: &mut Vec<Entry>,
) {
    let mut pending: Vec<String> = Vec::new();
    let mut hidden = false;
    let mut in_impl: Option<String> = None;
    for raw in source.lines() {
        let line = raw.trim_end();
        let indent = line.len() - line.trim_start().len();
        let text = line.trim_start();
        if let Some(rest) = text.strip_prefix("///") {
            if !rest.starts_with('/') {
                pending.push(rest.strip_prefix(' ').unwrap_or(rest).to_owned());
            }
            continue;
        }
        if text.starts_with("#[") {
            if text.contains("doc(hidden)") {
                hidden = true;
            }
            continue;
        }
        if text.starts_with("//") {
            continue;
        }
        if indent == 0 {
            if text.starts_with('}') {
                in_impl = None;
            } else if text.starts_with("impl") {
                in_impl = inherent_impl_type(text);
            }
        }
        let wanted_indent = if in_impl.is_some() { 4 } else { 0 };
        if indent == wanted_indent
            && let Some(name) = public_fn_name(text)
            && !pending.is_empty()
            && !hidden
        {
            let path = match &in_impl {
                | Some(ty) => format!("{module}::{ty}::{name}"),
                | None => format!("{module}::{name}"),
            };
            let doc = clean(&pending);
            if !doc.is_empty() {
                entries.push(Entry {
                    path,
                    name: name.to_owned(),
                    doc,
                });
            }
        }
        pending.clear();
        hidden = false;
    }
}

/// The type of `impl Type {` / `impl<T> Type<T> {`; `None` for a trait impl,
/// whose methods are documented by the trait.
fn inherent_impl_type(text: &str) -> Option<String> {
    if !text.trim_end().ends_with('{') && !text.contains(" where") {
        return None;
    }
    let mut rest = text.strip_prefix("impl")?;
    if rest.starts_with('<') {
        let mut depth = 0_usize;
        let mut end = None;
        for (index, ch) in rest.char_indices() {
            match ch {
                | '<' => depth += 1,
                | '>' => {
                    depth -= 1;
                    if depth == 0 {
                        end = Some(index + 1);
                        break;
                    }
                },
                | _ => {},
            }
        }
        rest = &rest[end?..];
    }
    let rest = rest.trim_start();
    if rest.contains(" for ") {
        return None;
    }
    let name: String = rest
        .chars()
        .take_while(|ch| ch.is_alphanumeric() || *ch == '_')
        .collect();
    (!name.is_empty()).then_some(name)
}

/// The name of a `pub fn` (also `pub const fn`, `pub unsafe fn`, `pub async fn`).
fn public_fn_name(text: &str) -> Option<&str> {
    let mut rest = text.strip_prefix("pub ")?;
    loop {
        if let Some(next) = rest.strip_prefix("fn ") {
            rest = next;
            break;
        }
        let word_end = rest.find(' ')?;
        if !matches!(&rest[..word_end], "const" | "unsafe" | "async") {
            return None;
        }
        rest = &rest[word_end + 1..];
    }
    let end = rest
        .find(|ch: char| !(ch.is_alphanumeric() || ch == '_'))
        .unwrap_or(rest.len());
    (end > 0).then(|| &rest[..end])
}

/// Drops hidden doctest lines and trims blank edges.
fn clean(lines: &[String]) -> String {
    let mut out: Vec<&str> = Vec::new();
    let mut fenced = false;
    for line in lines {
        let trimmed = line.trim_start();
        if trimmed.starts_with("```") {
            fenced = !fenced;
        } else if fenced && (trimmed == "#" || trimmed.starts_with("# ")) {
            continue;
        }
        out.push(line);
    }
    let text = out.join("\n");
    let text = text.trim();
    if text.len() <= MAX_DOC {
        return text.to_owned();
    }
    let mut cut = MAX_DOC;
    while !text.is_char_boundary(cut) {
        cut -= 1;
    }
    let head = &text[..cut];
    let head = head.rfind("\n\n").map_or(head, |at| &head[..at]);
    head.trim_end().to_owned()
}
