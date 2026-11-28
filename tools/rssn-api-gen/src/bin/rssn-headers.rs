//! Regenerate or verify the C and C++ headers without compiling libRSSN.

#[path = "../../../../buildsrc/headers.rs"]
mod headers;

use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let action = args
        .next()
        .ok_or("usage: rssn-headers (--write|--check) [RSSN_ROOT]")?;
    let root = args.next().map_or_else(
        || PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../.."),
        PathBuf::from,
    );
    if !matches!(action.as_str(), "--write" | "--check") || args.next().is_some() {
        return Err("usage: rssn-headers (--write|--check) [RSSN_ROOT]".into());
    }
    for (name, bindings) in headers::generate(&root)? {
        let path = root.join(name);
        let generated = headers::header_bytes(name, &bindings)?;
        if action == "--write" {
            std::fs::write(path, generated)?;
        } else {
            let expected = String::from_utf8(generated)?.replace("\r\n", "\n");
            if std::fs::read_to_string(&path)?.replace("\r\n", "\n") != expected {
                return Err(format!(
                    "{} is stale; run bazel run //:refresh_headers",
                    path.display()
                )
                .into());
            }
        }
    }
    Ok(())
}
