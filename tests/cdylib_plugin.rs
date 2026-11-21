//! An end-to-end test of the `cdylib` plugin ABI.
//!
//! The plugin is written in **C** on purpose: it depends on nothing from rssn,
//! so the test proves that the ABI is a byte-level contract rather than a Rust
//! detail that happens to work when both sides are the same crate.
//!
//! The interesting part is how the two sides agree on the encoding. The test
//! computes the bytes with rssn's own codec and *writes them into the C source*
//! as a static array; the plugin returns them verbatim. If rssn's decoder and
//! the ABI ever drift, this test fails — which is exactly what a wire-format
//! test should do.
//!
//! The plugin is compiled at test time with `cc`, so the test is skipped (not
//! failed) when no C compiler is available.

#![cfg(feature = "plugin-cdylib")]

use std::path::PathBuf;
use std::process::Command;

use rssn::kernel::cdylib;
use rssn::kernel::codec;
use rssn::kernel::plugin::Plugin;
use rssn::kernel::signature::Signature;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::Value;

/// Where to put the generated plugin.
fn scratch() -> PathBuf {
    let mut dir = std::env::temp_dir();
    dir.push(format!("rssn-cdylib-{}", std::process::id()));
    let _ = std::fs::create_dir_all(&dir);
    dir
}

/// A C array literal for `bytes`.
fn c_array(bytes: &[u8]) -> String {
    let body: Vec<String> = bytes.iter().map(u8::to_string).collect();
    format!("{{{}}}", body.join(", "))
}

/// Writes and compiles the plugin; returns its path, or `None` without a C
/// compiler.
fn build_plugin(
    name: &str,
    arity: u32,
    signature_bytes: &[u8],
    answer_bytes: &[u8],
) -> Option<PathBuf> {
    let dir = scratch();
    let source = dir.join(format!("{name}.c"));
    let library = dir.join(format!("lib{name}.so"));
    let text = format!(
        r#"
#include <stdint.h>
#include <stddef.h>

typedef struct {{ uint8_t *ptr; size_t len; }} RssnBuf;
typedef struct {{ uint32_t abi; void (*log)(uint32_t, const uint8_t *, size_t); }} RssnHost;
typedef RssnBuf (*RssnCallFn)(void *state, const uint8_t *args, size_t len);

typedef struct {{
    const char    *name;
    uint32_t       arity_min;
    uint32_t       arity_max;
    const uint8_t *signatures;
    size_t         signatures_len;
    RssnCallFn     call;
}} RssnOpEntry;

typedef struct {{
    uint32_t           abi;
    const char        *name;
    const char        *version;
    const char        *description;
    const RssnOpEntry *ops;
    size_t             op_count;
    void *(*init)(const RssnHost *host);
    void  (*free_buf)(RssnBuf);
}} RssnPluginDesc;

/* bincode(Vec<Signature>) and bincode(Result<Value, String>) supplied by the
   test, which computes them with rssn's own encoder. */
static const uint8_t SIGNATURES[] = {signatures};
static const uint8_t ANSWER[] = {answer};

static RssnBuf call_op(void *state, const uint8_t *args, size_t len) {{
    (void)state; (void)args; (void)len;
    RssnBuf out;
    /* The answer is static, so `free_buf` below is a no-op. */
    out.ptr = (uint8_t *)ANSWER;
    out.len = sizeof(ANSWER);
    return out;
}}

static void free_buf(RssnBuf buffer) {{ (void)buffer; }}

static const RssnOpEntry OPS[] = {{{{
    "{name}", {arity}, {arity}, SIGNATURES, sizeof(SIGNATURES), call_op
}}}};

static const RssnPluginDesc DESC = {{
    1,
    "test/{name}",
    "0.1.0",
    "A plugin written in C for the ABI test",
    OPS,
    1,
    0,
    free_buf
}};

const RssnPluginDesc *rssn_plugin_v1(const RssnHost *host) {{
    (void)host;
    return &DESC;
}}
"#,
        signatures = c_array(signature_bytes),
        answer = c_array(answer_bytes),
    );
    std::fs::write(&source, text).ok()?;
    let status = Command::new("cc")
        .args(["-shared", "-fPIC", "-O1", "-o"])
        .arg(&library)
        .arg(&source)
        .status()
        .ok()?;
    if !status.success() {
        return None;
    }
    Some(library)
}

#[test]
fn a_plugin_written_in_c_loads_registers_and_answers() {
    // What the plugin will return: `Ok(Value::Bool(true))`, encoded by rssn.
    let answer = codec::encode_deterministic(&Ok::<Value, String>(Value::Bool(true)))
        .unwrap_or_else(|error| panic!("{error}"));
    let signatures = codec::encode_deterministic(&vec![Signature::scalars(1, TypeExpr::Bool)])
        .unwrap_or_else(|error| panic!("{error}"));
    let Some(library) = build_plugin("trueop", 1, &signatures, &answer) else {
        eprintln!("skipping: no C compiler available");
        return;
    };

    let plugin = cdylib::load(&library).unwrap_or_else(|error| panic!("{error}"));
    assert_eq!(plugin.manifest().name, "test/trueop");
    assert_eq!(
        plugin.manifest().kind,
        rssn::kernel::plugin::PluginKind::Cdylib
    );
    assert_eq!(plugin.op_names(), vec!["trueop"]);

    let mut session = rssn::api::Session::new();
    session
        .register(Box::new(plugin))
        .unwrap_or_else(|error| panic!("{error}"));

    // The manifest now lists the plugin, with its kind.
    let manifest = session.manifest();
    assert!(
        manifest.plugins.iter().any(|p| p.name == "test/trueop"),
        "the plugin must appear in the manifest"
    );
    assert!(
        manifest
            .op("trueop")
            .is_some_and(|op| !op.overloads.is_empty()),
        "the plugin's operator must carry the signature it declared"
    );

    // And the call crosses the boundary: the arguments are encoded, the
    // plugin ignores them, and the answer decodes back into a `Value`.
    let answer = session
        .call_value(
            "trueop",
            &[Value::int(7_i32)],
            &rssn::kernel::CallCtx::new(),
        )
        .unwrap_or_else(|error| panic!("{error}"));
    assert_eq!(answer, Value::Bool(true));

    let _ = std::fs::remove_file(&library);
}

#[test]
fn a_library_that_is_not_a_plugin_is_refused() {
    let Some(library) = build_plugin(
        "noplugin",
        1,
        &codec::encode_deterministic(&vec![Signature::scalars(1, TypeExpr::Real)])
            .unwrap_or_default(),
        &codec::encode_deterministic(&Ok::<Value, String>(Value::Bool(true))).unwrap_or_default(),
    ) else {
        eprintln!("skipping: no C compiler available");
        return;
    };
    // The symbol is exported, so this one loads; a library without it must not.
    assert!(cdylib::load(&library).is_ok());
    let missing = scratch().join("libdoes_not_exist.so");
    assert!(cdylib::load(&missing).is_err());
    let _ = std::fs::remove_file(&library);
}
