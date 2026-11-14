---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

### 📦 Bug Description

Please clearly and concisely describe the bug you encountered. Specify which module it occurs in (e.g., `rules::calculus`, `kernels::ode`, `sim::scenario`, `graph::dim`, `ffi`).

### 🔢 Steps to Reproduce

A minimal, reproducible example is key to fixing bugs quickly.

1.  `rssn` Version Used: [e.g., `3.0.0-rc.1`]
2.  Your `Cargo.toml` dependency configuration (especially the features enabled):
    ```toml
    # [Paste relevant dependency configuration here]
    ```
3.  Minimal, runnable Rust code that reproduces the bug:
    ```rust
    // [Paste your runnable Rust code snippet here]
    
    fn main() {
        // ...
    }
    ```
4.  The command used to execute the code (e.g., `cargo run --release --features jit`)

### 📝 Expected Behavior

What did you expect the code to do or the result to be?

### Observed Behavior

What actually happened? Please provide the full error message, Panic Trace, or the unexpected calculated result.

### ⚙️ Your Environment

* **Operating System:** [e.g., Windows 10 / macOS 14 / Ubuntu 22.04]
* **Rust Version:** [Obtained via `rustc --version`]
* **CPU Architecture:** [e.g., x86_64, aarch64]

### 💡 Suggested Fix (Optional)

If you have any theories about the cause of the bug or a potential fix, please describe it here.
