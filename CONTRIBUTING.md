# Contributing to rssn

First of all, thank you for considering contributing to **rssn**!  
This project aims to become a next-generation scientific computing ecosystem in Rust, and your help is highly appreciated.

---

## 🔧 Development Workflow

1. **Fork and clone** the repository.  
   ```bash
   git clone https://github.com/Apich-Organization/rssn.git
   cd rssn
````

2. **Set up the environment**:

   * Rust (latest stable)
   * Cargo

3. **Code style & formatting**:

   * All code must compile with **zero warnings** on the latest stable Rust.
   * Additional **lint rules** are configured in `lib.rs` and must be respected.
   * To ensure maximum readability and long-term maintainability, we require all contributions to follow a strict language standard.
   * Please avoid using abbreviations for variable names, function names, and comments.
   * Always use full words and complete phrases to clearly describe your intent (e.g., use message instead of msg, initialization instead of init).
   * The only exception is for abbreviations that are widely recognized and unambiguous industry standards (e.g., HTTP, JSON, API).
   * This helps new contributors quickly understand the codebase and significantly reduces cognitive load during code reviews.
   * Before pushing, always run:

     ```bash
     cargo fmt --all
     cargo clippy --all-targets -- -D warnings
     cargo test --all
     ```

4. **AI reviewer**:
   Every pull request is automatically reviewed by an AI-assisted reviewer.
   Please write clear commit messages and PR descriptions to help the review process.

---

## 🧪 Testing

* Add unit tests for new features in the corresponding module.
* Integration tests should be placed in the `tests/` directory.
* Benchmarks can be added under `benches/`.

We follow the principle: **new features require tests, bug fixes require regression tests.**

---

## 📖 Documentation

* All public functions, structs, and traits must include `///` doc comments.
* Use `cargo doc --open` to locally verify documentation.
* Examples should be concise and runnable.

---

## ✅ Contribution Areas

* **Bug fixes**: Help us improve stability.
* **Performance improvements**: Optimize algorithms and solvers.
* **New functionality**: Expand symbolic, numerical, physics, or output modules.
* **Testing**: Improve coverage and add benchmarks.
* **Documentation**: Enhance clarity and usability.

---

## 📬 Submitting Changes

1. Create a feature branch:

   ```bash
   git checkout -b feature/my-new-feature
   ```

2. Commit changes with clear messages:

   ```bash
   git commit -m "feat(symbolic): add new integration method"
   ```

3. Push and open a Pull Request:

   ```bash
   git push origin feature/my-new-feature
   ```

4. Ensure CI checks (format, lint, tests) all pass.

---

## 🙏 Acknowledgements

Contributors are credited in release notes and on the GitHub page.
We value every contribution, from fixing typos to implementing new solvers.

Thank you for making **rssn** better!
