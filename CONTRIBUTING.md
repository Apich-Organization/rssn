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

## 🚀 Roadmap: C++ Adapter for FFI

We aim to provide a first-class experience for C++ users. While the C-style FFI is functional, it is not idiomatic for C++ developers. We are looking for contributors to help create a modern, header-only C++ wrapper library.

### Goal

The goal is to create a `rssn.hpp` that provides a clean, object-oriented C++ interface over the raw C FFI.

### Core Design: `RssnExpr` Class

The central piece of the adapter would be a `RssnExpr` class that wraps the `*mut Expr` handle.

```cpp
#include <string>
#include <memory>
#include <optional> // C++17
#include "nlohmann/json.hpp" // Example JSON library

// Forward declarations of the C FFI functions
extern "C" {
    struct Expr;
    Expr* expr_from_json(const char* json_ptr);
    void expr_free(Expr* handle);
    char* expr_to_string(Expr* handle);
    // ... other functions
}

class RssnExpr {
private:
    Expr* handle_ = nullptr;

public:
    // Constructor is private to force creation via factory methods
    explicit RssnExpr(Expr* handle) : handle_(handle) {}

    // RAII: Destructor to automatically free the Rust object
    ~RssnExpr() {
        if (handle_) {
            expr_free(handle_);
        }
    }

    // Disable copy constructor and assignment to prevent double-freeing
    RssnExpr(const RssnExpr&) = delete;
    RssnExpr& operator=(const RssnExpr&) = delete;

    // Enable move semantics
    RssnExpr(RssnExpr&& other) noexcept : handle_(other.handle_) {
        other.handle_ = nullptr; // Prevent the moved-from object from freeing the handle
    }
    RssnExpr& operator=(RssnExpr&& other) noexcept {
        if (this != &other) {
            if (handle_) {
                expr_free(handle_);
            }
            handle_ = other.handle_;
            other.handle_ = nullptr;
        }
        return *this;
    }

    // Static factory method
    static RssnExpr fromJson(const std::string& json_str) {
        Expr* handle = expr_from_json(json_str.c_str());
        if (!handle) {
            throw std::runtime_error("Failed to create Expr from JSON");
        }
        return RssnExpr(handle);
    }

    // Method to wrap an FFI function
    RssnExpr simplify() {
        Expr* new_handle = expr_simplify(handle_);
        if (!new_handle) {
            throw std::runtime_error("Failed to simplify expression");
        }
        return RssnExpr(new_handle);
    }

    // Method to wrap a function that returns a string
    std::string toString() {
        char* c_str = expr_to_string(handle_);
        std::string str(c_str);
        free_string(c_str); // Remember to free the string from Rust
        return str;
    }
};
```

### Key Responsibilities for the Contributor

1.  **RAII and Memory Management**: Implement robust RAII to ensure that no memory is leaked. The C++ class should handle all calls to `_free` functions automatically.
2.  **JSON Integration**: The C++ adapter will need a dependency on a JSON library (like `nlohmann/json`, `RapidJSON`, etc.) to construct the JSON strings required by the FFI and parse the JSON returned by it.
3.  **Error Handling**: For FFI functions that return a result via JSON (like `expr_unify_expression`), the C++ wrapper should parse the JSON and translate it into idiomatic C++ error handling, such as returning a `std::optional` or throwing an exception.
4.  **API Design**: Design an intuitive C++ API that hides the complexity of the underlying C FFI. This includes wrapping functions, overloading operators (e.g., `+`, `*` for expressions), and providing clear documentation.
5.  **Build System Integration**: Provide instructions or a simple CMake/Makefile example to show how a C++ project can easily include and link against the `rssn` dynamic library and use the header-only adapter.

If you are interested in leading this effort, please open an issue on GitHub to discuss the design further!

---

## 🙏 Acknowledgements

Contributors are credited in release notes and on the GitHub page.
We value every contribution, from fixing typos to implementing new solvers.

Thank you for making **rssn** better!

