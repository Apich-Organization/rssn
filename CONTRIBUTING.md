# Contributing to rssn

First of all, thank you for considering contributing to **rssn**! We are thrilled to have you. This project aims to become a next-generation scientific computing ecosystem in Rust, and every contribution, no matter how small, helps us get there.

We are building a friendly, open community and are committed to helping you get started.

---

## How Can I Contribute?

**We welcome contributions across all areas of the project!** Whether you are a seasoned Rust developer, a numerical methods expert, a student learning symbolic math, or just someone who wants to fix a typo, there is a place for you here.

Good places to start include:

*   **Symbolic Engine**: Help expand our core CAS capabilities in `src/symbolic`. This could mean adding new simplification rules, improving polynomial algebra, or implementing new calculus functions.
*   **Numerical Methods**: Add new algorithms or improve existing ones in `src/numerical`, such as optimizers, ODE solvers, or statistical tools.
*   **Physics Simulations**: Contribute new models or solvers to the `src/physics` module.
*   **Documentation**: Good documentation is as important as good code. Improving explanations, adding examples, and clarifying concepts is invaluable.
*   **Testing**: Help us improve test coverage, add new regression tests for bugs, or create new benchmarks.

If you have an idea, open an issue and we would be happy to discuss it!

---

## Getting Started: It's Easy!

We have designed the project to be as easy to set up as possible. You do not need any complex dependencies or special environment configuration.

1.  **Install Rust**: If you don't have it already, install the Rust toolchain via [rustup](https://rustup.rs/).

2.  **Fork and Clone**: Fork the repository on GitHub and clone it to your local machine.
    ```bash
    git clone https://github.com/YOUR-USERNAME/rssn.git
    cd rssn
    ```

3.  **Build and Test**: Check that everything is working correctly by running the test suite.
    ```bash
    cargo test --all
    ```

That's it! You are now ready to start contributing.

---

## 🔧 Development Workflow

1. **Create a feature branch**:

   ```bash
   git checkout -b feature/my-new-feature
   ```

2. **Write your code**: Make your changes, and please add tests for any new functionality or bug fixes!

3. **Code style & formatting**:

   * All code must compile with **zero warnings** on the latest stable Rust.
   * Additional **lint rules** are configured in `lib.rs` and must be respected.
   * Before creating a pull request, please always run:

     ```bash
     cargo fmt --all
     cargo clippy --all-targets -- -D warnings
     cargo test --all
     ```

4. **Commit and Push**:

   Commit your changes with a clear message and push them to your fork.

   ```bash
   git commit -m "feat(symbolic): add new integration method"
   ```

5. **Open a Pull Request**: Open a PR against the `main` branch of the `Apich-Organization/rssn` repository. Please provide a clear description of your changes.

---

## 🚀 Roadmap: C++ Adapter for FFI

We aim to provide a first-class experience for C++ users. While the C-style FFI is functional, it is not idiomatic for C++ developers. We are looking for contributors to help create a modern, header-only C++ wrapper library.

### Goal

The goal is to create a `rssn.hpp` that provides a clean, object-oriented C++ interface over the raw C FFI.

<details>
<summary>Click for C++ Adapter Design Details</summary>

### Core Design: `RssnExpr` Class

The central piece of the adapter would be a `RssnExpr` class that wraps the `*mut Expr` handle.

```cpp
#include <string>
#include <memory>
#include <stdexcept> 

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

1.  **RAII and Memory Management**: Implement robust RAII to ensure that no memory is leaked.
2.  **JSON Integration**: The C++ adapter will need a dependency on a JSON library (like `nlohmann/json`) to construct and parse JSON strings.
3.  **Error Handling**: Translate FFI error-reporting (e.g., via JSON) into idiomatic C++ exceptions or error codes.
4.  **API Design**: Design an intuitive C++ API that hides the complexity of the underlying C FFI.
5.  **Build System Integration**: Provide a simple CMake/Makefile example to show how a C++ project can use the adapter.

If you are interested in leading this effort, please open an issue on GitHub to discuss the design further!

</details>

---

## 🙏 Acknowledgements

Contributors are credited in release notes and on the GitHub page.
We value every contribution, from fixing typos to implementing new solvers.

Thank you for making **rssn** better!
