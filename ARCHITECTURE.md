# RSSN Project Architecture

## 1. High-Level Overview

The `rssn` library is designed as a modular, multi-layer system that separates the core mathematical engine from its various interfaces and applications. This architecture prioritizes performance, safety, and extensibility.

```mermaid
graph TD
    A[External Languages <br>(C, C++, Python, Fortran)] --> B{FFI Layer w/ Blinding};
    B --> C[High-Level APIs <br> symbolic, numerical, physics];
    C --> D{Core Symbolic System <br> (DAG, Primitives)};
    C --> E[Core Numerical System <br> (Traits, Algorithms)];
    F[Plugin System] -.-> C;
    C --> G[Output Module <br> (LaTeX, Typst, Pretty-Print)];
```

## 2. Key Architectural Pillars

### 2.1. The Core Symbolic System

The heart of `rssn` is its symbolic computation engine, designed around performance and canonical representation.

- **Dual Representation (`Expr` and `DagNode`)**: 
  -   `Expr`: A public-facing Abstract Syntax Tree (AST) that is easy to construct and pattern match. It represents the logical structure of a mathematical expression.
  -   `DagNode`: A private, internal **Directed Acyclic Graph (DAG)** representation. All expressions are transparently converted into a content-addressed DAG managed by a central `DagManager`.
- **`DagManager`**: This singleton service ensures that any structurally identical subexpression is represented by a single node in memory. This provides automatic canonicalization (e.g., `a+b` and `b+a` can resolve to the same node) and massive memory savings and performance gains for large, complex expressions.
- **Iterative Simplification Engine**: The `simplify_dag` module provides a powerful, stack-safe simplification engine. It operates directly on the DAG and uses a **bottom-up, iterative, fixpoint** algorithm to apply simplification rules until the expression reaches a stable state. This design avoids the stack overflows common in recursive CAS engines and ensures that rules are applied exhaustively.
- **Advanced Algebraic Modules**: High-level mathematical concepts are built on top of the core. The `grobner` module, for example, provides functionality to compute Gröbner bases, which is then used by `cas_foundations` to offer powerful simplification and normalization of expressions with respect to a set of polynomial side-relations.

### 2.2. The Flexible Numerical System

The `numerical` module provides a suite of algorithms for scientific computing. Its architecture is based on Rust's trait system to remain flexible and generic.

- **Trait-Based**: Algorithms are implemented against traits (e.g., `Fn`), allowing them to operate on a wide variety of numeric types and functions.
- **Comprehensive**: It includes a rich collection of common numerical methods for integration, optimization, linear algebra, and statistics.

### 2.3. The Parallel-Ready Physics Simulation System

The `physics` module is a high-level framework for building scientific simulations (e.g., FEM, FDM, FVM solvers).

- **Data-Oriented**: The design favors data-oriented patterns that are amenable to high-performance computing.
- **Extensible for Parallelism**: While the initial implementation may be single-threaded, the architecture is designed to be extended for parallel computation using libraries like `rayon` or `MPI`, a crucial requirement for large-scale simulations.

### 2.4. Pluggable Subsystems

A key design goal is extensibility without modification of the core library. The `plugins` system is envisioned to allow third parties to dynamically extend the functionality of the CAS.

- **Dynamic Dispatch**: A `PluginManager` can load and manage plugins that register new functions or simplification rules.
- **Trait-Based Interface**: Plugins would implement a `Plugin` trait, providing a clear and stable API for defining their behavior.

### 2.5. The FFI Layer and Extensible Blinding

The FFI (Foreign Function Interface) layer is a critical component for interoperability, designed for safety, stability, and high performance.

- **Modern Handle-Based System**: Instead of exposing raw Rust pointers, `rssn` uses a `HANDLE_MANAGER` that gives out `usize` integers as handles. This abstracts away the memory layout and lifetime of Rust objects, preventing the calling language from causing memory corruption.
- **FFI Blinding**: This handle-based approach is a form of "blinding." The internal data structures are completely opaque to the foreign caller. This is a deliberate design for both stability and high-performance computing (HPC).
- **Extensible for HPC**: The blinding mechanism is extensible. A `usize` handle could, in a future implementation, refer to data residing on a separate device (like a GPU). The FFI functions in Rust would manage the data locality and dispatch computations to the appropriate device, all while the foreign caller remains unaware of these details, simply holding an integer ID.
- **Data Serialization**: Complex data that needs to be passed by value is serialized to **JSON**, providing a language-agnostic and robust data exchange format.
- **Error Handling**: A thread-local `LAST_ERROR` mechanism is used to provide clear, retrievable error messages without complicating the C-compatible function signatures.

### 2.6. The Output Module

The `output` module is responsible for rendering expressions into human-readable formats. It is decoupled from the core symbolic engine.

- **Multiple Backends**: It supports multiple output formats, including pretty-printed console text, LaTeX, and Typst.
- **Extensible**: The design allows for new formatters to be added easily by implementing a common rendering trait for the `Expr` type.
