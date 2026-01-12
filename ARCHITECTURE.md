# RSSN Project Architecture

## 1. High-Level Overview

The `rssn` library is a modular, multi-layer system designed for high-performance scientific computing. Its architecture is centered around an asynchronous **`ComputeEngine`** that acts as the primary user entry point and orchestrates the various specialized subsystems. This design prioritizes performance, safety, and extensibility.

```mermaid
graph TD
    subgraph User Interaction Layer
        A[External Languages <br>(C, C++, Python)] --> B{FFI Layer};
        U[Rust Library Users] --> C{ComputeEngine <br> (Async Task Orchestrator)};
        B --> C;
    end

    subgraph Core Systems
        C --> D{Core Symbolic System <br> (DAG, Primitives, Simplifier)};
        C --> E[Core Numerical System <br> (Traits, Algorithms)];
        C --> P[Physics Simulation System];
        C --> J[JIT Compilation Engine];
    end

    subgraph Utility Modules
        I[Input & Parsing] --> C;
        D --> O[Output Module <br> (LaTeX, Typst)];
        F[Plugin System] -.-> C;
    end
```

## 2. Key Architectural Pillars

### 2.1. The `ComputeEngine`: Central Orchestrator

The primary interface to the `rssn` library is the `ComputeEngine`. This asynchronous component manages the entire lifecycle of a computation.

- **Task Management**: It accepts tasks (e.g., parsing, simplification, numerical evaluation) and manages their execution.
- **Caching**: The engine features an intelligent caching layer for both parsed expressions and computation results, avoiding redundant work and ensuring high performance.
- **Orchestration**: It acts as the "brain" of the library, delegating tasks to the appropriate subsystem (symbolic, numerical, etc.) and composing their results.

### 2.2. The Core Symbolic System

The heart of `rssn` is its symbolic computation engine, which is utilized by the `ComputeEngine` for all mathematical manipulations.

- **Dual Representation (`Expr` and `DagNode`)**:
  -   `Expr`: A public-facing Abstract Syntax Tree (AST) that is easy to construct and pattern match.
  -   `DagNode`: A private, internal **Directed Acyclic Graph (DAG)** representation. All expressions are transparently converted into a content-addressed DAG managed by a central `DagManager`.
- **`DagManager`**: This singleton service ensures that any structurally identical subexpression is represented by a single node in memory. This provides automatic canonicalization (e.g., `a+b` and `b+a` resolve to the same node) and massive memory savings for complex expressions.
- **Iterative Simplification Engine**: The `simplify_dag` module provides a powerful, stack-safe simplification engine. It operates directly on the DAG using a **bottom-up, iterative, fixpoint** algorithm to apply simplification rules until the expression reaches a stable state.

### 2.3. The Flexible Numerical System

The `numerical` module provides a suite of algorithms for scientific computing, implemented against generic traits to remain flexible. It includes methods for integration, optimization, and solving differential equations.

### 2.4. The FFI Layer and Extensible Blinding

The Foreign Function Interface (FFI) is a critical component for interoperability, designed for safety and performance.

- **Handle-Based System**: Instead of exposing raw Rust pointers, `rssn` uses a `HANDLE_MANAGER` that gives out opaque `usize` integers as handles. This "blinding" technique abstracts away the memory layout and lifetime of Rust objects, preventing memory corruption from the calling language.
- **Extensible for HPC**: The blinding mechanism is extensible. A handle could, in a future implementation, refer to data on a separate device (like a GPU), with the FFI layer managing data locality transparently.
- **Data Serialization**: Complex data is serialized to **JSON** for robust, language-agnostic data exchange.

### 2.5. Other Key Modules

- **`physics`**: A high-level framework for building scientific simulations (FEM, FDM, FVM), designed with data-oriented patterns for performance.
- **`jit`**: A Just-In-Time (JIT) compilation engine that can compile symbolic expressions to native machine code at runtime using backends like `cranelift`.
- **`output`**: A decoupled module for rendering expressions into formats like pretty-printed text, LaTeX, and Typst.
- **`plugins`**: An extensible system allowing third parties to register new functions or simplification rules.
