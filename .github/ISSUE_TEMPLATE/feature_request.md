---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

### 🚀 Is your feature request related to a problem?

A clear and concise description of the problem or limitation you are currently facing without this feature.

*Example: "I am trying to solve large systems of linear equations but the current solver lacks support for multi-threading, leading to slow computation times."*

### 💡 Describe the solution you'd like

A clear and concise description of what you want to happen. Ideally, describe how the new API would look in Rust code.

*Which module would this feature belong to? (e.g., `kernels::krylov`, `rules::poly`, `sim`, `io`)*

```rust
// Example of how the new API/functionality would be used:

use rssn::prelude::*;

fn main() {
    let result = new_function(input_data);
    // ...
}
```

### 📋 Describe alternatives you've considered

A clear and concise description of any alternative solutions or workarounds you have already tried or considered.

*Example: "I currently use the `nalgebra` crate for this specific task, but integrating it adds complexity and overhead that `rssn` could eliminate if the feature were internal."*

### ➕ Does this require new dependencies or features?

If this feature requires adding a new dependency to `Cargo.toml` or enabling a new feature flag, please list them here.

### ❓ Additional context

Add any other context, links to research papers, or images about the feature request here.
