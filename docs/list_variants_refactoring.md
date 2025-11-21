# List Variants Refactoring Summary

## Overview
This document summarizes the refactoring work done to add list variants to the `Expr` enum in `src/symbolic/core.rs` as part of the AST to DAG migration strategy.

## Changes Made

### 1. New Expr Variants Added

#### N-ary Operation Variants
- **`AddList(Vec<Expr>)`**: N-ary addition for efficient sum operations
- **`MulList(Vec<Expr>)`**: N-ary multiplication for efficient product operations

These variants allow representing `a + b + c + d` or `a * b * c * d` as single operations instead of nested binary operations, improving efficiency and simplifying the DAG representation.

#### Generic/Dynamic Operation Variants
- **`UnaryList(String, Arc<Expr>)`**: Generic unary operations identified by name
- **`BinaryList(String, Arc<Expr>, Arc<Expr>)`**: Generic binary operations identified by name
- **`NaryList(String, Vec<Expr>)`**: Generic n-ary operations identified by name

These variants enable runtime registration of custom operations without modifying the core `Expr` enum, supporting plugin systems and domain-specific extensions.

### 2. Dynamic Operation Registry

Added a global registry system for managing dynamic operations:

```rust
pub struct DynamicOpProperties {
    pub name: String,
    pub description: String,
    pub is_associative: bool,
    pub is_commutative: bool,
}

pub static ref DYNAMIC_OP_REGISTRY: RwLock<HashMap<String, DynamicOpProperties>>;
```

Functions:
- `register_dynamic_op(name: &str, props: DynamicOpProperties)`: Register a new operation
- `get_dynamic_op_properties(name: &str) -> Option<DynamicOpProperties>`: Retrieve operation properties

### 3. Updated Trait Implementations

All `Expr` trait implementations were updated to handle the new variants:

- **`Clone`**: Properly clones all list variants
- **`fmt::Display`**: Provides readable string representation
- **`variant_order`**: Assigns consistent ordering for canonicalization
- **`pre_order_walk`**: Traverses children in pre-order
- **`post_order_walk`**: Traverses children in post-order
- **`in_order_walk`**: Traverses children in in-order
- **`get_children_internal`**: Returns child expressions
- **`normalize`**: Flattens nested lists and sorts children (respecting commutativity)
- **`to_dag_op_internal`**: Maps to corresponding `DagOp` variants

### 4. DagOp Enum Updates

Added corresponding `DagOp` variants:
- `DagOp::UnaryList(String)`
- `DagOp::BinaryList(String)`
- `DagOp::NaryList(String)`

Note: `AddList` and `MulList` map to existing `DagOp::Add` and `DagOp::Mul`.

### 5. DagNode::to_expr Updates

Updated to handle n-ary operations:
- Returns `Expr::AddList` when `DagOp::Add` has != 2 children
- Returns `Expr::MulList` when `DagOp::Mul` has != 2 children
- Properly converts generic list variants

### 6. Simplify.rs Enhancements

Added simplification rules for list variants:

**AddList simplification**:
- Flattens nested `AddList` expressions
- Filters out zeros
- Combines constant terms
- Returns simplified form (single term, constant, or reduced list)

**MulList simplification**:
- Flattens nested `MulList` expressions
- Short-circuits on zero
- Filters out ones
- Combines constant factors
- Returns simplified form

**Generic list variants**:
- Recursively simplifies children
- Preserves operation name

### 7. Documentation

Added comprehensive documentation including:
- Detailed docstrings for all new variants
- Usage examples for each variant
- Explanation of the design rationale
- Integration with the DAG migration strategy

## Benefits

1. **Performance**: N-ary operations reduce tree depth and improve cache locality
2. **Extensibility**: Dynamic operations allow plugins without core modifications
3. **Backward Compatibility**: Existing binary operations remain unchanged
4. **Future-Proofing**: Supports ongoing AST to DAG migration
5. **Flexibility**: Operations can be registered at runtime with custom properties

## Testing

All changes compile successfully and maintain backward compatibility with existing code.

## Next Steps

1. Add comprehensive unit tests for list variants
2. Add property-based tests for simplification rules
3. Create examples demonstrating dynamic operation registration
4. Add benchmarks comparing binary vs n-ary operations
5. Update simplify_dag.rs to leverage n-ary operations
6. Consider adding more sophisticated simplification rules for list variants
