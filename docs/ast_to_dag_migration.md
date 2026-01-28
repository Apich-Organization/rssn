# AST to DAG Migration Plan for symbolic/core.rs

## Problem Statement

Currently, `Expr` has a dual representation:
1. **Old AST variants**: `Expr::Add`, `Expr::Sin`, `Expr::Variable`, etc. (100+ variants)
2. **New DAG variant**: `Expr::Dag(Arc<DagNode>)` with `DagOp` enum

This creates several issues:
- Code duplication in implementations (`Clone`, `PartialEq`, `Hash`, `Display`, etc.)
- Serialization problems (`Expr::Dag` is marked `skip_serializing`)
- Inconsistent API (some functions return AST, others return DAG)
- Performance issues (no automatic expression sharing)

## Goal

Transition to a DAG-first architecture where:
- All `Expr` values are internally `Expr::Dag`
- Constructor functions (`new_add`, `new_sin`, etc.) create DAG nodes
- Old AST variants are deprecated but kept for backward compatibility
- Serialization works correctly
- No breaking changes to public API

## Migration Strategy

### Phase 1: Add Compatibility Layer (Non-Breaking)

1. **Mark old variants as deprecated**:
   ```rust
   #[deprecated(since = "0.2.0", note = "Use Expr::new_add() instead")]
   Add(Arc<Expr>, Arc<Expr>),
   ```

2. **Ensure all constructors use DAG**:
   - Already done for most: `new_variable`, `new_add`, etc.
   - Verify all return `Expr::Dag`

3. **Add conversion methods**:
   ```rust
   impl Expr {
       /// Converts any Expr variant to DAG representation
       pub fn to_dag(&self) -> Expr {
           match self {
               Expr::Dag(_) => self.clone(),
               _ => {
                   // Convert AST to DAG
                   DAG_MANAGER.get_or_create(self).unwrap()
               }
           }
       }
       
       /// Ensures this expression is in DAG form
       pub fn normalize(&mut self) {
           *self = self.to_dag();
       }
   }
   ```

### Phase 2: Update Serialization (Breaking for Dag variant)

1. **Make `DagNode` serializable**:
   ```rust
   #[derive(Debug, Clone, Serialize, Deserialize)]
   pub struct DagNode {
       pub op: DagOp,
       #[serde(serialize_with = "serialize_dag_children")]
       #[serde(deserialize_with = "deserialize_dag_children")]
       pub children: Vec<Arc<DagNode>>,
       #[serde(skip)]
       pub hash: u64,
   }
   ```

2. **Custom serialization for Dag variant**:
   ```rust
   // Remove skip_serializing from Expr::Dag
   #[derive(serde::Serialize, serde::Deserialize)]
   pub enum Expr {
       // ...
       Dag(Arc<DagNode>),  // Now serializable
       // ...
   }
   ```

3. **Backward compatibility**:
   - Deserialize old AST format
   - Convert to DAG on load
   - Always serialize as DAG

### Phase 3: Deprecate AST Variants (Breaking)

1. **Move AST variants to separate enum**:
   ```rust
   #[deprecated]
   pub enum ExprAst {
       Add(Arc<Expr>, Arc<Expr>),
       Sin(Arc<Expr>),
       // ... all old variants
   }
   
   pub enum Expr {
       // Only DAG variant remains
       Dag(Arc<DagNode>),
   }
   ```

2. **Provide conversion**:
   ```rust
   impl From<ExprAst> for Expr {
       fn from(ast: ExprAst) -> Self {
           // Convert to DAG
       }
   }
   ```

### Phase 4: Remove AST Variants (Major Breaking)

1. **Remove `ExprAst` enum entirely**
2. **Only `Expr::Dag` remains**
3. **Update all dependent code**

## Implementation Plan (Phase 1 - Non-Breaking)

### Step 1: Audit Constructor Functions

Verify all public constructors return `Expr::Dag`:
- [x] `new_variable` ✓
- [x] `new_add`, `new_sub`, `new_mul`, `new_div` ✓
- [ ] Check all ~100 constructors

### Step 2: Add Conversion Utilities

```rust
impl Expr {
    /// Converts this expression to DAG form if not already
    pub fn to_dag(&self) -> Result<Expr, String> {
        match self {
            Expr::Dag(_) => Ok(self.clone()),
            _ => {
                let dag_node = DAG_MANAGER.get_or_create(self)?;
                Ok(Expr::Dag(dag_node))
            }
        }
    }
    
    /// Checks if this expression is in DAG form
    pub fn is_dag(&self) -> bool {
        matches!(self, Expr::Dag(_))
    }
}
```

### Step 3: Update DagManager

Ensure `DagManager::get_or_create` handles all `Expr` variants:
```rust
impl DagManager {
    pub fn get_or_create(&self, expr: &Expr) -> Result<Arc<DagNode>, String> {
        match expr {
            Expr::Dag(node) => Ok(node.clone()),
            _ => {
                let op = expr.to_dag_op_internal()?;
                let children = expr.children()
                    .into_iter()
                    .map(|child| self.get_or_create(&child))
                    .collect::<Result<Vec<_>, _>>()?;
                self.get_or_create_normalized(op, children)
            }
        }
    }
}
```

### Step 4: Make Serialization Work

Option A: Serialize DAG as AST (backward compatible)
```rust
impl Serialize for Expr {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        match self {
            Expr::Dag(node) => {
                // Convert to AST for serialization
                let ast = node.to_expr().map_err(serde::ser::Error::custom)?;
                ast.serialize(serializer)
            }
            _ => {
                // Serialize AST directly
                // ... existing code
            }
        }
    }
}
```

Option B: Serialize DAG natively (new format)
```rust
// Make DagNode serializable
#[derive(Serialize, Deserialize)]
pub struct DagNode {
    pub op: DagOp,
    pub children: Vec<Arc<DagNode>>,
    #[serde(skip)]
    pub hash: u64,
}
```

**Recommendation**: Use Option A for now (backward compatible), migrate to Option B in Phase 2.

### Step 5: Add Tests

```rust
#[cfg(test)]
mod dag_migration_tests {
    use super::*;
    
    #[test]
    fn test_ast_to_dag_conversion() {
        let ast = Expr::Add(
            Arc::new(Expr::Variable("x".to_string())),
            Arc::new(Expr::new_constant(1.0))
        );
        let dag = ast.to_dag().unwrap();
        assert!(dag.is_dag());
    }
    
    #[test]
    fn test_dag_serialization() {
        let expr = Expr::new_add(
            Expr::new_variable("x"),
            Expr::new_constant(1.0)
        );
        let json = serde_json::to_string(&expr).unwrap();
        let deserialized: Expr = serde_json::from_str(&json).unwrap();
        assert_eq!(expr, deserialized);
    }
}
```

## Risks and Mitigation

### Risk 1: Breaking Existing Code
**Mitigation**: Phase 1 is non-breaking, only adds deprecation warnings

### Risk 2: Serialization Incompatibility
**Mitigation**: Support both formats during transition

### Risk 3: Performance Regression
**Mitigation**: Benchmark before/after, DAG should be faster

### Risk 4: Complex Codebase
**Mitigation**: Incremental migration, extensive testing

## Timeline

- **Phase 1** (Current): 1-2 weeks - Add compatibility layer
- **Phase 2** (v0.2.0): 2-3 weeks - Update serialization
- **Phase 3** (v0.3.0): 1-2 weeks - Deprecate AST
- **Phase 4** (v1.0.0): 1 week - Remove AST

## Success Criteria

- [ ] All constructors return `Expr::Dag`
- [ ] Serialization works for both AST and DAG
- [ ] All tests pass
- [ ] No performance regression
- [ ] Documentation updated
- [ ] Migration guide written

## Current Status

**Phase 1 - In Progress**
- Constructors mostly use DAG (need audit)
- Serialization broken for DAG variant
- Need to add conversion utilities
- Need to update DagManager

**Next Steps**:
1. Add `to_dag()` and `is_dag()` methods
2. Fix serialization (Option A)
3. Audit all constructors
4. Add comprehensive tests
5. Update documentation
