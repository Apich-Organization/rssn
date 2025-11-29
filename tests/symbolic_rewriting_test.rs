use rssn::symbolic::core::Expr;
use rssn::symbolic::rewriting::{apply_rules_to_normal_form, knuth_bendix, RewriteRule};
use std::sync::Arc;

#[test]
fn test_apply_rules_simple() {
    // Rule: x + 0 -> x
    let p_x = Expr::Pattern("x".to_string());
    let zero = Expr::Constant(0.0);
    
    // Rule: x + 0 -> x
    let rule = RewriteRule {
        lhs: Expr::new_add(p_x.clone(), zero.clone()),
        rhs: p_x.clone(),
    };

    // Expr: (y + 0) + 0
    let y = Expr::new_variable("y");
    let expr = Expr::new_add(Expr::new_add(y.clone(), zero.clone()), zero.clone());

    let result = apply_rules_to_normal_form(&expr, &[rule]);
    
    // Should reduce to y
    assert_eq!(format!("{}", result), "y");
}

#[test]
fn test_apply_rules_associativity() {
    // Rule: op(op(x, y), z) -> op(x, op(y, z))
    // We use NaryList for "op"
    let px = Expr::Pattern("x".to_string());
    let py = Expr::Pattern("y".to_string());
    let pz = Expr::Pattern("z".to_string());
    
    let op_name = "op".to_string();
    
    // lhs: op(op(x, y), z)
    let inner_lhs = Expr::NaryList(op_name.clone(), vec![px.clone(), py.clone()]);
    let lhs = Expr::NaryList(op_name.clone(), vec![inner_lhs, pz.clone()]);
    
    // rhs: op(x, op(y, z))
    let inner_rhs = Expr::NaryList(op_name.clone(), vec![py.clone(), pz.clone()]);
    let rhs = Expr::NaryList(op_name.clone(), vec![px.clone(), inner_rhs]);
    
    let rule = RewriteRule { lhs, rhs };
    
    // Expr: op(op(a, b), c)
    let a = Expr::new_variable("a");
    let b = Expr::new_variable("b");
    let c = Expr::new_variable("c");
    
    let inner_expr = Expr::NaryList(op_name.clone(), vec![a.clone(), b.clone()]);
    let expr = Expr::NaryList(op_name.clone(), vec![inner_expr, c.clone()]);
    
    let result = apply_rules_to_normal_form(&expr, &[rule]);
    
    // Should be op(a, op(b, c))
    let expected_inner = Expr::NaryList(op_name.clone(), vec![b.clone(), c.clone()]);
    let expected = Expr::NaryList(op_name.clone(), vec![a.clone(), expected_inner]);
    
    assert_eq!(format!("{}", result), format!("{}", expected));
}

#[test]
fn test_knuth_bendix_simple() {
    // Simple system:
    // f(x) = g(x)
    // g(x) = h(x)
    // Should derive f(x) = h(x) or similar rules.
    
    let px = Expr::Pattern("x".to_string());
    
    // f(x)
    let fx = Expr::UnaryList("f".to_string(), Arc::new(px.clone()));
    // g(x)
    let gx = Expr::UnaryList("g".to_string(), Arc::new(px.clone()));
    // h(x)
    let hx = Expr::UnaryList("h".to_string(), Arc::new(px.clone()));
    
    let eq1 = Expr::Eq(Arc::new(fx.clone()), Arc::new(gx.clone()));
    let eq2 = Expr::Eq(Arc::new(gx.clone()), Arc::new(hx.clone()));
    
    let equations = vec![eq1, eq2];
    
    let result = knuth_bendix(&equations);
    
    match result {
        Ok(rules) => {
            println!("Generated {} rules", rules.len());
            for r in &rules {
                println!("{} -> {}", r.lhs, r.rhs);
            }
            // We expect at least the input rules oriented, and maybe a critical pair resolved.
            assert!(!rules.is_empty());
        }
        Err(e) => panic!("Knuth-Bendix failed: {}", e),
    }
}

#[test]
fn test_rewrite_rule_serialization() {
    let rule = RewriteRule {
        lhs: Expr::new_variable("a"),
        rhs: Expr::new_variable("b"),
    };
    
    let json = serde_json::to_string(&rule).expect("Serialize");
    let deserialized: RewriteRule = serde_json::from_str(&json).expect("Deserialize");
    
    assert_eq!(format!("{}", rule.lhs), format!("{}", deserialized.lhs));
    assert_eq!(format!("{}", rule.rhs), format!("{}", deserialized.rhs));
}
