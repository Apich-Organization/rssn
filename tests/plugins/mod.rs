use std::collections::HashMap;
use rssn::plugins::manager::PluginManager;
use rssn::plugins::plugin_c::{Plugin, PluginError, PluginHealth, PluginErrorKind};
use rssn::symbolic::core::{Expr, DagOp};
use std::sync::Arc;

struct MockPlugin;

impl Plugin for MockPlugin {
    fn name(&self) -> &'static str {
        "mock_plugin"
    }

    fn description(&self) -> &'static str {
        "A mock plugin for testing"
    }

    fn api_version(&self) -> &'static str {
        "0.1.19"
    }

    fn on_load(&self) -> Result<(), PluginError> {
        Ok(())
    }

    fn execute(&self, command: &str, args: &Expr) -> Result<Expr, PluginError> {
        match command {
            "identity" => Ok(args.clone()),
            "double" => {
                if let DagOp::Constant(c) = args.op() {
                    Ok(Expr::new_constant(c.into_inner() * 2.0))
                } else {
                    Err(PluginError::new(PluginErrorKind::ExecutionFailed, "Expected constant"))
                }
            }
            _ => Err(PluginError::new(PluginErrorKind::NotFound, "Command not found")),
        }
    }

    fn health_check(&self) -> PluginHealth {
        PluginHealth::Ok
    }

    fn metadata(&self) -> HashMap<String, String> {
        let mut meta = HashMap::new();
        meta.insert("author".to_string(), "tester".to_string());
        meta
    }
}

#[test]
fn test_plugin_manager_registration() {
    let manager = PluginManager::empty();

    manager.register_plugin(Box::new(MockPlugin));
    
    let names = manager.get_loaded_plugin_names();
    assert!(names.contains(&"mock_plugin".to_string()));
}

#[test]
fn test_plugin_execution() {
    let manager = PluginManager::empty();

    manager.register_plugin(Box::new(MockPlugin));
    
    let args = Expr::new_constant(5.0);
    let res = manager.execute_plugin("mock_plugin", "double", &args).unwrap();
    
    if let DagOp::Constant(c) = res.op() {
        assert_eq!(c.into_inner(), 10.0);
    } else {
        panic!("Execution failed");
    }
}

#[test]
fn test_plugin_metadata() {
    let manager = PluginManager::empty();

    manager.register_plugin(Box::new(MockPlugin));
    
    let metadata = manager.get_plugin_metadata();
    assert_eq!(metadata["mock_plugin"]["author"], "tester");
}

#[test]
fn test_plugin_unload() {
    let manager = PluginManager::empty();

    manager.register_plugin(Box::new(MockPlugin));
    assert!(manager.get_loaded_plugin_names().contains(&"mock_plugin".to_string()));
    
    assert!(manager.unload_plugin("mock_plugin"));
    assert!(!manager.get_loaded_plugin_names().contains(&"mock_plugin".to_string()));
}
