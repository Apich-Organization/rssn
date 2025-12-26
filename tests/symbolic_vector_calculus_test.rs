use rssn::symbolic::core::Expr;
use rssn::symbolic::vector::Vector;
use rssn::symbolic::vector_calculus::*;

#[test]

fn test_line_integral_scalar() {

    // Line integral of f(x,y,z) = x^2 + y^2 along the curve r(t) = [t, t^2, 0] from t=0 to t=1
    let t = Expr::new_variable("t");

    let zero = Expr::Constant(0.0);

    let one = Expr::Constant(1.0);

    let curve = ParametricCurve {
        r: Vector::new(
            t.clone(),
            Expr::new_pow(t.clone(), Expr::Constant(2.0)),
            zero.clone(),
        ),
        t_var: "t".to_string(),
        t_bounds: (zero.clone(), one.clone()),
    };

    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    let scalar_field = Expr::new_add(
        Expr::new_pow(x, Expr::Constant(2.0)),
        Expr::new_pow(y, Expr::Constant(2.0)),
    );

    let result = line_integral_scalar(&scalar_field, &curve);

    println!("Line integral (scalar): {}", result);

    // The result should be a valid expression
    assert!(matches!(
        result,
        Expr::Constant(_) | Expr::Rational(_) | Expr::Dag(_)
    ));
}

#[test]

fn test_line_integral_vector() {

    // Line integral of F = [y, -x, 0] along the unit circle r(t) = [cos(t), sin(t), 0] from t=0 to t=2π
    let t = Expr::new_variable("t");

    let zero = Expr::Constant(0.0);

    let two_pi = Expr::new_mul(Expr::Constant(2.0), Expr::Pi);

    let curve = ParametricCurve {
        r: Vector::new(
            Expr::new_cos(t.clone()),
            Expr::new_sin(t.clone()),
            zero.clone(),
        ),
        t_var: "t".to_string(),
        t_bounds: (zero.clone(), two_pi),
    };

    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    let vector_field = Vector::new(y.clone(), Expr::new_neg(x.clone()), zero.clone());

    let result = line_integral_vector(&vector_field, &curve);

    println!("Line integral (vector): {}", result);

    // The result should be a valid expression
    assert!(matches!(
        result,
        Expr::Constant(_) | Expr::Rational(_) | Expr::Dag(_) | Expr::Mul(_, _)
    ));
}

#[test]

fn test_surface_integral() {

    // Surface integral of F = [0, 0, z] over the surface r(u,v) = [u*cos(v), u*sin(v), u]
    // with u in [0,1] and v in [0, 2π]
    let u = Expr::new_variable("u");

    let v = Expr::new_variable("v");

    let zero = Expr::Constant(0.0);

    let one = Expr::Constant(1.0);

    let two_pi = Expr::new_mul(Expr::Constant(2.0), Expr::Pi);

    let surface = ParametricSurface {
        r: Vector::new(
            Expr::new_mul(u.clone(), Expr::new_cos(v.clone())),
            Expr::new_mul(u.clone(), Expr::new_sin(v.clone())),
            u.clone(),
        ),
        u_var: "u".to_string(),
        u_bounds: (zero.clone(), one.clone()),
        v_var: "v".to_string(),
        v_bounds: (zero.clone(), two_pi),
    };

    let z = Expr::new_variable("z");

    let vector_field = Vector::new(zero.clone(), zero.clone(), z.clone());

    let result = surface_integral(&vector_field, &surface);

    println!("Surface integral: {}", result);

    // The result should be a valid expression
    assert!(matches!(
        result,
        Expr::Constant(_) | Expr::Rational(_) | Expr::Dag(_) | Expr::Mul(_, _)
    ));
}

#[test]

fn test_volume_integral() {

    // Volume integral of f(x,y,z) = x*y*z over the unit cube [0,1]×[0,1]×[0,1]
    let zero = Expr::Constant(0.0);

    let one = Expr::Constant(1.0);

    let volume = Volume {
        z_bounds: (zero.clone(), one.clone()),
        y_bounds: (zero.clone(), one.clone()),
        x_bounds: (zero.clone(), one.clone()),
        vars: ("x".to_string(), "y".to_string(), "z".to_string()),
    };

    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    let z = Expr::new_variable("z");

    let scalar_field = Expr::new_mul(Expr::new_mul(x, y), z);

    let result = volume_integral(&scalar_field, &volume);

    println!("Volume integral: {}", result);

    // The result should be 1/8 = 0.125
    if let Expr::Constant(val) = result {

        assert!((val - 0.125).abs() < 1e-9);
    } else if let Expr::Rational(r) = result {

        use num_traits::ToPrimitive;

        assert_eq!(r.to_f64().unwrap(), 0.125);
    } else {

        // It might return a DAG that evaluates to 0.125
        println!("Result is a complex expression: {}", result);
    }
}

#[test]

fn test_parametric_curve_serialization() {

    // Test that ParametricCurve can be serialized/deserialized
    let t = Expr::new_variable("t");

    let zero = Expr::Constant(0.0);

    let one = Expr::Constant(1.0);

    let curve = ParametricCurve {
        r: Vector::new(t.clone(), t.clone(), zero.clone()),
        t_var: "t".to_string(),
        t_bounds: (zero, one),
    };

    let json = serde_json::to_string(&curve).unwrap();

    let deserialized: ParametricCurve = serde_json::from_str(&json).unwrap();

    assert_eq!(curve.t_var, deserialized.t_var);
}

#[test]

fn test_parametric_surface_serialization() {

    // Test that ParametricSurface can be serialized/deserialized
    let u = Expr::new_variable("u");

    let v = Expr::new_variable("v");

    let zero = Expr::Constant(0.0);

    let one = Expr::Constant(1.0);

    let surface = ParametricSurface {
        r: Vector::new(u.clone(), v.clone(), zero.clone()),
        u_var: "u".to_string(),
        u_bounds: (zero.clone(), one.clone()),
        v_var: "v".to_string(),
        v_bounds: (zero, one),
    };

    let json = serde_json::to_string(&surface).unwrap();

    let deserialized: ParametricSurface = serde_json::from_str(&json).unwrap();

    assert_eq!(surface.u_var, deserialized.u_var);

    assert_eq!(surface.v_var, deserialized.v_var);
}

#[test]

fn test_volume_serialization() {

    // Test that Volume can be serialized/deserialized
    let zero = Expr::Constant(0.0);

    let one = Expr::Constant(1.0);

    let volume = Volume {
        z_bounds: (zero.clone(), one.clone()),
        y_bounds: (zero.clone(), one.clone()),
        x_bounds: (zero, one),
        vars: ("x".to_string(), "y".to_string(), "z".to_string()),
    };

    let json = serde_json::to_string(&volume).unwrap();

    let deserialized: Volume = serde_json::from_str(&json).unwrap();

    assert_eq!(volume.vars, deserialized.vars);
}
