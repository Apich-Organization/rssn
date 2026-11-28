target triple = "x86_64-unknown-linux-gnu"

declare double @llvm.exp.f64(double)
declare double @llvm.fabs.f64(double)
declare i1 @llvm.is.fpclass.f64(double, i32 immarg)
declare void @observe()

define double @endpoint_node(double %weight, double %arg, i1 %left, i32 %n) {
entry:
  br label %loop

loop:
  %i = phi i32 [ 0, %entry ], [ %next, %join ]
  %e = call double @llvm.exp.f64(double %arg)
  %plus = fadd double %e, 1.000000e+00
  %distance = fdiv double 2.000000e+00, %plus
  %positive = fcmp ugt double %distance, 0.000000e+00
  %abs = call double @llvm.fabs.f64(double %weight)
  %finite = fcmp one double %abs, +inf
  %valid = and i1 %positive, %finite
  br i1 %valid, label %lower, label %join

lower:
  %scaled = fmul double %distance, 5.000000e-01
  %reflected = fsub double 1.000000e+00, %scaled
  %node = select i1 %left, double %scaled, double %reflected
  %above = fcmp ugt double %node, 0.000000e+00
  br i1 %above, label %upper, label %join

upper:
  %outside = fcmp oge double %node, 1.000000e+00
  %zero = fcmp oeq double %weight, 0.000000e+00
  %reject = or i1 %outside, %zero
  br i1 %reject, label %join, label %interior

interior:
  %term = fmul double %weight, %node
  br label %join

join:
  %result = phi double [ 0.000000e+00, %loop ], [ 0.000000e+00, %lower ],
                       [ 0.000000e+00, %upper ], [ %term, %interior ]
  %next = add i32 %i, 1
  %more = icmp slt i32 %next, %n
  br i1 %more, label %loop, label %exit

exit:
  ret double %result
}

define i32 @mapped_weight(double %factor, double %seed) {
entry:
  br label %loop

loop:
  %i = phi i32 [ 0, %entry ], [ %next, %join ]
  %count = phi i32 [ 0, %entry ], [ %result, %join ]
  %index = uitofp i32 %i to double
  %arg = fadd double %seed, %index
  %node = call double @llvm.exp.f64(double %arg)
  %weight = fmul double %node, %factor
  %invalid = call i1 @llvm.is.fpclass.f64(double %weight, i32 615)
  br i1 %invalid, label %join, label %valid

valid:
  %count.valid = add i32 %count, 1
  br label %join

join:
  %result = phi i32 [ %count, %loop ], [ %count.valid, %valid ]
  %next = add i32 %i, 1
  %again = icmp ult i32 %next, 3
  br i1 %again, label %loop, label %exit

exit:
  ret i32 %result
}

define i32 @carried_magnitude(double %first, double %candidate, i32 %n) {
entry:
  %initial = call double @llvm.fabs.f64(double %first)
  br label %header

header:
  %i = phi i32 [ 0, %entry ], [ %next, %latch ]
  %previous = phi double [ %initial, %entry ], [ %magnitude, %latch ]
  %abs = call double @llvm.fabs.f64(double %candidate)
  %take = fcmp ogt double %abs, %previous
  %magnitude = select i1 %take, double %abs, double %previous
  %invalid.class = call i1 @llvm.is.fpclass.f64(double %magnitude, i32 579)
  br i1 %invalid.class, label %invalid, label %latch

latch:
  %next = add i32 %i, 1
  %again = icmp ult i32 %next, %n
  br i1 %again, label %header, label %exit

invalid:
  ret i32 0

exit:
  ret i32 %next
}
