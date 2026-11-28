target triple = "x86_64-pc-windows-msvc"

declare double @llvm.pow.f64(double, double)
declare double @consume(double, double, double, double, double, double)

define double @share_pow_budget(double %base, double %e0, double %e1,
                                double %e2, double %e3, i32 %n,
                                double %a0, double %a1, double %a2,
                                double %a3, double %a4) {
entry:
  br label %header

header:
  %i = phi i32 [ 0, %entry ], [ %next, %join ]
  %parity = and i32 %i, 1
  %even = icmp eq i32 %parity, 0
  br i1 %even, label %left, label %right

left:
  %l0 = call double @llvm.pow.f64(double %base, double %e0)
  %l1 = call double @llvm.pow.f64(double %base, double %e1)
  %l01 = fmul double %l0, %l1
  %l2 = call double @llvm.pow.f64(double %base, double %e2)
  %l012 = fmul double %l01, %l2
  %l3 = call double @llvm.pow.f64(double %base, double %e3)
  %left_value = fmul double %l012, %l3
  br label %join

right:
  %r0 = call double @llvm.pow.f64(double %base, double %e0)
  %r1 = call double @llvm.pow.f64(double %base, double %e1)
  %r01 = fmul double %r0, %r1
  %r2 = call double @llvm.pow.f64(double %base, double %e2)
  %r012 = fmul double %r01, %r2
  %r3 = call double @llvm.pow.f64(double %base, double %e3)
  %right_value = fmul double %r012, %r3
  br label %join

join:
  %v = phi double [ %left_value, %left ], [ %right_value, %right ]
  %out = call double @consume(double %v, double %a0, double %a1,
                             double %a2, double %a3, double %a4)
  %next = add i32 %i, 1
  %more = icmp slt i32 %next, %n
  br i1 %more, label %header, label %exit, !prof !0

exit:
  ret double %out
}

!0 = !{!"branch_weights", i32 100, i32 1}
