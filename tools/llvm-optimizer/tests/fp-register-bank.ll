target triple = "x86_64-pc-windows-msvc"

declare double @llvm.pow.f64(double, double)
declare double @observe(double, double, double, double, double, double,
                        double, double, double, double, double, double)

define double @bank_exponents(double %base, double %e0, double %a0, i32 %n,
  double noundef %e1, double %a1, double %a2, double %a3, double %a4,
  double %a5, double %a6, double %a7, double %a8, double %a9) {
entry:
  %exponent = fneg double %e0
  br label %header

header:
  %i = phi i32 [ 0, %entry ], [ %next, %join ]
  %parity = and i32 %i, 1
  %even = icmp eq i32 %parity, 0
  br i1 %even, label %positive, label %negative

positive:
  %a = call double @llvm.pow.f64(double %base, double %exponent)
  %b = call double @llvm.pow.f64(double %base, double %e1)
  br label %join

negative:
  %neg = fneg double %base
  %c = call double @llvm.pow.f64(double %neg, double %exponent)
  %d = call double @llvm.pow.f64(double %neg, double %e1)
  br label %join

join:
  %low = phi double [ %a, %positive ], [ %c, %negative ]
  %high = phi double [ %b, %positive ], [ %d, %negative ]
  %sum = call double @observe(double %low, double %high, double %a0,
    double %a1, double %a2, double %a3, double %a4, double %a5, double %a6,
    double %a7, double %a8, double %a9)
  %next = add i32 %i, 1
  %more = icmp slt i32 %next, %n
  br i1 %more, label %header, label %exit, !prof !0

exit:
  ret double %sum
}

!0 = !{!"branch_weights", i32 100, i32 1}
