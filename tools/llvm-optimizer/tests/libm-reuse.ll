target triple = "x86_64-unknown-linux-gnu"

declare double @sinh(double)
declare double @llvm.fabs.f64(double)
declare double @llvm.exp.f64(double)

define double @mirrored_pair(double %x) {
entry:
  br label %loop

loop:
  %k = phi double [ %x, %entry ], [ %next, %right ]
  %a = call double @sinh(double %k)
  %ax = fmul double %a, 2.000000e+00
  %aa = call double @llvm.fabs.f64(double %ax)
  %e = call double @llvm.exp.f64(double %aa)
  %p = fadd double %e, 1.000000e+00
  %d = fdiv double 2.000000e+00, %p
  br label %right

right:
  %neg = fneg double %k
  %b = call double @sinh(double %neg)
  %bx = fmul double %b, 2.000000e+00
  %ba = call double @llvm.fabs.f64(double %bx)
  %f = call double @llvm.exp.f64(double %ba)
  %q = fadd double %f, 1.000000e+00
  %r = fdiv double 2.000000e+00, %q
  %sum = fadd double %d, %r
  %next = fadd double %k, 1.000000e+00
  %again = fcmp olt double %next, %x
  br i1 %again, label %loop, label %exit

exit:
  ret double %sum
}
