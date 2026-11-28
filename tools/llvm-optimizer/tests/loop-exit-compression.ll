target triple = "x86_64-pc-windows-msvc"

declare double @seed(double)
declare double @step(double, double, double, double, double, double, double,
                     double, double, double, double)
declare void @touch(double, double, double)

define double @compress_exit(double %x, double %r0, double %r1, double %r2,
                            double %r3, double %r4, double %r5, double %r6,
                            double %r7, double %r8, double %r9, i32 %n) {
entry:
  %a = call double @seed(double %x)
  %b = call double @seed(double %r0)
  %c = call double @seed(double %r1)
  br label %loop

loop:
  %index = phi i32 [ 0, %entry ], [ %next, %loop ]
  %value = phi double [ %x, %entry ], [ %y, %loop ]
  %t = call double @step(double %value, double %r0, double %r1, double %r2,
                        double %r3, double %r4, double %r5, double %r6,
                        double %r7, double %r8, double %r9)
  %y = call double @step(double %t, double %r0, double %r1, double %r2,
                        double %r3, double %r4, double %r5, double %r6,
                        double %r7, double %r8, double %r9)
  %next = add i32 %index, 1
  %more = icmp slt i32 %next, %n
  br i1 %more, label %loop, label %exit, !prof !0

exit:
  %product = fmul double %b, %c
  %coefficient = fdiv double %a, %product
  %result = fmul double %y, %coefficient
  ret double %result
}

!0 = !{!"branch_weights", i32 100, i32 1}
