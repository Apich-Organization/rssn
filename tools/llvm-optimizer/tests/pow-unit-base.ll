target triple = "x86_64-unknown-linux-gnu"

declare double @llvm.pow.f64(double, double)
declare void @observe()

define double @unit_bases(double %node, double %x, double %y, double %e0,
                          double %e1, double %e2, double %initial, i32 %rounds) {
entry:
  br label %header

header:
  %t = phi double [ %node, %entry ], [ %next.t, %join ]
  %i = phi i32 [ 0, %entry ], [ %next, %join ]
  %odd = trunc i32 %i to i1
  br i1 %odd, label %right, label %left

left:
  %b0 = fsub double 1.0, %t
  %p0 = call double @llvm.pow.f64(double %b0, double %e0)
  %r0 = fmul double %initial, %p0
  %sx = fmul double %t, %x
  %bx = fsub double 1.0, %sx
  %px = call double @llvm.pow.f64(double %bx, double %e1)
  %rx = fmul double %r0, %px
  %sy = fmul double %y, %t
  %by = fsub double 1.0, %sy
  %py = call double @llvm.pow.f64(double %by, double %e2)
  %ry = fmul double %py, %rx
  br label %join

right:
  %m = fmul double %t, 2.0
  %c0 = fsub double 1.0, %m
  %q0 = call double @llvm.pow.f64(double %c0, double %e0)
  %s0 = fmul double %initial, %q0
  %tx = fmul double %m, %x
  %cx = fsub double 1.0, %tx
  %qx = call double @llvm.pow.f64(double %cx, double %e1)
  %ux = fmul double %s0, %qx
  %ty = fmul double %m, %y
  %cy = fsub double 1.0, %ty
  %qy = call double @llvm.pow.f64(double %cy, double %e2)
  %uy = fmul double %ux, %qy
  br label %join

join:
  %ret = phi double [ %ry, %left ], [ %uy, %right ]
  %next = add i32 %i, 1
  %next.t = fmul double %t, 0.5
  %more = icmp ult i32 %next, %rounds
  br i1 %more, label %header, label %exit, !prof !0

exit:
  ret double %ret
}

!0 = !{!"branch_weights", i32 100, i32 1}
