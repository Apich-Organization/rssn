target triple = "x86_64-pc-windows-msvc"

declare double @llvm.exp.f64(double)
declare void @consume(<4 x i64>, i64)

define double @dispatch(i64 noundef %choice, i64 noundef %v,
                        double noundef %x) {
entry:
  %lane = insertelement <4 x i64> poison, i64 %v, i64 0
  %wide = shufflevector <4 x i64> %lane, <4 x i64> poison, <4 x i32> zeroinitializer
  br label %header

header:
  %i = phi i64 [ 0, %entry ], [ %next, %latch ]
  %pick = icmp eq i64 %i, %choice
  br i1 %pick, label %arm, label %call, !prof !0

arm:
  call void @consume(<4 x i64> %wide, i64 %i)
  br label %latch

call:
  %r = call double @llvm.exp.f64(double %x)
  br label %latch

latch:
  %out = phi double [ 0.0, %arm ], [ %r, %call ]
  %next = add i64 %i, 1
  %done = icmp eq i64 %next, 32
  br i1 %done, label %exit, label %header, !prof !1

exit:
  ret double %out
}

!0 = !{!"branch_weights", i32 1, i32 1000}
!1 = !{!"branch_weights", i32 1, i32 31}
