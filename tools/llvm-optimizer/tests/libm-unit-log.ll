target triple = "x86_64-unknown-linux-gnu"

declare double @llvm.log.f64(double)
declare double @llvm.exp.f64(double)
declare void @observe(double)

define double @unit_log_chain(double %seed, double %c, i32 %rounds) {
entry:
  br label %header

header:
  %i = phi i32 [ 0, %entry ], [ %next, %join ]
  %t = phi double [ %seed, %entry ], [ %next.t, %join ]
  %odd = trunc i32 %i to i1
  br i1 %odd, label %right, label %left

left:
  %a = call noundef double @llvm.exp.f64(double %t)
  %x = fadd double %a, 1.0
  %lx = call double @llvm.log.f64(double %x)
  %cx = fmul double %c, %lx
  %dx = fsub double %cx, %x
  %ex = call double @llvm.exp.f64(double %dx)
  br label %join

right:
  %b = call noundef double @llvm.exp.f64(double %t)
  %y = fadd double 1.0, %b
  %ly = call double @llvm.log.f64(double %y)
  %cy = fmul double %ly, %c
  %dy = fsub double %cy, %y
  %ey = call double @llvm.exp.f64(double %dy)
  br label %join

join:
  %result = phi double [ %ex, %left ], [ %ey, %right ]
  %next = add i32 %i, 1
  %next.t = fadd double %t, 0.25
  %more = icmp ult i32 %next, %rounds
  br i1 %more, label %header, label %exit, !prof !0

exit:
  ret double %result
}

!0 = !{!"branch_weights", i32 100, i32 1}

; The matcher sees these ordinary operations before recognizing a full chain.
define double @unit_log_load_input(ptr %data, i32 %count) {
entry:
  br label %loop

loop:
  %i = phi i32 [ 0, %entry ], [ %next, %loop ]
  %cell = getelementptr double, ptr %data, i32 %i
  %value = load double, ptr %cell
  %log = call double @llvm.log.f64(double %value)
  %next = add i32 %i, 1
  %more = icmp ult i32 %next, %count
  br i1 %more, label %loop, label %exit

exit:
  ret double %log
}

define i64 @unit_log_bitcast_after_log(i32 %count) {
entry:
  br label %loop

loop:
  %i = phi i32 [ 0, %entry ], [ %next, %loop ]
  %seed = phi double [ -900.0, %entry ], [ %step, %loop ]
  %exp = call double @llvm.exp.f64(double %seed)
  %sum = fadd double %exp, 1.0
  %log = call double @llvm.log.f64(double %sum)
  %bits = bitcast double %log to i64
  %step = fadd double %seed, 1.0
  %next = add i32 %i, 1
  %more = icmp ult i32 %next, %count
  br i1 %more, label %loop, label %exit

exit:
  ret i64 %bits
}

define i64 @unit_log_bitcast_after_multiply(i32 %count, double %coefficient) {
entry:
  br label %loop

loop:
  %i = phi i32 [ 0, %entry ], [ %next, %loop ]
  %seed = phi double [ -900.0, %entry ], [ %step, %loop ]
  %exp = call double @llvm.exp.f64(double %seed)
  %sum = fadd double %exp, 1.0
  %log = call double @llvm.log.f64(double %sum)
  %product = fmul double %coefficient, %log
  %bits = bitcast double %product to i64
  %step = fadd double %seed, 1.0
  %next = add i32 %i, 1
  %more = icmp ult i32 %next, %count
  br i1 %more, label %loop, label %exit

exit:
  ret i64 %bits
}
