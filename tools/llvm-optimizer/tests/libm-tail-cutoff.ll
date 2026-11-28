target triple = "x86_64-unknown-linux-gnu"

declare double @llvm.log.f64(double)
declare double @llvm.exp.f64(double)
declare i1 @llvm.is.fpclass.f64(double, i32 immarg)
declare void @observe()

define double @tail_existing(double %seed, double %c, double %scale, i32 %rounds) {
entry:
  br label %header

header:
  %i = phi i32 [ 0, %entry ], [ %next, %join ]
  %t = phi double [ %seed, %entry ], [ %next.t, %join ]
  %odd = trunc i32 %i to i1
  %x = call noundef double @llvm.exp.f64(double %t)
  %weight = fmul double %x, 2.0
  %bad = call i1 @llvm.is.fpclass.f64(double %weight, i32 615)
  br i1 %bad, label %join, label %valid

valid:
  br i1 %odd, label %right, label %left

left:
  %pos = fcmp ugt double %x, 0.0
  br i1 %pos, label %left.score, label %join

left.score:
  %log = call double @llvm.log.f64(double %x)
  %product = fmul double %c, %log
  %decay = fmul double %scale, %x
  %score = fsub double %product, %decay
  %cut = fcmp olt double %score, -745.0
  br i1 %cut, label %join, label %left.exp

left.exp:
  %out = call double @llvm.exp.f64(double %score)
  br label %join

right:
  %pos.r = fcmp ugt double %x, 0.0
  br i1 %pos.r, label %right.score, label %join

right.score:
  %log.r = call double @llvm.log.f64(double %x)
  %product.r = fmul double %log.r, %c
  %decay.r = fmul double %x, %scale
  %score.r = fsub double %product.r, %decay.r
  %cut.r = fcmp ogt double -745.0, %score.r
  br i1 %cut.r, label %join, label %right.exp

right.exp:
  %out.r = call double @llvm.exp.f64(double %score.r)
  br label %join

join:
  %answer = phi double [ 0.0, %header ], [ 0.0, %left ], [ 0.0, %left.score ],
                       [ %out, %left.exp ], [ 0.0, %right ], [ 0.0, %right.score ],
                       [ %out.r, %right.exp ]
  %next.t = fadd double %t, 0.25
  %next = add i32 %i, 1
  %more = icmp ult i32 %next, %rounds
  br i1 %more, label %header, label %exit, !prof !0

exit:
  ret double %answer
}

define double @tail_zero_numerator(double %seed, double %c, double %den, i32 %rounds) {
entry:
  br label %header

header:
  %i = phi i32 [ 0, %entry ], [ %next, %join ]
  %t = phi double [ %seed, %entry ], [ %next.t, %join ]
  %odd = trunc i32 %i to i1
  %raw = call noundef double @llvm.exp.f64(double %t)
  %x = fadd double 1.0, %raw
  %weight = fmul double %raw, 2.0
  %bad = call i1 @llvm.is.fpclass.f64(double %weight, i32 615)
  br i1 %bad, label %join, label %valid

valid:
  br i1 %odd, label %right, label %left

left:
  %log = call double @llvm.log.f64(double %x)
  %product = fmul double %c, %log
  %score = fsub double %product, %x
  %out = call double @llvm.exp.f64(double %score)
  %ratio = fdiv double %out, %den
  br label %join

right:
  %log.r = call double @llvm.log.f64(double %x)
  %product.r = fmul double %log.r, %c
  %score.r = fsub double %product.r, %x
  %out.r = call double @llvm.exp.f64(double %score.r)
  %ratio.r = fdiv double %out.r, %den
  br label %join

join:
  %answer = phi double [ 0.0, %header ], [ %ratio, %left ], [ %ratio.r, %right ]
  %next.t = fadd double %t, 0.25
  %next = add i32 %i, 1
  %more = icmp ult i32 %next, %rounds
  br i1 %more, label %header, label %exit, !prof !0

exit:
  ret double %answer
}

!0 = !{!"branch_weights", i32 100, i32 1}

; A score with a loaded log argument does not establish a bounded tail.
define double @tail_load_input(ptr %data, double %coefficient, i32 %count) {
entry:
  br label %loop

loop:
  %i = phi i32 [ 0, %entry ], [ %next, %loop ]
  %cell = getelementptr double, ptr %data, i32 %i
  %node = load double, ptr %cell
  %log = call double @llvm.log.f64(double %node)
  %product = fmul double %coefficient, %log
  %score = fsub double %product, %node
  %output = call double @llvm.exp.f64(double %score)
  %next = add i32 %i, 1
  %more = icmp ult i32 %next, %count
  br i1 %more, label %loop, label %exit

exit:
  ret double %output
}
