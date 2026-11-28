target triple = "x86_64-pc-windows-msvc"
target datalayout = "e-m:w-p270:32:32-p271:32:32-p272:64:64-i64:64-i128:128-f80:128-n8:16:32:64-S128"

declare win64cc i64 @opaque(double)
declare x86_64_sysvcc void @foreign()
declare double @llvm.exp.f64(double)

define double @win_pressure(ptr %p, <4 x double> %wide, <2 x double> %narrow,
                            double %explicit_arg, double %math_arg) {
entry:
  %saved = load i64, ptr %p, align 8
  %r = call win64cc i64 @opaque(double %explicit_arg)
  call x86_64_sysvcc void @foreign()
  %math = call double @llvm.exp.f64(double %math_arg)
  %bits = add i64 %saved, %r
  %integer = uitofp i64 %bits to double
  %x = extractelement <4 x double> %wide, i32 3
  %y = extractelement <2 x double> %narrow, i32 1
  %xy = fadd double %x, %y
  %a = fadd double %xy, %integer
  %ret = fadd double %a, %math
  ret double %ret
}
