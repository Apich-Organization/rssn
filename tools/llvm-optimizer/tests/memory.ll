target triple = "x86_64-unknown-linux-gnu"

declare void @llvm.lifetime.end.p0(ptr)
declare void @quiet(double) nounwind nofree nosync willreturn memory(none)
declare void @unknown()

define double @merge_disjoint(ptr noundef noalias %p, ptr noundef noalias %q) {
  %a = load double, ptr %p
  store double %a, ptr %q
  %b = load double, ptr %p
  %s = fadd double %a, %b
  ret double %s
}

define double @keep_alias(ptr %p, ptr %q) {
  %a = load double, ptr %p
  store double %a, ptr %q
  %b = load double, ptr %p
  %s = fadd double %a, %b
  ret double %s
}

define double @keep_atomic(ptr %p, ptr noalias %q) {
  %a = load double, ptr %p
  store double %a, ptr %q
  fence acquire
  %b = load double, ptr %p
  %s = fadd double %a, %b
  ret double %s
}

define double @keep_volatile(ptr %p, ptr noalias %q) {
  %a = load double, ptr %p
  store volatile double %a, ptr %q
  %b = load double, ptr %p
  %s = fadd double %a, %b
  ret double %s
}

define double @keep_lifetime(ptr noalias %q) {
  %p = alloca double
  store double 1.0, ptr %p
  %a = load double, ptr %p
  store double %a, ptr %q
  call void @llvm.lifetime.end.p0(ptr %p)
  %b = load double, ptr %p
  ret double %b
}

define double @keep_call_pressure(ptr %p) {
  %a = load double, ptr %p
  call void @quiet(double %a)
  %b = load double, ptr %p
  %s = fadd double %b, 1.0
  ret double %s
}

define double @keep_unknown(ptr %p, ptr noalias %q) {
  %a = load double, ptr %p
  store double %a, ptr %q
  call void @unknown()
  %b = load double, ptr %p
  %s = fadd double %a, %b
  ret double %s
}

define ptr @keep_metadata(ptr %p, ptr noalias %q) {
  %a = load ptr, ptr %p, !nonnull !0
  store ptr %a, ptr %q
  %b = load ptr, ptr %p
  ret ptr %b
}

!0 = !{}
