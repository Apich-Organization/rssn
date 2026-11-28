target triple = "x86_64-unknown-linux-gnu"

define <2 x double> @pack_contiguous(ptr %p, <2 x double> %v) {
  %q = getelementptr inbounds double, ptr %p, i64 1
  %a = load double, ptr %p, align 8, !noundef !0
  %b = load double, ptr %q, align 8, !noundef !0
  %lo = insertelement <2 x double> poison, double %a, i32 0
  %wide = insertelement <2 x double> %lo, double %b, i32 1
  %r = fadd <2 x double> %v, %wide
  ret <2 x double> %r
}

define <4 x double> @concat_halves(ptr %p, <4 x double> %v) {
  %q = getelementptr inbounds double, ptr %p, i64 2
  %a = load <2 x double>, ptr %p, align 8
  %b = load <2 x double>, ptr %q, align 8
  %wide = shufflevector <2 x double> %a, <2 x double> %b, <4 x i32> <i32 0, i32 1, i32 2, i32 3>
  %r = fadd <4 x double> %v, %wide
  ret <4 x double> %r
}

define <2 x double> @keep_gap(ptr %p, <2 x double> %v) {
  %q = getelementptr inbounds double, ptr %p, i64 2
  %a = load double, ptr %p, align 8
  %b = load double, ptr %q, align 8
  %lo = insertelement <2 x double> poison, double %a, i32 0
  %wide = insertelement <2 x double> %lo, double %b, i32 1
  %r = fadd <2 x double> %v, %wide
  ret <2 x double> %r
}

define <2 x double> @keep_swapped(ptr %p, <2 x double> %v) {
  %q = getelementptr inbounds double, ptr %p, i64 1
  %a = load double, ptr %p, align 8
  %b = load double, ptr %q, align 8
  %lo = insertelement <2 x double> poison, double %b, i32 0
  %wide = insertelement <2 x double> %lo, double %a, i32 1
  %r = fadd <2 x double> %v, %wide
  ret <2 x double> %r
}

define <2 x double> @keep_clobbered(ptr %p, ptr %possibly_aliasing, <2 x double> %v) {
  %q = getelementptr inbounds double, ptr %p, i64 1
  %a = load double, ptr %p, align 8
  store double 1.0, ptr %possibly_aliasing
  %b = load double, ptr %q, align 8
  %lo = insertelement <2 x double> poison, double %a, i32 0
  %wide = insertelement <2 x double> %lo, double %b, i32 1
  %r = fadd <2 x double> %v, %wide
  ret <2 x double> %r
}

define <2 x double> @keep_volatile_lane(ptr %p, <2 x double> %v) {
  %q = getelementptr inbounds double, ptr %p, i64 1
  %a = load volatile double, ptr %p, align 8
  %b = load double, ptr %q, align 8
  %lo = insertelement <2 x double> poison, double %a, i32 0
  %wide = insertelement <2 x double> %lo, double %b, i32 1
  %r = fadd <2 x double> %v, %wide
  ret <2 x double> %r
}

define double @merge_equivalent_addresses(ptr noalias %p, ptr noalias %out) {
  %physical1 = getelementptr inbounds double, ptr %p, i64 2
  %a = load double, ptr %physical1, align 8
  store double %a, ptr %out, align 8
  %next = getelementptr inbounds double, ptr %p, i64 1
  %physical2 = getelementptr inbounds double, ptr %next, i64 1
  %b = load double, ptr %physical2, align 8
  %r = fadd double %a, %b
  ret double %r
}

!0 = !{}

declare double @opaque_scale(ptr)

define double @pack_live_scalars(ptr %p) {
  %q = getelementptr inbounds double, ptr %p, i64 1
  %a = load double, ptr %p, align 8
  %b = load double, ptr %q, align 8
  %coef = call double @opaque_scale(ptr %p)
  %x = fmul double %a, %coef
  %y = fmul double %b, %coef
  %r = fdiv double %x, %y
  ret double %r
}

define double @keep_lanes_clobbered(ptr %p, ptr %possibly_aliasing) {
  %q = getelementptr inbounds double, ptr %p, i64 1
  %a = load double, ptr %p, align 8
  store double 1.0, ptr %possibly_aliasing
  %b = load double, ptr %q, align 8
  %coef = call double @opaque_scale(ptr %p)
  %x = fmul double %a, %coef
  %y = fmul double %b, %coef
  %r = fdiv double %x, %y
  ret double %r
}


; Independent lane projections let SelectionDAG scalarize the combined read.
; Mixed scalar operations stay scalar until there is a packed lowering.
define double @keep_mixed_scalar_ops(ptr %p) {
  %q = getelementptr inbounds double, ptr %p, i64 1
  %a = load double, ptr %p, align 8
  %b = load double, ptr %q, align 8
  %coef = call double @opaque_scale(ptr %p)
  %x = fmul double %a, %coef
  %y = fadd double %b, %coef
  %r = fdiv double %x, %y
  ret double %r
}
