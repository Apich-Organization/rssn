target triple = "x86_64-unknown-linux-gnu"

declare double @quiet_gain(i64) nounwind nofree nosync willreturn memory(none)

; Exactly five GP values already cross the call. Moving either read alone
; adds a sixth, while moving %key saves one. The two reads share the extension
; of %p. Only the joint plan relieves pressure without exceeding the GP budget.
define double @shared_joint(ptr %p, i64 %a, i64 %b, i64 %c, i64 %d) {
  %q = getelementptr inbounds double, ptr %p, i64 1
  %x = load double, ptr %p, align 8
  %y = load double, ptr %q, align 8
  %key = xor i64 %a, 17
  %gain = call double @quiet_gain(i64 0)
  %xa = fmul double %x, %gain
  %ya = fadd double %y, %gain
  %both = fadd double %xa, %ya
  %ab = add i64 %a, %b
  %cd = add i64 %c, %d
  %abcd = add i64 %ab, %cd
  %mixed = add i64 %abcd, %key
  %fp = sitofp i64 %mixed to double
  %total = fadd double %both, %fp
  ret double %total
}

; Six GP values survive the call. The seventh is only a call operand.
define double @gp_call_excess(i64 %a, i64 %b, i64 %c, i64 %d, i64 %e, i64 %f, i64 %key) {
  %gain = call double @quiet_gain(i64 %key)
  %ab = add i64 %a, %b
  %cd = add i64 %c, %d
  %ef = add i64 %e, %f
  %abcd = add i64 %ab, %cd
  %sum = add i64 %abcd, %ef
  %fp = sitofp i64 %sum to double
  %total = fadd double %gain, %fp
  ret double %total
}

; TTI must count a 256-bit vector as one register on AVX, two on SSE2.
define <4 x double> @width_probe(<4 x double> %data, i64 %key) {
  %c = call double @quiet_gain(i64 %key)
  %first = insertelement <4 x double> poison, double %c, i32 0
  %wide = shufflevector <4 x double> %first, <4 x double> poison, <4 x i32> zeroinitializer
  %r = fadd <4 x double> %data, %wide
  ret <4 x double> %r
}
