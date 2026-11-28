target triple = "x86_64-unknown-linux-gnu"

declare <2 x double> @llvm.exp.v2f64(<2 x double>)
declare <4 x double> @llvm.log.v4f64(<4 x double>)

define <2 x double> @uniform(i1 %c, <2 x double> %x, <2 x double> %other) {
  %one = insertelement <2 x i1> poison, i1 %c, i64 0
  %mask = shufflevector <2 x i1> %one, <2 x i1> poison, <2 x i32> zeroinitializer
  %y = call <2 x double> @llvm.exp.v2f64(<2 x double> %x)
  %r = select <2 x i1> %mask, <2 x double> %y, <2 x double> %other
  ret <2 x double> %r
}

define <4 x double> @false_edge(i1 %c, <4 x double> %x, <4 x double> %other) {
  %y = call <4 x double> @llvm.log.v4f64(<4 x double> %x)
  %r = select i1 %c, <4 x double> %other, <4 x double> %y, !prof !0
  ret <4 x double> %r
}

define <2 x double> @keep_nonuniform(<2 x i1> %c, <2 x double> %x, <2 x double> %other) {
  %y = call <2 x double> @llvm.exp.v2f64(<2 x double> %x)
  %r = select <2 x i1> %c, <2 x double> %y, <2 x double> %other
  ret <2 x double> %r
}

define <2 x double> @keep_shared(i1 %c, <2 x double> %x, <2 x double> %other) {
  %y = call <2 x double> @llvm.exp.v2f64(<2 x double> %x)
  %r = select i1 %c, <2 x double> %y, <2 x double> %other
  %sum = fadd <2 x double> %y, %r
  ret <2 x double> %sum
}

define <2 x double> @keep_unprofitable(i1 %c, <2 x double> %x, <2 x double> %other) {
  %y = call <2 x double> @llvm.exp.v2f64(<2 x double> %x)
  %r = select i1 %c, <2 x double> %y, <2 x double> %other, !prof !1
  ret <2 x double> %r
}

define <2 x double> @keep_strict(i1 %c, <2 x double> %x, <2 x double> %other) strictfp {
  %y = call <2 x double> @llvm.exp.v2f64(<2 x double> %x)
  %r = select i1 %c, <2 x double> %y, <2 x double> %other
  ret <2 x double> %r
}

define <2 x double> @shared_branch(i1 %c, <2 x double> %x, <2 x double> %z, <2 x double> %other) {
  %a = call <2 x double> @llvm.exp.v2f64(<2 x double> %x)
  %s = select i1 %c, <2 x double> %a, <2 x double> %other
  %b = call <2 x double> @llvm.exp.v2f64(<2 x double> %z)
  %t = select i1 %c, <2 x double> %b, <2 x double> %other
  %r = fadd <2 x double> %s, %t
  ret <2 x double> %r
}

define <2 x double> @dependent_branches(i1 %c, <2 x double> %x, <2 x double> %other) {
  %a = call <2 x double> @llvm.exp.v2f64(<2 x double> %x)
  %s = select i1 %c, <2 x double> %a, <2 x double> %other
  %b = call <2 x double> @llvm.exp.v2f64(<2 x double> %s)
  %t = select i1 %c, <2 x double> %b, <2 x double> %other
  ret <2 x double> %t
}

!0 = !{!"branch_weights", i32 20, i32 80}
!1 = !{!"branch_weights", i32 99, i32 1}
