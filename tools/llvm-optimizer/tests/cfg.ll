target triple = "x86_64-unknown-linux-gnu"

declare i64 @gate(i64)

define double @split_defined(i1 %branch, <2 x double> noundef %source, i64 %key) {
entry:
  %neg = fneg <2 x double> %source
  %call = call i64 @gate(i64 %key)
  %odd = trunc i64 %call to i1
  br i1 %branch, label %left, label %right

left:
  %sum = fadd <2 x double> %neg, %source
  %leftval = extractelement <2 x double> %sum, i1 %odd
  ret double %leftval

right:
  %diff = fsub <2 x double> %neg, %source
  %rightval = extractelement <2 x double> %diff, i1 %odd
  ret double %rightval
}

define double @keep_undef_operand(i1 %branch, <2 x double> %source, i64 %key) {
entry:
  %neg = fneg <2 x double> %source
  %call = call i64 @gate(i64 %key)
  %odd = trunc i64 %call to i1
  br i1 %branch, label %left, label %right

left:
  %sum = fadd <2 x double> %neg, %source
  %leftval = extractelement <2 x double> %sum, i1 %odd
  ret double %leftval

right:
  %diff = fsub <2 x double> %neg, %source
  %rightval = extractelement <2 x double> %diff, i1 %odd
  ret double %rightval
}

define double @keep_phi_edge(i1 %branch, <2 x double> noundef %source, i64 %key) {
entry:
  %neg = fneg <2 x double> %source
  %call = call i64 @gate(i64 %key)
  %odd = trunc i64 %call to i1
  br i1 %branch, label %left, label %join

left:
  br label %join

join:
  %either = phi <2 x double> [ %neg, %entry ], [ %source, %left ]
  %val = extractelement <2 x double> %either, i1 %odd
  ret double %val
}

define double @sink_single_successor(i1 %branch, <2 x double> %source, i64 %key) {
entry:
  %neg = fneg <2 x double> %source
  %call = call i64 @gate(i64 %key)
  %odd = trunc i64 %call to i1
  br i1 %branch, label %left, label %right

left:
  %sum = fadd <2 x double> %neg, %source
  %leftval = extractelement <2 x double> %sum, i1 %odd
  ret double %leftval

right:
  ret double 0.0
}
