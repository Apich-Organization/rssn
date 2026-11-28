target triple = "x86_64-unknown-linux-gnu"

%Result = type { i64, i64, i64 }
@exposed = global ptr @address_taken
declare i64 @opaque(i64) nounwind

define internal fastcc void @leaf_result(ptr noalias noundef sret(%Result) align 8 %out, ptr nocapture noundef readonly %a, i64 noundef range(i64 0, 1152921504606846976) %na, ptr nocapture noundef readonly %b, i64 noundef %nb, ptr nocapture noundef readonly %c, i64 noundef %nc, i1 noundef zeroext %inverse) nounwind noinline {
  %x = load i64, ptr %a, align 8
  %y = load i64, ptr %b, align 8
  %z = load i64, ptr %c, align 8
  %one = zext i1 %inverse to i64
  %lengths = add i64 %na, %nb
  %total = add i64 %lengths, %nc
  %m = xor i64 %x, %y
  %n = xor i64 %m, %z
  %r = add i64 %n, %total
  %v = add i64 %r, %one
  store i64 %v, ptr %out, align 8
  ret void
}

define i64 @drive(ptr %a, ptr %b, ptr %c, i64 %rounds) nounwind {
entry:
  %out = alloca %Result, align 8
  br label %loop
loop:
  %i = phi i64 [ 0, %entry ], [ %next, %loop ]
  %acc = phi i64 [ 0, %entry ], [ %sum, %loop ]
  %even = trunc i64 %i to i1
  call fastcc void @leaf_result(ptr noundef sret(%Result) align 8 %out, ptr readonly %a, i64 noundef 9, ptr readonly %b, i64 13, ptr readonly %c, i64 17, i1 zeroext %even)
  %v = load i64, ptr %out, align 8
  %sum = add i64 %acc, %v
  %next = add i64 %i, 1
  %done = icmp uge i64 %next, %rounds
  br i1 %done, label %exit, label %loop
exit:
  ret i64 %sum
}

define internal i64 @address_taken(i64 %a, i64 %b, i64 %c, i64 %d, i64 %e, i64 %f, i64 %stack) nounwind noinline {
  %r = add i64 %a, %stack
  ret i64 %r
}

define internal i64 @nonleaf(i64 %a, i64 %b, i64 %c, i64 %d, i64 %e, i64 %f, i64 %stack) nounwind noinline {
  %v = call i64 @opaque(i64 %a)
  %r = add i64 %v, %stack
  ret i64 %r
}

define i64 @drive_nonleaf(i64 %v) nounwind {
  %a = add i64 %v, 1
  %b = xor i64 %a, 2
  %c = call i64 @nonleaf(i64 %a, i64 %a, i64 %a, i64 %a, i64 %b, i64 %v, i64 %b)
  %r = add i64 %v, %c
  ret i64 %r
}

define internal i64 @leaf_forward(i64 %a, i64 %b, i64 %c, i64 %d, i64 %e, i64 %f, i64 %stack) nounwind noinline {
  %ab = add i64 %a, %b
  %cd = xor i64 %c, %d
  %ef = add i64 %e, %f
  %p = add i64 %stack, %ab
  %q = xor i64 %cd, %ef
  %r = add i64 %p, %q
  ret i64 %r
}

define i64 @forward(i64 %a, i64 %b, i64 %c, i64 %d, i64 %e, i64 %f, i64 %stack) nounwind {
  %v = call i64 @leaf_forward(i64 %a, i64 %b, i64 %c, i64 %d, i64 %e, i64 %f, i64 %stack)
  %x = add i64 %v, %a
  %y = xor i64 %x, %b
  ret i64 %y
}
