target triple = "x86_64-unknown-linux-gnu"
target datalayout = "e-m:e-p270:32:32-p271:32:32-p272:64:64-i64:64-f80:128-n8:16:32:64-S128"

declare void @consume_fields(i16, i8, i8, i32, i16, i8, i8, i32, i16, i8, i8, i32)

define void @field_loop(ptr %a, ptr %b, ptr %c, i64 %size) {
entry:
  br label %loop

loop:
  %i = phi i64 [ 0, %entry ], [ %next, %loop ]
  %ra = getelementptr inbounds [8 x i8], ptr %a, i64 %i
  %rb = getelementptr inbounds [8 x i8], ptr %b, i64 %i
  %rc = getelementptr inbounds [8 x i8], ptr %c, i64 %i
  %ra.high = getelementptr inbounds i8, ptr %ra, i64 4
  %rb.high = getelementptr inbounds i8, ptr %rb, i64 4
  %rc.high = getelementptr inbounds i8, ptr %rc, i64 4
  %a3 = load i32, ptr %ra.high, align 4
  %b3 = load i32, ptr %rb.high, align 4
  %c3 = load i32, ptr %rc.high, align 4
  %a0 = load i16, ptr %ra, align 4
  %b0 = load i16, ptr %rb, align 4
  %c0 = load i16, ptr %rc, align 4
  %ra.middle = getelementptr inbounds i8, ptr %ra, i64 2
  %rb.middle = getelementptr inbounds i8, ptr %rb, i64 2
  %rc.middle = getelementptr inbounds i8, ptr %rc, i64 2
  %a1 = load i8, ptr %ra.middle, align 2
  %b1 = load i8, ptr %rb.middle, align 2
  %c1 = load i8, ptr %rc.middle, align 2
  %ra.byte = getelementptr inbounds i8, ptr %ra, i64 3
  %rb.byte = getelementptr inbounds i8, ptr %rb, i64 3
  %rc.byte = getelementptr inbounds i8, ptr %rc, i64 3
  %a2 = load i8, ptr %ra.byte, align 1
  %b2 = load i8, ptr %rb.byte, align 1
  %c2 = load i8, ptr %rc.byte, align 1
  call void @consume_fields(i16 %a0, i8 %a1, i8 %a2, i32 %a3,
      i16 %b0, i8 %b1, i8 %b2, i32 %b3, i16 %c0, i8 %c1, i8 %c2, i32 %c3)
  %next = add nuw i64 %i, 1
  %more = icmp ult i64 %next, %size
  br i1 %more, label %loop, label %exit

exit:
  ret void
}
