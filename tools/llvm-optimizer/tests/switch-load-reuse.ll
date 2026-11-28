target triple = "x86_64-unknown-linux-gnu"

define i64 @hash_mode(i32 %mode, ptr %source, ptr %out) {
entry:
  switch i32 %mode, label %low [
    i32 1, label %one
    i32 2, label %two
  ]

low:
  %part = load i32, ptr %source, align 4
  %h0 = mul i32 %part, 2654435761
  %h0.ext = zext i32 %h0 to i64
  %w0 = load i64, ptr %source, align 4
  br label %join

one:
  %w1 = load i64, ptr %source, align 8
  %h1 = lshr i64 %w1, 13
  br label %join

two:
  %w2 = load i64, ptr %source, align 1
  %h2 = lshr i64 %w2, 19
  br label %join

join:
  %word = phi i64 [ %w0, %low ], [ %w1, %one ], [ %w2, %two ]
  %hash = phi i64 [ %h0.ext, %low ], [ %h1, %one ], [ %h2, %two ]
  store i64 %hash, ptr %out
  ret i64 %word
}
