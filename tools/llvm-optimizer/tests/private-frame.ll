target triple = "x86_64-unknown-linux-gnu"

@escape = global ptr null

declare void @fill_local(ptr noundef captures(none) writeonly)
declare void @fill_unknown(ptr noundef writeonly)
declare double @observer()
declare double @touch(ptr noundef captures(none))
declare void @llvm.lifetime.end.p0(ptr captures(none))

define double @private_frame() {
  %frame = alloca double, align 8
  call void @fill_local(ptr %frame)
  %old = load double, ptr %frame, align 8
  %seen = call double @observer()
  %result = fmul double %old, %seen
  ret double %result
}

define double @unknown_capture() {
  %frame = alloca double, align 8
  call void @fill_unknown(ptr %frame)
  %old = load double, ptr %frame, align 8
  %seen = call double @observer()
  %result = fmul double %old, %seen
  ret double %result
}

define double @escaped_frame() {
  %frame = alloca double, align 8
  call void @fill_local(ptr %frame)
  store ptr %frame, ptr @escape
  %old = load double, ptr %frame, align 8
  %seen = call double @observer()
  %result = fmul double %old, %seen
  ret double %result
}

define double @frame_visible_to_call() {
  %frame = alloca double, align 8
  call void @fill_local(ptr %frame)
  %old = load double, ptr %frame, align 8
  %seen = call double @touch(ptr %frame)
  %result = fmul double %old, %seen
  ret double %result
}

define double @expired_frame() {
  %frame = alloca double, align 8
  call void @fill_local(ptr %frame)
  %old = load double, ptr %frame, align 8
  call void @llvm.lifetime.end.p0(ptr %frame)
  %seen = call double @observer()
  %result = fmul double %old, %seen
  ret double %result
}

define double @private_frame_gep() {
  %frame = alloca [2 x double], align 8
  %field = getelementptr inbounds [2 x double], ptr %frame, i64 0, i64 1
  call void @fill_local(ptr %field)
  %old = load double, ptr %field, align 8
  %seen = call double @observer()
  %result = fmul double %old, %seen
  ret double %result
}

declare void @llvm.localescape(...)

define double @frame_exposed_implicitly() {
  %frame = alloca double, align 8
  call void (...) @llvm.localescape(ptr %frame)
  call void @fill_local(ptr %frame)
  %old = load double, ptr %frame, align 8
  %seen = call double @observer()
  %result = fmul double %old, %seen
  ret double %result
}
