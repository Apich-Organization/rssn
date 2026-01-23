//! The macros used by the symbolic_core_ffi module.
//! These macros are not intended to be used by users.

/// Macro for generating handle-based FFI unary functions.

macro_rules! handle_unary_ffi {
    ($ffi_name:ident, $rs_func:ident) => {
        #[unsafe(no_mangle)]
        /// # Safety
        ///
        /// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
        /// The caller must ensure:
        /// 1. All pointer arguments are valid and point to initialized memory.
        /// 2. The memory layout of passed structures matches the expected C-ABI layout.
        /// 3. Any pointers returned by this function are managed according to the API's ownership rules.
        #[doc = concat!("FFI wrapper for `", stringify!($rs_func), "`.")]
        #[doc = ""]
        #[doc = concat!("This function takes a handle, parses it into `Expr`, calls `Expr::", stringify!($rs_func), "`, and returns the result as a handle.")]
        pub unsafe extern "C" fn $ffi_name(
            arg: *const Expr
        ) -> *mut Expr {
            unsafe {
                if arg.is_null() {
                    return std::ptr::null_mut();
                }
                let arg_expr = &*arg;
                Box::into_raw(Box::new(Expr::$rs_func(arg_expr)))
            }
        }
    };
}

/// Macro for generating handle-based FFI binary functions.

macro_rules! handle_binary_ffi {
    ($ffi_name:ident, $rs_func:ident) => {
        #[unsafe(no_mangle)]
        /// # Safety
        ///
        /// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
        /// The caller must ensure:
        /// 1. All pointer arguments are valid and point to initialized memory.
        /// 2. The memory layout of passed structures matches the expected C-ABI layout.
        /// 3. Any pointers returned by this function are managed according to the API's ownership rules.
        #[doc = concat!("FFI wrapper for `", stringify!($rs_func), "`.")]
        #[doc = ""]
        #[doc = concat!("This function takes two handles, parses them into `Expr`, calls `Expr::", stringify!($rs_func), "`, and returns the result as a handle.")]
        pub unsafe extern "C" fn $ffi_name(
            lhs: *const Expr,
            rhs: *const Expr,
        ) -> *mut Expr {
            unsafe {
                if lhs.is_null() || rhs.is_null() {
                    return std::ptr::null_mut();
                }
                let lhs_expr = &*lhs;
                let rhs_expr = &*rhs;
                Box::into_raw(Box::new(Expr::$rs_func(lhs_expr, rhs_expr)))
            }
        }
    };
}

/// Macro for generating handle-based FFI n-ary functions (from C array).

macro_rules! handle_nary_ffi {
    ($ffi_name:ident, $rs_func:ident) => {
        #[unsafe(no_mangle)]
        /// # Safety
        ///
        /// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
        /// The caller must ensure:
        /// 1. All pointer arguments are valid and point to initialized memory.
        /// 2. The memory layout of passed structures matches the expected C-ABI layout.
        /// 3. Any pointers returned by this function are managed according to the API's ownership rules.
        #[doc = concat!("FFI wrapper for `", stringify!($rs_func), "`.")]
        #[doc = ""]
        #[doc = concat!("This function takes a vector of handles, parses them into `Expr`, calls `Expr::", stringify!($rs_func), "`, and returns the result as a handle.")]
        pub unsafe extern "C" fn $ffi_name(
            args: *const *const Expr,
            len: usize
        ) -> *mut Expr {
            unsafe {
                if args.is_null() {
                    return std::ptr::null_mut();
                }
                let mut exprs = Vec::with_capacity(len);
                for i in 0..len {
                    let ptr = *args.add(i);
                    if ptr.is_null() {
                        return std::ptr::null_mut();
                    }
                    exprs.push((*ptr).clone());
                }
                Box::into_raw(Box::new(Expr::$rs_func(exprs)))
            }
        }
    };
}

/// Macro for generating Bincode-based FFI unary functions.

macro_rules! bincode_unary_ffi {
    ($ffi_name:ident, $rs_func:ident) => {
        #[unsafe(no_mangle)]
        /// # Safety
        ///
        /// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
        /// The caller must ensure:
        /// 1. All pointer arguments are valid and point to initialized memory.
        /// 2. The memory layout of passed structures matches the expected C-ABI layout.
        /// 3. Any pointers returned by this function are managed according to the API's ownership rules.
        #[doc = concat!("FFI wrapper for `", stringify!($rs_func), "`.")]
        #[doc = ""]
        #[doc = concat!("This function takes a bincode buffer, parses it into `Expr`, calls `Expr::", stringify!($rs_func), "`, and returns the result as a bincode buffer.")]
        pub extern "C" fn $ffi_name(
            arg_buf: BincodeBuffer
        ) -> BincodeBuffer {
            let arg: Option<Expr> = from_bincode_buffer(&arg_buf);
            if let Some(a) = arg {
                let expr = Expr::$rs_func(&a);
                to_bincode_buffer(&expr)
            } else {
                BincodeBuffer::empty()
            }
        }
    };
}

/// Macro for generating Bincode-based FFI binary functions.

macro_rules! bincode_binary_ffi {
    ($ffi_name:ident, $rs_func:ident) => {
        #[unsafe(no_mangle)]
        /// # Safety
        ///
        /// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
        /// The caller must ensure:
        /// 1. All pointer arguments are valid and point to initialized memory.
        /// 2. The memory layout of passed structures matches the expected C-ABI layout.
        /// 3. Any pointers returned by this function are managed according to the API's ownership rules.
        #[doc = concat!("FFI wrapper for `", stringify!($rs_func), "`.")]
        #[doc = ""]
        #[doc = concat!("This function takes two bincode buffers, parses them into `Expr`, calls `Expr::", stringify!($rs_func), "`, and returns the result as a bincode buffer.")]
        pub extern "C" fn $ffi_name(
            lhs_buf: BincodeBuffer,
            rhs_buf: BincodeBuffer
        ) -> BincodeBuffer {
            let lhs: Option<Expr> = from_bincode_buffer(&lhs_buf);
            let rhs: Option<Expr> = from_bincode_buffer(&rhs_buf);
            match (lhs, rhs) {
                (Some(l), Some(r)) => {
                    let expr = Expr::$rs_func(&l, &r);
                    to_bincode_buffer(&expr)
                },
                _ => BincodeBuffer::empty(),
            }
        }
    };
}

/// Macro for generating JSON-based FFI unary functions.

macro_rules! json_unary_ffi {
    ($ffi_name:ident, $rs_func:ident) => {
        #[unsafe(no_mangle)]
        /// # Safety
        ///
        /// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
        /// The caller must ensure:
        /// 1. All pointer arguments are valid and point to initialized memory.
        /// 2. The memory layout of passed structures matches the expected C-ABI layout.
        /// 3. Any pointers returned by this function are managed according to the API's ownership rules.
        #[doc = concat!("FFI wrapper for `", stringify!($rs_func), "`.")]
        #[doc = ""]
        #[doc = concat!("This function takes a JSON string, parses it into `Expr`, calls `Expr::", stringify!($rs_func), "`, and returns the result as a JSON string.")]
        pub extern "C" fn $ffi_name(
            arg_json: *const c_char
        ) -> *mut c_char {
            let arg: Option<Expr> = from_json_string(arg_json);
            if let Some(a) = arg {
                let expr = Expr::$rs_func(&a);
                to_json_string(&expr)
            } else {
                std::ptr::null_mut()
            }
        }
    };
}

/// Macro for generating JSON-based FFI binary functions.

macro_rules! json_binary_ffi {
    ($ffi_name:ident, $rs_func:ident) => {
        #[unsafe(no_mangle)]
        /// # Safety
        ///
        /// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
        /// The caller must ensure:
        /// 1. All pointer arguments are valid and point to initialized memory.
        /// 2. The memory layout of passed structures matches the expected C-ABI layout.
        /// 3. Any pointers returned by this function are managed according to the API's ownership rules.
        #[doc = concat!("FFI wrapper for `", stringify!($rs_func), "`.")]
        #[doc = ""]
        #[doc = concat!("This function takes two JSON strings, parses them into `Expr`, calls `Expr::", stringify!($rs_func), "`, and returns the result as a JSON string.")]
        pub extern "C" fn $ffi_name(
            lhs_json: *const c_char,
            rhs_json: *const c_char
        ) -> *mut c_char {
            let lhs: Option<Expr> = from_json_string(lhs_json);
            let rhs: Option<Expr> = from_json_string(rhs_json);
            match (lhs, rhs) {
                (Some(l), Some(r)) => {
                    let expr = Expr::$rs_func(&l, &r);
                    to_json_string(&expr)
                },
                _ => std::ptr::null_mut(),
            }
        }
    };
}

/// Macro for generating Bincode-based FFI n-ary functions.

macro_rules! bincode_nary_ffi {
    ($ffi_name:ident, $rs_func:ident) => {
        #[unsafe(no_mangle)]
        /// # Safety
        ///
        /// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
        /// The caller must ensure:
        /// 1. All pointer arguments are valid and point to initialized memory.
        /// 2. The memory layout of passed structures matches the expected C-ABI layout.
        /// 3. Any pointers returned by this function are managed according to the API's ownership rules.
        #[doc = concat!("FFI wrapper for `", stringify!($rs_func), "`.")]
        #[doc = ""]
        #[doc = concat!("This function takes a vector of bincode buffers, parses them into `Expr`, calls `Expr::", stringify!($rs_func), "`, and returns the result as a bincode buffer.")]
        pub extern "C" fn $ffi_name(
            args_buf: BincodeBuffer
        ) -> BincodeBuffer {
            let args: Option<Vec<Expr>> = from_bincode_buffer(&args_buf);
            if let Some(a) = args {
                let expr = Expr::$rs_func(a);
                to_bincode_buffer(&expr)
            } else {
                BincodeBuffer::empty()
            }
        }
    };
}

/// Macro for generating JSON-based FFI n-ary functions.

macro_rules! json_nary_ffi {
    ($ffi_name:ident, $rs_func:ident) => {
        #[unsafe(no_mangle)]
        /// # Safety
        ///
        /// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
        /// The caller must ensure:
        /// 1. All pointer arguments are valid and point to initialized memory.
        /// 2. The memory layout of passed structures matches the expected C-ABI layout.
        /// 3. Any pointers returned by this function are managed according to the API's ownership rules.
        #[doc = concat!("FFI wrapper for `", stringify!($rs_func), "`.")]
        #[doc = ""]
        #[doc = concat!("This function takes a vector of JSON strings, parses them into `Expr`, calls `Expr::", stringify!($rs_func), "`, and returns the result as a JSON string.")]
        pub extern "C" fn $ffi_name(
            args_json: *const c_char
        ) -> *mut c_char {
            let args: Option<Vec<Expr>> = from_json_string(args_json);
            if let Some(a) = args {
                let expr = Expr::$rs_func(a);
                to_json_string(&expr)
            } else {
                std::ptr::null_mut()
            }
        }
    };
}
