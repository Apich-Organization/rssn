#[macro_export]

macro_rules! json_ffi_unary {
    ($name:ident, $input_type:ty, | $arg:ident | $body:expr) => {
        #[no_mangle]

        pub extern "C" fn $name(input_json: *const std::ffi::c_char) -> *mut std::ffi::c_char {

            let input: Option<$input_type> = $crate::ffi_apis::common::from_json_string(input_json);

            if let Some($arg) = input {

                let result = $body;

                $crate::ffi_apis::common::to_json_string(&result)
            } else {

                std::ptr::null_mut()
            }
        }
    };
}

#[macro_export]

macro_rules! json_ffi_binary {
    ($name:ident, $input1_type:ty, $input2_type:ty, | $arg1:ident, $arg2:ident | $body:expr) => {
        #[no_mangle]

        pub extern "C" fn $name(
            input1_json: *const std::ffi::c_char,
            input2_json: *const std::ffi::c_char,
        ) -> *mut std::ffi::c_char {

            let input1: Option<$input1_type> =
                $crate::ffi_apis::common::from_json_string(input1_json);

            let input2: Option<$input2_type> =
                $crate::ffi_apis::common::from_json_string(input2_json);

            if let (Some($arg1), Some($arg2)) = (input1, input2) {

                let result = $body;

                $crate::ffi_apis::common::to_json_string(&result)
            } else {

                std::ptr::null_mut()
            }
        }
    };
}

#[macro_export]

macro_rules! handle_ffi_unary {
    ($name:ident, $input_type:ty, | $arg:ident | $body:expr) => {
        #[no_mangle]

        pub extern "C" fn $name(input: *const $input_type) -> *mut $crate::symbolic::core::Expr {

            let $arg = unsafe {

                &*input
            };

            let result = $body;

            Box::into_raw(Box::new(result))
        }
    };
    // Generic return type version
    ($name:ident, $input_type:ty, $ret_type:ty, | $arg:ident | $body:expr) => {
        #[no_mangle]

        pub extern "C" fn $name(input: *const $input_type) -> *mut $ret_type {

            let $arg = unsafe {

                &*input
            };

            let result = $body;

            Box::into_raw(Box::new(result))
        }
    };
}

#[macro_export]

macro_rules! handle_ffi_binary {
    (
        $name:ident,
        $input1_type:ty,
        $input2_type:ty,
        $ret_type:ty, |
        $arg1:ident,
        $arg2:ident |
        $body:expr
    ) => {
        #[no_mangle]

        pub extern "C" fn $name(
            input1: *const $input1_type,
            input2: *const $input2_type,
        ) -> *mut $ret_type {

            let $arg1 = unsafe {

                &*input1
            };

            let $arg2 = unsafe {

                &*input2
            };

            let result = $body;

            Box::into_raw(Box::new(result))
        }
    };
}

#[macro_export]

macro_rules! bincode_ffi_unary {
    ($name:ident, $input_type:ty, | $arg:ident | $body:expr) => {
        #[no_mangle]

        pub extern "C" fn $name(
            input_buf: $crate::ffi_apis::common::BincodeBuffer
        ) -> $crate::ffi_apis::common::BincodeBuffer {

            let input: Option<$input_type> =
                $crate::ffi_apis::common::from_bincode_buffer(&input_buf);

            if let Some($arg) = input {

                let result = $body;

                $crate::ffi_apis::common::to_bincode_buffer(&result)
            } else {

                $crate::ffi_apis::common::BincodeBuffer::empty()
            }
        }
    };
}

#[macro_export]

macro_rules! bincode_ffi_binary {
    ($name:ident, $input1_type:ty, $input2_type:ty, | $arg1:ident, $arg2:ident | $body:expr) => {
        #[no_mangle]

        pub extern "C" fn $name(
            input1_buf: $crate::ffi_apis::common::BincodeBuffer,
            input2_buf: $crate::ffi_apis::common::BincodeBuffer,
        ) -> $crate::ffi_apis::common::BincodeBuffer {

            let input1: Option<$input1_type> =
                $crate::ffi_apis::common::from_bincode_buffer(&input1_buf);

            let input2: Option<$input2_type> =
                $crate::ffi_apis::common::from_bincode_buffer(&input2_buf);

            if let (Some($arg1), Some($arg2)) = (input1, input2) {

                let result = $body;

                $crate::ffi_apis::common::to_bincode_buffer(&result)
            } else {

                $crate::ffi_apis::common::BincodeBuffer::empty()
            }
        }
    };
}
