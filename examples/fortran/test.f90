! FFI Example for the rssn library (Final Robust Version)
! ======================================================================
! This version safely handles C-String input/output, memory management,
! and adheres strictly to gfortran's C-binding syntax requirements.
!
! Compilation command (for gfortran, assuming lib is in ./target/release):
! gfortran -o test examples/fortran/test.f90 -L./target/release -lrssn
! ======================================================================

program test_rssn_ffi
    use iso_c_binding
    implicit none

    ! --- 1. FFI Interface Definitions ---
    interface
        ! FFI Function 1: Constructor (JSON Input -> Handle)
        function expr_from_json(json_ptr) bind(C, name='expr_from_json')
            import c_ptr
            type(c_ptr), value :: json_ptr  ! C address passed by value
            type(c_ptr) :: expr_from_json
        end function expr_from_json

        ! FFI Function 2: Serialization (Handle -> C-String Pointer)
        function expr_to_string(handle) bind(C, name='expr_to_string')
            import c_ptr
            type(c_ptr), value :: handle
            type(c_ptr) :: expr_to_string
        end function expr_to_string
        
        ! FFI Function 3: Simplification (Handle -> New Handle)
        function expr_simplify(handle) bind(C, name='expr_simplify')
            import c_ptr
            type(c_ptr), value :: handle
            type(c_ptr) :: expr_simplify
        end function expr_simplify

        ! FFI Function 4: Destructor (Free Handle Memory)
        subroutine expr_free(handle) bind(C, name='expr_free')
            import c_ptr
            type(c_ptr), value :: handle
        end subroutine expr_free

        ! FFI Function 5: Destructor for C-Strings (Free C-String Memory)
        subroutine free_string(s) bind(C, name='free_string')
            import c_ptr
            type(c_ptr), value :: s
        end subroutine free_string

        ! FFI Function 6: Unify Expression (Handle Input，JSON C Pointer Output)
        function expr_unify_expression(handle) bind(C, name='expr_unify_expression')
            import c_ptr
            type(c_ptr), value :: handle
            type(c_ptr) :: expr_unify_expression
        end function expr_unify_expression
    end interface

    ! --- 2. Variables and Pointers ---
    
    type(c_ptr) :: expr_handle = c_null_ptr
    type(c_ptr) :: simplified_expr_handle = c_null_ptr
    type(c_ptr) :: c_str_ptr = c_null_ptr 

    ! Fortran pointer for reading C-String output
    character(kind=c_char), pointer :: fortran_str 

    character(kind=c_char, len=100), target :: json_buffer 

    ! Temporary C string buffer (TARGET required for C-string output handling)
    character(kind=c_char, len=200), target :: temp_c_buffer 
    
    ! --- 3. Program Logic ---

    ! 1. Construct the JSON C-String.
    json_buffer = '{"Add":[{"Variable":"x"},{"Mul":[{"Constant":2.0},{"Variable":"x"}]}]}' // c_null_char
    
    ! Create expression handle by passing the C location (address) of the buffer.
    expr_handle = expr_from_json(C_LOC(json_buffer))
    
    if (.not. c_associated(expr_handle)) then
        print *, "Error: Failed to create expression from JSON."
        stop
    end if

    ! 2. Get and print the string representation of the original expression.
    c_str_ptr = expr_to_string(expr_handle)
    ! Associate C string pointer with Fortran pointer
    call c_f_pointer(c_str_ptr, fortran_str) 
    ! Copy C string to Fortran buffer for safe printing
    call c_copy_string(fortran_str, temp_c_buffer)

    print *, "Original expression: ", trim(temp_c_buffer)
    flush(unit=6) 
    call free_string(c_str_ptr)

    ! 3. Simplify the expression, returning a NEW handle.
    simplified_expr_handle = expr_simplify(expr_handle)
    if (.not. c_associated(simplified_expr_handle)) then
        print *, "Error: Failed to simplify expression."
        call expr_free(expr_handle)
        stop
    end if

    ! 4. Get and print the string representation of the simplified expression.
    c_str_ptr = expr_to_string(simplified_expr_handle)
    
    ! Associate C string pointer with Fortran pointer
    call c_f_pointer(c_str_ptr, fortran_str) 
    ! Copy C string to Fortran buffer for safe printing
    call c_copy_string(fortran_str, temp_c_buffer)

    print *, "Simplified expression: ", trim(temp_c_buffer)
    flush(unit=6) 
    call free_string(c_str_ptr)

    ! --- 5. Unify Units in the Simplified Expression (Core Feature 1) ---
    print *, " "
    print *, "--- Core Feature Validation: Unit Unification ---"
    
    ! 5a. Call the unification function, which returns a FfiResult JSON string
    c_str_ptr = expr_unify_expression(simplified_expr_handle)
    
    if (.not. c_associated(c_str_ptr)) then
        print *, "Error: expr_unify_expression returned NULL pointer."
        call expr_free(expr_handle)
        call expr_free(simplified_expr_handle)
        stop
    end if

    ! 5b. Safely get the JSON output (FfiResult)
    call c_f_pointer(c_str_ptr, fortran_str) 
    call c_copy_string(fortran_str, temp_c_buffer)

    print *, "FfiResult JSON from Unification: ", trim(temp_c_buffer)
    print *, " (User would parse this JSON to get the result or error.)"
    flush(unit=6) 
    
    ! CRITICAL: Free the C-String memory allocated by Rust
    call free_string(c_str_ptr)

    ! 6. Clean up all handles created by the Rust library.
    call expr_free(expr_handle)
    call expr_free(simplified_expr_handle)

    print *, "Fortran example finished successfully."
    flush(unit=6)

contains
    subroutine c_copy_string(c_str, f_str)
        use iso_c_binding
        implicit none
        character(kind=c_char, len=1), pointer :: c_str 
        character(kind=c_char, len=*) :: f_str
        integer :: i, len_c
        
        len_c = 0
        do i = 1, len(f_str)
            if (c_str(i:i) == c_null_char) exit
            len_c = len_c + 1
        end do
        
        f_str = ' '  
        if (len_c > 0) then
            f_str(1:len_c) = c_str(1:len_c)
        end if
    end subroutine c_copy_string

end program test_rssn_ffi