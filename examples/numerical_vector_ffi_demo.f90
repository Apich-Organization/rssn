! FFI Example for Numerical Vector Operations
! ======================================================================
! Compilation command (for gfortran):
! gfortran -o vector_demo examples/numerical_vector_ffi_demo.f90 -L./target/debug -lrssn
! ======================================================================

program vector_demo
    use iso_c_binding
    implicit none

    ! --- FFI Interfaces ---
    interface
        function rssn_num_vec_create(data_ptr, n_len) bind(C, name='rssn_num_vec_create')
            import c_ptr, c_double, c_size_t
            type(c_ptr), value :: data_ptr
            integer(c_size_t), value :: n_len
            type(c_ptr) :: rssn_num_vec_create
        end function rssn_num_vec_create

        subroutine rssn_num_vec_free(v_ptr) bind(C, name='rssn_num_vec_free')
            import c_ptr
            type(c_ptr), value :: v_ptr
        end subroutine rssn_num_vec_free

        function rssn_num_vec_len(v_ptr) bind(C, name='rssn_num_vec_len')
            import c_ptr, c_size_t
            type(c_ptr), value :: v_ptr
            integer(c_size_t) :: rssn_num_vec_len
        end function rssn_num_vec_len

        function rssn_num_vec_data(v_ptr) bind(C, name='rssn_num_vec_data')
            import c_ptr
            type(c_ptr), value :: v_ptr
            type(c_ptr) :: rssn_num_vec_data
        end function rssn_num_vec_data

        function rssn_num_vec_add(v1, v2) bind(C, name='rssn_num_vec_add')
            import c_ptr
            type(c_ptr), value :: v1, v2
            type(c_ptr) :: rssn_num_vec_add
        end function rssn_num_vec_add

        function rssn_num_vec_dot_product(v1, v2, res) bind(C, name='rssn_num_vec_dot_product')
            import c_ptr, c_double, c_int
            type(c_ptr), value :: v1, v2
            type(c_ptr), value :: res
            integer(c_int) :: rssn_num_vec_dot_product
        end function rssn_num_vec_dot_product
    end interface

    ! --- Variables ---
    real(c_double), target :: d1(3) = [1.0_c_double, 2.0_c_double, 3.0_c_double]
    real(c_double), target :: d2(3) = [4.0_c_double, 5.0_c_double, 6.0_c_double]
    real(c_double), target :: dot_res
    type(c_ptr) :: v1, v2, v_sum
    integer(c_size_t) :: res_len
    type(c_ptr) :: res_data_ptr
    real(c_double), pointer :: res_array(:)
    integer :: i

    print *, "--- Fortran FFI Demo: Numerical Vectors ---"

    ! Create vectors
    v1 = rssn_num_vec_create(c_loc(d1), 3_c_size_t)
    v2 = rssn_num_vec_create(c_loc(d2), 3_c_size_t)

    ! Add vectors
    v_sum = rssn_num_vec_add(v1, v2)
    
    if (c_associated(v_sum)) then
        res_len = rssn_num_vec_len(v_sum)
        res_data_ptr = rssn_num_vec_data(v_sum)
        call c_f_pointer(res_data_ptr, res_array, [res_len])
        
        print *, "v1 + v2 = [", (res_array(i), i=1, int(res_len)), "]"
    end if

    ! Dot product
    i = rssn_num_vec_dot_product(v1, v2, c_loc(dot_res))
    if (i == 0) then
        print *, "v1 . v2 = ", dot_res
    end if

    ! Cleanup
    call rssn_num_vec_free(v1)
    call rssn_num_vec_free(v2)
    call rssn_num_vec_free(v_sum)

end program vector_demo
