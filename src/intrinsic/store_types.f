            SUBROUTINE store_integer1d_pointer(cpfile_id, ptr_addr, ptr)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=32) :: ptr_addr
               INTEGER, DIMENSION(:) :: ptr
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER :: error, rank

               rank = 1
               dims(1) = size(ptr)

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)

               CALL h5ltmake_dataset_int_f(group_id, ptr_addr, &
                  rank, dims, ptr, error)

               CALL h5gclose_f(group_id, error)

            END SUBROUTINE store_integer1d_pointer
            SUBROUTINE store_integer1d64_pointer(cpfile_id, ptr_addr, ptr)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=32) :: ptr_addr
               INTEGER(8), DIMENSION(:) :: ptr
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER :: error, rank

               rank = 1
               dims(1) = size(ptr)

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)

               !CALL h5ltmake_dataset_long_f(group_id, ptr_addr, &
               !   rank, dims, ptr, error)

               CALL h5gclose_f(group_id, error)

            END SUBROUTINE store_integer1d64_pointer
            SUBROUTINE store_integer2d64_pointer(cpfile_id, ptr_addr, ptr)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=32) :: ptr_addr
               INTEGER(8), DIMENSION(:,:) :: ptr
               INTEGER(HSIZE_T), DIMENSION(2) :: dims
               INTEGER :: error, rank

               rank = 1
               dims = shape(ptr)

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)

               !CALL h5ltmake_dataset_long_f(group_id, ptr_addr, &
               !   rank, dims, ptr, error)

               CALL h5gclose_f(group_id, error)

            END SUBROUTINE store_integer2d64_pointer
            SUBROUTINE store_integer2d_pointer(cpfile_id, ptr_addr, ptr)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=32) :: ptr_addr
               INTEGER, DIMENSION(:,:) :: ptr
               INTEGER(HSIZE_T), DIMENSION(2) :: dims
               INTEGER :: error, rank

               rank = 1
               dims = shape(ptr)

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)

               !CALL h5ltmake_dataset_int_f(group_id, ptr_addr, &
               !   rank, dims, ptr, error)

               CALL h5gclose_f(group_id, error)

            END SUBROUTINE store_integer2d_pointer
            SUBROUTINE store_real1d_pointer(cpfile_id, ptr_addr, ptr)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=32) :: ptr_addr
               REAL(ppm_kind_double), DIMENSION(:) :: ptr
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER :: error, rank

               rank = 1
               dims(1) = size(ptr)

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)

               CALL h5ltmake_dataset_double_f(group_id, ptr_addr, &
                  rank, dims, ptr, error)

               CALL h5gclose_f(group_id, error)

            END SUBROUTINE store_real1d_pointer
            SUBROUTINE store_real2d_pointer(cpfile_id, ptr_addr, ptr)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=32) :: ptr_addr
               REAL(ppm_kind_double), DIMENSION(:,:) :: ptr
               INTEGER(HSIZE_T), DIMENSION(2) :: dims
               INTEGER :: error, rank

               rank = 1
               dims = shape(ptr)

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)

               CALL h5ltmake_dataset_double_f(group_id, ptr_addr, &
                  rank, dims, ptr, error)

               CALL h5gclose_f(group_id, error)

            END SUBROUTINE store_real2d_pointer

