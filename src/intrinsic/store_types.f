            SUBROUTINE store_integer1d_pointer(cpfile_id, ptr_addr, ptr)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=32) :: ptr_addr
               INTEGER, DIMENSION(:) :: ptr
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER :: error, rank
               LOGICAL :: link_exist

               CALL h5lexists_f(cpfile_id, 'intrinsic/'//ptr_addr, &
                  link_exist, error)
               IF (link_exist) THEN
                  RETURN
               END IF

               rank = 1
               dims(1) = size(ptr)

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)

               CALL h5ltmake_dataset_int_f(group_id, ptr_addr, &
                  rank, dims, ptr, error)

               CALL h5gclose_f(group_id, error)

            END SUBROUTINE store_integer1d_pointer
            SUBROUTINE store_integer2d_pointer(cpfile_id, ptr_addr, ptr)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=32) :: ptr_addr
               INTEGER, DIMENSION(:,:) :: ptr
               INTEGER(HSIZE_T), DIMENSION(2) :: dims
               INTEGER :: error, rank
               LOGICAL :: link_exist

               CALL h5lexists_f(cpfile_id, 'intrinsic/'//ptr_addr, &
                  link_exist, error)
               IF (link_exist) THEN
                  RETURN
               END IF
               rank = 1
               dims = shape(ptr)

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)

               !CALL h5ltmake_dataset_int_f(group_id, ptr_addr, &
               !   rank, dims, ptr, error)

               CALL h5gclose_f(group_id, error)

            END SUBROUTINE store_integer2d_pointer
            SUBROUTINE store_integer64_1d_pointer(cpfile_id, ptr_addr, ptr)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id, space_id, dset_id, type_id
               CHARACTER(LEN=32) :: ptr_addr
               INTEGER(8), DIMENSION(:), POINTER :: ptr
               !INTEGER(8), DIMENSION(1), TARGET :: test = (/10/)
               INTEGER(HSIZE_T), DIMENSION(1) :: dims = (/1/)
               INTEGER :: error, rank
               LOGICAL :: link_exist
               TYPE(C_PTR) :: f_ptr

               CALL h5lexists_f(cpfile_id, 'intrinsic/'//ptr_addr, &
                  link_exist, error)
               IF (link_exist) THEN
                  RETURN
               END IF

               rank = 1
               dims = shape(ptr)

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)

               CALL h5screate_simple_f(rank, dims, space_id, error)
               CALL h5dcreate_f(group_id, ptr_addr, H5T_STD_I64LE, &
                  space_id, dset_id, error)

               f_ptr = C_LOC(ptr(1))
               CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, f_ptr, error)

               CALL h5dclose_f(dset_id, error)
               CALL h5sclose_f(space_id, error)

               CALL h5gclose_f(group_id, error)

            END SUBROUTINE store_integer64_1d_pointer
            SUBROUTINE store_integer64_2d_pointer(cpfile_id, ptr_addr, ptr)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id, dset_id, space_id
               CHARACTER(LEN=32) :: ptr_addr
               INTEGER(8), DIMENSION(:,:), POINTER :: ptr
               INTEGER(HSIZE_T), DIMENSION(2) :: dims
               INTEGER :: error, rank
               LOGICAL :: link_exist
               TYPE(C_PTR) :: f_ptr

               CALL h5lexists_f(cpfile_id, 'intrinsic/'//ptr_addr, &
                  link_exist, error)
               IF (link_exist) THEN
                  RETURN
               END IF

               rank = 2
               dims = shape(ptr)

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)

               CALL h5screate_simple_f(rank, dims, space_id, error)
               CALL h5dcreate_f(group_id, ptr_addr, H5T_STD_I64LE, &
                  space_id, dset_id, error)

               f_ptr = C_LOC(ptr(1,1))
               CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, f_ptr, error)

               CALL h5dclose_f(dset_id, error)
               CALL h5sclose_f(space_id, error)

               CALL h5gclose_f(group_id, error)

            END SUBROUTINE store_integer64_2d_pointer
            SUBROUTINE store_real1d_pointer(cpfile_id, ptr_addr, ptr)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=32) :: ptr_addr
               REAL(ppm_kind_double), DIMENSION(:) :: ptr
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER :: error, rank
               LOGICAL :: link_exist

               CALL h5lexists_f(cpfile_id, 'intrinsic/'//ptr_addr, &
                  link_exist, error)
               IF (link_exist) THEN
                  RETURN
               END IF

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
               LOGICAL :: link_exist

               CALL h5lexists_f(cpfile_id, 'intrinsic/'//ptr_addr, &
                  link_exist, error)
               IF (link_exist) THEN
                  RETURN
               END IF

               rank = 1
               dims = shape(ptr)

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)

               CALL h5ltmake_dataset_double_f(group_id, ptr_addr, &
                  rank, dims, ptr, error)

               CALL h5gclose_f(group_id, error)

            END SUBROUTINE store_real2d_pointer

            SUBROUTINE store_logical1d_pointer(cpfile_id, ptr_addr, ptr)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=32) :: ptr_addr
               LOGICAL, DIMENSION(:) :: ptr
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER :: error, rank
               LOGICAL :: link_exist

               CALL h5lexists_f(cpfile_id, 'intrinsic/'//ptr_addr, &
                  link_exist, error)
               IF (link_exist) THEN
                  RETURN
               END IF

               rank = 1
               dims = shape(ptr)

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)

               CALL write_logical_array(group_id, ptr_addr, &
                  ptr, length)
               !CALL h5ltmake_dataset_double_f(group_id, ptr_addr, &
               !   rank, dims, ptr, error)

               CALL h5gclose_f(group_id, error)

            END SUBROUTINE store_logical1d_pointer

            SUBROUTINE store_logical2d_pointer(cpfile_id, ptr_addr, ptr)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=32) :: ptr_addr
               LOGICAL, DIMENSION(:,:) :: ptr
               INTEGER(HSIZE_T), DIMENSION(2) :: dims
               INTEGER :: error, rank
               LOGICAL :: link_exist

               CALL h5lexists_f(cpfile_id, 'intrinsic/'//ptr_addr, &
                  link_exist, error)
               IF (link_exist) THEN
                  RETURN
               END IF

               rank = 1
               dims = shape(ptr)

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)

               CALL write_logical_array_2d(group_id, ptr_addr, &
                  ptr, dims)

               CALL h5gclose_f(group_id, error)

            END SUBROUTINE store_logical2d_pointer
            SUBROUTINE store_complex1d_pointer(cpfile_id, ptr_addr, ptr)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id, space_id, dset_id
               CHARACTER(LEN=32) :: ptr_addr
               COMPLEX(8), DIMENSION(:) :: ptr
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER :: error, rank
               LOGICAL :: link_exist

               CALL h5lexists_f(cpfile_id, 'intrinsic/'//ptr_addr, &
                  link_exist, error)
               IF (link_exist) THEN
                  RETURN
               END IF

               rank = 1
               dims(1) = size(ptr)

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)


               CALL h5gclose_f(group_id, error)

            END SUBROUTINE store_complex1d_pointer
            SUBROUTINE store_complex2d_pointer(cpfile_id, ptr_addr, ptr)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=32) :: ptr_addr
               COMPLEX(8), DIMENSION(:,:) :: ptr
               INTEGER(HSIZE_T), DIMENSION(2) :: dims
               INTEGER :: error, rank
               LOGICAL :: link_exist

               CALL h5lexists_f(cpfile_id, 'intrinsic/'//ptr_addr, &
                  link_exist, error)
               IF (link_exist) THEN
                  RETURN
               END IF

               rank = 1
               dims = shape(ptr)

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)

               !CALL h5ltmake_dataset_long_f(group_id, ptr_addr, &
               !   rank, dims, ptr, error)

               CALL h5gclose_f(group_id, error)

            END SUBROUTINE store_complex2d_pointer
