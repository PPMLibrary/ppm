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

            END SUBROUTINE

