            SUBROUTINE read_integer1d_pointer(cpfile_id, ptr_addr, input)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id, dset_id, space_id
               CHARACTER(LEN=32) :: ptr_addr
               INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: input
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER :: error, rank, tclass
               INTEGER(SIZE_T) :: tsize
               LOGICAL :: link_exist
               CLASS(intrinsic_tree), POINTER :: ptr_tree

               rank = 1

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)
               CALL h5ltget_dataset_info_f(group_id, ptr_addr, dims, &
                  tclass, tsize, error)
               ALLOCATE (input(int(dims(1))))
               CALL h5dopen_f(group_id, ptr_addr, dset_id, error)
               CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, &
                  input, dims, error)
               CALL h5dclose_f(dset_id, error)
               !CALL h5ltread_dataset_f(group_id, "test", &
               !   H5T_NATIVE_INTEGER, &
               !   dims, input, error)

               ALLOCATE(integer1d_tree::ptr_tree)
               SELECT TYPE (ptr_tree)
               TYPE is (integer1d_tree)
                  ptr_tree%val => input
                  ptr_tree%hash = FNVHash(input, int(dims(1)))
                  ptr_tree%key = ptr_addr
               END SELECT

               !CALL pointer_insert(pointer_data, ptr_tree, ptr_addr)
               CALL insert_intrinsic(pointer_data%itree, ptr_tree)

               CALL h5gclose_f(group_id, error)

            END SUBROUTINE read_integer1d_pointer
            SUBROUTINE read_integer2d_pointer(cpfile_id, ptr_addr, input)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id, dset_id, space_id
               CHARACTER(LEN=32) :: ptr_addr
               INTEGER, DIMENSION(:), ALLOCATABLE :: int_buf
               INTEGER, DIMENSION(:,:), ALLOCATABLE, TARGET :: input
               INTEGER(HSIZE_T), DIMENSION(2) :: dims
               INTEGER :: error, rank, tclass, bsize
               INTEGER(SIZE_T) :: tsize
               LOGICAL :: link_exist
               CLASS(intrinsic_tree), POINTER :: ptr_tree

               rank = 2

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)
               CALL h5ltget_dataset_info_f(group_id, ptr_addr, dims, &
                  tclass, tsize, error)
               ALLOCATE (input(int(dims(1)),int(dims(2))))
               CALL h5dopen_f(group_id, ptr_addr, dset_id, error)
               CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, &
                  input, dims, error)
               CALL h5dclose_f(dset_id, error)
               !CALL h5ltread_dataset_f(group_id, "test", &
               !   H5T_NATIVE_INTEGER, &
               !   dims, input, error)

               bsize = size(transfer(input, int_buf))
               ALLOCATE (int_buf(bsize))
               int_buf = transfer(input,int_buf)
               ALLOCATE(integer2d_tree::ptr_tree)
               SELECT TYPE (ptr_tree)
               TYPE is (integer2d_tree)
                  ptr_tree%val => input
                  ptr_tree%hash = FNVHash(int_buf, bsize)
                  ptr_tree%key = ptr_addr
               END SELECT
               DEALLOCATE(int_buf)

               !CALL pointer_insert(pointer_data, ptr_tree, ptr_addr)
               CALL insert_intrinsic(pointer_data%itree, ptr_tree)

               CALL h5gclose_f(group_id, error)

            END SUBROUTINE read_integer2d_pointer
            SUBROUTINE read_REAL1d_pointer(cpfile_id, ptr_addr, input)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id, dset_id, space_id
               CHARACTER(LEN=32) :: ptr_addr
               REAL(ppm_kind_double), DIMENSION(:), ALLOCATABLE, TARGET :: input
               INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: int_buf
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER :: error, rank, tclass, bsize
               INTEGER(SIZE_T) :: tsize
               LOGICAL :: link_exist
               CLASS(intrinsic_tree), POINTER :: ptr_tree

               rank = 1

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)
               CALL h5ltget_dataset_info_f(group_id, ptr_addr, dims, &
                  tclass, tsize, error)
               ALLOCATE (input(int(dims(1))))
               CALL h5dopen_f(group_id, ptr_addr, dset_id, error)
               CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, &
                  input, dims, error)
               CALL h5dclose_f(dset_id, error)
               !CALL h5ltread_dataset_f(group_id, "test", &
               !   H5T_NATIVE_INTEGER, &
               !   dims, input, error)

               bsize = size(transfer(input, int_buf))
               ALLOCATE (int_buf(bsize))
               int_buf = transfer(input,int_buf)
               ALLOCATE(real1d_tree::ptr_tree)
               SELECT TYPE (ptr_tree)
               TYPE is (real1d_tree)
                  ptr_tree%val => input
                  ptr_tree%hash = FNVHash(int_buf, int(dims(1)))
                  ptr_tree%key = ptr_addr
               END SELECT

               !CALL pointer_insert(pointer_data, ptr_tree, ptr_addr)
               CALL insert_intrinsic(pointer_data%itree, ptr_tree)

               CALL h5gclose_f(group_id, error)

            END SUBROUTINE read_real1d_pointer
            SUBROUTINE read_real2d_pointer(cpfile_id, ptr_addr, input)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id, dset_id, space_id
               CHARACTER(LEN=32) :: ptr_addr
               INTEGER, DIMENSION(:), ALLOCATABLE :: int_buf
               REAL(ppm_kind_double), DIMENSION(:,:), ALLOCATABLE, TARGET :: input
               INTEGER(HSIZE_T), DIMENSION(2) :: dims
               INTEGER :: error, rank, tclass, bsize
               INTEGER(SIZE_T) :: tsize
               LOGICAL :: link_exist
               CLASS(intrinsic_tree), POINTER :: ptr_tree

               rank = 2

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)
               CALL h5ltget_dataset_info_f(group_id, ptr_addr, dims, &
                  tclass, tsize, error)
               ALLOCATE (input(int(dims(1)),int(dims(2))))
               CALL h5dopen_f(group_id, ptr_addr, dset_id, error)
               CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, &
                  input, dims, error)
               CALL h5dclose_f(dset_id, error)
               !CALL h5ltread_dataset_f(group_id, "test", &
               !   H5T_NATIVE_INTEGER, &
               !   dims, input, error)

               bsize = size(transfer(input, int_buf))
               ALLOCATE (int_buf(bsize))
               int_buf = transfer(input,int_buf)
               ALLOCATE(real2d_tree::ptr_tree)
               SELECT TYPE (ptr_tree)
               TYPE is (real2d_tree)
                  ptr_tree%val => input
                  ptr_tree%hash = FNVHash(int_buf, bsize)
                  ptr_tree%key = ptr_addr
               END SELECT
               DEALLOCATE(int_buf)

               !CALL pointer_insert(pointer_data, ptr_tree, ptr_addr)
               CALL insert_intrinsic(pointer_data%itree, ptr_tree)

               CALL h5gclose_f(group_id, error)

            END SUBROUTINE read_real2d_pointer
            SUBROUTINE READ_integer64_1d_pointer(cpfile_id, ptr_addr, ptr)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id, space_id, dset_id, type_id
               CHARACTER(LEN=32) :: ptr_addr
               INTEGER(8), DIMENSION(:), POINTER :: ptr
               INTEGER, DIMENSION(:), POINTER :: int_buf
               INTEGER(HSIZE_T) :: tsize
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER :: error, rank, bsize, tclass
               TYPE(C_PTR) :: f_ptr
               CLASS(intrinsic_tree), POINTER :: ptr_tree

               !lookup if pointer exists in tree

               rank = 1
               dims = shape(ptr)

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)

               CALL h5ltget_dataset_info_f(group_id, ptr_addr, dims, &
                  tclass, tsize, error)
               ALLOCATE (ptr(int(dims(1))))

               CALL h5dOPEN_f(group_id, ptr_addr, &
                  dset_id, error)

               f_ptr = C_LOC(ptr(1))
               CALL h5dREAD_f(dset_id, H5T_NATIVE_INTEGER, f_ptr, error)

               CALL h5dclose_f(dset_id, error)

               CALL h5gclose_f(group_id, error)

               bsize = size(transfer(ptr, int_buf))
               ALLOCATE (int_buf(bsize))
               int_buf = transfer(ptr,int_buf)
               ALLOCATE(integer64_1d_tree::ptr_tree)
               SELECT TYPE (ptr_tree)
               TYPE is (integer64_1d_tree)
                  ptr_tree%val => ptr
                  ptr_tree%hash = FNVHash(int_buf, bsize)
                  ptr_tree%key = ptr_addr
               END SELECT
               DEALLOCATE(int_buf)

               CALL insert_intrinsic(pointer_data%itree, ptr_tree)
               ! record the pointer in the tree

            END SUBROUTINE READ_integer64_1d_pointer
            SUBROUTINE READ_integer64_2d_pointer(cpfile_id, ptr_addr, ptr)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id, space_id, dset_id, type_id
               CHARACTER(LEN=32) :: ptr_addr
               INTEGER(8), DIMENSION(:,:), POINTER :: ptr
               INTEGER, DIMENSION(:), POINTER :: int_buf
               INTEGER(HSIZE_T) :: tsize
               INTEGER(HSIZE_T), DIMENSION(2) :: dims
               INTEGER :: error, rank, bsize, tclass
               TYPE(C_PTR) :: f_ptr
               CLASS(intrinsic_tree), POINTER :: ptr_tree

               !lookup if pointer exists in tree

               rank = 2
               dims = shape(ptr)

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)

               CALL h5ltget_dataset_info_f(group_id, ptr_addr, dims, &
                  tclass, tsize, error)
               ALLOCATE (ptr(int(dims(1)),int(dims(2))))

               CALL h5dOPEN_f(group_id, ptr_addr, &
                  dset_id, error)

               f_ptr = C_LOC(ptr(1,1))
               CALL h5dREAD_f(dset_id, H5T_NATIVE_INTEGER, f_ptr, error)

               CALL h5dclose_f(dset_id, error)

               CALL h5gclose_f(group_id, error)

               bsize = size(transfer(ptr, int_buf))
               ALLOCATE (int_buf(bsize))
               int_buf = transfer(ptr,int_buf)
               ALLOCATE(integer64_2d_tree::ptr_tree)
               SELECT TYPE (ptr_tree)
               TYPE is (integer64_2d_tree)
                  ptr_tree%val => ptr
                  ptr_tree%hash = FNVHash(int_buf, bsize)
                  ptr_tree%key = ptr_addr
               END SELECT
               DEALLOCATE(int_buf)

               CALL insert_intrinsic(pointer_data%itree, ptr_tree)
               ! record the pointer in the tree

            END SUBROUTINE READ_integer64_2d_pointer
            SUBROUTINE READ_logical1d_pointer(cpfile_id, ptr_addr, ptr)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=32) :: ptr_addr
               LOGICAL, DIMENSION(:), POINTER :: ptr
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER :: error, rank, bsize, tclass
               INTEGER, DIMENSION(:), POINTER :: int_buf
               INTEGER(HSIZE_T) :: tsize
               CLASS(intrinsic_tree), POINTER :: ptr_tree

               rank = 1
               dims = shape(ptr)

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)

               CALL h5ltget_dataset_info_f(group_id, ptr_addr, dims, &
                  tclass, tsize, error)
               ALLOCATE (ptr(int(dims(1))))

               CALL READ_logical_array(group_id, ptr_addr, &
                  ptr, int(dims(1)))

               CALL h5gclose_f(group_id, error)

               bsize = size(transfer(ptr, int_buf))
               ALLOCATE (int_buf(bsize))
               int_buf = transfer(ptr,int_buf)
               ALLOCATE(logical1d_tree::ptr_tree)
               SELECT TYPE (ptr_tree)
               TYPE is (logical1d_tree)
                  ptr_tree%val => ptr
                  ptr_tree%hash = FNVHash(int_buf, bsize)
                  ptr_tree%key = ptr_addr
               END SELECT
               DEALLOCATE(int_buf)

               CALL insert_intrinsic(pointer_data%itree, ptr_tree)

            END SUBROUTINE READ_logical1d_pointer
            SUBROUTINE READ_logical2d_pointer(cpfile_id, ptr_addr, ptr)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=32) :: ptr_addr
               LOGICAL, DIMENSION(:,:), POINTER :: ptr
               INTEGER(HSIZE_T), DIMENSION(2) :: dims
               INTEGER :: error, rank, bsize, tclass
               INTEGER, DIMENSION(:), POINTER :: int_buf
               INTEGER(HSIZE_T) :: tsize
               CLASS(intrinsic_tree), POINTER :: ptr_tree

               rank = 2
               dims = shape(ptr)

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)

               CALL h5ltget_dataset_info_f(group_id, ptr_addr, dims, &
                  tclass, tsize, error)
               ALLOCATE (ptr(int(dims(1)),int(dims(2))))

               CALL READ_logical_array_2d(group_id, ptr_addr, &
                  ptr, dims)

               CALL h5gclose_f(group_id, error)

               bsize = size(transfer(ptr, int_buf))
               ALLOCATE (int_buf(bsize))
               int_buf = transfer(ptr,int_buf)
               ALLOCATE(logical2d_tree::ptr_tree)
               SELECT TYPE (ptr_tree)
               TYPE is (logical2d_tree)
                  ptr_tree%val => ptr
                  ptr_tree%hash = FNVHash(int_buf, bsize)
                  ptr_tree%key = ptr_addr
               END SELECT
               DEALLOCATE(int_buf)

               CALL insert_intrinsic(pointer_data%itree, ptr_tree)

            END SUBROUTINE READ_logical2d_pointer

!            SUBROUTINE READ_complex1d_pointer(cpfile_id, ptr_addr, ptr)
!               IMPLICIT NONE
!               INTEGER(HID_T), INTENT(IN) :: cpfile_id
!               INTEGER(HID_T) :: group_id
!               CHARACTER(LEN=32) :: ptr_addr
!               COMPLEX(8), DIMENSION(:) :: ptr
!               INTEGER(HSIZE_T), DIMENSION(2) :: dims
!               INTEGER :: error, rank, bsize, tclass, i
!               INTEGER, DIMENSION(:), POINTER :: int_buf
!               INTEGER(HSIZE_T) :: tsize
!               CLASS(intrinsic_tree), POINTER :: ptr_tree
!
!               INTEGER(HID_T) :: group_id, space_id, dset_id, type_id, &
!                  sub_id, array_id
!               INTEGER(HSIZE_T) :: sizef, offset
!               REAL(8), DIMENSION(:), ALLOCATABLE :: buffer
!
!               rank = 1
!               dims(1) = size(ptr)
!
!               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)
!
!               offset = 0
!               CALL h5tcreate_f(H5T_COMPOUND_F, sizef*2, type_id, error)
!               CALL h5tinsert_f(type_id, "real", offset, &
!                  H5T_NATIVE_REAL, error)
!               offset = offset + sizef
!               CALL h5tinsert_f(type_id, "im", sizef, &
!                  H5T_NATIVE_REAL, error)
!
!               CALL h5tarray_create_f(type_id, rank, dims, array_id, &
!                     error)
!
!               ! store the real parts
!               ALLOCATE(buffer(int(dims(1))))
!               DO i=1, int(dims(1))
!                  buffer(i) = real(ptr(i))
!               ENDDO
!               offset = 0
!               CALL h5tcreate_f(H5T_COMPOUND_F, sizef, sub_id, error)
!               CALL h5tinsert_f(sub_id, "real", sizef, H5T_NATIVE_REAL, error)
!               CALL h5dwrite_f(dset_id, sub_id, buffer, dims, error)
!               CALL h5tclose_f(sub_id, error)
!
!               ! Now the imaginary parts
!               DO i=1, int(dims(1))
!                  buffer(i) = aimag(ptr(i))
!               ENDDO
!               offset = 0
!               CALL h5tcreate_f(H5T_COMPOUND_F, sizef, sub_id, error)
!               CALL h5tinsert_f(sub_id, "im", sizef, H5T_NATIVE_REAL, error)
!               CALL h5dwrite_f(dset_id, sub_id, buffer, dims, error)
!               CALL h5tclose_f(sub_id, error)
!               DEALLOCATE (buffer)
!
!               CALL h5tclose_f(sub_id, error)
!               CALL h5tclose_f(type_id, error)
!
!               CALL h5gclose_f(group_id, error)
!
!            END SUBROUTINE store_complex1d_pointer
!            SUBROUTINE store_complex2d_pointer(cpfile_id, ptr_addr, ptr)
!               INTEGER(HID_T), INTENT(IN) :: cpfile_id
!               INTEGER(HID_T) :: group_id, space_id, dset_id, type_id, &
!                  sub_id, array_id
!               CHARACTER(LEN=32) :: ptr_addr
!               COMPLEX(8), DIMENSION(:,:) :: ptr
!               INTEGER(HSIZE_T), DIMENSION(2) :: dims
!               INTEGER(HSIZE_T) :: sizef, offset
!               INTEGER :: error, rank, i,j
!               LOGICAL :: link_exist
!               REAL(8), DIMENSION(:,:), ALLOCATABLE :: buffer
!
!               CALL h5lexists_f(cpfile_id, 'intrinsic/'//ptr_addr, &
!                  link_exist, error)
!               IF (link_exist) THEN
!                  RETURN
!               END IF
!
!               rank = 2
!               dims = shape(ptr)
!
!               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)
!
!               offset = 0
!               CALL h5tcreate_f(H5T_COMPOUND_F, sizef*2, type_id, error)
!               CALL h5tinsert_f(type_id, "real", offset, &
!                  H5T_NATIVE_REAL, error)
!               offset = offset + sizef
!               CALL h5tinsert_f(type_id, "im", sizef, &
!                  H5T_NATIVE_REAL, error)
!
!               CALL h5tarray_create_f(type_id, rank, dims, array_id, &
!                     error)
!
!               ! store the real parts
!               ALLOCATE(buffer(int(dims(1)),int(dims(2))))
!               DO i=1, int(dims(1))
!                  DO j=1, int(dims(2))
!                     buffer(i,j) = real(ptr(i,j))
!                  ENDDO
!               ENDDO
!               offset = 0
!               CALL h5tcreate_f(H5T_COMPOUND_F, sizef, sub_id, error)
!               CALL h5tinsert_f(sub_id, "real", sizef, H5T_NATIVE_REAL, error)
!               CALL h5dwrite_f(dset_id, sub_id, buffer, dims, error)
!               CALL h5tclose_f(sub_id, error)
!
!               ! Now the imaginary parts
!               DO i=1, int(dims(1))
!                  DO j=1, int(dims(2))
!                     buffer(i,j) = aimag(ptr(i,j))
!                  ENDDO
!               ENDDO
!               offset = 0
!               CALL h5tcreate_f(H5T_COMPOUND_F, sizef, sub_id, error)
!               CALL h5tinsert_f(sub_id, "im", sizef, H5T_NATIVE_REAL, error)
!               CALL h5dwrite_f(dset_id, sub_id, buffer, dims, error)
!               CALL h5tclose_f(sub_id, error)
!               DEALLOCATE (buffer)
!
!               CALL h5tclose_f(sub_id, error)
!               CALL h5tclose_f(type_id, error)
!
!               CALL h5gclose_f(group_id, error)
!            END SUBROUTINE store_complex2d_pointer
