            FUNCTION recover_integer1d_pointer(cpfile_id, ptr_addr, &
                  & ptr) RESULT(ptr2)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id, dset_id
               CHARACTER(LEN=32) :: ptr_addr
               INTEGER, DIMENSION(:), POINTER :: ptr, ptr2
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER :: error, rank, tclass
               INTEGER(SIZE_T) :: tsize
               CLASS(intrinsic_tree), POINTER :: ptr_tree

               CALL lookup_pointer(ptr_addr, pointer_data%itree, ptr_tree)
               IF (associated(ptr_tree)) THEN
                  SELECT TYPE(ptr_tree)
                  CLASS is (integer1d_tree)
                     ptr2 => ptr_tree%val
                  END SELECT
                  RETURN
               ENDIF
               rank = 1

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)
               CALL h5ltget_dataset_info_f(group_id, ptr_addr, dims, &
                  tclass, tsize, error)
               ALLOCATE (ptr(int(dims(1))))
               CALL h5dopen_f(group_id, ptr_addr, dset_id, error)
               CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, &
                  ptr, dims, error)
               CALL h5dclose_f(dset_id, error)
               !CALL h5ltread_dataset_f(group_id, "test", &
               !   H5T_NATIVE_INTEGER, &
               !   dims, ptr, error)

               ALLOCATE(integer1d_tree::ptr_tree)
               SELECT TYPE (ptr_tree)
               TYPE is (integer1d_tree)
                  ptr_tree%val => ptr
                  ptr_tree%hash = FNVHash(ptr, int(dims(1)))
                  ptr_tree%key = ptr_addr
               END SELECT

               !CALL pointer_insert(pointer_data, ptr_tree, ptr_addr)
               CALL insert_intrinsic(pointer_data%itree, ptr_tree)

               CALL h5gclose_f(group_id, error)

               ptr2 => ptr
            END FUNCTION recover_integer1d_pointer
            FUNCTION recover_integer2d_pointer(cpfile_id, ptr_addr,&
                & ptr) RESULT(ptr2)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id, dset_id
               CHARACTER(LEN=32) :: ptr_addr
               INTEGER, DIMENSION(:), ALLOCATABLE :: int_buf
               INTEGER, DIMENSION(:,:), POINTER :: ptr, ptr2
               INTEGER(HSIZE_T), DIMENSION(2) :: dims
               INTEGER :: error, rank, tclass, bsize
               INTEGER(SIZE_T) :: tsize
               CLASS(intrinsic_tree), POINTER :: ptr_tree

               CALL lookup_pointer(ptr_addr, pointer_data%itree, ptr_tree)
               IF (associated(ptr_tree)) THEN
                  SELECT TYPE(ptr_tree)
                  CLASS is (integer2d_tree)
                     ptr2 => ptr_tree%val
                  END SELECT
                  RETURN
               ENDIF
               rank = 2

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)
               CALL h5ltget_dataset_info_f(group_id, ptr_addr, dims, &
                  tclass, tsize, error)
               ALLOCATE (ptr(int(dims(1)),int(dims(2))))
               CALL h5dopen_f(group_id, ptr_addr, dset_id, error)
               CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, &
                  ptr, dims, error)
               CALL h5dclose_f(dset_id, error)
               !CALL h5ltread_dataset_f(group_id, "test", &
               !   H5T_NATIVE_INTEGER, &
               !   dims, ptr, error)

               bsize = size(transfer(ptr, int_buf))
               ALLOCATE (int_buf(bsize))
               int_buf = transfer(ptr,int_buf)
               ALLOCATE(integer2d_tree::ptr_tree)
               SELECT TYPE (ptr_tree)
               TYPE is (integer2d_tree)
                  ptr_tree%val => ptr
                  ptr_tree%hash = FNVHash(int_buf, bsize)
                  ptr_tree%key = ptr_addr
               END SELECT
               DEALLOCATE(int_buf)

               !CALL pointer_insert(pointer_data, ptr_tree, ptr_addr)
               CALL insert_intrinsic(pointer_data%itree, ptr_tree)

               CALL h5gclose_f(group_id, error)

               ptr2 => ptr

            END FUNCTION recover_integer2d_pointer
            FUNCTION recover_REAL1d_pointer(cpfile_id, ptr_addr, &
                   & ptr) RESULT(ptr2)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id, dset_id
               CHARACTER(LEN=32) :: ptr_addr
               REAL(ppm_kind_double), DIMENSION(:), POINTER :: ptr, ptr2
               INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: int_buf
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER :: error, rank, tclass, bsize
               INTEGER(SIZE_T) :: tsize
               CLASS(intrinsic_tree), POINTER :: ptr_tree

               CALL lookup_pointer(ptr_addr, pointer_data%itree, ptr_tree)
               IF (associated(ptr_tree)) THEN
                  SELECT TYPE(ptr_tree)
                  CLASS is (real1d_tree)
                     ptr2 => ptr_tree%val
                  END SELECT
                  RETURN
               ENDIF
               rank = 1

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)
               CALL h5ltget_dataset_info_f(group_id, ptr_addr, dims, &
                  tclass, tsize, error)
               ALLOCATE (ptr(int(dims(1))))
               CALL h5dopen_f(group_id, ptr_addr, dset_id, error)
               CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, &
                  ptr, dims, error)
               CALL h5dclose_f(dset_id, error)
               !CALL h5ltread_dataset_f(group_id, "test", &
               !   H5T_NATIVE_INTEGER, &
               !   dims, ptr, error)

               bsize = size(transfer(ptr, int_buf))
               ALLOCATE (int_buf(bsize))
               int_buf = transfer(ptr,int_buf)
               ALLOCATE(real1d_tree::ptr_tree)
               SELECT TYPE (ptr_tree)
               TYPE is (real1d_tree)
                  ptr_tree%val => ptr
                  ptr_tree%hash = FNVHash(int_buf, int(dims(1)))
                  ptr_tree%key = ptr_addr
               END SELECT

               !CALL pointer_insert(pointer_data, ptr_tree, ptr_addr)
               CALL insert_intrinsic(pointer_data%itree, ptr_tree)

               CALL h5gclose_f(group_id, error)
               ptr2 => ptr

            END FUNCTION recover_real1d_pointer
            FUNCTION recover_real2d_pointer(cpfile_id, ptr_addr, &
                     & ptr) RESULT(ptr2)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id, dset_id
               CHARACTER(LEN=32) :: ptr_addr
               INTEGER, DIMENSION(:), ALLOCATABLE :: int_buf
               REAL(ppm_kind_double), DIMENSION(:,:), POINTER :: ptr, ptr2
               INTEGER(HSIZE_T), DIMENSION(2) :: dims
               INTEGER :: error, rank, tclass, bsize
               INTEGER(SIZE_T) :: tsize
               CLASS(intrinsic_tree), POINTER :: ptr_tree

               CALL lookup_pointer(ptr_addr, pointer_data%itree, ptr_tree)
               IF (associated(ptr_tree)) THEN
                  SELECT TYPE(ptr_tree)
                  CLASS is (real2d_tree)
                     ptr2 => ptr_tree%val
                  END SELECT
                  RETURN
               ENDIF
               rank = 2

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)
               CALL h5ltget_dataset_info_f(group_id, ptr_addr, dims, &
                  tclass, tsize, error)
               ALLOCATE (ptr(int(dims(1)),int(dims(2))))
               CALL h5dopen_f(group_id, ptr_addr, dset_id, error)
               CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, &
                  ptr, dims, error)
               CALL h5dclose_f(dset_id, error)
               !CALL h5ltread_dataset_f(group_id, "test", &
               !   H5T_NATIVE_INTEGER, &
               !   dims, ptr, error)

               bsize = size(transfer(ptr, int_buf))
               ALLOCATE (int_buf(bsize))
               int_buf = transfer(ptr,int_buf)
               ALLOCATE(real2d_tree::ptr_tree)
               SELECT TYPE (ptr_tree)
               TYPE is (real2d_tree)
                  ptr_tree%val => ptr
                  ptr_tree%hash = FNVHash(int_buf, bsize)
                  ptr_tree%key = ptr_addr
               END SELECT
               DEALLOCATE(int_buf)

               !CALL pointer_insert(pointer_data, ptr_tree, ptr_addr)
               CALL insert_intrinsic(pointer_data%itree, ptr_tree)

               CALL h5gclose_f(group_id, error)

               ptr2 => ptr
            END FUNCTION recover_real2d_pointer
            FUNCTION recover_integer64_1d_pointer(cpfile_id, ptr_addr, &
                     ptr) RESULT(ptr2)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id, dset_id
               CHARACTER(LEN=32) :: ptr_addr
               INTEGER(8), DIMENSION(:), POINTER :: ptr, ptr2
               INTEGER, DIMENSION(:), POINTER :: int_buf
               INTEGER(HSIZE_T) :: tsize
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER :: error, rank, bsize, tclass
               TYPE(C_PTR) :: f_ptr
               CLASS(intrinsic_tree), POINTER :: ptr_tree

               CALL lookup_pointer(ptr_addr, pointer_data%itree, ptr_tree)
               IF (associated(ptr_tree)) THEN
                  SELECT TYPE(ptr_tree)
                  CLASS is (integer64_1d_tree)
                     ptr2 => ptr_tree%val
                  END SELECT
                  RETURN
               ENDIF

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

            END FUNCTION recover_integer64_1d_pointer
            FUNCTION recover_integer64_2d_pointer(cpfile_id, ptr_addr, &
                     ptr) RESULT(ptr2)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id, dset_id
               CHARACTER(LEN=32) :: ptr_addr
               INTEGER(8), DIMENSION(:,:), POINTER :: ptr, ptr2
               INTEGER, DIMENSION(:), POINTER :: int_buf
               INTEGER(HSIZE_T) :: tsize
               INTEGER(HSIZE_T), DIMENSION(2) :: dims
               INTEGER :: error, rank, bsize, tclass
               TYPE(C_PTR) :: f_ptr
               CLASS(intrinsic_tree), POINTER :: ptr_tree

               CALL lookup_pointer(ptr_addr, pointer_data%itree, ptr_tree)
               IF (associated(ptr_tree)) THEN
                  SELECT TYPE(ptr_tree)
                  CLASS is (integer64_2d_tree)
                     ptr2 => ptr_tree%val
                  END SELECT
                  RETURN
               ENDIF

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

               ptr2 => ptr
            END FUNCTION recover_integer64_2d_pointer
            FUNCTION recover_logical1d_pointer(cpfile_id, ptr_addr, &
                     ptr) RESULT(ptr2)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=32) :: ptr_addr
               LOGICAL, DIMENSION(:), POINTER :: ptr, ptr2
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER :: error, rank, bsize, tclass
               INTEGER, DIMENSION(:), POINTER :: int_buf
               INTEGER(HSIZE_T) :: tsize
               CLASS(intrinsic_tree), POINTER :: ptr_tree

               CALL lookup_pointer(ptr_addr, pointer_data%itree, ptr_tree)
               IF (associated(ptr_tree)) THEN
                  SELECT TYPE(ptr_tree)
                  CLASS is (logical1d_tree)
                     ptr2 => ptr_tree%val
                  END SELECT
                  RETURN
               ENDIF
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

               ptr2 => ptr
            END FUNCTION recover_logical1d_pointer
            FUNCTION recover_logical2d_pointer(cpfile_id, ptr_addr, &
                     ptr) RESULT(ptr2)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=32) :: ptr_addr
               LOGICAL, DIMENSION(:,:), POINTER :: ptr, ptr2
               INTEGER(HSIZE_T), DIMENSION(2) :: dims
               INTEGER :: error, rank, bsize, tclass
               INTEGER, DIMENSION(:), POINTER :: int_buf
               INTEGER(HSIZE_T) :: tsize
               CLASS(intrinsic_tree), POINTER :: ptr_tree

               CALL lookup_pointer(ptr_addr, pointer_data%itree, ptr_tree)
               IF (associated(ptr_tree)) THEN
                  SELECT TYPE(ptr_tree)
                  CLASS is (logical2d_tree)
                     ptr2 => ptr_tree%val
                  END SELECT
                  RETURN
               ENDIF
               rank = 2

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

               ptr2 => ptr
            END FUNCTION recover_logical2d_pointer

            FUNCTION recover_complex1d_pointer(cpfile_id, ptr_addr, &
                     ptr) RESULT(ptr2)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               CHARACTER(LEN=32) :: ptr_addr
               COMPLEX(8), DIMENSION(:), POINTER :: ptr, ptr2
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER :: error, rank, bsize, tclass, i
               INTEGER, DIMENSION(:), POINTER :: int_buf
               INTEGER(HSIZE_T) :: tsize, sizef, offset
               CLASS(intrinsic_tree), POINTER :: ptr_tree
               INTEGER(HID_T) :: group_id, dset_id, type_id, &
                  sub_id, array_id
               REAL(8), DIMENSION(:), ALLOCATABLE :: buffer

               CALL lookup_pointer(ptr_addr, pointer_data%itree, ptr_tree)
               IF (associated(ptr_tree)) THEN
                  SELECT TYPE(ptr_tree)
                  CLASS is (complex1d_tree)
                     ptr2 => ptr_tree%val
                  END SELECT
                  RETURN
               ENDIF
               rank = 1

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)

               CALL h5ltget_dataset_info_f(group_id, ptr_addr, dims, &
                  tclass, tsize, error)
               ALLOCATE (ptr(int(dims(1))))

               offset = 0
               CALL h5tcreate_f(H5T_COMPOUND_F, sizef*2, type_id, error)
               CALL h5tinsert_f(type_id, "real", offset, &
                  H5T_NATIVE_REAL, error)
               offset = offset + sizef
               CALL h5tinsert_f(type_id, "im", sizef, &
                  H5T_NATIVE_REAL, error)

               CALL h5tarray_create_f(type_id, rank, dims, array_id, &
                     error)

               ! store the real parts
               ALLOCATE(buffer(int(dims(1))))
               offset = 0
               CALL h5tcreate_f(H5T_COMPOUND_F, sizef, sub_id, error)
               CALL h5tinsert_f(sub_id, "real", sizef, H5T_NATIVE_REAL, error)
               CALL h5dread_f(dset_id, sub_id, buffer, dims, error)
               CALL h5tclose_f(sub_id, error)
               DO i=1, int(dims(1))
                 ptr(i) = cmplx(buffer(i), 0)
               ENDDO

               ! Now the imaginary parts
               offset = 0
               CALL h5tcreate_f(H5T_COMPOUND_F, sizef, sub_id, error)
               CALL h5tinsert_f(sub_id, "im", sizef, H5T_NATIVE_REAL, error)
               CALL h5dread_f(dset_id, sub_id, buffer, dims, error)
               CALL h5tclose_f(sub_id, error)
               DO i=1, int(dims(1))
                  ptr(i) = cmplx(REAL(ptr(i)), buffer(i))
               ENDDO
               DEALLOCATE (buffer)

               CALL h5tclose_f(sub_id, error)
               CALL h5tclose_f(type_id, error)

               CALL h5gclose_f(group_id, error)

               bsize = size(transfer(ptr, int_buf))
               ALLOCATE (int_buf(bsize))
               int_buf = transfer(ptr,int_buf)
               ALLOCATE(complex1d_tree::ptr_tree)
               SELECT TYPE (ptr_tree)
               TYPE is (complex1d_tree)
                  ptr_tree%val => ptr
                  ptr_tree%hash = FNVHash(int_buf, bsize)
                  ptr_tree%key = ptr_addr
               END SELECT
               DEALLOCATE(int_buf)

               CALL insert_intrinsic(pointer_data%itree, ptr_tree)

               ptr2 => ptr
            END FUNCTION recover_complex1d_pointer
            FUNCTION recover_complex2d_pointer(cpfile_id, ptr_addr, &
                     ptr) RESULT(ptr2)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               CHARACTER(LEN=32) :: ptr_addr
               COMPLEX(8), DIMENSION(:,:), POINTER :: ptr, ptr2
               INTEGER(HSIZE_T), DIMENSION(2) :: dims
               INTEGER :: error, rank, bsize, tclass, i ,j
               INTEGER, DIMENSION(:), POINTER :: int_buf
               INTEGER(HSIZE_T) :: tsize, sizef, offset
               CLASS(intrinsic_tree), POINTER :: ptr_tree
               INTEGER(HID_T) :: group_id, dset_id, type_id, &
                  sub_id, array_id
               REAL(8), DIMENSION(:,:), ALLOCATABLE :: buffer

               CALL lookup_pointer(ptr_addr, pointer_data%itree, ptr_tree)
               IF (associated(ptr_tree)) THEN
                  SELECT TYPE(ptr_tree)
                  CLASS is (complex2d_tree)
                     ptr2 => ptr_tree%val
                  END SELECT
                  RETURN
               ENDIF
               rank = 2

               CALL h5gopen_f(cpfile_id, "intrinsic", group_id, error)

               CALL h5ltget_dataset_info_f(group_id, ptr_addr, dims, &
                  tclass, tsize, error)
               ALLOCATE (ptr(int(dims(1)),int(dims(2))))

               offset = 0
               CALL h5tcreate_f(H5T_COMPOUND_F, sizef*2, type_id, error)
               CALL h5tinsert_f(type_id, "real", offset, &
                  H5T_NATIVE_REAL, error)
               offset = offset + sizef
               CALL h5tinsert_f(type_id, "im", sizef, &
                  H5T_NATIVE_REAL, error)

               CALL h5tarray_create_f(type_id, rank, dims, array_id, &
                     error)

               ! store the real parts
               ALLOCATE(buffer(int(dims(1)), int(dims(2))))
               offset = 0
               CALL h5tcreate_f(H5T_COMPOUND_F, sizef, sub_id, error)
               CALL h5tinsert_f(sub_id, "real", sizef, H5T_NATIVE_REAL, error)
               CALL h5dread_f(dset_id, sub_id, buffer, dims, error)
               CALL h5tclose_f(sub_id, error)
               DO i=1, int(dims(1))
                 DO j=1, int(dims(2))
                    ptr(i,j) = cmplx(buffer(i,j), 0)
                 ENDDO
               ENDDO

               ! Now the imaginary parts
               offset = 0
               CALL h5tcreate_f(H5T_COMPOUND_F, sizef, sub_id, error)
               CALL h5tinsert_f(sub_id, "im", sizef, H5T_NATIVE_REAL, error)
               CALL h5dread_f(dset_id, sub_id, buffer, dims, error)
               CALL h5tclose_f(sub_id, error)
               DO i=1, int(dims(1))
                  DO j=1, int(dims(2))
                     ptr(i,j) = cmplx(REAL(ptr(i,j)), buffer(i,j))
                  ENDDO
               ENDDO
               DEALLOCATE (buffer)

               CALL h5tclose_f(sub_id, error)
               CALL h5tclose_f(type_id, error)

               CALL h5gclose_f(group_id, error)

               bsize = size(transfer(ptr, int_buf))
               ALLOCATE (int_buf(bsize))
               int_buf = transfer(ptr,int_buf)
               ALLOCATE(complex2d_tree::ptr_tree)
               SELECT TYPE (ptr_tree)
               TYPE is (complex2d_tree)
                  ptr_tree%val => ptr
                  ptr_tree%hash = FNVHash(int_buf, bsize)
                  ptr_tree%key = ptr_addr
               END SELECT
               DEALLOCATE(int_buf)

               CALL insert_intrinsic(pointer_data%itree, ptr_tree)

               ptr2 => ptr
            END FUNCTION recover_complex2d_pointer
