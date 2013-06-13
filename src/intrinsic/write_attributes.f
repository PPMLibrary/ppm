            ! This function will be removed for the genericized version
            SUBROUTINE store_logical_dim(group_id, dname, buffer,&
                   length)
               INTEGER(HID_T), INTENT(IN) :: group_id
               CHARACTER(LEN=*), INTENT(IN) :: dname
               LOGICAL, DIMENSION(length) :: buffer
               INTEGER, INTENT(in) :: length
               CHARACTER, DIMENSION(length) :: charbuf

               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER(HID_T) :: space_id, dset_id
               INTEGER :: error, i, rank
               rank = 1
               dims = (/length/)

               DO i=1, length
                  IF (buffer(i)) THEN
                     charbuf(i) = 'T'
                  ELSE
                     charbuf(i) = 'F'
                  ENDIF
               ENDDO
               CALL h5screate_simple_f(rank, dims, space_id, error)
               CALL h5dcreate_f(group_id, dname, H5T_NATIVE_CHARACTER, &
                  space_id, dset_id, error)
               CALL h5dwrite_f(dset_id, H5T_NATIVE_CHARACTER, charbuf, &
                  dims, error)
               CALL h5dclose_f(dset_id, error)
               CALL h5sclose_f(space_id, error)
            END SUBROUTINE store_logical_dim

            ! Write a BOOL array attribute to a dataset
            SUBROUTINE write_logical_array(dset_id, &
                  dname, buffer, length)
               INTEGER(HID_T), INTENT(IN) :: dset_id
               CHARACTER(LEN=*), INTENT(IN) :: dname
               LOGICAL, DIMENSION(:) :: buffer
               INTEGER, INTENT(in) :: length
               CHARACTER, DIMENSION(length) :: charbuf
               !ALLOCATE(charbuf(length))

               INTEGER :: i

               DO i=1, length
                  IF (buffer(i)) THEN
                     charbuf(i) = 'T'
                  ELSE
                     charbuf(i) = 'F'
                  ENDIF
               ENDDO
               CALL write_character_array(dset_id, dname, &
                  charbuf, length)

            END SUBROUTINE write_logical_array

            SUBROUTINE write_character_array(dset_id, dname, &
                  charbuf, length)
               INTEGER(HID_T), INTENT(IN) :: dset_id
               CHARACTER(LEN=*), INTENT(IN) :: dname
               CHARACTER, DIMENSION(length) :: charbuf
               INTEGER, INTENT(in) :: length

               INTEGER(HID_T) :: attr_id, plist_id, bool_id
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER(HSIZE_T) :: csize, arr_size, offset
               INTEGER :: error, rank
               rank = 1
               dims = (/length/)
               offset = 0

               ! Create the type properties
               ! preserve partially initialized fields
               CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
               CALL h5pset_preserve_f(plist_id, .TRUE., error)

               ! Get the data size
               CALL h5tget_size_f(H5T_NATIVE_CHARACTER, csize, error)
               arr_size = length*csize

               ! Create the array type for the attribute
               CALL h5tarray_create_f(H5T_NATIVE_CHARACTER, rank, dims,&
                  bool_id, error)
               CALL h5tget_size_f(bool_id, arr_size, error)
               CALL h5tcreate_f(H5T_COMPOUND_F, arr_size, attr_id, &
                  error)
               CALL h5tinsert_f(attr_id, dname, offset, bool_id, error)

               CALL h5dwrite_f(dset_id, attr_id , charbuf, &
                  dims, error)
               CALL h5tclose_f(attr_id, error)
            END SUBROUTINE write_character_array

            SUBROUTINE write_string_attribute(dset_id, dname, charbuf, &
               length)
               INTEGER(HID_T), INTENT(IN) :: dset_id
               CHARACTER(LEN=*), INTENT(IN) :: dname
               CHARACTER(LEN=length) :: charbuf
               INTEGER, INTENT(in) :: length

               INTEGER(HID_T) :: attr_id, plist_id, bool_id
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER(HSIZE_T) :: csize, arr_size, offset
               INTEGER :: error, rank
               rank = 1
               dims = (/0/)
               offset = 0

               ! Create the type properties
               ! preserve partially initialized fields
               CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
               CALL h5pset_preserve_f(plist_id, .TRUE., error)

               ! Get the data size
               CALL h5tget_size_f(H5T_NATIVE_CHARACTER, csize, error)
               arr_size = length*csize

               ! Create the array type for the attribute
               CALL h5tcreate_f(H5T_STRING_F, arr_size, &
                  bool_id, error)
               CALL h5tget_size_f(bool_id, arr_size, error)
               CALL h5tcreate_f(H5T_COMPOUND_F, arr_size, attr_id, &
                  error)
               WRITE(*,*) "insert create error = " , error
               CALL h5tinsert_f(attr_id, dname, offset, bool_id, error)
               WRITE(*,*) "insert type error = " , error

               CALL h5dwrite_f(dset_id, attr_id , charbuf, &
                  dims, error)
               WRITE(*,*) "write string error = " , error
               CALL h5tclose_f(attr_id, error)
            END SUBROUTINE write_string_attribute

            ! Write integer attribute to a dataset
            SUBROUTINE write_integer_attribute(dset_id, &
                  vname, val)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: dset_id
               CHARACTER(LEN=*), INTENT(IN) :: vname
               INTEGER :: val

               INTEGER error
               INTEGER(HID_T) :: plist_id, attr_id
               INTEGER(HSIZE_T) :: isize, offset
               INTEGER(HSIZE_T), DIMENSION(1) :: dims = (/0/)
               offset = 0

               ! Create the type properties
               ! preserve partially initialized fields
               CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
               CALL h5pset_preserve_f(plist_id, .TRUE., error)

               ! Get the data size
               CALL h5tget_size_f(H5T_NATIVE_INTEGER, isize, error)

               ! Create the compound type for the attribute
               CALL h5tcreate_f(H5T_COMPOUND_F, isize, attr_id, error)
               CALL h5tinsert_f(attr_id, vname, offset, &
                   H5T_NATIVE_INTEGER, error)
               CALL h5dwrite_f(dset_id, attr_id, val, dims, &
                     error, xfer_prp = plist_id)

               CALL h5tclose_f(attr_id, error)

            END SUBROUTINE write_integer_attribute

            ! Write double attribute to a dataset
            SUBROUTINE write_double_attribute(dset_id, vname, &
                  val)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: dset_id
               CHARACTER(LEN=*), INTENT(IN) :: vname
               REAL(ppm_kind_double) :: val

               INTEGER error
               INTEGER(HID_T) :: plist_id, attr_id
               INTEGER(HSIZE_T) :: dsize, offset
               INTEGER(HSIZE_T), DIMENSION(1) :: dims = (/0/)
               offset = 0

               ! Create the type properties
               ! preserve partially initialized fields
               CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
               CALL h5pset_preserve_f(plist_id, .TRUE., error)

               ! Get the data size
               CALL h5tget_size_f(H5T_NATIVE_DOUBLE, dsize, error)

               ! Create the compound type for the attribute
               CALL h5tcreate_f(H5T_COMPOUND_F, dsize, attr_id, error)
               CALL h5tinsert_f(attr_id, vname, offset, &
                   H5T_NATIVE_DOUBLE, error)
               CALL h5dwrite_f(dset_id, attr_id, val, dims, &
                     error, xfer_prp = plist_id)

               CALL h5tclose_f(attr_id, error)

            END SUBROUTINE write_double_attribute

            SUBROUTINE write_logical_attribute(dset_id, &
                  vname, val)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: dset_id
               CHARACTER(LEN=*), INTENT(IN) :: vname
               LOGICAL :: val
               CHARACTER :: c
               IF (val) THEN
                  c = 'T'
               ELSE
                  c = 'F'
               ENDIF
               CALL write_character_attribute(dset_id, vname, c)
            END SUBROUTINE write_logical_attribute

            SUBROUTINE write_character_attribute(dset_id, &
                  vname, val)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: dset_id
               CHARACTER(LEN=*), INTENT(IN) :: vname
               CHARACTER :: val

               INTEGER error
               INTEGER(HID_T) :: plist_id, attr_id
               INTEGER(HSIZE_T) :: csize, offset
               INTEGER(HSIZE_T), DIMENSION(1) :: dims = (/0/)
               offset = 0

               ! Create the type properties
               ! preserve partially initialized fields
               CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
               CALL h5pset_preserve_f(plist_id, .TRUE., error)

               ! Get the data size
               CALL h5tget_size_f(H5T_NATIVE_CHARACTER, csize, error)

               ! Create the compound type for the attribute
               CALL h5tcreate_f(H5T_COMPOUND_F, csize, attr_id, error)
               CALL h5tinsert_f(attr_id, vname, offset, &
                   H5T_NATIVE_CHARACTER, error)
               CALL h5dwrite_f(dset_id, attr_id, val, dims, &
                     error, xfer_prp = plist_id)

               CALL h5tclose_f(attr_id, error)

            END SUBROUTINE write_character_attribute
