            SUBROUTINE READ_logical_array_2d(group_id, dname, buffer,&
                   dims)
               INTEGER(HID_T), INTENT(IN) :: group_id
               CHARACTER(LEN=*), INTENT(IN) :: dname
               LOGICAL, DIMENSION(:,:) :: buffer
               INTEGER(HSIZE_T), DIMENSION(2) :: dims
               CHARACTER, DIMENSION(dims(1), dims(2)) :: charbuf

               INTEGER(HID_T) :: space_id, dset_id
               INTEGER :: error, i, v, rank
               rank = 2

               CALL h5screate_simple_f(rank, dims, space_id, error)
               CALL h5dcreate_f(group_id, dname, H5T_NATIVE_CHARACTER, &
                  space_id, dset_id, error)
               CALL h5dread_f(dset_id, H5T_NATIVE_CHARACTER, charbuf, &
                  dims, error)
               CALL h5dclose_f(dset_id, error)
               CALL h5sclose_f(space_id, error)

               DO i=1, int(dims(1))
                  DO v=1, int(dims(2))
                     IF (charbuf(i,v) == 'T') THEN
                        buffer(i,v) = .TRUE.
                     ELSE
                        buffer(i,v) = .FALSE.
                     ENDIF
                  ENDDO
               ENDDO
            END SUBROUTINE READ_logical_array_2d

            ! Write a BOOL array attribute to a dataset
            SUBROUTINE READ_logical_array(dset_id, &
                  dname, buffer, length)
               INTEGER(HID_T), INTENT(IN) :: dset_id
               CHARACTER(LEN=*), INTENT(IN) :: dname
               LOGICAL, DIMENSION(:) :: buffer
               INTEGER, INTENT(in) :: length
               CHARACTER, DIMENSION(length) :: charbuf
               !ALLOCATE(charbuf(length))

               INTEGER :: i

               CALL read_character_array(dset_id, dname, &
                  charbuf, length)

               DO i=1, length
                  IF (charbuf(i) == 'T') THEN
                     buffer(i) = .TRUE.
                  ELSE
                     buffer(i) = .FALSE.
                  ENDIF
               ENDDO

            END SUBROUTINE read_logical_array

            SUBROUTINE read_character_array(dset_id, dname, &
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

               CALL h5dread_f(dset_id, attr_id , charbuf, &
                  dims, error)
               CALL h5tclose_f(attr_id, error)
            END SUBROUTINE read_character_array

            SUBROUTINE read_string_attribute(dset_id, dname, charbuf, &
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
               CALL h5tinsert_f(attr_id, dname, offset, bool_id, error)

               CALL h5dread_f(dset_id, attr_id , charbuf, &
                  dims, error)
               CALL h5tclose_f(attr_id, error)
            END SUBROUTINE read_string_attribute

            ! Write integer attribute to a dataset
            SUBROUTINE READ_integer_attribute(dset_id, &
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
               CALL h5dREAD_f(dset_id, attr_id, val, dims, &
                     error, xfer_prp = plist_id)

               CALL h5tclose_f(attr_id, error)

            END SUBROUTINE READ_integer_attribute

            ! Write double attribute to a dataset
            SUBROUTINE read_double_attribute(dset_id, vname, &
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
               CALL h5dread_f(dset_id, attr_id, val, dims, &
                     error, xfer_prp = plist_id)

               CALL h5tclose_f(attr_id, error)

            END SUBROUTINE read_double_attribute

            SUBROUTINE read_logical_attribute(dset_id, &
                  vname, val)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: dset_id
               CHARACTER(LEN=*), INTENT(IN) :: vname
               LOGICAL :: val
               CHARACTER :: c
               CALL read_character_attribute(dset_id, vname, c)
               IF (c == 'T') THEN
                  val = .TRUE.
               ELSE
                  val = .FALSE.
               ENDIF
            END SUBROUTINE read_logical_attribute

            SUBROUTINE read_character_attribute(dset_id, &
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
               CALL h5dread_f(dset_id, attr_id, val, dims, &
                     error, xfer_prp = plist_id)

               CALL h5tclose_f(attr_id, error)

            END SUBROUTINE read_character_attribute
