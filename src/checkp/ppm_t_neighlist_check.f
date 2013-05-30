      ! Store routines for ppm_t_neighlist_d_
            SUBROUTINE make_type_ppm_t_neighlist_d_(dtype_id)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(OUT) :: dtype_id
               INTEGER(HID_T) :: array_id
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER rank
               INTEGER error
               INTEGER(HSIZE_T) :: tsize, offset, subsize
               INTEGER(HSIZE_T) :: isize
               INTEGER(HSIZE_T) :: dsize
               INTEGER(HSIZE_T) :: csize

               offset = 0
               tsize = 0


               ! Calculate datatype size
               CALL h5tget_size_f(dtype_id, tsize, error) ! initial size

               ! Calculate datatype size
               CALL h5tget_size_f(H5T_NATIVE_INTEGER, isize, error)
               CALL h5tget_size_f(H5T_NATIVE_DOUBLE, dsize, error)
               CALL h5tget_size_f(H5T_NATIVE_CHARACTER, csize, error)


               tsize = tsize + isize*(3) + dsize*(2) + csize*(ppm_char+0)



               ! Create/Expand the datatype
               CALL h5tset_size_f(dtype_id, tsize, error)
               CALL h5tcreate_f(H5T_COMPOUND_F, tsize, dtype_id, error)

               ! Insert the members
               ! Insert the members
               ! Integer members
               CALL h5tinsert_f(dtype_id, "nneighmin", offset, &
                     H5T_NATIVE_INTEGER, error)
               offset = offset + isize
               CALL h5tinsert_f(dtype_id, "isymm", offset, &
                     H5T_NATIVE_INTEGER, error)
               offset = offset + isize
               CALL h5tinsert_f(dtype_id, "nneighmax", offset, &
                     H5T_NATIVE_INTEGER, error)
               offset = offset + isize

               ! Real members
               CALL h5tinsert_f(dtype_id, "skin", offset, &
                     H5T_NATIVE_DOUBLE, error)
               offset = offset + dsize
               CALL h5tinsert_f(dtype_id, "cutoff", offset, &
                     H5T_NATIVE_DOUBLE, error)
               offset = offset + dsize

               ! Character members
               dims = (/csize*ppm_char/)
               CALL h5tarray_create_f(H5T_NATIVE_CHARACTER, rank, &
                   dims, array_id, error)
               CALL h5tinsert_f(dtype_id, "name", offset, &
                   array_id, error)
               offset = offset + (csize*ppm_char)

               ! Logical members
               CALL h5tinsert_f(dtype_id, "uptodate", offset, &
                     H5T_NATIVE_CHARACTER, error)
               offset = offset + csize



            END SUBROUTINE make_type_ppm_t_neighlist_d_

            SUBROUTINE store_ppm_t_neighlist_d_(cpfile_id, &
                  type_ptr_id, type_ptr)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id, type_id, dset_id, &
                  dspace_id
               CHARACTER(LEN=*), INTENT(IN) :: type_ptr_id
               CLASS(ppm_t_neighlist_d_), POINTER :: type_ptr
               INTEGER error

               CALL h5gopen_f(cpfile_id, 'ppm_t_neighlist_d_', &
                  group_id, error)


               ! Make our dataset
               CALL make_type_ppm_t_neighlist_d_(type_id) ! get type
               CALL h5screate_f(H5S_SCALAR_F, dspace_id, error)  ! get space
               CALL h5dcreate_f(group_id, type_ptr_id, type_id, &
                   dspace_id, dset_id, error)

               CALL write_ppm_t_neighlist_d_(dset_id, type_ptr)

               CALL h5dclose_f(dset_id, error)
               CALL h5sclose_f(dspace_id, error)

               CALL h5gclose_f(group_id, error)
            END SUBROUTINE store_ppm_t_neighlist_d_

            SUBROUTINE write_ppm_t_neighlist_d_(dset_id, type_ptr)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: dset_id
               CLASS(ppm_t_neighlist_d_), POINTER :: type_ptr

               !WRITE_STUB

            END SUBROUTINE write_ppm_t_neighlist_d_
