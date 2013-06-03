      ! Store routines for ppm_t_mapping_d_
            SUBROUTINE make_type_ppm_t_mapping_d_(dtype_id)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(OUT) :: dtype_id
               INTEGER error
               INTEGER(HSIZE_T) :: tsize, offset
               INTEGER(HSIZE_T) :: isize

               offset = 0
               tsize = 0


               ! Calculate datatype size

               CALL h5tget_size_f(H5T_NATIVE_INTEGER, isize, error)


               tsize = tsize + isize*(12)



               ! Create/Expand the datatype
               CALL h5tcreate_f(H5T_COMPOUND_F, tsize, dtype_id, error)

               ! Insert the members
               ! Insert the members
               ! Integer members
               CALL h5tinsert_f(dtype_id, "source_topoid", offset, &
                     H5T_NATIVE_INTEGER, error)
               offset = offset + isize
               CALL h5tinsert_f(dtype_id, "ppm_nrecvlist", offset, &
                     H5T_NATIVE_INTEGER, error)
               offset = offset + isize
               CALL h5tinsert_f(dtype_id, "ppm_nsendbuffer", offset, &
                     H5T_NATIVE_INTEGER, error)
               offset = offset + isize
               CALL h5tinsert_f(dtype_id, "old_buffer_set", offset, &
                     H5T_NATIVE_INTEGER, error)
               offset = offset + isize
               CALL h5tinsert_f(dtype_id, "old_nsendlist", offset, &
                     H5T_NATIVE_INTEGER, error)
               offset = offset + isize
               CALL h5tinsert_f(dtype_id, "ppm_map_type", offset, &
                     H5T_NATIVE_INTEGER, error)
               offset = offset + isize
               CALL h5tinsert_f(dtype_id, "target_topoid", offset, &
                     H5T_NATIVE_INTEGER, error)
               offset = offset + isize
               CALL h5tinsert_f(dtype_id, "ppm_buffer_set", offset, &
                     H5T_NATIVE_INTEGER, error)
               offset = offset + isize
               CALL h5tinsert_f(dtype_id, "ppm_nrecvbuffer", offset, &
                     H5T_NATIVE_INTEGER, error)
               offset = offset + isize
               CALL h5tinsert_f(dtype_id, "ppm_recvbufsize", offset, &
                     H5T_NATIVE_INTEGER, error)
               offset = offset + isize
               CALL h5tinsert_f(dtype_id, "ppm_nsendlist", offset, &
                     H5T_NATIVE_INTEGER, error)
               offset = offset + isize
               CALL h5tinsert_f(dtype_id, "ppm_sendbufsize", offset, &
                     H5T_NATIVE_INTEGER, error)
               offset = offset + isize

               ! Real members

               ! Character members

               ! Logical members



            END SUBROUTINE make_type_ppm_t_mapping_d_

            SUBROUTINE store_ppm_t_mapping_d_(cpfile_id, &
                  type_ptr_id, type_ptr)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id, type_id, dset_id, &
                  dspace_id
               CHARACTER(LEN=*), INTENT(IN) :: type_ptr_id
               CLASS(ppm_t_mapping_d_), POINTER :: type_ptr
               INTEGER error

               CALL h5gopen_f(cpfile_id, 'ppm_t_mapping_d_', &
                  group_id, error)


               ! Make our dataset
               CALL make_type_ppm_t_mapping_d_(type_id) ! get type
               CALL h5screate_f(H5S_SCALAR_F, dspace_id, error)  ! get space
               CALL h5dcreate_f(group_id, type_ptr_id, type_id, &
                   dspace_id, dset_id, error)

               CALL write_ppm_t_mapping_d_(dset_id, type_ptr)

               CALL h5dclose_f(dset_id, error)
               CALL h5sclose_f(dspace_id, error)

               CALL h5gclose_f(group_id, error)
            END SUBROUTINE store_ppm_t_mapping_d_

            SUBROUTINE write_ppm_t_mapping_d_(dset_id, type_ptr)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: dset_id
               CLASS(ppm_t_mapping_d_), POINTER :: type_ptr

               CALL write_attribute(dset_id, 'source_topoid', &
                   type_ptr%source_topoid)
               CALL write_attribute(dset_id, 'ppm_nrecvlist', &
                   type_ptr%ppm_nrecvlist)
               CALL write_attribute(dset_id, 'ppm_nsendbuffer', &
                   type_ptr%ppm_nsendbuffer)
               CALL write_attribute(dset_id, 'old_buffer_set', &
                   type_ptr%old_buffer_set)
               CALL write_attribute(dset_id, 'old_nsendlist', &
                   type_ptr%old_nsendlist)
               CALL write_attribute(dset_id, 'ppm_map_type', &
                   type_ptr%ppm_map_type)
               CALL write_attribute(dset_id, 'target_topoid', &
                   type_ptr%target_topoid)
               CALL write_attribute(dset_id, 'ppm_buffer_set', &
                   type_ptr%ppm_buffer_set)
               CALL write_attribute(dset_id, 'ppm_nrecvbuffer', &
                   type_ptr%ppm_nrecvbuffer)
               CALL write_attribute(dset_id, 'ppm_recvbufsize', &
                   type_ptr%ppm_recvbufsize)
               CALL write_attribute(dset_id, 'ppm_nsendlist', &
                   type_ptr%ppm_nsendlist)
               CALL write_attribute(dset_id, 'ppm_sendbufsize', &
                   type_ptr%ppm_sendbufsize)


            END SUBROUTINE write_ppm_t_mapping_d_
