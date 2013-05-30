            ! This type def needs to be expanded
            SUBROUTINE make_type_ppm_t_mapping_d_(dtype_id)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(OUT) :: dtype_id
               INTEGER error
               INTEGER(HSIZE_T) :: typei, type_size, offset, subsize
               offset = 0
               type_size = 0

               CALL h5tget_size_f(H5T_NATIVE_INTEGER, typei, error)
               type_size = type_size + typei*10

               ! Create the data type
               CALL h5tcreate_f(H5T_COMPOUND_F, type_size, dtype_id, &
                     error)

               ! Insert integer members
               CALL h5tinsert_f(dtype_id, "source_topoid", offset, &
                  H5T_NATIVE_INTEGER, error)
               offset = offset + typei
               CALL h5tinsert_f(dtype_id, "target_topoid", offset, &
                  H5T_NATIVE_INTEGER, error)
               offset = offset + typei
               CALL h5tinsert_f(dtype_id, "ppm_nsendbuffer", offset, &
                  H5T_NATIVE_INTEGER, error)
               offset = offset + typei
               CALL h5tinsert_f(dtype_id, "ppm_nrecvbuffer", offset, &
                  H5T_NATIVE_INTEGER, error)
               offset = offset + typei
               CALL h5tinsert_f(dtype_id, "ppm_sendbufsize", offset, &
                  H5T_NATIVE_INTEGER, error)
               offset = offset + typei
               CALL h5tinsert_f(dtype_id, "ppm_recvbuffsize", offset, &
                  H5T_NATIVE_INTEGER, error)
               offset = offset + typei
               CALL h5tinsert_f(dtype_id, "ppm_buffer_set", offset, &
                  H5T_NATIVE_INTEGER, error)
               offset = offset + typei
               CALL h5tinsert_f(dtype_id, "ppm_map_type", offset, &
                  H5T_NATIVE_INTEGER, error)
               offset = offset + typei
               CALL h5tinsert_f(dtype_id, "ppm_nsendlist", offset, &
                  H5T_NATIVE_INTEGER, error)
               offset = offset + typei
               CALL h5tinsert_f(dtype_id, "ppm_nrecvlist", offset, &
                  H5T_NATIVE_INTEGER, error)
               offset = offset + typei

            END SUBROUTINE make_type_ppm_t_mapping_d_

            ! Add store subroutine
            ! This is possibly unecessary
            SUBROUTINE store_ppm_t_mapping_d_(cpfile_id, mapping_id, &
                  mapping)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               CHARACTER(LEN=*), INTENT(in) :: mapping_id
               CLASS(ppm_t_mapping_d_), POINTER, INTENT(IN) :: mapping

               INTEGER(HID_T) :: group_id, dset_id, space_id, dtype_id
               INTEGER :: error

               CALL h5gopen_f(cpfile_id, 'ppm_t_mapping', group_id, &
                  error)

               CALL make_type_ppm_t_mapping_d_(dtype_id)
               CALL h5screate_f(H5S_SCALAR_F, space_id, error)
               CALL h5dcreate_f(group_id, mapping_id, dtype_id, &
                  space_id, dset_id, error)

               CALL write_ppm_t_mapping_d_(dset_id, mapping)

               CALL h5dclose_f(dset_id, error)
               CALL h5sclose_f(space_id, error)

               CALL h5gclose_f(group_id, error)

            END SUBROUTINE store_ppm_t_mapping_d_

            SUBROUTINE write_ppm_t_mapping_d_(dset_id, map)
               INTEGER(HID_T), INTENT(IN) :: dset_id
               CLASS(ppm_t_mapping_d_), POINTER, INTENT(IN) :: map

               CALL write_attribute(dset_id, 'source_topoid', &
                  map%source_topoid)
               CALL write_attribute(dset_id, 'target_topoid', &
                  map%target_topoid)
               CALL write_attribute(dset_id, 'ppm_nsendbuffer', &
                  map%ppm_nsendbuffer)
               CALL write_attribute(dset_id, 'ppm_nrecvbuffer', &
                  map%ppm_nrecvbuffer)
               CALL write_attribute(dset_id, 'ppm_sendbufsize', &
                  map%ppm_sendbufsize)
               CALL write_attribute(dset_id, 'ppm_recvbufsize', &
                  map%ppm_recvbufsize)
               CALL write_attribute(dset_id, 'ppm_buffer_set', &
                  map%ppm_buffer_set)
               CALL write_attribute(dset_id, 'ppm_map_type', &
                  map%ppm_map_type)
               CALL write_attribute(dset_id, 'ppm_nsendlist', &
                  map%ppm_nsendlist)
               CALL write_attribute(dset_id, 'ppm_nrecvlist', &
                  map%ppm_nrecvlist)

            END SUBROUTINE write_ppm_t_mapping_d_


