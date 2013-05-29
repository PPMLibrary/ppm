            ! This type def needs to be expanded
            SUBROUTINE make_ppm_t_mapping_d_type(dtype_id)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(OUT) :: dtype_id
               INTEGER error
               INTEGER(HSIZE_T) :: typei, type_size, offset
               offset = 0

               CALL h5tget_size_f(H5T_NATIVE_INTEGER, typei, error)

               ! Just two integers for now
               type_size = typei*2

               ! Create the data type
               CALL h5tcreate_f(H5T_COMPOUND_F, type_size, dtype_id, &
                     error)

               CALL h5tinsert_f(dtype_id, "source_topoid", offset, &
                  H5T_NATIVE_INTEGER, error)
               offset = offset + typei
               CALL h5tinsert_f(dtype_id, "target_topoid", offset, &
                  H5T_NATIVE_INTEGER, error)
            END SUBROUTINE make_ppm_t_mapping_d_type
            SUBROUTINE write_ppm_t_mapping_d_(dset_id, map)
               INTEGER(HID_T), INTENT(IN) :: dset_id
               CLASS(ppm_t_mapping_d_), POINTER, INTENT(IN) :: map

               INTEGER (HID_T) :: type_id, space_id
               INTEGER error

               CALL write_attribute(dset_id, 'source_topoid', &
                  map%source_topoid)
               CALL write_attribute(dset_id, 'target_topoid', &
                  map%TARGET_topoid)

            END SUBROUTINE


