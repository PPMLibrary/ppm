            SUBROUTINE make_type_ppm_t_part_mapping_d_(dtype_id)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(OUT) :: dtype_id
               INTEGER(HID_T) :: parent_id
               INTEGER error 
               INTEGER(HSIZE_T) ::  typei, tsize, offset

               ! Start with the parent Type
               CALL make_type_ppm_t_mapping_d_(parent_id)
               CALL h5tcopy_f(parent_id, dtype_id, error)
               CALL h5tget_size_f(dtype_id, tsize, error) ! initial size
               offset = tsize ! initial offset

               ! Calculate datatype size
               CALL h5tget_size_f(H5T_NATIVE_INTEGER, typei, error)
               tsize = tsize + typei * 2

               ! Expand the datatype
               CALL h5tset_size_f(dtype_id, tsize, error)

               ! Insert the members
               CALL h5tinsert_f(dtype_id, "oldNpart", offset, &
                  H5T_NATIVE_INTEGER, error)
               offset = offset + typei

               CALL h5tinsert_f(dtype_id, "newNpart", offset, &
                  H5T_NATIVE_INTEGER, error)
               offset = offset + typei

            END SUBROUTINE make_type_ppm_t_part_mapping_d_

            SUBROUTINE store_ppm_t_part_mapping_d_(cpfile_id, &
                  mapping_id, mapping)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id, type_id, dset_id, &
                  dspace_id
               CHARACTER(LEN=*), INTENT(IN) :: mapping_id
               CLASS(ppm_t_part_mapping_d_), POINTER :: mapping
               INTEGER error

               CALL h5gopen_f(cpfile_id, 'ppm_t_part_mapping', &
                  group_id, error)


               ! Make our dataset
               CALL make_type_ppm_t_part_mapping_d_(type_id) ! get type
               CALL h5screate_f(H5S_SCALAR_F, dspace_id, error)  ! get space
               CALL h5dcreate_f(group_id, mapping_id, type_id, &
                   dspace_id, dset_id, error)

               CALL write_ppm_t_part_mapping_d_(dset_id, mapping)

               CALL h5dclose_f(dset_id, error)
               CALL h5sclose_f(dspace_id, error)

               CALL h5gclose_f(group_id, error)
            END SUBROUTINE store_ppm_t_part_mapping_d_

            SUBROUTINE write_ppm_t_part_mapping_d_(dset_id, mapping)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: dset_id
               CLASS(ppm_t_part_mapping_d_), POINTER :: mapping
               CLASS(ppm_t_mapping_d_), POINTER :: mparent
               CALL write_attribute(dset_id, 'oldNpart', &
                  mapping%oldNpart)
               CALL write_attribute(dset_id, 'newPpart', &
                  mapping%newNpart)

               mparent => mapping
               CALL write_ppm_t_mapping_d_(dset_id, &
                      mparent)
            END SUBROUTINE write_ppm_t_part_mapping_d_
