            SUBROUTINE store_ppm_c_part_mapping_d(cpfile_id,&
                         mapping_id, mapping)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=*), INTENT(IN) :: mapping_id
               CLASS(ppm_c_part_mapping_d_), POINTER :: mapping
               INTEGER error

               CALL h5gcreate_f(cpfile_id, 'maps/'//mapping_id, &
                     group_id, error)


               CALL h5gclose_f(group_id, error)
            END SUBROUTINE store_ppm_c_part_mapping_d
