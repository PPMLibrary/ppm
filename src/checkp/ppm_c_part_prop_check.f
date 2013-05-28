            SUBROUTINE store_ppm_c_part_prop_d_(cpfile_id, cprop_id, &
                         cprop)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=*), INTENT(IN) :: cprop_id
               CLASS(ppm_c_part_prop_d_), POINTER :: cprop
               INTEGER error

               CALL h5gcreate_f(cpfile_id, 'props/'//cprop_id, &
                  group_id, error)


               CALL h5gclose_f(group_id, error)
            END SUBROUTINE store_ppm_c_part_prop_d_
