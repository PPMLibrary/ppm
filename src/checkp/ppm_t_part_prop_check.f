            SUBROUTINE store_ppm_t_part_prop_d_(cpfile_id, tprop_id, &
                  tprop)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               CHARACTER(LEN=20), INTENT(IN) :: tprop_id
               CLASS(ppm_t_part_prop_d_), POINTER :: tprop
               INTEGER(HID_T) :: group_id
               INTEGER :: error

               CALL h5gcreate_f(cpfile_id, 'gi/'//tprop_id, group_id, &
                  error)

               ! TODO
               ! All pointers will work on later
               ! Integer dimensions
               ! Real Dimensions
               ! Complex Dimensions
               ! Logical Dimensions

               CALL h5gclose_f(group_id, error)
            END SUBROUTINE store_ppm_t_part_prop_d_
