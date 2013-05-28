            SUBROUTINE store_ppm_v_main_abstr(cpfile_id, abstr_id,&
                         abstr)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=*), INTENT(IN) :: abstr_id
               CLASS(ppm_v_main_abstr), POINTER :: abstr
               INTEGER error

               CALL h5gcreate_f(cpfile_id, 'field_ptr/'//abstr_id, &
                      group_id, error)


               CALL h5gclose_f(group_id, error)
            END SUBROUTINE store_ppm_v_main_abstr
