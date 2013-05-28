            SUBROUTINE store_ppm_c_operator_discr_(cpfile_id, discr_id,&
                         operator_discr)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=*), INTENT(IN) :: discr_id
               CLASS(ppm_c_operator_discr_), POINTER :: discr
               INTEGER error

               CALL h5gcreate_f(cpfile_id, 'ops/'//discr_id, group_id,&
                      error)


               CALL h5gclose_f(group_id, error)
            END SUBROUTINE store_ppm_c_operator_discr_
