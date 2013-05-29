         SUBROUTINE store_ppm_t_operator_discr_(cpfile_id, discr_id,&
                   discr)
            INTEGER(HID_T), INTENT(in) :: cpfile_id
            CHARACTER(LEN=*), INTENT(IN) :: discr_id
            CLASS(ppm_t_operator_discr_), POINTER :: discr

            INTEGER(HID_T) :: group_id
            INTEGER error

            CALL h5gcreate_f(cpfile_id, 'ops_t'//discr_id, group_id,&
                   error)

            ! integer dimension pointer order

            ! classes
            ! ppm_t_discr_kind discr_src disc_to
            ! ppm_t_operator_ op_ptr

            ! logical dimension 
            CALL store_logical_dim(group_id, 'flags', &
               discr%flags, ppm_param_length_opsflags)

            CALL h5gclose_f(group_id, error)
         END SUBROUTINE store_ppm_t_operator_discr_
