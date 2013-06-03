            SUBROUTINE make_type_ppm_t_container(dtype_id)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(OUT) :: dtype_id
               INTEGER error 
               INTEGER(HSIZE_T) ::  typei, tsize, offset

               offset = 0
               tsize = 0

               ! Calculate datatype size
               CALL h5tget_size_f(H5T_NATIVE_INTEGER, typei, error)
               tsize = tsize + typei * 5

               ! Create the datatype
               CALL h5tcreate_f(H5T_COMPOUND_F, tsize, dtype_id, &
                  error)

               ! Insert the members
               CALL h5tinsert_f(dtype_id, "nb", offset, &
                  H5T_NATIVE_INTEGER, error)
               offset = offset + typei
               CALL h5tinsert_f(dtype_id, "iter_id", offset, &
                  H5T_NATIVE_INTEGER, error)
               offset = offset + typei
               CALL h5tinsert_f(dtype_id, "min_id", offset, &
                  H5T_NATIVE_INTEGER, error)
               offset = offset + typei
               CALL h5tinsert_f(dtype_id, "max_id", offset, &
                  H5T_NATIVE_INTEGER, error)
               offset = offset + typei
               CALL h5tinsert_f(dtype_id, "vec_size", offset, &
                  H5T_NATIVE_INTEGER, error)
               offset = offset + typei

            END SUBROUTINE make_type_ppm_t_container

            ! Abstract class this is likely unecessary
            SUBROUTINE store_ppm_t_container(cpfile_id, &
                  cont_id, cont)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id, type_id, dset_id, &
                  dspace_id
               CHARACTER(LEN=*), INTENT(IN) :: cont_id
               CLASS(ppm_t_container), POINTER :: cont
               INTEGER error

               CALL h5gopen_f(cpfile_id, 'ppm_t_container', &
                  group_id, error)


               ! Make our dataset
               CALL make_type_ppm_t_container(type_id) ! get type
               CALL h5screate_f(H5S_SCALAR_F, dspace_id, error)  ! get space
               CALL h5dcreate_f(group_id, cont_id, type_id, &
                   dspace_id, dset_id, error)

               CALL write_ppm_t_container(dset_id, cont)

               CALL h5dclose_f(dset_id, error)
               CALL h5sclose_f(dspace_id, error)

               CALL h5gclose_f(group_id, error)
            END SUBROUTINE store_ppm_t_container

            SUBROUTINE write_ppm_t_container(dset_id, cont)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: dset_id
               CLASS(ppm_t_container), POINTER :: cont

               CALL write_attribute(dset_id, "nb", &
                  cont%nb)
               CALL write_attribute(dset_id, "iter_id", &
                  cont%iter_id)
               CALL write_attribute(dset_id, "min_id", &
                  cont%min_id)
               CALL write_attribute(dset_id, "max_id", &
                  cont%max_id)
               CALL write_attribute(dset_id, "vec_size", &
                  cont%vec_size)


            END SUBROUTINE write_ppm_t_container
