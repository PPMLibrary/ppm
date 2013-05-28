            SUBROUTINE store_ppm_c_neighlist_d_(cpfile_id, &
                neighlist_id, neighlist)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=*), INTENT(IN) :: neighlist_id
               CLASS(ppm_c_neighlist_d_), POINTER :: neighlist
               CLASS(ppm_t_neighlist_d_), POINTER :: nl
               INTEGER :: error, rank = 1
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               CHARACTER(LEN=ppm_char) :: buffer
               !CHARACTER(LEN=ppm_char), DIMENSION(neighlist%nb) :: buff
               INTEGER, DIMENSION(neighlist%nb) :: idlist
               INTEGER :: id

               CALL h5gcreate_f(cpfile_id, 'neighs/'//neighlist_id, &
                     group_id, error)

               ! TODO
               ! iterate the collection
               ! store each ppm_t_neighlist
               nl => neighlist%begin()
               id = 1 ! This will get fixed later
               DO WHILE(associated(nl))
                  WRITE(buffer,*) id
                  buffer = adjustl(buffer)
                  CALL store_ppm_t_neighlist(cpfile_id, &
                     buffer, nl)
                  idlist(id) = id
                  id = id + 1
                  nl => neighlist%next()
               ENDDO
               ! Store the integer ids of the neighborlists
               dims = (/neighlist%nb/)
               CALL h5ltmake_dataset_int_f(group_id, 'idlist', rank, &
                     dims, idlist, error)

               CALL h5gclose_f(group_id, error)
            END SUBROUTINE store_ppm_c_neighlist_d_
