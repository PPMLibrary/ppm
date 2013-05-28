      ! TODO find source of segfaults here
            SUBROUTINE store_ppm_c_operator_discr_(cpfile_id, discr_id,&
                         discr)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=*), INTENT(IN) :: discr_id
               CLASS(ppm_c_operator_discr_), POINTER :: discr
               CLASS(ppm_t_operator_discr_), POINTER :: dl => null()
               INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE :: dims
               CHARACTER(LEN=ppm_char) :: buffer
               INTEGER, DIMENSION(1) :: idlist
               INTEGER :: error, id, rank = 1

               CALL h5gcreate_f(cpfile_id, 'ops/'//discr_id, group_id,&
                      error)

               id = 1
               dl => discr%begin()
               DO WHILE(associated(dl))
                  WRITE(buffer,*) id
                  buffer = adjustl(buffer)
                  CALL store_ppm_t_operator_discr_(cpfile_id, &
                     buffer, dl)
                  idlist(id) = id
                  id = id +1
                  dl => discr%next()
               ENDDO
               dims = (/discr%nb/)
               CALL h5ltmake_dataset_int_f(group_id, 'idlist', rank, &
                  dims, idlist, error)


               CALL h5gclose_f(group_id, error)
            END SUBROUTINE store_ppm_c_operator_discr_
