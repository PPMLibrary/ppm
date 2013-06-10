            SUBROUTINE store_ppm_c_part_mapping_d(cpfile_id,&
                         mapping_id, mapping)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=*), INTENT(IN) :: mapping_id
               CLASS(ppm_c_part_mapping_d_), POINTER :: mapping
               CLASS(ppm_t_part_mapping_d_), POINTER :: mapl
               CHARACTER(LEN=ppm_char) :: buffer
               INTEGER, DIMENSION(mapping%nb) :: idlist


               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER error, id, rank

               CALL h5gcreate_f(cpfile_id, 'maps/'//mapping_id, &
                     group_id, error)

               id = 1
               rank = 1
               mapl => mapping%begin()
               DO WHILE(associated(mapl))
                  WRITE(buffer,*) id
                  buffer = adjustl(buffer)
                  CALL store_ppm_t_mapping_d_(cpfile_id, &
                     buffer, mapl)
                  idlist(id) = id
                  id = id +1
                  mapl => mapping%next()
               ENDDO
               dims = (/mapping%nb/)
               CALL h5ltmake_dataset_int_f(group_id, 'idlist', rank, &
                  dims, idlist, error)

               CALL h5gclose_f(group_id, error)
            END SUBROUTINE store_ppm_c_part_mapping_d
