            SUBROUTINE store_ppm_t_neighlist(cpfile_id, &
                  neighlist_id, neighlist)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=*), INTENT(IN) :: neighlist_id
               CLASS(ppm_t_neighlist_d_), POINTER :: neighlist
               INTEGER error, rank
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               CHARACTER :: logbuf

               CALL h5gcreate_f(cpfile_id, 'neighlists/'//neighlist_id,&
                     group_id, error)

               ! Store the character values
               CALL h5ltmake_dataset_string_f(group_id, 'name', &
                  neighlist%name, error)

               ! Store the Integer Values next
               dims = (/0/)
               rank = 0
               CALL h5ltmake_dataset_int_f(group_id, 'isymm', &
                  rank, dims, (/neighlist%isymm/), error)
               CALL h5ltmake_dataset_int_f(group_id, 'nneighmin', &
                  rank, dims, (/neighlist%nneighmin/), error)
               CALL h5ltmake_dataset_int_f(group_id, 'nneighmax', &
                  rank, dims, (/neighlist%nneighmax/), error)
               ! TODO
               ! Rest are dimension pointers will do later
               ! vlist
               ! nvlist

               ! Store the Real arguments
               CALL h5ltmake_dataset_double_f(group_id, 'cutoff', &
                  rank, dims, (/neighlist%cutoff/), error)
               CALL h5ltmake_dataset_double_f(group_id, 'skin', &
                  rank, dims, (/neighlist%skin/), error)

               ! Convert and store the logical
               IF (neighlist%uptodate) THEN
                  logbuf = 'T'
               ELSE
                  logbuf = 'F'
               END IF
               CALL h5ltmake_dataset_string_f(group_id, 'uptodate', &
                     logbuf, error)

               CALL h5gclose_f(group_id, error)
            END SUBROUTINE store_ppm_t_neighlist
