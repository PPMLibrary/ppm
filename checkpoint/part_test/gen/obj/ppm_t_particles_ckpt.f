            SUBROUTINE store_ppm_t_particles_d(cpfile_id, particle_id, &
                        particle)
               IMPLICIT NONE

               INTEGER(HID_T) :: cpfile_id
               INTEGER(HID_T) :: group_id
               !CHARACTER (LEN=20) :: particle_id ! look at boundary case
               INTEGER :: error, rank = 0
               INTEGER(HID_T) :: space_id, dset_id
               INTEGER(HSIZE_T), DIMENSION(1) :: dims = (/0/)
               CHARACTER(LEN=*), INTENT(IN):: particle_id
               TYPE(ppm_t_particles_d), POINTER :: particle
               INTEGER(HSIZE_T), DIMENSION(2) :: dims2
               CHARACTER, DIMENSION(ppm_param_length_partflags) :: chfl
               INTEGER :: I ! loop variable

               ! Identifiers for class pointers
               CHARACTER(LEN=20) :: tprop, cprop, pstat, abstr, discr, &
                  mapping, neighborlist


               ! Change later
               !CHARACTER (LEN=*), intent(IN) :: particle

               ! Create a group for our particle
               CALL h5gcreate_f(cpfile_id, particle_id, group_id, error)

               ! Store the classes
               pstat = particle_id//'_stat'
               CALL store_particles_stats_d_(cpfile_id, pstat, &
                        particle%stats)
               CALL h5ltmake_dataset_string_f(group_id, 'stat', pstat, &
                        error)

               ! Store the integer values of the particle
               CALL h5ltmake_dataset_int_f(group_id, 'Npart', rank, &
                  dims, (/particle%Npart/), error)
               CALL h5ltmake_dataset_int_f(group_id, 'Mpart', rank, &
                  dims, (/particle%Mpart/), error)
               CALL h5ltmake_dataset_int_f(group_id, 'active_topoid', &
                  rank, dims, (/particle%active_topoid/), error)
               CALL h5ltmake_dataset_int_f(group_id, 'isymm', &
                  rank, dims, (/particle%isymm/), error)
               CALL h5ltmake_dataset_int_f(group_id, 'NewNpart', &
                  rank, dims, (/particle%NewNpart/), error)
               CALL h5ltmake_dataset_int_f(group_id, 'itime', &
                  rank, dims, (/particle%itime/), error)


               ! Store the Real arguements of the particle
               CALL h5ltmake_dataset_double_f(group_id, 'time', &
                  rank, dims, (/particle%time/), error)
               CALL h5ltmake_dataset_double_f(group_id, 'ghostlayer', &
                  rank, dims, (/particle%ghostlayer/), error)
               CALL h5ltmake_dataset_double_f(group_id, 'h_avg', &
                  rank, dims, (/particle%h_avg/), error)
               CALL h5ltmake_dataset_double_f(group_id, 'h_min', &
                  rank, dims, (/particle%h_min/), error)

               ! Store the subtypes
               dims2 = shape(particle%xp)
               rank = 2
               CALL h5ltmake_dataset_double_f(group_id, 'xp', &
                  rank, dims2, particle%xp, error)
               rank = 1
               dims = shape(particle%pcost)
               CALL h5ltmake_dataset_double_f(group_id, 'pcost', &
                  rank, dims, particle%pcost, error)

               ! Convert flags to char and store
               DO i=1, ppm_param_length_partflags
                  IF (particle%flags(i)) THEN
                     chfl(i) = 'T'
                  ELSE
                     chfl(i) = 'F'
                  END IF
               ENDDO
               dims = (/ppm_param_length_partflags/)
               CALL h5screate_simple_f(rank, dims, space_id, error)
               CALL h5dcreate_f(group_id, "flags", H5T_NATIVE_CHARACTER&
                  , space_id, dset_id, error)
               CALL h5dwrite_f(dset_id, H5T_NATIVE_CHARACTER, chfl, &
                     dims, error)
               CALL h5dclose_f(dset_id, error)
               CALL h5sclose_f(space_id, error)

               ! Close the group and file
               CALL h5gclose_f(group_id, error)
            END SUBROUTINE
