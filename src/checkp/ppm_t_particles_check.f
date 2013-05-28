            SUBROUTINE store_ppm_t_particles_d(cpfile_id, particle_id, &
                        particle)
               IMPLICIT NONE

               INTEGER(HID_T) :: cpfile_id
               INTEGER(HID_T) :: group_id
               !CHARACTER (LEN=20) :: particle_id ! look at boundary case
               INTEGER :: error, rank = 0
               INTEGER(HSIZE_T), DIMENSION(1) :: dims = (/0/)
               CHARACTER(LEN=*), INTENT(IN):: particle_id
               TYPE(ppm_t_particles_d), POINTER :: particle
               INTEGER(HSIZE_T), DIMENSION(2) :: dims2

               ! Identifiers for class pointers
               CHARACTER(LEN=20) :: tprop, cprop, pstat, abstr, discr, &
                  mapping, neighborlist


               ! Create a group for our particle
               CALL h5gcreate_f(cpfile_id, 'particles/'//particle_id, &
                  group_id, error)

               ! Store the classes
               ! Replace the ID assign line when we have reference
               ! tracking
               WRITE(pstat,'(I20)') 1
               pstat = adjustl(pstat)
               IF (associated(particle%stats)) THEN
                  CALL store_particles_stats_d_(cpfile_id, pstat, &
                           particle%stats)
                  CALL h5ltmake_dataset_string_f(group_id, 'stat', &
                           pstat, error)
               ELSE
                  CALL h5ltmake_dataset_string_f(group_id, 'stat', &
                           '0', error)
               ENDIF

               !tprop = particle_id//'_tprop'
               tprop = pstat
               IF (associated(particle%gi)) THEN
                  CALL store_ppm_t_part_prop_d_(cpfile_id, tprop, &
                     particle%gi)
                  CALL h5ltmake_dataset_string_f(group_id, 'gi', tprop,&
                     error)
               ELSE
                  CALL h5ltmake_dataset_string_f(group_id, 'gi', '0',&
                     error)
               ENDIF

               !cprop = particle_id//'_stat'
               cprop = pstat
               IF (associated(particle%props)) THEN
                  CALL store_ppm_c_part_prop_d_(cpfile_id, cprop, &
                     particle%props)
                  CALL h5ltmake_dataset_string_f(group_id, 'props', &
                     cprop, error)
               ELSE
                  CALL h5ltmake_dataset_string_f(group_id, 'props', '0'&
                     , error)
               ENDIF

               !abstr = particle_id//'_field_ptr'
               abstr = pstat
               IF (associated(particle%field_ptr)) THEN
                  CALL store_ppm_v_main_abstr(cpfile_id, abstr, &
                     particle%field_ptr)
                  CALL h5ltmake_dataset_string_f(group_id, 'field_ptr',&
                     abstr, error)
               ELSE
                  CALL h5ltmake_dataset_string_f(group_id, 'field_ptr',&
                     '0', error)
               ENDIF

               !neighborlist = particle_id//'neighs'
               neighborlist = pstat
               IF (associated(particle%neighs)) THEN
                  CALL store_ppm_c_neighlist_d_(cpfile_id, abstr, &
                     particle%neighs)
                  CALL h5ltmake_dataset_string_f(group_id, 'neighs', &
                     neighborlist, error)
               ELSE
                  CALL h5ltmake_dataset_string_f(group_id, 'neighs', &
                     '0', error)
               ENDIF

               !discr = particle_id//'ops'
               discr = pstat
               IF (associated(particle%ops)) THEN
                  CALL store_ppm_c_operator_discr_(cpfile_id, discr, &
                     particle%ops)
                  CALL h5ltmake_dataset_string_f(group_id, 'ops', &
                     discr, error)
               ELSE
                  CALL h5ltmake_dataset_string_f(group_id, 'ops', &
                     '0', error)
               ENDIF

               !mapping = particle_id//'maps'
               mapping = pstat
               IF (associated(particle%maps)) THEN
                  CALL store_ppm_c_part_mapping_d(cpfile_id, mapping, &
                     particle%maps)
                  CALL h5ltmake_dataset_string_f(group_id, 'maps', &
                     mapping, error)
               ELSE
                  CALL h5ltmake_dataset_string_f(group_id, 'maps', &
                     '0', error)
               ENDIF

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

               CALL store_logical_dim(group_id, 'flags', &
                     particle%flags, ppm_param_length_partflags)

               ! Close the group and file
               CALL h5gclose_f(group_id, error)
            END SUBROUTINE store_ppm_t_particles_d
