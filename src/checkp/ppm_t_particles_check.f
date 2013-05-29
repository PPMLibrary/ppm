            SUBROUTINE make_ppm_t_particles_d_type(dtype_id)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(OUT) :: dtype_id
               INTEGER(HID_T) :: bool_id
               INTEGER :: error
               INTEGER(HSIZE_T) :: int_size, double_size, type_size, &
                  offset, char_size, bool_size
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER :: rank
               rank = 1
               dims = (/ppm_param_length_partflags/)

               CALL h5tarray_create_f(H5T_NATIVE_CHARACTER, rank, &
                  dims, bool_id, error)

               type_size = 0
               CALL h5tget_size_f(H5T_NATIVE_INTEGER, int_size, error)
               CALL h5tget_size_f(H5T_NATIVE_DOUBLE, double_size, error)
               CALL h5tget_size_f(H5T_NATIVE_CHARACTER, char_size, &
                     error)
               CALL h5tget_size_f(bool_id, bool_size, error)

               type_size = (6*int_size) + (4*double_size) + &
                  bool_size

               CALL h5tcreate_f(H5T_COMPOUND_F, type_size, dtype_id, &
                  error)

               CALL h5tinsert_f(dtype_id, "Npart", offset, &
                  H5T_NATIVE_INTEGER, error)
               offset = offset + int_size
               CALL h5tinsert_f(dtype_id, "Mpart", offset, &
                  H5T_NATIVE_INTEGER, error)
               offset = offset + int_size
               CALL h5tinsert_f(dtype_id, "active_topoid", offset, &
                  H5T_NATIVE_INTEGER, error)
               offset = offset + int_size
               CALL h5tinsert_f(dtype_id, "isymm", offset, &
                  H5T_NATIVE_INTEGER, error)
               offset = offset + int_size
               CALL h5tinsert_f(dtype_id, "NewNpart", offset, &
                  H5T_NATIVE_INTEGER, error)
               offset = offset + int_size
               CALL h5tinsert_f(dtype_id, "itime", offset, &
                  H5T_NATIVE_INTEGER, error)
               offset = offset + int_size

               ! Insert the doubles
               CALL h5tinsert_f(dtype_id, "time", offset, &
                  H5T_NATIVE_DOUBLE, error)
               offset = offset + double_size
               CALL h5tinsert_f(dtype_id, "ghostlayer", offset, &
                  H5T_NATIVE_DOUBLE, error)
               offset = offset + double_size
               CALL h5tinsert_f(dtype_id, "h_avg", offset, &
                  H5T_NATIVE_DOUBLE, error)
               offset = offset + double_size
               CALL h5tinsert_f(dtype_id, "h_min", offset, &
                  H5T_NATIVE_DOUBLE, error)
               offset = offset + double_size

               ! Finally the logical flags
               CALL h5tinsert_f(dtype_id, "flags", offset, &
                  bool_id, error)
               offset = offset + bool_size

            END SUBROUTINE
            SUBROUTINE store_ppm_t_particles_d(cpfile_id, particle_id, &
                        particle)
               IMPLICIT NONE

               INTEGER(HID_T) :: cpfile_id
               INTEGER(HID_T) :: group_id, space_id, dset_id, type_id
               !CHARACTER (LEN=20) :: particle_id ! look at boundary case
               INTEGER :: error, rank = 0
               INTEGER(HSIZE_T), DIMENSION(1) :: dims = (/0/)
               CHARACTER(LEN=*), INTENT(IN):: particle_id
               CLASS(ppm_t_particles_d), POINTER :: particle
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


               CALL make_ppm_t_particles_d_type(type_id)
               CALL h5screate_f(H5S_SCALAR_F, space_id, error)
               CALL h5dcreate_f(group_id, particle_id, type_id, &
                  space_id, dset_id, error)
               CALL write_ppm_t_particles_d(dset_id, particle)
               CALL h5dclose_f(dset_id, error)
               CALL h5sclose_f(space_id, error)
               ! Close the group and file
               CALL h5gclose_f(group_id, error)
            END SUBROUTINE store_ppm_t_particles_d

            SUBROUTINE write_ppm_t_particles_d(dset_id, particle)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: dset_id
               CLASS(ppm_t_particles_d), POINTER :: particle

               ! Integers
               CALL write_attribute(dset_id, 'Npart', &
                  particle%Npart)
               CALL write_attribute(dset_id, 'Mpart', &
                  particle%Mpart)
               CALL write_attribute(dset_id, 'active_topoid', &
                  particle%active_topoid)
               CALL write_attribute(dset_id, 'isymm', &
                  particle%isymm)
               CALL write_attribute(dset_id, 'NewNpart', &
                  particle%NewNpart)
               CALL write_attribute(dset_id, 'itime', &
                  particle%itime)

               ! Doubles
               CALL write_attribute(dset_id, 'time', &
                  particle%time)
               CALL write_attribute(dset_id, 'ghostlayer', &
                  particle%ghostlayer)
               CALL write_attribute(dset_id, 'h_avg', &
                  particle%h_avg)
               CALL write_attribute(dset_id, 'h_min', &
                  particle%h_min)

               ! Logicals
               CALL write_attribute(dset_id, "flags", &
                  particle%flags, ppm_param_length_partflags)
            END SUBROUTINE


