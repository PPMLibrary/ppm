      PROGRAM test_cases
         USE checkpoint_hdf5
         USE ppm_module_init

         CLASS(ppm_t_part_mapping_d_), POINTER :: map
         CLASS(ppm_t_particles_d), POINTER :: parts
         INTEGER(HID_T) :: cpfile_id
         INTEGER :: error

         ALLOCATE(ppm_t_part_mapping_d::map)
         map%oldNpart = 5
         map%target_topoid = 10

         ALLOCATE(ppm_t_particles_d::parts)

         CALL make_checkpoint_file('test.h5', cpfile_id)
         CALL store_ppm_t_part_mapping_d_(cpfile_id, '10', map)
         CALL store_ppm_t_particles_d(cpfile_id, '10', parts)
         CALL close_checkpoint_file(cpfile_id, error)

      END PROGRAM test_cases
