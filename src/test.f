      PROGRAM test_cases
         USE checkpoint_hdf5
         USE ppm_module_init
         IMPLICIT NONE


         CLASS(ppm_t_part_mapping_d_), POINTER :: map
         CLASS(ppm_t_particles_d), POINTER :: parts
         CLASS(particles_stats_d_), POINTER :: stats
         CLASS(ppm_t_neighlist_d_), POINTER :: neigh
         CHARACTER(32) :: pointer_addr
         INTEGER(HID_T) :: cpfile_id
         INTEGER :: error
         INTEGER, DIMENSION(:), POINTER :: datablock
         ALLOCATE(datablock(10))

         ALLOCATE(ppm_t_part_mapping_d::map)
         map%oldNpart = 5
         map%target_topoid = 10
         ALLOCATE(particles_stats_d::stats)
         ALLOCATE(ppm_t_neighlist_d::neigh)
         neigh%nvlist => datablock
         ALLOCATE(ppm_t_particles_d::parts)

         CALL make_checkpoint_file('test.h5', cpfile_id)
         !CALL store_ppm_t_part_mapping_d_(cpfile_id, '10', map)
         parts%stats => stats
         pointer_addr = get_pointer(parts)
         WRITE(*,*) pointer_addr
         CALL store_ppm_t_particles_d(cpfile_id, '10', parts)
         !CALL store_ppm_t_neighlist_d_(cpfile_id, '10', neigh)
         CALL close_checkpoint_file(cpfile_id, error)

         !pointer_addr = get_pointer(map)
         !WRITE(*,*) pointer_addr
      END PROGRAM test_cases
