      PROGRAM test_cases
         USE checkpoint_hdf5
         USE ppm_module_init
         IMPLICIT NONE


         CLASS(ppm_t_part_mapping_d_), POINTER :: map
         CLASS(ppm_c_part_mapping_d_), POINTER :: maplist
         CLASS(ppm_t_particles_d_), POINTER :: parts
         CLASS(particles_stats_d_), POINTER :: stats
         CLASS(ppm_t_neighlist_d_), POINTER :: neigh
         CHARACTER(32) :: pointer_addr
         INTEGER(HID_T) :: cpfile_id
         INTEGER :: error
         INTEGER(8), DIMENSION(:), POINTER :: datablock
         ALLOCATE(datablock(10))
         datablock(3) = 50
         !ALLOCATE(ppm_c_part_mapping_d::maplist)

         ALLOCATE(ppm_t_part_mapping_d::map)
         map%oldNpart = 5
         map%target_topoid = 10
         ALLOCATE(particles_stats_d::stats)
         ALLOCATE(ppm_t_neighlist_d::neigh)
         ALLOCATE(ppm_t_particles_d::parts)

         CALL make_checkpoint_file('test.h5', cpfile_id)
         parts%id = 10
         parts%name = "kevin"
         parts%maps => maplist
         CALL store_TYPE(cpfile_id, '10', parts)
         CALL store_TYPE(cpfile_id, '10', neigh)

         deallocate(parts)
         parts => recover_ppm_t_particles_d_(cpfile_id, "10", parts)
         WRITE(*,*) parts%name
         CALL close_checkpoint_file(cpfile_id, error)

         CALL open_checkpoint_file(cpfile_id, 'checkpoint.h5', error)
         DEALLOCATE (parts)
         parts => recover_ppm_t_particles_d_(cpfile_id, "p1", parts)
         WRITE(*,*) parts%isymm
         CALL close_checkpoint_file(cpfile_id, error)

         !pointer_addr = get_pointer(map)
         !WRITE(*,*) pointer_addr
         DEALLOCATE(map, stats, parts, neigh)
      END PROGRAM test_cases
