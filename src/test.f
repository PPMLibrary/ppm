      PROGRAM test_cases
         USE checkpoint_hdf5
         USE ppm_module_init
         IMPLICIT NONE


         TYPE(ppm_t_particles_d), POINTER :: dumb
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
         pointer_addr = get_pointer(parts)
         CALL store_TYPE(cpfile_id, pointer_addr, parts)
         CALL save_pointer(cpfile_id, 'parts', pointer_addr)
         !CALL store_TYPE(cpfile_id, '10', neigh)
         deallocate(parts)
         CALL close_checkpoint_file(cpfile_id, error)

         CALL open_checkpoint_file(cpfile_id, 'test.h5', error)
         CALL get_saved_pointer(cpfile_id, 'parts', pointer_addr)
         parts => recover_ppm_t_particles_d_(cpfile_id, pointer_addr,
     C      parts)
         IF (associated(parts)) THEN
            WRITE(*,*) parts%name
            DEALLOCATE (parts)
         ELSE
            WRITE(*,*) "Recovery test 1 failed"
         ENDIF
         CALL close_checkpoint_file(cpfile_id, error)

         !CALL open_checkpoint_file(cpfile_id, 'checkpoint.h5', error)
         !pointer_addr = "30788A02000000004096700000000000"
         !parts => recover_ppm_t_particles_d_(cpfile_id,
        !   pointer_addr, parts)
         !IF (associated(parts)) THEN
         !   WRITE(*,*) parts%isymm
         !   WRITE(*,*) parts%stats%nb_ghost_push
         !   WRITE (*,*) parts%neighs%vec_size
         !   WRITE(*,*) parts%neighs%vec(1)%t%isymm
         !   SELECT TYPE (parts)
         !   TYPE is(ppm_t_particles_d)
         !      neigh => parts%neighs%begin()
         !      WRITE(*,*) neigh%nneighmax
         !      WRITE(*,*) (parts%xp(1,error), error=1,20)
         !      !neigh => parts%get_neighlist()
         !   END SELECT
         !   DEALLOCATE(parts)
         !ELSE
         !   WRITE(*,*) "Recovery test 2 failed"
         !ENDIF
         !CALL close_checkpoint_file(cpfile_id, error)

         !pointer_addr = get_pointer(map)
         DEALLOCATE(map, stats, neigh)
      END PROGRAM test_cases
