            SUBROUTINE store_particles_stats_d_(cpfile_id, &
                     pstat_id, stats)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=20) :: pstat_id
               CLASS(particles_stats_d_), INTENT(IN), POINTER :: stats
               INTEGER :: rank = 0, error
               INTEGER(HSIZE_T) , DIMENSION(1) :: dims = (/0/)

               ! Create a group for our particle stats
               !CALL h5gcreate_f(cpfile_id, 'stats', stat_group, error)
               CALL h5gcreate_f(cpfile_id, 'stats/'//pstat_id,group_id,&
                     error)

               ! Store the integer arguments
               CALL h5ltmake_dataset_int_f(group_id, 'nb_nl', &
                  rank, dims, (/stats%nb_nl/), error)
               CALL h5ltmake_dataset_int_f(group_id, 'nb_inl', &
                  rank, dims, (/stats%nb_inl/), error)
               CALL h5ltmake_dataset_int_f(group_id, 'nb_cinl', &
                  rank, dims, (/stats%nb_cinl/), error)
               CALL h5ltmake_dataset_int_f(group_id, 'nb_xset_inl', &
                  rank, dims, (/stats%nb_xset_inl/), error)
               CALL h5ltmake_dataset_int_f(group_id, 'nb_xset_nl', &
                  rank, dims, (/stats%nb_xset_nl/), error)
               CALL h5ltmake_dataset_int_f(group_id, 'nb_dc_comp', &
                  rank, dims, (/stats%nb_dc_comp/), error)
               CALL h5ltmake_dataset_int_f(group_id, 'nb_dc_apply', &
                  rank, dims, (/stats%nb_dc_apply/), error)
               CALL h5ltmake_dataset_int_f(group_id, 'nb_kdtree', &
                  rank, dims, (/stats%nb_kdtree/), error)
               CALL h5ltmake_dataset_int_f(group_id, 'nb_global_map', &
                  rank, dims, (/stats%nb_global_map/), error)
               CALL h5ltmake_dataset_int_f(group_id, 'nb_part_map', &
                  rank, dims, (/stats%nb_part_map/), error)
               CALL h5ltmake_dataset_int_f(group_id, 'nb_ghost_get', &
                  rank, dims, (/stats%nb_ghost_get/), error)
               CALL h5ltmake_dataset_int_f(group_id, 'nb_ghost_push', &
                  rank, dims, (/stats%nb_ghost_push/), error)
               CALL h5ltmake_dataset_int_f(group_id, 'nb_ls', &
                  rank, dims, (/stats%nb_ls/), error)

               ! Next we store the Real arguments
               CALL h5ltmake_dataset_double_f(group_id, 't_nl', &
                  rank, dims, (/stats%t_nl/), error)
               CALL h5ltmake_dataset_double_f(group_id, 't_inl', &
                  rank, dims, (/stats%t_inl/), error)
               CALL h5ltmake_dataset_double_f(group_id, 't_cinl', &
                  rank, dims, (/stats%t_cinl/), error)
               CALL h5ltmake_dataset_double_f(group_id, 't_xset_inl', &
                  rank, dims, (/stats%t_xset_inl/), error)
               CALL h5ltmake_dataset_double_f(group_id, 't_xset_nl', &
                  rank, dims, (/stats%t_xset_nl/), error)
               CALL h5ltmake_dataset_double_f(group_id, 't_dc_comp', &
                  rank, dims, (/stats%t_dc_comp/), error)
               CALL h5ltmake_dataset_double_f(group_id, 't_dc_apply', &
                  rank, dims, (/stats%t_dc_apply/), error)
               CALL h5ltmake_dataset_double_f(group_id, 't_kdtree', &
                  rank, dims, (/stats%t_kdtree/), error)
               CALL h5ltmake_dataset_double_f(group_id, 't_global_map',&
                  rank, dims, (/stats%t_global_map/), error)
               CALL h5ltmake_dataset_double_f(group_id, 't_part_map', &
                  rank, dims, (/stats%t_part_map/), error)
               CALL h5ltmake_dataset_double_f(group_id, 't_ghost_get', &
                  rank, dims, (/stats%t_ghost_get/), error)
               CALL h5ltmake_dataset_double_f(group_id, 't_ghost_push',&
                  rank, dims, (/stats%t_ghost_push/), error)
               CALL h5ltmake_dataset_double_f(group_id, 't_ls', &
                  rank, dims, (/stats%t_ls/), error)
               CALL h5ltmake_dataset_double_f(group_id, 't_add', &
                  rank, dims, (/stats%t_add/), error)
               CALL h5ltmake_dataset_double_f(group_id, 't_del', &
                  rank, dims, (/stats%t_del/), error)
               CALL h5ltmake_dataset_double_f(group_id, 't_compD', &
                  rank, dims, (/stats%t_compD/), error)

               ! close the group
               CALL h5gclose_f(group_id, error)
            END SUBROUTINE store_particles_stats_d_
