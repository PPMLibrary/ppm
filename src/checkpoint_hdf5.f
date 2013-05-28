      MODULE checkpoint_hdf5
         USE hdf5
         USE h5lt
         USE ppm_module_core
         CONTAINS
            SUBROUTINE make_checkpoint_file(filename, file_id)
               IMPLICIT NONE

               ! The filename of the checkpoint
               CHARACTER(LEN=*), INTENT(IN) :: filename

               INTEGER(HID_T), INTENT(out):: file_id
               INTEGER(HID_T) :: group_id
               INTEGER :: error

               ! Open the interface and create the file
               CALL h5open_f(error)
               CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, &
                        error)
               ! Create the required storage groups
               ! Use a group to encapsulate all objects of a given class

               ! ppm_t_particles
               CALL h5gcreate_f(file_id, 'particles', group_id, error)
               CALL h5gclose_f(group_id, error)

               ! ppm_particle_stats
               CALL h5gcreate_f(file_id, 'stats', group_id, error)
               CALL h5gclose_f(group_id, error)

               ! ppm_c_part_mapping
               CALL h5gcreate_f(file_id, 'maps', group_id, error)
               CALL h5gclose_f(group_id, error)

               ! ppm_v_main_abstr
               CALL h5gcreate_f(file_id, 'field_ptr', group_id, error)
               CALL h5gclose_f(group_id, error)

               ! ppm_c_operator_discr_
               CALL h5gcreate_f(file_id, 'ops', group_id, error)
               CALL h5gclose_f(group_id, error)

               ! ppm_c_neighlist (collections of neighlists)
               CALL h5gcreate_f(file_id, 'neighs', group_id, error)
               CALL h5gclose_f(group_id, error)

               ! ppm_t_neighlist (the actual neighlists)
               CALL h5gcreate_f(file_id, 'neighlists', group_id, error)
               CALL h5gclose_f(group_id, error)

               ! ppm_c_part_prop
               CALL h5gcreate_f(file_id, 'props', group_id, error)
               CALL h5gclose_f(group_id, error)

               ! ppm_t_part_prop
               CALL h5gcreate_f(file_id, 'gi', group_id, error)
               CALL h5gclose_f(group_id, error)
            END SUBROUTINE make_checkpoint_file

            SUBROUTINE close_checkpoint_file(file_id, error)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: file_id
               INTEGER error
               CALL h5fclose_f(file_id, error)
               CALL h5close_f(error)
            END SUBROUTINE close_checkpoint_file

            INCLUDE 'checkp/ppm_t_particles_check.f'

            INCLUDE 'checkp/particles_stats_check.f'

            INCLUDE 'checkp/ppm_t_part_prop_check.f'

            INCLUDE 'checkp/ppm_c_part_prop_check.f'

            INCLUDE 'checkp/ppm_c_neighlist_check.f'

            INCLUDE 'checkp/ppm_t_neighlist_check.f'

            INCLUDE 'checkp/ppm_c_part_mapping_check.f'

            INCLUDE 'checkp/ppm_v_main_abstr_check.f'

            INCLUDE 'checkp/ppm_c_operator_discr_check.f'
      END MODULE
