      MODULE checkpoint_hdf5
         USE hdf5
         USE h5lt
         USE ppm_module_core
         !INTERFACE store_type
            !MODULE PROCEDURE store_ppm_t_part_prop_d_, &
               !store_particles_stats_d_
         !END INTERFACE store_type
         INTERFACE store_type
         END INTERFACE
         INTERFACE store_collection
            MODULE PROCEDURE store_ppm_c_neighlist_d_
         END INTERFACE store_collection
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

               ! ppm_c_operator_discr_ (Collection)
               CALL h5gcreate_f(file_id, 'ops', group_id, error)
               CALL h5gclose_f(group_id, error)

               ! ppm_t_operator_discr_ (actual ops)
               CALL h5gcreate_f(file_id, 'ops_t', group_id, error)
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

            SUBROUTINE store_logical_dim(group_id, dname, buffer,&
                   length)
               INTEGER(HID_T), INTENT(IN) :: group_id
               CHARACTER(LEN=*), INTENT(IN) :: dname
               LOGICAL, DIMENSION(length) :: buffer
               INTEGER, INTENT(in) :: length
               CHARACTER, DIMENSION(length) :: charbuf

               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER(HID_T) :: space_id, dset_id
               INTEGER :: error, i, rank
               rank = 1
               dims = (/length/)

               DO i=1, length
                  IF (buffer(i)) THEN
                     charbuf(i) = 'T'
                  ELSE
                     charbuf(i) = 'F'
                  ENDIF
               ENDDO
               CALL h5screate_simple_f(rank, dims, space_id, error)
               CALL h5dcreate_f(group_id, dname, H5T_NATIVE_CHARACTER, &
                  space_id, dset_id, error)
               CALL h5dwrite_f(dset_id, H5T_NATIVE_CHARACTER, charbuf, &
                  dims, error)
               CALL h5dclose_f(dset_id, error)
               CALL h5sclose_f(space_id, error)
            END SUBROUTINE store_logical_dim

            SUBROUTINE checkpoint_type(data_pointer)
               CLASS(ppm_t_main_abstr), POINTER, INTENT(IN) :: &
                  data_pointer
            END SUBROUTINE checkpoint_type
            SUBROUTINE checkpoint_container(container_ptr)
               CLASS(ppm_t_container), POINTER, INTENT(IN) :: &
                  container_ptr
            END SUBROUTINE checkpoint_container
            !SUBROUTINE store_pointer(cpfid, loc, ident, vname, ptr)
               !INTEGER(HID_T), INTENT(IN) :: cpfid, loc
               !CHARACTER(LEN=*), INTENT(IN) :: ident, vname
               !POINTER :: ptr
               !INTEGER :: error
               !IF (associated(ptr)) THEN
                  !CALL store_type(loc, ident, ptr)
                  !CALL h5ltmake_dataset_string_f(loc, vname, ident, &
                     !error)
               !ELSE
                  !CALL h5ltmake_dataset_string_f(loc, vname, '0', error)
               !ENDIF
            !END SUBROUTINE


            ! Done except for pointers
            INCLUDE 'checkp/ppm_t_particles_check.f'

            ! Done except for pointers
            INCLUDE 'checkp/particles_stats_check.f'

            INCLUDE 'checkp/ppm_t_part_prop_check.f'

            INCLUDE 'checkp/ppm_c_part_prop_check.f'

            ! Done except for pointers
            INCLUDE 'checkp/ppm_c_neighlist_check.f'
            INCLUDE 'checkp/ppm_t_neighlist_check.f'

            ! Need to test these, default tests have a null
            ! collection
            INCLUDE 'checkp/ppm_c_part_mapping_check.f'
            INCLUDE 'checkp/ppm_t_part_mapping_check.f'

            INCLUDE 'checkp/ppm_v_main_abstr_check.f'

            ! Done except for pointers
            ! Need to test these, default tests have a null
            ! collection
            INCLUDE 'checkp/ppm_c_operator_discr_check.f'
            INCLUDE 'checkp/ppm_t_operator_discr_check.f'
      END MODULE
