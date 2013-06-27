      MODULE checkpoint_hdf5
         USE hdf5
         USE h5lt
         USE ppm_module_core
         USE ppm_module_mapping_typedef
         USE pointer_tracker
         USE ISO_C_BINDING

         ! Generic interface definitions for abstract storage functions
         INCLUDE 'intrinsic/interfaces.f'

         ! Store a defined Datatype
         ! These interfaces have issues with arguments that
         ! derive one another
         INTERFACE store_type
            MODULE PROCEDURE store_ppm_t_main_abstr, &
                  store_particles_stats_d_, &
                  store_ppm_t_mapping_d_, &
                  store_ppm_t_container, &
                  store_ppm_t_operator_discr_, &
                  store_logical1d_pointer, &
                  store_logical2d_pointer, &
                  store_integer1d_pointer, &
                  store_integer2d_pointer, &
                  store_complex1d_pointer, &
                  store_complex2d_pointer, &
                  store_integer64_1d_pointer, &
                  store_integer64_2d_pointer, &
                  store_real1d_pointer, &
                  store_real2d_pointer, &
                  store_ppm_t_neighlist_d_, &
                  store_ppm_t_ptr_main_abstr, &
                  store_ppm_t_ptr_neighlist_d, &
                  store_ppm_t_ptr_part_mapping_d, &
                  store_ppm_t_ptr_part_prop_d, &
                  store_ppm_t_ptr_operator_discr
         END INTERFACE

!         INTERFACE recover_type
!            MODULE PROCEDURE recover_ppm_t_ptr_main_abstr, &
!                  recover_ppm_t_ptr_neighlist_d, &
!                  recover_ppm_t_ptr_part_mapping_d, &
!                  recover_ppm_t_ptr_part_prop_d, &
!                  recover_ppm_t_ptr_operator_discr, &
!                  recover_integer1d_pointer, &
!                  recover_integer2d_pointer, &
!                  recover_logical1d_pointer, &
!                  recover_logical2d_pointer, &
!                  recover_complex1d_pointer, &
!                  recover_complex2d_pointer, &
!                  recover_integer64_1d_pointer, &
!                  recover_integer64_2d_pointer, &
!                  recover_real1d_pointer, &
!                  recover_real2d_pointer, &
!                  INCLUDE 'types/read_derived.f'
!!                  read_ppm_t_main_abstr, &
!!                  read_particles_stats_d_, &
!!                  read_ppm_t_mapping_d_, &
!!                  read_ppm_t_container, &
!!                  read_ppm_t_operator_discr_, &
!         END INTERFACE recover_type
!
         CONTAINS
            ! Creates and initializes the checkpoint file
            ! Creates necessary groups in the file for types
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

               CALL h5gcreate_f(file_id, 'ppm_t_main_abstr', group_id,&
                  error)
               CALL h5gclose_f(group_id, error)

               CALL h5gcreate_f(file_id, 'ppm_t_operator_discr_', group_id,&
                  error)
               CALL h5gclose_f(group_id, error)

               CALL h5gcreate_f(file_id, 'particles_stats_d_', group_id, error)
               CALL h5gclose_f(group_id, error)

               CALL h5gcreate_f(file_id, 'ppm_t_mapping_d_', &
                  group_id, error)
               CALL h5gclose_f(group_id, error)

               ! ppm_c_neighlist (collections of neighlists)
               CALL h5gcreate_f(file_id, 'ppm_t_container', group_id, error)
               CALL h5gclose_f(group_id, error)

               ! ppm_t_neighlist (the actual neighlists)
               CALL h5gcreate_f(file_id, 'ppm_t_neighlist_d_', group_id, error)
               CALL h5gclose_f(group_id, error)

               ! for the intrinsic types
               CALL h5gcreate_f(file_id, 'intrinsic', group_id, error)
               CALL h5gclose_f(group_id, error)

               ! for the ptr types
               CALL h5gcreate_f(file_id, 'ptr_lists', group_id, error)
               CALL h5gclose_f(group_id, error)

            END SUBROUTINE make_checkpoint_file

            ! Opens an existing checkpoint file
            SUBROUTINE open_checkpoint_file(file_id, filename, error)
               IMPLICIT NONE

               ! The filename of the checkpoint
               CHARACTER(LEN=*), INTENT(IN) :: filename

               INTEGER(HID_T), INTENT(out):: file_id
               INTEGER :: error

               ! Open the interface and create the file
               CALL h5open_f(error)
               CALL h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, &
                        error)
            END SUBROUTINE open_checkpoint_file

            ! Close the checkpoint file
            ! Called at the end of the checkpoint
            SUBROUTINE close_checkpoint_file(file_id, error)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: file_id
               INTEGER error
               CALL h5fclose_f(file_id, error)
               CALL h5close_f(error)
               CALL delete_tree(pointer_data%itree)
               CALL delete_tree(pointer_data%dtree)
            END SUBROUTINE close_checkpoint_file

            ! pointer definitions for aquiring address of derived type
            INCLUDE 'pointers/ppm_t_ptr.f'

            ! subroutines for writing individual attributes of intrinsic
            ! types
            INCLUDE 'intrinsic/write_attributes.f'
            INCLUDE 'intrinsic/read_attributes.f'

            ! subroutines for reading intrinsic type pointers
            !INCLUDE 'intrinsic/read_types.f'
            INCLUDE 'intrinsic/recover_types.f'

            ! subroutines for the storage of pointer references to
            ! intrinsic types
            INCLUDE 'intrinsic/store_types.f'

            ! derived type definitions for hdf5 storage along with the
            ! routines to store the types
            INCLUDE 'types/ppm_t_operator_discr_check.f'
            INCLUDE 'types/ppm_t_main_abstr_check.f'
            INCLUDE 'types/ppm_t_mapping_check.f'
            INCLUDE 'types/ppm_t_container_check.f'
            INCLUDE 'types/ppm_t_neighlist_check.f'
            INCLUDE 'types/particles_stats_check.f'
      END MODULE
