      MODULE checkpoint_hdf5
         USE hdf5
         USE h5lt
         USE ppm_module_core
         USE pointer_tracker

         ! Generic interface definitions for abstract storage functions

         ! Write an attribute to a dataset
         INTERFACE write_attribute
            module procedure write_integer_attribute, &
                  write_double_attribute, &
                  write_logical_attribute, &
                  write_character_attribute, &
                  write_logical_array, &
                  write_character_array, &
                  write_string_attribute
         END INTERFACE write_attribute
         INTERFACE store_integer
            MODULE PROCEDURE store_integer_ptr, &
                  store_integer1d, &
                  store_integer2d
         END INTERFACE store_integer

         ! Store a defined Datatype
         ! These interfaces have issues with arguments that
         ! derive one another
         INTERFACE store_type
            MODULE PROCEDURE store_ppm_t_main_abstr, &
                  store_ppm_t_mapping_d_, &
                  !store_ppm_c_neighlist_d_, &
                  !store_ppm_c_part_prop_d_, &
                  store_ppm_t_container, &
                  store_ppm_t_operator_discr_, &
         !         store_ppm_t_discr_data, & is abstract
         !         store_ppm_t_discr_kind, & is abstract
         !         store_ppm_t_mapping_d_, & is abstract
                  store_ppm_t_neighlist_d_
         !         store_ppm_t_part_prop_d_, &
         !         store_ppm_v_main_abstr is abstract
         END INTERFACE

         INTERFACE write_type
            MODULE PROCEDURE write_ppm_t_main_abstr, &
                  write_ppm_t_mapping_d_
                  !write_ppm_t_mapping_d_
         END INTERFACE


         ! Generic abstraction in order to store generic containers
         !INTERFACE store_collection
         !   MODULE PROCEDURE store_ppm_c_neighlist_d_
         !END INTERFACE store_collection

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

               ! ppm_t_particles
               !CALL h5gcreate_f(file_id, 'ppm_t_particles_d', group_id,&
               !   error)
               !CALL h5gclose_f(group_id, error)
               CALL h5gcreate_f(file_id, 'ppm_t_main_abstr', group_id,&
                  error)
               CALL h5gclose_f(group_id, error)

               ! ppm_particle_stats
               !CALL h5gcreate_f(file_id, 'stats', group_id, error)
               CALL h5gcreate_f(file_id, 'particles_stats_d_', group_id, error)
               CALL h5gclose_f(group_id, error)

               ! ppm_c_part_mapping
               CALL h5gcreate_f(file_id, 'maps', group_id, error)
               CALL h5gclose_f(group_id, error)
               ! ppm_t_part_mapping
               CALL h5gcreate_f(file_id, 'ppm_t_mapping_d_', &
                  group_id, error)
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
               CALL h5gcreate_f(file_id, 'ppm_t_container', group_id, error)
               CALL h5gclose_f(group_id, error)

               ! ppm_c_neighlist (collections of neighlists)
               CALL h5gcreate_f(file_id, 'ppm_c_neighlist', group_id, error)
               CALL h5gclose_f(group_id, error)

               ! ppm_t_neighlist (the actual neighlists)
               CALL h5gcreate_f(file_id, 'ppm_t_neighlist_d_', group_id, error)
               CALL h5gclose_f(group_id, error)

               ! ppm_c_part_prop
               CALL h5gcreate_f(file_id, 'props', group_id, error)
               CALL h5gclose_f(group_id, error)

               ! ppm_t_part_prop
               CALL h5gcreate_f(file_id, 'gi', group_id, error)
               CALL h5gclose_f(group_id, error)

               ! other group for miscellaneous and fortran type pointers
               CALL h5gcreate_f(file_id, 'other', group_id, error)
               CALL h5gclose_f(group_id, error)

               ! Abstract groups, should remove these later if possible
               CALL h5gcreate_f(file_id, 'ppm_t_discr_data', group_id, error)
               CALL h5gclose_f(group_id, error)

            END SUBROUTINE make_checkpoint_file

            ! Opens an existing checkpoint file
            SUBROUTINE open_checkpoint_file(filename, file_id)
               IMPLICIT NONE

               ! The filename of the checkpoint
               CHARACTER(LEN=*), INTENT(IN) :: filename

               INTEGER(HID_T), INTENT(out):: file_id
               INTEGER :: error

               ! Open the interface and create the file
               CALL h5open_f(error)
               CALL h5fcreate_f(filename, H5F_ACC_RDWR_F, file_id, &
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
            END SUBROUTINE close_checkpoint_file

            ! This function will be removed for the genericized version
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

            ! Write a BOOL array attribute to a dataset
            SUBROUTINE write_logical_array(dset_id, &
                  dname, buffer, length)
               INTEGER(HID_T), INTENT(IN) :: dset_id
               CHARACTER(LEN=*), INTENT(IN) :: dname
               LOGICAL, DIMENSION(:) :: buffer
               INTEGER, INTENT(in) :: length
               CHARACTER, DIMENSION(length) :: charbuf
               !ALLOCATE(charbuf(length))

               INTEGER :: i

               DO i=1, length
                  IF (buffer(i)) THEN
                     charbuf(i) = 'T'
                  ELSE
                     charbuf(i) = 'F'
                  ENDIF
               ENDDO
               CALL write_character_array(dset_id, dname, &
                  charbuf, length)

            END SUBROUTINE write_logical_array

            SUBROUTINE write_character_array(dset_id, dname, &
                  charbuf, length)
               INTEGER(HID_T), INTENT(IN) :: dset_id
               CHARACTER(LEN=*), INTENT(IN) :: dname
               CHARACTER, DIMENSION(length) :: charbuf
               INTEGER, INTENT(in) :: length

               INTEGER(HID_T) :: attr_id, plist_id, bool_id
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER(HSIZE_T) :: csize, arr_size, offset
               INTEGER :: error, rank
               rank = 1
               dims = (/length/)
               offset = 0

               ! Create the type properties
               ! preserve partially initialized fields
               CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
               CALL h5pset_preserve_f(plist_id, .TRUE., error)

               ! Get the data size
               CALL h5tget_size_f(H5T_NATIVE_CHARACTER, csize, error)
               arr_size = length*csize

               ! Create the array type for the attribute
               CALL h5tarray_create_f(H5T_NATIVE_CHARACTER, rank, dims,&
                  bool_id, error)
               CALL h5tget_size_f(bool_id, arr_size, error)
               CALL h5tcreate_f(H5T_COMPOUND_F, arr_size, attr_id, &
                  error)
               CALL h5tinsert_f(attr_id, dname, offset, bool_id, error)

               CALL h5dwrite_f(dset_id, attr_id , charbuf, &
                  dims, error)
               CALL h5tclose_f(attr_id, error)
            END SUBROUTINE write_character_array

            SUBROUTINE write_string_attribute(dset_id, dname, charbuf, &
               length)
               INTEGER(HID_T), INTENT(IN) :: dset_id
               CHARACTER(LEN=*), INTENT(IN) :: dname
               CHARACTER(LEN=length) :: charbuf
               INTEGER, INTENT(in) :: length

               INTEGER(HID_T) :: attr_id, plist_id, bool_id
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER(HSIZE_T) :: csize, arr_size, offset
               INTEGER :: error, rank
               rank = 1
               dims = (/0/)
               offset = 0

               ! Create the type properties
               ! preserve partially initialized fields
               CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
               CALL h5pset_preserve_f(plist_id, .TRUE., error)

               ! Get the data size
               CALL h5tget_size_f(H5T_NATIVE_CHARACTER, csize, error)
               arr_size = length*csize

               ! Create the array type for the attribute
               CALL h5tcreate_f(H5T_STRING_F, arr_size, &
                  bool_id, error)
               CALL h5tget_size_f(bool_id, arr_size, error)
               CALL h5tcreate_f(H5T_COMPOUND_F, arr_size, attr_id, &
                  error)
               CALL h5tinsert_f(attr_id, dname, offset, bool_id, error)

               CALL h5dwrite_f(dset_id, attr_id , charbuf, &
                  dims, error)
               CALL h5tclose_f(attr_id, error)
            END SUBROUTINE write_string_attribute

            !SUBROUTINE checkpoint_type(data_pointer)
            !   CLASS(ppm_t_main_abstr), POINTER, INTENT(IN) :: &
            !      data_pointer
            !END SUBROUTINE checkpoint_type
            !SUBROUTINE checkpoint_container(container_ptr)
            !   CLASS(ppm_t_container), POINTER, INTENT(IN) :: &
            !      container_ptr
            !END SUBROUTINE checkpoint_container

            ! Write integer attribute to a dataset
            SUBROUTINE write_integer_attribute(dset_id, &
                  vname, val)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: dset_id
               CHARACTER(LEN=*), INTENT(IN) :: vname
               INTEGER :: val

               INTEGER error
               INTEGER(HID_T) :: plist_id, attr_id
               INTEGER(HSIZE_T) :: isize, offset
               INTEGER(HSIZE_T), DIMENSION(1) :: dims = (/0/)
               offset = 0

               ! Create the type properties
               ! preserve partially initialized fields
               CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
               CALL h5pset_preserve_f(plist_id, .TRUE., error)

               ! Get the data size
               CALL h5tget_size_f(H5T_NATIVE_INTEGER, isize, error)

               ! Create the compound type for the attribute
               CALL h5tcreate_f(H5T_COMPOUND_F, isize, attr_id, error)
               CALL h5tinsert_f(attr_id, vname, offset, &
                   H5T_NATIVE_INTEGER, error)
               CALL h5dwrite_f(dset_id, attr_id, val, dims, &
                     error, xfer_prp = plist_id)

               CALL h5tclose_f(attr_id, error)

            END SUBROUTINE write_integer_attribute

            ! Write double attribute to a dataset
            SUBROUTINE write_double_attribute(dset_id, vname, &
                  val)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: dset_id
               CHARACTER(LEN=*), INTENT(IN) :: vname
               REAL(ppm_kind_double) :: val

               INTEGER error
               INTEGER(HID_T) :: plist_id, attr_id
               INTEGER(HSIZE_T) :: dsize, offset
               INTEGER(HSIZE_T), DIMENSION(1) :: dims = (/0/)
               offset = 0

               ! Create the type properties
               ! preserve partially initialized fields
               CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
               CALL h5pset_preserve_f(plist_id, .TRUE., error)

               ! Get the data size
               CALL h5tget_size_f(H5T_NATIVE_DOUBLE, dsize, error)

               ! Create the compound type for the attribute
               CALL h5tcreate_f(H5T_COMPOUND_F, dsize, attr_id, error)
               CALL h5tinsert_f(attr_id, vname, offset, &
                   H5T_NATIVE_DOUBLE, error)
               CALL h5dwrite_f(dset_id, attr_id, val, dims, &
                     error, xfer_prp = plist_id)

               CALL h5tclose_f(attr_id, error)

            END SUBROUTINE write_double_attribute

            SUBROUTINE write_logical_attribute(dset_id, &
                  vname, val)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: dset_id
               CHARACTER(LEN=*), INTENT(IN) :: vname
               LOGICAL :: val
               CHARACTER :: c
               IF (val) THEN
                  c = 'T'
               ELSE
                  c = 'F'
               ENDIF
               CALL write_character_attribute(dset_id, vname, c)
            END SUBROUTINE write_logical_attribute

            SUBROUTINE write_character_attribute(dset_id, &
                  vname, val)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: dset_id
               CHARACTER(LEN=*), INTENT(IN) :: vname
               CHARACTER :: val

               INTEGER error
               INTEGER(HID_T) :: plist_id, attr_id
               INTEGER(HSIZE_T) :: csize, offset
               INTEGER(HSIZE_T), DIMENSION(1) :: dims = (/0/)
               offset = 0

               ! Create the type properties
               ! preserve partially initialized fields
               CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
               CALL h5pset_preserve_f(plist_id, .TRUE., error)

               ! Get the data size
               CALL h5tget_size_f(H5T_NATIVE_CHARACTER, csize, error)

               ! Create the compound type for the attribute
               CALL h5tcreate_f(H5T_COMPOUND_F, csize, attr_id, error)
               CALL h5tinsert_f(attr_id, vname, offset, &
                   H5T_NATIVE_CHARACTER, error)
               CALL h5dwrite_f(dset_id, attr_id, val, dims, &
                     error, xfer_prp = plist_id)

               CALL h5tclose_f(attr_id, error)

            END SUBROUTINE write_character_attribute
            SUBROUTINE store_integer_ptr(cpfile_id, mem_id, &
                  buffer)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=*), INTENT(IN) :: mem_id
               INTEGER :: rank
               INTEGER, INTENT(in) :: buffer
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER :: error

               rank = 0
               dims = (/0/)
               CALL h5gopen_f(cpfile_id, "other", group_id, error)
               CALL h5ltmake_dataset_int_f(group_id, mem_id, rank, &
                  dims, (/buffer/), error)

               CALL h5gclose_f(group_id, error)
            END SUBROUTINE store_integer_ptr
            SUBROUTINE store_integer1d(cpfile_id, mem_id, &
                  buffer)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=*), INTENT(IN) :: mem_id
               INTEGER :: rank
               INTEGER, DIMENSION(:), INTENT(in) :: buffer
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER :: error

               rank = 1
               dims = shape(buffer)
               CALL h5gopen_f(cpfile_id, "other", group_id, error)
               CALL h5ltmake_dataset_int_f(group_id, mem_id, rank, &
                  dims, buffer, error)

               CALL h5gclose_f(group_id, error)
            END SUBROUTINE store_integer1d
            SUBROUTINE store_integer2d(cpfile_id, mem_id, &
                  buffer)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id
               CHARACTER(LEN=*), INTENT(IN) :: mem_id
               INTEGER :: rank
               INTEGER, DIMENSION(:,:), INTENT(in) :: buffer
               INTEGER(HSIZE_T), DIMENSION(2) :: dims
               INTEGER :: error

               rank = 2
               dims = shape(buffer)
               CALL h5gopen_f(cpfile_id, "other", group_id, error)
               CALL h5ltmake_dataset_int_f(group_id, mem_id, rank, &
                  dims, buffer, error)

               CALL h5gclose_f(group_id, error)
            END SUBROUTINE store_integer2d
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


            ! Basic abstraction layer, no pointers
            ! Still in sub group for compatibility
            !INCLUDE 'checkp/ppm_t_particles_check.f'

            INCLUDE 'checkp/particles_stats_check.f'

            !INCLUDE 'checkp/ppm_t_part_prop_check.f'

            !INCLUDE 'checkp/ppm_c_part_prop_check.f'

            !INCLUDE 'checkp/ppm_c_neighlist_check.f'
            !INCLUDE 'checkp/ppm_t_neighlist_check.f'

            ! Need to test these, default tests have a null
            ! collection
            ! Needs container abstraction layer
            !INCLUDE 'checkp/ppm_c_part_mapping_check.f'
            ! Simple abstraction layer, no pointers
            !INCLUDE 'checkp/ppm_t_part_mapping_check.f'
            !INCLUDE 'checkp/ppm_t_mapping_check.f'

            !INCLUDE 'checkp/ppm_v_main_abstr_check.f'

            ! Need to test these, default tests have a null
            ! collection
            INCLUDE 'checkp/ppm_c_operator_discr_check.f'
            !INCLUDE 'checkp/ppm_t_operator_discr_check.f'
            INCLUDE 'poly/ppm_t_operator_discr_check.f'

            !INCLUDE 'checkp/ppm_t_container_check.f'

            !INCLUDE 'checkp/ppm_t_discr_kind_check.f'
            !INCLUDE 'checkp/ppm_t_discr_data_check.f'
            !INCLUDE 'checkp/ppm_t_main_abstr_check.f'
            INCLUDE 'poly/ppm_t_main_abstr_check.f'
            INCLUDE 'poly/ppm_t_mapping_check.f'
            INCLUDE 'poly/ppm_t_container_check.f'
            INCLUDE 'poly/ppm_t_neighlist_check.f'
      END MODULE
