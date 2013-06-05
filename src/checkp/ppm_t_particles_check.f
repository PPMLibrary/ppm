      ! Store routines for ppm_t_particles_d
            SUBROUTINE make_type_ppm_t_particles_d(dtype_id)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(OUT) :: dtype_id
               INTEGER(HID_T) :: array_id
               INTEGER(HSIZE_T), DIMENSION(1) :: dims
               INTEGER rank
               INTEGER error
               INTEGER(HSIZE_T) :: tsize, offset
               INTEGER(HSIZE_T) :: isize
               INTEGER(HSIZE_T) :: dsize
               INTEGER(HSIZE_T) :: csize

               offset = 0
               tsize = 0
               rank = 1


               ! Calculate datatype size

               CALL h5tget_size_f(H5T_NATIVE_INTEGER, isize, error)
               CALL h5tget_size_f(H5T_NATIVE_DOUBLE, dsize, error)
               CALL h5tget_size_f(H5T_NATIVE_CHARACTER, csize, error)


               tsize = tsize + isize*(5) + dsize*(4) + csize*(ppm_param_length_partflags+0) + csize*32*1



               ! Create/Expand the datatype
               CALL h5tcreate_f(H5T_COMPOUND_F, tsize, dtype_id, error)

               ! Insert the members
               ! Insert the members
               ! Integer members
               CALL h5tinsert_f(dtype_id, "NewNpart", offset, &
                     H5T_NATIVE_INTEGER, error)
               offset = offset + isize
               CALL h5tinsert_f(dtype_id, "isymm", offset, &
                     H5T_NATIVE_INTEGER, error)
               offset = offset + isize
               CALL h5tinsert_f(dtype_id, "Npart", offset, &
                     H5T_NATIVE_INTEGER, error)
               offset = offset + isize
               CALL h5tinsert_f(dtype_id, "itime", offset, &
                     H5T_NATIVE_INTEGER, error)
               offset = offset + isize
               CALL h5tinsert_f(dtype_id, "active_topoid", offset, &
                     H5T_NATIVE_INTEGER, error)
               offset = offset + isize

               ! Real members
               CALL h5tinsert_f(dtype_id, "h_min", offset, &
                     H5T_NATIVE_DOUBLE, error)
               offset = offset + dsize
               CALL h5tinsert_f(dtype_id, "h_avg", offset, &
                     H5T_NATIVE_DOUBLE, error)
               offset = offset + dsize
               CALL h5tinsert_f(dtype_id, "ghostlayer", offset, &
                     H5T_NATIVE_DOUBLE, error)
               offset = offset + dsize
               CALL h5tinsert_f(dtype_id, "time", offset, &
                     H5T_NATIVE_DOUBLE, error)
               offset = offset + dsize

               ! Character members

               ! Logical members
               dims = (/csize*ppm_param_length_partflags/)
               CALL h5tarray_create_f(H5T_NATIVE_CHARACTER, rank, &
                   dims, array_id, error)
               CALL h5tinsert_f(dtype_id, "flags", offset, &
                   array_id, error)
               offset = offset + (csize*ppm_param_length_partflags)

               ! Pointer members
               CALL h5tcreate_f(H5T_STRING_F, 32*csize, &
                   array_id, error)
               CALL h5tinsert_f(dtype_id, "stats", offset, &
                   array_id, error)
               offset = offset + (32*csize)



            END SUBROUTINE make_type_ppm_t_particles_d

            SUBROUTINE store_ppm_t_particles_d(cpfile_id, &
                  type_ptr_id, type_ptr)
               INTEGER(HID_T), INTENT(IN) :: cpfile_id
               INTEGER(HID_T) :: group_id, type_id, dset_id, &
                  dspace_id
               CHARACTER(LEN=*), INTENT(IN) :: type_ptr_id
               CLASS(ppm_t_particles_d), POINTER :: type_ptr
               INTEGER error

               CALL h5gopen_f(cpfile_id, 'ppm_t_particles_d', &
                  group_id, error)


               ! Make our dataset
               CALL make_type_ppm_t_particles_d(type_id) ! get type
               CALL h5screate_f(H5S_SCALAR_F, dspace_id, error)  ! get space
               CALL h5dcreate_f(group_id, type_ptr_id, type_id, &
                   dspace_id, dset_id, error)

               CALL write_ppm_t_particles_d(cpfile_id, dset_id, type_ptr)
               !CALL write_ppm_t_particles_d(dset_id, type_ptr)

               CALL h5dclose_f(dset_id, error)
               CALL h5sclose_f(dspace_id, error)

               CALL h5gclose_f(group_id, error)
            END SUBROUTINE store_ppm_t_particles_d

            SUBROUTINE write_ppm_t_particles_d(cpfile_id, dset_id, type_ptr)
            !SUBROUTINE write_ppm_t_particles_d(cpfile_id, dset_id, type_ptr)
               IMPLICIT NONE
               INTEGER(HID_T), INTENT(IN) :: dset_id
               INTEGER(HID_T), INTENT(in) :: cpfile_id
               CHARACTER(LEN=32) :: pointer_addr
               CLASS(ppm_t_particles_d), POINTER :: type_ptr

               CALL write_attribute(dset_id, 'NewNpart', &
                   type_ptr%NewNpart)
               CALL write_attribute(dset_id, 'isymm', &
                   type_ptr%isymm)
               CALL write_attribute(dset_id, 'Npart', &
                   type_ptr%Npart)
               CALL write_attribute(dset_id, 'itime', &
                   type_ptr%itime)
               CALL write_attribute(dset_id, 'active_topoid', &
                   type_ptr%active_topoid)
               CALL write_attribute(dset_id, 'h_min', &
                   type_ptr%h_min)
               CALL write_attribute(dset_id, 'h_avg', &
                   type_ptr%h_avg)
               CALL write_attribute(dset_id, 'ghostlayer', &
                   type_ptr%ghostlayer)
               CALL write_attribute(dset_id, 'time', &
                   type_ptr%time)
               CALL write_attribute(dset_id, 'flags', &
                   type_ptr%flags, ppm_param_length_partflags)
               IF (associated(type_ptr%stats)) THEN
                  pointer_addr = get_pointer(type_ptr%stats)
                  CALL store_particles_stats_d_(cpfile_id, &
                      pointer_addr, type_ptr%stats)
               ELSE
                  pointer_addr = "00000000000000000000000000000000"
               ENDIF
               CALL write_attribute(dset_id, "stats", pointer_addr, 32)


            END SUBROUTINE write_ppm_t_particles_d
