      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_alloc_1d
      !-------------------------------------------------------------------------
      ! Copyright (c) 2012 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      !                    Center for Fluid Dynamics (DTU)
      !
      !
      ! This file is part of the Parallel Particle Mesh Library (PPM).
      !
      ! PPM is free software: you can redistribute it and/or modify
      ! it under the terms of the GNU Lesser General Public License
      ! as published by the Free Software Foundation, either
      ! version 3 of the License, or (at your option) any later
      ! version.
      !
      ! PPM is distributed in the hope that it will be useful,
      ! but WITHOUT ANY WARRANTY; without even the implied warranty of
      ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      ! GNU General Public License for more details.
      !
      ! You should have received a copy of the GNU General Public License
      ! and the GNU Lesser General Public License along with PPM. If not,
      ! see <http://www.gnu.org/licenses/>.
      !
      ! Parallel Particle Mesh Library (PPM)
      ! ETH Zurich
      ! CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------


      SUBROUTINE ppm_mesh_alloc(mesh,iopt,info)
      !!! (Re)allocates the memory of one-dimensional arrays
      !!! (pointers) based on the number of elements.
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_interface

      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_t_equi_mesh_), POINTER       :: mesh
      INTEGER,                 INTENT(IN   ) :: iopt
      !!! Allocation mode. One of:
      !!!
      !!! * ppm_param_alloc_fit
      !!! * ppm_param_dealloc
      INTEGER,                 INTENT(  OUT) :: info
      !!! Returns status, 0 upon success.
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      CHARACTER(LEN=ppm_char) :: caller='ppm_mesh_alloc'

      LOGICAL :: lalloc,ldealloc

      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      ! maybe add some sanity checks later

      !-------------------------------------------------------------------------
      !  Check the allocation type
      !-------------------------------------------------------------------------
      lalloc   = .FALSE.
      ldealloc = .FALSE.
      SELECT CASE (iopt)
      CASE (ppm_param_alloc_fit)
         !----------------------------------------------------------------------
         !  fit memory but skip the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(mesh)) ldealloc = .TRUE.
         lalloc   = .TRUE.

      CASE (ppm_param_dealloc)
         ldealloc = .TRUE.

      CASE DEFAULT
         !----------------------------------------------------------------------
         !  Unknown iopt
         !----------------------------------------------------------------------
         info = ppm_error_error
         CALL ppm_error(ppm_err_argument,caller, &
         &    'unknown iopt',__LINE__,info)
         GOTO 9999
      END SELECT

      !-------------------------------------------------------------------------
      !  If reallocating, deallocate old data first
      !-------------------------------------------------------------------------
      IF (ldealloc) THEN
         !----------------------------------------------------------------------
         !  deallocate
         !----------------------------------------------------------------------
         IF (ASSOCIATED(mesh)) THEN
            ! first deallocate all content of mesh
            IF (ASSOCIATED(mesh%Nm))                 DEALLOCATE(mesh%Nm,                 STAT=info)
            IF (ASSOCIATED(mesh%Offset))             DEALLOCATE(mesh%Offset,             STAT=info)
            IF (ASSOCIATED(mesh%h))                  DEALLOCATE(mesh%h,                  STAT=info)
            IF (ASSOCIATED(mesh%ghostsize))          DEALLOCATE(mesh%ghostsize,          STAT=info)
            IF (ASSOCIATED(mesh%nnodes))             DEALLOCATE(mesh%nnodes,             STAT=info)
            IF (ASSOCIATED(mesh%istart))             DEALLOCATE(mesh%istart,             STAT=info)
            IF (ASSOCIATED(mesh%iend))               DEALLOCATE(mesh%iend,               STAT=info)
            IF (ASSOCIATED(mesh%subpatch))           DEALLOCATE(mesh%subpatch,           STAT=info)
            IF (ASSOCIATED(mesh%mdata))              DEALLOCATE(mesh%mdata,              STAT=info)
            IF (ASSOCIATED(mesh%patch))              DEALLOCATE(mesh%patch,              STAT=info)
            IF (ASSOCIATED(mesh%subpatch_by_sub))    DEALLOCATE(mesh%subpatch_by_sub,    STAT=info)
            IF (ASSOCIATED(mesh%field_ptr))          DEALLOCATE(mesh%field_ptr,          STAT=info)
            IF (ASSOCIATED(mesh%mapping_s))          DEALLOCATE(mesh%mapping_s,          STAT=info)
            IF (ASSOCIATED(mesh%mapping_d))          DEALLOCATE(mesh%mapping_d,          STAT=info)
            IF (ASSOCIATED(mesh%ghost_fromsub))      DEALLOCATE(mesh%ghost_fromsub,      STAT=info)
            IF (ASSOCIATED(mesh%ghost_tosub))        DEALLOCATE(mesh%ghost_tosub,        STAT=info)
            IF (ASSOCIATED(mesh%ghost_patchid))      DEALLOCATE(mesh%ghost_patchid,      STAT=info)
            IF (ASSOCIATED(mesh%ghost_blkstart))     DEALLOCATE(mesh%ghost_blkstart,     STAT=info)
            IF (ASSOCIATED(mesh%ghost_blksize))      DEALLOCATE(mesh%ghost_blksize,      STAT=info)
            IF (ASSOCIATED(mesh%ghost_blk))          DEALLOCATE(mesh%ghost_blk,          STAT=info)
            IF (ASSOCIATED(mesh%ghost_recvtosub))    DEALLOCATE(mesh%ghost_recvtosub,    STAT=info)
            IF (ASSOCIATED(mesh%ghost_recvpatchid))  DEALLOCATE(mesh%ghost_recvpatchid,  STAT=info)
            IF (ASSOCIATED(mesh%ghost_recvblkstart)) DEALLOCATE(mesh%ghost_recvblkstart, STAT=info)
            IF (ASSOCIATED(mesh%ghost_recvblksize))  DEALLOCATE(mesh%ghost_recvblksize,  STAT=info)
            IF (ASSOCIATED(mesh%ghost_recvblk))      DEALLOCATE(mesh%ghost_recvblk,      STAT=info)
            IF (ASSOCIATED(mesh%mapping))            DEALLOCATE(mesh%mapping,            STAT=info)

            DEALLOCATE(mesh, STAT=info)
            NULLIFY(mesh)
            or_fail('Deallocating mesh')
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate new memory
      !-------------------------------------------------------------------------
      IF (lalloc) THEN
         ALLOCATE(mesh, STAT=info)
         or_fail('Allocating mesh')

         NULLIFY(mesh%Nm)
         NULLIFY(mesh%Offset)
         NULLIFY(mesh%h)
         NULLIFY(mesh%ghostsize)
         NULLIFY(mesh%nnodes)
         NULLIFY(mesh%istart)
         NULLIFY(mesh%iend)
         NULLIFY(mesh%subpatch)
         NULLIFY(mesh%mdata)
         NULLIFY(mesh%patch)
         NULLIFY(mesh%subpatch_by_sub)
         NULLIFY(mesh%field_ptr)
         NULLIFY(mesh%mapping_s)
         NULLIFY(mesh%mapping_d)
         NULLIFY(mesh%ghost_fromsub)
         NULLIFY(mesh%ghost_tosub)
         NULLIFY(mesh%ghost_patchid)
         NULLIFY(mesh%ghost_blkstart)
         NULLIFY(mesh%ghost_blksize)
         NULLIFY(mesh%ghost_blk)
         NULLIFY(mesh%ghost_recvtosub)
         NULLIFY(mesh%ghost_recvpatchid)
         NULLIFY(mesh%ghost_recvblkstart)
         NULLIFY(mesh%ghost_recvblksize)
         NULLIFY(mesh%ghost_recvblk)
         NULLIFY(mesh%mapping)

         mesh%ID = 0
         mesh%name = "default_mesh_name"
         mesh%topoid = 0
         mesh%npatch = 0
         mesh%ghost_initialized = .FALSE.
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE ppm_mesh_alloc

