      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_map_connect_prune
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

      SUBROUTINE ppm_map_connect_prune(topoid,cd,lda,Ncon,id,Npart,lsymm,info)
      !!! This routine prunes all connections that are not needed
      !!! on this processor according to the current topology
      !!! and whether symmetry will be used or not.

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_util_invert_list
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER                , INTENT(IN   ) :: topoid
      !!! ID of current topology
      INTEGER, DIMENSION(:,:), POINTER       :: cd
      !!! Connection data
      INTEGER                , INTENT(IN   ) :: lda
      !!! Leading dimension of cd(:,:)
      INTEGER                , INTENT(INOUT) :: Ncon
      !!! Number of connections
      INTEGER, DIMENSION(:)  , TARGET        :: id
      !!! Local to global particle number mapping
      INTEGER                , INTENT(IN   ) :: Npart
      !!! Number of particles
      LOGICAL                , INTENT(IN   ) :: lsymm
      !!! Is symmetric?
      INTEGER                , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      INTEGER, DIMENSION(2) :: ldu
      INTEGER               :: i,j,k,iopt,icon,part,prunemode
      INTEGER               :: ncons_local
      INTEGER               :: min_,max_

      CHARACTER(LEN=ppm_char) :: caller='ppm_map_connect_prune'

      LOGICAL :: valid
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         CALL check
         IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Check if we have a ring topology ...
      !-------------------------------------------------------------------------
      prunemode = 0
      IF (topoid .EQ. ppm_param_topo_undefined) THEN
         prunemode = 1
      !-------------------------------------------------------------------------
      !  ... or any other valid topology
      !-------------------------------------------------------------------------
      ELSE IF (topoid .NE. ppm_param_topo_undefined) THEN
         IF (lsymm) THEN
            ! if all particles are there
            prunemode = 2
         ELSE
            ! lowest id
            prunemode = 1
         ENDIF
      ELSE
         ! no topology defined
         prunemode = 0
      ENDIF

      IF (prunemode .NE. 0) THEN
          !---------------------------------------------------------------------
          !  Allocate
          !---------------------------------------------------------------------
          iopt = ppm_param_alloc_fit
          ldu(1) = lda
          ldu(2) = Ncon
          CALL ppm_alloc(cd_local,ldu,iopt,info)
          or_fail_alloc('allocating connection list CD_LOCAL',ppm_error=ppm_error_fatal)

          cd_local(:,:) = cd(:,:)

          ncons_local = Ncon

          !---------------------------------------------------------------------
          !  Get the inverse of the local to global particle id mapping
          !---------------------------------------------------------------------
          id_temp => id(1:Npart)

          CALL ppm_util_invert_list(id_temp,id_inv,info)
          or_fail("ppm_util_invert_list",ppm_error=ppm_error_fatal)

          NULLIFY(id_temp)
      ENDIF

      min_ = LBOUND(id_inv,1)
      max_ = UBOUND(id_inv,1)

      !-------------------------------------------------------------------------
      !  prunemode 1: remove all connections where we do not have the particle
      !               with the lowest ID of the connection
      !-------------------------------------------------------------------------
      IF (prunemode .EQ. 1) THEN
         !---------------------------------------------------------------------
         !  Loop over all connections
         !---------------------------------------------------------------------
         icon = 0
         DO j = 1,ncons_local
            !------------------------------------------------------------------
            !  Get the particle with the lowest ID in this connection
            !------------------------------------------------------------------
            part = MINVAL(cd_local(:,j))

            !------------------------------------------------------------------
            !  Check if we have this particle
            !------------------------------------------------------------------
            IF ((part.GE.min_).AND.(part.LE.max_)) THEN
               IF (id_inv(part).EQ.-ppm_big_i) THEN
                  cd_local(1,j) = -1
                  icon = icon + 1
               ENDIF
            ELSE
               cd_local(1,j) = -1
               icon = icon + 1
            ENDIF
         ENDDO
      !-------------------------------------------------------------------------
      !  prunemode 2: remove all connections we do not have all particles
      !               for
      !-------------------------------------------------------------------------
      ELSE IF (prunemode .EQ. 2) THEN
         icon = 0
         DO j = 1,ncons_local
            k = 0
            !------------------------------------------------------------------
            !  Loop over all particles in the connection and count how
            !  many are missing.
            !------------------------------------------------------------------
            DO i = 1,lda
               IF ((cd_local(i,j) .GE. min_) .AND. &
               &   (cd_local(i,j) .LE. max_)) THEN
                  IF (id_inv(cd_local(i,j)).EQ.-ppm_big_i) THEN
                     k = k + 1
                  ENDIF
               ELSE
                  k = k + 1
               ENDIF
            ENDDO

            !------------------------------------------------------------------
            !  If at least one particle is missing we prune this connection
            !------------------------------------------------------------------
            IF (k .GT. 0) THEN
               cd_local(1,j) = -1
               icon = icon + 1
            ENDIF
         ENDDO
      ENDIF

      IF (prunemode .NE. 0) THEN
         !---------------------------------------------------------------------
         !  Prepare the final connection array
         !---------------------------------------------------------------------
         ldu(1) = lda
         ldu(2) = Ncon - icon
         iopt = ppm_param_alloc_fit
         CALL ppm_alloc(cd,ldu,iopt,info)
         or_fail_alloc('allocating connection array CD',ppm_error=ppm_error_fatal)

         !---------------------------------------------------------------------
         !  Go through all connections and save all that are not marked
         !---------------------------------------------------------------------
         Ncon = 0
         DO j = 1,ncons_local
            IF (cd_local(1,j) .NE. -1) THEN
               Ncon = Ncon + 1
               cd(:,Ncon) = cd_local(:,j)
            ENDIF
         ENDDO

         iopt = ppm_param_dealloc
         CALL ppm_alloc(cd_local,ldu,iopt,info)
         or_fail_dealloc('deallocating connection list CD_LOCAL', &
         & ppm_error=ppm_error_fatal)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
        IF (Ncon .LT. 0) THEN
           fail('Number of connections must be >=0', &
           & exit_point=8888,ppm_error=ppm_error_fatal)
        ENDIF
        IF (topoid .NE. ppm_param_topo_undefined) THEN
           CALL ppm_check_topoid(topoid,valid,info)
           IF (.NOT. valid) THEN
              fail('topoid out of range',exit_point=8888)
           ENDIF
        ENDIF
      8888 CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_map_connect_prune
