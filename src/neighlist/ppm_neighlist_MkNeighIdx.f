      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_neighlist_MkNeighIdx
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

      SUBROUTINE ppm_neighlist_MkNeighIdx(lsymm,ind,jnd,nnd,info)
      !!! Creates the index offset list of cell interactions.
      !!! The interaction for the cell with itself is always included as
      !!! the first entry.
      !!!
      !!! [NOTE]
      !!! If the loops do not vectorize, maybe we need to
      !!! duplicate them and put IF(ppm_dim...) statements
      !!! around!
      !!!
      !!! [WARNING]
      !!! `ind` and `jnd` are allocated inside this routine! The
      !!! lists are always (3 x nnd) since then no IF
      !!! statements in the inner loop are needed to
      !!! distinguish between 2D and 3D case. Since `nnd` is
      !!! <28 anyway, this should not be too much of a memory
      !!! waste.

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
      USE ppm_module_alloc
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      LOGICAL,                 INTENT(IN   ) :: lsymm
      !!! T for using symmetry, F for full list
      INTEGER, DIMENSION(:,:), POINTER       :: ind
      !!! First interaction partner (box which *interacts*).
      !!!
      !!! 1st index: 1...3 (x,y,[z]) index shift.                              +
      !!! 2nd index: interaction number 1...nnd.
      INTEGER, DIMENSION(:,:), POINTER       :: jnd
      !!! Second interaction partner (box which *is interacted with*).
      !!!
      !!! 1st index: 1...3 (x,y,[z]) index shift.                              +
      !!! 2nd index: interaction number 1...nnd.
      INTEGER,                 INTENT(  OUT) :: nnd
      !!! Number of box-box interactions to be performed.
      INTEGER,                 INTENT(  OUT) :: info
      !!! Returns status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      ! alloc
      INTEGER, DIMENSION(2) :: lda
      INTEGER               :: iopt
      INTEGER               :: i,j,k,l,nz

      CHARACTER(LEN=ppm_char) :: caller='ppm_neighlist_MkNeighIdx'

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         IF (.NOT. ppm_initialized) THEN
            fail('Please call ppm_init first!',ppm_err_ppm_noinit)
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Determine number of box-box interactions needed
      !-------------------------------------------------------------------------
      SELECT CASE (ppm_dim)
      CASE (2)
         IF (lsymm) THEN
            ! 2D using symmetry
            nnd = 5
         ELSE
            ! 2D NOT using symmetry
            nnd = 9
         ENDIF
      CASE (3)
         IF (lsymm) THEN
            ! 3D using symmetry
            nnd = 14
         ELSE
            ! 3D NOT using symmetry
            nnd = 27
         ENDIF
      END SELECT

      !-------------------------------------------------------------------------
      !  Allocate memory for interaction lists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      lda(1) = 3
      lda(2) = nnd

      CALL ppm_alloc(ind,lda,iopt,info)
      or_fail_alloc('Interaction list IND',ppm_error=ppm_error_fatal)

      CALL ppm_alloc(jnd,lda,iopt,info)
      or_fail_alloc('Interaction list JND',ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Initialize ind and jnd
      !-------------------------------------------------------------------------
      ind = 0
      jnd = 0

      !---------------------------------------------------------------------
      !  Compute neighbour indices
      !---------------------------------------------------------------------
      IF (lsymm) THEN
         !------------------------------------------------------------------
         !  Using symmetry
         !------------------------------------------------------------------
                          ! interaction   0 -- 0 (as initialized)

         jnd(1,2) = 1     ! interaction   0 -- 1

         ind(2,3) = 1     ! interaction   0 -- 3

         jnd(1,4) = 1     ! interaction   0 -- 4
         jnd(2,4) = 1

         ind(1,5) = 1     ! interaction   1 -- 3
         jnd(2,5) = 1

         IF (ppm_dim .EQ. 3) THEN
             jnd(3,5)   = 1   ! interaction   0 -- 9
             ind(1,5)   = 0   ! reset to zero (1-3 will be further down)
             jnd(2,5)   = 0

             jnd(1,6)   = 1   ! interaction   0 -- 10
             jnd(3,6)   = 1

             jnd(2,7)   = 1   ! interaction   0 -- 12
             jnd(3,7)   = 1

             jnd(1,8)   = 1   ! interaction   0 -- 13
             jnd(2,8)   = 1
             jnd(3,8)   = 1

             ind(1,9)   = 1   ! interaction   1 -- 3
             jnd(2,9)   = 1

             ind(1,10)  = 1   ! interaction   1 -- 9
             jnd(3,10)  = 1

             ind(1,11)  = 1   ! interaction   1 -- 12
             jnd(2,11)  = 1
             jnd(3,11)  = 1

             ind(2,12)  = 1   ! interaction   3 -- 9
             jnd(3,12)  = 1

             ind(2,13)  = 1   ! interaction   3 -- 10
             jnd(1,13)  = 1
             jnd(3,13)  = 1

             ind(1,14)  = 1   ! interaction   4 -- 9
             ind(2,14)  = 1
             jnd(3,14)  = 1
         ENDIF
      ELSE
         !-------------------------------------------------------------------------
         !  Set z direction according to dimensionality
         !-------------------------------------------------------------------------
         nz = MERGE(1,0,ppm_dim.EQ.3)

         !------------------------------------------------------------------
         !  Full list
         !------------------------------------------------------------------
         l = 0
         DO k=-nz,nz
            DO j=-1,1
               DO i=-1,1
                  !---------------------------------------------------------
                  !  Add to the interaction lists
                  !---------------------------------------------------------
                  l = l + 1
                  ! center box interacts WITH all boxes around it
                  jnd(1,l) = i
                  jnd(2,l) = j
                  jnd(3,l) = k     ! will always be 0 for 2d case
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE ppm_neighlist_MkNeighIdx
