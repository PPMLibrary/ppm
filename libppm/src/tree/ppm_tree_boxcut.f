      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_tree_boxcut
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Lab (ETH Zurich), 
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

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_tree_boxcut_s(xp,cutbox,min_box,max_box,ncut,cutdir,  &
     &    cutpos,mincut,maxcut,lhbx,lpdx,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_tree_boxcut_d(xp,cutbox,min_box,max_box,ncut,cutdir,  &
     &    cutpos,mincut,maxcut,lhbx,lpdx,info)
#endif
      !!! This routine cuts a box into 2, 4 or 8 pieces using
      !!! the given cut directions and cut positions.
      !!!
      !!! [NOTE]
      !!! The routine only works for ncut=1, 2, or 3.
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_tree
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xp
      !!! Positions of the points.
      REAL(MK), DIMENSION(:  ), INTENT(IN   ) :: min_box
      !!! Minimum coordinates of the box to be cut.
      REAL(MK), DIMENSION(:  ), INTENT(IN   ) :: max_box
      !!! Maximum coordinates of the box to be cut.
      REAL(MK), DIMENSION(:  ), INTENT(IN   ) :: cutpos
      !!! Positions of the cuts along the given axes.
      INTEGER                 , INTENT(IN   ) :: ncut
      !!! Number of cuts to apply
      INTEGER                 , INTENT(IN   ) :: cutbox
      !!! ID of the box to be cut.
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: cutdir
      !!! Directions of the cuts. cutdir(1)=2 means that the first
      !!! cut is perpendicular to the second coordinate axis (y axis).
      REAL(MK), DIMENSION(:,:), POINTER       :: mincut
      !!! Minimum coordinates of the cut boxes.
      !!!
      !!! 1st: x,y[,z]                                                         +
      !!! 2nd: ibox.
      REAL(MK), DIMENSION(:,:), POINTER       :: maxcut
      !!! Maximum coordinates of the cut boxes.
      !!!
      !!! 1st: x,y[,z]                                                         +
      !!! 2nd: ibox.
      INTEGER , DIMENSION(:  ), POINTER       :: lhbx
      !!! Pointer to the first point (in lpdx) in each sub-box. This is
      !!! only allocated and returned if there are particles at all.
      INTEGER , DIMENSION(:  ), POINTER       :: lpdx
      !!! Index of particles in sub-boxes. This is only allocated and returned
      !!! if there are particles at all.
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                                :: t0
      INTEGER                                 :: nnbox,cd,i,j,iopt,Np,ip
      INTEGER , DIMENSION(2)                  :: ldc
      INTEGER , DIMENSION(:), POINTER         :: pbox
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_tree_boxcut',t0,info)

      !-------------------------------------------------------------------------
      !  Check input arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Exit if there is nothing to be cut
      !-------------------------------------------------------------------------
      IF (ncut .LT. 1) THEN
          IF (ppm_debug .GT. 0) THEN
              CALL ppm_write(ppm_rank,'ppm_tree_boxcut',    &
     &            'Nothing to cut. Exiting.',info)
          ENDIF
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Number of new boxes to be created
      !-------------------------------------------------------------------------
      nnbox = 2**ncut
      iopt  = ppm_param_alloc_grow
      IF (have_particles) THEN
          Np = - tree_lhbx(1,cutbox) + tree_lhbx(2,cutbox) + 1
          ldc(1) = Np
          CALL ppm_alloc(lpdx,ldc,iopt,info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_tree_boxcut',          &
     &            'particle list pointers LPDX',__LINE__,info)
              GOTO 9999
          ENDIF 
      ENDIF 

      !-------------------------------------------------------------------------
      !  Split box into 2 pieces
      !-------------------------------------------------------------------------
      IF (ncut .EQ. 1) THEN
          IF (ppm_dim .GT. 2) THEN
              DO i=1,nnbox
                  mincut(1,i) = min_box(1)
                  mincut(2,i) = min_box(2)
                  mincut(3,i) = min_box(3)
                  maxcut(1,i) = max_box(1)
                  maxcut(2,i) = max_box(2)
                  maxcut(3,i) = max_box(3)
              ENDDO
          ELSE
              DO i=1,nnbox
                  mincut(1,i) = min_box(1)
                  mincut(2,i) = min_box(2)
                  maxcut(1,i) = max_box(1)
                  maxcut(2,i) = max_box(2)
              ENDDO
          ENDIF

          cd  = cutdir(1)
          mincut(cd,1) = min_box(cd)
          mincut(cd,2) = cutpos(1)

          maxcut(cd,1) = cutpos(1)
          maxcut(cd,2) = max_box(cd)
          
      !-------------------------------------------------------------------------
      !  Split box into 4 pieces
      !-------------------------------------------------------------------------
      ELSEIF (ncut .EQ. 2) THEN
          IF (ppm_dim .GT. 2) THEN
              DO i=1,nnbox
                  mincut(1,i) = min_box(1)
                  mincut(2,i) = min_box(2)
                  mincut(3,i) = min_box(3)
                  maxcut(1,i) = max_box(1)
                  maxcut(2,i) = max_box(2)
                  maxcut(3,i) = max_box(3)
              ENDDO
          ELSE
              DO i=1,nnbox
                  mincut(1,i) = min_box(1)
                  mincut(2,i) = min_box(2)
                  maxcut(1,i) = max_box(1)
                  maxcut(2,i) = max_box(2)
              ENDDO
          ENDIF

          cd = cutdir(1)
          mincut(cd,1) = min_box(cd)
          mincut(cd,2) = min_box(cd)
          mincut(cd,3) = cutpos(1)
          mincut(cd,4) = cutpos(1)

          maxcut(cd,1) = cutpos(1)
          maxcut(cd,2) = cutpos(1)
          maxcut(cd,3) = max_box(cd)
          maxcut(cd,4) = max_box(cd)

          cd = cutdir(2)
          mincut(cd,1) = min_box(cd)
          mincut(cd,2) = cutpos(2)
          mincut(cd,3) = min_box(cd)
          mincut(cd,4) = cutpos(2)

          maxcut(cd,1) = cutpos(2)
          maxcut(cd,2) = max_box(cd)
          maxcut(cd,3) = cutpos(2)
          maxcut(cd,4) = max_box(cd)
      
      !-------------------------------------------------------------------------
      !  Split box into 8 pieces
      !-------------------------------------------------------------------------
      ELSEIF (ncut .EQ. 3) THEN
          IF (ppm_dim .GT. 2) THEN
              DO i=1,nnbox
                  mincut(1,i) = min_box(1)
                  mincut(2,i) = min_box(2)
                  mincut(3,i) = min_box(3)
                  maxcut(1,i) = max_box(1)
                  maxcut(2,i) = max_box(2)
                  maxcut(3,i) = max_box(3)
              ENDDO
          ELSE
              DO i=1,nnbox
                  mincut(1,i) = min_box(1)
                  mincut(2,i) = min_box(2)
                  maxcut(1,i) = max_box(1)
                  maxcut(2,i) = max_box(2)
              ENDDO
          ENDIF

          cd = cutdir(1)
          mincut(cd,1) = min_box(cd)
          mincut(cd,2) = min_box(cd)
          mincut(cd,3) = cutpos(1)
          mincut(cd,4) = cutpos(1)

          maxcut(cd,1) = cutpos(1)
          maxcut(cd,2) = cutpos(1)
          maxcut(cd,3) = max_box(cd)
          maxcut(cd,4) = max_box(cd)

          mincut(cd,5) = min_box(cd)
          mincut(cd,6) = min_box(cd)
          mincut(cd,7) = cutpos(1)
          mincut(cd,8) = cutpos(1)

          maxcut(cd,5) = cutpos(1)
          maxcut(cd,6) = cutpos(1)
          maxcut(cd,7) = max_box(cd)
          maxcut(cd,8) = max_box(cd)

          cd = cutdir(2)
          mincut(cd,1) = min_box(cd)
          mincut(cd,2) = cutpos(2)
          mincut(cd,3) = min_box(cd)
          mincut(cd,4) = cutpos(2)

          maxcut(cd,1) = cutpos(2)
          maxcut(cd,2) = max_box(cd)
          maxcut(cd,3) = cutpos(2)
          maxcut(cd,4) = max_box(cd)

          mincut(cd,5) = min_box(cd)
          mincut(cd,6) = cutpos(2)
          mincut(cd,7) = min_box(cd)
          mincut(cd,8) = cutpos(2)

          maxcut(cd,5) = cutpos(2)
          maxcut(cd,6) = max_box(cd)
          maxcut(cd,7) = cutpos(2)
          maxcut(cd,8) = max_box(cd)

          cd = cutdir(3)
          mincut(cd,1) = min_box(cd)
          mincut(cd,2) = min_box(cd)
          mincut(cd,3) = min_box(cd)
          mincut(cd,4) = min_box(cd)

          maxcut(cd,1) = cutpos(3)
          maxcut(cd,2) = cutpos(3)
          maxcut(cd,3) = cutpos(3)
          maxcut(cd,4) = cutpos(3)

          mincut(cd,5) = cutpos(3)
          mincut(cd,6) = cutpos(3)
          mincut(cd,7) = cutpos(3)
          mincut(cd,8) = cutpos(3)

          maxcut(cd,5) = max_box(cd)
          maxcut(cd,6) = max_box(cd)
          maxcut(cd,7) = max_box(cd)
          maxcut(cd,8) = max_box(cd)
      
      !-------------------------------------------------------------------------
      !  Unknown cutting mode
      !-------------------------------------------------------------------------
      ELSE
          info = ppm_error_error
          CALL ppm_error(ppm_err_argument,'ppm_tree_boxcut',     &
     &        'ncut must be < 4.',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Rank particles into new boxes. We cannot use ppm_util_rank
      !  here, since the child boxes need not be all of the same size!
      !-------------------------------------------------------------------------
      IF (have_particles) THEN
          !---------------------------------------------------------------------
          !  Allocate particle box ID list
          !---------------------------------------------------------------------
          iopt = ppm_param_alloc_fit
          ldc(1) = Np
          CALL ppm_alloc(pbox,ldc,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_tree_boxcut',   &
     &            'list of particle box IDs PBOX',__LINE__,info)
              GOTO 9999
          ENDIF
          pbox = -1
          npbx = 0

          !---------------------------------------------------------------------
          !  Determine sub-box for each particle. TODO: use the snipplet
          !  here... but avoid alloc!
          !---------------------------------------------------------------------
          IF (ppm_dim .EQ. 3) THEN
              DO i=1,nnbox
                  DO j=1,Np
                      ip = tree_lpdx(tree_lhbx(1,cutbox)+j-1)
                      IF (.NOT.(xp(1,ip) .LT. mincut(1,i)) .AND.    &
     &                    (xp(1,ip) .LT. maxcut(1,i)) .AND.    &
     &                    .NOT.(xp(2,ip) .LT. mincut(2,i)) .AND.    &
     &                    (xp(2,ip) .LT. maxcut(2,i)) .AND.    &
     &                    .NOT.(xp(3,ip) .LT. mincut(3,i)) .AND.    &
     &                    (xp(3,ip) .LT. maxcut(3,i))) THEN
                          pbox(j) = i
                      ENDIF
                  ENDDO
              ENDDO
          ELSE
              DO i=1,nnbox
                  DO j=1,Np
                      ip = tree_lpdx(tree_lhbx(1,cutbox)+j-1)
                      IF (.NOT.(xp(1,ip) .LT. mincut(1,i)) .AND.    &
     &                    (xp(1,ip) .LT. maxcut(1,i)) .AND.    &
     &                    .NOT.(xp(2,ip) .LT. mincut(2,i)) .AND.    &
     &                    (xp(2,ip) .LT. maxcut(2,i))) THEN
                          pbox(j) = i
                      ENDIF
                  ENDDO
              ENDDO
          ENDIF

          !---------------------------------------------------------------------
          !  Determine number of particles in each sub-box
          !---------------------------------------------------------------------
          DO j=1,Np
              i = pbox(j)
              IF (i .GT. 0) npbx(i) = npbx(i) + 1
          ENDDO

          !---------------------------------------------------------------------
          !  Build the box particle pointers
          !---------------------------------------------------------------------
          cbox(1) = 1
          lhbx(1) = 1
          DO i=2,nnbox
              cbox(i) = cbox(i-1) + npbx(i-1)
              lhbx(i) = cbox(i)
          ENDDO
          lhbx(nnbox+1) = lhbx(nnbox) + npbx(nnbox)

          !---------------------------------------------------------------------
          !  Rank the particles in the sub-boxes
          !---------------------------------------------------------------------
          DO j=1,Np
              ip = tree_lpdx(tree_lhbx(1,cutbox)+j-1)
              i = pbox(j)
              IF (i .GT. 0) THEN
                  lpdx(cbox(i)) = ip
                  cbox(i)       = cbox(i) + 1
              ENDIF
          ENDDO

          !---------------------------------------------------------------------
          !  Deallocate memory
          !---------------------------------------------------------------------
          iopt = ppm_param_dealloc
          CALL ppm_alloc(pbox,ldc,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_tree_boxcut',     &
     &            'list of particle box IDs PBOX',__LINE__,info)
          ENDIF
      ENDIF           ! have_particles

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_tree_boxcut',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
         IF (cutbox .LT. 1) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_tree_boxcut',     &
     &          'cutbox must be >0 !',__LINE__,info)
            GOTO 8888
         ENDIF
         IF ((ppm_dim .LT. 3) .AND. (ncut .GT. 2)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_tree_boxcut',     &
     &          'Cannot cut more than 2 times in 2D!',__LINE__,info)
            GOTO 8888
         ENDIF
         DO i=1,ppm_dim
            IF (min_box(i) .GT. max_box(i)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_tree_boxcut',     &
     &             'min_box must be <= max_box !',__LINE__,info)
               GOTO 8888
            ENDIF
         ENDDO
         DO i=1,ncut-1
            DO j=i+1,ncut
               IF (cutdir(i) .EQ. cutdir(j)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_tree_boxcut',     &
     &                'cut directions must be distinct !',__LINE__,info)
                  GOTO 8888
               ENDIF
            ENDDO
         ENDDO
 8888    CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_tree_boxcut_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_tree_boxcut_d
#endif
