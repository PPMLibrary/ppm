      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_tree_divcheck
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

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_tree_divcheck_s(min_box,max_box,nbox,minboxsize,   &
     &    fixed,boxcost,ndiv,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_tree_divcheck_d(min_box,max_box,nbox,minboxsize,   &
     &    fixed,boxcost,ndiv,info)
#endif
      !!! This routine checks how many dimensions of a box are divisible.
      !!!
      !!! NOTE: A box with cost 0 is counted as non-divisible.
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
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
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: min_box
      !!! Lower coordinates of the boxes.
      !!!
      !!! 1st index: x,y[,z]                                                   +
      !!! 2nd: box ID
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: max_box
      !!! Upper coordinates of the boxes.
      !!!
      !!! 1st index: x,y[,z]                                                   +
      !!! 2nd: box ID
      INTEGER                 , INTENT(IN   ) :: nbox
      !!! Number of boxes to check
      REAL(MK), DIMENSION(:  ), INTENT(IN   ) :: minboxsize
      !!! Minimum box size in all dimensions.
      REAL(MK), DIMENSION(:  ), INTENT(IN   ) :: boxcost
      !!! Costs associated with boxes
      LOGICAL , DIMENSION(:  ), INTENT(IN   ) :: fixed
      !!! Flags telling which dimensions are fixed (i.e. must not be divided).
      INTEGER , DIMENSION(:  ), POINTER       :: ndiv
      !!! Number of divisible directions of each box (1..nbox).
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                                :: t0,lmyeps,boxlen
      REAL(MK), DIMENSION(ppm_dim)            :: ms2
      INTEGER                                 :: iopt,i,j
      INTEGER, DIMENSION(2)                   :: ldc
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_tree_divcheck',t0,info)
#if   __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif

      !-------------------------------------------------------------------------
      !  Check input arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF 

      !-------------------------------------------------------------------------
      !  If there are no boxes to check, we quit
      !-------------------------------------------------------------------------
      IF (nbox .LT. 1) THEN
          CALL ppm_write(ppm_rank,'ppm_tree_divcheck',   &
     &        'No boxes to be checked. Exiting.',info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Count and determine divisible dimensions
      !-------------------------------------------------------------------------
      DO i=1,ppm_dim
          ms2(i) = 2.0_MK*minboxsize(i)
      ENDDO

      DO i=1,nbox
          ndiv(i) = 0
          IF (boxcost(i) .GT. lmyeps) THEN
              DO j=1,ppm_dim
                  boxlen = max_box(j,i)-min_box(j,i)
                  IF (((boxlen-ms2(j)).GT.lmyeps*boxlen).AND.(.NOT.fixed(j))) THEN
                      ndiv(i) = ndiv(i) + 1
                  ENDIF
              ENDDO
          ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_tree_divcheck',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
         IF (nbox .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_tree_divcheck',     &
     &          'Number of boxes must be >= 0',__LINE__,info)
            GOTO 8888
         ENDIF
         DO i=1,ppm_dim
            IF (minboxsize(i) .LT. 0.0_MK) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_tree_divcheck',     &
     &             'the minimum box size must be > 0 !',__LINE__,info)
               GOTO 8888
            ENDIF
            DO j=1,nbox
               IF (min_box(i,j) .GT. max_box(i,j)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_tree_divcheck',   &
     &                'min_box must be <= max_box !',__LINE__,info)
                  GOTO 8888
               ENDIF
            ENDDO
         ENDDO
 8888    CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_tree_divcheck_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_tree_divcheck_d
#endif
