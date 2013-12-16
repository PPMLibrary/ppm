      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_tree_done
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
      SUBROUTINE ppm_tree_done_s(minboxes,nsubs,boxcost,iboxlist,nboxlist, &
      &          nlevel,maxvariance,maxboxcost,maxlevels,lcontinue,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_tree_done_d(minboxes,nsubs,boxcost,iboxlist,nboxlist, &
      &          nlevel,maxvariance,maxboxcost,maxlevels,lcontinue,info)
#endif
      !!! This routine decides if a decomposition is done or
      !!! one or more boxes need further refinement.
      !!!
      !!! [NOTE]
      !!! Decomposition is considered done if more boxes than
      !!! minboxes are present and the variance of the costs
      !!! of all divisible (i.e. childless and non-empty)
      !!! boxes is below maxvariance. Alternatively, a
      !!! decomposition is always done when no more divisible
      !!! boxes exist or the specified maximum number of levels
      !!! has been reached.
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
      INTEGER                 , INTENT(IN   ) :: minboxes
      !!! Minimum number of childless boxes of non-zero cost to create.
      INTEGER                 , INTENT(IN   ) :: nsubs
      !!! Number of childless boxes with non-zero cost. This is not the
      !!! same as number of divisible boxes as they need not be larger than
      !!! 2*minboxsize in this case.
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: boxcost
      !!! Computational cost associated with each box
      INTEGER , DIMENSION(:)  , INTENT(IN   ) :: iboxlist
      !!! List of divisible boxes.
      INTEGER                 , INTENT(IN   ) :: nboxlist
      !!! Number of boxes which could be divided further (i.e. length of
      !!! iboxlist).
      INTEGER                 , INTENT(IN   ) :: nlevel
      !!! Number if tree levels (tree depth) so far.
      REAL(MK)                , INTENT(IN   ) :: maxvariance
      !!! Maximum variance of cost allowed   between boxes. Set to .LE. 0 to
      !!! disable this.
      REAL(MK)                , INTENT(IN   ) :: maxboxcost
      !!! Maximum allowed cost of a box. Tree will stop if all boxes are
      !!! below this cost. If .LE. 0, cost is unlimited.
      INTEGER                 , INTENT(IN   ) :: maxlevels
      !!! Maximum number of levels to be created. Tree stops as soon as
      !!! this is reached. If .LE. 0, levels are unlimited.
      LOGICAL                 , INTENT(  OUT) :: lcontinue
      !!! `FALSE` if no further subdivision is needed, `TRUE` otherwise.
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK) :: t0,meancost,diffcost,varcost
      REAL(MK) :: maxcost,dm

      INTEGER :: i,j

      CHARACTER(LEN=ppm_char) :: caller='ppm_tree_done'
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)
      lcontinue = .TRUE.

      !-------------------------------------------------------------------------
      !  Check input arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         CALL check
         IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  If there are no boxes, we quit
      !-------------------------------------------------------------------------
      IF (nsubs .LT. 1) THEN
         lcontinue = .FALSE. ! no boxes = nothing to subdivide !
         CALL ppm_write(ppm_rank,caller,   &
         &    'No non-empty boxes present. Done.',info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  If there are no more boxes that could potentially be divided, we
      !  have to stop.
      !-------------------------------------------------------------------------
      IF (nboxlist .LT. 1) THEN
         lcontinue = .FALSE.
         !---------------------------------------------------------------------
         !  If there are less boxes than processors, we have a problem
         !  THIS SHOULD GO TO THE DECOMP ROUTINE AND NOT INTO THE GENERIC
         !  TREE !!!
         !---------------------------------------------------------------------
         IF (nsubs .LT. minboxes) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_few_subs,caller, &
            &    'Could not create the minimum number of non-empty boxes!', &
            &    __LINE__,info)
         ENDIF
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  If the max tree depth has been reached, we stop
      !-------------------------------------------------------------------------
      IF (maxlevels .GT. 0) THEN
          IF (nlevel .GE. maxlevels) THEN
              lcontinue = .FALSE.
              IF (ppm_debug .GT. 1) THEN
                  CALL ppm_write(ppm_rank,caller,     &
     &                'Max number of levels reached. Done.',info)
              ENDIF
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  If there are less boxes than processors, decomposition is not OK
      !-------------------------------------------------------------------------
      IF (nsubs .LT. minboxes) GOTO 9999

      !-------------------------------------------------------------------------
      !  Compute the variance of divisible box costs
      !-------------------------------------------------------------------------
      varcost  = 0.0_MK
      meancost = 0.0_MK
      maxcost  = -HUGE(maxcost)
      DO i=1,nboxlist
          j = iboxlist(i)
          dm = boxcost(j)
          meancost = meancost + dm
          IF (dm .GT. maxcost) THEN
              maxcost = dm
          ENDIF
      ENDDO
      IF (nboxlist .GT. 1) THEN
          meancost = meancost/REAL(nboxlist,MK)
          DO i=1,nboxlist
              j = iboxlist(i)
              diffcost = boxcost(j) - meancost
              varcost  = varcost + (diffcost*diffcost)
          ENDDO
          varcost  = varcost/REAL(nboxlist-1,MK)
      ENDIF

      !-------------------------------------------------------------------------
      !  If variance of costs is below threshold, decomposition is OK
      !-------------------------------------------------------------------------
      IF (varcost .LT. maxvariance) lcontinue = .FALSE.

      !-------------------------------------------------------------------------
      !  If all boxes are below the maxcost we are done.
      !-------------------------------------------------------------------------
      IF (maxcost .LT. maxboxcost) lcontinue = .FALSE.

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
         IF (nsubs .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,     &
     &          'Number of non-empty boxes must be >= 0',__LINE__,info)
            GOTO 8888
         ENDIF
         IF (nlevel .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,     &
     &          'Number of levels must be >= 0',__LINE__,info)
            GOTO 8888
         ENDIF
         IF (nboxlist .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,     &
     &          'Number of boxes in list must be >= 0',__LINE__,info)
            GOTO 8888
         ENDIF
 8888    CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_tree_done_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_tree_done_d
#endif
