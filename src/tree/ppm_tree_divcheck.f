      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_tree_divcheck
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich), 
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
     &    fixed,boxcost,neigh_constraints,num_constr,ndiv,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_tree_divcheck_d(min_box,max_box,nbox,minboxsize,   &
     &    fixed,boxcost,neigh_constraints,num_constr,ndiv,info)
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
      USE ppm_module_alloc
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
      REAL(MK), DIMENSION(:,:,:,:), INTENT(IN   ) :: neigh_constraints
      ! access: boxid, dimension, constraintid, 1 (from) 2 (to)
      INTEGER, DIMENSION(:,:), INTENT(IN )      :: num_constr
      ! number of constraints in a box and dimension
      INTEGER , DIMENSION(:  ), POINTER       :: ndiv
      !!! Number of divisible directions of each box (1..nbox).
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                                :: t0,lmyeps,boxlen
      REAL(MK), DIMENSION(ppm_dim)            :: ms2
      INTEGER                                 :: iopt,i,j,k,temp_r
      INTEGER, DIMENSION(2)                   :: ldc
      LOGICAL, DIMENSION(ppm_dim)             :: is_possible
      REAL(MK), DIMENSION(:,:,:), POINTER       :: result_array
      ! access: dimension, constraintid, 1 (from) 2 (to)
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
      DO i=1,nbox

!          IF (ppm_rank .EQ. 0)THEN
!             print *, ' '
!             print *, 'box: ', min_box(1,i),min_box(2,i),max_box(1,i),max_box(2,i)
!          ENDIF

         DO j=1,ppm_dim
            ms2(j) = 2.0_MK*minboxsize(j)

            !--------------------------------------------------------------------
            ! We need to check if path of non cutable regions is equal to box
            !--------------------------------------------------------------------
            is_possible(j) = .TRUE.
            
             IF(num_constr(i,j) .GT. 0) THEN
               CALL get_subarray(neigh_constraints(i,:,:,:),num_constr(i,:),result_array)
               k = 1
               IF (result_array(j,k,1)+lmyeps .LT. min_box(j,i)+minboxsize(j)) THEN
               
                  ! Find the first constraint in the minboxsize position
                  temp_r = k
                  DO WHILE (.NOT.((min_box(j,i)+minboxsize(j)-lmyeps).GT.result_array(j,k,1) &
            &                         .AND. (min_box(j,i)+minboxsize(j)+lmyeps).LT.result_array(j,k,2)))
                     IF (result_array(j,k,2) .GT. result_array(j,temp_r,2)) THEN
                        temp_r = k
                     ENDIF  
                     IF (k .EQ. num_constr(i,j)) THEN
                        EXIT
                     ENDIF
                     k = k+1
                  ENDDO
               


               ! check if we have found one
               IF (k .LE. num_constr(i,j)) THEN
                  ! we have found one
                   IF (result_array(j,k,2) .GT. result_array(j,temp_r,2)) THEN
                        temp_r = k
                     ENDIF  
                  ! start to travel until rightest position
                  k = k+1
                  IF (k .LE. num_constr(i,j)) THEN
                     DO WHILE(result_array(j,k,1)+lmyeps .LT. result_array(j,temp_r,2) &
         &                      .OR. result_array(j,k,2)+lmyeps .LT. min_box(j,i)+minboxsize(j))
                        IF (result_array(j,k,2) .GT. result_array(j,temp_r,2)) THEN
                           temp_r = k
                        ENDIF
                        IF (k .EQ. num_constr(i,j)) THEN
                           EXIT
                        ENDIF
                        k = k+1
                     ENDDO
                  ELSE
                     k = k-1
                  ENDIF 
                  
!                IF (result_array(j,k,1)+lmyeps .LT. min_box(j,i)+minboxsizes(j,i)) THEN
!                   !     until distance is >= 0
!                   IF (num_constr(i,j) .GT. 1) THEN
!                      k = k+1
!                      DO WHILE(result_array(j,k,1)+lmyeps .LT. result_array(j,k-1,2) &
!          &                      .OR. result_array(j,k,1)+lmyeps .LT. min_box(j,i)+minboxsizes(j,i))
!                         temp_r = k-1
!                         DO WHILE(result_array(j,k,2)+lmyeps .LT. result_array(j,temp_r,2))
!                            k = k+1
!                            IF (k .GT. num_constr(i,j)) THEN
!                               k = k-1
!                               EXIT
!                            ENDIF
!                         ENDDO
!                         IF (result_array(j,k,1)+lmyeps .GT. result_array(j,temp_r,2) .AND. &
!          &                      .NOT. (result_array(j,k,1)-lmyeps .LT. min_box(j,i)+minboxsizes(j,i))) THEN
!                            ! temp_right is the first possible
!                            k = temp_r+1
!                            EXIT
!                         ENDIF
!                         k = k+1
!                         IF (k .GT. num_constr(i,j)) THEN
!                            IF (result_array(j,k-1,2) .LT. result_array(j,temp_r,2) .AND. &
!             &                      .NOT. (result_array(j,k,1)-lmyeps .LT. min_box(j,i)+minboxsizes(j,i))) THEN
!                               ! temp_right is the first possible
!                               k = temp_r+1
!                            ENDIF
!                            EXIT
!                         ENDIF
!                      ENDDO
!                      k = k-1
!                   ENDIF
                  ! if it is inside max - ghostsize -> ok
                  IF (result_array(j,temp_r,2)-lmyeps .GT. max_box(j,i)-minboxsize(j)) THEN
                     is_possible(j) = .FALSE.
!                      IF (ppm_rank .EQ. 0)THEN
!                          print *, ' CASE IN divcheck ' , i
!                       ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDIF



         ENDDO
      
        

          ndiv(i) = 0
          IF (boxcost(i) .GT. lmyeps) THEN
              DO j=1,ppm_dim
                  boxlen = max_box(j,i)-min_box(j,i)
                  IF (((boxlen-ms2(j)).GT.lmyeps*boxlen).AND.(.NOT.fixed(j)) .AND. is_possible(j)) THEN
                      ndiv(i) = ndiv(i) + 1
                  ENDIF
!                IF (ppm_rank .EQ. 0)THEN
!                   boxlen = max_box(j,i)-min_box(j,i)
!                   print *, j, ndiv(i), is_possible(j), (((boxlen-ms2(j)).GT.lmyeps*boxlen).AND.(.NOT.fixed(j)))
!                ENDIF
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

      SUBROUTINE get_subarray(neigh_ghost,n,res)
         ! This subroutine is used to get the sublist for the neighboring constraints
         INTEGER                                    :: iopt, info, i, j, k
         INTEGER , DIMENSION(3)                     :: ldc
         INTEGER, DIMENSION(ppm_dim), INTENT(IN  )  :: n
         REAL(MK), DIMENSION(:,:,:), INTENT(IN   )  :: neigh_ghost
         REAL(MK), DIMENSION(:,:,:), POINTER        :: res
         REAL(MK)                                   :: temp1, temp2

         ! Res is first deallocated then allocated here
         iopt = ppm_param_dealloc
         ldc(1) = ppm_dim
         ldc(2) = n(1)
         ! take maximum of length of constraints
         DO i = 2,ppm_dim
            IF (ldc(2) < n(i)) THEN
               ldc(2) = n(i)
            ENDIF
         ENDDO
         ldc(3) = 2
         CALL ppm_alloc(res,ldc,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_dealloc,'ppm_tree',     &
      &        'list of sorted neighboring constraints in get_subarray',__LINE__,info)
         ENDIF
         iopt = ppm_param_alloc_fit
         CALL ppm_alloc(res,ldc,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_dealloc,'ppm_tree',     &
      &        'list of sorted neighboring constraints in get_subarray',__LINE__,info)
         ENDIF

         ! Store the sorted array in res, Insertion sort
         ! sort in all dimensions
         DO k = 1,ppm_dim

            IF(n(k) .GT. 0) THEN
               res(k,1,1) = neigh_ghost(k,1,1)
               res(k,1,2) = neigh_ghost(k,1,2)
               DO i = 2,n(k)
                  j = i - 1
                  temp1 = neigh_ghost(k,i,1)
                  temp2 = neigh_ghost(k,i,2)

                  DO WHILE (j .GE. 2 .AND. res(k,j,1) .GT. temp1)
                     res(k,j+1,1) = res(k,j,1)
                     res(k,j+1,2) = res(k,j,2)
                     j = j - 1
                  END DO
                  IF (j .EQ. 1) THEN
                     IF (res(k,j,1) .GT. temp1) THEN
                        res(k,j+1,1) = res(k,j,1)
                        res(k,j+1,2) = res(k,j,2)
                        j = j - 1
                     ENDIF
                  ENDIF

                  res(k,j+1,1) = temp1
                  res(k,j+1,2) = temp2
               END DO
            ENDIF

         ENDDO

      END SUBROUTINE get_subarray

#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_tree_divcheck_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_tree_divcheck_d
#endif
