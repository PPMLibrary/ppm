      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_tree_alloc
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

#if   __TYPE == __TREE
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_tree_alloc_ts(iopt,nbox,nbpd,nlevel,min_box,     &
     &    max_box,boxcost,parent,nchld,child,blevel,nbpl,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_tree_alloc_td(iopt,nbox,nbpd,nlevel,min_box,     &
     &    max_box,boxcost,parent,nchld,child,blevel,nbpl,info)
#endif
#elif __TYPE == __DECOMP
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_tree_alloc_ds(iopt,nbox,nbpd,min_box,max_box,    &
      &    boxcost,nchld,blevel,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_tree_alloc_dd(iopt,nbox,nbpd,min_box,max_box,    &
      &    boxcost,nchld,blevel,info)
#endif
#endif
      !!! This routine (re)allocates the tree data structures.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_tree
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
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
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), POINTER       :: min_box
      !!! Lower coordinates of the box.
      !!!
      !!! 1st index: x,y[,z]                                                   +
      !!! 2nd: box ID
      REAL(MK), DIMENSION(:,:), POINTER       :: max_box
      !!! Upper coordinates of the box.
      !!!
      !!! 1st index: x,y[,z]                                                   +
      !!! 2nd: box ID
      REAL(MK), DIMENSION(:  ), POINTER       :: boxcost
      !!! Cost of all the boxes.
      INTEGER , DIMENSION(:  ), POINTER       :: nchld
      !!! Number of children of each box.
      INTEGER                 , INTENT(IN   ) :: iopt
      !!! Allocation mode (passed on to ppm_alloc)
      INTEGER                 , INTENT(IN   ) :: nbox
      !!! New number of boxes to allocate
      INTEGER                 , INTENT(IN   ) :: nbpd
      !!! Number of children per parent
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      INTEGER , DIMENSION(:  ), POINTER       :: blevel
      !!! Tree level of each box
#if   __TYPE == __TREE
      INTEGER                 , INTENT(IN   ) :: nlevel
      !!! Number of levels in the tree
      INTEGER , DIMENSION(:  ), POINTER       :: parent
      !!! Index of the parent box of each box. `ppm_param_undefined` if no
      !!! parent (i.e. root box)
      INTEGER , DIMENSION(:  ), POINTER       :: nbpl
      !!! The number of boxes per level
      INTEGER , DIMENSION(:,:), POINTER       :: child
      !!! Indices of all children of a box. 1st index: child ID, 2nd: box ID.
#endif
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK) :: t0

      INTEGER, DIMENSION(2) :: ldc

      CHARACTER(LEN=ppm_char) :: caller = 'ppm_tree_alloc'
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      !-------------------------------------------------------------------------
      !  Check input arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         CALL check
         IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate
      !-------------------------------------------------------------------------
      IF (have_particles) THEN
          ldc(1) = 2
          ldc(2) = nbox
          CALL ppm_alloc(tree_lhbx,ldc,iopt,info)
          or_fail_alloc('pointer to headers TREE_LHBX',ppm_error=ppm_error_fatal)
      ENDIF

      ldc(1) = ppm_dim
      ldc(2) = nbox
      CALL ppm_alloc(min_box,ldc,iopt,info)
      or_fail_alloc('lower box boundaries MIN_BOX',ppm_error=ppm_error_fatal)

      CALL ppm_alloc(max_box,ldc,iopt,info)
      or_fail_alloc('upper box boundaries MAX_BOX',ppm_error=ppm_error_fatal)

      IF (have_mesh) THEN
         CALL ppm_alloc(Nm_box,ldc,iopt,info)
         or_fail_alloc('box grid size NM_BOX',ppm_error=ppm_error_fatal)
      ENDIF

      ldc(1) = nbox
      CALL ppm_alloc(ndiv,ldc,iopt,info)
      or_fail_alloc('number of divisible directions NDIV',ppm_error=ppm_error_fatal)

      CALL ppm_alloc(blevel,ldc,iopt,info)
      or_fail_alloc('tree levels of boxes BLEVEL',ppm_error=ppm_error_fatal)

      CALL ppm_alloc(boxcost,ldc,iopt,info)
      or_fail_alloc('box costs BOXCOST',ppm_error=ppm_error_fatal)

      CALL ppm_alloc(nchld,ldc,iopt,info)
      or_fail_alloc('number of children NCHLD',ppm_error=ppm_error_fatal)

#if   __TYPE == __TREE
      ldc(1) = nbpd
      ldc(2) = nbox
      CALL ppm_alloc(child,ldc,iopt,info)
      or_fail_alloc('list of children CHILD',ppm_error=ppm_error_fatal)

      ldc(1) = nbox
      CALL ppm_alloc(parent,ldc,iopt,info)
      or_fail_alloc('parent pointer PARENT',ppm_error=ppm_error_fatal)

      ldc(1) = nlevel
      CALL ppm_alloc(nbpl,ldc,iopt,info)
      or_fail_alloc('number of boxes per level NBPL',ppm_error=ppm_error_fatal)
#endif

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
         IF (nbox .LT. 0) THEN
            fail('Number of boxes must be >= 0',exit_point=8888)
         ENDIF
         IF (nbpd .LT. 0) THEN
            fail('Number of boxes per step must be >= 0',exit_point=8888)
         ENDIF
#if   __TYPE == __TREE
         IF (nlevel .LT. 0) THEN
            fail('Number of levels must be >= 0',exit_point=8888)
         ENDIF
#endif
      8888 CONTINUE
      END SUBROUTINE check
#if   __TYPE == __TREE
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_tree_alloc_ts
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_tree_alloc_td
#endif
#elif __TYPE == __DECOMP
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_tree_alloc_ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_tree_alloc_dd
#endif
#endif
