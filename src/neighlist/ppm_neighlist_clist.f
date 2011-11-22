      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_neighlist_clist
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
      SUBROUTINE ppm_neighlist_clist_s(topoid,xp,np,cutoff,lsymm,clist, &
     &                                 info,pidx)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_neighlist_clist_d(topoid,xp,np,cutoff,lsymm,clist, &
     &                                 info,pidx)
      !!! Create cell lists for all subs of this processor.
      !!!
      !!! [NOTE]
      !!! ====================================================
      !!! Symmetry is used as follows:
      !!!
      !!! ----------------------------------------------
      !!!                                   - 2 1
      !!!                                   - 0 2
      !!!                                   - - -
      !!! ----------------------------------------------
      !!! cell 0 interacts with all cells >0
      !!! cells 2 interact with each other.
      !!!
      !!! In 3D, the top Layer (larger z) above the numbered
      !!! 2 x 2 block is also included. (Plus the appropriate
      !!! diagonal interactions. see MkNeighIdx).
      !!!
      !!! The cell list is created for the current topology
      !!! topoid, which is passed by the user
      !!!
      !!! Particles in cell icell of sub isub are:                             +
      !!!   `LET a = clist(isub)%lhbx(icell)`                                  +
      !!!   `LET b = clist(isub)%lhbx(icell+1)                                 +
      !!!   `clist(isub)%lpdx(a:b-1)`
      !!! ====================================================
      !!!
      !!! [WARNING]
      !!! the sub index isub is NOT the global sub
      !!! number, but just linear from 1 to ppm_nsublist!
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_typedef
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_util_rank
      USE ppm_module_check_id
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:),          POINTER :: xp
      !!! Particle co-ordinates
      INTEGER                 , INTENT(IN   )    :: np
      !!! Number of particles
      INTEGER                 , INTENT(IN   )    :: topoid
      !!! ID of current topology
      REAL(MK), DIMENSION(:)  , INTENT(IN   )    :: cutoff
      !!! Cutoff in all (2,3) space directions. Actual cell size may differ,
      !!! due to round-off, but it always >= cutoff.
      LOGICAL                 , INTENT(IN   )    :: lsymm
      !!! Use symmetry?
      TYPE(ppm_t_clist), DIMENSION(:), POINTER   :: clist
      !!! Cell list data structure
      INTEGER                 , INTENT(  OUT)    :: info
      !!! Returns 0 upon success
      INTEGER, DIMENSION(:)   , OPTIONAL         :: pidx
      !!! Indices of those particles that are
      !!! to be ranked. By default, all particles are ranked. If given,
      !!! particle indices in cell list are relative to `xp(:,pidx(:))` and not
      !!! `xp(:,:)`
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! timer
      REAL(MK)                                :: t0
      ! counters
      REAL(MK), DIMENSION(:,:), POINTER       :: wxp => NULL()
      ! work array, must be allocated of pidx is passed
      INTEGER                                 :: idom,jdom,i,npidx,wnp
      ! extent of cell mesh
      REAL(MK), DIMENSION(ppm_dim)            :: xmin,xmax
      ! actual cell size
      REAL(MK), DIMENSION(ppm_dim)            :: min_phys,max_phys
      ! domain extents
      REAL(MK), DIMENSION(ppm_dim)            :: cellsize
      ! number of ghostlayers
      INTEGER, DIMENSION(2*ppm_dim)           :: ngl
      ! parameter for alloc
      INTEGER                                 :: nsbc
      LOGICAL, DIMENSION(2*ppm_dim)           :: isbc
      INTEGER                                 :: iopt
      INTEGER, DIMENSION(2)                   :: ldc
      LOGICAL                                 :: valid
      TYPE(ppm_t_topo)          , POINTER     :: topo => NULL()
      REAL(MK)                                :: eps
      LOGICAL                                 :: lpidx
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_neighlist_clist',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF
#if   __KIND == __DOUBLE_PRECISION
      eps = ppm_myepsd
#elif __KIND == __SINGLE_PRECISION
      eps = ppm_myepss
#endif

      topo => ppm_topo(topoid)%t

      !-------------------------------------------------------------------------
      !  Allocate or set wxp
      !-------------------------------------------------------------------------
      IF (PRESENT(pidx)) THEN
          lpidx = .TRUE.
          npidx = SIZE(pidx,1)
          iopt = ppm_param_alloc_fit
          ldc(1) = ppm_dim
          ldc(2) = npidx
          CALL ppm_alloc(wxp,ldc,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_neighlist_clist',  &
 &                   'work xp array wxp',__LINE__,info)
              GOTO 9999
          ENDIF
          FORALL(i=1:npidx) wxp(:,i) = xp(:,pidx(i))
          wnp = npidx
      ELSE
          lpidx = .FALSE.
          wxp => xp
          wnp = np
      ENDIF
      !-------------------------------------------------------------------------
      !  Allocate clist to the number of subs this processor has
      !-------------------------------------------------------------------------
      IF (ASSOCIATED(clist)) THEN
          IF (SIZE(clist,1) .NE. topo%nsublist) THEN
              CALL ppm_clist_destroy(clist,info)
              IF (info .NE. 0) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_sub_failed,'ppm_neighlist_clist',  &
     &                'Could not destroy old cell list',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDIF
      ENDIF
      IF (.NOT. ASSOCIATED(clist)) THEN
          ALLOCATE(clist(topo%nsublist), STAT=info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_neighlist_clist',  &
     &                'cell list array CLIST',__LINE__,info)
              GOTO 9999
          ENDIF
          DO i=1,topo%nsublist
              NULLIFY(clist(i)%nm)
              NULLIFY(clist(i)%lpdx)
              NULLIFY(clist(i)%lhbx)
          ENDDO
      ENDIF

#if   __KIND == __DOUBLE_PRECISION
      min_phys(:) = topo%min_physd
      max_phys(:) = topo%max_physd
#elif __KIND == __SINGLE_PRECISION
      min_phys(:) = topo%min_physs
      max_phys(:) = topo%max_physs
#endif

      !-------------------------------------------------------------------------
      ! Determine if there are any (non-)symmetric boundary conditions
      !-------------------------------------------------------------------------
      nsbc = 0
      isbc(:) = .FALSE. 
      DO i=1,2*ppm_dim
          SELECT CASE (topo%bcdef(i))
          CASE (ppm_param_bcdef_symmetry)
              nsbc = nsbc + 1
              isbc(i) = .TRUE.
          CASE (ppm_param_bcdef_antisymmetry)
              nsbc = nsbc + 1
              isbc(i) = .TRUE.
          CASE (ppm_param_bcdef_neumann)
              nsbc = nsbc + 1
              isbc(i) = .TRUE.
          CASE (ppm_param_bcdef_dirichlet)
              nsbc = nsbc + 1
              isbc(i) = .TRUE.
          END SELECT
      ENDDO

      !-------------------------------------------------------------------------
      !  Loop over all subs of this processor and create cell lists for all
      !  of them. The extent of the cell list mesh is larger than the sub
      !  (but not on all sides if we are using symmetry). 
      !-------------------------------------------------------------------------
      DO idom=1,topo%nsublist
          jdom = topo%isublist(idom)

          !---------------------------------------------------------------------
          !  Determine extent of mesh around this sub plus tolerance for
          !  outer boundaries (to avoid real particles in ghost cells)
          !---------------------------------------------------------------------
#if   __KIND == __DOUBLE_PRECISION
              !-----------------------------------------------------------------
              !  Get sub extent -- 2D double precision
              !-----------------------------------------------------------------
              xmin(1) = topo%min_subd(1,jdom)
              xmax(1) = topo%max_subd(1,jdom)
              xmin(2) = topo%min_subd(2,jdom)
              xmax(2) = topo%max_subd(2,jdom)
          IF (ppm_dim .EQ. 3) THEN
              !-----------------------------------------------------------------
              !  Get sub extent -- 3D double precision
              !-----------------------------------------------------------------
              xmin(3) = topo%min_subd(3,jdom)
              xmax(3) = topo%max_subd(3,jdom)
          ENDIF
#elif __KIND == __SINGLE_PRECISION
              !-----------------------------------------------------------------
              !  Get sub extent -- 2D single precision
              !-----------------------------------------------------------------
              xmin(1) = topo%min_subs(1,jdom)
              xmax(1) = topo%max_subs(1,jdom)
              xmin(2) = topo%min_subs(2,jdom)
              xmax(2) = topo%max_subs(2,jdom)
          IF (ppm_dim .EQ. 3) THEN
              !-----------------------------------------------------------------
              !  Get sub extent -- 3D single precision
              !-----------------------------------------------------------------
              xmin(3) = topo%min_subs(3,jdom)
              xmax(3) = topo%max_subs(3,jdom)
          ENDIF
#endif

          !---------------------------------------------------------------------
          !  Allocate cell number array
          !---------------------------------------------------------------------
          iopt = ppm_param_alloc_fit
          ldc(1) = ppm_dim
          CALL ppm_alloc(clist(idom)%nm,ldc,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_neighlist_clist',  &
 &                   'Numbers of cells NM',__LINE__,info)
              GOTO 9999
          ENDIF
          !---------------------------------------------------------------------
          !  Determine number of cell boxes and effective cell size.
          !---------------------------------------------------------------------
          DO i=1,ppm_dim
              ! number of cells based on a cellsize = cutoff 
              clist(idom)%nm(i) = INT((xmax(i) - xmin(i))/cutoff(i))
              ! make at least one box
              IF (clist(idom)%nm(i) .LT. 1) clist(idom)%nm(i) = 1
              cellsize(i) = (xmax(i) - xmin(i))/REAL(clist(idom)%nm(i),MK)
          ENDDO

          !---------------------------------------------------------------------
          !  Find out how many ghost layers are needed:
          !  Do no longer add ghost layers to the domain since this would result
          !  in errors in the ranking due to nummerical errors
          !---------------------------------------------------------------------
          ngl(:) = 0
          IF (lsymm) THEN                ! EXPLOIT SYMMETRY
              DO i=1,ppm_dim
                  ! if we are at are the phys_dom border and have (non-)symmetirc
                  ! BCs then add a ghost layer
                  IF ((ABS(xmin(i)-min_phys(i)).LT.eps).AND.isbc((i-1)*2+1)) THEN
                      ngl(i) = 1
                  ENDIF
              ENDDO
              DO i=ppm_dim+1,2*ppm_dim  ! layers on upper-right side ngl(i) = 1
                  ngl(i) = 1
              ENDDO
          ELSE                       ! DO NOT EXPLOIT SYMMETRY => ghost layers 
              DO i=1,2*ppm_dim             ! all around
                 ngl(i) = 1
              ENDDO
          ENDIF

          !---------------------------------------------------------------------
          !  Rank the particles in this extended sub
          !---------------------------------------------------------------------
          IF (ppm_dim .EQ. 2) THEN
              CALL ppm_util_rank2d(wxp,wnp,xmin(1:2),xmax(1:2),&
     &                 clist(idom)%nm(1:2),ngl(1:4),clist(idom)%lpdx,&
     &                 clist(idom)%lhbx,info)
              !-----------------------------------------------------------------
              !  We have to increase nm by the ghost layers to provide the same
              !  behaviour as before the change of interface of ppm_util_rank
              !-----------------------------------------------------------------
              clist(idom)%nm(1) = clist(idom)%nm(1) + ngl(1) + ngl(3)
              clist(idom)%nm(2) = clist(idom)%nm(2) + ngl(2) + ngl(4)
          ELSEIF (ppm_dim .EQ. 3) THEN
              CALL ppm_util_rank3d(wxp,wnp,xmin(1:3),xmax(1:3),&
     &                 clist(idom)%nm(1:3),ngl(1:6),clist(idom)%lpdx,&
     &                 clist(idom)%lhbx,info)
              !-----------------------------------------------------------------
              !  We have to increase nm by the ghost layers to provide the same
              !  behaviour as before the change of interface of ppm_util_rank
              !-----------------------------------------------------------------
              clist(idom)%nm(1) = clist(idom)%nm(1) + ngl(1) + ngl(4)
              clist(idom)%nm(2) = clist(idom)%nm(2) + ngl(2) + ngl(5)
              clist(idom)%nm(3) = clist(idom)%nm(3) + ngl(3) + ngl(6)
          ENDIF
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_sub_failed,'ppm_neighlist_clist',   &
     &                  'ranking of particles failed!',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDDO
      IF (lpidx) THEN
          iopt = ppm_param_dealloc
          CALL ppm_alloc(wxp,ldc,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_neighlist_clist',  &
 &                   'work xp array wxp',__LINE__,info)
              GOTO 9999
          ENDIF
      ELSE
          NULLIFY(wxp)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_neighlist_clist',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_neighlist_clist',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (SIZE(cutoff,1) .LT. ppm_dim) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_neighlist_clist',  &
     &            'cutoff must be given in all dimensions',__LINE__,info)
              GOTO 8888
          ENDIF
          DO i=1,ppm_dim
              IF (cutoff(i) .LE. 0.0_MK) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_neighlist_clist',  &
     &                'cutoff must be >0 in all dimensions',__LINE__,info)
                  GOTO 8888
              ENDIF
          ENDDO
          IF (np .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_neighlist_clist',  &
     &            'np must be >0',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (topoid .EQ. ppm_param_topo_undefined) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_neighlist_clist',  &
     &            'A defined topology must be given',__LINE__,info)
              GOTO 8888
          ENDIF
          CALL ppm_check_topoid(topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_neighlist_clist',  &
     &               'topoid out of range',__LINE__,info)
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_neighlist_clist_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_neighlist_clist_d
#endif
