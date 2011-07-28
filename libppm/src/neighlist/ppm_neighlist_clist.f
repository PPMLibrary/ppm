      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_neighlist_clist
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_neighlist_clist_s(topoid,xp,Np,cutoff,lsymm,clist,Nm, &
     &                                 info,pidx)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_neighlist_clist_d(topoid,xp,Np,cutoff,lsymm,clist,Nm, &
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
      !!! In 3D, the TOP LAYER (larger z) above the numbered
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
#include "ppm_define.h"

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_neighlist
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
      REAL(MK), DIMENSION(:,:), INTENT(IN   )    :: xp
      !!! Particle co-ordinates
      INTEGER                 , INTENT(IN   )    :: Np
      !!! Number of particles
      INTEGER                 , INTENT(IN   )    :: topoid
      !!! ID of current topology
      REAL(MK), DIMENSION(:)  , INTENT(IN   )    :: cutoff
      !!! Cutoff in all (2,3) space directions. Actual cell size may differ,
      !!! due to round-off, but it always >= cutoff.
      LOGICAL                 , INTENT(IN   )    :: lsymm
      !!! Use symmetry?
      TYPE(ppm_type_ptr_to_clist), DIMENSION(:), POINTER :: clist
      !!! Number of cells in each space direction
      !!! clist(isub)%lhbx(ibox)
      INTEGER, DIMENSION(:,:) , POINTER          :: Nm
      !!! Number of cells in x,y,(z) direction (including the ghosts cells)
      !!! in each subdomain. 1st index: direction. second index: subid.
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
      INTEGER                                 :: idom,jdom,i,Npdx
      ! extent of cell mesh
      REAL(MK), DIMENSION(3)                  :: xmin,xmax
      ! actual cell size
      REAL(MK), DIMENSION(3)                  :: cellsize
      ! number of ghostlayers
      INTEGER, DIMENSION(6)                   :: Ngl
      ! parameter for alloc
      INTEGER                                 :: lda,iopt
      INTEGER, DIMENSION(2)                   :: ldc
      LOGICAL                                 :: valid
      TYPE(ppm_t_topo)          , POINTER     :: topo
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

      topo => ppm_topo(topoid)%t

      !-------------------------------------------------------------------------
      !  Allocate Nm
      !-------------------------------------------------------------------------
      lda = topo%nsublist
      iopt = ppm_param_alloc_fit
      ldc(1) = ppm_dim
      ldc(2) = lda
      CALL ppm_alloc(Nm,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_neighlist_clist',  &
     &            'Numbers of cells NM',__LINE__,info)
          GOTO 9999
      ENDIF
      Nm = 0

      !-------------------------------------------------------------------------
      !  Allocate clist to the number of subs this processor has
      !-------------------------------------------------------------------------
      IF (ASSOCIATED(clist)) THEN
          IF (SIZE(clist,1) .NE. lda) THEN
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
          ALLOCATE(clist(lda), STAT=info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_neighlist_clist',  &
     &                'cell list array CLIST',__LINE__,info)
              GOTO 9999
          ENDIF
          DO i=1,lda
              NULLIFY(clist(i)%lpdx)
              NULLIFY(clist(i)%lhbx)
          ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Loop over all subs of this processor and create cell lists for all
      !  of them. The extent of the cell list mesh is larger than the sub
      !  (but not on all sides if we are using symmetry). 
      !-------------------------------------------------------------------------
      DO idom=1,lda
          jdom = topo%isublist(idom)

          !---------------------------------------------------------------------
          !  Determine extent of mesh around this sub plus tolerance for
          !  outer boundaries (to avoid real particles in ghost cells)
          !---------------------------------------------------------------------
#if   __KIND == __DOUBLE_PRECISION
          IF (ppm_dim .EQ. 2) THEN
              !-----------------------------------------------------------------
              !  Get sub extent -- 2D double precision
              !-----------------------------------------------------------------
              xmin(1) = topo%min_subd(1,jdom)
              xmax(1) = topo%max_subd(1,jdom)
              xmin(2) = topo%min_subd(2,jdom)
              xmax(2) = topo%max_subd(2,jdom)
          ELSE
              !-----------------------------------------------------------------
              !  Get sub extent -- 3D double precision
              !-----------------------------------------------------------------
              xmin(1) = topo%min_subd(1,jdom)
              xmax(1) = topo%max_subd(1,jdom)
              xmin(2) = topo%min_subd(2,jdom)
              xmax(2) = topo%max_subd(2,jdom)
              xmin(3) = topo%min_subd(3,jdom)
              xmax(3) = topo%max_subd(3,jdom)
          ENDIF
#elif __KIND == __SINGLE_PRECISION
          IF (ppm_dim .EQ. 2) THEN
              !-----------------------------------------------------------------
              !  Get sub extent -- 2D single precision
              !-----------------------------------------------------------------
              xmin(1) = topo%min_subs(1,jdom)
              xmax(1) = topo%max_subs(1,jdom)
              xmin(2) = topo%min_subs(2,jdom)
              xmax(2) = topo%max_subs(2,jdom)
          ELSE
              !-----------------------------------------------------------------
              !  Get sub extent -- 3D single precision
              !-----------------------------------------------------------------
              xmin(1) = topo%min_subs(1,jdom)
              xmax(1) = topo%max_subs(1,jdom)
              xmin(2) = topo%min_subs(2,jdom)
              xmax(2) = topo%max_subs(2,jdom)
              xmin(3) = topo%min_subs(3,jdom)
              xmax(3) = topo%max_subs(3,jdom)
          ENDIF
#endif

          !---------------------------------------------------------------------
          !  Determine number of cell boxes and effective cell size.
          !---------------------------------------------------------------------
          DO i=1,ppm_dim
              ! number of cells based on a cellsize = cutoff 
              Nm(i,idom) = INT((xmax(i) - xmin(i))/cutoff(i))
              ! make at least one box
              IF (Nm(i,idom) .LT. 1) Nm(i,idom) = 1
              cellsize(i) = (xmax(i) - xmin(i))/REAL(Nm(i,idom),MK)
          ENDDO

          !---------------------------------------------------------------------
          !  Find out how many ghost layers are needed:
          !  Do no longer add ghost layers to the domain since this would result
          !  in errors in the ranking due to nummerical errors
          !---------------------------------------------------------------------
          Ngl(1:6) = 0
          IF (lsymm) THEN                ! EXPLOIT SYMMETRY => only need ghost 
              DO i=1,ppm_dim             ! layers on one side
                  Ngl(ppm_dim + i) = 1
              ENDDO
          ELSE                       ! DO NOT EXPLOIT SYMMETRY => ghost layers 
              DO i=1,ppm_dim             ! all around
                 Ngl(i) = 1
                 Ngl(ppm_dim + i) = 1
              ENDDO
          ENDIF

          !---------------------------------------------------------------------
          !  Rank the particles in this extended sub
          !---------------------------------------------------------------------
          IF (ppm_dim .EQ. 2) THEN
              IF (PRESENT(pidx)) THEN
                  Npdx = SIZE(pidx,1)
                  CALL ppm_util_rank2d(xp(1:2,pidx),Npdx,xmin(1:2),xmax(1:2),&
     &                    Nm(1:2,idom),Ngl(1:4),clist(idom)%lpdx,&
     &                    clist(idom)%lhbx,info)
              ELSE
                  CALL ppm_util_rank2d(xp,Np,xmin(1:2),xmax(1:2),Nm(1:2,idom),&
     &                    Ngl(1:4),clist(idom)%lpdx,clist(idom)%lhbx,info)
              ENDIF
              !-----------------------------------------------------------------
              !  We have to increase Nm by the ghost layers to provide the same
              !  behaviour as before the change of interface of ppm_util_rank
              !-----------------------------------------------------------------
              Nm(1,idom) = Nm(1,idom) + Ngl(1) + Ngl(3)
              Nm(2,idom) = Nm(2,idom) + Ngl(2) + Ngl(4)
          ELSEIF (ppm_dim .EQ. 3) THEN
              IF (PRESENT(pidx)) THEN
                  Npdx = SIZE(pidx,1)
                  CALL ppm_util_rank3d(xp(1:3,pidx),Npdx,xmin(1:3),xmax(1:3),&
     &                    Nm(1:3,idom),Ngl(1:6),clist(idom)%lpdx,&
     &                    clist(idom)%lhbx,info)
              ELSE
                  CALL ppm_util_rank3d(xp,Np,xmin(1:3),xmax(1:3),Nm(1:3,idom),&
     &                    Ngl(1:6),clist(idom)%lpdx,clist(idom)%lhbx,info)
              ENDIF
              !-----------------------------------------------------------------
              !  We have to increase Nm by the ghost layers to provide the same
              !  behaviour as before the change of interface of ppm_util_rank
              !-----------------------------------------------------------------
              Nm(1,idom) = Nm(1,idom) + Ngl(1) + Ngl(4)
              Nm(2,idom) = Nm(2,idom) + Ngl(2) + Ngl(5)
              Nm(3,idom) = Nm(3,idom) + Ngl(3) + Ngl(6)
          ENDIF
          IF (info .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_sub_failed,'ppm_neighlist_clist',   &
     &                  'ranking of particles failed!',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDDO

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
          IF (Np .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_neighlist_clist',  &
     &            'Np must be >0',__LINE__,info)
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
