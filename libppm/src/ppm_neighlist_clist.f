      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_neighlist_clist
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Create cell lists for all subs of this processor.
      !
      !  Input        : xp(:,:)    (F) particle co-ordinates
      !                 Np         (I) number of particles
      !                 lsymm      (L) T to use symmetry, else do not use
      !                                symmetry
      !                 pidx(:)    (I) indices of those particles that are
      !                                to be ranked (OPTIONAL). By default,
      !                                all particles are ranked. If given,
      !                                particle indices in cell list are
      !                                relative to xp(:,pidx(:)) and not
      !                                xp(:,:)!
      !                 cutoff(:)  (F) cutoff in all (2,3) space
      !                                directions. Actual cell size may differ,
      !                                due to round-off, but it always .GE.
      !                                cutoff.
      !
      !  Input/output : 
      !
      !  Output       : info       (I) return status. =0 if no error.
      !                 clist(:)       Cell list as a list of ptr_to_clist.
      !                                particle index list of isub: 
      !                                        clist(isub)%lpdx(:)
      !                                pointer to first particle in ibox of
      !                                isub: 
      !                                        clist(isub)%lhbx(ibox)
      !                 Nm(:,:)    (I) number of cells in x,y,(z)
      !                                direction (including the ghosts cells) 
      !                                in each subdomain. 1st index:
      !                                direction. second index: subid.
      !
      !  Remarks      : Symmetry is used as follows:      - 2 1
      !                                                   - 0 2
      !                                                   - - -
      !
      !                      cell 0 interacts with all cells >0
      !                      cells 2 interact with each other.
      !
      !                 In 3D, the TOP LAYER (larger z) above the numbered
      !                 2 x 2 block is also included. (Plus the appropriate
      !                 diagonal interactions. see MkNeighIdx).
      !
      !                 The cell list is created for the CURRENT TOPOLOGY
      !                 ppm_topoid.
      !
      !                 Particles in cell icell of sub isub are:
      !                   LET a = clist(isub)%lhbx(icell)
      !                   LET b = clist(isub)%lhbx(icell+1)
      !                   clist(isub)%lpdx(a:b-1)
      !
      !                 ATTENTION: the sub index isub is NOT the global sub
      !                 number, but just linear from 1 to ppm_nsublist!
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_neighlist_clist.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.25  2006/07/04 15:51:35  ivos
      !  Added explicit array boundaries for initialization of Ngl.
      !
      !  Revision 1.24  2006/04/10 12:42:21  pchatela
      !  Comment at the end of a line containing the USE word:
      !  quite disorienting for the Makefile. Replaced it with EXPLOIT.
      !  This will happen every time a comment containing the word USE does
      !  not begin on a new line... Will think about a fix..
      !
      !  Revision 1.23  2006/03/28 19:37:02  ivos
      !  Added comments.
      !
      !  Revision 1.22  2004/10/28 12:38:19  davidch
      !  Fixed numerical bug in cell lists that resulted in real particles 
      !  being treated as ghosts
      !  and vice versa. The new ranking and cell list routines are 
      !  supposed to be exact. All
      !  epsilons that were added to the domains in order to prevent the 
      !  mentioned problems were
      !  removed since they are no longer needed.
      !  Modified Files:
      !      ppm_util_rank2d.f ppm_util_rank3d.f ppm_util_sort2d.f
      !      ppm_util_sort3d.f ppm_find_duplicates.f ppm_neighlist_clist.f
      !      ppm_error.h
      !
      !  Revision 1.21  2004/10/14 11:28:18  ivos
      !  bugfix: eps should only be added to OUTER faces and not to internal
      !  sub boundaries. Otherwise we might end up having the same particle
      !  assigned to more than one cell list!
      !
      !  Revision 1.20  2004/10/13 07:36:19  hiebers
      !  enlarged domain of the cells by the tolerance ppm_myeps to ensure that
      !  particles on the boundary are sorted in as well
      !
      !  Revision 1.19  2004/10/01 16:09:11  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.18  2004/09/28 15:22:50  hiebers
      !  BUG FIX: eliminated a typo in SINGLE_PRECISION
      !
      !  Revision 1.17  2004/07/26 15:38:50  ivos
      !  Inserted missing USE statements to resolve undefined references
      !  at link stage.
      !
      !  Revision 1.16  2004/07/26 11:59:39  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.15  2004/07/26 07:42:49  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.14  2004/07/16 14:47:21  ivos
      !  Added check for ppm_initialized.
      !
      !  Revision 1.13  2004/06/15 08:44:29  ivos
      !  bugfixes: (1) changed cutoff to INTENT(IN) since the actual cell size
      !  can be different on each sub. (2) changed Nm to be a 2d array since the
      !  number of cells depends on the subdomain id.
      !
      !  Revision 1.12  2004/06/11 10:20:57  ivos
      !  bugfix: problem if domain .LT. cutoff resolved.
      !
      !  Revision 1.11  2004/02/24 11:36:38  ivos
      !  Changed imode (INTEGER symmetry flag) argument to lsymm (LOGICAL) in
      !  order to have the same interface as the ghost routines.
      !
      !  Revision 1.10  2004/02/24 11:20:11  ivos
      !  The same as gonnetp. CVS comflict resolved.
      !
      !  Revision 1.9  2004/02/24 09:14:49  gonnetp
      !  nullify pointers in clist after allocation.
      !
      !  Revision 1.8  2004/02/04 17:20:48  ivos
      !  renamed TYPE ptr_to_clist to ppm_type_ptr_to_clist.
      !
      !  Revision 1.7  2004/01/27 14:51:40  ivos
      !  Added success check after CALL to ppm_neighlist_clist_destroy.
      !
      !  Revision 1.6  2004/01/23 17:24:17  ivos
      !  Now includes ppm_define.h for the cpp defines and no longer sets them
      !  in the Makefile.
      !
      !  Revision 1.5  2004/01/22 14:17:52  ivos
      !  Bugfix: cutoff vector is now checked element-wise.
      !
      !  Revision 1.4  2004/01/22 13:27:56  ivos
      !  Did (1) update of the header, (2) replaced pwrite with ppm_write or
      !  ppm_error calls, (3) inserted validity check of INTENT(IN) arguments
      !  where needed.
      !
      !  Revision 1.3  2004/01/09 16:34:23  ivos
      !  Updated comment header.
      !
      !  Revision 1.2  2004/01/09 16:30:26  ivos
      !  Corrected ghost layer size and location (new interaction handling) and
      !  added switch to toggle symmetry use.
      !
      !  Revision 1.1  2004/01/08 17:51:54  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_neighlist_clist_s(xp,Np,cutoff,lsymm,clist,Nm,info,pidx)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_neighlist_clist_d(xp,Np,cutoff,lsymm,clist,Nm,info,pidx)
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
      USE ppm_module_clist_destroy
      USE ppm_module_util_rank
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
      INTEGER                 , INTENT(IN   )    :: Np
      REAL(MK), DIMENSION(:)  , INTENT(IN   )    :: cutoff
      LOGICAL                 , INTENT(IN   )    :: lsymm
      TYPE(ppm_type_ptr_to_clist), DIMENSION(:), POINTER :: clist
      ! number of cells in each space direction
      INTEGER, DIMENSION(:,:) , POINTER          :: Nm
      INTEGER                 , INTENT(  OUT)    :: info
      INTEGER, DIMENSION(:)   , OPTIONAL         :: pidx
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
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_neighlist_clist',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (SIZE(cutoff,1) .LT. ppm_dim) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_neighlist_clist',  &
     &            'cutoff must be given in all dimensions',__LINE__,info)
              GOTO 9999
          ENDIF
          DO i=1,ppm_dim
              IF (cutoff(i) .LE. 0.0_MK) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_neighlist_clist',  &
     &                'cutoff must be >0 in all dimensions',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDDO
          IF (Np .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_neighlist_clist',  &
     &            'Np must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate Nm
      !-------------------------------------------------------------------------
      lda = ppm_nsublist(ppm_topoid)
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
          jdom = ppm_isublist(idom,ppm_topoid)

          !---------------------------------------------------------------------
          !  Determine extent of mesh around this sub plus tolerance for
          !  outer boundaries (to avoid real particles in ghost cells)
          !---------------------------------------------------------------------
#if   __KIND == __DOUBLE_PRECISION
          IF (ppm_dim .EQ. 2) THEN
              !-----------------------------------------------------------------
              !  Get sub extent -- 2D double precision
              !-----------------------------------------------------------------
              xmin(1) = ppm_min_subd(1,jdom,ppm_topoid)
              xmax(1) = ppm_max_subd(1,jdom,ppm_topoid)
              xmin(2) = ppm_min_subd(2,jdom,ppm_topoid)
              xmax(2) = ppm_max_subd(2,jdom,ppm_topoid)
          ELSE
              !-----------------------------------------------------------------
              !  Get sub extent -- 3D double precision
              !-----------------------------------------------------------------
              xmin(1) = ppm_min_subd(1,jdom,ppm_topoid)
              xmax(1) = ppm_max_subd(1,jdom,ppm_topoid)
              xmin(2) = ppm_min_subd(2,jdom,ppm_topoid)
              xmax(2) = ppm_max_subd(2,jdom,ppm_topoid)
              xmin(3) = ppm_min_subd(3,jdom,ppm_topoid)
              xmax(3) = ppm_max_subd(3,jdom,ppm_topoid)
          ENDIF
#elif __KIND == __SINGLE_PRECISION
          IF (ppm_dim .EQ. 2) THEN
              !-----------------------------------------------------------------
              !  Get sub extent -- 2D single precision
              !-----------------------------------------------------------------
              xmin(1) = ppm_min_subs(1,jdom,ppm_topoid)
              xmax(1) = ppm_max_subs(1,jdom,ppm_topoid)
              xmin(2) = ppm_min_subs(2,jdom,ppm_topoid)
              xmax(2) = ppm_max_subs(2,jdom,ppm_topoid)
          ELSE
              !-----------------------------------------------------------------
              !  Get sub extent -- 3D single precision
              !-----------------------------------------------------------------
              xmin(1) = ppm_min_subs(1,jdom,ppm_topoid)
              xmax(1) = ppm_max_subs(1,jdom,ppm_topoid)
              xmin(2) = ppm_min_subs(2,jdom,ppm_topoid)
              xmax(2) = ppm_max_subs(2,jdom,ppm_topoid)
              xmin(3) = ppm_min_subs(3,jdom,ppm_topoid)
              xmax(3) = ppm_max_subs(3,jdom,ppm_topoid)
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
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_neighlist_clist_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_neighlist_clist_d
#endif
