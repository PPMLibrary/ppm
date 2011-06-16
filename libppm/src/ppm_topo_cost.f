      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_topo_cost
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine calculates the computational cost of
      !                 each subdomain based on either the particles (if
      !                 there are any) or the mesh points. The per-particle
      !                 costs are summed up if present. If not, a unity
      !                 cost is assumed for each particle.
      !
      !  Input        : xp(:,:)      (F) paricle positions
      !                 Np           (I) number of particles. Set to .LE. 0
      !                                  if mesh-based costs are desired.
      !                 min_sub(:,:) (F) min. extent of the subdomains
      !                 max_sub(:,:) (F) max. extent of the subdomains
      !                 nsubs        (I) total number of subdomains
      !                 nnodes(:,:)  (I) number of mesh nodes in each
      !                                  direction (first index) on each
      !                                  sub (second index). If Np is .LE.
      !                                  0, this is taken to compute the
      !                                  cost. If SIZE(nnodes,2).LT.nsubs,
      !                                  and there are no particles, cost
      !                                  is computed as geometric volume.
      !                 pcost(:)     (F) per-particle costs. OPTIONAL. If
      !                                  not present, a cost of 1 per
      !                                  particle is assumed.
      !
      !  Input/output : 
      !
      !  Output       : cost(:)      (F) Aggregate cost for each subdomain
      !                 info         (I) return status
      !
      !  Remarks      : If Np.GT.0, cost is computed based on particles. If
      !                 SIZE(nnodes,2).GE.nsubs, cost based on mesh points
      !                 is computed and added to the particle cost. If
      !                 neither particles nor mesh are given, cost of a
      !                 subdomain is equal to its volume.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_topo_cost.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.11  2006/09/04 18:34:55  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.10  2004/09/24 15:00:55  ivos
      !  Extended such that particles and meshes can both be present and
      !  their costs added. Also, if neither particles nor mesh are present,
      !  the cost is now goven by the geometric volume of a sub.
      !
      !  Revision 1.9  2004/07/26 07:42:35  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.8  2004/04/20 11:33:13  ivos
      !  Comments corrected
      !
      !  Revision 1.7  2004/04/20 11:32:29  ivos
      !  bugfix: costsum must be copied to cost and not pointer-set for
      !  AllReduce to work properly.
      !
      !  Revision 1.6  2004/04/15 15:45:29  ivos
      !  Added proper handling for the ring topology 0.
      !
      !  Revision 1.5  2004/03/22 09:50:51  walther
      !  Optimized the initialization of the the particle list and bug fixed
      !  the pointer reallocation - later calls to MPI_AllReduce() would fail
      !  because of the pointer would point to the same memory and the IN/OUT
      !  of MPI_All=Reduce() must be disjoined.
      !
      !  Revision 1.4  2004/03/04 14:29:13  ivos
      !  bugfix: in argument check changed .LE. to .LT. when comparing 
      !  SIZE(nnodes) with nsubs.
      !
      !  Revision 1.3  2004/03/03 17:37:30  ivos
      !  bugfix: Added MPI_AllReduce in parallel version if the particles are
      !  used for computing the cost. Needed because each processor only
      !  knows part of the particles.
      !
      !  Revision 1.2  2004/03/03 09:13:05  ivos
      !  Removed topoid from the argument list (since the topology is not yet
      !  stored by the time this routine is called!!) and explicitly added
      !  min_sub,max_sub,nsubs,etc. instead.
      !
      !  Revision 1.1  2004/03/02 16:24:44  ivos
      !  Initial implementation. Not tested yet.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_topo_cost_s(xp,Np,min_sub,max_sub,nsubs,nnodes,cost,   &
     &              info,pcost)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_topo_cost_d(xp,Np,min_sub,max_sub,nsubs,nnodes,cost,   &
     &              info,pcost)
#endif
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
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
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xp
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: min_sub,max_sub
      REAL(MK), DIMENSION(:  ), POINTER       :: cost
      REAL(MK), DIMENSION(:  ), OPTIONAL, INTENT(IN) :: pcost
      INTEGER , DIMENSION(:,:), INTENT(IN   ) :: nnodes
      INTEGER                 , INTENT(IN   ) :: Np,nsubs
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                          :: t0
#ifdef __MPI
      REAL(MK), DIMENSION(:), POINTER   :: costsum
#endif
      INTEGER, DIMENSION(ppm_dim)       :: ldu,Nc
      REAL(MK),DIMENSION(ppm_dim)       :: len_sub
      INTEGER                           :: i,ipart,j
      INTEGER                           :: iopt,nlist1,nlist2,idom
      INTEGER, DIMENSION(:), POINTER    :: ilist1,ilist2,ilist3
      CHARACTER(LEN=ppm_char)           :: mesg
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_topo_cost',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (nsubs .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_cost',  &
     &            'nsubs must be > 0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (SIZE(min_sub,2) .LT. nsubs .OR. SIZE(max_sub,2) .LT. nsubs) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_cost',  &
     &            'not enough subs specified in min_sub, max_sub',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (SIZE(min_sub,1).LT.ppm_dim.OR.SIZE(max_sub,1).LT.ppm_dim) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_cost',  &
     &            'leading dimension in min_sub, max_sub wrong',__LINE__,info)
              GOTO 9999
          ENDIF
          DO i=1,nsubs
              DO j=1,ppm_dim
                  IF (max_sub(j,i) .LT. min_sub(j,i)) THEN
                      info = ppm_error_error
                      CALL ppm_error(ppm_err_argument,'ppm_topo_cost',  &
     &                    'max_sub must always be >= min_sub',__LINE__,info)
                      GOTO 9999
                  ENDIF
              ENDDO
          ENDDO
          IF (Np .GT. 0) THEN
              IF (SIZE(xp,2) .LT. Np) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_topo_cost',  &
     &                'not enough particles contained in xp',__LINE__,info)
                  GOTO 9999
              ENDIF
              IF (SIZE(xp,1) .LT. ppm_dim) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_topo_cost',  &
     &                'leading dimension of xp insufficient',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDIF
          IF (PRESENT(pcost)) THEN
              IF (SIZE(pcost,1) .LT. Np) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_topo_cost',  &
     &                'pcost must be of at least length Np',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for the costs
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldu(1) = nsubs
      CALL ppm_alloc(cost,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_cost',    &
     &        'costs per sub COST',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  And initialize to zero
      !-------------------------------------------------------------------------
      cost = 0.0_MK

      !-------------------------------------------------------------------------
      !  Determine the total cost of each sub either based on particles or
      !  based on mesh points
      !-------------------------------------------------------------------------
      IF (Np .GT. 0) THEN
          IF (ppm_debug .GT. 0) THEN
              CALL ppm_write(ppm_rank,'ppm_topo_cost',   &
     &            'Computing costs based on particles',info)
          ENDIF
          !---------------------------------------------------------------------
          !  Allocate memory for the particle lists
          !---------------------------------------------------------------------
          iopt = ppm_param_alloc_fit
          ldu(1) = Np
          CALL ppm_alloc(ilist1,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_topo_cost',    &
     &            'particle list 1 ILIST1',__LINE__,info)
              GOTO 9999
          ENDIF
          CALL ppm_alloc(ilist2,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_topo_cost',    &
     &            'particle list 2 ILIST2',__LINE__,info)
              GOTO 9999
          ENDIF

          !---------------------------------------------------------------------
          !  Initialize the particle list
          !---------------------------------------------------------------------
          DO ipart=1,Np
              ilist1(ipart) = ipart
          ENDDO
          nlist1 = Np

          !---------------------------------------------------------------------
          !  Loop over all subdomains
          !---------------------------------------------------------------------
          DO idom=nsubs,1,-1
              nlist2 = 0
              DO i=1,nlist1
                  ipart = ilist1(i)
                  !-------------------------------------------------------------
                  !  If the particle is inside the current subdomain, assign it
                  !-------------------------------------------------------------
                  IF (ppm_dim .GT. 2) THEN
                      IF (xp(1,ipart).GE.min_sub(1,idom).AND.   &
     &                    xp(1,ipart).LT.max_sub(1,idom).AND.   &
     &                    xp(2,ipart).GE.min_sub(2,idom).AND.   & 
     &                    xp(2,ipart).LT.max_sub(2,idom).AND.   &
     &                    xp(3,ipart).GE.min_sub(3,idom).AND.   & 
     &                    xp(3,ipart).LT.max_sub(3,idom)) THEN
                          !-----------------------------------------------------
                          !  Cost based on particles
                          !-----------------------------------------------------
                          IF (PRESENT(pcost)) THEN
                              cost(idom) = cost(idom) + pcost(ipart)
                          ELSE
                              cost(idom) = cost(idom) + 1.0_MK
                          ENDIF
                      ELSE
                         !------------------------------------------------------
                         !  if not, add it the list of particle to search
                         !------------------------------------------------------
                         nlist2         = nlist2 + 1
                         ilist2(nlist2) = ipart
                      ENDIF
                  ELSE
                      IF (xp(1,ipart).GE.min_sub(1,idom).AND.   &
     &                    xp(1,ipart).LT.max_sub(1,idom).AND.   &
     &                    xp(2,ipart).GE.min_sub(2,idom).AND.   & 
     &                    xp(2,ipart).LT.max_sub(2,idom)) THEN
                          !-----------------------------------------------------
                          !  Cost based on particles
                          !-----------------------------------------------------
                          IF (PRESENT(pcost)) THEN
                              cost(idom) = cost(idom) + pcost(ipart)
                          ELSE
                              cost(idom) = cost(idom) + 1.0_MK
                          ENDIF
                      ELSE
                         !------------------------------------------------------
                         !  if not, add it the list of particle to search
                         !------------------------------------------------------
                         nlist2         = nlist2 + 1
                         ilist2(nlist2) = ipart
                      ENDIF
                  ENDIF
              ENDDO   ! loop over particles

              !-----------------------------------------------------------------
              !  Copy the lists (well, only if nlist2 changed - decreased)
              !-----------------------------------------------------------------
              IF (nlist2.NE.nlist1) THEN
                 nlist1 = nlist2
                 ilist3 => ilist1
                 ilist1 => ilist2
                 ilist2 => ilist3
              ENDIF 

              !-----------------------------------------------------------------
              !  Exit if the list is empty
              !-----------------------------------------------------------------
              IF (nlist1.EQ.0) EXIT
          ENDDO    ! loop over subdomains
              
          !---------------------------------------------------------------------
          !  Check that we did not miss a particle
          !---------------------------------------------------------------------
          IF (nlist2 .GT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_part_unass,'ppm_topo_cost',        &
     &            'Beware: computed costs are wrong!',__LINE__,info)
          ENDIF

          !---------------------------------------------------------------------
          !  Deallocate the memory for the lists
          !---------------------------------------------------------------------
          iopt = ppm_param_dealloc
          CALL ppm_alloc(ilist2,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_topo_cost',    &
     &            'particle list 2 ILIST2',__LINE__,info)
          ENDIF
          CALL ppm_alloc(ilist1,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_topo_cost',    &
     &            'particle list 1 ILIST1',__LINE__,info)
          ENDIF

          !---------------------------------------------------------------------
          !  Do an AllReduce of the costs as one processor only has ITS
          !  paricles, but ALL particles contribute to the cost (the
          !  particles are not yet mapped!)
          !---------------------------------------------------------------------
#ifdef __MPI
          !---------------------------------------------------------------------
          !  Skip this if we only have one sub, because that means that
          !  either there is only one processor, or we have a ring
          !  topology. In the latter case the AllReduce would lead to wrong
          !  results.
          !---------------------------------------------------------------------
          IF (nsubs .GT. 1) THEN
             iopt = ppm_param_alloc_fit
             ldu(1) = nsubs
             CALL ppm_alloc(costsum,ldu,iopt,info)
             IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_topo_cost',    &
     &              'sum of costs of all processors COSTSUM',__LINE__,info)
                 GOTO 9999
             ENDIF
             costsum = 0.0_MK
#if   __KIND == __SINGLE_PRECISION
             CALL MPI_AllReduce(cost,costsum,nsubs,MPI_REAL,MPI_SUM,   &
     &          ppm_comm,info)
#elif __KIND == __DOUBLE_PRECISION
             CALL MPI_AllReduce(cost,costsum,nsubs,MPI_DOUBLE_PRECISION, &
     &          MPI_SUM,ppm_comm,info)
#endif
             !------------------------------------------------------------------
             !  Copy data back
             !------------------------------------------------------------------
             cost = costsum

             !------------------------------------------------------------------
             !  Deallocate the temporary costsum array
             !------------------------------------------------------------------
             iopt = ppm_param_dealloc
             CALL ppm_alloc(costsum,ldu,iopt,info)
             IF (info .NE. 0) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_dealloc,'ppm_topo_cost',    &
     &              'sum of costs costs COSTSUM',__LINE__,info)
             ENDIF
          ENDIF
#endif
      ENDIF

      IF (SIZE(nnodes,2) .GE. nsubs) THEN
          !---------------------------------------------------------------------
          !  Cost based on mesh points
          !---------------------------------------------------------------------
          IF (ppm_debug .GT. 0) THEN
              CALL ppm_write(ppm_rank,'ppm_topo_cost',   &
     &            'Computing costs based on mesh points',info)
          ENDIF
          IF (ppm_dim .EQ. 3) THEN
              DO i=1,nsubs
                  Nc(1) = nnodes(1,i)
                  Nc(2) = nnodes(2,i)
                  Nc(3) = nnodes(3,i)
                  cost(i) = cost(i) + (REAL(Nc(1),MK)*REAL(Nc(2),MK)*   &
     &                REAL(Nc(3),MK))
              ENDDO
          ELSE
              DO i=1,nsubs
                  Nc(1) = nnodes(1,i)
                  Nc(2) = nnodes(2,i)
                  cost(i) = cost(i) + (REAL(Nc(1),MK)*REAL(Nc(2),MK))
              ENDDO
          ENDIF
      ELSEIF (Np .LT. 1) THEN
          !---------------------------------------------------------------------
          !  Cost based on geometry if we have no particles and no mesh
          !---------------------------------------------------------------------
          IF (ppm_debug .GT. 0) THEN
              CALL ppm_write(ppm_rank,'ppm_topo_cost',   &
     &            'Computing costs based on geometry',info)
          ENDIF
          IF (ppm_dim .EQ. 3) THEN
              DO i=1,nsubs
                  len_sub(1) = max_sub(1,i)-min_sub(1,i)
                  len_sub(2) = max_sub(2,i)-min_sub(2,i)
                  len_sub(3) = max_sub(3,i)-min_sub(3,i)
                  cost(i) = len_sub(1)*len_sub(2)*len_sub(3)
              ENDDO
          ELSE
              DO i=1,nsubs
                  len_sub(1) = max_sub(1,i)-min_sub(1,i)
                  len_sub(2) = max_sub(2,i)-min_sub(2,i)
                  cost(i) = len_sub(1)*len_sub(2)
              ENDDO
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Some diagnostics
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
          CALL ppm_write(ppm_rank,'ppm_topo_cost',  &
     &           '----------------------------------------------',info)
          DO i=1,nsubs
              WRITE(mesg,'(I4,A,F15.3)') i,' sub cost: ',cost(i)
              CALL ppm_write(ppm_rank,'ppm_topo_cost',  &
     &               mesg,info)
          ENDDO
          CALL ppm_write(ppm_rank,'ppm_topo_cost',  &
     &           '----------------------------------------------',info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_topo_cost',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_topo_cost_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_topo_cost_d
#endif
