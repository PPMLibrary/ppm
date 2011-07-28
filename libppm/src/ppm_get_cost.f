      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_get_cost
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine calculates the computational cost of
      !                 each subdomain and each processor.
      !
      !  Input        : xp(:,:)      (F) paricle positions
      !                 Np           (I) number of particles. Set to .LE. 0
      !                                  if mesh-based costs are desired.
      !                 topo_id      (I) Topology ID for which to compute
      !                                  the cost (user numbering)
      !                 mesh_id      (I) mesh ID for which to compute the
      !                                  cost (user numbering). If -1 is
      !                                  passed, only particles are
      !                                  considered and the mesh is
      !                                  ignored. If there are no particles
      !                                  and mesh_id is -1, the costs are
      !                                  computed based on the geometry:
      !                                  cost(sub) = volume(sub).
      !                 pcost(:)     (F) per-particle costs. OPTIONAL. If
      !                                  not present, a cost of 1 per
      !                                  particle is assumed.
      !
      !  Input/output : 
      !
      !  Output       : cost(:)      (F) Aggregate cost for each subdomain
      !                 proc_cost(:) (F) Aggregate cost for each processor
      !                 info         (I) return status
      !
      !  Remarks      : If Np.GT.0, cost is computed based on particles. If
      !                 mesh_id.GT.-1, cost based on mesh points
      !                 is computed and added to the particle cost. If
      !                 neither particles nor mesh are given, cost of a
      !                 subdomain is equal to its volume.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_get_cost.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.8  2004/09/24 15:00:55  ivos
      !  Extended such that particles and meshes can both be present and
      !  their costs added. Also, if neither particles nor mesh are present,
      !  the cost is now goven by the geometric volume of a sub.
      !
      !  Revision 1.7  2004/08/31 13:29:57  ivos
      !  changed argument checks to use ppm_check_topoid and ppm_check_meshid.
      !
      !  Revision 1.6  2004/07/26 07:45:29  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.5  2004/07/16 14:46:25  ivos
      !  Added check for ppm_initialized.
      !
      !  Revision 1.4  2004/04/22 10:35:54  ivos
      !  bugfix: use of SIZE(ppm_internal_topid) caused problems in argument
      !  checking when first user-topoid was > 1. Resolved by use of UBOUND.
      !
      !  Revision 1.3  2004/04/15 15:45:29  ivos
      !  Added proper handling for the ring topology 0.
      !
      !  Revision 1.2  2004/03/04 14:26:23  ivos
      !  bugfix: added dummy mesh argument for the particles-only case.
      !  Routine now tested with both particles and meshes on 1,2,3 processors.
      !
      !  Revision 1.1  2004/03/03 17:36:11  ivos
      !  Initial implementation. Not yet tested.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_get_cost_s(xp,Np,topo_id,mesh_id,cost,proc_cost,info,pcost)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_get_cost_d(xp,Np,topo_id,mesh_id,cost,proc_cost,info,pcost)
#endif
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh, ONLY: ppm_cart_mesh, ppm_meshid
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      USE ppm_module_check_topoid
      USE ppm_module_check_meshid
      USE ppm_module_topo_cost
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
      REAL(MK), DIMENSION(:  ), POINTER       :: cost
      REAL(MK), DIMENSION(:  ), POINTER       :: proc_cost
      REAL(MK), DIMENSION(:  ), OPTIONAL, INTENT(IN) :: pcost
      INTEGER                 , INTENT(IN   ) :: Np,topo_id,mesh_id
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                          :: t0,mincost,maxcost
      INTEGER, DIMENSION(ppm_dim)       :: ldl,ldu
      INTEGER, DIMENSION(3,1), TARGET   :: ndummy
      INTEGER, DIMENSION(:,:), POINTER  :: nnodes
      INTEGER                           :: i,meshid,topoid,proc
      INTEGER                           :: iopt,minproc,maxproc
      LOGICAL                           :: valid
      CHARACTER(LEN=ppm_char)           :: mesg
#ifdef __MPI
      REAL(MK), DIMENSION(:), POINTER   :: sendcost, proccost
#endif
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_get_cost',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_get_cost',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (Np .GT. 0) THEN
              IF (SIZE(xp,2) .LT. Np) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_get_cost',  &
     &                'not enough particles contained in xp',__LINE__,info)
                  GOTO 9999
              ENDIF
              IF (SIZE(xp,1) .LT. ppm_dim) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_get_cost',  &
     &                'leading dimension of xp insufficient',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDIF
          IF (PRESENT(pcost)) THEN
              IF (SIZE(pcost,1) .LT. Np) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_get_cost',  &
     &                'pcost must be of at least length Np',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDIF
          IF ((Np .LE. 0) .AND. (topo_id .EQ. 0)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_get_cost',  &
     &            'null topology must use particles',__LINE__,info)
              GOTO 9999
          ENDIF
          CALL ppm_check_topoid(ppm_param_id_user,topo_id,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_get_cost',  &
     &            'topo_id is invalid',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Translate topo_id to internal numbering
      !-------------------------------------------------------------------------
      topoid = ppm_internal_topoid(topo_id)

      !-------------------------------------------------------------------------
      !  Check mesh_id validity for this topology
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (mesh_id .GT. 0) THEN
              CALL ppm_check_meshid(ppm_param_id_user,mesh_id,topoid,valid,info)
              IF (.NOT. valid) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_get_cost',  &
     &                'mesh_id is invalid',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Translate mesh_id to internal numbering and set pointer
      !-------------------------------------------------------------------------
      ndummy(1:3,1) = 0
      IF (mesh_id .GT. 0) THEN
          meshid = ppm_meshid(topoid)%internal(mesh_id)
          nnodes => ppm_cart_mesh(meshid,topoid)%nnodes
      ELSE
          meshid = -1
          nnodes => ndummy
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for the proc_costs
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldl(1) = 0
      ldu(1) = ppm_nproc-1
      CALL ppm_alloc(proc_cost,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_get_cost',    &
     &        'costs per processor PROC_COST',__LINE__,info)
          GOTO 9999
      ENDIF
      proc_cost = 0.0_MK

      !-------------------------------------------------------------------------
      !  Determine the subdomain costs either based on particles or mesh
      !  points
      !-------------------------------------------------------------------------
      IF (PRESENT(pcost)) THEN
#if   __KIND == __SINGLE_PRECISION
          CALL ppm_topo_cost(xp,Np,ppm_min_subs(:,:,topoid),       &
     &        ppm_max_subs(:,:,topoid),ppm_nsubs(topoid),          &
     &        nnodes,cost,info,pcost)
#elif __KIND == __DOUBLE_PRECISION
          CALL ppm_topo_cost(xp,Np,ppm_min_subd(:,:,topoid),       &
     &        ppm_max_subd(:,:,topoid),ppm_nsubs(topoid),          &
     &        nnodes,cost,info,pcost)
#endif
      ELSE
#if   __KIND == __SINGLE_PRECISION
          CALL ppm_topo_cost(xp,Np,ppm_min_subs(:,:,topoid),       &
     &        ppm_max_subs(:,:,topoid),ppm_nsubs(topoid),          &
     &        nnodes,cost,info)
#elif __KIND == __DOUBLE_PRECISION
          CALL ppm_topo_cost(xp,Np,ppm_min_subd(:,:,topoid),       &
     &        ppm_max_subd(:,:,topoid),ppm_nsubs(topoid),          &
     &        nnodes,cost,info)
#endif
      ENDIF
      IF (info.NE.0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_get_cost',  &
     &            'Computing costs of subdomains failed',__LINE__,info)
          GOTO 9999
      ENDIF

      NULLIFY(nnodes)

      !-------------------------------------------------------------------------
      !  Sum up the costs of all subs of the processors
      !-------------------------------------------------------------------------
      IF (topoid .EQ. 0) THEN
          !---------------------------------------------------------------------
          !  In the ring topology, each processor now knows ITS cost. These
          !  have to be communicated to all in order for all to know the
          !  cost of all processors.
          !---------------------------------------------------------------------
#ifdef __MPI
          !---------------------------------------------------------------------
          !  Broadcast costs to all processors
          !---------------------------------------------------------------------
          iopt = ppm_param_alloc_fit
          ldu(1) = ppm_nproc
          CALL ppm_alloc(sendcost,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_get_cost',    &
     &            'send buffer for costs SENDCOST',__LINE__,info)
              GOTO 9999
          ENDIF
          ! For MPI we also need an array starting at 1...
          CALL ppm_alloc(proccost,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_get_cost',    &
     &            'recv buffer for costs PROCCOST',__LINE__,info)
              GOTO 9999
          ENDIF
          sendcost(1:ppm_nproc) = cost(1)
          proccost(1:ppm_nproc) = 0.0_MK
#if   __KIND == __SINGLE_PRECISION
          CALL MPI_Alltoall(sendcost,1,MPI_REAL,proccost,1,MPI_REAL,      &
     &        ppm_comm,info)
#elif __KIND == __DOUBLE_PRECISION
          CALL MPI_Alltoall(sendcost,1,MPI_DOUBLE_PRECISION,proccost,1,   &
     &                      MPI_DOUBLE_PRECISION,ppm_comm,info)
#endif
          ! fill in the final array starting from 0
          proc_cost(0:(ppm_nproc-1)) = proccost(1:ppm_nproc)
          iopt = ppm_param_dealloc
          CALL ppm_alloc(proccost,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_get_cost',     &
     &            'recv buffer for costs PROCCOST',__LINE__,info)
          ENDIF
          CALL ppm_alloc(sendcost,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_get_cost',     &
     &            'send buffer for costs SENDCOST',__LINE__,info)
          ENDIF
#else
          proc_cost(0) = cost(1)
#endif
      ELSE
          !---------------------------------------------------------------------
          !  If this is not a ring, all processors know the costs of all the
          !  subs and we simply sum them up according to subs2proc
          !---------------------------------------------------------------------
          DO i=1,ppm_nsubs(topoid)
              proc = ppm_subs2proc(i,topoid)
              proc_cost(proc) = proc_cost(proc) + cost(i)
          ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Diagnostics
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          mincost = HUGE(cost(1))
          maxcost = 0.0_MK
          minproc = -1
          maxproc = -1
          DO i=0,ppm_nproc-1
              IF (proc_cost(i) .LT. mincost) THEN
                  mincost = proc_cost(i)
                  minproc = i
              ENDIF
              IF (proc_cost(i) .GT. maxcost) THEN
                  maxcost = proc_cost(i)
                  maxproc = i
              ENDIF
          ENDDO
          WRITE(mesg,'(A,F15.3,A,I6)') 'min cost: ',mincost,   &
     &        ' on processor ',minproc
          CALL ppm_write(ppm_rank,'ppm_get_cost',mesg,info)
          WRITE(mesg,'(A,F15.3,A,I6)') 'max cost: ',maxcost,   &
     &        ' on processor ',maxproc
          CALL ppm_write(ppm_rank,'ppm_get_cost',mesg,info)
      ENDIF
      IF (ppm_debug .GT. 1) THEN
          DO i=0,ppm_nproc-1
              WRITE(mesg,'(I4,A,F15.3)') i,' proc cost: ',proc_cost(i)
              CALL ppm_write(ppm_rank,'ppm_get_cost',  &
     &               mesg,info)
          ENDDO
          CALL ppm_write(ppm_rank,'ppm_get_cost',  &
     &           '----------------------------------------------',info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_get_cost',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_get_cost_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_get_cost_d
#endif
