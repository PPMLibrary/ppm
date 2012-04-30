      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_get_cost
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
      SUBROUTINE ppm_get_cost_s(topoid,meshid,xp,Np,cost,proc_cost,info,pcost)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_get_cost_d(topoid,meshid,xp,Np,cost,proc_cost,info,pcost)
#endif
      !!! This routine calculates the computational cost of each subdomain
      !!! and each processor.
      !!!
      !!! [NOTE]
      !!! If Np > 0, cost is computed based on particles. If mesh_id > -1,
      !!! cost based on mesh points is computed and added to the particle
      !!! cost. If neither particles nor mesh are given, cost of a
      !!! subdomain is equal to its volume.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      USE ppm_module_check_id
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
#ifdef __MPI
       INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xp
      !!! Particle positions
      REAL(MK), DIMENSION(:  ), POINTER       :: cost
      !!! Aggregate cost for each subdomain
      REAL(MK), DIMENSION(:  ), POINTER       :: proc_cost
      !!! Aggregate cost for each processor
      INTEGER                 , INTENT(IN   ) :: Np
      !!! Number of particles. Set to <= 0 if mesh-based costs are desired.
      INTEGER                 , INTENT(IN   ) :: topoid
      !!! Topology ID for which to compute the cost.
      !!! If `ppm_param_topo_undefined` all particle positions on the processor
      !!! are considered to compute the cost
      INTEGER                 , INTENT(IN   ) :: meshid
      !!! mesh ID for which to compute the cost. 
      !!! If -1 is passed, only particles are considered and the mesh is
      !!! ignored.                                                             +
      !!! If there are no particles and mesh_id is -1, the costs are
      !!! computed based on the geometry: cost(sub) = volume(sub).
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      REAL(MK), DIMENSION(:  ), OPTIONAL, INTENT(IN) :: pcost
      !!! per-particle costs.
      !!! If not present, a cost of 1.0 per particle is assumed.
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK)                          :: t0,mincost,maxcost
      INTEGER, DIMENSION(ppm_dim)       :: ldl,ldu
      INTEGER, DIMENSION(3,1), TARGET   :: ndummy
      INTEGER, DIMENSION(:,:), POINTER  :: nnodes => NULL()
      INTEGER                           :: i,proc
      INTEGER                           :: iopt,minproc,maxproc
      LOGICAL                           :: valid
      CHARACTER(LEN=ppm_char)           :: mesg
#ifdef __MPI
      REAL(MK), DIMENSION(:), POINTER   :: sendcost => NULL()
      REAL(MK), DIMENSION(:), POINTER   :: proccost => NULL()
#endif
      TYPE(ppm_t_topo)      , POINTER   :: topo     => NULL()
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
        CALL check
        IF (info .NE. 0) GOTO 9999
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
          CALL ppm_error(ppm_err_alloc,'ppm_get_cost',              &
     &        'costs per processor PROC_COST',__LINE__,info)
          GOTO 9999
      ENDIF
      proc_cost = 0.0_MK



      !-------------------------------------------------------------------------
      !  Sum up the costs of all subs of the processors
      !-------------------------------------------------------------------------
      IF (topoid .EQ. ppm_param_topo_undefined) THEN

        !-----------------------------------------------------------------------
        !  In this case we setup the cost vector ourself. As there are no
        !  subdomains, the size is 1
        !-----------------------------------------------------------------------
        iopt = ppm_param_alloc_fit
        ldu(1) = 1
        CALL ppm_alloc(cost,ldu,iopt,info)
        IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_topo_cost',    &
     &          'costs per sub COST',__LINE__,info)
            GOTO 9999
        ENDIF

        IF (PRESENT(pcost)) THEN
            DO i = 1, Np
                cost(1) = cost(1) + pcost(i)
            ENDDO
        ELSE
            DO i = 1, Np
                cost(1) = cost(1) + 1.0_MK
            ENDDO
        ENDIF

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
     &          'send buffer for costs SENDCOST',__LINE__,info)
            GOTO 9999
        ENDIF
        ! For MPI we also need an array starting at 1...
        CALL ppm_alloc(proccost,ldu,iopt,info)
        IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_get_cost',    &
     &          'recv buffer for costs PROCCOST',__LINE__,info)
            GOTO 9999
        ENDIF
        sendcost(1:ppm_nproc) = cost(1)
        proccost(1:ppm_nproc) = 0.0_MK
#if   __KIND == __SINGLE_PRECISION
        CALL MPI_Alltoall(sendcost,1,MPI_REAL,proccost,1,MPI_REAL,      &
     &      ppm_comm,info)
#elif __KIND == __DOUBLE_PRECISION
        CALL MPI_Alltoall(sendcost,1,MPI_DOUBLE_PRECISION,proccost,1,   &
     &                    MPI_DOUBLE_PRECISION,ppm_comm,info)
#endif
        ! fill in the final array starting from 0
        proc_cost(0:(ppm_nproc-1)) = proccost(1:ppm_nproc)
        iopt = ppm_param_dealloc
        CALL ppm_alloc(proccost,ldu,iopt,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_dealloc,'ppm_get_cost',     &
     &          'recv buffer for costs PROCCOST',__LINE__,info)
        ENDIF
        CALL ppm_alloc(sendcost,ldu,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_error
            CALL ppm_error(ppm_err_dealloc,'ppm_get_cost',     &
     &          'send buffer for costs SENDCOST',__LINE__,info)
        ENDIF
#else
        proc_cost(0) = cost(1)
#endif
      ELSE

        topo => ppm_topo(topoid)%t
        !-------------------------------------------------------------------------
        !  set pointer for mesh
        !-------------------------------------------------------------------------
        !FIXME: we have to account for the data stored on the mesh, etc...
        ! should be doable with the new DS
        ndummy(1:3,1) = 0
        !IF (meshid .GT. 0) THEN
            !nnodes => topo%mesh(meshid)%nnodes
        !ELSE
            nnodes => ndummy
        !ENDIF
        !-------------------------------------------------------------------------
        !  Determine the subdomain costs either based on particles or mesh
        !  points
        !-------------------------------------------------------------------------
        IF (PRESENT(pcost)) THEN
#if   __KIND == __SINGLE_PRECISION
            CALL ppm_topo_cost(xp,Np,topo%min_subs(:,:),       &
     &          topo%max_subs(:,:),topo%nsubs,         &
     &          nnodes,cost,info,pcost)
#elif __KIND == __DOUBLE_PRECISION
            CALL ppm_topo_cost(xp,Np,topo%min_subd(:,:),       &
     &          topo%max_subd(:,:),topo%nsubs,         &
     &          nnodes,cost,info,pcost)
#endif
        ELSE
#if   __KIND == __SINGLE_PRECISION
            CALL ppm_topo_cost(xp,Np,topo%min_subs(:,:),       &
     &          topo%max_subs(:,:),topo%nsubs,         &
     &          nnodes,cost,info)
#elif __KIND == __DOUBLE_PRECISION
            CALL ppm_topo_cost(xp,Np,topo%min_subd(:,:),       &
     &          topo%max_subd(:,:),topo%nsubs,         &
     &          nnodes,cost,info)
#endif
        ENDIF
        IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,'ppm_get_cost',         &
     &              'Computing costs of subdomains failed',__LINE__,info)
            GOTO 9999
        ENDIF

        NULLIFY(nnodes)
        !---------------------------------------------------------------------
        !  All processors know the costs of all the
        !  subs and we simply sum them up according to subs2proc
        !---------------------------------------------------------------------
        DO i=1,topo%nsubs
            proc = topo%sub2proc(i)
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
      CONTAINS
      SUBROUTINE check
        IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_get_cost',    &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 8888
        ENDIF
        IF (Np .GT. 0) THEN
            IF (SIZE(xp,2) .LT. Np) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_get_cost',  &
     &              'not enough particles contained in xp',__LINE__,info)
                GOTO 8888
            ENDIF
            IF (SIZE(xp,1) .LT. ppm_dim) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_get_cost',  &
     &              'leading dimension of xp insufficient',__LINE__,info)
                GOTO 8888
            ENDIF
        ENDIF
        IF (PRESENT(pcost)) THEN
          IF (SIZE(pcost,1) .LT. Np) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_get_cost',  &
     &            'pcost must be of at least length Np',__LINE__,info)
              GOTO 8888
          ENDIF
        ENDIF
        IF ((Np .LE. 0) .AND. (topoid .EQ. ppm_param_topo_undefined)) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_get_cost',  &
     &            'null topology must use particles',__LINE__,info)
             GOTO 8888
        ENDIF
        IF (topoid .NE. ppm_param_topo_undefined) THEN
           CALL ppm_check_meshid(topoid,meshid,valid,info)
           IF (.NOT. valid) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_get_cost',  &
     &             'meshid is invalid',__LINE__,info)
               GOTO 8888
           ENDIF
        ENDIF
 8888   CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_get_cost_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_get_cost_d
#endif
