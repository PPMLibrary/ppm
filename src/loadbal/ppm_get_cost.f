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
      SUBROUTINE ppm_get_cost_s(topoid,meshid,xp,Np,cost,proc_cost,info,&
      &                         vlist,nvlist,pcost)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_get_cost_d(topoid,meshid,xp,Np,cost,proc_cost,info,&
      &                         vlist,nvlist,pcost)
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
      USE ppm_module_mpi
      USE ppm_module_check_id
      USE ppm_module_topo_cost
      USE ppm_module_data_loadbal
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
      INTEGER, DIMENSION(:)  , OPTIONAL, INTENT(IN   ) :: nvlist
      !!! Number of particles with which ip has to interact. Index: ip.
      INTEGER, DIMENSION(:,:), OPTIONAL, INTENT(IN   ) :: vlist
      !!! Verlet list. First index: particles with which particle ip interacts.
      !!! Second index: ip. The second index only runs up to the
      !!! largest ip with non-zero nvlist. This is to save memory since the
      !!! last particles are the ghosts and they do not have a verlet list.
      !!! This is only allocated and returned if lstore is .TRUE.
      REAL(MK), DIMENSION(:  ), OPTIONAL, INTENT(IN) :: pcost
      !!! per-particle costs.
      !!! If it is not present but the nvlist and vlist are given,
      !!! pcost is computed by PPM using the relative communication cost
      !!! If it is not present AND nvlist and vlist are NOT given,
      !!! pcost is 1.0 for each particle
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      TYPE(ppm_t_topo), POINTER :: topo

#ifdef __MPI
      REAL(MK), DIMENSION(:), ALLOCATABLE :: sendcost
      REAL(MK), DIMENSION(:), ALLOCATABLE :: proccost
#endif
      REAL(MK)                          :: t0,mincost,maxcost,comp_cost
      REAL(MK)                          :: comm_cost

      INTEGER, DIMENSION(ppm_dim)       :: ldl,ldu
      INTEGER, DIMENSION(3,1), TARGET   :: ndummy
      INTEGER, DIMENSION(:,:), POINTER  :: nnodes
      INTEGER                           :: i,j,k,proc
      INTEGER                           :: iopt,minproc,maxproc

      CHARACTER(LEN=ppm_char) :: mesg
      CHARACTER(LEN=ppm_char) :: caller='ppm_get_cost'

      LOGICAL :: valid
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

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
      or_fail_alloc('costs per processor PROC_COST',ppm_error=ppm_error_fatal)

      proc_cost = 0.0_MK
      comp_cost = 1._MK
      comm_cost = comp_cost * 5._MK

      !-------------------------------------------------------------------------
      !  Sum up the costs of all subs of the processors
      !-------------------------------------------------------------------------
      SELECT CASE (topoid)
      CASE (ppm_param_topo_undefined)
        !-----------------------------------------------------------------------
        !  In this case we setup the cost vector ourself. As there are no
        !  subdomains, the size is 1
        !-----------------------------------------------------------------------
        iopt = ppm_param_alloc_fit
        ldu(1) = 1
        CALL ppm_alloc(cost,ldu,iopt,info)
        or_fail_alloc('costs per sub COST',ppm_error=ppm_error_fatal)

        !---------------------------------------------------------------------
        !  Advanced communication costs.
        !  One particle cost: comp_cost
        !  One neighbor on the SAME proc cost: comp_cost (aka interaction cost)
        !  One neighbor on ANOTHER proc cost:  comm_cost
        !---------------------------------------------------------------------
        IF (PRESENT(pcost) .AND. PRESENT(nvlist) .AND. PRESENT(vlist)) THEN
            DO i = 1, Np
                cost(1) = cost(1) + pcost(i)
                DO j=1,nvlist(i)
           !--------------------------------------------------------------------
           !  Computational cost of a neighboring particle on THIS proc is 1
           !  on another processor, it's comm_cost.
           !  If the ID of a particle is GREATER THAN Np, it's a ghost particle
           !  thus it belongs to another processor.
           !--------------------------------------------------------------------
                    IF (vlist(j,i) .LE. Np) THEN
                        cost(1) = cost(1) + comp_cost
                    ELSE ! ghost particle --> add comm_cost
                        cost(1) = cost(1) + comm_cost
                    ENDIF
                ENDDO
            ENDDO
        ELSEIF (PRESENT(pcost)) THEN
            DO i = 1, Np
                cost(1) = cost(1) + pcost(i)
            ENDDO
        ELSE
            DO i = 1, Np
                cost(1) = cost(1) + comp_cost
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
        !  For MPI we also need an array starting at 1...
        !---------------------------------------------------------------------
        ALLOCATE(sendcost(ppm_nproc),proccost(ppm_nproc),STAT=info)
        or_fail_alloc('send & recv buffer for costs SENDCOST',ppm_error=ppm_error_fatal)

        sendcost=cost(1)

#if   __KIND == __SINGLE_PRECISION
        CALL MPI_Alltoall(sendcost,1,MPI_REAL,proccost,1,MPI_REAL,      &
        &    ppm_comm,info)
#elif __KIND == __DOUBLE_PRECISION
        CALL MPI_Alltoall(sendcost,1,MPI_DOUBLE_PRECISION,proccost,1,   &
        &    MPI_DOUBLE_PRECISION,ppm_comm,info)
#endif
        ! fill in the final array starting from 0
        proc_cost(0:(ppm_nproc-1)) = proccost(1:ppm_nproc)

        DEALLOCATE(sendcost,proccost,STAT=info)
        or_fail_dealloc('send & recv buffer for costs PROCCOST')
#else
        proc_cost(0) = cost(1)
#endif

      CASE DEFAULT

        topo => ppm_topo(topoid)%t
        !-------------------------------------------------------------------------
        !  set pointer for mesh
        !-------------------------------------------------------------------------
        !TODO
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
            CALL ppm_topo_cost(xp,Np,topo%min_subs(:,:), &
            &    topo%max_subs(:,:),topo%nsubs,          &
            &    nnodes,cost,info,pcost)
#elif __KIND == __DOUBLE_PRECISION
            CALL ppm_topo_cost(xp,Np,topo%min_subd(:,:), &
            &    topo%max_subd(:,:),topo%nsubs,          &
            &    nnodes,cost,info,pcost)
#endif
        ELSE
#if   __KIND == __SINGLE_PRECISION
            CALL ppm_topo_cost(xp,Np,topo%min_subs(:,:), &
            &    topo%max_subs(:,:),topo%nsubs,          &
            &    nnodes,cost,info)
#elif __KIND == __DOUBLE_PRECISION
            CALL ppm_topo_cost(xp,Np,topo%min_subd(:,:), &
            &    topo%max_subd(:,:),topo%nsubs,          &
            &    nnodes,cost,info)
#endif
        ENDIF
        or_fail('Computing costs of subdomains failed')

        NULLIFY(nnodes)
        !---------------------------------------------------------------------
        !  Advanced communication costs.
        !  One particle cost: comp_cost
        !  One neighbor on the SAME proc cost: comp_cost (aka interaction cost)
        !  One neighbor on ANOTHER proc cost:  comm_cost
        !---------------------------------------------------------------------
        IF (PRESENT(pcost) .AND. PRESENT(nvlist) .AND. PRESENT(vlist)) THEN
            DO k = 1, topo%nsubs
                DO i = 1, Np
                    cost(k) = cost(k) + pcost(i)
                    DO j=1,nvlist(i)
               !--------------------------------------------------------------------
               !  Computational cost of a neighboring particle on THIS proc is 1
               !  on another processor, it's comm_cost.
               !  If the ID of a particle is GREATER THAN Np, it's a ghost particle
               !  thus it belongs to another processor.
               !--------------------------------------------------------------------
                        IF (vlist(j,i) .LE. Np) THEN
                            cost(k) = cost(k) + comp_cost
                        ELSE ! ghost particle --> add comm_cost
                            cost(k) = cost(k) + comm_cost
                        ENDIF
                    ENDDO
                ENDDO
            ENDDO
        ELSEIF (PRESENT(pcost)) THEN
            DO k = 1, topo%nsubs
                DO i = 1, Np
                    cost(k) = cost(k) + pcost(i)
                ENDDO
            ENDDO
        ELSE
            DO k = 1, topo%nsubs
                DO i = 1, Np
                    cost(k) = cost(k) + comp_cost
                ENDDO
            ENDDO
        ENDIF
        !---------------------------------------------------------------------
        !  All processors know the costs of all the
        !  subs and we simply sum them up according to subs2proc
        !---------------------------------------------------------------------
        DO i=1,topo%nsubs
            proc = topo%sub2proc(i)
            proc_cost(proc) = proc_cost(proc) + cost(i)
        ENDDO

      END SELECT

#if    __KIND == __SINGLE_PRECISION
      ppm_loadbal_subcosts => cost
      ppm_loadbal_proccosts= proc_cost(ppm_rank)
      print*,ppm_rank,ppm_loadbal_subcosts,ppm_loadbal_proccosts
#else
      ppm_loadbal_subcostd => cost
      ppm_loadbal_proccostd= proc_cost(ppm_rank)
      print*,ppm_rank,ppm_loadbal_subcostd,ppm_loadbal_proccostd
#endif

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
          WRITE(mesg,'(A,F15.3,A,I6)') 'min cost: ',mincost,' on processor ',minproc
          CALL ppm_write(ppm_rank,caller,mesg,info)
          WRITE(mesg,'(A,F15.3,A,I6)') 'max cost: ',maxcost,' on processor ',maxproc
          CALL ppm_write(ppm_rank,caller,mesg,info)
      ENDIF
      IF (ppm_debug .GT. 1) THEN
          DO i=0,ppm_nproc-1
              WRITE(mesg,'(I4,A,F15.3)') i,' proc cost: ',proc_cost(i)
              CALL ppm_write(ppm_rank,caller,mesg,info)
          ENDDO
          CALL ppm_write(ppm_rank,caller,  &
          &    '----------------------------------------------',info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
        IF (.NOT. ppm_initialized) THEN
           fail('Please call ppm_init first!',ppm_err_ppm_noinit,exit_point=8888)
        ENDIF
        IF (Np .GT. 0) THEN
            IF (SIZE(xp,2) .LT. Np) THEN
               fail('not enough particles contained in xp',exit_point=8888)
            ENDIF
            IF (SIZE(xp,1) .LT. ppm_dim) THEN
               fail('leading dimension of xp insufficient',exit_point=8888)
            ENDIF
        ENDIF
        IF (PRESENT(pcost)) THEN
          IF (SIZE(pcost,1) .LT. Np) THEN
             fail('pcost must be of at least length Np',exit_point=8888)
          ENDIF
        ENDIF
        IF ((Np .LE. 0) .AND. (topoid .EQ. ppm_param_topo_undefined)) THEN
           fail('null topology must use particles',exit_point=8888)
        ENDIF
        IF (topoid .NE. ppm_param_topo_undefined) THEN
           CALL ppm_check_meshid(topoid,meshid,valid,info)
           IF (.NOT. valid) THEN
              fail('meshid is invalid',exit_point=8888)
           ENDIF
        ENDIF
      8888 CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_get_cost_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_get_cost_d
#endif
