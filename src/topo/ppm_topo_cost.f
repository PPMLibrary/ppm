      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_topo_cost
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
      SUBROUTINE ppm_topo_cost_s(xp,Np,min_sub,max_sub,nsubs,nnodes,cost,   &
     &              info,pcost)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_topo_cost_d(xp,Np,min_sub,max_sub,nsubs,nnodes,cost,   &
     &              info,pcost)
#endif
      !!! This routine calculates the computational cost of
      !!! each subdomain based on either the particles (if there are any)
      !!! or the mesh points. The per-particle costs are summed up if present.
      !!! If not, a unity cost is assumed for each particle.
      !!!
      !!! [NOTE]
      !!! If `Np < 0`, cost is computed based on particles. If
      !!! `SIZE(nnodes,2) >= nsubs`, cost based on mesh points
      !!! is computed and added to the particle cost. If  neither particles
      !!! nor mesh are given, cost of a subdomain is equal to its volume.
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
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xp
      !!! Particle positions
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: min_sub
      !!! Min. extent of the subdomains
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: max_sub
      !!! Max. extent of the subdomains
      REAL(MK), DIMENSION(:  ), POINTER       :: cost
      !!! Aggregate cost for each subdomain
      REAL(MK), DIMENSION(:  ), OPTIONAL, INTENT(IN) :: pcost
      !!! Per-particle costs. OPTIONAL. If  not present, a cost of 1 per
      !!! particle is assumed.
      INTEGER , DIMENSION(:,:), INTENT(IN   ) :: nnodes
      !!! Number of mesh nodes in each direction (first index) on each
      !!! sub (second index). If Np is .LE. 0, this is taken to compute the
      !!! cost. If SIZE(nnodes,2).LT.nsubs, and there are no particles, cost
      !!! is computed as geometric volume.
      INTEGER                 , INTENT(IN   ) :: Np
      !!! Number of particles. Set to .LE. 0 if mesh-based costs are desired.
      INTEGER                 , INTENT(IN   ) :: nsubs
      !!! Total number of subdomains
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                          :: t0
#ifdef __MPI
      REAL(MK), DIMENSION(:), POINTER   :: costsum => NULL()
#endif
      INTEGER, DIMENSION(ppm_dim)       :: ldu,Nc
      REAL(MK),DIMENSION(ppm_dim)       :: len_sub
      INTEGER                           :: i,ipart,j
      INTEGER                           :: iopt,nlist1,nlist2,idom
      INTEGER, DIMENSION(:), POINTER    :: ilist1 => NULL()
      INTEGER, DIMENSION(:), POINTER    :: ilist2 => NULL()
      INTEGER, DIMENSION(:), POINTER    :: ilist3 => NULL()
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
        CALL check
        IF (info .NE. 0) GOTO 9999
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
             CALL MPI_AllReduce(cost,costsum,nsubs,MPI_REAL,MPI_SUM,    &
     &          ppm_comm,info)
#elif __KIND == __DOUBLE_PRECISION
             CALL MPI_AllReduce(cost,costsum,nsubs,MPI_DOUBLE_PRECISION,&
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
                 CALL ppm_error(ppm_err_dealloc,'ppm_topo_cost',  &
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
      CONTAINS
      SUBROUTINE check
          IF (nsubs .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_cost',  &
     &            'nsubs must be > 0',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (SIZE(min_sub,2) .LT. nsubs .OR. SIZE(max_sub,2) .LT. nsubs) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_cost',  &
     &            'not enough subs specified in min_sub, max_sub',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (SIZE(min_sub,1).LT.ppm_dim.OR.SIZE(max_sub,1).LT.ppm_dim) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_cost',  &
     &            'leading dimension in min_sub, max_sub wrong',__LINE__,info)
              GOTO 8888
          ENDIF
          DO i=1,nsubs
              DO j=1,ppm_dim
                  IF (max_sub(j,i) .LT. min_sub(j,i)) THEN
                      info = ppm_error_error
                      CALL ppm_error(ppm_err_argument,'ppm_topo_cost',  &
     &                    'max_sub must always be >= min_sub',__LINE__,info)
                      GOTO 8888
                  ENDIF
              ENDDO
          ENDDO
          IF (Np .GT. 0) THEN
              IF (SIZE(xp,2) .LT. Np) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_topo_cost',  &
     &                'not enough particles contained in xp',__LINE__,info)
                  GOTO 8888
              ENDIF
              IF (SIZE(xp,1) .LT. ppm_dim) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_topo_cost',  &
     &                'leading dimension of xp insufficient',__LINE__,info)
                  GOTO 8888
              ENDIF
          ENDIF
          IF (PRESENT(pcost)) THEN
              IF (SIZE(pcost,1) .LT. Np) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_topo_cost',  &
     &                'pcost must be of at least length Np',__LINE__,info)
                  GOTO 8888
              ENDIF
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_topo_cost_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_topo_cost_d
#endif
