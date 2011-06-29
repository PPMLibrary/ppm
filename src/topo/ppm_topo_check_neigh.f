      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_topo_check
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
      SUBROUTINE ppm_topo_check_neigh_s(topoid,xp,ghost_req,Npart,has_one_way,topo_ok,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_topo_check_neigh_d(topoid,xp,ghost_req,Npart,has_one_way,topo_ok,info)
#endif
      !!! Checks if all neighbors are correct
      !!! For all particles find particles to interact with
      !!! Check if this box is neighbor or itself
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_check_id
      use ppm_module_alloc

      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xp
      !!! Particle locations
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: ghost_req
      !!! Particle ghost requirements
      !!! 1st: dim, 2nd: particleid
      INTEGER                 , INTENT(IN   ) :: Npart
      !!! Number of particles on this processor
      LOGICAL                 , INTENT(  OUT) :: topo_ok
      !!! Is the topology consistent
      LOGICAL                 , INTENT(  IN) ::  has_one_way
      !!! Is one way interaction between particles
      INTEGER                 , INTENT(IN   ) :: topoid
      !!! The topology to be checked
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK)                       :: t0
      INTEGER                        :: ipart,idom,i,j,k,l,ison, in, iid, iproc,ip,ip2,iproc2,k2,box1,box2
      LOGICAL                        :: valid
      INTEGER, DIMENSION(:), POINTER :: bcdef => NULL()
      TYPE(ppm_t_topo)     , POINTER :: topo  => NULL()

      ! Needed the collect data
      INTEGER, DIMENSION(3)             :: lda
      INTEGER, DIMENSION(:),    POINTER :: allnp => NULL()
      REAL(MK), DIMENSION(:,:,:),POINTER:: allxp => NULL()
      REAL(MK), DIMENSION(:,:,:),POINTER:: allghost_req => NULL()
      INTEGER, DIMENSION(:,:),  POINTER :: buf   => NULL()
      INTEGER, DIMENSION(:),  POINTER   :: req => NULL()
      INTEGER                           :: maxnp
      REAL(MK), DIMENSION(ppm_dim)      :: len_phys
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_topo_check_neigh',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Set topoid
      !-------------------------------------------------------------------------
      topo => ppm_topo(topoid)%t
      bcdef => topo%bcdef

      !--------------------------------------------------------------------
      ! Send all data to rank 0
      !--------------------------------------------------------------------
       ! first allocate the size info arrays
       lda(1) = ppm_nproc
       CALL ppm_alloc(allnp,lda,ppm_param_alloc_fit,info)
       IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_check_neigh',     &
 &            'failed to allocate allnp or allmp',__LINE__,info)
          GOTO 9999
       ENDIF
!       ! gather the np and mp at the root
       CALL mpi_gather(Npart,1,MPI_INTEGER,allnp,1,MPI_INTEGER,0,ppm_comm,info)
       IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_mpi_fail,'ppm_check_neigh',     &
 &            'failed to gather allnp or allmp',__LINE__,info)
          GOTO 9999
       ENDIF
       
       IF (ppm_rank.eq.0) THEN
          ! allocate allxp array
          maxnp = maxval(allnp)
          lda(1) = ppm_dim
          lda(2) = maxnp
          lda(3) = ppm_nproc
          CALL ppm_alloc(allxp,lda,ppm_param_alloc_fit,info)
          CALL ppm_alloc(allghost_req,lda,ppm_param_alloc_fit,info)
          CALL ppm_alloc(buf,lda,ppm_param_alloc_fit,info)
          IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_check_neigh',     &
    &            'failed to allocate allxp or buf',__LINE__,info)
             GOTO 9999
          ENDIF
          
          lda(1) = ppm_nproc
          CALL ppm_alloc(req,lda,ppm_param_alloc_fit,info)
          IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_check_neigh',     &
    &            'failed to allocate req',__LINE__,info)
             GOTO 9999
          ENDIF
          
          DO i=1,allnp(1)
             allxp(:,i,1) = xp(:,i)
             allghost_req(:,i,1) = ghost_req(:,i)
          ENDDO

          ! now let all procs communicate with rank 0
          DO iproc=1,ppm_nproc-1

             CALL mpi_recv(allxp(:,:,iproc+1),allnp(iproc+1)*ppm_dim,ppm_mpi_kind,iproc,  &
 &                              0,ppm_comm,MPI_STATUS_IGNORE,info)
             CALL mpi_recv(allghost_req(:,:,iproc+1),allnp(iproc+1)*ppm_dim,ppm_mpi_kind,iproc,  &
 &                              0,ppm_comm,MPI_STATUS_IGNORE,info)
             IF (info .NE. 0) THEN
                   info = ppm_error_fatal
                   CALL ppm_error(ppm_err_mpi_fail,'ppm_check_neigh',   &
 &                'failed to sendrecv xp',__LINE__,info)
                   GOTO 9999
             ENDIF
          ENDDO

       ELSE
          CALL mpi_send(xp,Npart*ppm_dim,ppm_mpi_kind,0,0,ppm_comm,info)
          CALL mpi_send(ghost_req,Npart*ppm_dim,ppm_mpi_kind,0,0,ppm_comm,info)
       ENDIF
      
      ! We have all particles now in allxp with maxindex allnp (dim,index,proc)
      topo_ok = .TRUE.

      ! include periodicity!!!!!!!!!!!!
   
      IF (ppm_rank.eq.0) THEN

         IF (ppm_dim .EQ. 2) THEN

#if    __KIND == __SINGLE_PRECISION
            len_phys(1) = topo%max_physs(1) - topo%min_physs(1)
            len_phys(2) = topo%max_physs(2) - topo%min_physs(2)
#elif  __KIND == __DOUBLE_PRECISION
            len_phys(1) = topo%max_physd(1) - topo%min_physd(1)
            len_phys(2) = topo%max_physd(2) - topo%min_physd(2)
#endif

            ! For all particles ip
            DO iproc = 1, ppm_nproc
               DO ip = 1, allnp(iproc)
                  
                  ! Find box of ip
                  DO j=1,topo%nsubs
                     idom = j
#if    __KIND == __SINGLE_PRECISION
                     IF (allxp(1,ip,iproc).GE.topo%min_subs(1,idom).AND.  &
      &                 allxp(1,ip,iproc).LE.topo%max_subs(1,idom).AND.  &
      &                 allxp(2,ip,iproc).GE.topo%min_subs(2,idom).AND.  &
      &                 allxp(2,ip,iproc).LE.topo%max_subs(2,idom)) THEN
                        !------------------------------------------------------------
                        !  In the non-periodic case, allow particles that are
                        !  exactly ON an upper EXTERNAL boundary.
                        !------------------------------------------------------------
                        IF((allxp(1,ip,iproc).LT.topo%max_subs(1,idom) .OR.  &
         &                (topo%subs_bc(2,idom).EQ.1           .AND.  &
         &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND.  &
         &                (allxp(2,ip,iproc).LT.topo%max_subs(2,idom) .OR.  &
         &                (topo%subs_bc(4,idom).EQ.1           .AND.  &
         &                bcdef(4).NE. ppm_param_bcdef_periodic))) THEN

                           !box found
                           box1 = idom
                           EXIT
                        ENDIF
                     ENDIF
#elif  __KIND == __DOUBLE_PRECISION
                     IF (allxp(1,ip,iproc).GE.topo%min_subd(1,idom).AND.  &
      &                 allxp(1,ip,iproc).LE.topo%max_subd(1,idom).AND.  &
      &                 allxp(2,ip,iproc).GE.topo%min_subd(2,idom).AND.  &
      &                 allxp(2,ip,iproc).LE.topo%max_subd(2,idom)) THEN
                        !------------------------------------------------------------
                        !  In the non-periodic case, allow particles that are
                        !  exactly ON an upper EXTERNAL boundary.
                        !------------------------------------------------------------
                        IF((allxp(1,ip,iproc).LT.topo%max_subd(1,idom) .OR.  &
         &                (topo%subs_bc(2,idom).EQ.1           .AND.  &
         &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND.  &
         &                (allxp(2,ip,iproc).LT.topo%max_subd(2,idom) .OR.  &
         &                (topo%subs_bc(4,idom).EQ.1           .AND.  &
         &                bcdef(4).NE. ppm_param_bcdef_periodic))) THEN

                           !box found

                           box1 = idom
                           EXIT
                        ENDIF
                     ENDIF
#endif
                  ENDDO

                  ! For all particles jp that are closer than ghost_req
                  DO iproc2 = 1, ppm_nproc
                        DO ip2 = 1, allnp(iproc2)

                           ! Check if they enclose each other
                           IF ( (( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
              &                         allghost_req(1,ip,iproc)) .AND. &
              &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
              &                         allghost_req(2,ip,iproc)) ) .AND. &
              &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
              &                         allghost_req(1,ip2,iproc2)) .AND. &
              &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
              &                         allghost_req(2,ip2,iproc2)) )) &
              &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
              &                         allghost_req(1,ip,iproc)) .AND. & 
              &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
              &                         allghost_req(2,ip,iproc)))) ) THEN

                              
            
                              ! Find box of ip2
                              DO j=1,topo%nsubs
                                 idom = j
#if    __KIND == __SINGLE_PRECISION
                                 IF (allxp(1,ip2,iproc2).GE.topo%min_subs(1,idom).AND.  &
                  &                 allxp(1,ip2,iproc2).LE.topo%max_subs(1,idom).AND.  &
                  &                 allxp(2,ip2,iproc2).GE.topo%min_subs(2,idom).AND.  &
                  &                 allxp(2,ip2,iproc2).LE.topo%max_subs(2,idom)) THEN
                                    !------------------------------------------------------------
                                    !  In the non-periodic case, allow particles that are
                                    !  exactly ON an upper EXTERNAL boundary.
                                    !------------------------------------------------------------
                                    IF((allxp(1,ip2,iproc2).LT.topo%max_subs(1,idom) .OR.  &
                     &                (topo%subs_bc(2,idom).EQ.1           .AND.  &
                     &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND.  &
                     &                (allxp(2,ip2,iproc2).LT.topo%max_subs(2,idom) .OR.  &
                     &                (topo%subs_bc(4,idom).EQ.1           .AND.  &
                     &                bcdef(4).NE. ppm_param_bcdef_periodic))) THEN

                                    !box found
                                    box2 = idom
                                    EXIT
                                 ENDIF
                              ENDIF
#elif  __KIND == __DOUBLE_PRECISION
                              IF (allxp(1,ip2,iproc2).GE.topo%min_subd(1,idom).AND.  &
               &                 allxp(1,ip2,iproc2).LE.topo%max_subd(1,idom).AND.  &
               &                 allxp(2,ip2,iproc2).GE.topo%min_subd(2,idom).AND.  &
               &                 allxp(2,ip2,iproc2).LE.topo%max_subd(2,idom)) THEN
                                 !------------------------------------------------------------
                                 !  In the non-periodic case, allow particles that are
                                 !  exactly ON an upper EXTERNAL boundary.
                                 !------------------------------------------------------------
                                 IF((allxp(1,ip2,iproc2).LT.topo%max_subd(1,idom) .OR.  &
                  &                (topo%subs_bc(2,idom).EQ.1           .AND.  &
                  &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND.  &
                  &                (allxp(2,ip2,iproc2).LT.topo%max_subd(2,idom) .OR.  &
                  &                (topo%subs_bc(4,idom).EQ.1           .AND.  &
                  &                bcdef(4).NE. ppm_param_bcdef_periodic))) THEN

                                    !box found
                                    box2 = idom
                                    EXIT
                                 ENDIF
                              ENDIF
#endif
                              ENDDO
                              ! Check if subid of box of jp is in array of boxes of ip

                              IF (.NOT. (box1 .EQ. box2)) THEN
               
                                 DO in = 1, topo%nneigh(box1)
                                    iid = topo%ineigh(in,box1)
               
                                    IF(box2 .EQ. iid) THEN
                                       ! we have the right neighbor
                                       EXIT
                                    ELSEIF (in .EQ. topo%nneigh(box1)) THEN
                                       ! no neighbor found until end of list
                                       ! found a particle that violates the topology
               print *, box1, box2, allxp(1,ip,iproc), allxp(2,ip,iproc), allxp(1,ip2,iproc2), allxp(2,ip2,iproc2)

               print *,' ', box1, topo%min_subd(1,box1),topo%min_subd(2,box1),topo%max_subd(1,box1),topo%max_subd(2,box1)
               print *,' ', box2, topo%min_subd(1,box2),topo%min_subd(2,box2),topo%max_subd(1,box2),topo%max_subd(2,box2)

                                       topo_ok = .FALSE.
                                       ! no need to search any further

                                      
                                       GOTO 9999
                                    ENDIF
                                 ENDDO

                              ENDIF

                           ENDIF
                           
                           ! Check if they enclose each other x dim
                           IF (topo%bcdef(1).EQ.ppm_param_bcdef_periodic) THEN

                               IF (( (( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
              &                         allghost_req(1,ip,iproc)) .AND. &
              &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
              &                         allghost_req(2,ip,iproc)) ) .AND. &
              &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
              &                         allghost_req(1,ip2,iproc2)) .AND. &
              &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
              &                         allghost_req(2,ip2,iproc2)) )) &
              &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
              &                         allghost_req(1,ip,iproc)) .AND. & 
              &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
              &                         allghost_req(2,ip,iproc)))) ) .OR. &
              &                 ( (( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
              &                         allghost_req(1,ip,iproc)) .AND. &
              &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
              &                         allghost_req(2,ip,iproc)) ) .AND. &
              &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
              &                         allghost_req(1,ip2,iproc2)) .AND. &
              &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
              &                         allghost_req(2,ip2,iproc2)) )) &
              &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
              &                         allghost_req(1,ip,iproc)) .AND. & 
              &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
              &                         allghost_req(2,ip,iproc))))))THEN

                              
            
                              ! Find box of ip2
                              DO j=1,topo%nsubs
                                 idom = j
#if    __KIND == __SINGLE_PRECISION
                                 IF (allxp(1,ip2,iproc2).GE.topo%min_subs(1,idom).AND.  &
                  &                 allxp(1,ip2,iproc2).LE.topo%max_subs(1,idom).AND.  &
                  &                 allxp(2,ip2,iproc2).GE.topo%min_subs(2,idom).AND.  &
                  &                 allxp(2,ip2,iproc2).LE.topo%max_subs(2,idom)) THEN
                                    !------------------------------------------------------------
                                    !  In the non-periodic case, allow particles that are
                                    !  exactly ON an upper EXTERNAL boundary.
                                    !------------------------------------------------------------
                                    IF((allxp(1,ip2,iproc2).LT.topo%max_subs(1,idom) .OR.  &
                     &                (topo%subs_bc(2,idom).EQ.1           .AND.  &
                     &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND.  &
                     &                (allxp(2,ip2,iproc2).LT.topo%max_subs(2,idom) .OR.  &
                     &                (topo%subs_bc(4,idom).EQ.1           .AND.  &
                     &                bcdef(4).NE. ppm_param_bcdef_periodic))) THEN

                                    !box found
                                    box2 = idom
                                    EXIT
                                 ENDIF
                              ENDIF
#elif  __KIND == __DOUBLE_PRECISION
                              IF (allxp(1,ip2,iproc2).GE.topo%min_subd(1,idom).AND.  &
               &                 allxp(1,ip2,iproc2).LE.topo%max_subd(1,idom).AND.  &
               &                 allxp(2,ip2,iproc2).GE.topo%min_subd(2,idom).AND.  &
               &                 allxp(2,ip2,iproc2).LE.topo%max_subd(2,idom)) THEN
                                 !------------------------------------------------------------
                                 !  In the non-periodic case, allow particles that are
                                 !  exactly ON an upper EXTERNAL boundary.
                                 !------------------------------------------------------------
                                 IF((allxp(1,ip2,iproc2).LT.topo%max_subd(1,idom) .OR.  &
                  &                (topo%subs_bc(2,idom).EQ.1           .AND.  &
                  &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND.  &
                  &                (allxp(2,ip2,iproc2).LT.topo%max_subd(2,idom) .OR.  &
                  &                (topo%subs_bc(4,idom).EQ.1           .AND.  &
                  &                bcdef(4).NE. ppm_param_bcdef_periodic))) THEN

                                    !box found
                                    box2 = idom
                                    EXIT
                                 ENDIF
                              ENDIF
#endif
                              ENDDO
                              ! Check if subid of box of jp is in array of boxes of ip

                              IF (.NOT. (box1 .EQ. box2)) THEN
               
                                 DO in = 1, topo%nneigh(box1)
                                    iid = topo%ineigh(in,box1)
               
                                    IF(box2 .EQ. iid) THEN
                                       ! we have the right neighbor
                                       EXIT
                                    ELSEIF (in .EQ. topo%nneigh(box1)) THEN
                                       ! no neighbor found until end of list
                                       ! found a particle that violates the topology
               print *, 'x periodic'
               print *, box1, box2, allxp(1,ip,iproc), allxp(2,ip,iproc), allxp(1,ip2,iproc2), allxp(2,ip2,iproc2)

               print *,' ', box1, topo%min_subd(1,box1),topo%min_subd(2,box1),topo%max_subd(1,box1),topo%max_subd(2,box1)
               print *,' ', box2, topo%min_subd(1,box2),topo%min_subd(2,box2),topo%max_subd(1,box2),topo%max_subd(2,box2)


                                       topo_ok = .FALSE.
                                       ! no need to search any further

                                      
                                       GOTO 9999
                                    ENDIF
                                 ENDDO

                              ENDIF

                           ENDIF


                           ENDIF

                           IF (topo%bcdef(3).EQ.ppm_param_bcdef_periodic) THEN

                              IF (( (( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
              &                         allghost_req(1,ip,iproc)) .AND. &
              &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
              &                         allghost_req(2,ip,iproc)) ) .AND. &
              &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
              &                         allghost_req(1,ip2,iproc2)) .AND. &
              &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
              &                         allghost_req(2,ip2,iproc2)) )) &
              &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
              &                         allghost_req(1,ip,iproc)) .AND. & 
              &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
              &                         allghost_req(2,ip,iproc)))) ) .OR. &
              &                 ( (( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
              &                         allghost_req(1,ip,iproc)) .AND. &
              &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
              &                         allghost_req(2,ip,iproc)) ) .AND. &
              &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
              &                         allghost_req(1,ip2,iproc2)) .AND. &
              &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
              &                         allghost_req(2,ip2,iproc2)) )) &
              &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
              &                         allghost_req(1,ip,iproc)) .AND. & 
              &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
              &                         allghost_req(2,ip,iproc))))))THEN

                              
            
                              ! Find box of ip2
                              DO j=1,topo%nsubs
                                 idom = j
#if    __KIND == __SINGLE_PRECISION
                                 IF (allxp(1,ip2,iproc2).GE.topo%min_subs(1,idom).AND.  &
                  &                 allxp(1,ip2,iproc2).LE.topo%max_subs(1,idom).AND.  &
                  &                 allxp(2,ip2,iproc2).GE.topo%min_subs(2,idom).AND.  &
                  &                 allxp(2,ip2,iproc2).LE.topo%max_subs(2,idom)) THEN
                                    !------------------------------------------------------------
                                    !  In the non-periodic case, allow particles that are
                                    !  exactly ON an upper EXTERNAL boundary.
                                    !------------------------------------------------------------
                                    IF((allxp(1,ip2,iproc2).LT.topo%max_subs(1,idom) .OR.  &
                     &                (topo%subs_bc(2,idom).EQ.1           .AND.  &
                     &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND.  &
                     &                (allxp(2,ip2,iproc2).LT.topo%max_subs(2,idom) .OR.  &
                     &                (topo%subs_bc(4,idom).EQ.1           .AND.  &
                     &                bcdef(4).NE. ppm_param_bcdef_periodic))) THEN

                                    !box found
                                    box2 = idom
                                    EXIT
                                 ENDIF
                              ENDIF
#elif  __KIND == __DOUBLE_PRECISION
                              IF (allxp(1,ip2,iproc2).GE.topo%min_subd(1,idom).AND.  &
               &                 allxp(1,ip2,iproc2).LE.topo%max_subd(1,idom).AND.  &
               &                 allxp(2,ip2,iproc2).GE.topo%min_subd(2,idom).AND.  &
               &                 allxp(2,ip2,iproc2).LE.topo%max_subd(2,idom)) THEN
                                 !------------------------------------------------------------
                                 !  In the non-periodic case, allow particles that are
                                 !  exactly ON an upper EXTERNAL boundary.
                                 !------------------------------------------------------------
                                 IF((allxp(1,ip2,iproc2).LT.topo%max_subd(1,idom) .OR.  &
                  &                (topo%subs_bc(2,idom).EQ.1           .AND.  &
                  &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND.  &
                  &                (allxp(2,ip2,iproc2).LT.topo%max_subd(2,idom) .OR.  &
                  &                (topo%subs_bc(4,idom).EQ.1           .AND.  &
                  &                bcdef(4).NE. ppm_param_bcdef_periodic))) THEN

                                    !box found
                                    box2 = idom
                                    EXIT
                                 ENDIF
                              ENDIF
#endif
                              ENDDO
                              ! Check if subid of box of jp is in array of boxes of ip

                              IF (.NOT. (box1 .EQ. box2)) THEN
               
                                 DO in = 1, topo%nneigh(box1)
                                    iid = topo%ineigh(in,box1)
               
                                    IF(box2 .EQ. iid) THEN
                                       ! we have the right neighbor
                                       EXIT
                                    ELSEIF (in .EQ. topo%nneigh(box1)) THEN
                                       ! no neighbor found until end of list
                                       ! found a particle that violates the topology
               print *, 'y periodic'
               print *, box1, box2, allxp(1,ip,iproc), allxp(2,ip,iproc), allxp(1,ip2,iproc2), allxp(2,ip2,iproc2)

               print *,' ', box1, topo%min_subd(1,box1),topo%min_subd(2,box1),topo%max_subd(1,box1),topo%max_subd(2,box1)
               print *,' ', box2, topo%min_subd(1,box2),topo%min_subd(2,box2),topo%max_subd(1,box2),topo%max_subd(2,box2)
     

                                       topo_ok = .FALSE.
                                       ! no need to search any further

                                      
                                       GOTO 9999
                                    ENDIF
                                 ENDDO

                              ENDIF

                           ENDIF

                           ENDIF

                           IF (topo%bcdef(1).EQ.ppm_param_bcdef_periodic .AND. topo%bcdef(3).EQ.ppm_param_bcdef_periodic) THEN

                              

                              IF (( (( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
              &                         allghost_req(1,ip,iproc)) .AND. &
              &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
              &                         allghost_req(2,ip,iproc)) ) .AND. &
              &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
              &                         allghost_req(1,ip2,iproc2)) .AND. &
              &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
              &                         allghost_req(2,ip2,iproc2)) )) &
              &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
              &                         allghost_req(1,ip,iproc)) .AND. & 
              &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
              &                         allghost_req(2,ip,iproc)))) ) .OR. &
              &                 ( (( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
              &                         allghost_req(1,ip,iproc)) .AND. &
              &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
              &                         allghost_req(2,ip,iproc)) ) .AND. &
              &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
              &                         allghost_req(1,ip2,iproc2)) .AND. &
              &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
              &                         allghost_req(2,ip2,iproc2)) )) &
              &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
              &                         allghost_req(1,ip,iproc)) .AND. & 
              &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
              &                         allghost_req(2,ip,iproc))))))THEN

                              ! Find box of ip2
                              DO j=1,topo%nsubs
                                 idom = j
#if    __KIND == __SINGLE_PRECISION
                                 IF (allxp(1,ip2,iproc2).GE.topo%min_subs(1,idom).AND.  &
                  &                 allxp(1,ip2,iproc2).LE.topo%max_subs(1,idom).AND.  &
                  &                 allxp(2,ip2,iproc2).GE.topo%min_subs(2,idom).AND.  &
                  &                 allxp(2,ip2,iproc2).LE.topo%max_subs(2,idom)) THEN
                                    !------------------------------------------------------------
                                    !  In the non-periodic case, allow particles that are
                                    !  exactly ON an upper EXTERNAL boundary.
                                    !------------------------------------------------------------
                                    IF((allxp(1,ip2,iproc2).LT.topo%max_subs(1,idom) .OR.  &
                     &                (topo%subs_bc(2,idom).EQ.1           .AND.  &
                     &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND.  &
                     &                (allxp(2,ip2,iproc2).LT.topo%max_subs(2,idom) .OR.  &
                     &                (topo%subs_bc(4,idom).EQ.1           .AND.  &
                     &                bcdef(4).NE. ppm_param_bcdef_periodic))) THEN

                                    !box found
                                    box2 = idom
                                    EXIT
                                 ENDIF
                              ENDIF
#elif  __KIND == __DOUBLE_PRECISION
                              IF (allxp(1,ip2,iproc2).GE.topo%min_subd(1,idom).AND.  &
               &                 allxp(1,ip2,iproc2).LE.topo%max_subd(1,idom).AND.  &
               &                 allxp(2,ip2,iproc2).GE.topo%min_subd(2,idom).AND.  &
               &                 allxp(2,ip2,iproc2).LE.topo%max_subd(2,idom)) THEN
                                 !------------------------------------------------------------
                                 !  In the non-periodic case, allow particles that are
                                 !  exactly ON an upper EXTERNAL boundary.
                                 !------------------------------------------------------------
                                 IF((allxp(1,ip2,iproc2).LT.topo%max_subd(1,idom) .OR.  &
                  &                (topo%subs_bc(2,idom).EQ.1           .AND.  &
                  &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND.  &
                  &                (allxp(2,ip2,iproc2).LT.topo%max_subd(2,idom) .OR.  &
                  &                (topo%subs_bc(4,idom).EQ.1           .AND.  &
                  &                bcdef(4).NE. ppm_param_bcdef_periodic))) THEN

                                    !box found
                                    box2 = idom
                                    EXIT
                                 ENDIF
                              ENDIF
#endif
                              ENDDO
                              ! Check if subid of box of jp is in array of boxes of ip

                              IF (.NOT. (box1 .EQ. box2)) THEN
               
                                 DO in = 1, topo%nneigh(box1)
                                    iid = topo%ineigh(in,box1)
               
                                    IF(box2 .EQ. iid) THEN
                                       ! we have the right neighbor
                                       EXIT
                                    ELSEIF (in .EQ. topo%nneigh(box1)) THEN
                                       ! no neighbor found until end of list
                                       ! found a particle that violates the topology
               print *, 'xy periodic'
               print *, box1, box2, allxp(1,ip,iproc), allxp(2,ip,iproc), allxp(1,ip2,iproc2), allxp(2,ip2,iproc2)

               print *,' ', box1, topo%min_subd(1,box1),topo%min_subd(2,box1),topo%max_subd(1,box1),topo%max_subd(2,box1)
               print *,' ', box2, topo%min_subd(1,box2),topo%min_subd(2,box2),topo%max_subd(1,box2),topo%max_subd(2,box2)
     

                                       topo_ok = .FALSE.
                                       ! no need to search any further

                                      
                                       GOTO 9999
                                    ENDIF
                                 ENDDO

                              ENDIF

                           ENDIF

                           ENDIF

                          
                     ENDDO
                  ENDDO

               ENDDO
            ENDDO

         ELSEIF(ppm_dim .EQ. 3) THEN


         ENDIF

         CALL ppm_alloc(allxp,lda,ppm_param_dealloc,info)

      ENDIF
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_topo_check_neigh',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
         IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_topo_check_neigh',       &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (Npart .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_check_neigh',         &
     &            'Npart must be >= 0',__LINE__,info)
              GOTO 8888
          ENDIF
          CALL ppm_check_topoid(topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_check_neigh', &
     &             'topoid invalid',&
     &              __LINE__, info)
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_topo_check_neigh_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_topo_check_neigh_d
#endif
