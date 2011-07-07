      !-------------------------------------------------------------------------
      !  Subroutine   :                   map_check_ghosts
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
      SUBROUTINE ppm_map_check_ghosts_s(topoid,xp,ghost_req,Npart,mp,has_one_way,map_ok,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_check_ghosts_d(topoid,xp,ghost_req,Npart,mp,has_one_way,map_ok,info)
#endif
      !!! Checks if all ghosts are present that we need
      !!! For all particles find particles to interact with, which are on other proc
      !!! Check if this particle we have as ghost
! TODO: symmetry and 3d
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
      INTEGER                 , INTENT(IN   ) :: mp
      !!! Number of particles on this processor plus ghost particles
      LOGICAL                 , INTENT(  OUT) :: map_ok
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
      REAL(MK)                       :: t0,lmyeps
      INTEGER                        :: ipart,idom,i,j,k,l,ison, in, iid, iproc,ip,ip2,iproc2,k2,box1,box2
      LOGICAL                        :: valid
      INTEGER, DIMENSION(:), POINTER :: bcdef => NULL()
      TYPE(ppm_t_topo)     , POINTER :: topo  => NULL()

      ! Needed the collect data
      INTEGER, DIMENSION(3)             :: lda
      LOGICAL, DIMENSION(:),    POINTER :: have_xp => NULL()
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
#if   __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_check_ghosts',t0,info)

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
      ! Send all data to all procs, first to 0 and then to all others
      !--------------------------------------------------------------------
       ! first allocate the size info arrays
       lda(1) = ppm_nproc
       CALL ppm_alloc(allnp,lda,ppm_param_alloc_fit,info)
       IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_check_ghosts',     &
 &            'failed to allocate allnp or allmp',__LINE__,info)
          GOTO 9999
       ENDIF
!       ! gather the np and mp at the root and send back to all
       CALL mpi_allgather(Npart,1,MPI_INTEGER,allnp,1,MPI_INTEGER,ppm_comm,info)
       IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_mpi_fail,'ppm_check_ghosts',     &
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
      
      !--------------------------------------------------------------------
      ! Send this array holding all particles to all other procs
      !--------------------------------------------------------------------
      IF (ppm_rank .EQ. 0) THEN

         ! send to all other procs
         DO i=1,ppm_nproc-1
            DO j = 1,ppm_nproc
               CALL mpi_send(allxp(:,:,j),allnp(j)*ppm_dim,ppm_mpi_kind,i,0,ppm_comm,info)
               CALL mpi_send(allghost_req(:,:,j),allnp(j)*ppm_dim,ppm_mpi_kind,i,0,ppm_comm,info)
            ENDDO
         ENDDO

      ELSE

         ! allocate and receive
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

         DO j = 1,ppm_nproc

            CALL mpi_recv(allxp(:,:,j),allnp(j)*ppm_dim,ppm_mpi_kind,0,  &
   &                              0,ppm_comm,MPI_STATUS_IGNORE,info)
            CALL mpi_recv(allghost_req(:,:,j),allnp(j)*ppm_dim,ppm_mpi_kind,0,  &
   &                              0,ppm_comm,MPI_STATUS_IGNORE,info)

               
         ENDDO

         ENDIF

!       print *, ppm_rank, maxnp, mp, mp-Npart

      ! We have all particles now in allxp with maxindex allnp (dim,index,proc)
      map_ok = .TRUE.

      ! allocate a have already checked
      lda(1) = mp
      CALL ppm_alloc(have_xp,lda,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_check_ghosts',     &
 &           'failed to allocate allnp or allmp',__LINE__,info)
         GOTO 9999
      ENDIF

      ! init have_xp
      DO i=1,mp
         have_xp(i) = .FALSE.
      ENDDO
   
      IF (ppm_dim .EQ. 2) THEN

#if    __KIND == __SINGLE_PRECISION
            len_phys(1) = topo%max_physs(1) - topo%min_physs(1)
            len_phys(2) = topo%max_physs(2) - topo%min_physs(2)
#elif  __KIND == __DOUBLE_PRECISION
            len_phys(1) = topo%max_physd(1) - topo%min_physd(1)
            len_phys(2) = topo%max_physd(2) - topo%min_physd(2)
#endif

         ! For all particles ip on this proc
         iproc = ppm_rank+1
            DO ip = 1, allnp(iproc)
               
               ! For all particles on other procs check if we have it in our ghosts particles
               DO iproc2 = 1, ppm_nproc
                     !IF (.NOT. (iproc2 .EQ. iproc)) THEN

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

                              
                                 ! We need to find ip2 in the ghost of iproc
                                 DO in = 1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
!                                        print *, ' error ', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
!          &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO

                           ENDIF

                           !if we have periodicity, also check that
                           IF (topo%bcdef(1).EQ.ppm_param_bcdef_periodic) THEN

                              IF ( (( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
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
               &                         allghost_req(2,ip,iproc))))) THEN
 
                                 ! We need to find ip2 in the ghost of iproc
                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)+len_phys(1)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
!                                        print *, ' error x', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
!          &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO


                              ENDIF

                              IF ( (( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
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
               &                         allghost_req(2,ip,iproc))))) THEN

                              
                                 ! We need to find ip2 in the ghost of iproc
                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)-len_phys(1)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
!                                        print *, ' error x', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
!          &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO

                              ENDIF



                           ENDIF

                           IF (topo%bcdef(3).EQ.ppm_param_bcdef_periodic) THEN

                              IF ( (( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
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
               &                         allghost_req(2,ip,iproc))))) THEN
                                 ! We need to find ip2 in the ghost of iproc
                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)+len_phys(2)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
!                                        print *, ' error y', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
!          &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO

                              ENDIF  
                              IF ( (( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
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
               &                         allghost_req(2,ip,iproc)))) ) THEN
                                 ! We need to find ip2 in the ghost of iproc
                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)-len_phys(2)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
!                                        print *, ' error y', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
!          &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO

                              ENDIF

                           ENDIF

                           IF (topo%bcdef(1).EQ.ppm_param_bcdef_periodic .AND. topo%bcdef(3).EQ.ppm_param_bcdef_periodic) THEN
                              IF  ((( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
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
               &                         allghost_req(2,ip,iproc))))) THEN
                                 ! We need to find ip2 in the ghost of iproc

                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)+len_phys(1)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)+len_phys(2)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
!                                        print *, ' error xy', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
!          &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO
                              ENDIF

                           IF  ((( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
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
               &                         allghost_req(2,ip,iproc))))) THEN
                              
                                 ! We need to find ip2 in the ghost of iproc
                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)-len_phys(1)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)-len_phys(2)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
!                                        print *, ' error xy', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
!          &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO

                              ENDIF

                              IF  ((( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc)) ) .AND. &
               &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
               &                         allghost_req(1,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
               &                         allghost_req(2,ip2,iproc2)) )) &
               &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. & 
               &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc))))) THEN
                                 ! We need to find ip2 in the ghost of iproc

                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)+len_phys(1)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)-len_phys(2)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
!                                        print *, ' error xy', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
!          &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO
                              ENDIF

                           IF  ((( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc)) ) .AND. &
               &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
               &                         allghost_req(1,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
               &                         allghost_req(2,ip2,iproc2)) )) &
               &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. & 
               &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc))))) THEN

                              
                                 ! We need to find ip2 in the ghost of iproc
                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)-len_phys(1)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)+len_phys(2)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
!                                        print *, ' errmak(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO

                              ENDIF

                           ENDIF

                           
                        ENDDO

                  !ENDIF

               ENDDO

            ENDDO

         ! plot the number of particles we do not need, but have in the ghosts
         IF(ppm_debug .GT. 0) THEN
            j = 0
            DO i = 1,Npart
               IF(have_xp(i) .EQV. .FALSE.) THEN
                  j = j+1
               ENDIF
            ENDDO

            IF (mp-Npart .GT. 0) THEN
               print *,'[', ppm_rank,']', int(100*float(j)/float((mp-Npart))) , '% of ghost particles not needed!'
            ELSE
               print *,'[', ppm_rank,']', 'has no ghost particles and does not need any!'
            ENDIF
         ENDIF
         
      ELSEIF(ppm_dim .EQ. 3) THEN

#if    __KIND == __SINGLE_PRECISION
            len_phys(1) = topo%max_physs(1) - topo%min_physs(1)
            len_phys(2) = topo%max_physs(2) - topo%min_physs(2)
            len_phys(3) = topo%max_physs(3) - topo%min_physs(3)
#elif  __KIND == __DOUBLE_PRECISION
            len_phys(1) = topo%max_physd(1) - topo%min_physd(1)
            len_phys(2) = topo%max_physd(2) - topo%min_physd(2)
            len_phys(3) = topo%max_physd(3) - topo%min_physd(3)
#endif

         ! For all particles ip on this proc
         iproc = ppm_rank+1
            DO ip = 1, allnp(iproc)
               
               ! For all particles on other procs check if we have it in our ghosts particles
               DO iproc2 = 1, ppm_nproc
                     !IF (.NOT. (iproc2 .EQ. iproc)) THEN

                        DO ip2 = 1, allnp(iproc2)

                           ! Check if they enclose each other
                           IF ( (( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip,iproc)) ) .AND. &
               &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip2,iproc2)) ) .AND. &
               &                  ( (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip,iproc)) .AND. &
               &                    (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) )) &
               &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. & 
               &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip,iproc)) .AND. & 
               &                       (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip,iproc)))) ) THEN

                              
                                 ! We need to find ip2 in the ghost of iproc
                                 DO in = 1, mp
                                    
                                    IF(      (ABS(allxp(1,ip2,iproc2)-xp(1,in)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)) .LT. lmyeps &
               &                        .AND. ABS(allxp(3,ip2,iproc2)-xp(3,in)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
                                       print *, ' error ', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
         &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO

                           ENDIF

                        !if we have periodicity, also check that
                        IF (topo%bcdef(1).EQ.ppm_param_bcdef_periodic) THEN

                              IF ( (( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip,iproc)) ) .AND. &
               &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
               &                         allghost_req(1,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip2,iproc2)) ) .AND. &
               &                  ( (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip,iproc)) .AND. &
               &                    (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) )) &
               &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. & 
               &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip,iproc)) .AND. & 
               &                       (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip,iproc))))) THEN
 
                                 ! We need to find ip2 in the ghost of iproc
                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)+len_phys(1)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)) .LT. lmyeps &
               &                        .AND. ABS(allxp(3,ip2,iproc2)-xp(3,in)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
                                       print *, ' error x', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
         &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO


                              ENDIF

                              IF ( (( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip,iproc)) ) .AND. &
               &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
               &                         allghost_req(1,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip2,iproc2)) ) .AND. &
               &                  ( (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip,iproc)) .AND. &
               &                    (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) )) &
               &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. & 
               &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip,iproc)) .AND. & 
               &                       (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip,iproc))))) THEN

                              
                                 ! We need to find ip2 in the ghost of iproc
                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)-len_phys(1)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)) .LT. lmyeps &
               &                        .AND. ABS(allxp(3,ip2,iproc2)-xp(3,in)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
                                       print *, ' error x', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
         &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO

                              ENDIF



                           ENDIF

                        IF (topo%bcdef(3).EQ.ppm_param_bcdef_periodic) THEN

                              IF ( (( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc)) ) .AND. &
               &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
               &                         allghost_req(2,ip2,iproc2)) ) .AND. &
               &                  ( (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip,iproc)) .AND. &
               &                    (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) )) &
               &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. & 
               &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc)) .AND. & 
               &                       (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip,iproc))))) THEN
                                 ! We need to find ip2 in the ghost of iproc
                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)+len_phys(2)) .LT. lmyeps &
               &                        .AND. ABS(allxp(3,ip2,iproc2)-xp(3,in)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
                                       print *, ' error y', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
         &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO

                              ENDIF  
                              IF ( (( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc)) ) .AND. &
               &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
               &                         allghost_req(2,ip2,iproc2)) ) .AND. &
               &                  ( (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip,iproc)) .AND. &
               &                    (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) )) &
               &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. & 
               &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc)) .AND. & 
               &                       (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip,iproc)))) ) THEN
                                 ! We need to find ip2 in the ghost of iproc
                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)-len_phys(2)) .LT. lmyeps&
               &                        .AND. ABS(allxp(3,ip2,iproc2)-xp(3,in)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
                                       print *, ' error y', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
         &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO

                              ENDIF

                           ENDIF

!
                        IF (topo%bcdef(5).EQ.ppm_param_bcdef_periodic) THEN

                              IF ( (( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip,iproc)) ) .AND. &
               &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip2,iproc2)) ) .AND. &
               &                  ( (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)-len_phys(3)) .LT. &
               &                         allghost_req(3,ip,iproc)) .AND. &
               &                    (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)-len_phys(3)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) )) &
               &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. & 
               &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip,iproc)) .AND. & 
               &                       (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)-len_phys(3)) .LT. &
               &                         allghost_req(3,ip,iproc))))) THEN
                                 ! We need to find ip2 in the ghost of iproc
                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)) .LT. lmyeps &
               &                        .AND. ABS(allxp(3,ip2,iproc2)-xp(3,in)+len_phys(3)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
                                       print *, ' error z', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
         &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO

                              ENDIF  
                              IF ( (( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip,iproc)) ) .AND. &
               &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip2,iproc2)) ) .AND. &
               &                  ( (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)+len_phys(3)) .LT. &
               &                         allghost_req(3,ip,iproc)) .AND. &
               &                    (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)+len_phys(3)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) )) &
               &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. & 
               &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip,iproc)) .AND. & 
               &                       (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)+len_phys(3)) .LT. &
               &                         allghost_req(3,ip,iproc)))) ) THEN
                                 ! We need to find ip2 in the ghost of iproc
                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)) .LT. lmyeps&
               &                        .AND. ABS(allxp(3,ip2,iproc2)-xp(3,in)-len_phys(3)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
                                       print *, ' error z', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
         &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO

                              ENDIF

                           ENDIF


                        IF (topo%bcdef(1).EQ.ppm_param_bcdef_periodic .AND. topo%bcdef(3).EQ.ppm_param_bcdef_periodic &
                        & .AND. topo%bcdef(5).EQ.ppm_param_bcdef_periodic) THEN
                              
                             IF  ((( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc)) ) .AND. &
               &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
               &                         allghost_req(1,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
               &                         allghost_req(2,ip2,iproc2)) )  .AND. &
               &                  ( (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)-len_phys(3)) .LT. &
               &                         allghost_req(3,ip,iproc)) .AND. &
               &                    (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)-len_phys(3)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) )) &
               &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. & 
               &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc)) .AND. & 
               &                       (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)-len_phys(3)) .LT. &
               &                         allghost_req(3,ip,iproc))))) THEN
                                 ! We need to find ip2 in the ghost of iproc

                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)+len_phys(1)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)+len_phys(2)) .LT. lmyeps&
               &                        .AND. ABS(allxp(3,ip2,iproc2)-xp(3,in)+len_phys(3)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
                                       print *, ' error xyz', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
         &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO
                              ENDIF

                           IF  ((( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc)) ) .AND. &
               &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
               &                         allghost_req(1,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
               &                         allghost_req(2,ip2,iproc2)) ) .AND. &
               &                  ( (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)+len_phys(3)) .LT. &
               &                         allghost_req(3,ip,iproc)) .AND. &
               &                    (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)+len_phys(3)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) )) &
               &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. & 
               &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc)) .AND. & 
               &                       (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)+len_phys(3)) .LT. &
               &                         allghost_req(3,ip,iproc))))) THEN
                              
                                 ! We need to find ip2 in the ghost of iproc
                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)-len_phys(1)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)-len_phys(2)) .LT. lmyeps&
               &                        .AND. ABS(allxp(3,ip2,iproc2)-xp(3,in)-len_phys(3)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
                                       print *, ' error xyz', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
         &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO

                              ENDIF


                              IF  ((( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc)) ) .AND. &
               &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
               &                         allghost_req(1,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
               &                         allghost_req(2,ip2,iproc2)) ) .AND. &
               &                  ( (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip,iproc)) .AND. &
               &                    (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) )) &
               &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. & 
               &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc)) .AND. & 
               &                       (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip,iproc))))) THEN
                                 ! We need to find ip2 in the ghost of iproc

                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)-len_phys(1)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)-len_phys(2)) .LT. lmyeps&
               &                        .AND. ABS(allxp(3,ip2,iproc2)-xp(3,in)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
                                       print *, ' error xy', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
         &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO
                              ENDIF

                            IF  ((( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc)) ) .AND. &
               &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
               &                         allghost_req(1,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
               &                         allghost_req(2,ip2,iproc2)) ) .AND. &
               &                  ( (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) )) &
               &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. & 
               &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc)) .AND. & 
               &                       (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip,iproc))))) THEN

                              
                                 ! We need to find ip2 in the ghost of iproc
                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)+len_phys(1)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)+len_phys(2)) .LT. lmyeps&
               &                        .AND. ABS(allxp(3,ip2,iproc2)-xp(3,in)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
                                       print *, ' error xy', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
         &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO

                              ENDIF


                            IF  ((( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc)) ) .AND. &
               &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
               &                         allghost_req(1,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
               &                         allghost_req(2,ip2,iproc2)) ) .AND. &
               &                  ( (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip,iproc)) .AND. &
               &                    (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) )) &
               &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. & 
               &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc)) .AND. & 
               &                       (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip,iproc))))) THEN
                                 ! We need to find ip2 in the ghost of iproc

                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)+len_phys(1)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)-len_phys(2)) .LT. lmyeps&
               &                        .AND. ABS(allxp(3,ip2,iproc2)-xp(3,in)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
                                       print *, ' error xy', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
         &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO
                              ENDIF

                           IF  ((( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc)) ) .AND. &
               &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
               &                         allghost_req(1,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
               &                         allghost_req(2,ip2,iproc2)) ) .AND. &
               &                  ( (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) )) &
               &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. & 
               &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc)) .AND. & 
               &                       (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)) .LT. &
               &                         allghost_req(3,ip,iproc))))) THEN

                              
                                 ! We need to find ip2 in the ghost of iproc
                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)-len_phys(1)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)+len_phys(2)) .LT. lmyeps&
               &                        .AND. ABS(allxp(3,ip2,iproc2)-xp(3,in)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
                              print *, ' error xy', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
                                       &      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO

                              ENDIF

                            IF  ((( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip,iproc)) ) .AND. &
               &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
               &                         allghost_req(1,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip2,iproc2)) ) .AND. &
               &                  ( (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)+len_phys(3)) .LT. &
               &                         allghost_req(3,ip,iproc)) .AND. &
               &                    (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)+len_phys(3)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) )) &
               &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. & 
               &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip,iproc)) .AND. & 
               &                       (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)+len_phys(3)) .LT. &
               &                         allghost_req(3,ip,iproc))))) THEN
                                 ! We need to find ip2 in the ghost of iproc

                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)-len_phys(1)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)) .LT. lmyeps&
               &                        .AND. ABS(allxp(3,ip2,iproc2)-xp(3,in)-len_phys(3)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
                                       print *, ' error xz', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
         &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO
                              ENDIF

                           IF  ((( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip,iproc)) ) .AND. &
               &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
               &                         allghost_req(1,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip2,iproc2)) ) .AND. &
               &                  ( (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)-len_phys(3)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)-len_phys(3)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) )) &
               &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. & 
               &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip,iproc)) .AND. & 
               &                       (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)-len_phys(3)) .LT. &
               &                         allghost_req(3,ip,iproc))))) THEN

                              
                                 ! We need to find ip2 in the ghost of iproc
                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)+len_phys(1)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)) .LT. lmyeps&
               &                        .AND. ABS(allxp(3,ip2,iproc2)-xp(3,in)+len_phys(3)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
                                       print *, ' error xz', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
         &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO

                              ENDIF

                            IF  ((( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip,iproc)) ) .AND. &
               &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
               &                         allghost_req(1,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip2,iproc2)) ) .AND. &
               &                  ( (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)+len_phys(3)) .LT. &
               &                         allghost_req(3,ip,iproc)) .AND. &
               &                    (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)+len_phys(3)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) )) &
               &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)-len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. & 
               &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip,iproc)) .AND. & 
               &                       (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)+len_phys(3)) .LT. &
               &                         allghost_req(3,ip,iproc))))) THEN
                                 ! We need to find ip2 in the ghost of iproc

                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)+len_phys(1)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)) .LT. lmyeps&
               &                        .AND. ABS(allxp(3,ip2,iproc2)-xp(3,in)-len_phys(3)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
                                       print *, ' error xz', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
         &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO
                              ENDIF

                           IF  ((( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip,iproc)) ) .AND. &
               &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
               &                         allghost_req(1,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip2,iproc2)) ) .AND. &
               &                  ( (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)-len_phys(3)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)-len_phys(3)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) )) &
               &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)+len_phys(1)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. & 
               &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)) .LT. &
               &                         allghost_req(2,ip,iproc)) .AND. & 
               &                       (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)-len_phys(3)) .LT. &
               &                         allghost_req(3,ip,iproc))))) THEN

                              
                                 ! We need to find ip2 in the ghost of iproc
                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)-len_phys(1)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)) .LT. lmyeps&
               &                        .AND. ABS(allxp(3,ip2,iproc2)-xp(3,in)+len_phys(3)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
                                       print *, ' error xz', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
         &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO

                              ENDIF
                              

                              IF  ((( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc)) ) .AND. &
               &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
               &                         allghost_req(2,ip2,iproc2)) ) .AND. &
               &                  ( (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)+len_phys(3)) .LT. &
               &                         allghost_req(3,ip,iproc)) .AND. &
               &                    (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)+len_phys(3)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) )) &
               &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. & 
               &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc)) .AND. & 
               &                       (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)+len_phys(3)) .LT. &
               &                         allghost_req(3,ip,iproc))))) THEN
                                 ! We need to find ip2 in the ghost of iproc

                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)+len_phys(2)) .LT. lmyeps&
               &                        .AND. ABS(allxp(3,ip2,iproc2)-xp(3,in)-len_phys(3)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
                                       print *, ' error yz', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
         &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO
                              ENDIF

                           IF  ((( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc)) ) .AND. &
               &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
               &                         allghost_req(2,ip2,iproc2)) ) .AND. &
               &                  ( (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)-len_phys(3)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)-len_phys(3)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) )) &
               &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. & 
               &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc)) .AND. & 
               &                       (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)-len_phys(3)) .LT. &
               &                         allghost_req(3,ip,iproc))))) THEN

                              
                                 ! We need to find ip2 in the ghost of iproc
                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)-len_phys(2)) .LT. lmyeps&
               &                        .AND. ABS(allxp(3,ip2,iproc2)-xp(3,in)+len_phys(3)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
!                                        print *, ' errmak(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO

                              ENDIF

                           IF  ((( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc)) ) .AND. &
               &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
               &                         allghost_req(2,ip2,iproc2)) ) .AND. &
               &                  ( (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)-len_phys(3)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)-len_phys(3)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) )) &
               &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. & 
               &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)-len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc)) .AND. & 
               &                       (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)-len_phys(3)) .LT. &
               &                         allghost_req(3,ip,iproc))))) THEN

                              
                                 ! We need to find ip2 in the ghost of iproc
                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)+len_phys(2)) .LT. lmyeps&
               &                        .AND. ABS(allxp(3,ip2,iproc2)-xp(3,in)+len_phys(3)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
                                       print *, ' error yz', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
         &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO

                              ENDIF
                              
                              IF  ((( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc)) ) .AND. &
               &                  ( (ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip2,iproc2)) .AND. &
               &                    (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
               &                         allghost_req(2,ip2,iproc2)) ) .AND. &
               &                  ( (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)+len_phys(3)) .LT. &
               &                         allghost_req(3,ip,iproc)) .AND. &
               &                    (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)+len_phys(3)) .LT. &
               &                         allghost_req(3,ip2,iproc2)) )) &
               &                  .OR. (has_one_way .AND.((ABS(allxp(1,ip,iproc)-allxp(1,ip2,iproc2)) .LT. &
               &                         allghost_req(1,ip,iproc)) .AND. & 
               &                       (ABS(allxp(2,ip,iproc)-allxp(2,ip2,iproc2)+len_phys(2)) .LT. &
               &                         allghost_req(2,ip,iproc)) .AND. & 
               &                       (ABS(allxp(3,ip,iproc)-allxp(3,ip2,iproc2)+len_phys(3)) .LT. &
               &                         allghost_req(3,ip,iproc))))) THEN
                                 ! We need to find ip2 in the ghost of iproc

                                 DO in = Npart+1, mp
                                    
                                    IF((ABS(allxp(1,ip2,iproc2)-xp(1,in)) .LT. lmyeps &
               &                        .AND. ABS(allxp(2,ip2,iproc2)-xp(2,in)-len_phys(2)) .LT. lmyeps&
               &                        .AND. ABS(allxp(3,ip2,iproc2)-xp(3,in)-len_phys(3)) .LT. lmyeps)) THEN
                                       !we found it
                                       !print *, 'found'
!                                        IF (ppm_rank .eq. 0) THEN
!                                           print *, ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), xp(1,in), xp(2,in)
!                                        ENDIF
                                       have_xp(in) = .TRUE.
                                       EXIT
                                    ENDIF
                                    IF (in .EQ. mp) THEN
                                       !We have not found any until the last one -> error
                                       print *, ' error yz', ppm_rank, allxp(1,ip2,iproc2), allxp(2,ip2,iproc2), &
         &                                      allxp(1,ip,iproc),allxp(2,ip,iproc)
                                        map_ok = .FALSE.
                                       GOTO 9999
                                    ENDIF
                                 ENDDO
                              ENDIF


                           ENDIF
                           
                     ENDDO

                  !ENDIF

               ENDDO

            ENDDO

         ! plot the number of particles we do not need, but have in the ghosts
         IF(ppm_debug .GT. 0) THEN
            j = 0
            DO i = 1,mp-Npart
               IF(have_xp(i) .EQV. .FALSE.) THEN
                  j = j+1
               ENDIF
            ENDDO

            IF (mp-Npart .GT. 0) THEN
               print *,'[', ppm_rank,']', int(100*float(j)/float((mp-Npart))) , '% of ghost particles not needed!'
            ELSE
               print *,'[', ppm_rank,']', 'has no ghost particles and does not need any!'
            ENDIF
         ENDIF

      ENDIF

      CALL ppm_alloc(allxp,lda,ppm_param_dealloc,info)
      CALL ppm_alloc(allghost_req,lda,ppm_param_dealloc,info)
      CALL ppm_alloc(have_xp,lda,ppm_param_dealloc,info)

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
   
      !MPI allreduce for success

      CALL substop('ppm_map_check_ghosts',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
         IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_map_check_ghosts',       &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (Npart .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_check_ghosts',         &
     &            'Npart must be >= 0',__LINE__,info)
              GOTO 8888
          ENDIF
          CALL ppm_check_topoid(topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_check_ghosts', &
     &             'topoid invalid',&
     &              __LINE__, info)
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_check_ghosts_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_check_ghosts_d
#endif
