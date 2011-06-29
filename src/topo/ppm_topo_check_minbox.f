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
      SUBROUTINE ppm_topo_check_minbox_s(topoid,xp,ghost_req,Npart,has_one_way,topo_ok,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_topo_check_minbox_d(topoid,xp,ghost_req,Npart,has_one_way,topo_ok,info)
#endif
      !!! Checks if all boxes fulfill minimum sizes for all particle inside this box
      !!! Assumes neighbors are correct
      !!! (on the local processor)
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
      USE ppm_module_find_neigh

      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
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
      LOGICAL                 , INTENT(  IN) ::  has_one_way
      !!! Is one way interaction between particles
      LOGICAL                 , INTENT(  OUT) :: topo_ok
      !!! Is the topology consistent
      INTEGER                 , INTENT(IN   ) :: topoid
      !!! The topology to be checked
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK)                       :: t0, lmyeps
      INTEGER                        :: ipart,idom,j,ison, in, iid
      LOGICAL                        :: valid
      INTEGER, DIMENSION(:), POINTER :: bcdef => NULL()
      TYPE(ppm_t_topo)     , POINTER :: topo  => NULL()

      ! Needed for has_one_way neighbor checking
      INTEGER , DIMENSION(  :), POINTER :: nneigh  => NULL()
      INTEGER , DIMENSION(:,:), POINTER :: ineigh  => NULL()
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
      CALL substart('ppm_topo_check_minbox',t0,info)
      topo_ok = .TRUE.

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

      !-------------------------------------------------------------------------
      !  Check all particles
      !-------------------------------------------------------------------------
      IF (ppm_dim .EQ. 2) THEN
          
            DO j=1,topo%nsublist
               idom = topo%isublist(j)
               DO ipart=1,Npart
                  ison = 0
#if    __KIND == __SINGLE_PRECISION
                   IF (xp(1,ipart).GE.topo%min_subs(1,idom).AND.  &
     &                 xp(1,ipart).LE.topo%max_subs(1,idom).AND.  &
     &                 xp(2,ipart).GE.topo%min_subs(2,idom).AND.  &
     &                 xp(2,ipart).LE.topo%max_subs(2,idom)) THEN
                     !------------------------------------------------------------
                     !  In the non-periodic case, allow particles that are
                     !  exactly ON an upper EXTERNAL boundary.
                     !------------------------------------------------------------
                     IF((xp(1,ipart).LT.topo%max_subs(1,idom) .OR.  &
      &                (topo%subs_bc(2,idom).EQ.1           .AND.  &
      &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND.  &
      &                (xp(2,ipart).LT.topo%max_subs(2,idom) .OR.  &
      &                (topo%subs_bc(4,idom).EQ.1           .AND.  &
      &                bcdef(4).NE. ppm_param_bcdef_periodic))) THEN
                        !------------------------------------------------------------
                        !  Box has to be at least as large as the ghost_req
                        !------------------------------------------------------------
                           IF((topo%max_subs(1,idom)-topo%min_subs(1,idom)) .LT. ghost_req(1,ipart)-lmyeps .OR. &
         &                    (topo%max_subs(2,idom)-topo%min_subs(2,idom)) .LT. ghost_req(2,ipart)-lmyeps) THEN

                   
#elif  __KIND == __DOUBLE_PRECISION
                   IF (xp(1,ipart).GE.topo%min_subd(1,idom).AND.  &
     &                 xp(1,ipart).LE.topo%max_subd(1,idom).AND.  &
     &                 xp(2,ipart).GE.topo%min_subd(2,idom).AND.  &
     &                 xp(2,ipart).LE.topo%max_subd(2,idom)) THEN
                     !------------------------------------------------------------
                     !  In the non-periodic case, allow particles that are
                     !  exactly ON an upper EXTERNAL boundary.
                     !------------------------------------------------------------
                     IF((xp(1,ipart).LT.topo%max_subd(1,idom) .OR.  &
      &                (topo%subs_bc(2,idom).EQ.1           .AND.  &
      &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND.  &
      &                (xp(2,ipart).LT.topo%max_subd(2,idom) .OR.  &
      &                (topo%subs_bc(4,idom).EQ.1           .AND.  &
      &                bcdef(4).NE. ppm_param_bcdef_periodic))) THEN
                        !------------------------------------------------------------
                        !  Box has to be at least as large as the ghost_req
                        !------------------------------------------------------------
                        IF((topo%max_subd(1,idom)-topo%min_subd(1,idom)) .LT. ghost_req(1,ipart)-lmyeps .OR. &
      &                    (topo%max_subd(2,idom)-topo%min_subd(2,idom)) .LT. ghost_req(2,ipart)-lmyeps) THEN
                           !print *, topo%min_subd(1,idom), topo%min_subd(2,idom), topo%max_subd(1,idom), topo%max_subd(2,idom)
                           print *, topo%min_subd(1,idom)+ghost_req(1,ipart)-topo%max_subd(1,idom), lmyeps
                           print *, topo%min_subd(2,idom)+ghost_req(2,ipart)-topo%max_subd(2,idom)
#endif
                           ! particle is in one of my subs
                           !print *, xp(1,ipart), xp(2,ipart), ghost_req(1,ipart), ghost_req(2,ipart)
                           !print *, xp(1,ipart)-ghost_req(1,ipart), xp(1,ipart)+ghost_req(1,ipart)
                           !print *, xp(2,ipart)-ghost_req(2,ipart), xp(2,ipart)+ghost_req(2,ipart)
                           ison = 1
                           EXIT
                        ENDIF
                     ENDIF
                   ENDIF
              ENDDO
              IF (ison .EQ. 1) THEN
                  ! found a particle that violates the topology
                  topo_ok = .FALSE.
                  ! no need to search any further
                  GOTO 9999
              ENDIF
          ENDDO
      ELSEIF (ppm_dim .EQ. 3) THEN
          
            DO j=1,topo%nsublist
               idom = topo%isublist(j)
               DO ipart=1,Npart
                  ison = 0
#if    __KIND == __SINGLE_PRECISION
                   IF (xp(1,ipart).GE.topo%min_subs(1,idom).AND.  &
     &                 xp(1,ipart).LE.topo%max_subs(1,idom).AND.  &
     &                 xp(2,ipart).GE.topo%min_subs(2,idom).AND.  &
     &                 xp(2,ipart).LE.topo%max_subs(2,idom).AND.  &
     &                 xp(3,ipart).GE.topo%min_subs(3,idom).AND.  &
     &                 xp(3,ipart).LE.topo%max_subs(3,idom)) THEN
                        !------------------------------------------------------------
                        !  In the non-periodic case, allow particles that are
                        !  exactly ON an upper EXTERNAL boundary.
                        !------------------------------------------------------------
                        IF((xp(1,ipart).LT.topo%max_subs(1,idom) .OR.  &
         &                (topo%subs_bc(2,idom).EQ.1           .AND.  &
         &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND.  &
         &                (xp(2,ipart).LT.topo%max_subs(2,idom) .OR.  &
         &                (topo%subs_bc(4,idom).EQ.1           .AND.  &
         &                bcdef(4).NE. ppm_param_bcdef_periodic))    .AND.  &
         &                (xp(3,ipart).LT.topo%max_subs(3,idom) .OR.  &
         &                (topo%subs_bc(6,idom).EQ.1           .AND.  &
         &                bcdef(6).NE. ppm_param_bcdef_periodic))   ) THEN
                        !------------------------------------------------------------
                        !  Box has to be at least as large as the ghost_req
                        !------------------------------------------------------------
                        IF((topo%max_subs(1,idom)-topo%min_subs(1,idom)) .LT. ghost_req(1,ipart)-lmyeps .OR. &
         &                    (topo%max_subs(2,idom)-topo%min_subs(2,idom)) .LT. ghost_req(2,ipart)-lmyeps .OR. &
      &                         (topo%max_subs(3,idom)-topo%min_subs(3,idom)) .LT. ghost_req(3,ipart)-lmyeps) THEN
#elif  __KIND == __DOUBLE_PRECISION
                   IF (xp(1,ipart).GE.topo%min_subd(1,idom).AND.  &
     &                 xp(1,ipart).LE.topo%max_subd(1,idom).AND.  &
     &                 xp(2,ipart).GE.topo%min_subd(2,idom).AND.  &
     &                 xp(2,ipart).LE.topo%max_subd(2,idom).AND.  &
     &                 xp(3,ipart).GE.topo%min_subd(3,idom).AND.  &
     &                 xp(3,ipart).LE.topo%max_subd(3,idom)) THEN
                        !------------------------------------------------------------
                        !  In the non-periodic case, allow particles that are
                        !  exactly ON an upper EXTERNAL boundary.
                        !------------------------------------------------------------
                        IF((xp(1,ipart).LT.topo%max_subd(1,idom) .OR.  &
         &                (topo%subs_bc(2,idom).EQ.1           .AND.  &
         &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND.  &
         &                (xp(2,ipart).LT.topo%max_subd(2,idom) .OR.  &
         &                (topo%subs_bc(4,idom).EQ.1           .AND.  &
         &                bcdef(4).NE. ppm_param_bcdef_periodic))    .AND.  &
         &                (xp(3,ipart).LT.topo%max_subd(3,idom) .OR.  &
         &                (topo%subs_bc(6,idom).EQ.1           .AND.  &
         &                bcdef(6).NE. ppm_param_bcdef_periodic))   ) THEN
                        !------------------------------------------------------------
                        !  Box has to be at least as large as the ghost_req
                        !------------------------------------------------------------
                        IF((topo%max_subd(1,idom)-topo%min_subd(1,idom)) .LT. ghost_req(1,ipart)-lmyeps .OR. &
         &                    (topo%max_subd(2,idom)-topo%min_subd(2,idom)) .LT. ghost_req(2,ipart)-lmyeps .OR. &
      &                         (topo%max_subd(3,idom)-topo%min_subd(3,idom)) .LT. ghost_req(3,ipart)-lmyeps) THEN
#endif
                           ! particle is in one of my subs
                           ison = 1
                           EXIT
                        ENDIF
                      ENDIF
                   ENDIF
              ENDDO
              IF (ison .EQ. 1) THEN
                  ! found a particle that violates the topology
                  topo_ok = .FALSE.
                  ! no need to search any further
                  GOTO 9999
              ENDIF
          ENDDO
      END IF

      

      IF (has_one_way) THEN
         ! If has one way, then we check if neighbors are at least the size of particle size

         ! Get neighbors
#if    __KIND == __SINGLE_PRECISION
         CALL ppm_find_neigh(topo%min_physs,topo%max_physs,topo%bcdef, &
     &                    topo%min_subs,topo%max_subs,topo%nsubs,nneigh,ineigh,info)
#elif  __KIND == __DOUBLE_PRECISION
         CALL ppm_find_neigh(topo%min_physd,topo%max_physd,topo%bcdef, &
     &                    topo%min_subd,topo%max_subd,topo%nsubs,nneigh,ineigh,info)
#endif

         IF (ppm_dim .EQ. 2) THEN
            
               DO j=1,topo%nsublist
                  idom = topo%isublist(j)
                  DO ipart=1,Npart
                     ison = 0
#if    __KIND == __SINGLE_PRECISION
                     IF (xp(1,ipart).GE.topo%min_subs(1,idom).AND.  &
      &                 xp(1,ipart).LE.topo%max_subs(1,idom).AND.  &
      &                 xp(2,ipart).GE.topo%min_subs(2,idom).AND.  &
      &                 xp(2,ipart).LE.topo%max_subs(2,idom)) THEN
                        !------------------------------------------------------------
                        !  In the non-periodic case, allow particles that are
                        !  exactly ON an upper EXTERNAL boundary.
                        !------------------------------------------------------------
                        IF((xp(1,ipart).LT.topo%max_subs(1,idom) .OR.  &
         &                (topo%subs_bc(2,idom).EQ.1           .AND.  &
         &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND.  &
         &                (xp(2,ipart).LT.topo%max_subs(2,idom) .OR.  &
         &                (topo%subs_bc(4,idom).EQ.1           .AND.  &
         &                bcdef(4).NE. ppm_param_bcdef_periodic))) THEN
                           !------------------------------------------------------------
                           !  Neighbor box has to be at least as large as the ghost_req
                           !------------------------------------------------------------
                           ! Iterate through all neighbors
                              DO in = 1, nneigh(idom)
                                 iid = ineigh(in,idom)
                                 IF((topo%max_subs(1,iid)-topo%min_subs(1,iid)) .LT. ghost_req(1,ipart)-lmyeps .OR. &
               &                    (topo%max_subs(2,iid)-topo%min_subs(2,iid)) .LT. ghost_req(2,ipart)-lmyeps) THEN
                                    ison = 1
                                    EXIT
                                 ENDIF
                              ENDDO
                        ENDIF
                     ENDIF

#elif  __KIND == __DOUBLE_PRECISION
                     IF (xp(1,ipart).GE.topo%min_subd(1,idom).AND.  &
      &                 xp(1,ipart).LE.topo%max_subd(1,idom).AND.  &
      &                 xp(2,ipart).GE.topo%min_subd(2,idom).AND.  &
      &                 xp(2,ipart).LE.topo%max_subd(2,idom)) THEN
                        !------------------------------------------------------------
                        !  In the non-periodic case, allow particles that are
                        !  exactly ON an upper EXTERNAL boundary.
                        !------------------------------------------------------------
                        IF((xp(1,ipart).LT.topo%max_subd(1,idom) .OR.  &
         &                (topo%subs_bc(2,idom).EQ.1           .AND.  &
         &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND.  &
         &                (xp(2,ipart).LT.topo%max_subd(2,idom) .OR.  &
         &                (topo%subs_bc(4,idom).EQ.1           .AND.  &
         &                bcdef(4).NE. ppm_param_bcdef_periodic))) THEN
                           !------------------------------------------------------------
                           !  Neighbor box has to be at least as large as the ghost_req
                           !------------------------------------------------------------
                           ! Iterate through all neighbors
                           DO in = 1, nneigh(idom)
                              iid = ineigh(in,idom)
                              IF((topo%max_subd(1,iid)-topo%min_subd(1,iid)) .LT. ghost_req(1,ipart)-lmyeps .OR. &
            &                    (topo%max_subd(2,iid)-topo%min_subd(2,iid)) .LT. ghost_req(2,ipart)-lmyeps) THEN
                                 ison = 1
                                 EXIT
                              ENDIF
                           ENDDO
                        ENDIF
                     ENDIF
#endif
               ENDDO
               IF (ison .EQ. 1) THEN
                     ! found a particle that violates the topology
                     topo_ok = .FALSE.
                     ! no need to search any further
                     GOTO 9999
               ENDIF
            ENDDO
         ELSEIF (ppm_dim .EQ. 3) THEN
            DO j=1,topo%nsublist
                  idom = topo%isublist(j)
                  DO ipart=1,Npart
                     ison = 0
#if    __KIND == __SINGLE_PRECISION
                     IF (xp(1,ipart).GE.topo%min_subs(1,idom).AND.  &
      &                 xp(1,ipart).LE.topo%max_subs(1,idom).AND.  &
      &                 xp(2,ipart).GE.topo%min_subs(2,idom).AND.  &
      &                 xp(2,ipart).LE.topo%max_subs(2,idom).AND.  &
      &                 xp(3,ipart).GE.topo%min_subs(3,idom).AND.  &
      &                 xp(3,ipart).LE.topo%max_subs(3,idom)) THEN
                     !------------------------------------------------------------
                     !  In the non-periodic case, allow particles that are
                     !  exactly ON an upper EXTERNAL boundary.
                     !------------------------------------------------------------
                     IF((xp(1,ipart).LT.topo%max_subs(1,idom) .OR.  &
      &                (topo%subs_bc(2,idom).EQ.1           .AND.  &
      &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND.  &
      &                (xp(2,ipart).LT.topo%max_subs(2,idom) .OR.  &
      &                (topo%subs_bc(4,idom).EQ.1           .AND.  &
      &                bcdef(4).NE. ppm_param_bcdef_periodic))    .AND.  &
      &                (xp(3,ipart).LT.topo%max_subs(3,idom) .OR.  &
      &                (topo%subs_bc(6,idom).EQ.1           .AND.  &
      &                bcdef(6).NE. ppm_param_bcdef_periodic))   ) THEN
                        !------------------------------------------------------------
                        !  Neighbor box has to be at least as large as the ghost_req
                        !------------------------------------------------------------
                        ! Iterate through all neighbors
                           DO in = 1, nneigh(idom)
                              iid = ineigh(in,idom)
                              IF((topo%max_subs(1,iid)-topo%min_subs(1,iid)) .LT. ghost_req(1,ipart)-lmyeps .OR. &
            &                    (topo%max_subs(2,iid)-topo%min_subs(2,iid)) .LT. ghost_req(2,ipart)-lmyeps.OR. &
            &                    (topo%max_subs(3,iid)-topo%min_subs(3,iid)) .LT. ghost_req(3,ipart)-lmyeps) THEN
                                 ison = 1
                                 EXIT
                              ENDIF
                           ENDDO
                        ENDIF
                     ENDIF

#elif  __KIND == __DOUBLE_PRECISION
                     IF (xp(1,ipart).GE.topo%min_subd(1,idom).AND.  &
      &                 xp(1,ipart).LE.topo%max_subd(1,idom).AND.  &
      &                 xp(2,ipart).GE.topo%min_subd(2,idom).AND.  &
      &                 xp(2,ipart).LE.topo%max_subd(2,idom).AND.  &
      &                 xp(3,ipart).GE.topo%min_subd(3,idom).AND.  &
      &                 xp(3,ipart).LE.topo%max_subd(3,idom)) THEN
                        !------------------------------------------------------------
                        !  In the non-periodic case, allow particles that are
                        !  exactly ON an upper EXTERNAL boundary.
                        !------------------------------------------------------------
                        IF((xp(1,ipart).LT.topo%max_subd(1,idom) .OR.  &
         &                (topo%subs_bc(2,idom).EQ.1           .AND.  &
         &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND.  &
         &                (xp(2,ipart).LT.topo%max_subd(2,idom) .OR.  &
         &                (topo%subs_bc(4,idom).EQ.1           .AND.  &
         &                bcdef(4).NE. ppm_param_bcdef_periodic))    .AND.  &
         &                (xp(3,ipart).LT.topo%max_subd(3,idom) .OR.  &
         &                (topo%subs_bc(6,idom).EQ.1           .AND.  &
         &                bcdef(6).NE. ppm_param_bcdef_periodic))   ) THEN
                           !------------------------------------------------------------
                           !  Neighbor box has to be at least as large as the ghost_req
                           !------------------------------------------------------------
                           ! Iterate through all neighbors
                              DO in = 1, nneigh(idom)
                                 iid = ineigh(in,idom)
                                 IF((topo%max_subd(1,iid)-topo%min_subd(1,iid)) .LT. ghost_req(1,ipart)-lmyeps .OR. &
               &                    (topo%max_subd(2,iid)-topo%min_subd(2,iid)) .LT. ghost_req(2,ipart)-lmyeps .OR. &
               &                    (topo%max_subd(3,iid)-topo%min_subd(3,iid)) .LT. ghost_req(3,ipart)-lmyeps) THEN
                                    ison = 1
                                    EXIT
                                 ENDIF
                              ENDDO
                        ENDIF
                     ENDIF
#endif
               ENDDO
               IF (ison .EQ. 1) THEN
                     ! found a particle that violates the topology
                     topo_ok = .FALSE.
                     ! no need to search any further
                     GOTO 9999
               ENDIF
            ENDDO

         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_topo_check_minbox',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
         IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_topo_check_minbox',       &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (Npart .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_check_minbox',         &
     &            'Npart must be >= 0',__LINE__,info)
              GOTO 8888
          ENDIF
          CALL ppm_check_topoid(topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_check_minbox', &
     &             'topoid invalid',&
     &              __LINE__, info)
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_topo_check_minbox_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_topo_check_minbox_d
#endif
