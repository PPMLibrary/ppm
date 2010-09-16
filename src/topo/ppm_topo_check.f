      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_topo_check
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_topo_check_s(topoid,xp,Npart,topo_ok,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_topo_check_d(topoid,xp,Npart,topo_ok,info)
#endif
      !!! Checks if all particles are on the current topology.
      !!! (on the ocal processor)
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
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
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xp
      !!! Particle locations
      INTEGER                 , INTENT(IN   ) :: Npart
      !!! Number of particles on this processor
      LOGICAL                 , INTENT(  OUT) :: topo_ok
      !!! Is the topology consistent
      INTEGER                 , INTENT(IN   ) :: topoid
      !!! The topology to be checked
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK)                       :: t0
      INTEGER                        :: ipart,idom,j,ison
      LOGICAL                        :: valid
      INTEGER, DIMENSION(:), POINTER :: bcdef
      TYPE(ppm_t_topo)     , POINTER :: topo
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_topo_check',t0,info)
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
          DO ipart=1,Npart
              ison = 0
              DO j=1,topo%nsublist
                  idom = topo%isublist(j)
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
#endif
                      ! particle is in one of my subs
                      ison = 1
                      EXIT
                   ENDIF
                   ENDIF
              ENDDO
              IF (ison .EQ. 0) THEN
                  ! found a particle that violates the topology
                  topo_ok = .FALSE.
                  ! no need to search any further
                  GOTO 9999
              ENDIF
          ENDDO
      ELSEIF (ppm_dim .EQ. 3) THEN
          DO ipart=1,Npart
              ison = 0
              DO j=1,topo%nsublist
                  idom = topo%isublist(j)
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
#endif
                      ! particle is in one of my subs
                      ison = 1
                      EXIT
                   ENDIF
                   ENDIF
              ENDDO
              IF (ison .EQ. 0) THEN
                  ! found a particle that violates the topology
                  topo_ok = .FALSE.
                  ! no need to search any further
                  GOTO 9999
              ENDIF
          ENDDO
      END IF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_topo_check',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
         IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_topo_check',       &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (Npart .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_check',         &
     &            'Npart must be >= 0',__LINE__,info)
              GOTO 8888
          ENDIF
          CALL ppm_check_topoid(topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_check', &
     &             'topoid invalid',&
     &              __LINE__, info)
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_topo_check_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_topo_check_d
#endif
