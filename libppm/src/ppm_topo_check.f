      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_topo_check
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Checks if all particles are on the current
      !                 topology. (on the local processor)
      !
      !  Input        : xp         (F) particle positions
      !                 Npart      (I) number of particles (local)
      !                 topo_id    (I) OPTIONAL. User topoid specifying
      !                                which topo the particles are to be
      !                                checked against. If not present, the
      !                                current topology is assumed.
      !
      !  Output       : topo_ok    (L) .TRUE. if the present topo is OK.
      !                 info       (I) return status. 
      !
      !  Remarks      : 
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_topo_check.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.12  2006/12/06 13:39:55  pchatela
      !  Bugfix: Moved the bcdef assignment down, after Set topoid,
      !   so that we used an initialized value for topoid... How could this run?
      !
      !  Revision 1.11  2006/08/09 16:12:35  ivos
      !  Npart = 0 is now allowed.
      !
      !  Revision 1.10  2005/12/01 08:27:20  ivos
      !  For non-periodic cases, particles ON the upper boundary of the
      !  computational domain are now properly accounted for.
      !
      !  Revision 1.9  2005/05/23 17:34:36  ivos
      !  Added OPTIONAL argument topo_id to allow checking against
      !  arbitrary topology.
      !
      !  Revision 1.8  2004/10/01 16:09:12  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.7  2004/07/26 07:42:35  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.6  2004/07/16 14:47:22  ivos
      !  Added check for ppm_initialized.
      !
      !  Revision 1.5  2004/06/10 16:20:04  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.4  2004/01/23 17:24:18  ivos
      !  Now includes ppm_define.h for the cpp defines and no longer sets them
      !  in the Makefile.
      !
      !  Revision 1.3  2004/01/22 13:27:57  ivos
      !  Did (1) update of the header, (2) replaced pwrite with ppm_write or
      !  ppm_error calls, (3) inserted validity check of INTENT(IN) arguments
      !  where needed.
      !
      !  Revision 1.2  2004/01/06 14:04:15  ivos
      !  Removed unnecessary (since already done in substart) initialization of
      !  info to 0.
      !
      !  Revision 1.1  2003/12/16 12:24:21  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_topo_check_s(xp,Npart,topo_ok,info,topo_id)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_topo_check_d(xp,Npart,topo_ok,info,topo_id)
#endif

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
      USE ppm_module_check_topoid
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      ! Particle locations
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xp
      ! number of particles on this processor
      INTEGER                 , INTENT(IN   ) :: Npart
      ! Is the topology consitent ?
      LOGICAL                 , INTENT(  OUT) :: topo_ok
      ! return status
      INTEGER                 , INTENT(  OUT) :: info
      ! what topology to check against
      INTEGER , OPTIONAL      , INTENT(IN   ) :: topo_id
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                       :: t0
      INTEGER                        :: ipart,idom,j,ison,topoid
      LOGICAL                        :: valid
      INTEGER, DIMENSION(:), POINTER :: bcdef
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
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_topo_check',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (Npart .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_check',  &
     &            'Npart must be >= 0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (PRESENT(topo_id)) THEN
              CALL ppm_check_topoid(ppm_param_id_user,topo_id,valid,info)
              IF (.NOT. valid) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_topo_check',  &
     &                'Topology ID (topo_id) is invalid!',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Set topoid
      !-------------------------------------------------------------------------
      IF (PRESENT(topo_id)) THEN
          topoid = ppm_internal_topoid(topo_id)
      ELSE
          topoid = ppm_topoid
      ENDIF
      bcdef => ppm_bcdef(:,topoid)

      !-------------------------------------------------------------------------
      !  Check all particles
      !-------------------------------------------------------------------------
      IF (ppm_dim .EQ. 2) THEN
          DO ipart=1,Npart
              ison = 0
              DO j=1,ppm_nsublist(topoid)
                  idom = ppm_isublist(j,topoid)
#if    __KIND == __SINGLE_PRECISION
                   IF (xp(1,ipart).GE.ppm_min_subs(1,idom,topoid).AND.   &
     &                 xp(1,ipart).LE.ppm_max_subs(1,idom,topoid).AND.   &
     &                 xp(2,ipart).GE.ppm_min_subs(2,idom,topoid).AND.   & 
     &                 xp(2,ipart).LE.ppm_max_subs(2,idom,topoid)) THEN
                   !------------------------------------------------------------
                   !  In the non-periodic case, allow particles that are
                   !  exactly ON an upper EXTERNAL boundary.
                   !------------------------------------------------------------
                   IF((xp(1,ipart).LT.ppm_max_subs(1,idom,topoid) .OR.  &
     &                (ppm_subs_bc(2,idom,topoid).EQ.1           .AND.  &
     &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND.  &
     &                (xp(2,ipart).LT.ppm_max_subs(2,idom,topoid) .OR.  &
     &                (ppm_subs_bc(4,idom,topoid).EQ.1           .AND.  &
     &                bcdef(4).NE. ppm_param_bcdef_periodic))) THEN
#elif  __KIND == __DOUBLE_PRECISION
                   IF (xp(1,ipart).GE.ppm_min_subd(1,idom,topoid).AND.   &
     &                 xp(1,ipart).LE.ppm_max_subd(1,idom,topoid).AND.   &
     &                 xp(2,ipart).GE.ppm_min_subd(2,idom,topoid).AND.   & 
     &                 xp(2,ipart).LE.ppm_max_subd(2,idom,topoid)) THEN
                   !------------------------------------------------------------
                   !  In the non-periodic case, allow particles that are
                   !  exactly ON an upper EXTERNAL boundary.
                   !------------------------------------------------------------
                   IF((xp(1,ipart).LT.ppm_max_subd(1,idom,topoid) .OR.  &
     &                (ppm_subs_bc(2,idom,topoid).EQ.1           .AND.  &
     &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND.  &
     &                (xp(2,ipart).LT.ppm_max_subd(2,idom,topoid) .OR.  &
     &                (ppm_subs_bc(4,idom,topoid).EQ.1           .AND.  &
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
              DO j=1,ppm_nsublist(topoid)
                  idom = ppm_isublist(j,topoid)
#if    __KIND == __SINGLE_PRECISION
                   IF (xp(1,ipart).GE.ppm_min_subs(1,idom,topoid).AND.   &
     &                 xp(1,ipart).LE.ppm_max_subs(1,idom,topoid).AND.   &
     &                 xp(2,ipart).GE.ppm_min_subs(2,idom,topoid).AND.   & 
     &                 xp(2,ipart).LE.ppm_max_subs(2,idom,topoid).AND.   &
     &                 xp(3,ipart).GE.ppm_min_subs(3,idom,topoid).AND.   & 
     &                 xp(3,ipart).LE.ppm_max_subs(3,idom,topoid)) THEN
                   !------------------------------------------------------------
                   !  In the non-periodic case, allow particles that are
                   !  exactly ON an upper EXTERNAL boundary.
                   !------------------------------------------------------------
                   IF((xp(1,ipart).LT.ppm_max_subs(1,idom,topoid) .OR. &
     &                (ppm_subs_bc(2,idom,topoid).EQ.1           .AND. &
     &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND. &
     &                (xp(2,ipart).LT.ppm_max_subs(2,idom,topoid) .OR. &
     &                (ppm_subs_bc(4,idom,topoid).EQ.1           .AND. &
     &                bcdef(4).NE. ppm_param_bcdef_periodic))    .AND. &
     &                (xp(3,ipart).LT.ppm_max_subs(3,idom,topoid) .OR. &
     &                (ppm_subs_bc(6,idom,topoid).EQ.1           .AND. &
     &                bcdef(6).NE. ppm_param_bcdef_periodic))   ) THEN
#elif  __KIND == __DOUBLE_PRECISION
                   IF (xp(1,ipart).GE.ppm_min_subd(1,idom,topoid).AND.   &
     &                 xp(1,ipart).LE.ppm_max_subd(1,idom,topoid).AND.   &
     &                 xp(2,ipart).GE.ppm_min_subd(2,idom,topoid).AND.   & 
     &                 xp(2,ipart).LE.ppm_max_subd(2,idom,topoid).AND.   &
     &                 xp(3,ipart).GE.ppm_min_subd(3,idom,topoid).AND.   & 
     &                 xp(3,ipart).LE.ppm_max_subd(3,idom,topoid)) THEN
                   !------------------------------------------------------------
                   !  In the non-periodic case, allow particles that are
                   !  exactly ON an upper EXTERNAL boundary.
                   !------------------------------------------------------------
                   IF((xp(1,ipart).LT.ppm_max_subd(1,idom,topoid) .OR. &
     &                (ppm_subs_bc(2,idom,topoid).EQ.1           .AND. &
     &                bcdef(2).NE. ppm_param_bcdef_periodic))    .AND. &
     &                (xp(2,ipart).LT.ppm_max_subd(2,idom,topoid) .OR. &
     &                (ppm_subs_bc(4,idom,topoid).EQ.1           .AND. &
     &                bcdef(4).NE. ppm_param_bcdef_periodic))    .AND. &
     &                (xp(3,ipart).LT.ppm_max_subd(3,idom,topoid) .OR. &
     &                (ppm_subs_bc(6,idom,topoid).EQ.1           .AND. &
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
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_topo_check_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_topo_check_d
#endif
