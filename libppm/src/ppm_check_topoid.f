      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_check_topoid
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Checks a given topology ID for its validity. Can
      !                 handle both user and internal numberings.
      !
      !  Input        : typ      (I) Type of topoid to be checked. On of:
      !                                 ppm_param_id_internal
      !                                 ppm_param_id_user
      !                 topoid   (I) Topo ID to be checked.
      !
      !  Input/output : 
      !
      !  Output       : valid    (L) Returns .TRUE. if the given topoid is
      !                              valid and defined, .FALSE. otherwise.
      !                 info     (I) Return code. 0 on success.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_check_topoid.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.4  2004/10/01 16:33:32  ivos
      !  cosmetics.
      !
      !  Revision 1.3  2004/10/01 16:08:56  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.2  2004/08/31 13:30:22  ivos
      !  Also included check of association status of arrays.
      !
      !  Revision 1.1  2004/08/31 12:13:47  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_check_topoid(typ,topoid,valid,info)
      
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_write
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER , INTENT(IN   )        :: typ,topoid
      LOGICAL , INTENT(  OUT)        :: valid
      INTEGER , INTENT(  OUT)        :: info

      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)          :: t0

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_check_topoid',t0,info)
      valid = .TRUE.

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_check_topoid',  &
     &            'Please call ppm_init first!',__LINE__,info)
              valid = .FALSE.
              GOTO 9999
          ENDIF
          IF (typ .NE. ppm_param_id_internal .AND. typ .NE.      &
     &        ppm_param_id_user) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_check_topoid',  &
     &            'Invalid type specified!',__LINE__,info)
              valid = .FALSE.
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Validity check for internal topoid
      !-------------------------------------------------------------------------
      IF (typ .EQ. ppm_param_id_internal) THEN
          IF (ASSOCIATED(ppm_user_topoid)) THEN
              IF ((topoid .LT. LBOUND(ppm_user_topoid,1)) .OR.      &
     &            (topoid .GT. ppm_max_topoid)) THEN
                  valid = .FALSE.
                  IF (ppm_debug .GT. 0) THEN
                      CALL ppm_write(ppm_rank,'ppm_check_topoid',   &
     &                    'Internal topoid is out of bounds!',info)
                  ENDIF
              ENDIF
          ELSE
              valid = .FALSE.
              IF (ppm_debug .GT. 0) THEN
                  CALL ppm_write(ppm_rank,'ppm_check_topoid',   &
     &                'No topologies are defined so far!',info)
              ENDIF
          ENDIF
      
      !-------------------------------------------------------------------------
      !  Validity check for user topoid
      !-------------------------------------------------------------------------
      ELSEIF (typ .EQ. ppm_param_id_user) THEN
          IF (ASSOCIATED(ppm_internal_topoid)) THEN
              IF ((topoid .GT. UBOUND(ppm_internal_topoid,1)) .OR.   &
     &            (topoid .LT. LBOUND(ppm_internal_topoid,1))) THEN
                  valid = .FALSE.
                  IF (ppm_debug .GT. 0) THEN
                      CALL ppm_write(ppm_rank,'ppm_check_topoid',    &
     &                    'User topoid is out of bounds!',info)
                  ENDIF
                  GOTO 9999
              ENDIF
              IF (ppm_internal_topoid(topoid) .EQ. -HUGE(topoid)) THEN
                  valid = .FALSE.
                  IF (ppm_debug .GT. 0) THEN
                      CALL ppm_write(ppm_rank,'ppm_check_topoid',  &
     &                    'Desired topology is not defined!',info)
                  ENDIF
              ENDIF
          ELSE
              valid = .FALSE.
              IF (ppm_debug .GT. 0) THEN
                  CALL ppm_write(ppm_rank,'ppm_check_topoid',   &
     &                'No topologies are defined so far!',info)
              ENDIF
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_check_topoid',t0,info)
      RETURN
      END SUBROUTINE ppm_check_topoid
