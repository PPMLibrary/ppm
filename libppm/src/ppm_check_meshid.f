      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_check_meshid
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Checks a given mesh ID for its validity. Can
      !                 handle both user and internal numberings.
      !
      !  Input        : typ      (I) Type of mesh to be checked. On of:
      !                                 ppm_param_id_internal
      !                                 ppm_param_id_user
      !                 meshid   (I) Mesh ID to be checked.
      !                 topoid   (I) Topology ID on which the mesh is
      !                              defined. Internal numbering.
      !
      !  Input/output : 
      !
      !  Output       : valid    (L) Returns .TRUE. if the given meshid is
      !                              valid and defined, .FALSE. otherwise.
      !                 info     (I) Return code. 0 on success.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_check_meshid.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.4  2004/10/01 16:33:32  ivos
      !  cosmetics.
      !
      !  Revision 1.3  2004/10/01 16:08:56  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.2  2004/08/31 13:30:23  ivos
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

      SUBROUTINE ppm_check_meshid(typ,meshid,topoid,valid,info)
      
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_write
      USE ppm_module_check_topoid
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER , INTENT(IN   )        :: typ,meshid,topoid
      LOGICAL , INTENT(  OUT)        :: valid
      INTEGER , INTENT(  OUT)        :: info

      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)          :: t0
      LOGICAL                        :: topo_ok

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_check_meshid',t0,info)
      valid = .TRUE.

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_check_meshid',  &
     &            'Please call ppm_init first!',__LINE__,info)
              valid = .FALSE.
              GOTO 9999
          ENDIF
          IF (typ .NE. ppm_param_id_internal .AND. typ .NE.      &
     &        ppm_param_id_user) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_check_meshid',  &
     &            'Invalid type specified!',__LINE__,info)
              valid = .FALSE.
              GOTO 9999
          ENDIF
          CALL ppm_check_topoid(ppm_param_id_internal,topoid,topo_ok,info)
          IF (.NOT. topo_ok) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_check_meshid',  &
     &            'Topoid out of range!',__LINE__,info)
              valid = .FALSE.
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Validity check for internal meshid
      !-------------------------------------------------------------------------
      IF (typ .EQ. ppm_param_id_internal) THEN
          IF (ASSOCIATED(ppm_meshid(topoid)%user)) THEN
              IF ((meshid .LT. LBOUND(ppm_meshid(topoid)%user,1)) .OR.     &
     &            (meshid .GT. ppm_max_meshid(topoid))) THEN
                  valid = .FALSE.
                  IF (ppm_debug .GT. 0) THEN
                      CALL ppm_write(ppm_rank,'ppm_check_meshid',   &
     &                    'Internal meshid is out of bounds!',info)
                  ENDIF
              ENDIF
          ELSE
              valid = .FALSE.
              IF (ppm_debug .GT. 0) THEN
                  CALL ppm_write(ppm_rank,'ppm_check_meshid',   &
     &                'No meshes are defined so far!',info)
              ENDIF
          ENDIF
      
      !-------------------------------------------------------------------------
      !  Validity check for user meshid
      !-------------------------------------------------------------------------
      ELSEIF (typ .EQ. ppm_param_id_user) THEN
          IF (ASSOCIATED(ppm_meshid(topoid)%internal)) THEN
              IF ((meshid .GT. UBOUND(ppm_meshid(topoid)%internal,1)) .OR.   &
     &            (meshid .LT. LBOUND(ppm_meshid(topoid)%internal,1))) THEN
                  valid = .FALSE.
                  IF (ppm_debug .GT. 0) THEN
                      CALL ppm_write(ppm_rank,'ppm_check_meshid',    &
     &                    'User meshid is out of bounds!',info)
                  ENDIF
                  GOTO 9999
              ENDIF
              IF (meshid .NE. 0) THEN
                  IF (ppm_meshid(topoid)%internal(meshid) .EQ.      &
     &                -HUGE(topoid)) THEN
                      valid = .FALSE.
                      IF (ppm_debug .GT. 0) THEN
                          CALL ppm_write(ppm_rank,'ppm_check_meshid',  &
     &                        'Desired mesh is not defined!',info)
                      ENDIF
                  ENDIF
              ENDIF
          ELSE
              valid = .FALSE.
              IF (ppm_debug .GT. 0) THEN
                  CALL ppm_write(ppm_rank,'ppm_check_meshid',   &
     &                'No meshes are defined so far!',info)
              ENDIF
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_check_meshid',t0,info)
      RETURN
      END SUBROUTINE ppm_check_meshid
