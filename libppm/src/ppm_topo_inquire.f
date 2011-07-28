      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_topo_inquire
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Returnes the current topology ID (in user
      !                 numbering).
      !
      !  Input        : (INT) what        Choose whether to inquire the
      !                                   particle topology ID or the field
      !                                   topology ID. One of:
      !                                        ppm_param_topo_part
      !                                        ppm_param_topo_field
      !                                   Default is ppm_param_topo_part.
      !
      !  Input/output : 
      !
      !  Output       : (INT) topo_id     user-ID of the current particle
      !                                   or field topology
      !                 (INT) info        0 if everything is OK
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_topo_inquire.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.9  2004/10/01 16:33:39  ivos
      !  cosmetics.
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
      !  Revision 1.5  2004/07/16 11:42:45  ivos
      !  Added possibility to inquire field topology.
      !
      !  Revision 1.4  2004/04/14 15:06:36  ivos
      !  removed info=0 after the substart.
      !
      !  Revision 1.3  2004/01/23 17:24:18  ivos
      !  Now includes ppm_define.h for the cpp defines and no longer sets them
      !  in the Makefile.
      !
      !  Revision 1.2  2004/01/22 13:27:57  ivos
      !  Did (1) update of the header, (2) replaced pwrite with ppm_write or
      !  ppm_error calls, (3) inserted validity check of INTENT(IN) arguments
      !  where needed.
      !
      !  Revision 1.1  2003/12/12 16:22:48  ivos
      !  Initial version. returnes the user-ID of the current topology.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_topo_inquire(what,topo_id,info)
      
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER , INTENT(IN   )        :: what
      INTEGER , INTENT(  OUT)        :: topo_id
      INTEGER , INTENT(  OUT)        :: info

      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)          :: t0

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_topo_inquire',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_topo_inquire',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Read the user ID of the current topology 
      !-------------------------------------------------------------------------
      IF (what .EQ. ppm_param_topo_field) THEN
          topo_id = ppm_user_topoid(ppm_field_topoid)
      ELSE
          topo_id = ppm_user_topoid(ppm_topoid)
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_topo_inquire',t0,info)
      RETURN
      END SUBROUTINE ppm_topo_inquire
