      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_is_initialized
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Return global initialization status of ppm.
      !
      !  Input        :
      !
      !  Input/output : 
      !
      !  Output       : initd         (L) .TRUE. if ppm is initialized
      !                                   (i.e. ppm_init has been called,
      !                                   but not yet ppm_finalize).
      !                                   Otherwise .FALSE.
      !                 info          (I) error status. 0 on success.
      !
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_is_initialized.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.4  2004/10/01 16:33:34  ivos
      !  cosmetics.
      !
      !  Revision 1.3  2004/10/01 16:09:02  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.2  2004/07/26 07:46:37  ivos
      !  Changed to use single-interface modules. Updated all USE statements.
      !
      !  Revision 1.1  2004/07/16 11:05:48  ivos
      !  Initial implementation. Tested.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_is_initialized(initd,info)
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data, ONLY: ppm_kind_double, ppm_initialized
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      LOGICAL         , INTENT(  OUT) :: initd
      INTEGER         , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)           :: t0
      
      !-------------------------------------------------------------------------
      !  Initialize 
      !-------------------------------------------------------------------------
      CALL substart('ppm_is_initialized',t0,info)

      !-------------------------------------------------------------------------
      !  Read global status
      !-------------------------------------------------------------------------
      initd = ppm_initialized

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_is_initialized',t0,info)
      RETURN
      END SUBROUTINE ppm_is_initialized
