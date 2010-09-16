      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_is_initialized
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_is_initialized(initd,info)
      !!! Returns global initialization status of ppm
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
      !!! `TRUE` if ppm is initialized (i.e. `ppm_init` has been called,
      !!! but not yet ppm_finalize). Otherwise `FALSE`
      INTEGER         , INTENT(  OUT) :: info
      !!! Return status, 0 on success
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
