      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_io_inquire
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine inquires the state of a given ppm I/O
      !                 unit.
      !
      !  Input        : iUnit          (I) IO unit to be inquired
      !
      !  Input/output :
      !
      !  Output       : stat           (I) return code. 0 on success.
      !                                    OPTIONAL. 
      !                 isopen         (L) .TRUE. of the unit is open on
      !                                    the local processor. .FALSE.
      !                                    otherwise. OPTIONAL.
      !                 mode           (I) I/O mode. OPTIONAL. One of:
      !                                       ppm_param_io_centralized
      !                                       ppm_param_io_distributed
      !                                       ppm_param_undefined
      !                 pfmt           (I) I/O format. OPTIONAL. One of:
      !                                       ppm_param_io_ascii
      !                                       ppm_param_io_binary
      !                                       ppm_param_undefined
      !
      !  Remarks      : All optput arguments are OPTIONAL. It is adivised
      !                 to use explicit argument naming when using this
      !                 routine. Example:
      !                     CALL ppm_inquire(20,MODE=imode,STAT=info)
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_io_inquire.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.6  2004/10/01 16:33:34  ivos
      !  cosmetics.
      !
      !  Revision 1.5  2004/10/01 16:09:01  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.4  2004/07/26 07:45:26  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.3  2004/07/16 14:46:26  ivos
      !  Added check for ppm_initialized.
      !
      !  Revision 1.2  2004/05/14 14:47:25  ivos
      !  Tested. Added check that iUnit is not .LE. 0.
      !
      !  Revision 1.1  2004/05/06 10:44:59  ivos
      !  Preliminary check-in. Not tested yet.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_io_inquire(iUnit,isopen,mode,pfmt,stat)

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_io
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER             , INTENT(IN   ) :: iUnit
      INTEGER  , OPTIONAL , INTENT(  OUT) :: stat,pfmt,mode
      LOGICAL  , OPTIONAL , INTENT(  OUT) :: isopen
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)            :: t0
      LOGICAL                          :: lopen
      INTEGER                          :: info
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_io_inquire',t0,info)

      !-------------------------------------------------------------------------
      !  Check argument
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_io_inquire',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF
      IF (iUnit .LE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_argument,'ppm_io_inquire',    &
     &        'Unit number needs to be > 0',__LINE__,info)
          IF (PRESENT(isopen)) isopen = .FALSE.
          IF (PRESENT(pfmt)) pfmt = ppm_param_undefined
          IF (PRESENT(mode)) mode = ppm_param_undefined
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Check if the unit number is open
      !-------------------------------------------------------------------------
      lopen = .FALSE.
      IF(ASSOCIATED(ppm_io_unit)) THEN
          IF (SIZE(ppm_io_unit,1) .GE. iUnit) THEN
              IF (ppm_io_unit(iUnit) .GT. 0) lopen = .TRUE.
          ENDIF
      ENDIF
      IF (PRESENT(isopen)) isopen = lopen

      !-------------------------------------------------------------------------
      !  Get the other information if present
      !-------------------------------------------------------------------------
      IF (lopen) THEN
          IF (PRESENT(pfmt)) pfmt = ppm_io_format(ppm_io_unit(iUnit))
          IF (PRESENT(mode)) mode = ppm_io_mode(ppm_io_unit(iUnit))
      ELSE
          IF (PRESENT(pfmt)) pfmt = ppm_param_undefined
          IF (PRESENT(mode)) mode = ppm_param_undefined
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      IF (PRESENT(stat)) stat = info
      CALL substop('ppm_io_inquire',t0,info)
      RETURN
      END SUBROUTINE ppm_io_inquire
