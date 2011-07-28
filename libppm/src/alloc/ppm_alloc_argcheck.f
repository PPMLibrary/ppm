      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_alloc_argcheck
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_alloc_argcheck(caller,iopt,ldl,dimension,info,ldu)
      !!! Checks the arguments of the ppm_alloc routines

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop

      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CHARACTER(LEN=*)               , INTENT(IN)   :: caller
      !!! name of calling subroutine
      INTEGER, DIMENSION(:)          , INTENT(IN)   :: ldl
      !!! Lower index limit in leading dim.
      INTEGER, OPTIONAL, DIMENSION(:), INTENT(IN)   :: ldu
      !!! Upper index limit in leading dim. (>ldl(1)).
      !!! OPTIONAL in ppm_alloc_*dl.f
      INTEGER                        , INTENT(IN)   :: iopt
      !!! Allocation mode. One of:
      !!!
      !!! * ppm_param_alloc_fit
      !!! * ppm_param_alloc_fit_preserve
      !!! * ppm_param_alloc_grow
      !!! * ppm_param_alloc_grow_preserve
      !!! * ppm_param_dealloc
      INTEGER                        , INTENT(IN)   :: dimension
      !!! Helps to determine the dimension of the caller subroutine (1-5)
      INTEGER                        , INTENT(OUT)  :: info
      !!! Returns status, 0 upon success.
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0
      INTEGER               :: i

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_alloc_argcheck',t0,info)
      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (iopt .NE. ppm_param_alloc_fit           .AND.                       &
     &    iopt .NE. ppm_param_alloc_fit_preserve  .AND.                       &
     &    iopt .NE. ppm_param_alloc_grow          .AND.                       &
     &    iopt .NE. ppm_param_alloc_grow_preserve .AND.                       &
     &    iopt .NE. ppm_param_dealloc) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,'unknown iopt',__LINE__,info)
        GOTO 9999
      ENDIF
      IF (iopt .NE. ppm_param_dealloc) THEN
        IF (PRESENT(ldu)) THEN
            DO i=1,dimension
                IF (ldl(i) .GT. ldu(i)) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,caller,  &
     &              'ldu() must be >= ldl()',__LINE__,info)
                GOTO 9999
                ENDIF
            ENDDO
        ELSE
            DO i=1,dimension
                IF (ldl(i) .LT. 0) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,caller,  &
     &              'ldl() must be >= 0',__LINE__,info)
                GOTO 9999
                ENDIF
            ENDDO
        ENDIF
      ENDIF
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_alloc_argcheck',t0,info)
      RETURN
      END SUBROUTINE ppm_alloc_argcheck
