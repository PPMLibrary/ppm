
SUBROUTINE sop_init_opts(opts,info)
    !!! constructor for sop options derived type
    USE ppm_module_error

    IMPLICIT NONE

#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER  :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER  :: MK = ppm_kind_double
#endif
    TYPE(sop_t_opts), POINTER, INTENT(INOUT) :: opts
    INTEGER                  , INTENT(  OUT) :: info

    info = 0

    ALLOCATE(opts,STAT=info)
    IF (info.NE.0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,'sop_init_opts',       &
            &                  'allocation error',__LINE__,info)
        GOTO 9999
    ENDIF

    opts%level_set = .FALSE.
    opts%scale_D = 1._MK
    opts%nb_width = 1._MK
    opts%nb_width2 = 1._MK
    opts%adaptivity_criterion = 1._MK
    opts%fuse_radius = 0._MK
    opts%param_nb = 1._MK
    opts%diff_eq => NULL()
    opts%order_approx = 4
    opts%maximum_D = 1._MK
    opts%minimum_D = 0.01_MK     
    opts%nb_grad_desc_steps = 0

    9999 CONTINUE

END SUBROUTINE sop_init_opts
