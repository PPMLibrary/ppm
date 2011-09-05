
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
    opts%D_needs_gradients = .TRUE.
    opts%scale_D = 1._MK
    opts%rcp_over_D = 3._MK
    opts%nb_width = 1._MK
    opts%nb_width2 = 1._MK
    opts%adaptivity_criterion = 80.0_MK
    opts%attractive_radius0 = 0.5_MK
    opts%fuse_radius = 0.40_MK
    opts%param_nb = 1._MK
    opts%diff_eq => NULL()
    opts%order_approx = 4
    opts%maximum_D = 0.2_MK
    opts%minimum_D = 0.05_MK    
    opts%max_factor_aniso = 5.00_MK 
    opts%nb_grad_desc_steps = 0
    opts%nneigh_theo = 24
    opts%nneigh_critical = 20
    opts%nneigh_toobig = 2000
    opts%check_dcops = .FALSE.

    9999 CONTINUE

END SUBROUTINE sop_init_opts


SUBROUTINE sop_init_stats(stats,info)
    !!! constructor for sop statistics derived type
    USE ppm_module_error

    IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER  :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER  :: MK = ppm_kind_double
#endif
    TYPE(sop_t_stats), POINTER, INTENT(INOUT) :: stats
    INTEGER                  , INTENT(  OUT) :: info

    info = 0
    ALLOCATE(stats,STAT=info)
    IF (info.NE.0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,'sop_init_stats',       &
            &                  'allocation error',__LINE__,info)
        GOTO 9999
    ENDIF

    stats%nb_grad_desc_steps = 0
    stats%min_sv = HUGE(1._mk)

    9999 CONTINUE

END SUBROUTINE sop_init_stats
