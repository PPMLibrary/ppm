!emulating templating....
    !-------------------------------------------------------------------------
    ! local variables
    !-------------------------------------------------------------------------
    CHARACTER(LEN = ppm_char)               :: lcaller
    INTEGER                                 :: info


    IF (PRESENT(caller)) THEN
        lcaller = TRIM(ADJUSTL(caller))
    ELSE
        lcaller = 'ppm_DS_exists'
    ENDIF
    exists = .FALSE.
    !-------------------------------------------------------------------------
    ! Check arguments
    !-------------------------------------------------------------------------
    IF (id.LE.0 .OR. id.LT.this%min_id .OR. id.GT.this%max_id) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,lcaller,   &
            & 'Invalid id for this data structure, use create() first',&
            __LINE__,info)
        RETURN
    ENDIF


    info = ppm_error_error
    CALL ppm_error(ppm_err_argument,lcaller,   &
        & 'No data structure found, use create() first',&
        __LINE__,info)
    RETURN
