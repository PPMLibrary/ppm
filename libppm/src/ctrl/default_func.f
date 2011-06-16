#define WRAP(a) a
    DO i=1,WRAP(DTYPE)_args_i
       IF (ASSOCIATED(WRAP(DTYPE)_args(i)%default_func)) THEN
          IF (.NOT. WRAP(DTYPE)_args(i)%default_func(WRAP(DTYPE)_args(i)%variable)) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_argument, caller, &
                  'Error applying default functions!', __LINE__, info)
             GOTO 9999
          END IF
       END IF
    END DO
#undef DTYPE
