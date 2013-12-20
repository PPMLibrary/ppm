#define WRAP(a) a
    DO i=1,WRAP(DTYPE)_args_i
       IF (WRAP(DTYPE)_args(i)%default_func_set) THEN
          IF (.NOT. WRAP(DTYPE)_args(i)%default_func(WRAP(DTYPE)_args(i)%variable)) THEN
             IF (WRAP(DTYPE)_args(i)%default_set) THEN
                WRAP(DTYPE)_args(i)%variable = WRAP(DTYPE)_args(i)%default
             ELSE
                info = ppm_error_fatal
                CALL ppm_error(ppm_err_argument, caller, &
                &    'Error applying default functions!', __LINE__, info)
                GOTO 9999
             END IF
          END IF
       END IF
    END DO
#undef DTYPE
