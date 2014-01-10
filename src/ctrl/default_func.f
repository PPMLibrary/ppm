#define WRAP(a) a
          DO i=1,WRAP(DTYPE)_args_i
             IF (WRAP(DTYPE)_args(i)%default_func_set) THEN
                IF (.NOT. WRAP(DTYPE)_args(i)%default_func(WRAP(DTYPE)_args(i)%variable)) THEN
                   IF (WRAP(DTYPE)_args(i)%default_set) THEN
                      WRAP(DTYPE)_args(i)%variable = WRAP(DTYPE)_args(i)%default
                   ELSE
                      fail('Error applying default functions!',ppm_error=ppm_error_fatal)
                   END IF
                END IF
             END IF
          END DO
#undef DTYPE
