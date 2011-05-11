#define WRAP(a) a
    DO i=1,WRAP(DTYPE)_args_i
#ifdef ARRAY
       IF (ANY(WRAP(DTYPE)_args(i)%variable .GT. WRAP(DTYPE)_args(i)%max) .OR. &
           ANY(WRAP(DTYPE)_args(i)%variable .LT. WRAP(DTYPE)_args(i)%min)) THEN
#else
       IF (WRAP(DTYPE)_args(i)%variable .GT. WRAP(DTYPE)_args(i)%max .OR. &
            WRAP(DTYPE)_args(i)%variable .LT. WRAP(DTYPE)_args(i)%min) THEN
#endif
          WRITE (cvar,*) 'Argument ', WRAP(DTYPE)_args(i)&
               &%name(1:LEN_TRIM(WRAP(DTYPE)_args(i)%name)), &
               ' fails min max check!'
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_argument, caller, cvar, __LINE__, info)
          GOTO 9999
       END IF
    END DO
#undef DTYPE
