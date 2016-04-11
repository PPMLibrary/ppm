#define WRAP(a) a
          DO i=1,WRAP(DTYPE)_args_i
             IF (WRAP(DTYPE)_args(i)%validator_set) THEN
                IF (.NOT. WRAP(DTYPE)_args(i)%validator(WRAP(DTYPE)_args(i)%variable)) THEN
                   WRITE (cvar,*) 'Argument ', WRAP(DTYPE)_args(i)%name(1:LEN_TRIM(WRAP(DTYPE)_args(i)%name)), &
                   &  ' fails validator check!'

                   fail(cvar,ppm_error=ppm_error_fatal)
                ENDIF
             ENDIF
          ENDDO
#undef DTYPE
