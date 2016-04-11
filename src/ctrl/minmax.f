#define WRAP(a) a
          DO i=1,WRAP(DTYPE)_args_i
             IF (WRAP(DTYPE)_args(i)%min_set) THEN
#ifdef ARRAY
#ifdef __STRING
                IF (ANY(LEN_TRIM(WRAP(DTYPE)_args(i)%variable) .LT. WRAP(DTYPE)_args(i)%min)) THEN
#else
                IF (ANY(WRAP(DTYPE)_args(i)%variable .LT. WRAP(DTYPE)_args(i)%min)) THEN
#endif
#else
#ifdef __STRING
                IF (LEN_TRIM(WRAP(DTYPE)_args(i)%variable) .LT. WRAP(DTYPE)_args(i)%min) THEN
#else
                IF (WRAP(DTYPE)_args(i)%variable .LT. WRAP(DTYPE)_args(i)%min) THEN
#endif
#endif
                   WRITE (cvar,*) 'Argument ', WRAP(DTYPE)_args(i) &
                   & %name(1:LEN_TRIM(WRAP(DTYPE)_args(i)%name)),  &
                   & ' fails min check!'

                   fail(cvar,ppm_error=ppm_error_fatal)
                ENDIF
             ENDIF
             IF (WRAP(DTYPE)_args(i)%max_set) THEN
#ifdef ARRAY
#ifdef __STRING
                IF (ANY(LEN_TRIM(WRAP(DTYPE)_args(i)%variable) .GT. WRAP(DTYPE)_args(i)%max)) THEN
#else
                IF (ANY(WRAP(DTYPE)_args(i)%variable .GT. WRAP(DTYPE)_args(i)%max)) THEN
#endif
#else
#ifdef __STRING
                IF (LEN_TRIM(WRAP(DTYPE)_args(i)%variable) .GT. WRAP(DTYPE)_args(i)%max) THEN
#else
                IF (WRAP(DTYPE)_args(i)%variable .GT. WRAP(DTYPE)_args(i)%max) THEN
#endif
#endif
                   WRITE (cvar,*) 'Argument ', WRAP(DTYPE)_args(i)&
                   & %name(1:LEN_TRIM(WRAP(DTYPE)_args(i)%name)), &
                   & ' fails max check!'

                   fail(cvar,ppm_error=ppm_error_fatal)
                ENDIF
             ENDIF
          ENDDO
#undef DTYPE
#undef __STRING
