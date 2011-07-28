#define WRAP(a) a
          DO i=1,WRAP(DTYPE)_args_i
             IF (WRAP(DTYPE)_args(i)%ctrl_name_set) THEN
                cvar = WRAP(DTYPE)_args(i)%ctrl_name &
                     (1:LEN_TRIM(WRAP(DTYPE)_args(i)%ctrl_name))
                CALL UpperCase(cvar, LEN_TRIM(cvar), info)
                IF (carg .EQ. cvar) THEN
#if defined(STRING) && !defined(ARRAY)
                   WRAP(DTYPE)_args(i)%variable = ADJUSTL(cvalue)
#else
                   READ (cvalue, *, IOSTAT=ios, ERR=200) WRAP(DTYPE)_args(i)%variable
#endif
                   CYCLE var_loop
                END IF
             END IF
          END DO
