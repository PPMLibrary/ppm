#define WRAP(a) a
       DO i=1,WRAP(DTYPE)_args_i
          IF (ASSOCIATED(WRAP(DTYPE)_args(i)%validator)) THEN
             IF (.NOT. WRAP(DTYPE)_args(i)%validator(WRAP(DTYPE)_args(i)%variable)) THEN
                WRITE (*,*) 'Argument ', WRAP(DTYPE)_args(i)%name(1:LEN_TRIM(WRAP(DTYPE)_args(i)%name)), &
                     ' fails validator check!'
                info = 1
                GOTO 9999
             END IF
          END IF
       END DO
