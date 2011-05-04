#define WRAP(a) a
             ! find match
             DO i=1,WRAP(DTYPE)_args_i
                ! short flag
                ok = .FALSE.
                IF (WRAP(DTYPE)_args(i)%flag_set) THEN
#if defined(ARRAY) || !defined(BOOL)
                   CALL find_flag(WRAP(DTYPE)_args(i)%flag(1:2), ok, value)
#else
                   CALL find_flag(WRAP(DTYPE)_args(i)%flag(1:2), ok)
#endif
                END IF
                IF ((.NOT. ok) .AND. (WRAP(DTYPE)_args(i)%long_flag_set)) THEN
                   ! long flag
#if defined(ARRAY) || !defined(BOOL)
                   CALL find_flag(WRAP(DTYPE)_args(i) &
                        %long_flag(1:LEN_TRIM(WRAP(DTYPE)_args(i) &
                        %long_flag)), ok, value)
#else
                   CALL find_flag(WRAP(DTYPE)_args(i) &
                        %long_flag(1:LEN_TRIM(WRAP(DTYPE)_args(i) &
                        %long_flag)), ok)
#endif
                END IF
                IF (ok) THEN
#if defined(ARRAY) || !defined(BOOL)
                   ! match found - convert the arg
                   READ (value,*,IOSTAT=ios) WRAP(DTYPE)_args(i)%variable
                   ! check for failure
                   IF (ios .NE. 0) THEN
                      WRITE (*,*) 'Invalid argument (', value(1:LEN_TRIM(value)), ') for arg ', &
                           WRAP(DTYPE)_args(i)%name(1:LEN_TRIM(WRAP(DTYPE)_args(i)%name))
                      info = 1
                      GOTO 9999
                   END IF
#else
                   IF (WRAP(DTYPE)_args(i)%type .EQV. enabling_flag) THEN
                      WRAP(DTYPE)_args(i)%variable = .TRUE.
                   ELSE
                      WRAP(DTYPE)_args(i)%variable = .FALSE.
                   END IF
#endif
                END IF
             END DO
