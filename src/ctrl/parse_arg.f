#define WRAP(a) a
          ! find match
          DO i=1,WRAP(DTYPE)_args_i
             ! short flag
             ok = .FALSE.
             IF (WRAP(DTYPE)_args(i)%flag_set) THEN
#ifndef __LOGICAL
                CALL find_flag(WRAP(DTYPE)_args(i)%flag(1:2), ok, value, err)
#else
                CALL find_flag(WRAP(DTYPE)_args(i)%flag(1:2), ok)
#endif
             END IF
             IF ((.NOT.ok).AND.(WRAP(DTYPE)_args(i)%long_flag_set)) THEN
                ! long flag
#ifndef __LOGICAL
                CALL find_flag(WRAP(DTYPE)_args(i) &
                &    %long_flag(1:LEN_TRIM(WRAP(DTYPE)_args(i) &
                &    %long_flag)), ok, value, err)
#else
                CALL find_flag(WRAP(DTYPE)_args(i) &
                &    %long_flag(1:LEN_TRIM(WRAP(DTYPE)_args(i) &
                &    %long_flag)), ok)
#endif
             END IF
             IF (ok) THEN
#ifndef __LOGICAL
                ! match found - convert the arg
                READ (value,*,IOSTAT=ios) WRAP(DTYPE)_args(i)%variable
                ! check for failure
                IF (ios.NE.0) THEN
                   WRITE (cvar,*) 'Invalid argument (', value(1:LEN_TRIM(value)), ') for arg ', &
                   & WRAP(DTYPE)_args(i)%name(1:LEN_TRIM(WRAP(DTYPE)_args(i)%name))

                   fail(cvar,ppm_error=ppm_error_fatal)
                END IF
#else
                IF (WRAP(DTYPE)_args(i)%vtype .EQV. enabling_flag) THEN
                   WRAP(DTYPE)_args(i)%variable = .TRUE.
                ELSE
                   WRAP(DTYPE)_args(i)%variable = .FALSE.
                END IF
#endif
                WRAP(DTYPE)_args(i)%clf_supplied = .TRUE.
             ELSE IF (err) THEN
                WRITE (cvar,*) 'Flag ', value(1:LEN_TRIM(value)), &
                & ' has to be supplied with a value.'

                fail(cvar,ppm_error=ppm_error_fatal)
             END IF
          END DO
#undef DTYPE
#undef __LOGICAL
