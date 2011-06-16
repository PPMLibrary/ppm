#define WRAP(a) a
       DO j=1,WRAP(DTYPE)_args_i
          IF (WRAP(DTYPE)_args(j)%group .EQ. k .AND. &
              WRAP(DTYPE)_args(j)%group_i .EQ. i) THEN
             WRITE (*,'(/2A)',advance='no') "  ", WRAP(DTYPE)_args(j)%name(1:28)
             IF (WRAP(DTYPE)_args(j)%help_set) &
                  WRITE (*,*) "  ", &
                  WRAP(DTYPE)_args(j)%help(1:LEN_TRIM(WRAP(DTYPE)_args(j)%help))
             IF (WRAP(DTYPE)_args(j)%flag_set) &
                  WRITE (*,*) "                  short flag :  ", &
#if defined(ARRAY)
                  WRAP(DTYPE)_args(j)%flag(1:2), " <v1,v2,...>"
#elif !defined(BOOL)
                  WRAP(DTYPE)_args(j)%flag(1:2), " <value>"
#else
                  WRAP(DTYPE)_args(j)%flag(1:2)
#endif
             IF (WRAP(DTYPE)_args(j)%long_flag_set) &
                  WRITE (*,*) "                   long flag :  ", &
#if defined(ARRAY)
                  WRAP(DTYPE)_args(j)%long_flag(1:(LEN_TRIM(WRAP(DTYPE)_args(j)%long_flag))), &
                  " <v1,v2,...>"
#elif !defined(BOOL)
                  WRAP(DTYPE)_args(j)%long_flag(1:(LEN_TRIM(WRAP(DTYPE)_args(j)%long_flag))), &
                  " <value>"
#else
                  WRAP(DTYPE)_args(j)%long_flag(1:(LEN_TRIM(WRAP(DTYPE)_args(j)%long_flag)))
#endif
             IF (WRAP(DTYPE)_args(j)%ctrl_name_set) &
                  WRITE (*,*) "           control file name :  ", &
                  WRAP(DTYPE)_args(j)%ctrl_name(1:LEN_TRIM(WRAP(DTYPE)_args(j)%ctrl_name)), &
#ifdef ARRAY
                  " = <v1,v2,...>"
#else
                  " = <value>"
#endif
             IF (WRAP(DTYPE)_args(j)%default_set) THEN
                WRITE (*,'(A)',advance='no') "                default value : "
#ifdef STRING
#ifdef ARRAY
                WRITE (*,'(A)',advance='no') ' '
                DO l=LBOUND(WRAP(DTYPE)_args(j)%default,1), &
                     UBOUND(WRAP(DTYPE)_args(j)%default,1)
                   WRITE (*,'(2A)',advance='no') &
                        WRAP(DTYPE)_args(j)%default(l) &
                        (1:LEN_TRIM(WRAP(DTYPE)_args(j)%default(l))), &
                        ', '
                END DO
                WRITE (*,'(A)') ''
#else
                WRITE (*,*) WRAP(DTYPE)_args(j)%default &
                     (1:LEN_TRIM(WRAP(DTYPE)_args(j)%default))
#endif
#else
                WRITE(*,*) WRAP(DTYPE)_args(j)%default
#endif
             END IF
#ifndef STRING
#ifndef BOOL
             IF (WRAP(DTYPE)_args(j)%min_set) &
                  WRITE (*,*) "               minimum value :  ", &
                  WRAP(DTYPE)_args(j)%min
             IF (WRAP(DTYPE)_args(j)%max_set) &
                  WRITE (*,*) "               maximum value :  ", &
                  WRAP(DTYPE)_args(j)%max
#endif
#endif
             CYCLE group_loop
          END IF
       END DO
