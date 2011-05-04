#define WRAP(a) a
       DO j=1,WRAP(DTYPE)_args_i

          IF (WRAP(DTYPE)_args(j)%group   .EQ. k .AND. &
              WRAP(DTYPE)_args(j)%group_i .EQ. i) THEN

             IF (WRAP(DTYPE)_args(j)%ctrl_name_set) THEN

                WRITE (*,'(A,A15)',advance='no') "#  ", &
                     WRAP(DTYPE)_args(j)%name(1:LEN_TRIM(WRAP(DTYPE)_args(j)%name))

                IF (WRAP(DTYPE)_args(j)%help_set) &
                     WRITE (*,'(A,A)') ": ", &
                     WRAP(DTYPE)_args(j)%help(1:LEN_TRIM(WRAP(DTYPE)_args(j)%help))

                IF (WRAP(DTYPE)_args(j)%flag_set .OR. &
                     WRAP(DTYPE)_args(j)%long_flag_set) THEN
                   WRITE (*,'(A)',advance='no') "#                   "

                   IF (WRAP(DTYPE)_args(j)%flag_set) THEN
#if defined(ARRAY)
                      WRITE(*,'(A,A)',advance='no') &
                           WRAP(DTYPE)_args(j)%flag(1:2), " {v1,v2,...}"
#elif !defined(BOOL)
                      WRITE(*,'(A,A)',advance='no') &
                           WRAP(DTYPE)_args(j)%flag(1:2), " {value}"
#else
                      WRITE(*,'(A)',advance='no') &
                           WRAP(DTYPE)_args(j)%flag(1:2)
#endif
                      IF (WRAP(DTYPE)_args(j)%long_flag_set) &
                           WRITE(*,'(A)',advance='no') ', '
                   END IF

                   IF (WRAP(DTYPE)_args(j)%long_flag_set) &
                        WRITE(*,*) WRAP(DTYPE)_args(j)%long_flag &
#if defined(ARRAY)
                        (1:(LEN_TRIM(WRAP(DTYPE)_args(j)%long_flag))), " {v1,v2,...}"
#elif !defined(BOOL)
                        (1:(LEN_TRIM(WRAP(DTYPE)_args(j)%long_flag))), " {value}"
#else
                        (1:(LEN_TRIM(WRAP(DTYPE)_args(j)%long_flag)))
#endif
                END IF

#ifndef STRING
#ifndef BOOL
                IF (WRAP(DTYPE)_args(j)%min_set .OR. &
                     WRAP(DTYPE)_args(j)%max_set) THEN
                   WRITE (*,'(A)',advance='no') "#                   "
#ifdef ARRAY
                   WRITE (*,*) WRAP(DTYPE)_args(j)%min, ' < {v1,v2,...} < ', &
#else
                   WRITE (*,*) WRAP(DTYPE)_args(j)%min, ' < {value} < ', &
#endif
                        WRAP(DTYPE)_args(j)%max
                END IF
#endif
#endif
             END IF
             CYCLE comment_loop
          END IF
       END DO
