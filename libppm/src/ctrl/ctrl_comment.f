#define WRAP(a) a
       DO j=1,WRAP(DTYPE)_args_i

          IF (WRAP(DTYPE)_args(j)%group   .EQ. k .AND. &
              WRAP(DTYPE)_args(j)%group_i .EQ. i) THEN

             IF (WRAP(DTYPE)_args(j)%ctrl_name_set) THEN

                WRITE (*,'(A,A15)',advance='no') "#  ", &
                     WRAP(DTYPE)_args(j)%name(1:LEN_TRIM(WRAP(DTYPE)_args(j)%name))
                in_line = .TRUE.

                IF (WRAP(DTYPE)_args(j)%help_set) THEN
                   WRITE (*,'(A,A)') ": ", &
                        WRAP(DTYPE)_args(j)%help(1:LEN_TRIM(WRAP(DTYPE)_args(j)%help))
                   in_line = .FALSE.
                END IF

                IF (WRAP(DTYPE)_args(j)%flag_set .OR. &
                     WRAP(DTYPE)_args(j)%long_flag_set) THEN
                   IF (in_line) THEN
                      WRITE (*,'(A)',advance='no') ": "
                   ELSE
                      WRITE (*,'(A)',advance='no') "#                   "
                   END IF
                   in_line = .FALSE.
                   IF (WRAP(DTYPE)_args(j)%flag_set) THEN
#if defined(ARRAY)
                      WRITE(*,'(A,A)',advance='no') &
                           WRAP(DTYPE)_args(j)%flag(1:2), " {v1,v2,...}"
#elif !defined(__LOGICAL)
                      WRITE(*,'(A,A)',advance='no') &
                           WRAP(DTYPE)_args(j)%flag(1:2), " {value}"
#else
                      WRITE(*,'(A)',advance='no') &
                           WRAP(DTYPE)_args(j)%flag(1:2)
#endif
                      IF (WRAP(DTYPE)_args(j)%long_flag_set) THEN
                         WRITE (*,'(A)',advance='no') ', '
                      ELSE
                         WRITE (*,*) ''
                      END IF
                   END IF

                   IF (WRAP(DTYPE)_args(j)%long_flag_set) &
                        WRITE(*,*) WRAP(DTYPE)_args(j)%long_flag &
#if defined(ARRAY)
                        (1:(LEN_TRIM(WRAP(DTYPE)_args(j)%long_flag))), " {v1,v2,...}"
#elif !defined(__LOGICAL)
                        (1:(LEN_TRIM(WRAP(DTYPE)_args(j)%long_flag))), " {value}"
#else
                        (1:(LEN_TRIM(WRAP(DTYPE)_args(j)%long_flag)))
#endif
                END IF

#if defined(__INTEGER) || defined(__LONGINT) || defined(__SINGLE) || defined(__DOUBLE)
                IF (WRAP(DTYPE)_args(j)%min_set .OR. &
                     WRAP(DTYPE)_args(j)%max_set) THEN
                   IF (in_line) THEN
                      WRITE (*,'(A)',advance='no') ": "
                   ELSE
                      WRITE (*,'(A)',advance='no') "#                   "
                   END IF
                   in_line = .FALSE.
                   IF (WRAP(DTYPE)_args(j)%min_set) THEN
                      WRITE (scratch, *) WRAP(DTYPE)_args(j)%min
                      scratch = ADJUSTL(scratch)
                      WRITE (*,'(A,A)',advance='no') scratch(1:LEN_TRIM(scratch)) , ' <= '
                   END IF
#ifdef ARRAY
                   WRITE (*,'(A)',advance='no') '{v1,v2,...}'
#else
                   WRITE (*,'(A)',advance='no') '{value}'
#endif
                   IF (WRAP(DTYPE)_args(j)%max_set) THEN
                      WRITE (scratch, *) WRAP(DTYPE)_args(j)%max
                      scratch = ADJUSTL(scratch)
                      WRITE (*,'(A,A)',advance='no') ' <= ', scratch(1:LEN_TRIM(scratch))
                   END IF
                   WRITE(*,*) ''
                END IF
#endif
                IF (in_line) WRITE (*,*) ''
             END IF
             CYCLE comment_loop
          END IF
       END DO
#undef DTYPE
#undef __INTEGER
#undef __LONGINT
#undef __SINGLE
#undef __DOUBLE
#undef __LOGICAL
#undef __STRING
