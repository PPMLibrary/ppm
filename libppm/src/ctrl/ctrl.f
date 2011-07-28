#define WRAP(a) a
       DO j=1,WRAP(DTYPE)_args_i
          IF (WRAP(DTYPE)_args(j)%group   .EQ. k .AND. &
              WRAP(DTYPE)_args(j)%group_i .EQ. i) THEN

             IF (WRAP(DTYPE)_args(j)%ctrl_name_set) THEN

                WRITE (*,'(A)',advance='no') &
                     WRAP(DTYPE)_args(j)%ctrl_name &
                     (1:LEN_TRIM(WRAP(DTYPE)_args(j)%ctrl_name))

                IF (WRAP(DTYPE)_args(j)%default_set) THEN
#ifdef __STRING
#ifdef ARRAY
                   WRITE (*,'(A)',advance='no') '= '
                   DO l=LBOUND(WRAP(DTYPE)_args(j)%default,1), &
                        UBOUND(WRAP(DTYPE)_args(j)%default,1)
                      WRITE (*,'(2A)',advance='no') &
                           WRAP(DTYPE)_args(j)%default(l) &
                           (1:LEN_TRIM(WRAP(DTYPE)_args(j)%default(l))), &
                           ', '
                   END DO
                   WRITE (*,'(A)') ''
#else
                   WRITE (*,*) "= ", WRAP(DTYPE)_args(j)%default &
                        (1:LEN_TRIM(WRAP(DTYPE)_args(j)%default))
#endif
#else
                   WRITE (*,*) "=", WRAP(DTYPE)_args(j)%default
#endif
                END IF
             END IF

             CYCLE var_loop
          END IF
       END DO
#undef DTYPE
#undef __STRING
