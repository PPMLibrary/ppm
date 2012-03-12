#define WRAP(a) a
       DO j=1,WRAP(DTYPE)_args_i
          IF (WRAP(DTYPE)_args(j)%group   .EQ. k .AND. &
              WRAP(DTYPE)_args(j)%group_i .EQ. i) THEN

             IF (WRAP(DTYPE)_args(j)%ctrl_name_set) THEN
                WRITE (scratch, *) group_max_len(k)
                WRITE (*,'(A' // scratch(1:LEN_TRIM(scratch)) // ')',advance='no') &
                     WRAP(DTYPE)_args(j)%ctrl_name &
                     (1:LEN_TRIM(WRAP(DTYPE)_args(j)%ctrl_name))

!                 IF (WRAP(DTYPE)_args(j)%default_set) THEN
#ifdef __STRING
#ifdef ARRAY
                   WRITE (*,'(A)',advance='no') ' = '
                   DO l=LBOUND(WRAP(DTYPE)_args(j)%variable,1), &
                        UBOUND(WRAP(DTYPE)_args(j)%variable,1)
                      WRITE (*,'(2A)',advance='no') &
                           WRAP(DTYPE)_args(j)%variable(l) &
                           (1:LEN_TRIM(WRAP(DTYPE)_args(j)%variable(l))), &
                           ', '
                   END DO
                   WRITE (*,'(A)') ''
#else
                   WRITE (*,*) "= ", WRAP(DTYPE)_args(j)%variable &
                        (1:LEN_TRIM(WRAP(DTYPE)_args(j)%variable))
#endif
#else
#ifdef ARRAY
                   WRITE (*,'(A)',advance='no') ' = '
                   DO l=LBOUND(WRAP(DTYPE)_args(j)%variable,1), &
                        UBOUND(WRAP(DTYPE)_args(j)%variable,1)
                      WRITE(scratch,*) WRAP(DTYPE)_args(j)%variable(l)
                      WRITE (*,'(A)',advance='no') &
                           scratch(1:LEN_TRIM(scratch))
                      IF (l .EQ. UBOUND(WRAP(DTYPE)_args(j)%variable,1)) THEN
                         WRITE (*,*) ''
                      ELSE
                         WRITE (*,'(A)',advance='no') ', '
                      END IF
                   END DO
#else
                   WRITE (scratch, *) WRAP(DTYPE)_args(j)%variable
                   scratch = ADJUSTL(scratch)
                   WRITE (*,*) "= ", scratch(1:LEN_TRIM(scratch))
#endif
#endif
                END IF
!              END IF

             CYCLE var_loop
          END IF
       END DO
#undef DTYPE
#undef __STRING
