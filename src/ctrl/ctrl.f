#define WRAP(a) a
                DO j=1,WRAP(DTYPE)_args_i
                   IF (WRAP(DTYPE)_args(j)%group  .EQ. k .AND. &
                   &  WRAP(DTYPE)_args(j)%group_i .EQ. i) THEN

                     IF (WRAP(DTYPE)_args(j)%ctrl_name_set) THEN
                        WRITE (scratch, *) group_max_len(k)
                        WRITE (cbuf,'(A' // scratch(1:LEN_TRIM(scratch)) // ')') &
                        & WRAP(DTYPE)_args(j)%ctrl_name                          &
                        & (1:LEN_TRIM(WRAP(DTYPE)_args(j)%ctrl_name))

!                       IF (WRAP(DTYPE)_args(j)%default_set) THEN
#ifdef __STRING
#ifdef ARRAY
                        WRITE (cbuf1,'(2A)')TRIM(cbuf),' = '
                        DO l=LBOUND(WRAP(DTYPE)_args(j)%variable,1), &
                        &    UBOUND(WRAP(DTYPE)_args(j)%variable,1)
                           WRITE (scratch,'(4A)') TRIM(cbuf1), " ", &
                           & WRAP(DTYPE)_args(j)%variable(l) &
                           & (1:LEN_TRIM(WRAP(DTYPE)_args(j)%variable(l))), ', '
                           WRITE (cbuf1,'(A)') TRIM(scratch)
                        ENDDO
                        WRITE (cbuf,'(2A)')cbuf1(1:LEN_TRIM(cbuf1)-1), ''
                        stdout(cbuf)
                        cbuf=''
#else
                        WRITE (cbuf1,*) TRIM(cbuf)," = ", &
                        & WRAP(DTYPE)_args(j)%variable   &
                        & (1:LEN_TRIM(WRAP(DTYPE)_args(j)%variable))
                        WRITE (cbuf,'(2A)')TRIM(cbuf1), ''
                        stdout(cbuf)
                        cbuf=''
#endif
#else
#ifdef ARRAY
                        WRITE (cbuf1,'(2A)') TRIM(cbuf),' = '
                        DO l=LBOUND(WRAP(DTYPE)_args(j)%variable,1), &
                        &    UBOUND(WRAP(DTYPE)_args(j)%variable,1)
                           WRITE(scratch,*) TRIM(cbuf1)," ", &
                           & WRAP(DTYPE)_args(j)%variable(l)
                           IF (l .EQ. UBOUND(WRAP(DTYPE)_args(j)%variable,1)) THEN
                              WRITE (cbuf,'(2A)') scratch(1:LEN_TRIM(scratch)),''
                              stdout(cbuf)
                              cbuf1=''
                           ELSE
                              WRITE (cbuf1,'(2A)') scratch(1:LEN_TRIM(scratch)), ', '
                           ENDIF
                        ENDDO
#else
                        WRITE (scratch, *) WRAP(DTYPE)_args(j)%variable
                        scratch = ADJUSTL(scratch)
                        WRITE (cbuf1,'(3A)') TRIM(cbuf), " = ", &
                        & scratch(1:LEN_TRIM(scratch))
                        WRITE (cbuf,'(2A)')TRIM(cbuf1),''
                        stdout(cbuf)
                        cbuf=''
#endif
#endif
                     ENDIF
!                 ENDIF

                     CYCLE var_loop
                  ENDIF
               ENDDO
#undef DTYPE
#undef __STRING
