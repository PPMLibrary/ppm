#define WRAP(a) a
                DO j=1,WRAP(DTYPE)_args_i
                   IF (WRAP(DTYPE)_args(j)%group .EQ. k  .AND. &
                   &  WRAP(DTYPE)_args(j)%group_i .EQ. i .AND. &
                   &  (WRAP(DTYPE)_args(j)%flag_set .OR.       &
                   &  WRAP(DTYPE)_args(j)%long_flag_set)) THEN
                      IF (WRAP(DTYPE)_args(j)%help_set) THEN
                         WRITE (cbuf,'(3A)') "  ", WRAP(DTYPE)_args(j)%name(1:28),"  "
                         CALL break_help(WRAP(DTYPE)_args(j)%help  &
                         & (1:LEN_TRIM(WRAP(DTYPE)_args(j)%help)), &
                         & 50, "                                ", &
                         & cbuf(1:LEN_TRIM(cbuf)+2),caller)
                      ELSE
                         WRITE (cbuf,'(3A)') "  ", WRAP(DTYPE)_args(j)%name(1:28),""
                         stdout(cbuf)
                      ENDIF
                      cbuf=''

                      IF (WRAP(DTYPE)_args(j)%flag_set) THEN
                         WRITE (scratch,*) TRIM(cbuf),         &
                         & "                  short flag :  ", &
#if defined(ARRAY)
                         & WRAP(DTYPE)_args(j)%flag(1:2), " <v1,v2,...>"
#elif !defined(__LOGICAL)
                         & WRAP(DTYPE)_args(j)%flag(1:2), " <value>"
#else
                         & WRAP(DTYPE)_args(j)%flag(1:2)
#endif
                         WRITE (cbuf,'(2A)') TRIM(scratch),''
                         stdout(cbuf)
                         cbuf=''
                      ENDIF

                      IF (WRAP(DTYPE)_args(j)%long_flag_set) THEN
                         WRITE (scratch,*) TRIM(cbuf),          &
                         &  "                   long flag :  ", &
#if defined(ARRAY)
                         &  WRAP(DTYPE)_args(j)%long_flag(1:(LEN_TRIM(WRAP(DTYPE)_args(j)%long_flag))), &
                         &  " <v1,v2,...>"
#elif !defined(__LOGICAL)
                         &  WRAP(DTYPE)_args(j)%long_flag(1:(LEN_TRIM(WRAP(DTYPE)_args(j)%long_flag))), &
                         &  " <value>"
#else
                         &  WRAP(DTYPE)_args(j)%long_flag(1:(LEN_TRIM(WRAP(DTYPE)_args(j)%long_flag)))
#endif
                         WRITE (cbuf,'(2A)') TRIM(scratch),''
                         stdout(cbuf)
                         cbuf=''
                      ENDIF

                      IF (WRAP(DTYPE)_args(j)%ctrl_name_set) THEN
                         WRITE (scratch,*) TRIM(cbuf),          &
                         &  "           control file name :  ", &
                         &  WRAP(DTYPE)_args(j)%ctrl_name(1:LEN_TRIM(WRAP(DTYPE)_args(j)%ctrl_name)), &
#ifdef ARRAY
                         &  " = <v1,v2,...>"
#else
                         &  " = <value>"
#endif
                         WRITE (cbuf,'(2A)') TRIM(scratch),''
                         stdout(cbuf)
                         cbuf=''
                      ENDIF

                      IF (WRAP(DTYPE)_args(j)%default_set) THEN
#ifdef __STRING
#ifdef ARRAY
                         WRITE (scratch,'(2A)')TRIM(cbuf), &
                         & "               default value : "
                         WRITE (cbuf,'(A)') TRIM(scratch)
                         DO l=LBOUND(WRAP(DTYPE)_args(j)%default,1), &
                         &    UBOUND(WRAP(DTYPE)_args(j)%default,1)
                            WRITE (scratch,'(3A)') " ", TRIM(cbuf), &
                            & WRAP(DTYPE)_args(j)%default(l)        &
                            & (1:LEN_TRIM(WRAP(DTYPE)_args(j)%default(l)))
                            WRITE (cbuf,'(2A)') TRIM(scratch),','
                         END DO
                         WRITE (scratch,'(A)')cbuf(1:LEN_TRIM(cbuf)-1)
                         WRITE (cbuf,'(2A)') TRIM(scratch),''
                         stdout(cbuf)
                         cbuf=''
#else
                         WRITE (scratch,'(2A)')TRIM(cbuf), &
                         & "               default value : "
                         WRITE (cbuf,'(4A)') TRIM(scratch), ' ', &
                         & WRAP(DTYPE)_args(j)%default           &
                         & (1:LEN_TRIM(WRAP(DTYPE)_args(j)%default)),''
                         stdout(cbuf)
                         cbuf=''
#endif
#else
                         WRITE(scratch, *) WRAP(DTYPE)_args(j)%default
                         scratch = ADJUSTL(scratch)
                         WRITE(cbuf1,'(3A)')TRIM(cbuf),        &
                         & "                default value : ", &
                         & scratch(1:LEN_TRIM(scratch))
                         WRITE (cbuf,'(2A)') TRIM(cbuf1),''
                         stdout(cbuf)
                         cbuf=''
#endif
                      END IF
#if defined(__INTEGER) || defined(__LONGINT) || defined(__SINGLE) || defined(__DOUBLE)
                      IF (WRAP(DTYPE)_args(j)%min_set) THEN
                         WRITE(scratch, *) WRAP(DTYPE)_args(j)%min
                         scratch = ADJUSTL(scratch)
                         WRITE(cbuf1,'(3A)')TRIM(cbuf),        &
                         & "               minimum value :  ", &
                         & scratch(1:LEN_TRIM(scratch))
                         WRITE (cbuf,'(2A)') TRIM(cbuf1),''
                         stdout(cbuf)
                         cbuf=''
                      END IF
                      IF (WRAP(DTYPE)_args(j)%max_set) THEN
                         WRITE(scratch, *) WRAP(DTYPE)_args(j)%max
                         scratch = ADJUSTL(scratch)
                         WRITE(cbuf1,'(3A)')TRIM(cbuf),        &
                         & "               maximum value :  ", &
                         & scratch(1:LEN_TRIM(scratch))
                         WRITE (cbuf,'(2A)') TRIM(cbuf1),''
                         stdout(cbuf)
                         cbuf=''
                      END IF
#endif
                      WRITE(cbuf1,'(2A)') TRIM(cbuf),""
                      WRITE (cbuf,'(A)') TRIM(cbuf1)
                      stdout(cbuf)
                      cbuf=''
                      CYCLE group_loop
                   END IF
                END DO
#undef DTYPE
#undef __INTEGER
#undef __LONGINT
#undef __SINGLE
#undef __DOUBLE
#undef __LOGICAL
#undef __STRING
