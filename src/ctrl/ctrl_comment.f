#define WRAP(a) a
              DO j=1,WRAP(DTYPE)_args_i

                 IF (WRAP(DTYPE)_args(j)%group   .EQ. k .AND. &
                 &  WRAP(DTYPE)_args(j)%group_i .EQ. i) THEN

                    IF (WRAP(DTYPE)_args(j)%ctrl_name_set) THEN

                       WRITE (cbuf,'(A,A15)') "#  ", &
                       & WRAP(DTYPE)_args(j)%name(1:LEN_TRIM(WRAP(DTYPE)_args(j)%name))

                       in_line = .TRUE.

                       IF (WRAP(DTYPE)_args(j)%help_set) THEN
                          WRITE (scratch,'(2A)') TRIM(cbuf), ": "
                          CALL break_help(WRAP(DTYPE)_args(j)%help  &
                          & (1:LEN_TRIM(WRAP(DTYPE)_args(j)%help)), &
                          & 50, "#                   ",             &
                          & scratch(1:LEN_TRIM(scratch)+1),caller)
                          in_line = .FALSE.
                          cbuf=''
                       ENDIF

                       IF (WRAP(DTYPE)_args(j)%flag_set .OR.  &
                       &   WRAP(DTYPE)_args(j)%long_flag_set) THEN
                          IF (in_line) THEN
                             WRITE (scratch,'(2A)') TRIM(cbuf), ": "
                             lng=LEN_TRIM(cbuf)+2
                          ELSE
                             WRITE (scratch,'(2A)') TRIM(cbuf), "#                   "
                             lng=LEN_TRIM(cbuf)+20
                          ENDIF
                          WRITE (cbuf,'(A)') scratch
                          in_line = .FALSE.
                          IF (WRAP(DTYPE)_args(j)%flag_set) THEN
#if defined(ARRAY)
                             WRITE(scratch,'(3A)') cbuf(1:lng), &
                             & WRAP(DTYPE)_args(j)%flag(1:2), " {v1,v2,...}"
#elif !defined(__LOGICAL)
                             WRITE(scratch,'(3A)') cbuf(1:lng), &
                             & WRAP(DTYPE)_args(j)%flag(1:2), " {value}"
#else
                             WRITE(scratch,'(3A)') cbuf(1:lng), " ", &
                             & WRAP(DTYPE)_args(j)%flag(1:2)
#endif
                             IF (WRAP(DTYPE)_args(j)%long_flag_set) THEN
                                WRITE (cbuf,'(2A)') TRIM(scratch),', '
                                lng=LEN_TRIM(scratch)+2
                             ELSE
                                WRITE (cbuf,'(2A)') TRIM(scratch), ''
                                stdout(cbuf)
                                cbuf=''
                                lng=0
                             ENDIF
                          ENDIF

                          IF (WRAP(DTYPE)_args(j)%long_flag_set)  THEN
                             WRITE(cbuf1,'(2A)') cbuf(1:lng), " "
                             WRITE(scratch,*) cbuf1(1:lng+1),  &
                             & WRAP(DTYPE)_args(j)%long_flag   &
#if defined(ARRAY)
                             & (1:(LEN_TRIM(WRAP(DTYPE)_args(j)%long_flag))), " {v1,v2,...}"
#elif !defined(__LOGICAL)
                             & (1:(LEN_TRIM(WRAP(DTYPE)_args(j)%long_flag))), " {value}"
#else
                             & (1:(LEN_TRIM(WRAP(DTYPE)_args(j)%long_flag)))
#endif
                             WRITE (cbuf,'(2A)') TRIM(scratch), ''
                             stdout(cbuf)
                             cbuf=''
                             lng=0
                          ENDIF
                       ENDIF

#if defined(__INTEGER) || defined(__LONGINT) || defined(__SINGLE) || defined(__DOUBLE)
                       IF (WRAP(DTYPE)_args(j)%min_set .OR. &
                       &   WRAP(DTYPE)_args(j)%max_set) THEN
                          IF (in_line) THEN
                             WRITE(scratch,'(2A)') TRIM(cbuf), ": "
                             lng=LEN_TRIM(cbuf)+2
                          ELSE
                             WRITE(scratch,'(2A)') TRIM(cbuf), "#                   "
                             lng=LEN_TRIM(cbuf)+20
                          ENDIF
                          WRITE (cbuf,'(A)') scratch
                          in_line = .FALSE.
                          IF (WRAP(DTYPE)_args(j)%min_set) THEN
                             WRITE (scratch, *) WRAP(DTYPE)_args(j)%min
                             scratch = ADJUSTL(scratch)
                             WRITE (cbuf1,'(4A)') cbuf(1:lng),' ',scratch(1:LEN_TRIM(scratch)) , ' <= '
                             WRITE (cbuf,'(A)') cbuf1
                             lng=LEN_TRIM(cbuf)
                          ENDIF
#ifdef ARRAY
                          WRITE (scratch,'(3A)') cbuf(1:lng),' ','{v1,v2,...}'
#else
                          WRITE (scratch,'(3A)') cbuf(1:lng),' ','{value}'
#endif
                          WRITE (cbuf,'(A)') TRIM(scratch)

                          IF (WRAP(DTYPE)_args(j)%max_set) THEN
                             WRITE (scratch, *) WRAP(DTYPE)_args(j)%max
                             scratch = ADJUSTL(scratch)
                             WRITE (cbuf1,'(3A)') TRIM(cbuf),' <= ', scratch(1:LEN_TRIM(scratch))
                             WRITE (cbuf,'(A)') TRIM(cbuf1)
                          ENDIF
                          WRITE (scratch,'(A)')TRIM(cbuf)
                          WRITE (cbuf,'(2A)') TRIM(scratch), ''
                          stdout(cbuf)
                          cbuf=''
                       ENDIF
#endif
                       IF (in_line) THEN
                          stdout("")
                       ENDIF
                    ENDIF
                    CYCLE comment_loop
                 ENDIF
              ENDDO
#undef DTYPE
#undef __INTEGER
#undef __LONGINT
#undef __SINGLE
#undef __DOUBLE
#undef __LOGICAL
#undef __STRING
