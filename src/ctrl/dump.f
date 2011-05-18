
       WRITE (*,'(A/)') ' '
       WRITE (*,*) "name      : ", WRAP(DTYPE)_args(i)%name(1:LEN_TRIM(WRAP(DTYPE)_args(i)%name))
#ifdef __STRING
       WRITE (*,*) "value     : ", WRAP(DTYPE)_args(i)%variable(1:LEN_TRIM(WRAP(DTYPE)_args(i)%variable))
#else
       WRITE (*,*) "value     : ", WRAP(DTYPE)_args(i)%variable
#endif
       WRITE (*,*) "flag      : ", WRAP(DTYPE)_args(i)%flag(1:LEN_TRIM(WRAP(DTYPE)_args(i)%flag))
       WRITE (*,*) "long-flag : ", WRAP(DTYPE)_args(i)%long_flag(1:LEN_TRIM(WRAP(DTYPE)_args(i)%long_flag))
       WRITE (*,*) "ctrl-name : ", WRAP(DTYPE)_args(i)%ctrl_name(1:LEN_TRIM(WRAP(DTYPE)_args(i)%ctrl_name))
       WRITE (*,*) "default   : ", WRAP(DTYPE)_args(i)%default
#ifndef __LOGICAL
#ifndef __STRING
       WRITE (*,*) "min       : ", WRAP(DTYPE)_args(i)%min
       WRITE (*,*) "max       : ", WRAP(DTYPE)_args(i)%max
#endif
#endif
       WRITE (*,*) "help      : ", WRAP(DTYPE)_args(i)%help(1:LEN_TRIM(WRAP(DTYPE)_args(i)%help))
