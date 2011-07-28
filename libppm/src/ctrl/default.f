#define WRAP(a) a
    DO i=1,WRAP(DTYPE)_args_i
       IF (WRAP(DTYPE)_args(i)%default_set) &
            WRAP(DTYPE)_args(i)%variable = WRAP(DTYPE)_args(i)%default
    END DO
#undef DTYPE

