#define WRAP(a) a
        SUBROUTINE WRAP(DTYPE)_add_arg(variable, name,  &
        &          flag, long_flag, ctrl_name, default, &
#if defined(__INTEGER) || defined(__LONGINT) || defined(__SINGLE) || defined(__DOUBLE)
        &          min, max,                            &
#elif defined(__STRING)
        &          min, max,                            &
#endif
#if defined(__LOGICAL) && !defined(ARRAY)
        &          vtype,                               &
#endif
#ifdef __F2003
        &          default_func, validator,             &
#endif
        &          help)
        !!! Adds a new arg definition.

          IMPLICIT NONE
          !----------------------------------------------------------------------
          !  Arguments
          !----------------------------------------------------------------------
#ifdef ARRAY

#ifdef __INTEGER
          INTEGER,                  DIMENSION(:), TARGET,   INTENT(IN   ) :: variable
          !!! Global variable to bind to.
          INTEGER,                  DIMENSION(:), OPTIONAL, INTENT(IN   ) :: default
          !!! Default value.
#elif defined(__LONGINT)
          INTEGER(ppm_kind_int64),  DIMENSION(:), TARGET,   INTENT(IN   ) :: variable
          !!! Global variable to bind to.
          INTEGER(ppm_kind_int64),  DIMENSION(:), OPTIONAL, INTENT(IN   ) :: default
          !!! Default value.
#elif defined(__SINGLE)
          REAL(ppm_kind_single),    DIMENSION(:), TARGET,   INTENT(IN   ) :: variable
          !!! Global variable to bind to.
          REAL(ppm_kind_single),    DIMENSION(:), OPTIONAL, INTENT(IN   ) :: default
          !!! Default value.
#elif defined(__DOUBLE)
          REAL(ppm_kind_double),    DIMENSION(:), TARGET,   INTENT(IN   ) :: variable
          !!! Global variable to bind to.
          REAL(ppm_kind_double),    DIMENSION(:), OPTIONAL, INTENT(IN   ) :: default
          !!! Default value.
#elif defined(__LOGICAL)
          LOGICAL,                  DIMENSION(:), TARGET,   INTENT(IN   ) :: variable
          !!! Global variable to bind to.
          LOGICAL,                  DIMENSION(:), OPTIONAL, INTENT(IN   ) :: default
          !!! Default value.
#elif defined(__STRING)
          CHARACTER(LEN=*),         DIMENSION(:), TARGET,   INTENT(IN   ) :: variable
          !!! Global variable to bind to.
          CHARACTER(LEN=*),         DIMENSION(:), OPTIONAL, INTENT(IN   ) :: default
          !!! Default value.
#elif defined(__COMPLEX)
          COMPLEX(ppm_kind_single), DIMENSION(:), TARGET,   INTENT(IN   ) :: variable
          !!! Global variable to bind to.
          COMPLEX(ppm_kind_single), DIMENSION(:), OPTIONAL, INTENT(IN   ) :: default
          !!! Default value.
#elif defined(__DCOMPLEX)
          COMPLEX(ppm_kind_double), DIMENSION(:), TARGET,   INTENT(IN   ) :: variable
          !!! Global variable to bind to.
          COMPLEX(ppm_kind_double), DIMENSION(:), OPTIONAL, INTENT(IN   ) :: default
          !!! Default value.
#endif

#else

#ifdef __INTEGER
          INTEGER,                                TARGET,   INTENT(IN   ) :: variable
          !!! Global variable to bind to.
          INTEGER,                                OPTIONAL, INTENT(IN   ) :: default
          !!! Default value.
#elif defined(__LONGINT)
          INTEGER(ppm_kind_int64),                TARGET,   INTENT(IN   ) :: variable
          !!! Global variable to bind to.
          INTEGER(ppm_kind_int64),                OPTIONAL, INTENT(IN   ) :: default
          !!! Default value.
#elif defined(__SINGLE)
          REAL(ppm_kind_single),                  TARGET,   INTENT(IN   ) :: variable
          !!! Global variable to bind to.
          REAL(ppm_kind_single),                  OPTIONAL, INTENT(IN   ) :: default
          !!! Default value.
#elif defined(__DOUBLE)
          REAL(ppm_kind_double),                  TARGET,   INTENT(IN   ) :: variable
          !!! Global variable to bind to.
          REAL(ppm_kind_double),                  OPTIONAL, INTENT(IN   ) :: default
          !!! Default value.
#elif defined(__LOGICAL)
          LOGICAL,                                TARGET,   INTENT(IN   ) :: variable
          !!! Global variable to bind to.
          LOGICAL,                                OPTIONAL, INTENT(IN   ) :: default
          !!! Default value.
#elif defined(__STRING)
          CHARACTER(LEN=*),                       TARGET,   INTENT(IN   ) :: variable
          !!! Global variable to bind to.
          CHARACTER(LEN=*),                       OPTIONAL, INTENT(IN   ) :: default
          !!! Default value.
#elif defined(__COMPLEX)
          COMPLEX(ppm_kind_single),               TARGET,   INTENT(IN   ) :: variable
          !!! Global variable to bind to.
          COMPLEX(ppm_kind_single),               OPTIONAL, INTENT(IN   ) :: default
          !!! Default value.
#elif defined(__DCOMPLEX)
          COMPLEX(ppm_kind_double),               TARGET,   INTENT(IN   ) :: variable
          !!! Global variable to bind to.
          COMPLEX(ppm_kind_double),               OPTIONAL, INTENT(IN   ) :: default
          !!! Default value.
#endif

#endif

#ifdef __INTEGER
          INTEGER,                                OPTIONAL, INTENT(IN   ) :: min
          !!! Minimum value of arg.
          INTEGER,                                OPTIONAL, INTENT(IN   ) :: max
          !!! Maximum value of arg.
#elif defined(__LONGINT)
          INTEGER(ppm_kind_int64),                OPTIONAL, INTENT(IN   ) :: min
          !!! Minimum value of arg.
          INTEGER(ppm_kind_int64),                OPTIONAL, INTENT(IN   ) :: max
          !!! Maximum value of arg.
#elif defined(__SINGLE)
          REAL(ppm_kind_single),                  OPTIONAL, INTENT(IN   ) :: min
          !!! Minimum value of arg.
          REAL(ppm_kind_single),                  OPTIONAL, INTENT(IN   ) :: max
          !!! Maximum value of arg.
#elif defined(__DOUBLE)
          REAL(ppm_kind_double),                  OPTIONAL, INTENT(IN   ) :: min
          !!! Minimum value of arg.
          REAL(ppm_kind_double),                  OPTIONAL, INTENT(IN   ) :: max
          !!! Maximum value of arg.
#elif defined(__STRING)
          INTEGER,                                OPTIONAL, INTENT(IN   ) :: min
          !!! Minimum length of string.
          INTEGER,                                OPTIONAL, INTENT(IN   ) :: max
          !!! Maximum length of string.
#endif

#if defined(__LOGICAL) && !defined(ARRAY)
          LOGICAL,                                 OPTIONAL, INTENT(IN   ) :: vtype
          !!! Type of flag. Logical flags require no value to be supplied.
          !!! Instead the behavior depends on this argument. One of:
          !!! enabling_flag :: Presence of flag sets variable to .TRUE.
          !!! disabling_flag :: Presence of flag sets variable to .FALSE.
#endif
          CHARACTER(LEN=*),                                  INTENT(IN   ) :: name
          !!! Name of the arg for use in the auto generated usage message/ctrl file.
          CHARACTER(LEN=2),                        OPTIONAL, INTENT(IN   ) :: flag
          !!! Single character flag (eg. _'-f'_).
          CHARACTER(LEN=*),                        OPTIONAL, INTENT(IN   ) :: long_flag
          !!! Long flag (eg. _'--flag'_). Has to start with _'--'_!
          CHARACTER(LEN=*),                        OPTIONAL, INTENT(IN   ) :: ctrl_name
          !!! Control file variable name.
          CHARACTER(LEN=*),                        OPTIONAL, INTENT(IN   ) :: help
          !!! Help string for the auto generated usage message/ctrl file.
#ifdef __F2003
          PROCEDURE(WRAP(DTYPE)_func),             OPTIONAL                :: default_func
          !!! Default function.
          PROCEDURE(WRAP(DTYPE)_func),             OPTIONAL                :: validator
          !!! Validator function.
#endif
          !----------------------------------------------------------------------
          !  Local variables
          !----------------------------------------------------------------------
          TYPE(WRAP(DTYPE)_arg), DIMENSION(:), POINTER :: temp
          TYPE(WRAP(DTYPE)_arg)                        :: def

          INTEGER :: len

          !----------------------------------------------------------------------
          !  Body
          !----------------------------------------------------------------------
          ! create default group
          IF (groups_i .EQ. -1) CALL arg_group("General Options")
          ! allocate initial storage
          IF (.NOT.ASSOCIATED(WRAP(DTYPE)_args)) THEN
             ALLOCATE(WRAP(DTYPE)_args(1:di))
             WRAP(DTYPE)_args_i = 0
             len=di
          ELSE
             len=SIZE(WRAP(DTYPE)_args)
          ENDIF

          ! increment counter
          WRAP(DTYPE)_args_i = WRAP(DTYPE)_args_i + 1

          ! grow storage by di if needed
          IF (WRAP(DTYPE)_args_i .GT. len) THEN
             ALLOCATE(temp(1:len+di))

             temp(1:len)=WRAP(DTYPE)_args(1:len)

             DEALLOCATE(WRAP(DTYPE)_args)

             WRAP(DTYPE)_args => temp
          ENDIF
          ! populate structure
          def%variable  => variable
          def%name      =  name
          ! default values
#ifdef ARRAY
          ALLOCATE(def%default(LBOUND(def%variable,1):UBOUND(def%variable,1)))
#endif
          def%default   = def%variable
          def%flag      = ''
          def%long_flag = ''
          def%ctrl_name = ''
          def%help      = ''
#if defined(__INTEGER)
          def%min       = -ppm_big_i
          def%max       =  ppm_big_i
#elif defined(__LONGINT) || defined(__SINGLE) || defined(__DOUBLE)
          def%min       = -HUGE(min)
          def%max       =  HUGE(max)
#endif
          ! supplied values
          IF (PRESENT(flag)) THEN
             def%flag                 =  flag
             def%flag_set             =  .TRUE.
             def%settable             =  .TRUE.
             group_has_arg(groups_i)  =  .TRUE.
          ENDIF
          IF (PRESENT(long_flag)) THEN
             def%long_flag            =  long_flag
             def%long_flag_set        =  .TRUE.
             def%settable             =  .TRUE.
             group_has_arg(groups_i)  =  .TRUE.
          ENDIF
          IF (PRESENT(ctrl_name)) THEN
             def%ctrl_name            =  ctrl_name
             def%ctrl_name_set        =  .TRUE.
             def%settable             =  .TRUE.
             group_has_ctrl(groups_i) =  .TRUE.
             len                      =  LEN_TRIM(ctrl_name)
             IF (len .GT. group_max_len(groups_i)) group_max_len(groups_i) = len
          ENDIF
          IF (PRESENT(help)) THEN
             def%help                 =  help
             def%help_set             =  .TRUE.
          ENDIF
          IF (PRESENT(default)) THEN
             def%default              =  default
             def%default_set          =  .TRUE.
          ENDIF
#ifdef __F2003
          IF (PRESENT(default_func)) THEN
             def%default_func         => default_func
             def%default_func_set     =  .TRUE.
          ENDIF
          IF (PRESENT(validator)) THEN
             def%validator            => validator
             def%validator_set        =  .TRUE.
          ENDIF
#endif
#if defined(__INTEGER) || defined(__LONGINT) || defined(__SINGLE) || defined(__DOUBLE)
          IF (PRESENT(min)) THEN
             def%min                  =  min
             def%min_set              =  .TRUE.
          ENDIF
          IF (PRESENT(max)) THEN
             def%max                  =  max
             def%max_set              =  .TRUE.
          ENDIF
#elif defined(__STRING)
      !     IF (PRESENT(min)) THEN
      !        def%min                  =  min
      !        def%min_set              =  .TRUE.
      !     ENDIF
      !     IF (PRESENT(max)) THEN
      !        def%max                  =  max
      !        def%max_set              =  .TRUE.
      !     ENDIF
#endif
#if defined(__LOGICAL) && !defined(ARRAY)
          IF (PRESENT(vtype)) THEN
             def%vtype = vtype
          ENDIF
#endif
          ! group
          def%group            = groups_i
          group_size(groups_i) = group_size(groups_i) + 1
          def%group_i          = group_size(groups_i)
          WRAP(DTYPE)_args(WRAP(DTYPE)_args_i) = def

        END SUBROUTINE WRAP(DTYPE)_add_arg
#undef DTYPE
#undef __INTEGER
#undef __LONGINT
#undef __SINGLE
#undef __DOUBLE
#undef __LOGICAL
#undef __STRING
#undef __COMPLEX
#undef __DCOMPLEX
