#define WRAP(a) a
#ifdef ARRAY
  SUBROUTINE WRAP(DTYPE)_array_add_arg(variable, name, flag, long_flag, ctrl_name, &
#else
  SUBROUTINE WRAP(DTYPE)_add_arg(variable, name, flag, long_flag, ctrl_name, &
#endif
       &                         default, default_func,                      &
#ifndef BOOL
#ifndef STRING
       &                         min, max,                                   &
#endif
#endif
#ifdef BOOL
#ifndef ARRAY
                                 type, &
#endif
#endif
                                 validator, help)
    !----------------------------------------------------------------------
    !  Arguments
    !----------------------------------------------------------------------
#ifdef STRING
#ifdef ARRAY
    CHARACTER(LEN=*), DIMENSION(:), TARGET,   INTENT(IN   ) :: variable
    CHARACTER(LEN=*), DIMENSION(:), OPTIONAL, INTENT(IN   ) :: default    
#else
    CHARACTER(LEN=*), TARGET,              INTENT(IN   ) :: variable
    CHARACTER(LEN=*),            OPTIONAL, INTENT(IN   ) :: default
#endif
#else
#ifdef ARRAY
    DTYPE, DIMENSION(:), TARGET,           INTENT(IN   ) :: variable
    DTYPE, DIMENSION(:),         OPTIONAL, INTENT(IN   ) :: default    
#else
    DTYPE, TARGET,                         INTENT(IN   ) :: variable
    DTYPE,                       OPTIONAL, INTENT(IN   ) :: default
#endif
#ifndef BOOL
    DTYPE,                       OPTIONAL, INTENT(IN   ) :: min
    DTYPE,                       OPTIONAL, INTENT(IN   ) :: max
#endif
#endif
#ifdef BOOL
#ifndef ARRAY
    LOGICAL,                    OPTIONAL, INTENT(IN   ) :: type
#endif
#endif
    CHARACTER(LEN=*),                      INTENT(IN   ) :: name
    CHARACTER(LEN=2),            OPTIONAL, INTENT(IN   ) :: flag
    CHARACTER(LEN=*),            OPTIONAL, INTENT(IN   ) :: long_flag
    CHARACTER(LEN=*),            OPTIONAL, INTENT(IN   ) :: ctrl_name
    CHARACTER(LEN=*),            OPTIONAL, INTENT(IN   ) :: help
#ifdef ARRAY
    PROCEDURE(WRAP(DTYPE)_array_dflt), OPTIONAL         :: default_func
    PROCEDURE(WRAP(DTYPE)_array_vdtr), OPTIONAL         :: validator
#else
    PROCEDURE(WRAP(DTYPE)_dflt), OPTIONAL               :: default_func
    PROCEDURE(WRAP(DTYPE)_vdtr), OPTIONAL               :: validator
#endif
    !----------------------------------------------------------------------
    !  Local variables
    !----------------------------------------------------------------------
#ifdef ARRAY
    TYPE(WRAP(DTYPE)_array_arg), DIMENSION(:), POINTER      :: temp
    TYPE(WRAP(DTYPE)_array_arg)                             :: def
#else
    TYPE(WRAP(DTYPE)_arg), DIMENSION(:), POINTER            :: temp
    TYPE(WRAP(DTYPE)_arg)                                   :: def
#endif
    !----------------------------------------------------------------------
    !  Body
    !----------------------------------------------------------------------
    ! create default group
    IF (groups_i .EQ. -1) CALL arg_group("General Options")
#ifdef ARRAY
    ! allocate initial storage
    IF (.NOT. ASSOCIATED(WRAP(DTYPE)_array_args)) THEN
       ALLOCATE(WRAP(DTYPE)_array_args(1:di))
       WRAP(DTYPE)_array_args_i = 0
    END IF
    ! increment counter
    WRAP(DTYPE)_array_args_i = WRAP(DTYPE)_array_args_i + 1
    ! grow storage by di if needed
    IF (WRAP(DTYPE)_array_args_i .GT. SIZE(WRAP(DTYPE)_array_args)) THEN
       ALLOCATE(temp(1:SIZE(WRAP(DTYPE)_array_args)+di))
       temp(1:SIZE(WRAP(DTYPE)_array_args)) = WRAP(DTYPE)_array_args
       DEALLOCATE(WRAP(DTYPE)_array_args)
       WRAP(DTYPE)_array_args => temp
    END IF
#else
    ! allocate initial storage
    IF (.NOT. ASSOCIATED(WRAP(DTYPE)_args)) THEN
       ALLOCATE(WRAP(DTYPE)_args(1:di))
       WRAP(DTYPE)_args_i = 0
    END IF
    ! increment counter
    WRAP(DTYPE)_args_i = WRAP(DTYPE)_args_i + 1
    ! grow storage by di if needed
    IF (WRAP(DTYPE)_args_i .GT. SIZE(WRAP(DTYPE)_args)) THEN
       ALLOCATE(temp(1:SIZE(WRAP(DTYPE)_args)+di))
       temp(1:SIZE(WRAP(DTYPE)_args)) = WRAP(DTYPE)_args
       DEALLOCATE(WRAP(DTYPE)_args)
       WRAP(DTYPE)_args => temp
    END IF
#endif
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
#ifndef BOOL
#ifndef STRING
    def%min       = -HUGE(min)
    def%max       =  HUGE(max)
#endif
#endif
    ! supplied values
    IF (PRESENT(flag)) THEN
       def%flag                 =  flag
       def%flag_set             =  .TRUE.
       group_has_arg(groups_i)  = .TRUE.
    END IF
    IF (PRESENT(long_flag)) THEN
       def%long_flag            =  long_flag
       def%long_flag_set        =  .TRUE.
       group_has_arg(groups_i)  = .TRUE.
    END IF
    IF (PRESENT(ctrl_name)) THEN
       def%ctrl_name            =  ctrl_name
       def%ctrl_name_set        = .TRUE.
       group_has_ctrl(groups_i) = .TRUE.
    END IF
    IF (PRESENT(help)) THEN
       def%help                 =  help
       def%help_set             =  .TRUE.
    END IF
    IF (PRESENT(default)) THEN
       def%default              =  default
       def%default_set          =  .TRUE.
    END IF
    IF (PRESENT(default_func)) def%default_func => default_func
    IF (PRESENT(validator))    def%validator    => validator
#ifndef BOOL
#ifndef STRING
    IF (PRESENT(min)) THEN
       def%min                  =  min
       def%min_set              =  .TRUE.
    END IF
    IF (PRESENT(max)) THEN
       def%max                  =  max
       def%max_set              =  .TRUE.
    END IF
#endif
#endif
#ifdef BOOL
#ifndef ARRAY
    IF (PRESENT(type)) THEN
       def%type         =  type
    END IF
#endif
#endif
    ! group
    def%group            = groups_i
    group_size(groups_i) = group_size(groups_i) + 1
    def%group_i          = group_size(groups_i)
#ifdef ARRAY
    WRAP(DTYPE)_array_args(WRAP(DTYPE)_array_args_i) = def
  END SUBROUTINE WRAP(DTYPE)_array_add_arg
#else
    WRAP(DTYPE)_args(WRAP(DTYPE)_args_i) = def
  END SUBROUTINE WRAP(DTYPE)_add_arg
#endif
