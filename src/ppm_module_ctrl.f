      !----*- f90 -*-------------------------------------------------------------
      !  Module       :                    ppm_module_ctrl
      !--------------------------------------------------------------------------
      !
      !  Purpose      : Easy definition and parsing of command line and
      !                 control file arguments.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !--------------------------------------------------------------------------
      !
      !--------------------------------------------------------------------------
      !  Milan Mitrovic
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !--------------------------------------------------------------------------
      MODULE ppm_module_ctrl
      !!! This module provides an easy way to handle command line
      !!! arguments and control files.
      !!!
      !!! .Intended Usage
      !!!
      !!! This module assumes that you will create a module that contains
      !!! all your global variables (which we will call +client_global+) and
      !!! gives you an easy way of supplying values for these variables,
      !!! either through command line arguments or control files.
      !!!
      !!! To do this, you have to create a subroutine (which we will call
      !!! +define_args+) that will hold the initialization code. See the
      !!! following example.
      !!!
      !!! [source,fortran]
      !!! ----
      !!! MODULE client_global
      !!!   USE ppm_module_ctrl
      !!!   INTEGER :: example
      !!! CONTAINS
      !!!   SUBROUTINE define_args
      !!!      CALL arg(example, 'example',    &
      !!!               flag      = '-e',      &
      !!!               ctrl_name = 'example')
      !!!   END SUBROUTINE define_args
      !!! END MODULE client_global
      !!! ----
      !!!
      !!! This will create an argument of type +integer+ that can be set
      !!! either through the command line flag _'-e'_ or the control file
      !!! variable _'example'_. Please take note that the command line flag
      !!! will override the setting in the control file if both are
      !!! supplied.
      !!!
      !!! For this code to work you need to +call define_args+ somewhere in
      !!! your initialization code, followed by a +call parse_args(info)+.
      !!! After +parse_args+ completes your globals will be initialized
      !!! with the supplied values.
      !!!
      !!! The function +arg+ is overloaded to support +integer+, +real+,
      !!! +character+, +logical+ and +complex+ arguments or fixed sized
      !!! arrays of these types. The first two arguments are required and
      !!! they are: the variable itself (to which a pointer will be
      !!! stored) and the name of the variable (used for printing help
      !!! messages).
      !!!
      !!! All other arguments are optional and the details of how they work
      !!! can be found below. Here is just a short overview and reference
      !!! of the available options:
      !!! [horizontal]
      !!! flag :: [_'-f'_] single character command line flag
      !!! long_flag :: [_'--long-flag'_] long command line flag (starts with
      !!! --)
      !!! ctrl_name :: [_'name'_] name of the control file varible
      !!! default :: [_42_] default value
      !!! min :: [_4_] minimum value for numeric types
      !!! max :: [_30_] maximum value
      !!! help :: ['Description.'] help message to display in the auto
      !!! generated ctrl file and command line help. You can use +\n+ to
      !!! force a line break.
      !!! default_func :: [_external_func_] custom function to compute the
      !!! value of the variable after other globals have been set (only
      !!! available when compiled with F2003 support)
      !!! validator :: [_external_func_] custom function to validate the
      !!! variable value (only available when compiled with F2003 support)
      !!!
      !!! By default the module supports _-h_ and _--help_ flags for
      !!! printing the help message, and _--print-ctrl_ for printing a
      !!! sample control file. There is also an optional first positional
      !!! argument that is interpreted as the name of the control file.
      !!!
      !!! To make the output prettier you can add calls to
      !!! +arg_group(_'name'_)+ to your +define_args+ and all calls to
      !!! +arg+ following a group definition will be put into that group.
      !!!
      !!! If any one of the printing flags is present +parse_args+ will return
      !!! _exit_gracefully_ which you should check for and exit gracefully.
      !!!

        !------------------------------------------------------------------------
        !  Modules
        !------------------------------------------------------------------------
        USE ppm_module_data
        USE ppm_module_substart
        USE ppm_module_substop
        USE ppm_module_write
        USE ppm_module_error
        IMPLICIT NONE

        !------------------------------------------------------------------------
        !  Interface
        !------------------------------------------------------------------------
        PUBLIC :: arg, arg_group, parse_args, disable_help, disable_ctrl, &
             &    set_ctrl_name,                                          &
#ifdef __F2003
             &    integer_func, longint_func, single_func, double_func,   &
             &    logical_func, string_func, complex_func, dcomplex_func, &
             &    integer_array_func, longint_array_func,                 &
             &    single_array_func, double_array_func,                   &
             &    logical_array_func, string_array_func,                  &
             &    complex_array_func, dcomplex_array_func,                &
#endif
             &    reset, add_cmd, ctrl_file_name, break_help,             &
             &    find_arg, find_flag, arg_count,                         &
             &    enabling_flag, disabling_flag, exit_gracefully

        PRIVATE

#ifdef __MPI
        INCLUDE 'mpif.h'
#endif
        !------------------------------------------------------------------------
        !  Types
        !------------------------------------------------------------------------
        ! scalar

#define DTYPE INTEGER
#define __INTEGER
#include "ctrl/type.f"

#define DTYPE LONGINT
#define __LONGINT
#include "ctrl/type.f"

#define DTYPE SINGLE
#define __SINGLE
#include "ctrl/type.f"

#define DTYPE DOUBLE
#define __DOUBLE
#include "ctrl/type.f"

#define DTYPE LOGICAL
#define __LOGICAL
#include "ctrl/type.f"

#define DTYPE STRING
#define __STRING
#include "ctrl/type.f"

#define DTYPE COMPLEX
#define __COMPLEX
#include "ctrl/type.f"

#define DTYPE DCOMPLEX
#define __DCOMPLEX
#include "ctrl/type.f"

        ! array
#define ARRAY

#define DTYPE INTEGER_array
#define __INTEGER
#include "ctrl/type.f"

#define DTYPE LONGINT_array
#define __LONGINT
#include "ctrl/type.f"

#define DTYPE SINGLE_array
#define __SINGLE
#include "ctrl/type.f"

#define DTYPE DOUBLE_array
#define __DOUBLE
#include "ctrl/type.f"

#define DTYPE LOGICAL_array
#define __LOGICAL
#include "ctrl/type.f"

#define DTYPE STRING_array
#define __STRING
#include "ctrl/type.f"

#define DTYPE COMPLEX_array
#define __COMPLEX
#include "ctrl/type.f"

#define DTYPE DCOMPLEX_array
#define __DCOMPLEX
#include "ctrl/type.f"

#undef ARRAY

        !------------------------------------------------------------------------
        !  Interfaces
        !------------------------------------------------------------------------
        INTERFACE arg
           ! scalar
           MODULE PROCEDURE INTEGER_add_arg
           MODULE PROCEDURE LONGINT_add_arg
           MODULE PROCEDURE SINGLE_add_arg
           MODULE PROCEDURE DOUBLE_add_arg
           MODULE PROCEDURE LOGICAL_add_arg
           MODULE PROCEDURE STRING_add_arg
           MODULE PROCEDURE COMPLEX_add_arg
           MODULE PROCEDURE DCOMPLEX_add_arg
           ! array
           MODULE PROCEDURE INTEGER_array_add_arg
           MODULE PROCEDURE LONGINT_array_add_arg
           MODULE PROCEDURE SINGLE_array_add_arg
           MODULE PROCEDURE DOUBLE_array_add_arg
           MODULE PROCEDURE LOGICAL_array_add_arg
           MODULE PROCEDURE STRING_array_add_arg
           MODULE PROCEDURE COMPLEX_array_add_arg
           MODULE PROCEDURE DCOMPLEX_array_add_arg
        END INTERFACE

#ifdef __F2003
        ABSTRACT INTERFACE
           !---------------------------------------------------------------------
           !  Defaults and Validators
           !---------------------------------------------------------------------
           ! scalar
           LOGICAL FUNCTION INTEGER_func(variable)
             INTEGER, POINTER :: variable
           END FUNCTION INTEGER_func

           LOGICAL FUNCTION LONGINT_func(variable)
             INTEGER(8), POINTER :: variable
           END FUNCTION LONGINT_func

           LOGICAL FUNCTION SINGLE_func(variable)
             REAL(KIND(1.0E0)), POINTER :: variable
           END FUNCTION SINGLE_func

           LOGICAL FUNCTION DOUBLE_func(variable)
             REAL(KIND(1.0D0)), POINTER :: variable
           END FUNCTION DOUBLE_func

           LOGICAL FUNCTION LOGICAL_func(variable)
             LOGICAL, POINTER :: variable
           END FUNCTION LOGICAL_func

           LOGICAL FUNCTION STRING_func(variable)
             CHARACTER(LEN=*), POINTER :: variable
           END FUNCTION STRING_func

           LOGICAL FUNCTION COMPLEX_func(variable)
             COMPLEX(KIND(1.0E0)), POINTER :: variable
           END FUNCTION COMPLEX_func

           LOGICAL FUNCTION DCOMPLEX_func(variable)
             COMPLEX(KIND(1.0D0)), POINTER :: variable
           END FUNCTION DCOMPLEX_func

           ! array
           LOGICAL FUNCTION INTEGER_array_func(variable)
             INTEGER, DIMENSION(:), POINTER :: variable
           END FUNCTION INTEGER_array_func

           LOGICAL FUNCTION LONGINT_array_func(variable)
             INTEGER(8), DIMENSION(:), POINTER :: variable
           END FUNCTION LONGINT_array_func

           LOGICAL FUNCTION SINGLE_array_func(variable)
             REAL(KIND(1.0E0)), DIMENSION(:), POINTER :: variable
           END FUNCTION SINGLE_array_func

           LOGICAL FUNCTION DOUBLE_array_func(variable)
             REAL(KIND(1.0D0)), DIMENSION(:), POINTER :: variable
           END FUNCTION DOUBLE_array_func

           LOGICAL FUNCTION LOGICAL_array_func(variable)
             LOGICAL, DIMENSION(:), POINTER :: variable
           END FUNCTION LOGICAL_array_func

           LOGICAL FUNCTION STRING_array_func(variable)
             CHARACTER(LEN=*), DIMENSION(:), POINTER :: variable
           END FUNCTION STRING_array_func

           LOGICAL FUNCTION COMPLEX_array_func(variable)
             COMPLEX(KIND(1.0E0)), DIMENSION(:), POINTER :: variable
           END FUNCTION COMPLEX_array_func

           LOGICAL FUNCTION DCOMPLEX_array_func(variable)
             COMPLEX(KIND(1.0D0)), DIMENSION(:), POINTER :: variable
           END FUNCTION DCOMPLEX_array_func

        END INTERFACE
#endif
        !------------------------------------------------------------------------
        !  Constants
        !------------------------------------------------------------------------
        LOGICAL, PARAMETER :: enabling_flag   = .TRUE.
        !!! Value for type option of logical args. Presence of flag sets
        !!! variable to +.TRUE.+
        LOGICAL, PARAMETER :: disabling_flag  = .FALSE.
        !!! Value for type option of logical args. Presence of flag sets
        !!! variable to +.FALSE.+
        INTEGER, PARAMETER :: exit_gracefully = 42
        !!! Value of the +info+ argument of +parse_args+ that signals that
        !!! the program should exit without error.
        !------------------------------------------------------------------------
        !  Variables
        !------------------------------------------------------------------------
        ! by how much to grow storage
        INTEGER,                           PARAMETER    :: di = 10
        ! scalar
        TYPE(INTEGER_arg),        POINTER, DIMENSION(:) :: INTEGER_args    => NULL()
        INTEGER                                         :: INTEGER_args_i  = 0
        TYPE(LONGINT_arg),        POINTER, DIMENSION(:) :: LONGINT_args    => NULL()
        INTEGER                                         :: LONGINT_args_i  = 0
        TYPE(SINGLE_arg),         POINTER, DIMENSION(:) :: SINGLE_args     => NULL()
        INTEGER                                         :: SINGLE_args_i   = 0
        TYPE(DOUBLE_arg),         POINTER, DIMENSION(:) :: DOUBLE_args     => NULL()
        INTEGER                                         :: DOUBLE_args_i   = 0
        TYPE(LOGICAL_arg),        POINTER, DIMENSION(:) :: LOGICAL_args    => NULL()
        INTEGER                                         :: LOGICAL_args_i  = 0
        TYPE(STRING_arg),         POINTER, DIMENSION(:) :: STRING_args     => NULL()
        INTEGER                                         :: STRING_args_i   = 0
        TYPE(COMPLEX_arg),        POINTER, DIMENSION(:) :: COMPLEX_args    => NULL()
        INTEGER                                         :: COMPLEX_args_i  = 0
        TYPE(DCOMPLEX_arg),       POINTER, DIMENSION(:) :: DCOMPLEX_args   => NULL()
        INTEGER                                         :: DCOMPLEX_args_i = 0
        ! add arrays
        TYPE(INTEGER_array_arg),  POINTER, DIMENSION(:) :: INTEGER_array_args    => NULL()
        INTEGER                                         :: INTEGER_array_args_i  = 0
        TYPE(LONGINT_array_arg),  POINTER, DIMENSION(:) :: LONGINT_array_args    => NULL()
        INTEGER                                         :: LONGINT_array_args_i  = 0
        TYPE(SINGLE_array_arg),   POINTER, DIMENSION(:) :: SINGLE_array_args     => NULL()
        INTEGER                                         :: SINGLE_array_args_i   = 0
        TYPE(DOUBLE_array_arg),   POINTER, DIMENSION(:) :: DOUBLE_array_args     => NULL()
        INTEGER                                         :: DOUBLE_array_args_i   = 0
        TYPE(LOGICAL_array_arg),  POINTER, DIMENSION(:) :: LOGICAL_array_args    => NULL()
        INTEGER                                         :: LOGICAL_array_args_i  = 0
        TYPE(STRING_array_arg),   POINTER, DIMENSION(:) :: STRING_array_args     => NULL()
        INTEGER                                         :: STRING_array_args_i   = 0
        TYPE(COMPLEX_array_arg),  POINTER, DIMENSION(:) :: COMPLEX_array_args    => NULL()
        INTEGER                                         :: COMPLEX_array_args_i  = 0
        TYPE(DCOMPLEX_array_arg), POINTER, DIMENSION(:) :: DCOMPLEX_array_args   => NULL()
        INTEGER                                         :: DCOMPLEX_array_args_i = 0
        ! arg storage
        CHARACTER(LEN=ppm_char),  POINTER, DIMENSION(:) :: cmd_args       => NULL()
        INTEGER,                  POINTER, DIMENSION(:) :: cmd_args_len   => NULL()
        LOGICAL,                  POINTER, DIMENSION(:) :: cmd_args_used  => NULL()
        INTEGER                                         :: cmd_args_i     =  0
        INTEGER                                         :: cmd_i          =  0
        ! arg groups
        CHARACTER(LEN=ppm_char),  POINTER, DIMENSION(:) :: groups         => NULL()
        INTEGER,                  POINTER, DIMENSION(:) :: group_size     => NULL()
        INTEGER,                  POINTER, DIMENSION(:) :: group_max_len  => NULL()
        LOGICAL,                  POINTER, DIMENSION(:) :: group_has_ctrl => NULL()
        LOGICAL,                  POINTER, DIMENSION(:) :: group_has_arg  => NULL()
        INTEGER                                         :: groups_i       =  -1
        ! special args
        LOGICAL                                         :: help_enabled   = .TRUE.
        LOGICAL                                         :: ctrl_enabled   = .TRUE.
        CHARACTER(LEN=ppm_char)                         :: ctrl_file_name = 'Ctrl'
        ! test run
        LOGICAL                                         :: in_test        = .FALSE.

      CONTAINS
        !------------------------------------------------------------------------
        !  Master procedure
        !------------------------------------------------------------------------
        SUBROUTINE parse_args(info)

          IMPLICIT NONE

          !----------------------------------------------------------------------
          !  Arguments
          !----------------------------------------------------------------------
          INTEGER, INTENT(  OUT) :: info
          !----------------------------------------------------------------------
          !  Local variables
          !----------------------------------------------------------------------
          REAL(8) :: t0

          INTEGER :: info2
          INTEGER :: i
          INTEGER :: rank = 0

          CHARACTER(LEN=*), PARAMETER :: caller='parse_args'
          CHARACTER(LEN=ppm_char)     :: value

          LOGICAL :: ok
          LOGICAL :: printing_ctrl = .FALSE.

          !----------------------------------------------------------------------
          !  Externals
          !----------------------------------------------------------------------
          EXTERNAL iargc
          INTEGER  iargc
          !----------------------------------------------------------------------
          !  Initialize
          !----------------------------------------------------------------------
          CALL substart(caller, t0, info)
          !----------------------------------------------------------------------
          !  Do everything on rank 0 and bcast at the end
          !----------------------------------------------------------------------
#ifdef __MPI
          CALL MPI_Comm_Rank(MPI_COMM_WORLD, rank, info)
#endif
          IF (rank .EQ. 0) THEN
             !-------------------------------------------------------------------
             !  Copy default values into variables
             !-------------------------------------------------------------------
             CALL apply_defaults(info)
             or_fail('Applying defaults failed!',ppm_err_argument,exit_point=100,ppm_error=ppm_error_fatal)

             !-------------------------------------------------------------------
             !  Read in the command line
             !-------------------------------------------------------------------
             IF (.NOT. in_test) THEN
                CALL read_cmd_args(info)
                or_fail_alloc('Reading command line args failed!',exit_point=100,ppm_error=ppm_error_fatal)
             END IF
             !-------------------------------------------------------------------
             !  Parse help flag
             !-------------------------------------------------------------------
             IF (help_enabled) THEN
                CALL find_flag('-h', ok)
                IF (.NOT. ok) CALL find_flag('--help', ok)
                IF (ok) THEN
                   CALL print_help
                   info = exit_gracefully
                   GOTO 100
                END IF
             END IF
             !-------------------------------------------------------------------
             !  Print Control file
             !-------------------------------------------------------------------
             IF (ctrl_enabled) THEN
                CALL find_flag('--print-ctrl', ok)
                IF (ok) THEN
                   printing_ctrl = .TRUE.
                END IF
             END IF
             !-------------------------------------------------------------------
             !  Parse rest of the command line
             !-------------------------------------------------------------------
             CALL parse_cmd_line(info)
             or_fail('Parsing command line args failed!',ppm_err_argument,exit_point=100,ppm_error=ppm_error_fatal)

             !-------------------------------------------------------------------
             !  Parse Control file
             !-------------------------------------------------------------------
             IF (ctrl_enabled) THEN
                CALL find_arg(1, ok, ctrl_file_name)
                CALL parse_ctrl_file(info)
                or_fail('Parsing control file failed!',ppm_err_argument,exit_point=100,ppm_error=ppm_error_fatal)
             END IF
#ifdef __F2003
             !-------------------------------------------------------------------
             !  Call default funcs
             !-------------------------------------------------------------------
             CALL call_default_funcs(info)
             or_fail('Calling default functions failed!',ppm_err_argument,exit_point=100,ppm_error=ppm_error_fatal)
#endif
             !-------------------------------------------------------------------
             !  Check minmax
             !-------------------------------------------------------------------
             CALL check_minmax(info)
             or_fail('Min/max check failed!',ppm_err_argument,exit_point=100,ppm_error=ppm_error_fatal)
#ifdef __F2003
             !-------------------------------------------------------------------
             !  Run validators
             !-------------------------------------------------------------------
             CALL call_validator_funcs(info)
             or_fail('Calling validator functions failed!',ppm_err_argument,exit_point=100,ppm_error=ppm_error_fatal)
#endif
             !-------------------------------------------------------------------
             !  Print Control file
             !-------------------------------------------------------------------
             IF (printing_ctrl) THEN
                CALL print_ctrl
                info = exit_gracefully
                GOTO 100
             END IF
             !-------------------------------------------------------------------
             !  DONE!
             !-------------------------------------------------------------------
          END IF ! (ppm_rank .EQ. 0)
          !----------------------------------------------------------------------
          !  Exchange data
          !----------------------------------------------------------------------
      100 CONTINUE
#ifdef __MPI
      !   CALL MPI_BCast(what, length, MPI_TYPE, 0, comm, info)
          CALL MPI_BCast(info, 1, MPI_INTEGER, 0, ppm_comm, info2)
          IF (info .NE. 0) GOTO 9999

          ! scalar
          DO i=1,INTEGER_args_i
             CALL MPI_BCast(INTEGER_args(i)%variable, 1, MPI_INTEGER, 0, ppm_comm, info)
          END DO
          DO i=1,LONGINT_args_i
             CALL MPI_BCast(LONGINT_args(i)%variable, 1, MPI_INTEGER8, 0, ppm_comm, info)
          END DO
          DO i=1,SINGLE_args_i
             CALL MPI_BCast(SINGLE_args(i)%variable, 1, MPI_REAL, 0, ppm_comm, info)
          END DO
          DO i=1,DOUBLE_args_i
             CALL MPI_BCast(DOUBLE_args(i)%variable, 1, MPI_DOUBLE_PRECISION, 0, ppm_comm, info)
          END DO
          DO i=1,STRING_args_i
             CALL MPI_BCast(STRING_args(i)%variable, ppm_char, MPI_CHARACTER, 0, ppm_comm, info)
          END DO
          DO i=1,LOGICAL_args_i
             CALL MPI_BCast(LOGICAL_args(i)%variable, 1, MPI_LOGICAL, 0, ppm_comm, info)
          END DO
          DO i=1,COMPLEX_args_i
             CALL MPI_BCast(COMPLEX_args(i)%variable, 1, MPI_COMPLEX, 0, ppm_comm, info)
          END DO
          DO i=1,DCOMPLEX_args_i
             CALL MPI_BCast(DCOMPLEX_args(i)%variable, 1, MPI_DOUBLE_COMPLEX, 0, ppm_comm, info)
          END DO

          ! array
          DO i=1,INTEGER_array_args_i
             CALL MPI_BCast(INTEGER_array_args(i)%variable, &
             &    SIZE(INTEGER_array_args(i)%variable),     &
             &    MPI_INTEGER, 0, ppm_comm, info)
          END DO
          DO i=1,LONGINT_array_args_i
             CALL MPI_BCast(LONGINT_array_args(i)%variable, &
             &    SIZE(LONGINT_array_args(i)%variable),     &
             &    MPI_INTEGER8, 0, ppm_comm, info)
          END DO
          DO i=1,SINGLE_array_args_i
             CALL MPI_BCast(SINGLE_array_args(i)%variable, &
             &    SIZE(SINGLE_array_args(i)%variable),     &
             &    MPI_REAL, 0, ppm_comm, info)
          END DO
          DO i=1,DOUBLE_array_args_i
             CALL MPI_BCast(DOUBLE_array_args(i)%variable, &
             &    SIZE(DOUBLE_array_args(i)%variable),     &
             &    MPI_DOUBLE_PRECISION, 0, ppm_comm, info)
          END DO
          DO i=1,STRING_array_args_i
             CALL MPI_BCast(STRING_array_args(i)%variable,        &
             &    SIZE(STRING_array_args(i)%variable) * ppm_char, &
             &    MPI_CHARACTER, 0, ppm_comm, info)
          END DO
          DO i=1,LOGICAL_array_args_i
             CALL MPI_BCast(LOGICAL_array_args(i)%variable, &
             &    SIZE(LOGICAL_array_args(i)%variable),     &
             &    MPI_LOGICAL, 0, ppm_comm, info)
          END DO
          DO i=1,COMPLEX_array_args_i
             CALL MPI_BCast(COMPLEX_array_args(i)%variable, &
             &    SIZE(COMPLEX_array_args(i)%variable),     &
             &    MPI_COMPLEX, 0, ppm_comm, info)
          END DO
          DO i=1,DCOMPLEX_array_args_i
             CALL MPI_BCast(DCOMPLEX_array_args(i)%variable, &
             &    SIZE(DCOMPLEX_array_args(i)%variable),     &
             &    MPI_DOUBLE_COMPLEX, 0, ppm_comm, info)
          END DO

#endif
          !----------------------------------------------------------------------
          !  Error handling
          !----------------------------------------------------------------------
      9999 CONTINUE
          ! cleanup
          CALL deallocate_memory(info .NE. 0)
          CALL substop(caller, t0, info)
          RETURN
        END SUBROUTINE parse_args
        !------------------------------------------------------------------------
        !  Cleanup
        !------------------------------------------------------------------------
        SUBROUTINE deallocate_memory(all)

          IMPLICIT NONE

          LOGICAL, INTENT(IN   ) :: all

          IF (all) THEN
             IF (ASSOCIATED(cmd_args))         DEALLOCATE(cmd_args)
             IF (ASSOCIATED(cmd_args_len))     DEALLOCATE(cmd_args_len)
             IF (ASSOCIATED(cmd_args_used))    DEALLOCATE(cmd_args_used)
          END IF
          ! scalar
          IF (ASSOCIATED(INTEGER_args))        DEALLOCATE(INTEGER_args)
          IF (ASSOCIATED(LONGINT_args))        DEALLOCATE(LONGINT_args)
          IF (ASSOCIATED(SINGLE_args))         DEALLOCATE(SINGLE_args)
          IF (ASSOCIATED(DOUBLE_args))         DEALLOCATE(DOUBLE_args)
          IF (ASSOCIATED(LOGICAL_args))        DEALLOCATE(LOGICAL_args)
          IF (ASSOCIATED(STRING_args))         DEALLOCATE(STRING_args)
          IF (ASSOCIATED(COMPLEX_args))        DEALLOCATE(COMPLEX_args)
          IF (ASSOCIATED(DCOMPLEX_args))       DEALLOCATE(DCOMPLEX_args)
          ! array
          IF (ASSOCIATED(INTEGER_array_args))  DEALLOCATE(INTEGER_array_args)
          IF (ASSOCIATED(LONGINT_array_args))  DEALLOCATE(LONGINT_array_args)
          IF (ASSOCIATED(SINGLE_array_args))   DEALLOCATE(SINGLE_array_args)
          IF (ASSOCIATED(DOUBLE_array_args))   DEALLOCATE(DOUBLE_array_args)
          IF (ASSOCIATED(LOGICAL_array_args))  DEALLOCATE(LOGICAL_array_args)
          IF (ASSOCIATED(STRING_array_args))   DEALLOCATE(STRING_array_args)
          IF (ASSOCIATED(COMPLEX_array_args))  DEALLOCATE(COMPLEX_array_args)
          IF (ASSOCIATED(DCOMPLEX_array_args)) DEALLOCATE(DCOMPLEX_array_args)
          ! other
          IF (ASSOCIATED(groups))              DEALLOCATE(groups)
          IF (ASSOCIATED(group_size))          DEALLOCATE(group_size)
          IF (ASSOCIATED(group_has_ctrl))      DEALLOCATE(group_has_ctrl)
          IF (ASSOCIATED(group_has_arg))       DEALLOCATE(group_has_arg)

        END SUBROUTINE deallocate_memory

        SUBROUTINE reset
        !!! Debugging and testing routine. Resets all module variables.

          IMPLICIT NONE

          CALL deallocate_memory(.TRUE.)
          ! scalar
          INTEGER_args_i        = 0
          LONGINT_args_i        = 0
          SINGLE_args_i         = 0
          DOUBLE_args_i         = 0
          LOGICAL_args_i        = 0
          STRING_args_i         = 0
          COMPLEX_args_i        = 0
          DCOMPLEX_args_i       = 0
          ! array
          INTEGER_array_args_i  = 0
          LONGINT_array_args_i  = 0
          SINGLE_array_args_i   = 0
          DOUBLE_array_args_i   = 0
          LOGICAL_array_args_i  = 0
          STRING_array_args_i   = 0
          COMPLEX_array_args_i  = 0
          DCOMPLEX_array_args_i = 0
          ! other
          cmd_args_i            = 0
          cmd_i                 = 0
          groups_i              = -1
          help_enabled          = .TRUE.
          ctrl_enabled          = .TRUE.
          ctrl_file_name        = 'Ctrl'
          in_test               = .FALSE.

        END SUBROUTINE reset

        !-------------------------------------------------------------------------
        !  Apply defaults
        !-------------------------------------------------------------------------
        SUBROUTINE apply_defaults(info)

          IMPLICIT NONE

          INTEGER, INTENT(  OUT) :: info

          REAL(ppm_kind_double) :: t0

          INTEGER :: i

          CHARACTER(LEN=*), PARAMETER :: caller='apply_defaults'

          !-------------------------------------------------------------------------
          !  Initialise
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)
          ! scalar
#define DTYPE INTEGER
#include "ctrl/default.f"
#define DTYPE LONGINT
#include "ctrl/default.f"
#define DTYPE SINGLE
#include "ctrl/default.f"
#define DTYPE DOUBLE
#include "ctrl/default.f"
#define DTYPE LOGICAL
#include "ctrl/default.f"
#define DTYPE STRING
#include "ctrl/default.f"
#define DTYPE COMPLEX
#include "ctrl/default.f"
#define DTYPE DCOMPLEX
#include "ctrl/default.f"
        ! array
#define DTYPE INTEGER_array
#include "ctrl/default.f"
#define DTYPE LONGINT_array
#include "ctrl/default.f"
#define DTYPE SINGLE_array
#include "ctrl/default.f"
#define DTYPE DOUBLE_array
#include "ctrl/default.f"
#define DTYPE LOGICAL_array
#include "ctrl/default.f"
#define DTYPE STRING_array
#include "ctrl/default.f"
#define DTYPE COMPLEX_array
#include "ctrl/default.f"
#define DTYPE DCOMPLEX_array
#include "ctrl/default.f"
            !-------------------------------------------------------------------------
            !  Return
            !-------------------------------------------------------------------------
      9999 CONTINUE
           CALL substop(caller,t0,info)
        END SUBROUTINE apply_defaults
        !------------------------------------------------------------------------
        !  Read in command line args
        !------------------------------------------------------------------------
        SUBROUTINE read_cmd_args(info)

          IMPLICIT NONE

          INTEGER, INTENT(  OUT) :: info

          INTEGER :: i, start, nargc

          CHARACTER(LEN=ppm_char) :: cbuf
#ifdef __MPI
#ifdef __SUPER_UX
          nargc = iargc()
          start = 0
#elif defined __SunOS
          nargc = iargc() - 6
          start = 0
#elif defined __HP_UX
          nargc = iargc() + 1
          start = 1
#elif defined __IRIX64
          nargc = iargc() - 4
          start = 1
#elif defined __Linux
          nargc = COMMAND_ARGUMENT_COUNT()
          start = 0
#elif defined __MacOS
          nargc = iargc() - 4
          start = 0
#else
          nargc = COMMAND_ARGUMENT_COUNT()
          start = 0
#endif
#else
          nargc = COMMAND_ARGUMENT_COUNT()
          start = 0
#endif
          cmd_args_i = nargc-start
          ! allocate storage
          IF (ASSOCIATED(cmd_args))      DEALLOCATE(cmd_args)
          IF (ASSOCIATED(cmd_args_len))  DEALLOCATE(cmd_args_len)
          IF (ASSOCIATED(cmd_args_used)) DEALLOCATE(cmd_args_used)

          ALLOCATE(cmd_args(1:cmd_args_i),      STAT=info)
          ALLOCATE(cmd_args_len(1:cmd_args_i),  STAT=info)
          ALLOCATE(cmd_args_used(1:cmd_args_i), STAT=info)
          IF (info .NE. 0) RETURN

          cmd_args_used = .FALSE.
          ! read in all args
          DO i=1,cmd_args_i
             CALL GET_COMMAND_ARGUMENT(start+i, cbuf)
             cmd_args_len(i) = LEN_TRIM(cbuf)
             cmd_args(i)     = cbuf(1:cmd_args_len(i))
          END DO

        END SUBROUTINE read_cmd_args
        !-------------------------------------------------------------------------
        !  Parse command line
        !-------------------------------------------------------------------------
        SUBROUTINE parse_cmd_line(info)

          IMPLICIT NONE

          INTEGER, INTENT(  OUT) :: info

          INTEGER :: i, ios

          CHARACTER(LEN=*),       PARAMETER :: caller = 'parse_cmd_line'
          CHARACTER(LEN=ppm_char)           :: cvar
          CHARACTER(LEN=256)                :: value

          LOGICAL :: ok
          LOGICAL :: err = .FALSE.
          ! scalar
#define DTYPE INTEGER
#include "ctrl/parse_arg.f"
#define DTYPE LONGINT
#include "ctrl/parse_arg.f"
#define DTYPE SINGLE
#include "ctrl/parse_arg.f"
#define DTYPE DOUBLE
#include "ctrl/parse_arg.f"
#define DTYPE LOGICAL
#define __LOGICAL
#include "ctrl/parse_arg.f"
#define DTYPE STRING
#include "ctrl/parse_arg.f"
#define DTYPE COMPLEX
#include "ctrl/parse_arg.f"
#define DTYPE DCOMPLEX
#include "ctrl/parse_arg.f"
        ! array
#define DTYPE INTEGER_array
#include "ctrl/parse_arg.f"
#define DTYPE LONGINT_array
#include "ctrl/parse_arg.f"
#define DTYPE SINGLE_array
#include "ctrl/parse_arg.f"
#define DTYPE DOUBLE_array
#include "ctrl/parse_arg.f"
#define DTYPE LOGICAL_array
#include "ctrl/parse_arg.f"
#define DTYPE STRING_array
#include "ctrl/parse_arg.f"
#define DTYPE COMPLEX_array
#include "ctrl/parse_arg.f"
#define DTYPE DCOMPLEX_array
#include "ctrl/parse_arg.f"
      9999 CONTINUE

        END SUBROUTINE parse_cmd_line
        !------------------------------------------------------------------------
        !  Parse Control file
        !------------------------------------------------------------------------
        SUBROUTINE parse_ctrl_file(info)

          IMPLICIT NONE

          INTEGER, INTENT(  OUT) :: info

          INTEGER :: ilenctrl
          INTEGER :: iUnit, istat, ios
          INTEGER :: iline
          INTEGER :: ilen
          INTEGER :: i,j,idx

          CHARACTER(LEN=*), PARAMETER :: caller='parse_ctrl_file'
          CHARACTER(LEN=ppm_char)     :: cbuf, cvalue, carg, cvar
          CHARACTER(LEN=ppm_char)     :: current_var
          CHARACTER(LEN=10000)        :: errmsg

          LOGICAL :: lExist

          !-------------------------------------------------------------
          !  Check that the ctrl file exists
          !-------------------------------------------------------------
          ilenctrl = LEN_TRIM(ctrl_file_name)
          IF (ilenctrl .LT. 1) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_argument, caller, &
             &    'No ctrl file given', __LINE__, info)
             GOTO 9999
          END IF

          INQUIRE(FILE=ctrl_file_name, EXIST=lExist)
          IF (.NOT. lExist) THEN
             WRITE(cvar,'(2A)') 'No such file: ', ctrl_file_name(1:ilenctrl)

             fail(cvar,ppm_error=ppm_error_fatal)
          END IF
          !-------------------------------------------------------------
          !  Open the file
          !-------------------------------------------------------------
          iUnit = 19
          OPEN(iUnit, FILE=ctrl_file_name, IOSTAT=ios, ACTION='READ')
          IF (ios .NE. 0) THEN
             WRITE(cvar,'(2A)') 'Failed to open file: ', ctrl_file_name(1:ilenctrl)

             fail(cvar,ppm_error=ppm_error_fatal)
          END IF
          !-------------------------------------------------------------
          !  Scan file
          !-------------------------------------------------------------
          iline = 0
          var_loop: DO
             iline = iline + 1
             !----------------------------------------------------------
             !  Read line
             !----------------------------------------------------------
             READ(iUnit,'(A)', END=100, ERR=200) cbuf
             ilen = LEN_TRIM(cbuf)
             !----------------------------------------------------------
             !  Skip comment or empty lines
             !----------------------------------------------------------
             IF (ilen .GT. 0 .AND. cbuf(1:1) .NE. '#') THEN
                !-------------------------------------------------------
                !  Remove space
                !-------------------------------------------------------
                j = 0
                DO i=1,ilen
                   IF (cbuf(i:i) .NE. ' ' .AND. cbuf(i:i) .NE. ACHAR(9)) THEN
                      j = j + 1
                      cbuf(j:j) = cbuf(i:i)
                   END IF
                END DO
                !-------------------------------------------------------
                !  Update length of string
                !-------------------------------------------------------
                ilen = j
                !-------------------------------------------------------
                !  Find position of =
                !-------------------------------------------------------
                idx = INDEX(cbuf, '=')
                !-------------------------------------------------------
                !  Exit if missing
                !-------------------------------------------------------
                IF (idx .LT. 0) THEN
                   WRITE(cvar,'(A,I5)') 'Incorrect line: ', iline

                   fail(cvar,ppm_error=ppm_error_fatal)
                ENDIF
                !-------------------------------------------------------
                !  Get argument and value
                !-------------------------------------------------------
                carg   = ADJUSTL(cbuf(1:idx-1))
                cvalue = ADJUSTL(cbuf(idx+1:ilen))
                !-------------------------------------------------------
                !  Convert to upper case
                !-------------------------------------------------------
                CALL UpperCase(carg, idx-1, info)
                !-------------------------------------------------------
                !  Find the variable
                !-------------------------------------------------------
                ! scalar
#define DTYPE INTEGER
#include "ctrl/parse_ctrl.f"
#define DTYPE LONGINT
#include "ctrl/parse_ctrl.f"
#define DTYPE SINGLE
#include "ctrl/parse_ctrl.f"
#define DTYPE DOUBLE
#include "ctrl/parse_ctrl.f"
#define DTYPE LOGICAL
#include "ctrl/parse_ctrl.f"
#define DTYPE STRING
#define __STRING
#include "ctrl/parse_ctrl.f"
#define DTYPE COMPLEX
#include "ctrl/parse_ctrl.f"
#define DTYPE DCOMPLEX
#include "ctrl/parse_ctrl.f"
        ! array
#define DTYPE INTEGER_array
#include "ctrl/parse_ctrl.f"
#define DTYPE LONGINT_array
#include "ctrl/parse_ctrl.f"
#define DTYPE SINGLE_array
#include "ctrl/parse_ctrl.f"
#define DTYPE DOUBLE_array
#include "ctrl/parse_ctrl.f"
#define DTYPE LOGICAL_array
#include "ctrl/parse_ctrl.f"
#define DTYPE STRING_array
#include "ctrl/parse_ctrl.f"
#define DTYPE COMPLEX_array
#include "ctrl/parse_ctrl.f"
#define DTYPE DCOMPLEX_array
#include "ctrl/parse_ctrl.f"
             END IF
          END DO var_loop
          GOTO 9999
      300 CONTINUE
          WRITE(errmsg,'(5A,I5,2A)') 'Error reading variable: ', &
          & current_var(1:LEN_TRIM(current_var)),                &
          & ' from string: ', cvalue(1:LEN_TRIM(cvalue)),        &
          & ' on line: ', iline,                                 &
          & ' of file: ', ctrl_file_name(1:ilenctrl)

          fail(errmsg,ppm_error=ppm_error_fatal)
      200 CONTINUE

          WRITE(errmsg,'(A,I5,2A)') 'Error reading line: ', iline, &
          & ' of file: ', ctrl_file_name(1:ilenctrl)

          fail(errmsg,ppm_error=ppm_error_fatal)
      100 CONTINUE
          !----------------------------------------------------------------------
          !  Close file
          !----------------------------------------------------------------------
      9999 CONTINUE
          CLOSE(iUnit)
          RETURN
        END SUBROUTINE parse_ctrl_file
#ifdef __F2003
        !------------------------------------------------------------------------
        !  Call default funcs
        !------------------------------------------------------------------------
        SUBROUTINE call_default_funcs(info)

          IMPLICIT NONE

          INTEGER, INTENT(  OUT) :: info

          INTEGER :: i

          CHARACTER(LEN=*), PARAMETER :: caller = 'call_default_funcs'

#define DTYPE INTEGER
#include "ctrl/default_func.f"
#define DTYPE LONGINT
#include "ctrl/default_func.f"
#define DTYPE SINGLE
#include "ctrl/default_func.f"
#define DTYPE DOUBLE
#include "ctrl/default_func.f"
#define DTYPE LOGICAL
#include "ctrl/default_func.f"
#define DTYPE STRING
#include "ctrl/default_func.f"
#define DTYPE COMPLEX
#include "ctrl/default_func.f"
#define DTYPE DCOMPLEX
#include "ctrl/default_func.f"
        ! array
#define DTYPE INTEGER_array
#include "ctrl/default_func.f"
#define DTYPE LONGINT_array
#include "ctrl/default_func.f"
#define DTYPE SINGLE_array
#include "ctrl/default_func.f"
#define DTYPE DOUBLE_array
#include "ctrl/default_func.f"
#define DTYPE LOGICAL_array
#include "ctrl/default_func.f"
#define DTYPE STRING_array
#include "ctrl/default_func.f"
#define DTYPE COMPLEX_array
#include "ctrl/default_func.f"
#define DTYPE DCOMPLEX_array
#include "ctrl/default_func.f"
      9999 CONTINUE

        END SUBROUTINE call_default_funcs
#endif
        !------------------------------------------------------------------------
        !  Check min max
        !------------------------------------------------------------------------
        SUBROUTINE check_minmax(info)

          IMPLICIT NONE

          INTEGER, INTENT(  OUT) :: info

          INTEGER :: i

          CHARACTER(LEN=*), PARAMETER :: caller = 'check_minmax'
          CHARACTER(LEN=ppm_char)     :: cvar

#define DTYPE INTEGER
#include "ctrl/minmax.f"
#define DTYPE LONGINT
#include "ctrl/minmax.f"
#define DTYPE SINGLE
#include "ctrl/minmax.f"
#define DTYPE DOUBLE
#include "ctrl/minmax.f"
      ! #define DTYPE STRING
      ! #define __STRING
      ! #include "ctrl/minmax.f"
#define ARRAY
#define DTYPE INTEGER_array
#include "ctrl/minmax.f"
#define DTYPE LONGINT_array
#include "ctrl/minmax.f"
#define DTYPE SINGLE_array
#include "ctrl/minmax.f"
#define DTYPE DOUBLE_array
#include "ctrl/minmax.f"
      ! #define DTYPE STRING_array
      ! #define __STRING
      ! #include "ctrl/minmax.f"
#undef ARRAY
      9999 CONTINUE
        END SUBROUTINE check_minmax
#ifdef __F2003
        !------------------------------------------------------------------------
        !  Call validator functions
        !------------------------------------------------------------------------
        SUBROUTINE call_validator_funcs(info)

          IMPLICIT NONE

          INTEGER, INTENT(  OUT) :: info

          INTEGER :: i

          CHARACTER(LEN=*), PARAMETER :: caller = 'call_validator_funcs'
          CHARACTER(LEN=ppm_char)     :: cvar
          ! scalar
#define DTYPE INTEGER
#include "ctrl/validate.f"
#define DTYPE LONGINT
#include "ctrl/validate.f"
#define DTYPE SINGLE
#include "ctrl/validate.f"
#define DTYPE DOUBLE
#include "ctrl/validate.f"
#define DTYPE LOGICAL
#include "ctrl/validate.f"
#define DTYPE STRING
#include "ctrl/validate.f"
#define DTYPE COMPLEX
#include "ctrl/validate.f"
#define DTYPE DCOMPLEX
#include "ctrl/validate.f"
        ! array
#define DTYPE INTEGER_array
#include "ctrl/validate.f"
#define DTYPE LONGINT_array
#include "ctrl/validate.f"
#define DTYPE SINGLE_array
#include "ctrl/validate.f"
#define DTYPE DOUBLE_array
#include "ctrl/validate.f"
#define DTYPE LOGICAL_array
#include "ctrl/validate.f"
#define DTYPE STRING_array
#include "ctrl/validate.f"
#define DTYPE COMPLEX_array
#include "ctrl/validate.f"
#define DTYPE DCOMPLEX_array
#include "ctrl/validate.f"
      9999 CONTINUE
        END SUBROUTINE call_validator_funcs
#endif
        !------------------------------------------------------------------------
        !  Special args
        !------------------------------------------------------------------------
        SUBROUTINE disable_help
        !!! Turns of help flag parsing.
          IMPLICIT NONE

          help_enabled = .FALSE.
        END SUBROUTINE disable_help

        SUBROUTINE disable_ctrl
        !!! Turns of control file parsing.

          IMPLICIT NONE

          ctrl_enabled = .FALSE.
        END SUBROUTINE disable_ctrl

        SUBROUTINE set_ctrl_name(name)
        !!! Sets the control file name.

          IMPLICIT NONE

          CHARACTER(LEN=*), INTENT(IN   ) :: name
          !!! Name of the control file.

          ctrl_file_name = name
        END SUBROUTINE set_ctrl_name
        !------------------------------------------------------------------------
        !  Set command line args (usefull for running tests)
        !------------------------------------------------------------------------
        SUBROUTINE add_cmd(arg, value)
        !!! Debugging and testing procedure. First call replaces the actual
        !!! command arguments with the supplied values. Can be called many
        !!! times to incrementally build the fake argument list.

          IMPLICIT NONE

          CHARACTER(LEN=*),           INTENT(IN   ) :: arg
          !!! Argument string to be appended to the list.
          CHARACTER(LEN=*), OPTIONAL, INTENT(IN   ) :: value
          !!! Optional value to append.

          INTEGER, POINTER, DIMENSION(:) :: temp_l => NULL()
          INTEGER                        :: inc
          INTEGER                        :: l

          CHARACTER(LEN=256), POINTER, DIMENSION(:) :: temp_a => NULL()

          in_test = .TRUE.
          inc = 1
          IF (PRESENT(value)) inc=2
          IF (.NOT. ASSOCIATED(cmd_args)) THEN
             ALLOCATE(cmd_args(1:di))
             ALLOCATE(cmd_args_len(1:di))
             ALLOCATE(cmd_args_used(1:di))
             cmd_args_used = .FALSE.
          ELSE
             l = SIZE(cmd_args)
             IF (l .LT. cmd_args_i + inc) THEN
                ALLOCATE(temp_a(1:(l + di)))
                ALLOCATE(temp_l(1:(l + di)))
                temp_a(1:l) = cmd_args(1:l)
                temp_l(1:l) = cmd_args_len(1:l)
                DEALLOCATE(cmd_args)
                DEALLOCATE(cmd_args_len)
                DEALLOCATE(cmd_args_used)
                ALLOCATE(cmd_args_used(1:(l+di)))
                cmd_args     => temp_a
                cmd_args_len => temp_l
                cmd_args_used = .FALSE.
             END IF
          END IF
          cmd_i = cmd_i + 1
          cmd_args_i = cmd_args_i + 1
          cmd_args(cmd_args_i)     =     arg
          cmd_args_len(cmd_args_i) = LEN(arg)
          IF (inc .EQ. 2) THEN
             cmd_i = cmd_i + 1
             cmd_args_i = cmd_args_i + 1
             cmd_args(cmd_args_i)     =     value
             cmd_args_len(cmd_args_i) = LEN(value)
          END IF
        END SUBROUTINE add_cmd
        !-------------------------------------------------------------------------
        !  Look-up flags
        !-------------------------------------------------------------------------
        SUBROUTINE find_flag(flag, success, value, err)
        !!! Returns the value supplied with the flag _flag_.

          IMPLICIT NONE

          CHARACTER(LEN=*),           INTENT(IN   ) :: flag
          !!! Flag string (eg. _'-f'_ or _'--flag'_).
          LOGICAL,                    INTENT(  OUT) :: success
          !!! True on success. False if there is no such flag or the value is
          !!! not supplied.
          CHARACTER(LEN=*), OPTIONAL, INTENT(  OUT) :: value
          !!! The value of supplied for the flag (character string).
          LOGICAL,          OPTIONAL, INTENT(  OUT) :: err
          !!! Returns .TRUE. if the flag was supplied without a value.

          INTEGER :: i
          INTEGER :: info

          CHARACTER(LEN=*), PARAMETER :: caller = 'find_flag'
          CHARACTER(LEN=ppm_char)     :: cvar

          success = .FALSE.
          IF (PRESENT(err)) err = .FALSE.
          IF (cmd_args_i .EQ. 0) RETURN
          DO i=1,cmd_args_i
             IF (cmd_args(i)(1:cmd_args_len(i)) .EQ. flag) THEN
                success          = .TRUE.
                IF (.NOT. cmd_args_used(i)) cmd_i = cmd_i - 1
                cmd_args_used(i) = .TRUE.
                IF (PRESENT(value)) THEN
                   IF (i+1 .GT. cmd_args_i) THEN
                      success = .FALSE.
                      IF (PRESENT(err)) err = .TRUE.
                      WRITE (cvar,*) "Flag ", flag, " requires an argument."
                      info = ppm_error_warning
                      CALL ppm_error(ppm_err_argument, caller, cvar,__LINE__, info)
                   ELSE
                      value = cmd_args(i+1)
                      IF (.NOT. cmd_args_used(i+1)) cmd_i = cmd_i - 1
                      cmd_args_used(i+1) = .TRUE.
                   END IF
                END IF
                RETURN
             END IF
          END DO
        END SUBROUTINE find_flag
        !------------------------------------------------------------------------
        !  Get arg count
        !------------------------------------------------------------------------
        INTEGER FUNCTION arg_count()
        !!! Returns the number of positional arguments found.

          IMPLICIT NONE

          arg_count = cmd_i
        END FUNCTION arg_count
        !------------------------------------------------------------------------
        !  Look-up args (after all flags have been read!!!)
        !------------------------------------------------------------------------
        SUBROUTINE find_arg(position, success, value)
        !!! Returns a positional argument at _position_.

          IMPLICIT NONE

          INTEGER,          INTENT(IN   )  :: position
          !!! Index of the positional argument, 1 based.
          LOGICAL,          INTENT(  OUT)  :: success
          !!! True on success. False if there is no positional arg with the
          !!! supplied index.
          CHARACTER(LEN=*), INTENT(INOUT)  :: value
          !!! The value of the positional argument (character string).

          INTEGER :: i, j

          j = 0
          success = .FALSE.
          DO i=1,cmd_args_i
             IF (.NOT. cmd_args_used(i)) THEN
                j = j + 1
                IF (position .EQ. j) THEN
                   value   = cmd_args(i)
                   success = .TRUE.
                   RETURN
                END IF
             END IF
          END DO
        END SUBROUTINE find_arg
        !------------------------------------------------------------------------
        !  Argument groups
        !------------------------------------------------------------------------
        SUBROUTINE arg_group(name)
        !!! Defines an argument group.

          IMPLICIT NONE

          CHARACTER(LEN=*), INTENT(IN   ) :: name
          !!! Group name

          INTEGER, DIMENSION(:), POINTER :: temp_s => NULL()
          INTEGER, DIMENSION(:), POINTER :: temp_m => NULL()

          CHARACTER(LEN=256), DIMENSION(:), POINTER :: temp_g => NULL()

          LOGICAL, DIMENSION(:), POINTER :: temp_c => NULL()
          LOGICAL, DIMENSION(:), POINTER :: temp_a => NULL()

          groups_i = groups_i + 1
          IF (.NOT. ASSOCIATED(groups)) THEN
             ALLOCATE(groups(0:di))
             ALLOCATE(group_size(0:di))
             ALLOCATE(group_max_len(0:di))
             ALLOCATE(group_has_ctrl(0:di))
             ALLOCATE(group_has_arg(0:di))
          ELSE IF (groups_i+1 .GT. SIZE(groups)) THEN
             ALLOCATE(temp_g(0:(SIZE(groups)+di)))
             ALLOCATE(temp_s(0:(SIZE(groups)+di)))
             ALLOCATE(temp_m(0:(SIZE(groups)+di)))
             ALLOCATE(temp_c(0:(SIZE(groups)+di)))
             ALLOCATE(temp_a(0:(SIZE(groups)+di)))
             temp_g(0:SIZE(groups)) = groups(0:SIZE(groups))
             temp_s(0:SIZE(groups)) = group_size(0:SIZE(groups))
             temp_m(0:SIZE(groups)) = group_max_len(0:SIZE(groups))
             temp_c(0:SIZE(groups)) = group_has_ctrl(0:SIZE(groups))
             temp_a(0:SIZE(groups)) = group_has_arg(0:SIZE(groups))
             DEALLOCATE(groups)
             DEALLOCATE(group_size)
             DEALLOCATE(group_max_len)
             DEALLOCATE(group_has_ctrl)
             DEALLOCATE(group_has_arg)
             groups         => temp_g
             group_size     => temp_s
             group_max_len  => temp_m
             group_has_ctrl => temp_c
             group_has_arg  => temp_a
          END IF
          groups(groups_i)         = name
          group_size(groups_i)     = 0
          group_max_len(groups_i)  = 0
          group_has_ctrl(groups_i) = .FALSE.
          group_has_arg(groups_i)  = .FALSE.
          IF (groups_i .EQ. 0 .AND. (help_enabled .OR. ctrl_enabled)) THEN
             group_has_arg(0) = .TRUE.
          END IF
        END SUBROUTINE arg_group
        !------------------------------------------------------------------------
        !  Help string break
        !------------------------------------------------------------------------
        SUBROUTINE break_help(ht, w, prefix, cbuft, caller, skip_first_line)

          IMPLICIT NONE
          ! args
          CHARACTER(LEN=*),  INTENT(IN   ) :: ht
          INTEGER,           INTENT(IN   ) :: w
          CHARACTER(LEN=*),  INTENT(IN   ) :: prefix
          CHARACTER(LEN=*),  INTENT(INOUT) :: cbuft
          CHARACTER(LEN=*),  INTENT(IN   ) :: caller
          LOGICAL, OPTIONAL, INTENT(IN   ) :: skip_first_line

          ! vars
          INTEGER :: start
          INTEGER :: break
          INTEGER :: next
          INTEGER :: end
          INTEGER :: skip
          INTEGER :: lng,info

          CHARACTER(LEN=ppm_char) :: cbuf1=""

          LOGICAL :: fl

          IF (PRESENT(skip_first_line)) THEN
             fl = skip_first_line
          ELSE
             fl = .TRUE.
          END IF
          start = 1
          break = 1
          next  = 1
          end   = LEN(ht)
          ! line loop
          DO WHILE (start .LT. end)
             ! break < start + w < next
             DO WHILE (next .LT. end .AND. next .LT. start+w)
                ! find next blank or new line
                skip = 1
                DO WHILE (next .LT. end .AND. ht(next:next) .NE. ' ')
                   IF (next .LT. end) THEN
                      IF (ht(next:next+1) .EQ. '\n') THEN
                         ! line break
                         skip = 3
                         break = next - 1
                         next = next + 2
                         GOTO 1234
                      END IF
                   END IF
                   next = next + 1
                END DO
                break = next
                next = next + 1
             END DO
             ! print line
      1234   IF (fl) THEN
                fl = .FALSE.
                WRITE(cbuf1,'(3A)') cbuft,ht(start:break),''
             ELSE
                WRITE(cbuf1,'(4A)') cbuft,prefix,ht(start:break),''
             END IF
             cbuf1=ADJUSTL(cbuf1)
             stdout(cbuf1)
             cbuft=''
             ! iterate
             start = break + skip
             break = next
          END DO
        END SUBROUTINE break_help
        !------------------------------------------------------------------------
        !  Print command line help
        !------------------------------------------------------------------------
        SUBROUTINE print_help
        !!! Prints usage information to the standard output.

          IMPLICIT NONE

          INTEGER :: i, j, k, l, info

          CHARACTER(LEN=*), PARAMETER :: caller="print_help"
          CHARACTER(LEN=ppm_char)     :: scratch,cbuf1

          IF (ctrl_enabled) THEN
             stdout("Usage: progname {ctrl-file} [options] [args]")
          ELSE
             stdout("Usage: progname [options] [args]")
          END IF
          DO k=0,groups_i
             IF (.NOT. group_has_arg(k)) CYCLE
             stdout("")
             WRITE(cbuf,'(A)')groups(k)
             stdout(cbuf)

             IF (k .EQ. 0) THEN
                IF (help_enabled) THEN
                   stdout("  help                          Print this help message and exit.")
                   stdout("                  short flag :  -h")
                   stdout("                   long flag :  --help")
                   stdout("")
                END IF
                IF (ctrl_enabled) THEN
                   stdout("  control file                  Print sample control file and exit.")
                   stdout("                   long flag :  --print-ctrl")
                   stdout("")
                END IF
             END IF
             group_loop: DO i=1,group_size(k)
                ! scalar
#define DTYPE INTEGER
#define __INTEGER
#include "ctrl/help.f"
#define DTYPE LONGINT
#define __LONGINT
#include "ctrl/help.f"
#define DTYPE SINGLE
#define __SINGLE
#include "ctrl/help.f"
#define DTYPE DOUBLE
#define __DOUBLE
#include "ctrl/help.f"
#define DTYPE LOGICAL
#define __LOGICAL
#include "ctrl/help.f"
#define DTYPE STRING
#define __STRING
#include "ctrl/help.f"
#define DTYPE COMPLEX
#include "ctrl/help.f"
#define DTYPE DCOMPLEX
#include "ctrl/help.f"
                ! array
#define ARRAY
#define DTYPE INTEGER_array
#define __INTEGER
#include "ctrl/help.f"
#define DTYPE LONGINT_array
#define __LONGINT
#include "ctrl/help.f"
#define DTYPE SINGLE_array
#define __SINGLE
#include "ctrl/help.f"
#define DTYPE DOUBLE_array
#define __DOUBLE
#include "ctrl/help.f"
#define DTYPE LOGICAL_array
#define __LOGICAL
#include "ctrl/help.f"
#define DTYPE STRING_array
#define __STRING
#include "ctrl/help.f"
#define DTYPE COMPLEX_array
#include "ctrl/help.f"
#define DTYPE DCOMPLEX_array
#include "ctrl/help.f"
#undef ARRAY
             END DO group_loop
          END DO
        END SUBROUTINE print_help
        !------------------------------------------------------------------------
        !  Print sample control file
        !------------------------------------------------------------------------
        SUBROUTINE print_ctrl
        !!! Prints a sample control file to standard output.

          IMPLICIT NONE

          INTEGER :: i, j, k, l, lng, info

          CHARACTER(LEN=*), PARAMETER :: caller="print_ctrl"
          CHARACTER(LEN=ppm_char) :: scratch,cbuf1

          LOGICAL :: in_line

          IF (groups_i .EQ. -1) THEN
             stdout("No args have been defined...")
             RETURN
          END IF
          stdout("#-------------------------------------------------------------------------------")
          stdout("#  Sample control file for ppm_client")
          stdout("#")
          stdout("#  Edit the settings below")
          stdout("#")
          stdout("#-------------------------------------------------------------------------------")
          DO k=0,groups_i
             IF (.NOT. group_has_ctrl(k)) CYCLE
             stdout("")
             stdout("#-------------------------------------------------------------------------------")
             stdout("#  ", 'groups(k)(1:LEN_TRIM(groups(k)))')

             comment_loop: DO i=1,group_size(k)
                ! scalar
#define DTYPE INTEGER
#define __INTEGER
#include "ctrl/ctrl_comment.f"
#define DTYPE LONGINT
#define __LONGINT
#include "ctrl/ctrl_comment.f"
#define DTYPE SINGLE
#define __SINGLE
#include "ctrl/ctrl_comment.f"
#define DTYPE DOUBLE
#define __DOUBLE
#include "ctrl/ctrl_comment.f"
#define DTYPE LOGICAL
#define __LOGICAL
#include "ctrl/ctrl_comment.f"
#define DTYPE STRING
#define __STRING
#include "ctrl/ctrl_comment.f"
#define DTYPE COMPLEX
#include "ctrl/ctrl_comment.f"
#define DTYPE DCOMPLEX
#include "ctrl/ctrl_comment.f"
                ! array
#define ARRAY
#define DTYPE INTEGER_array
#define __INTEGER
#include "ctrl/ctrl_comment.f"
#define DTYPE LONGINT_array
#define __LONGINT
#include "ctrl/ctrl_comment.f"
#define DTYPE SINGLE_array
#define __SINGLE
#include "ctrl/ctrl_comment.f"
#define DTYPE DOUBLE_array
#define __DOUBLE
#include "ctrl/ctrl_comment.f"
#define DTYPE LOGICAL_array
#define __LOGICAL
#include "ctrl/ctrl_comment.f"
#define DTYPE STRING_array
#define __STRING
#include "ctrl/ctrl_comment.f"
#define DTYPE COMPLEX_array
#include "ctrl/ctrl_comment.f"
#define DTYPE DCOMPLEX_array
#include "ctrl/ctrl_comment.f"
#undef ARRAY
             END DO comment_loop

             stdout("#-------------------------------------------------------------------------------")

             var_loop: DO i=1,group_size(k)
                ! scalar
#define DTYPE INTEGER
#include "ctrl/ctrl.f"
#define DTYPE LONGINT
#include "ctrl/ctrl.f"
#define DTYPE SINGLE
#include "ctrl/ctrl.f"
#define DTYPE DOUBLE
#include "ctrl/ctrl.f"
#define DTYPE LOGICAL
#include "ctrl/ctrl.f"
#define DTYPE STRING
#define __STRING
#include "ctrl/ctrl.f"
#define DTYPE COMPLEX
#include "ctrl/ctrl.f"
#define DTYPE DCOMPLEX
#include "ctrl/ctrl.f"
                ! array
#define ARRAY
#define DTYPE INTEGER_array
#include "ctrl/ctrl.f"
#define DTYPE LONGINT_array
#include "ctrl/ctrl.f"
#define DTYPE SINGLE_array
#include "ctrl/ctrl.f"
#define DTYPE DOUBLE_array
#include "ctrl/ctrl.f"
#define DTYPE LOGICAL_array
#include "ctrl/ctrl.f"
#define DTYPE STRING_array
#define __STRING
#include "ctrl/ctrl.f"
#define DTYPE COMPLEX_array
#include "ctrl/ctrl.f"
#define DTYPE DCOMPLEX_array
#include "ctrl/ctrl.f"
#undef ARRAY
             END DO var_loop

          END DO
        END SUBROUTINE print_ctrl
        !------------------------------------------------------------------------
        !  Debug
        !------------------------------------------------------------------------
      !   SUBROUTINE dump_defines
      ! !!! Debugging procedure. Prints all defined args.
      !     INTEGER :: i
      !     WRITE (*,'(//A)') "DUMPING DEFINED ARGS"
      !     WRITE (*,'(//A/)') "Integers: "
      !     DO i=1,INTEGER_args_i
      ! #define DTYPE INTEGER
      ! #include "ctrl/dump.f"
      ! #undef DTYPE
      !     END DO
      !     WRITE (*,'(//A/)') "Reals: "
      !     DO i=1,REAL_args_i
      ! #define DTYPE REAL
      ! #include "ctrl/dump.f"
      ! #undef DTYPE
      !     END DO
      !     WRITE (*,'(//A/)') "Chars: "
      !     DO i=1,CHAR_args_i
      ! #define DTYPE CHAR
      ! #define STRING
      ! #include "ctrl/dump.f"
      ! #undef STRING
      ! #undef DTYPE
      !     END DO
      !     WRITE (*,'(//A/)') "Logicals: "
      !     DO i=1,LOGICAL_args_i
      ! #define DTYPE LOGICAL
      ! #define BOOL
      ! #include "ctrl/dump.f"
      ! #undef BOOL
      ! #undef DTYPE
      !     END DO
      !   END SUBROUTINE dump_defines
        SUBROUTINE dump_args
        !!! Debugging procedure. Prints all supplied command line args.

          IMPLICIT NONE

          INTEGER :: i,info

          CHARACTER(LEN=*), PARAMETER :: caller="dump_args"

          stdout("")
          stdout("")
          stdout("DUMPING COMMAND ARGS")
          stdout("")
          DO i=1,cmd_args_i
             WRITE (cbuf,*) 'arg : ', cmd_args(i)(1:cmd_args_len(i))
             stdout(cbuf)

             IF (cmd_args_used(i)) THEN
                stdout("used")
             ELSE
                stdout("not used")
             END IF
          END DO
        END SUBROUTINE dump_args
        !------------------------------------------------------------------------
        !  Adders - black or otherwise
        !------------------------------------------------------------------------

        ! scalar
#define DTYPE INTEGER
#define __INTEGER
#include "ctrl/adder.f"

#define DTYPE LONGINT
#define __LONGINT
#include "ctrl/adder.f"

#define DTYPE SINGLE
#define __SINGLE
#include "ctrl/adder.f"

#define DTYPE DOUBLE
#define __DOUBLE
#include "ctrl/adder.f"

#define DTYPE LOGICAL
#define __LOGICAL
#include "ctrl/adder.f"

#define DTYPE STRING
#define __STRING
#include "ctrl/adder.f"

#define DTYPE COMPLEX
#define __COMPLEX
#include "ctrl/adder.f"

#define DTYPE DCOMPLEX
#define __DCOMPLEX
#include "ctrl/adder.f"

        ! array
#define ARRAY

#define DTYPE INTEGER_array
#define __INTEGER
#include "ctrl/adder.f"

#define DTYPE LONGINT_array
#define __LONGINT
#include "ctrl/adder.f"

#define DTYPE SINGLE_array
#define __SINGLE
#include "ctrl/adder.f"

#define DTYPE DOUBLE_array
#define __DOUBLE
#include "ctrl/adder.f"

#define DTYPE LOGICAL_array
#define __LOGICAL
#include "ctrl/adder.f"

#define DTYPE STRING_array
#define __STRING
#include "ctrl/adder.f"

#define DTYPE COMPLEX_array
#define __COMPLEX
#include "ctrl/adder.f"

#define DTYPE DCOMPLEX_array
#define __DCOMPLEX
#include "ctrl/adder.f"

#undef ARRAY

        !-------------------------------------------------------------------------
        !  My own uppercase!!! YAY!!!
        !-------------------------------------------------------------------------
        SUBROUTINE UpperCase(string,ilen,info)
          IMPLICIT NONE
          CHARACTER(LEN=*), INTENT(INOUT) :: string
          INTEGER         , INTENT(IN   ) :: ilen
          INTEGER         , INTENT(  OUT) :: info
          INTEGER          :: i,j
          INTEGER          :: i1,i2,i3,iadd
          info = 0
          i1   = IACHAR('a') - 1
          i2   = IACHAR('z') + 1
          i3   = IACHAR('A')
          iadd = i3 - i1 - 1
          DO i=1,ilen
             j = IACHAR(string(i:i))
             IF (j.GT.i1.AND.j.LT.i2) THEN
                string(i:i) = CHAR(j+iadd)
             ENDIF
          ENDDO
          RETURN
        END SUBROUTINE UpperCase

!         !    %----------------------------------------------------------
!         !    | Converts the string "s1" to upper case
!         !    %----------------------------------------------------------
!         SUBROUTINE Uppercase(record)
!           !    %-----------------------------------------------------------------%
!           !    | author             - y.afshar         January     2010
!           !    %-----------------------------------------------------------------%
!           IMPLICIT NONE
!
!           CHARACTER(LEN = *), INTENT(INOUT) :: record
!
!           INTEGER :: length, LowerToUpper, i
!
!           CHARACTER(LEN = LEN(record)) :: temp
!
!           LowerToUpper = IACHAR("A") - IACHAR("a")
!
!           temp = TRIM(ADJUSTL(record))
!           length = LEN(temp)
!
!           record = ''
!           DO i=1, length
!              SELECT CASE (temp(i:i))
!              CASE ('a' : 'z')
!                 record(i:i) = ACHAR(IACHAR(temp(i:i)) + LowerToUpper)
!              CASE DEFAULT
!                 record(i:i) = temp(i:i)
!              END SELECT
!           ENDDO
!           !-------------------------------------------------------------------------
!           !  Return
!           !-------------------------------------------------------------------------
!  9999     CONTINUE
!           RETURN
!         END SUBROUTINE Uppercase

      END MODULE ppm_module_ctrl
