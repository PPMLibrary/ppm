!--------------------------------------------------------------------------
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
  !------------------------------------------------------------------------
  !  Modules
  !------------------------------------------------------------------------
  USE ppm_module_typedef
  USE ppm_module_data, ONLY: ppm_rank, ppm_comm
  IMPLICIT NONE
  !------------------------------------------------------------------------
  !  Interface
  !------------------------------------------------------------------------
  PUBLIC :: arg, arg_group, parse_args, disable_help, disable_ctrl, &
       &    set_ctrl_name,                                          &
       &    integer_dflt, real_dflt, char_dflt, logical_dflt,       & 
       &    integer_array_dflt, real_array_dflt,                    &
       &    char_array_dflt, logical_array_dflt,                    &
       &    integer_vdtr, real_vdtr, char_vdtr, logical_vdtr,       &
       &    integer_array_vdtr, real_array_vdtr,                    &
       &    char_array_vdtr, logical_array_vdtr,                    &
       &    reset, add_cmd, ctrl_file_test,                         &
       &    find_arg, find_flag, arg_count,                         &
       &    enabling_flag, disabling_flag

  PRIVATE
#ifdef __MPI
  INCLUDE 'mpif.h'
#endif
  !------------------------------------------------------------------------
  !  Types
  !------------------------------------------------------------------------

! scalars
#define DTYPE INTEGER
#include "ctrl/type.f"
#undef DTYPE
#define DTYPE REAL
#include "ctrl/type.f"
#undef DTYPE
#define DTYPE CHAR
#define STRING
#include "ctrl/type.f"
#undef STRING
#undef DTYPE
#define DTYPE LOGICAL
#define BOOL
#include "ctrl/type.f"
#undef BOOL
#undef DTYPE

! array versions
#define ARRAY
#define DTYPE INTEGER
#include "ctrl/type.f"
#undef DTYPE
#define DTYPE REAL
#include "ctrl/type.f"
#undef DTYPE
#define DTYPE CHAR
#define STRING
#include "ctrl/type.f"
#undef STRING
#undef DTYPE
#define DTYPE LOGICAL
#define BOOL
#include "ctrl/type.f"
#undef BOOL
#undef DTYPE
#undef ARRAY
  !------------------------------------------------------------------------
  !  Interfaces
  !------------------------------------------------------------------------
  INTERFACE arg
     ! scalar
     MODULE PROCEDURE INTEGER_add_arg
     MODULE PROCEDURE REAL_add_arg
     MODULE PROCEDURE LOGICAL_add_arg
     MODULE PROCEDURE CHAR_add_arg
     ! array
     MODULE PROCEDURE INTEGER_array_add_arg
     MODULE PROCEDURE REAL_array_add_arg
     MODULE PROCEDURE LOGICAL_array_add_arg
     MODULE PROCEDURE CHAR_array_add_arg
  END INTERFACE

  ABSTRACT INTERFACE
     !---------------------------------------------------------------------
     !  Default functions
     !---------------------------------------------------------------------
     ! scalar
     LOGICAL FUNCTION INTEGER_dflt(variable)
       INTEGER, POINTER :: variable
     END FUNCTION INTEGER_dflt
     LOGICAL FUNCTION REAL_dflt(variable)
       REAL, POINTER :: variable
     END FUNCTION REAL_dflt
     LOGICAL FUNCTION LOGICAL_dflt(variable)
       LOGICAL, POINTER :: variable
     END FUNCTION LOGICAL_dflt
     LOGICAL FUNCTION CHAR_dflt(variable)
       CHARACTER(LEN=256), POINTER :: variable
     END FUNCTION CHAR_dflt
     ! array
     LOGICAL FUNCTION INTEGER_array_dflt(variable)
       INTEGER, DIMENSION(:), POINTER :: variable
     END FUNCTION INTEGER_array_dflt
     LOGICAL FUNCTION REAL_array_dflt(variable)
       REAL, DIMENSION(:), POINTER :: variable
     END FUNCTION REAL_array_dflt
     LOGICAL FUNCTION LOGICAL_array_dflt(variable)
       LOGICAL, DIMENSION(:), POINTER :: variable
     END FUNCTION LOGICAL_array_dflt
     LOGICAL FUNCTION CHAR_array_dflt(variable)
       CHARACTER(LEN=*), DIMENSION(:), POINTER :: variable
     END FUNCTION CHAR_array_dflt
     !---------------------------------------------------------------------
     !  Validators
     !---------------------------------------------------------------------
     ! scalar
     LOGICAL FUNCTION INTEGER_vdtr(variable)
       INTEGER, POINTER           :: variable
     END FUNCTION INTEGER_vdtr
     LOGICAL FUNCTION REAL_vdtr(variable)
       REAL, POINTER              :: variable
     END FUNCTION REAL_vdtr
     LOGICAL FUNCTION LOGICAL_vdtr(variable)
       LOGICAL, POINTER           :: variable
     END FUNCTION LOGICAL_vdtr
     LOGICAL FUNCTION CHAR_vdtr(variable)
       CHARACTER(LEN=*), POINTER  :: variable
     END FUNCTION CHAR_vdtr
     ! arrays
     LOGICAL FUNCTION INTEGER_array_vdtr(variable)
       INTEGER, DIMENSION(:), POINTER           :: variable
     END FUNCTION INTEGER_array_vdtr
     LOGICAL FUNCTION REAL_array_vdtr(variable)
       REAL, DIMENSION(:), POINTER              :: variable
     END FUNCTION REAL_array_vdtr
     LOGICAL FUNCTION LOGICAL_array_vdtr(variable)
       LOGICAL, DIMENSION(:), POINTER           :: variable
     END FUNCTION LOGICAL_array_vdtr
     LOGICAL FUNCTION CHAR_array_vdtr(variable)
       CHARACTER(LEN=*), DIMENSION(:), POINTER  :: variable
     END FUNCTION CHAR_array_vdtr
  END INTERFACE
  !------------------------------------------------------------------------
  !  Constants
  !------------------------------------------------------------------------
  LOGICAL, PARAMETER       :: enabling_flag=.true.
  LOGICAL, PARAMETER       :: disabling_flag=.false.
  !------------------------------------------------------------------------
  !  Variables
  !------------------------------------------------------------------------
  ! by how much to grow storage
  INTEGER,             PARAMETER               :: di = 10
  ! scalar
  TYPE(INTEGER_arg),   POINTER, DIMENSION(:)   :: INTEGER_args   => NULL()
  INTEGER                                      :: INTEGER_args_i = 0
  TYPE(REAL_arg),      POINTER, DIMENSION(:)   :: REAL_args      => NULL()
  INTEGER                                      :: REAL_args_i    = 0
  TYPE(LOGICAL_arg),   POINTER, DIMENSION(:)   :: LOGICAL_args   => NULL()
  INTEGER                                      :: LOGICAL_args_i = 0
  TYPE(CHAR_arg),      POINTER, DIMENSION(:)   :: CHAR_args      => NULL()
  INTEGER                                      :: CHAR_args_i    = 0
  ! add arrays
  TYPE(INTEGER_array_arg), POINTER, DIMENSION(:)   :: INTEGER_array_args   => NULL()
  INTEGER                                      :: INTEGER_array_args_i = 0
  TYPE(REAL_array_arg),    POINTER, DIMENSION(:)   :: REAL_array_args      => NULL()
  INTEGER                                      :: REAL_array_args_i    = 0
  TYPE(LOGICAL_array_arg), POINTER, DIMENSION(:)   :: LOGICAL_array_args   => NULL()
  INTEGER                                      :: LOGICAL_array_args_i = 0
  TYPE(CHAR_array_arg),    POINTER, DIMENSION(:)   :: CHAR_array_args      => NULL()
  INTEGER                                      :: CHAR_array_args_i    = 0
  ! arg storage
  CHARACTER(LEN=256),  POINTER, DIMENSION(:)   :: cmd_args       => NULL()
  INTEGER,             POINTER, DIMENSION(:)   :: cmd_args_len   => NULL()
  LOGICAL,             POINTER, DIMENSION(:)   :: cmd_args_used  => NULL()
  INTEGER                                      :: cmd_args_i=0
  INTEGER                                      :: cmd_i=0
  ! arg groups
  CHARACTER(LEN=256),  POINTER, DIMENSION(:)   :: groups         => NULL()
  INTEGER,             POINTER, DIMENSION(:)   :: group_size     => NULL()
  LOGICAL,             POINTER, DIMENSION(:)   :: group_has_ctrl => NULL()
  LOGICAL,             POINTER, DIMENSION(:)   :: group_has_arg  => NULL()
  INTEGER                                      :: groups_i=-1
  ! special args
  LOGICAL                                      :: help_enabled   = .TRUE.
  LOGICAL                                      :: ctrl_enabled   = .TRUE.
  CHARACTER(LEN=256)                           :: ctrl_file_name = 'Ctrl'
  CHARACTER(LEN=ppm_char)                      :: ctrl_file_test = ''
  ! test run
  LOGICAL                                      :: in_test        = .FALSE.

CONTAINS
  !------------------------------------------------------------------------
  !  Master procedure
  !------------------------------------------------------------------------
  SUBROUTINE parse_args(info)
    !----------------------------------------------------------------------
    !  Arguments
    !----------------------------------------------------------------------
    INTEGER, INTENT(  OUT)                  :: info
    !----------------------------------------------------------------------
    !  Local variables
    !----------------------------------------------------------------------
    REAL(8)                                 :: t0
    CHARACTER(LEN=*), PARAMETER             :: caller='parse_args'
    CHARACTER(LEN=256)                      :: value
    LOGICAL                                 :: ok
    INTEGER                                 :: info2
    INTEGER                                 :: i
    !----------------------------------------------------------------------
    !  Externals
    !----------------------------------------------------------------------
    EXTERNAL iargc
    INTEGER  iargc
    !----------------------------------------------------------------------
    !  Initialize
    !----------------------------------------------------------------------
!    CALL substart(caller, t0, info)
    !----------------------------------------------------------------------
    !  Do everything on rank 0 and bcast at the end
    !----------------------------------------------------------------------
    IF (ppm_rank .EQ. 0) THEN
       !-------------------------------------------------------------------
       !  Copy default values into variables
       !-------------------------------------------------------------------
       CALL apply_defaults(info)
       !-------------------------------------------------------------------
       !  Read in the command line
       !-------------------------------------------------------------------
       IF (.NOT. in_test) CALL read_cmd_args
       !-------------------------------------------------------------------
       !  Parse help flag
       !-------------------------------------------------------------------
       IF (help_enabled) THEN
          CALL find_flag('-h', ok)
          IF (.NOT. ok) CALL find_flag('--help', ok)
          IF (ok) THEN
             CALL print_help
             info = 1
             GOTO 9999
          END IF
       END IF
       !-------------------------------------------------------------------
       !  Parse rest of the command line
       !-------------------------------------------------------------------
       CALL parse_cmd_line(info)
       IF (info .NE. 0) GOTO 9999
       !-------------------------------------------------------------------
       !  Control file
       !-------------------------------------------------------------------
       IF (ctrl_enabled) THEN
          CALL find_flag('--print-ctrl', ok)
          IF (ok) THEN
             CALL print_ctrl
             info = 1
             GOTO 9999
          END IF
          CALL find_arg(1, ok, ctrl_file_name)
          !----------------------------------------------------------------
          !  Parse control file
          !----------------------------------------------------------------
          CALL parse_ctrl_file(info)
       END IF
       !-------------------------------------------------------------------
       !  Call default funcs
       !-------------------------------------------------------------------
       CALL call_default_funcs(info)
       IF (info .NE. 0) GOTO 9999
       !-------------------------------------------------------------------
       !  Check minmax
       !-------------------------------------------------------------------
       CALL check_minmax(info)
       IF (info .NE. 0) GOTO 9999
       !-------------------------------------------------------------------
       !  Run validators
       !-------------------------------------------------------------------
       CALL call_validator_funcs(info)
       IF (info .NE. 0) GOTO 9999
       !-------------------------------------------------------------------
       !  DONE!
       !-------------------------------------------------------------------
    END IF ! (ppm_rank .EQ. 0)
    !----------------------------------------------------------------------
    !  Exchange data
    !----------------------------------------------------------------------
100 CONTINUE
#ifdef __MPI
!    CALL MPI_BCast(what, length, MPI_TYPE, 0, comm, info)
!     CALL MPI_BCast(info, 1, MPI_INTEGER, 0, ppm_comm, info2)
!     IF (info .NE. 0) GOTO 9999
!     DO i=1,INTEGER_args_i
!        CALL MPI_BCast(INTEGER_args(i)%variable, 1, MPI_INTEGER, 0, ppm_comm, info)
!     END DO
!     DO i=1,REAL_args_i
!        CALL MPI_BCast(REAL_args(i)%variable, 1, MPI_REAL, 0, ppm_comm, info)
!     END DO
!     DO i=1,CHAR_args_i
!        CALL MPI_BCast(CHAR_args(i)%variable, &
!             LEN_TRIM(CHAR_args(i)%variable), MPI_CHARACTER, 0, ppm_comm, info)
!     END DO
!     DO i=1,LOGICAL_args_i
!        CALL MPI_BCast(LOGICAL_args(i)%variable, 1, MPI_LOGICAL, 0, ppm_comm, info)
!     END DO
#endif
    !----------------------------------------------------------------------
    !  Error handling
    !----------------------------------------------------------------------
9999 CONTINUE
    ! cleanup
    CALL deallocate_memory
!    CALL substop(caller, t0, info)
    RETURN
  END SUBROUTINE parse_args
  !------------------------------------------------------------------------
  !  Cleanup
  !------------------------------------------------------------------------
  SUBROUTINE deallocate_memory
    IF (ASSOCIATED(INTEGER_args))       DEALLOCATE(INTEGER_args)
    IF (ASSOCIATED(REAL_args))          DEALLOCATE(REAL_args)
    IF (ASSOCIATED(LOGICAL_args))       DEALLOCATE(LOGICAL_args)
    IF (ASSOCIATED(CHAR_args))          DEALLOCATE(CHAR_args)
    IF (ASSOCIATED(INTEGER_array_args)) DEALLOCATE(INTEGER_array_args)
    IF (ASSOCIATED(REAL_array_args))    DEALLOCATE(REAL_array_args)
    IF (ASSOCIATED(LOGICAL_array_args)) DEALLOCATE(LOGICAL_array_args)
    IF (ASSOCIATED(CHAR_array_args))    DEALLOCATE(CHAR_array_args)
    IF (ASSOCIATED(groups))             DEALLOCATE(groups)
    IF (ASSOCIATED(group_size))         DEALLOCATE(group_size)
    IF (ASSOCIATED(group_has_ctrl))     DEALLOCATE(group_has_ctrl)
    IF (ASSOCIATED(group_has_arg))      DEALLOCATE(group_has_arg)
  END SUBROUTINE deallocate_memory

  SUBROUTINE reset
    CALL deallocate_memory
    IF (ASSOCIATED(cmd_args))           DEALLOCATE(cmd_args)
    IF (ASSOCIATED(cmd_args_len))       DEALLOCATE(cmd_args_len)
    IF (ASSOCIATED(cmd_args_used))      DEALLOCATE(cmd_args_used)
    INTEGER_args_i       = 0
    REAL_args_i          = 0
    LOGICAL_args_i       = 0
    CHAR_args_i          = 0
    INTEGER_array_args_i = 0
    REAL_array_args_i    = 0
    LOGICAL_array_args_i = 0
    CHAR_array_args_i    = 0
    cmd_args_i           = 0
    cmd_i                = 0
    groups_i             = -1
    help_enabled         = .TRUE.
    ctrl_enabled         = .TRUE.
    ctrl_file_name       = 'Ctrl'
    ctrl_file_test       = ''
    in_test              = .FALSE.
  END SUBROUTINE reset
  !-------------------------------------------------------------------------
  !  Apply defaults
  !-------------------------------------------------------------------------
  SUBROUTINE apply_defaults(info)
    INTEGER, INTENT(  OUT) :: info
    INTEGER                :: i
    ! scalar
    DO i=1,INTEGER_args_i
       IF (INTEGER_args(i)%default_set) &
            INTEGER_args(i)%variable = INTEGER_args(i)%default
    END DO
    DO i=1,REAL_args_i
       IF (REAL_args(i)%default_set) &
            REAL_args(i)%variable = REAL_args(i)%default
    END DO
    DO i=1,CHAR_args_i
       IF (CHAR_args(i)%default_set) &
            CHAR_args(i)%variable = CHAR_args(i)%default
    END DO
    DO i=1,LOGICAL_args_i
       IF (LOGICAL_args(i)%default_set) &
            LOGICAL_args(i)%variable = LOGICAL_args(i)%default
    END DO
    ! array
    DO i=1,INTEGER_array_args_i
       IF (INTEGER_array_args(i)%default_set) &
            INTEGER_array_args(i)%variable = INTEGER_array_args(i)%default
    END DO
    DO i=1,REAL_array_args_i
       IF (REAL_array_args(i)%default_set) &
            REAL_array_args(i)%variable = REAL_array_args(i)%default
    END DO
    DO i=1,CHAR_array_args_i
       IF (CHAR_array_args(i)%default_set) &
            CHAR_array_args(i)%variable = CHAR_array_args(i)%default
    END DO
    DO i=1,LOGICAL_array_args_i
       IF (LOGICAL_array_args(i)%default_set) &
            LOGICAL_array_args(i)%variable = LOGICAL_array_args(i)%default
    END DO
  END SUBROUTINE apply_defaults
  !------------------------------------------------------------------------
  !  Read in command line args
  !------------------------------------------------------------------------
  SUBROUTINE read_cmd_args
    CHARACTER(LEN=ppm_char)         :: cbuf
    INTEGER                         :: i, start, nargc
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
       ALLOCATE(cmd_args(1:cmd_args_i))
       ALLOCATE(cmd_args_len(1:cmd_args_i))
       ALLOCATE(cmd_args_used(1:cmd_args_i))
       cmd_args_used = .FALSE.
       ! read in all args
       DO i=1,cmd_args_i
          CALL getarg(start+i, cbuf)
          cmd_args_len(i) = LEN_TRIM(cbuf)
          cmd_args(i)     = cbuf(1:cmd_args_len(i))
       END DO
  END SUBROUTINE read_cmd_args
  !-------------------------------------------------------------------------
  !  Parse command line
  !-------------------------------------------------------------------------
  SUBROUTINE parse_cmd_line(info)
    INTEGER, INTENT(  OUT) :: info
    INTEGER                :: i, ios
    LOGICAL                :: ok
    CHARACTER(LEN=256)     :: value
    ! scalar
#define DTYPE INTEGER
#include "ctrl/parse_arg.f"
#undef DTYPE
#define DTYPE REAL
#include "ctrl/parse_arg.f"
#undef DTYPE
#define DTYPE CHAR
#include "ctrl/parse_arg.f"
#undef DTYPE
#define DTYPE LOGICAL
#define BOOL
#include "ctrl/parse_arg.f"
#undef BOOL
#undef DTYPE
    ! array
#define ARRAY
#define DTYPE INTEGER_array
#include "ctrl/parse_arg.f"
#undef DTYPE
#define DTYPE REAL_array
#include "ctrl/parse_arg.f"
#undef DTYPE
#define DTYPE CHAR_array
#include "ctrl/parse_arg.f"
#undef DTYPE
#define DTYPE LOGICAL_array
#define BOOL
#include "ctrl/parse_arg.f"
#undef BOOL
#undef DTYPE
#undef ARRAY
9999 CONTINUE
  END SUBROUTINE parse_cmd_line
  !------------------------------------------------------------------------
  !  Parse Control file
  !------------------------------------------------------------------------
  SUBROUTINE parse_ctrl_file(info)
    INTEGER, INTENT(  OUT)      :: info
    CHARACTER(LEN=*), PARAMETER :: caller='parse_ctrl_file'
    CHARACTER(LEN=ppm_char)     :: cbuf, cvalue, carg, cvar
    INTEGER                     :: ilenctrl
    INTEGER                     :: iUnit, istat, ios
    LOGICAL                     :: lExist
    INTEGER                     :: iline
    INTEGER                     :: ilen
    INTEGER                     :: i,j,idx
    !-------------------------------------------------------------
    !  Check that the ctrl file exists
    !-------------------------------------------------------------
    ilenctrl = LEN_TRIM(ctrl_file_name)
    IF (ilenctrl .LT. 1) THEN
       WRITE (*,*) "No ctrl file given!"
       info = -1
       GOTO 9999
    END IF
    INQUIRE(FILE=ctrl_file_name, EXIST=lExist)
    IF (.NOT. lExist) THEN
       WRITE(*,'(2A)') 'No such file: ', ctrl_file_name(1:ilenctrl)
       info = -1
       GOTO 9999
    END IF
    !-------------------------------------------------------------
    !  Open the file
    !-------------------------------------------------------------
    iUnit = 20
    OPEN(iUnit, FILE=ctrl_file_name, IOSTAT=ios, ACTION='READ')
    IF (ios .NE. 0) THEN
       WRITE(*,'(2A)') 'Failed to open file: ', ctrl_file_name(1:ilenctrl)
       info = -1
       GOTO 9999
    END IF
    !-------------------------------------------------------------
    !  Scan file
    !-------------------------------------------------------------
    iline = 0
    var_loop: DO
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
             IF (cbuf(i:i) .NE. ' ') THEN
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
             WRITE(*, '(A,I5)') 'Incorrect line: ', iline
             info = -1
             GOTO 9999
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
#undef DTYPE
#define DTYPE REAL
#include "ctrl/parse_ctrl.f"
#undef DTYPE
#define DTYPE CHAR
#include "ctrl/parse_ctrl.f"
#undef DTYPE
#define DTYPE LOGICAL
#define BOOL
#include "ctrl/parse_ctrl.f"
#undef BOOL
#undef DTYPE
          ! array
#define ARRAY
#define DTYPE INTEGER
#include "ctrl/parse_ctrl.f"
#undef DTYPE
#define DTYPE REAL
#include "ctrl/parse_ctrl.f"
#undef DTYPE
#define DTYPE CHAR
#include "ctrl/parse_ctrl.f"
#undef DTYPE
#define DTYPE LOGICAL
#define BOOL
#include "ctrl/parse_ctrl.f"
#undef BOOL
#undef DTYPE          
#undef ARRAY
       END IF
    END DO var_loop
200 CONTINUE
    WRITE(*,'(A,I5,2A)') 'Error reading line: ', iline,     &
         &               ' of file: ', ctrl_file_name(1:ilenctrl)
    info = -1
    GOTO 9999
100 CONTINUE
    !----------------------------------------------------------------------
    !  Close file
    !----------------------------------------------------------------------
    CLOSE(iUnit)
9999 CONTINUE
    RETURN
  END SUBROUTINE parse_ctrl_file
  !------------------------------------------------------------------------
  !  Call default funcs
  !------------------------------------------------------------------------
  SUBROUTINE call_default_funcs(info)
    INTEGER, INTENT(  OUT) :: info
    INTEGER  :: i
    DO i=1,INTEGER_args_i
       IF (ASSOCIATED(INTEGER_args(i)%default_func)) THEN
          IF (.NOT. INTEGER_args(i)%default_func(INTEGER_args(i)%variable)) THEN
             info = 1
             GOTO 9999
          END IF
       END IF
    END DO
    DO i=1,REAL_args_i
       IF (ASSOCIATED(REAL_args(i)%default_func)) THEN
          IF (.NOT. REAL_args(i)%default_func(REAL_args(i)%variable)) THEN
             info = 1
             GOTO 9999
          END IF
       END IF
    END DO
    DO i=1,CHAR_args_i
       IF (ASSOCIATED(CHAR_args(i)%default_func)) THEN
          IF (.NOT. CHAR_args(i)%default_func(CHAR_args(i)%variable)) THEN
             info = 1
             GOTO 9999
          END IF
       END IF
    END DO
    DO i=1,LOGICAL_args_i
       IF (ASSOCIATED(LOGICAL_args(i)%default_func)) THEN
           IF (.NOT. LOGICAL_args(i)%default_func(LOGICAL_args(i)%variable)) THEN
              info = 1
              GOTO 9999
           END IF
       END IF
    END DO
    DO i=1,INTEGER_array_args_i
       IF (ASSOCIATED(INTEGER_array_args(i)%default_func)) THEN
          IF (.NOT. INTEGER_array_args(i)%default_func(INTEGER_array_args(i)%variable)) THEN
             info = 1
             GOTO 9999
          END IF
       END IF
    END DO
    DO i=1,REAL_array_args_i
       IF (ASSOCIATED(REAL_array_args(i)%default_func)) THEN
          IF (.NOT. REAL_array_args(i)%default_func(REAL_array_args(i)%variable)) THEN
             info = 1
             GOTO 9999
          END IF
       END IF
    END DO
    DO i=1,CHAR_array_args_i
       IF (ASSOCIATED(CHAR_array_args(i)%default_func)) THEN
          IF (.NOT. CHAR_array_args(i)%default_func(CHAR_array_args(i)%variable)) THEN
             info = 1
             GOTO 9999
          END IF
       END IF
    END DO
    DO i=1,LOGICAL_array_args_i
       IF (ASSOCIATED(LOGICAL_array_args(i)%default_func)) THEN
           IF (.NOT. LOGICAL_array_args(i)%default_func(LOGICAL_array_args(i)%variable)) THEN
              info = 1
              GOTO 9999
           END IF
       END IF
    END DO
9999 CONTINUE
  END SUBROUTINE call_default_funcs
  !------------------------------------------------------------------------
  !  Check min max
  !------------------------------------------------------------------------
  SUBROUTINE check_minmax(info)
    INTEGER, INTENT(  OUT) :: info
    INTEGER  :: i
    DO i=1,INTEGER_args_i
       IF (INTEGER_args(i)%variable .GT. INTEGER_args(i)%max .OR. &
            INTEGER_args(i)%variable .LT. INTEGER_args(i)%min) THEN
          WRITE (*,*) 'Argument ', INTEGER_args(i)&
               &%name(1:LEN_TRIM(INTEGER_args(i)%name)), &
               ' fails min max check!'
          info = 1
          GOTO 9999
       END IF
    END DO
    DO i=1,REAL_args_i
       IF (REAL_args(i)%variable .GT. REAL_args(i)%max .OR. &
               REAL_args(i)%variable .LT. REAL_args(i)%min) THEN
          WRITE (*,*) 'Argument ', REAL_args(i)%name(1:LEN_TRIM(REAL_args(i)%name)), &
               ' fails min max check!'
          info = 1
          GOTO 9999
       END IF
    END DO
    DO i=1,INTEGER_array_args_i
       IF (ANY(INTEGER_array_args(i)%variable .GT. INTEGER_array_args(i)%max) .OR. &
           ANY(INTEGER_array_args(i)%variable .LT. INTEGER_array_args(i)%min)) THEN
          WRITE (*,*) 'Argument ', INTEGER_array_args(i)&
               &%name(1:LEN_TRIM(INTEGER_array_args(i)%name)), &
               ' fails min max check!'
          info = 1
          GOTO 9999
       END IF
    END DO
    DO i=1,REAL_array_args_i
       IF (ANY(REAL_array_args(i)%variable .GT. REAL_array_args(i)%max) .OR. &
           ANY(REAL_array_args(i)%variable .LT. REAL_array_args(i)%min)) THEN
          WRITE (*,*) 'Argument ', REAL_array_args(i)%name(1:LEN_TRIM(REAL_array_args(i)%name)), &
               ' fails min max check!'
          info = 1
          GOTO 9999
       END IF
    END DO
9999 CONTINUE
  END SUBROUTINE check_minmax
  !------------------------------------------------------------------------
  !  Call validator functions
  !------------------------------------------------------------------------
  SUBROUTINE call_validator_funcs(info)
    INTEGER, INTENT(  OUT) :: info
    INTEGER  :: i
    ! scalar
#define DTYPE INTEGER
#include "ctrl/validate.f"
#undef DTYPE
#define DTYPE REAL
#include "ctrl/validate.f"
#undef DTYPE
#define DTYPE CHAR
#include "ctrl/validate.f"
#undef DTYPE
#define DTYPE LOGICAL
#include "ctrl/validate.f"
#undef DTYPE
    ! array
#define ARRAY
#define DTYPE INTEGER_array
#include "ctrl/validate.f"
#undef DTYPE
#define DTYPE REAL_array
#include "ctrl/validate.f"
#undef DTYPE
#define DTYPE CHAR_array
#include "ctrl/validate.f"
#undef DTYPE
#define DTYPE LOGICAL_array
#include "ctrl/validate.f"
#undef DTYPE
#undef ARRAY
9999 CONTINUE
  END SUBROUTINE call_validator_funcs
  !------------------------------------------------------------------------
  !  Special args
  !------------------------------------------------------------------------
  SUBROUTINE disable_help
    help_enabled = .FALSE.
  END SUBROUTINE disable_help

  SUBROUTINE disable_ctrl
    ctrl_enabled = .FALSE.
  END SUBROUTINE disable_ctrl

  SUBROUTINE set_ctrl_name(name)
    CHARACTER(LEN=*), INTENT(IN   )  :: name
    ctrl_file_name = name
  END SUBROUTINE set_ctrl_name
  !------------------------------------------------------------------------
  !  Set command line args (usefull for running tests)
  !------------------------------------------------------------------------
  SUBROUTINE add_cmd(arg, value)
    CHARACTER(LEN=*), INTENT(IN   )           :: arg
    CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL :: value
    CHARACTER(LEN=256), POINTER, DIMENSION(:) :: temp_a => NULL()
    INTEGER,            POINTER, DIMENSION(:) :: temp_l => NULL()
    INTEGER                                   :: inc
    INTEGER                                   :: l
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
  SUBROUTINE find_flag(name, success, value)
    CHARACTER(LEN=*), INTENT(IN   )             :: name
    LOGICAL,          INTENT(  OUT)             :: success
    CHARACTER(LEN=*), INTENT(  OUT), OPTIONAL   :: value
    INTEGER                                     :: i
    success = .FALSE.
    DO i=1,cmd_args_i
       IF (cmd_args(i)(1:cmd_args_len(i)) .EQ. name) THEN
          success          = .TRUE.
          IF (.NOT. cmd_args_used(i)) cmd_i = cmd_i - 1
          cmd_args_used(i) = .TRUE.
          IF (PRESENT(value)) THEN
             IF (i+1 .GT. cmd_args_i) THEN
                success = .FALSE.
                WRITE (*,*) "Flag ", name, " requires an argument."
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
    arg_count = cmd_i
  END FUNCTION arg_count
  !------------------------------------------------------------------------
  !  Look-up args (after all flags have been read!!!)
  !------------------------------------------------------------------------
  SUBROUTINE find_arg(position, success, value)
    INTEGER,          INTENT(IN   )  :: position
    CHARACTER(LEN=*), INTENT(INOUT)  :: value
    LOGICAL,          INTENT(  OUT)  :: success
    INTEGER                          :: i, j
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
    CHARACTER(LEN=*), INTENT(IN   )            :: name
    CHARACTER(LEN=256), POINTER, DIMENSION(:)  :: temp_g
    INTEGER,            POINTER, DIMENSION(:)  :: temp_s
    LOGICAL,            POINTER, DIMENSION(:)  :: temp_c
    LOGICAL,            POINTER, DIMENSION(:)  :: temp_a
    groups_i = groups_i + 1
    IF (.NOT. ASSOCIATED(groups)) THEN
       ALLOCATE(groups(0:di))
       ALLOCATE(group_size(0:di))
       ALLOCATE(group_has_ctrl(0:di))
       ALLOCATE(group_has_arg(0:di))
    ELSE IF (groups_i+1 .GT. SIZE(groups)) THEN
       ALLOCATE(temp_g(0:(SIZE(groups)+di)))
       ALLOCATE(temp_s(0:(SIZE(groups)+di)))
       ALLOCATE(temp_c(0:(SIZE(groups)+di)))
       ALLOCATE(temp_a(0:(SIZE(groups)+di)))
       temp_g(0:SIZE(groups)) = groups(0:SIZE(groups))
       temp_s(0:SIZE(groups)) = group_size(0:SIZE(groups))
       temp_c(0:SIZE(groups)) = group_has_ctrl(0:SIZE(groups))
       temp_a(0:SIZE(groups)) = group_has_arg(0:SIZE(groups))
       DEALLOCATE(groups)
       DEALLOCATE(group_size)
       DEALLOCATE(group_has_ctrl)
       DEALLOCATE(group_has_arg)
       groups         => temp_g
       group_size     => temp_s
       group_has_ctrl => temp_c
       group_has_arg  => temp_a
    END IF
    groups(groups_i)         = name
    group_size(groups_i)     = 0
    group_has_ctrl(groups_i) = .false.
    group_has_arg(groups_i)  = .false.
  END SUBROUTINE arg_group
  !------------------------------------------------------------------------
  !  Print command line help
  !------------------------------------------------------------------------
  SUBROUTINE print_help
    INTEGER  :: i, j, k, l
    WRITE (*,'(A)') "Usage: progname {ctrl-file} [options]"
    DO k=0,groups_i
       WRITE (*,'(/A)') groups(k)
       IF (k .EQ. 0 .AND. i .EQ. 1) THEN
          IF (help_enabled) THEN
             WRITE (*,'(/A)') '   help                          Print &
                  &this help message and exit.'
             WRITE (*,*) '                  short flag :  -h'
             WRITE (*,*) '                   long flag :  --help'
          END IF
          IF (ctrl_enabled) THEN
             WRITE (*,'(/A)') '   control file                  Print &
                  &sample control file and exit.'
             WRITE (*,*) '                   long flag :  --print-ctrl'
          END IF  
       END IF
       group_loop: DO i=1,group_size(k)
          ! scalar
#define DTYPE INTEGER
#include "ctrl/help.f"
#undef DTYPE
#define DTYPE REAL
#include "ctrl/help.f"
#undef DTYPE
#define DTYPE LOGICAL
#define BOOL
#include "ctrl/help.f"
#undef BOOL
#undef DTYPE
#define DTYPE CHAR
#define STRING
#include "ctrl/help.f"
#undef STRING
#undef DTYPE
          ! array
#define ARRAY
#define DTYPE INTEGER_array
#include "ctrl/help.f"
#undef DTYPE
#define DTYPE REAL_array
#include "ctrl/help.f"
#undef DTYPE
#define DTYPE LOGICAL_array
#define BOOL
#include "ctrl/help.f"
#undef BOOL
#undef DTYPE
#define DTYPE CHAR_array
#define STRING
#include "ctrl/help.f"
#undef STRING
#undef DTYPE
#undef ARRAY
       END DO group_loop
    END DO
  END SUBROUTINE print_help
  !------------------------------------------------------------------------
  !  Print sample control file
  !------------------------------------------------------------------------
  SUBROUTINE print_ctrl
    INTEGER  :: i, j, k, l
    IF (groups_i .EQ. -1) THEN
       WRITE (*,*) "No args have been defined..."
       RETURN
    END IF
    WRITE (*,'(A)') "#--------------------------------------------------------------------------"
    WRITE (*,'(A)') "#  Sample control file for ppm_client"
    WRITE (*,'(A)') "#"
    WRITE (*,'(A)') "#  Edit the settings below"
    WRITE (*,'(A)') "#"
    WRITE (*,'(A)') "#--------------------------------------------------------------------------"
    DO k=0,groups_i
       WRITE (*,'(/A)') "#--------------------------------------------------------------------------"
       WRITE (*,'(2A)') "#  ", groups(k)(1:LEN_TRIM(groups(k)))

       comment_loop: DO i=1,group_size(k)
          ! scalar
#define DTYPE INTEGER
#include "ctrl/ctrl_comment.f"
#undef DTYPE
#define DTYPE REAL
#include "ctrl/ctrl_comment.f"
#undef DTYPE
#define DTYPE LOGICAL
#define BOOL
#include "ctrl/ctrl_comment.f"
#undef BOOL
#undef DTYPE
#define DTYPE CHAR
#define STRING
#include "ctrl/ctrl_comment.f"
#undef STRING
#undef DTYPE
          ! array
#define ARRAY
#define DTYPE INTEGER_array
#include "ctrl/ctrl_comment.f"
#undef DTYPE
#define DTYPE REAL_array
#include "ctrl/ctrl_comment.f"
#undef DTYPE
#define DTYPE LOGICAL_array
#define BOOL
#include "ctrl/ctrl_comment.f"
#undef BOOL
#undef DTYPE
#define DTYPE CHAR_array
#define STRING
#include "ctrl/ctrl_comment.f"
#undef STRING
#undef DTYPE
#undef ARRAY
       END DO comment_loop

       WRITE (*,'(A)') "#--------------------------------------------------------------------------"

       var_loop: DO i=1,group_size(k)
          ! scalar
#define DTYPE INTEGER
#include "ctrl/ctrl.f"
#undef DTYPE
#define DTYPE REAL
#include "ctrl/ctrl.f"
#undef DTYPE
#define DTYPE LOGICAL
#define BOOL
#include "ctrl/ctrl.f"
#undef BOOL
#undef DTYPE
#define DTYPE CHAR
#define STRING
#include "ctrl/ctrl.f"
#undef STRING
#undef DTYPE
          ! array
#define ARRAY
#define DTYPE INTEGER_array
#include "ctrl/ctrl.f"
#undef DTYPE
#define DTYPE REAL_array
#include "ctrl/ctrl.f"
#undef DTYPE
#define DTYPE LOGICAL_array
#define BOOL
#include "ctrl/ctrl.f"
#undef BOOL
#undef DTYPE
#define DTYPE CHAR_array
#define STRING
#include "ctrl/ctrl.f"
#undef STRING
#undef DTYPE
#undef ARRAY
       END DO var_loop

    END DO
  END SUBROUTINE print_ctrl
  !------------------------------------------------------------------------
  !  Debug
  !------------------------------------------------------------------------
  SUBROUTINE dump_defines
    INTEGER :: i
    WRITE (*,'(//A)') "DUMPING DEFINED ARGS"
    WRITE (*,'(//A/)') "Integers: "
    DO i=1,INTEGER_args_i
#define DTYPE INTEGER
#include "ctrl/dump.f"
#undef DTYPE
    END DO
    WRITE (*,'(//A/)') "Reals: "
    DO i=1,REAL_args_i
#define DTYPE REAL
#include "ctrl/dump.f"
#undef DTYPE
    END DO
    WRITE (*,'(//A/)') "Chars: "
    DO i=1,CHAR_args_i
#define DTYPE CHAR
#define STRING
#include "ctrl/dump.f"
#undef STRING
#undef DTYPE
    END DO
    WRITE (*,'(//A/)') "Logicals: "
    DO i=1,LOGICAL_args_i
#define DTYPE LOGICAL
#define BOOL
#include "ctrl/dump.f"
#undef BOOL
#undef DTYPE
    END DO
  END SUBROUTINE dump_defines
  SUBROUTINE dump_args
    INTEGER :: i
    WRITE (*,'(//A/)') 'DUMPING COMMAND ARGS'
    DO i=1,cmd_args_i
       WRITE (*,*) 'arg : ', cmd_args(i)(1:cmd_args_len(i))
       IF (cmd_args_used(i)) THEN
          WRITE (*,*) 'used'
       ELSE
          WRITE (*,*) 'not used'
       END IF
    END DO
  END SUBROUTINE dump_args
  !------------------------------------------------------------------------
  !  Adders - black or otherwise
  !------------------------------------------------------------------------

  ! scalar

#define DTYPE INTEGER
#include "ctrl/adder.f"
#undef DTYPE

#define DTYPE REAL
#include "ctrl/adder.f"
#undef DTYPE

#define DTYPE CHAR
#define STRING
#include "ctrl/adder.f"
#undef STRING
#undef DTYPE

#define DTYPE LOGICAL
#define BOOL
#include "ctrl/adder.f"
#undef BOOL
#undef DTYPE

  ! array

#define ARRAY

#define DTYPE INTEGER
#include "ctrl/adder.f"
#undef DTYPE

#define DTYPE REAL
#include "ctrl/adder.f"
#undef DTYPE

#define DTYPE CHAR
#define STRING
#include "ctrl/adder.f"
#undef STRING
#undef DTYPE

#define DTYPE LOGICAL
#define BOOL
#include "ctrl/adder.f"
#undef BOOL
#undef DTYPE

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

END MODULE ppm_module_ctrl
