#define WRAP(a) a

#ifdef ARRAY
  TYPE WRAP(DTYPE)_array_arg
#else
  TYPE WRAP(DTYPE)_arg
#endif

#ifdef ARRAY

#ifdef __INTEGER
     INTEGER,                  DIMENSION(:),     POINTER :: variable => NULL()
     INTEGER,                  DIMENSION(:), ALLOCATABLE :: default
#elif defined(__LONGINT)
     INTEGER(ppm_kind_int64),  DIMENSION(:),     POINTER :: variable => NULL()
     INTEGER(ppm_kind_int64),  DIMENSION(:), ALLOCATABLE :: default
#elif defined(__SINGLE)
     REAL(ppm_kind_single),    DIMENSION(:),     POINTER :: variable => NULL()
     REAL(ppm_kind_single),    DIMENSION(:), ALLOCATABLE :: default
#elif defined(__DOUBLE)
     REAL(ppm_kind_double),    DIMENSION(:),     POINTER :: variable => NULL()
     REAL(ppm_kind_double),    DIMENSION(:), ALLOCATABLE :: default
#elif defined(__LOGICAL)
     LOGICAL,                  DIMENSION(:),     POINTER :: variable => NULL()
     LOGICAL,                  DIMENSION(:), ALLOCATABLE :: default
#elif defined(__STRING)
     CHARACTER(LEN=ppm_char),  DIMENSION(:),     POINTER :: variable => NULL()
     CHARACTER(LEN=ppm_char),  DIMENSION(:), ALLOCATABLE :: default
#elif defined(__COMPLEX)
     COMPLEX(ppm_kind_single), DIMENSION(:),     POINTER :: variable => NULL()
     COMPLEX(ppm_kind_single), DIMENSION(:), ALLOCATABLE :: default
#elif defined(__DCOMPLEX)
     COMPLEX(ppm_kind_double), DIMENSION(:),     POINTER :: variable => NULL()
     COMPLEX(ppm_kind_double), DIMENSION(:), ALLOCATABLE :: default
#endif

#else

#ifdef __INTEGER
     INTEGER,                                    POINTER :: variable => NULL()
     INTEGER                                             :: default
#elif defined(__LONGINT)
     INTEGER(ppm_kind_int64),                    POINTER :: variable => NULL()
     INTEGER(ppm_kind_int64)                             :: default
#elif defined(__SINGLE)
     REAL(ppm_kind_single),                      POINTER :: variable => NULL()
     REAL(ppm_kind_single)                               :: default
#elif defined(__DOUBLE)
     REAL(ppm_kind_double),                      POINTER :: variable => NULL()
     REAL(ppm_kind_double)                               :: default
#elif defined(__LOGICAL)
     LOGICAL,                                    POINTER :: variable => NULL()
     LOGICAL                                             :: default
#elif defined(__STRING)
     CHARACTER(LEN=ppm_char),                    POINTER :: variable => NULL()
     CHARACTER(LEN=ppm_char)                             :: default
#elif defined(__COMPLEX)
     COMPLEX(ppm_kind_single),                   POINTER :: variable => NULL()
     COMPLEX(ppm_kind_single)                            :: default
#elif defined(__DCOMPLEX)
     COMPLEX(ppm_kind_double),                   POINTER :: variable => NULL()
     COMPLEX(ppm_kind_double)                            :: default
#endif

#endif

#ifdef __INTEGER
     INTEGER                                             :: min
     INTEGER                                             :: max
#elif defined(__LONGINT)
     INTEGER(ppm_kind_int64)                             :: min
     INTEGER(ppm_kind_int64)                             :: max
#elif defined(__SINGLE)
     REAL(ppm_kind_single)                               :: min
     REAL(ppm_kind_single)                               :: max
#elif defined(__DOUBLE)
     REAL(ppm_kind_double)                               :: min
     REAL(ppm_kind_double)                               :: max
#elif defined(__STRING)
!      INTEGER                                             :: min
!      INTEGER                                             :: max
#endif
     LOGICAL                                             :: min_set       = .FALSE.
     LOGICAL                                             :: max_set       = .FALSE.
#if defined(__LOGICAL) && !defined(ARRAY)
     LOGICAL                                             :: type          = .TRUE. ! enable
#endif
     LOGICAL                                             :: default_set   = .FALSE.
     CHARACTER(LEN=256)                                  :: name
     CHARACTER(LEN=2)                                    :: flag
     LOGICAL                                             :: flag_set      = .FALSE.
     CHARACTER(LEN=256)                                  :: long_flag
     LOGICAL                                             :: long_flag_set = .FALSE.
     CHARACTER(LEN=256)                                  :: ctrl_name
     LOGICAL                                             :: ctrl_name_set = .FALSE.
     LOGICAL                                             :: settable      = .FALSE.
     CHARACTER(LEN=ppm_char)                             :: help
     LOGICAL                                             :: help_set      = .FALSE.
     LOGICAL                                             :: clf_supplied  = .FALSE.
     INTEGER                                             :: group
     INTEGER                                             :: group_i
     LOGICAL                                             :: default_func_set = .FALSE.
     LOGICAL                                             :: validator_set = .FALSE.
#ifdef ARRAY
     PROCEDURE(WRAP(DTYPE)_array_func), POINTER, NOPASS :: default_func => NULL()
     PROCEDURE(WRAP(DTYPE)_array_func), POINTER, NOPASS :: validator    => NULL()     
  END TYPE WRAP(DTYPE)_array_arg
#else
     PROCEDURE(WRAP(DTYPE)_func), POINTER, NOPASS :: default_func => NULL()
     PROCEDURE(WRAP(DTYPE)_func), POINTER, NOPASS :: validator    => NULL()
  END TYPE WRAP(DTYPE)_arg
#endif
#undef __INTEGER
#undef __LONGINT
#undef __SINGLE
#undef __DOUBLE
#undef __LOGICAL
#undef __STRING
#undef __COMPLEX
#undef __DCOMPLEX
#undef DTYPE
