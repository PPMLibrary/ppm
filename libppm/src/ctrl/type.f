#define WRAP(a) a
#ifdef ARRAY
  TYPE WRAP(DTYPE)_array_arg
#else
  TYPE WRAP(DTYPE)_arg
#endif
#ifdef STRING
#ifdef ARRAY
     CHARACTER(LEN=ppm_char), DIMENSION(:), POINTER     :: variable => NULL()
     CHARACTER(LEN=ppm_char), DIMENSION(:), ALLOCATABLE :: default
#else
     CHARACTER(LEN=ppm_char), POINTER             :: variable => NULL()
     CHARACTER(LEN=ppm_char)                      :: default
#endif
#else
#ifdef ARRAY
     DTYPE, DIMENSION(:),     POINTER             :: variable => NULL()
     DTYPE, DIMENSION(:), ALLOCATABLE             :: default
#else
     DTYPE, POINTER                               :: variable => NULL()
     DTYPE                                        :: default
#endif
#ifndef BOOL
     DTYPE                                        :: min
     LOGICAL                                      :: min_set       = .FALSE.
     DTYPE                                        :: max
     LOGICAL                                      :: max_set       = .FALSE.
#endif
#endif
#ifdef BOOL
#ifndef ARRAY
     LOGICAL                                      :: type = .true. ! enable
#endif
#endif
     LOGICAL                                      :: default_set   = .FALSE.
     CHARACTER(LEN=256)                           :: name
     CHARACTER(LEN=2)                             :: flag
     LOGICAL                                      :: flag_set      = .FALSE.
     CHARACTER(LEN=256)                           :: long_flag
     LOGICAL                                      :: long_flag_set = .FALSE.
     CHARACTER(LEN=256)                           :: ctrl_name
     LOGICAL                                      :: ctrl_name_set = .FALSE.
     CHARACTER(LEN=ppm_char)                      :: help
     LOGICAL                                      :: help_set      = .FALSE.
     INTEGER                                      :: group
     INTEGER                                      :: group_i
#ifdef ARRAY
     PROCEDURE(WRAP(DTYPE)_array_dflt), POINTER, NOPASS :: default_func => NULL()
     PROCEDURE(WRAP(DTYPE)_array_vdtr), POINTER, NOPASS :: validator    => NULL()     
  END TYPE WRAP(DTYPE)_array_arg
#else
     PROCEDURE(WRAP(DTYPE)_dflt), POINTER, NOPASS :: default_func => NULL()
     PROCEDURE(WRAP(DTYPE)_vdtr), POINTER, NOPASS :: validator    => NULL()
  END TYPE WRAP(DTYPE)_arg
#endif
