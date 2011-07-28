test_suite ppm_module_ctrl

  integer            :: idefault, iflag, ilong_flag
  real               :: rdefault, rflag, rlong_flag
  character(len=256) :: cdefault, cflag, clong_flag
  logical            :: ldefault, leflag, lelong_flag
  logical            ::           ldflag, ldlong_flag

  integer,            dimension(3) :: iarray, iaflag, ialflag
  real,               dimension(3) :: rarray, raflag, ralflag
  character(len=256), dimension(3) :: carray, caflag, calflag
  logical,            dimension(3) :: larray, laflag, lalflag

  character(len=256) :: value
  logical            :: ok

  integer :: info, ios

  setup
    CALL disable_ctrl
  !  CALL arg_group("Integer Tests")
  !    CALL arg(idefault, 'idefault',     &
  !         flag      = '-d',             &
  !         long_flag = '--idefault',     &
  !         ctrl_name = 'idefault',       &
  !         default   = 42,               &
  !         help      = "Test if integer defaults work.")
  end setup

  teardown
    CALL reset
    idefault    = 0
    iflag       = 0
    ilong_flag  = 0
    rdefault    = 0
    rflag       = 0
    rlong_flag  = 0
    cdefault    = ''
    cflag       = ''
    clong_flag  = ''
    ldefault    = .false.
    leflag      = .false.
    lelong_flag = .false.
    ldflag      = .false.
    ldlong_flag = .false.
    value       = ''
    ok          = .false.
    info        = 0
  end teardown

  test printing
    CALL reset
    ! define some args
    CALL disable_ctrl
    CALL arg(idefault, 'some_int',     &
         flag      = '-i',             &
         long_flag = '--integer',      &
         ctrl_name = 'some_int',       &
         default   = 42,               &
         min       = 0,                &
         max       = 100,              &
         help      = "Test if integer printing works.")
    CALL arg(rdefault, 'some_real',    &
         flag      = '-r',             &
         long_flag = '--real',         &
         ctrl_name = 'some_real',      &
         default   = 0.1337,           &
         min       = 0.0,              &
         max       = 1.0,              &
         help      = "Test if real printing works.")
    CALL arg_group("Second Group")
    CALL arg(cdefault, 'some_char',    &
         flag      = '-c',             &
         long_flag = '--char',         &
         ctrl_name = 'some_char',      &
         default   = 'hrkljus',        &
         help      = "Test if char printing works.")
    CALL arg(ldefault, 'some_bool',    &
         flag      = '-b',             &
         long_flag = '--bool',         &
         ctrl_name = 'some_bool',      &
         default   = .true.,           &
         help      = "Test if bool printing works.")
    CALL arg(iarray, 'int_array',      &
         flag      = '-m',             &
         long_flag = '--int-array',    &
         ctrl_name = 'int_array',      &
         default   = (/1,2,3/),        &
         min       = 0,                &
         max       = 100,              &
         help      = "Test if integer array printing works.")
    CALL arg(rarray, 'real_array',     &
         flag      = '-n',             &
         long_flag = '--real-array',   &
         ctrl_name = 'real_array',     &
         default   = (/1.1,2.2,3.3/),  &
         min       = 0.0,              &
         max       = 3.5,              &
         help      = "Test if real array printing works.")
    CALL arg_group("Third Group")
    CALL arg(carray, 'char_array',     &
         flag      = '-o',             &
         long_flag = '--char-array',   &
         ctrl_name = 'char_array',     &
         default   = (/'foo','bar','baz'/), &
         help      = "Test if char array printing works.")
    CALL arg(larray, 'bool_array',     &
         flag      = '-p',             &
         long_flag = '--bool-array',   &
         ctrl_name = 'bool_array',     &
         default   = (/.true.,.false.,.true./), &
         help      = "Test if bool array printing works.")
    ! supply flag
!    CALL add_cmd('-h')
!    CALL add_cmd('--help')
!    CALL add_cmd('--print-ctrl')
    ! parse args
    CALL parse_args(info)
    Assert_Equal(info, 0)
  end test

  test arg_manipulation
    ! add args
    ! plain cmd arg
    CALL add_cmd('start')
    ! short flag
    ! with value
    CALL add_cmd('-a', '1')
    CALL add_cmd('-b', '2')
    ! plain cmd arg
    CALL add_cmd('intershort')
    CALL add_cmd('-c', '3')
    CALL add_cmd('-d', '4')
    ! no value
    CALL add_cmd('-e')
    CALL add_cmd('-f')
    CALL add_cmd('-g')
    ! plain cmd arg
    CALL add_cmd('middle')
    ! long flag
    ! with value
    CALL add_cmd('--long-1', '1')
    CALL add_cmd('--long-2', '2')
    ! plain cmd arg
    CALL add_cmd('interlong')
    CALL add_cmd('--long-3', '3')
    ! no value
    CALL add_cmd('--long-4')
    CALL add_cmd('--long-5')
    ! plain cmd arg
    CALL add_cmd('finalize')
    CALL add_cmd('stop')

    ! to see that args work after info
    CALL parse_args(info)

    ! check flag function (out of order but before arg)
    CALL find_flag('-b', ok, value)
    Assert_True(ok)
    Assert_Equal(value, '2')
    CALL find_flag('--long-3', ok, value)
    Assert_True(ok)
    Assert_Equal(value, '3')
    CALL find_flag('-e', ok)
    Assert_True(ok)
    CALL find_flag('-c', ok, value)
    Assert_True(ok)
    Assert_Equal(value, '3')
    CALL find_flag('--long-1', ok, value)
    Assert_True(ok)
    Assert_Equal(value, '1')
    CALL find_flag('--long-5', ok)
    Assert_True(ok)
    CALL find_flag('-a', ok, value)
    Assert_True(ok)
    Assert_Equal(value, '1')
    CALL find_flag('-d', ok, value)
    Assert_True(ok)
    Assert_Equal(value, '4')
    CALL find_flag('-f', ok)
    Assert_True(ok)
    CALL find_flag('--long-2', ok, value)
    Assert_True(ok)
    Assert_Equal(value, '2')
    CALL find_flag('-g', ok)
    Assert_True(ok)
    CALL find_flag('--long-4', ok)
    Assert_True(ok)

    ! check number of args
    Assert_Equal(6, arg_count())

    ! check arg function
    CALL find_arg(3, ok, value)
    Assert_True(ok)
    Assert_Equal(value, 'middle')
    CALL find_arg(2, ok, value)
    Assert_True(ok)
    Assert_Equal(value, 'intershort')
    CALL find_arg(6, ok, value)
    Assert_True(ok)
    Assert_Equal(value, 'stop')
    CALL find_arg(4, ok, value)
    Assert_True(ok)
    Assert_Equal(value, 'interlong')
    CALL find_arg(1, ok, value)
    Assert_True(ok)
    Assert_Equal(value, 'start')
    CALL find_arg(5, ok, value)
    Assert_True(ok)
    Assert_Equal(value, 'finalize')
    
  end test

  test integer_args
    ! define
    CALL arg(idefault,   'idefault',   default   = 42      )
    CALL arg(iflag,      'iflag',      flag      = '-f'    )
    CALL arg(ilong_flag, 'ilong_flag', long_flag = '--flag')
    ! supply
    CALL add_cmd('-f',     '42')
    CALL add_cmd('--flag', '42')
    ! parse
    CALL parse_args(info)
    ! test
    Assert_Equal(idefault,   42)
    Assert_Equal(iflag,      42)
    Assert_Equal(ilong_flag, 42)
  end test

  test integer_array_args
    ! define
    CALL arg(iarray,  'iarray',  default   = (/1,2,3/))
    CALL arg(iaflag,  'iaflag',  flag      = '-f'     )
    CALL arg(ialflag, 'ialflag', long_flag = '--flag' )
    ! supply
    CALL add_cmd('-f',     '1,2,3'  )
    CALL add_cmd('--flag', '1, 2, 3')
    ! parse
    CALL parse_args(info)
    ! test
    Assert_Array_Equal(iarray,  (/1, 2, 3/))
    Assert_Array_Equal(iaflag,  (/1, 2, 3/))
    Assert_Array_Equal(ialflag, (/1, 2, 3/))
  end test

  test real_args
    ! define
    CALL arg(rdefault,   'rdefault',   default   = 0.1337  )
    CALL arg(rflag,      'rflag',      flag      = '-f'    )
    CALL arg(rlong_flag, 'rlong_flag', long_flag = '--flag')
    ! supply
    CALL add_cmd('-f',     '0.1337')
    CALL add_cmd('--flag', '0.1337')
    ! parse
    CALL parse_args(info)
    ! test
    Assert_Equal(rdefault,   0.1337)
    Assert_Equal(rflag,      0.1337)
    Assert_Equal(rlong_flag, 0.1337)
  end test

  test real_array_args
    ! define
    CALL arg(rarray,  'rarray',  default   = (/1.1,2.2,3.3/))
    CALL arg(raflag,  'raflag',  flag      = '-f'     )
    CALL arg(ralflag, 'ralflag', long_flag = '--flag' )
    ! supply
    CALL add_cmd('-f',     '1.1,2.2,3.3'  )
    CALL add_cmd('--flag', '1.1, 2.2, 3.3')
    ! parse
    CALL parse_args(info)
    ! test
    Assert_Array_Equal(rarray,  (/1.1, 2.2, 3.3/))
    Assert_Array_Equal(raflag,  (/1.1, 2.2, 3.3/))
    Assert_Array_Equal(ralflag, (/1.1, 2.2, 3.3/))
  end test

  test char_args
    CALL arg(cdefault,   'cdefault',   default   = 'hrkljus')
    CALL arg(cflag,      'cflag',      flag      = '-f')
    CALL arg(clong_flag, 'clong_flag', long_flag = '--flag')
    ! supply
    CALL add_cmd('-f',     'hrkljus')
    CALL add_cmd('--flag', 'hrkljus')
    ! parse
    CALL parse_args(info)
    ! test
    Assert_Equal(cdefault,   'hrkljus')
    Assert_Equal(cflag,      'hrkljus')
    Assert_Equal(clong_flag, 'hrkljus')
  end test

  test char_array_args
    ! define
    CALL arg(carray,  'carray',  default   = (/'foo','bar','baz'/))
    CALL arg(caflag,  'caflag',  flag      = '-f'     )
    CALL arg(calflag, 'calflag', long_flag = '--flag' )
    ! supply
    CALL add_cmd('-f',     'foo,bar,baz'  )
    CALL add_cmd('--flag', 'foo, bar, baz')
    ! parse
    CALL parse_args(info)
    ! test
    Assert_Array_Equal(carray,  (/'foo','bar','baz'/))
    Assert_Array_Equal(caflag,  (/'foo','bar','baz'/))
    Assert_Array_Equal(calflag, (/'foo','bar','baz'/))
  end test

  test logical_args
    CALL arg(ldefault,    'ldefault',    default   = .true.)
    CALL arg(leflag,      'leflag',      default   = .false., &
         type = enabling_flag, flag      = '-e')
    CALL arg(lelong_flag, 'lelong_flag', default   = .false., &
         type = enabling_flag, long_flag = '--enable')
    CALL arg(ldflag,      'ldflag',      default   = .true.,  &
         type = disabling_flag, flag      = '-d')
    CALL arg(ldlong_flag, 'ldlong_flag', default   = .true.,  &
         type = disabling_flag, long_flag = '--disable')
    ! supply
    CALL add_cmd('-e'       )
    CALL add_cmd('--enable' )
    CALL add_cmd('-d'       )
    CALL add_cmd('--disable')
    ! parse
    CALL parse_args(info)
    ! test
    Assert_True( ldefault   )
    Assert_True( leflag     )
    Assert_True( lelong_flag)
    Assert_False(ldflag     )
    Assert_False(ldlong_flag)
  end test

  test logical_array_args
    ! define
    CALL arg(larray,  'larray',  default   = (/.true.,.false.,.true./))
    CALL arg(laflag,  'laflag',  flag      = '-f'     )
    CALL arg(lalflag, 'lalflag', long_flag = '--flag' )
    ! supply
    CALL add_cmd('-f',     'T,F,T'  )
    CALL add_cmd('--flag', 'T, F, T')
    ! parse
    CALL parse_args(info)
    ! test
    Assert_LArray_Equal(larray,  (/.true.,.false.,.true./))
    Assert_LArray_Equal(laflag,  (/.true.,.false.,.true./))
    Assert_LArray_Equal(lalflag, (/.true.,.false.,.true./))
  end test

  test ctrl_file
    ! auto enables ctrl_file
    CALL reset
    ! define
    CALL arg(idefault, 'idefault', ctrl_name = 'idefault')
    CALL arg(rdefault, 'rdefault', ctrl_name = 'rdefault')
    CALL arg(cdefault, 'cdefault', ctrl_name = 'cdefault')
    CALL arg(ldefault, 'ldefault', ctrl_name = 'ldefault')
    CALL arg(iarray,   'iarray',   ctrl_name = 'iarray')
    CALL arg(rarray,   'rarray',   ctrl_name = 'rarray')
    CALL arg(carray,   'carray',   ctrl_name = 'carray')
    CALL arg(larray,   'larray',   ctrl_name = 'larray')
    ! open file for writing
    OPEN(20, FILE='src/ctrl/test/__test_ctrl', IOSTAT=ios, ACTION='WRITE')
    Assert_Equal(ios, 0)
    WRITE(20,'(A)') 'idefault = 42'
    WRITE(20,'(A)') 'rdefault = 0.1337'
    WRITE(20,'(A)') 'cdefault = hrkljus'
    WRITE(20,'(A)') 'ldefault = T'
    WRITE(20,'(A)') 'iarray   = 1    2    3'
    WRITE(20,'(A)') 'rarray   = 1.1  2.2  3.3'
    WRITE(20,'(A)') 'carray   = foo, bar, baz'
    WRITE(20,'(A)') 'larray   = T    F    T'
    CLOSE(20)
    ! supply arg
    CALL add_cmd('src/ctrl/test/__test_ctrl')
    ! parse
    CALL parse_args(info)
    ! test
    Assert_Equal(idefault,   42)
    Assert_Equal(rdefault,   0.1337)
    Assert_Equal(cdefault,   'hrkljus')
    Assert_True(ldefault)
    Assert_Array_Equal(iarray, (/1, 2, 3/))
    Assert_Array_Equal(rarray, (/1.1, 2.2, 3.3/))
    Assert_Array_Equal(carray, (/'foo', 'bar', 'baz'/))
    Assert_LArray_Equal(larray, (/.true., .false., .true./))
    ! cleanup
    CALL SYSTEM('/bin/rm src/ctrl/test/__test_ctrl')
  end test

  test default_funcs
    ! procedure declarations
    procedure(INTEGER_dflt)       int_def
    procedure(REAL_dflt)          flt_def
    procedure(CHAR_dflt)          chr_def
    procedure(LOGICAL_dflt)       log_def
    procedure(INTEGER_array_dflt) inta_def
    procedure(REAL_array_dflt)    flta_def
    procedure(CHAR_array_dflt)    chra_def
    procedure(LOGICAL_array_dflt) loga_def
    ! define
    CALL arg(idefault, 'idefault', default_func = int_def)
    CALL arg(rdefault, 'rdefault', default_func = flt_def)
    CALL arg(cdefault, 'cdefault', default_func = chr_def)
    CALL arg(ldefault, 'ldefault', default_func = log_def)
    CALL arg(iarray,   'iarray',   default_func = inta_def)
    CALL arg(rarray,   'rarray',   default_func = flta_def)
    CALL arg(carray,   'carray',   default_func = chra_def)
    CALL arg(larray,   'larray',   default_func = loga_def)
    ! parse
    CALL parse_args(info)
    ! test
    Assert_Equal(idefault, 42       )
    Assert_Equal(rdefault, 0.1337   )
    Assert_Equal(cdefault, "hrkljus")
    Assert_True(ldefault)
    Assert_Array_Equal( iarray, (/1, 2, 3/))
    Assert_Array_Equal( rarray, (/1.1, 2.2, 3.3/))
    Assert_Array_Equal( carray, (/'foo', 'bar', 'baz'/))
    Assert_LArray_Equal(larray, (/.true., .false., .true./))
  end test

end test_suite

logical function int_def(var)
  integer, pointer, intent(in) :: var
  var = 42
  int_def = .true.
end function int_def

logical function flt_def(var)
  real, pointer, intent(in) :: var
  var = 0.1337
  flt_def = .true.
end function flt_def

logical function chr_def(var)
  character(len=*), pointer, intent(in) :: var
  var = "hrkljus"
  chr_def = .true.
end function chr_def

logical function log_def(var)
  logical, pointer, intent(in) :: var
  var = .true.
  log_def = .true.
end function log_def

logical function inta_def(var)
  integer, dimension(:), pointer, intent(in) :: var
  var = (/1,2,3/)
  inta_def = .true.
end function inta_def

logical function flta_def(var)
  real, dimension(:), pointer, intent(in) :: var
  var = (/1.1,2.2,3.3/)
  flta_def = .true.
end function flta_def

logical function chra_def(var)
  character(len=*), dimension(:), pointer, intent(in) :: var
  var = (/'foo','bar','baz'/)
  chra_def = .true.
end function chra_def

logical function loga_def(var)
  logical, dimension(:), pointer, intent(in) :: var
  var = (/.true.,.false.,.true./)
  loga_def = .true.
end function loga_def
