#-------------------------------------------------------------------------
#  File         :  configure.ac
#-------------------------------------------------------------------------
#
#  Purpose      :  Compilation
#
#  Remarks      :
#
#  References   :
#
#  Revisions    : This autoconf was written by Justin Park  8.17.2009
#                 Revised by Omar Awile
#-------------------------------------------------------------------------
#  Parallel Particle Mesh Library (PPM)
#  ETH Zurich
#  CH-8092 Zurich, Switzerland
#-------------------------------------------------------------------------

# Initialization
AC_PREREQ(2.61)
AC_INIT(PPM, 1.2.1)

# Reset existing preprocessing flags
#touch ./src/ppm_define.h
#rm ./src/ppm_define.h
#touch ./src/ppm_define.h
AC_SUBST(CPPFLAGS,[""])


# if MPI is enabled try grabbing those compiler wrappers first
use_mpi=no 
AC_ARG_ENABLE([mpi],
  [AS_HELP_STRING([--enable-mpi[[=impl.]]],
  [use MPI (default is no), If the MPI implementation of your choice provides
   compile wrappers that are in PATH, I can set them myself,  
   choose: guess (I will choose the first implementation I can find), openmpi, 
   lammpi, mpich, mpich2, intelmpi_gnu (Intel's MPI with GNU
   Compilers), intelmpi_intel (Intel's MPI with Intel Compilers), sun (Sun MPI), ibm
   (IBM AIX POE). Else, set this flag with no value and set CC, CXX and FC to
   the appropriate compiler wrappers (safest)])],
  [PPM_CHECK_MPI([$enableval])
    # compilers should be chosen
    use_mpi=yes],
  [AC_MSG_NOTICE([MPI was not enabled. Building non-parallel version of PPM.])])

: ${CXXFLAGS="-O3"}
: ${CFLAGS="-O3"}
: ${FCFLAGS="-O3"}

# Checks for basic compiler and flags
AC_PROG_CXX
AC_PROG_CPP
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_FC
AC_FC_FREEFORM
AC_FC_LIBRARY_LDFLAGS
# Check architecture
AC_CHECK_SIZEOF([long double])

AX_F90_MODULE_FLAG
AC_SUBST([ax_cv_f90_modflag])

# Check for Fortran2003 support
AC_MSG_CHECKING([whether your compiler properly supports Fortran 2003])
AC_LANG([Fortran])
use_f03=no 
AC_RUN_IFELSE([
      program procedure_test
        type container
           integer,         pointer         :: var
           procedure(func), pointer, nopass :: f
        end type container
      
        abstract interface
           logical function func(var)
             integer, pointer :: var
           end function func
        end interface
      
        type(container)            :: c
        integer                    :: var
        procedure(func)            :: test
      
        call add(var, test)
        if (c%f(c%var)) stop 0
        stop 1

      contains
        subroutine add(var,f)
          integer,         target    :: var
          procedure(func), optional  :: f
          c%var => var
          c%f   => f
        end subroutine add
      end program procedure_test
      
      logical function test(var)
        integer, pointer :: var
        test = .true.
      end function test],
               [use_f03=yes],
               [use_f03=no])
if test "x${use_f03}" = xyes; then
    AC_SUBST(CPPFLAGS,["$CPPFLAGS -D__F2003"])
#    echo "#define __F2003" >> ./src/ppm_define.h
    AC_MSG_RESULT([yes])
else
    AC_MSG_RESULT([no])
fi

if test "x${use_mpi}" = xyes; then
    AC_MSG_CHECKING([if your MPI is usable])
    AC_LINK_IFELSE([
                     program conftest_prog
                     include 'mpif.h'
                     call MPI_Init(ierr)
                     end program conftest_prog],
                     [AC_SUBST(CPPFLAGS,["$CPPFLAGS -D__MPI"])],
#                     [echo "#define __MPI" >> ./src/ppm_define.h],
                     [AC_MSG_ERROR([Cannot build PPM with MPI support like this])])
    AC_MSG_RESULT([yes])
#    AC_CHECK_LIB([],[MPI_Init],[
#                  echo "#define __MPI" >> ./src/ppm_define.h],
#                  [AC_MSG_ERROR([Cannot build PPM with MPI support like this])])
fi
AC_LANG([C])

if test "x${use_mpi}" = xyes; then
    AC_SUBST(FUNIT_FLAGS,["--with-mpi"])
fi

# Checks for libraries and set preprocessing flags

# this test should be omitted if the compiler is intel, they provide their own
# math implementation
#if test "x${ax_cv_c_compiler_vendor}" != xintel; then
#  AC_CHECK_LIB([m],[sqrt])
#fi

# Make sure that the libraries are compiled with Fortran
AC_LANG([Fortran])

# Check the compiler vendor.. this list should be extended of course
# but let's hope that this will work in most cases
AC_MSG_CHECKING([your Fortran compiler vendor])
vendor=none
for vendorkp in gnu:gfortran intel:Intel sun:Sun portland:Portland; do
    str=`echo $vendorkp | cut -d: -f2`
    if `$FC -V 2>&1 | grep -q $str`; then 
        vendor=`echo $vendorkp | cut -d: -f1`;
        AC_MSG_RESULT([$vendor])
        break;
    fi; 
done
if test "x${vendor}" = xnone; then
    AC_MSG_RESULT([unknown])
fi;

# Check for METIS
have_metis=no
AC_SEARCH_LIBS([METIS_PartGraphRecursive],[metis],[have_metis=yes])
if test "x${have_metis}" = xno; then
  echo "------------------------------------------"
  echo " The METIS library is required to build PPM. "
  echo " Either the library is missing or incorrectly compiled/linked."
  echo " Stopping...   Please download at:"
  echo " http://glaros.dtc.umn.edu/gkhome/metis/metis/download"
  echo "------------------------------------------"
  (exit 1); exit 1;
else
  AC_SUBST(CPPFLAGS,["$CPPFLAGS -D__METIS"])
#  echo "#define __METIS" >> ./src/ppm_define.h
fi

# Check to see if linux was enabled with -enable-linux
AC_ARG_ENABLE([linux],[AS_HELP_STRING([--enable-linux],
  [compile for linux (default is no)])],
  [AC_SUBST(CPPFLAGS,["$CPPFLAGS -D__Linux"])]
#  [echo "#define __Linux" >> ./src/ppm_define.h]
)

# Check to see if etime was enabled with -enable-etime
AC_ARG_ENABLE([etime],[AS_HELP_STRING([--enable-etime],
  [use etime (default is no)])],
  [AC_SUBST(CPPFLAGS,["$CPPFLAGS -D__ETIME"])]
#  [echo "#define __ETIME" >> ./src/ppm_define.h]
)

# Check to see if vector was enabled with -enable-vector
AC_ARG_ENABLE([vector],[AS_HELP_STRING([--enable-vector],
  [enable __VECTOR (default is no)])],
  [AC_SUBST(CPPFLAGS,["$CPPFLAGS -D__VECTOR"])]
#  [echo "#define __VECTOR" >> ./src/ppm_define.h]
)

AC_ARG_ENABLE([mathkeisan],[AS_HELP_STRING([--enable-mathkeisan],
  [enable __MATHKEISAN (default is no)])],
  [AC_SUBST(CPPFLAGS,["$CPPFLAGS -D__MATHKEISAN"])]
#  [echo "#define __MATHKEISAN" >> ./src/ppm_define.h]
)

AC_ARG_ENABLE([sxf90],[AS_HELP_STRING([--enable-sxf90],
  [enable __SXF90 (default is no)])],
  [AC_SUBST(CPPFLAGS,["$CPPFLAGS -D__SXF90"])]
#  [echo "#define __SXF90" >> ./src/ppm_define.h]
)

AC_ARG_ENABLE([crayfishpack],[AS_HELP_STRING([--enable-crayfishpack],
  [enable __CRAYFISHPACK (default is no)])],
  [AC_SUBST(CPPFLAGS,["$CPPFLAGS -D__CRAYFISHPACK"])]
#  [echo "#define __CRAYFISHPACK" >> ./src/ppm_define.h]
)

AC_ARG_ENABLE([hypre],[AS_HELP_STRING([--enable-hypre],
  [enable __HYPRE (default is no)])],
  [AC_SUBST(CPPFLAGS,["$CPPFLAGS -D__HYPRE"])]
#  [echo "#define __HYPRE" >> ./src/ppm_define.h]
)

AC_ARG_ENABLE([no-microinstr],[AS_HELP_STRING([--enable-no-microinstr],
  [enable __NOMICROINSTRUCTIONS (default is no)])],
  [AC_SUBST(CPPFLAGS,["$CPPFLAGS -D__NOMICROINSTRUCTIONS"])]
#  [echo "#define __NOMICROINSTRUCTIONS" >> ./src/ppm_define.h]
)

dev_compile=no
AC_ARG_ENABLE([dev],[AS_HELP_STRING([--enable-dev],
  [enable development mode (default is no)])],[dev_compile=yes])
if test "x${dev_compile}" = xyes; then
    if test "x${vendor}" = xintel; then
        AC_SUBST(DEBUG,["-g -O0 -check pointers -check bounds -traceback -diag-disable 10120"])
    elif test "x${vendor}" = xgnu; then
      AC_SUBST(DEBUG,["-g -O0 -fbounds-check -fbacktrace"])
    else
      AC_SUBST(DEBUG,["-g"])
    fi
fi

debug_compile=no
AC_ARG_ENABLE([debug],[AS_HELP_STRING([--enable-debug],
  [enable debug mode (default is no)])],[debug_compile=yes])
if test "x${debug_compile}" = xyes; then
    AC_SUBST(CPPFLAGS,["$CPPFLAGS -D__DEBUG"])
#    echo "#define __DEBUG" >> ./src/ppm_define.h
    if test "x${vendor}" = xintel; then
        AC_SUBST(DEBUG,["-O0 -g -debug all -check all -fp-stack-check -heap-arrays -ftrapuv -check pointers -check bounds -traceback -diag-disable 10120"])
    elif test "x${vendor}" = xgnu; then
      AC_SUBST(DEBUG,["-O0 -g -fbounds-check -fbacktrace -ftrapv -fno-automatic"])
    else
      AC_SUBST(DEBUG,["-O0 -g"])
    fi
else
    if test "x${vendor}" = xintel; then
        AC_SUBST(FCFLAGS,["$FCFLAGS -diag-disable 5462"])
    fi

fi

# set src as the default srcdir
if test "$srcdir" = "."; then
    srcdir=src
fi

  
# Generate Makefiles
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
