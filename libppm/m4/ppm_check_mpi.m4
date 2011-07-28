# SYNOPSIS
#
#   PPM_CHECK_MPI(mpi implementation)
#
#   the value can be one of
#   openmpi
#   lammpi
#   mpich
#   mpich2
#   intelmpi_gnu
#   intelmpi_intel
#   sun
#   ibm
# DESCRIPTION
#
#
# LICENSE
#

AC_DEFUN([PPM_CHECK_MPI],
[
if test "x$1" = xyes; then
    AC_MSG_NOTICE([*** Check for any MPI to compile and link with ***])
else
    AC_MSG_NOTICE([*** First check if I can use $1 to compile and link with ***])
fi
# check for one of the following MPI distributions
# OpenMPI
# LAMMPI
# MPICH
# MPICH2
# intel MPI
# Sun ClusterTools (tmcc,tmCC,tmf77,tmf90)
# IBM AIX POE (mpcc_r,mpCC_r,mpxlf_r)
# Many of thsoe implementation have same wrapper names
if test "x$1" = xopenmpi; then
    AC_CHECK_PROG(MPICC,[mpicc],[mpicc],[no])
    AS_VAR_IF(MPICC,[no],[AC_MSG_ERROR([Could not find mpicc for OpenMPI])],[CC="$MPICC"])
    AC_CHECK_PROG(MPICXX,[mpic++],[mpic++],[no])
    AS_VAR_IF(MPICXX,[no],[AC_MSG_ERROR([Could not find mpic++ for OpenMPI])],[CXX="$MPICXX"])
    AC_CHECK_PROG(MPIFC,[mpif90],[mpif90],[no])
    AS_VAR_IF(MPIFC,[no],[AC_MSG_ERROR([Could not find mpif90 for OpenMPI])],[FC="$MPIFC"])
elif test "x$1" = xlammpi; then
    AC_CHECK_PROG(MPICC,[mpicc],[mpicc],[no])
    AS_VAR_IF(MPICC,[no],[AC_MSG_ERROR([Could not find mpicc for LAM/MPI])],[CC="$MPICC"])
    AC_CHECK_PROG(MPICXX,[mpiCC],[mpiCC],[no])
    AS_VAR_IF(MPICXX,[no],[AC_MSG_ERROR([Could not find mpiCC for LAM/MPI])],[CXX="$MPICXX"])
    AC_CHECK_PROG(MPIFC,[mpif90],[mpif90],[no])
    AS_VAR_IF(MPIFC,[no],[AC_MSG_ERROR([Could not find mpif90 for LAM/MPI])],[FC="$MPIFC"])
elif test "x$1" = xmpich; then
    AC_CHECK_PROG(MPICC,[mpicc],[mpicc],[no])
    AS_VAR_IF(MPICC,[no],[AC_MSG_ERROR([Could not find mpicc for MPICH])],[CC="$MPICC"])
    AC_CHECK_PROG(MPICXX,[mpiCC],[mpiCC],[no])
    AS_VAR_IF(MPICXX,[no],[AC_MSG_ERROR([Could not find mpiCC for MPICH])],[CXX="$MPICXX"])
    AC_CHECK_PROG(MPIFC,[mpif90],[mpif90],[no])
    AS_VAR_IF(MPIFC,[no],[AC_MSG_ERROR([Could not find mpif90 for MPICH])],[FC="$MPIFC"])
elif test "x$1" = xmpich2; then
    AC_CHECK_PROG(MPICC,[mpicc],[mpicc],[no])
    AS_VAR_IF(MPICC,[no],[AC_MSG_ERROR([Could not find mpicc for MPICH 2])],[CC="$MPICC"])
    AC_CHECK_PROG(MPICXX,[mpicxx],[mpicxx],[no])
    AS_VAR_IF(MPICXX,[no],[AC_MSG_ERROR([Could not find mpicxx for MPICH 2])],[CXX="$MPICXX"])
    AC_CHECK_PROG(MPIFC,[mpif90],[mpif90],[no])
    AS_VAR_IF(MPIFC,[no],[AC_MSG_ERROR([Could not find mpif90 for MPICH 2])],[FC="$MPIFC"])
elif test "x$1" = xintelmpi_gnu; then
    AC_CHECK_PROG(MPICC,[mpigcc],[mpigcc],[no])
    AS_VAR_IF(MPICC,[no],[AC_MSG_ERROR([Could not find mpigcc for Intel MPI + GNU Compilers])],[CC="$MPICC"])
    AC_CHECK_PROG(MPICXX,[mpig++],[mpig++],[no])
    AS_VAR_IF(MPICXX,[no],[AC_MSG_ERROR([Could not find mpig++ for Intel MPI + GNU Compilers])],[CXX="$MPICXX"])
    AC_CHECK_PROG(MPIFC,[mpigfortran],[mpigfortran],[no])
    AS_VAR_IF(MPIFC,[no],[AC_MSG_ERROR([Could not find mpigfortran for Intel MPI + GNU Compilers])],[FC="$MPIFC"])
elif test "x$1" = xintelmpi_intel; then
    AC_CHECK_PROG(MPICC,[mpiicc],[mpiicc],[no])
    AS_VAR_IF(MPICC,[no],[AC_MSG_ERROR([Could not find mpigcc for Intel MPI + Intel Compilers])],[CC="$MPICC"])
    AC_CHECK_PROG(MPICXX,[mpiicpc],[mpicpc],[no])
    AS_VAR_IF(MPICXX,[no],[AC_MSG_ERROR([Could not find mpicpc for Intel MPI + Intel Compilers])],[CXX="$MPICXX"])
    AC_CHECK_PROG(MPIFC,[mpiifort],[mpiifort],[no])
    AS_VAR_IF(MPIFC,[no],[AC_MSG_ERROR([Could not find mpiifort for Intel MPI + Intel Compilers])],[FC="$MPIFC"])
elif test "x$1" = xsun; then
    AC_CHECK_PROG(MPICC,[tmcc],[tmcc],[no])
    AS_VAR_IF(MPICC,[no],[AC_MSG_ERROR([Could not find tmcc for Sun MPI])],[CC="$MPICC"])
    AC_CHECK_PROG(MPICXX,[tmCC],[tmCC],[no])
    AS_VAR_IF(MPICXX,[no],[AC_MSG_ERROR([Could not find tmCC for Sun MPI])],[CXX="$MPICXX"])
    AC_CHECK_PROG(MPIFC,[tmf90],[tmf90],[no])
    AS_VAR_IF(MPIFC,[no],[AC_MSG_ERROR([Could not find tmf90 for Sun MPI])],[FC="$MPIFC"])
elif test "x$1" = xibm; then
    AC_CHECK_PROG(MPICC,[mpcc_r],[mpcc_r],[no])
    AS_VAR_IF(MPICC,[no],[AC_MSG_ERROR([Could not find mpcc_r for IBM AIX MPI])],[CC="$MPICC"])
    AC_CHECK_PROG(MPICXX,[mpCC_r],[mpCC_r],[no])
    AS_VAR_IF(MPICXX,[no],[AC_MSG_ERROR([Could not find mpCC_r for IBM AIX MPI])],[CXX="$MPICXX"])
    AC_CHECK_PROG(MPIFC,[mpxf90_r],[mpxf90_r],[no])
    AS_VAR_IF(MPIFC,[no],[AC_MSG_ERROR([Could not find mpxf90_ for IBM AIX MPI])],[FC="$MPIFC"])
elif test "x$1" = xguess; then
    AC_CHECK_PROGS(MPICC,[mpicc mpigcc mpiicc tmcc mpcc_r],[no])
    AS_VAR_IF(MPICC,[no],[AC_MSG_WARN([Could not find any mpi wrapper for a C compiler])],[CC="$MPICC"])
    AC_CHECK_PROGS(MPICXX,[mpic++ mpiCC mpicxx mpig++ mpiicpc tmCC mpCC_r],[no])
    AS_VAR_IF(MPICXX,[no],[AC_MSG_WARN([Could not find any mpi wrapper for a C++ compiler])],[CXX="$MPICXX"])
    AC_CHECK_PROGS(MPIFC,[mpif90 mpifc mpigfortran mpiifort tmf90 mpxf90_r],[no])
    AS_VAR_IF(MPIFC,[no],[AC_MSG_WARN([Could not find any mpi wrapper for a Fortran compiler])],[FC="$MPIFC"])
    AC_MSG_NOTICE([found $CC $CXX $FC])
elif test "x$1" = xyes; then
    AC_MSG_NOTICE([You would like to use $CC $CXX $FC to compile with MPI support])
else
    AC_MSG_ERROR([Could not reocgnize this MPI implementation: $1])
fi
])

