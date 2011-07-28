      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_util_time
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_time_s(timing)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_time_d(timing)
#endif
      !!! Returns the current cpu time. Uses `ppm_util_time`,
      !!! which uses either `MPI_Wtime`, f90 `CPU_TIME` or `etime`,
      !!! based on how PPM was configured (see
      !!! ./configure --help):
      !!!
      !!! [NOTE]
      !!! etime is an C intrinsic - thus NOT standard fortran.
      !!! We therefore also write it in small letters.
      !!!
      !!! [WARNING]
      !!! This routine should not be called by the PPM client developer. Use
      !!! `ppm_time` instead.
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      IMPLICIT NONE
#include "ppm_define.h"
      INCLUDE 'ppm_param.h'
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Type kind
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      REAL(MK), INTENT(  OUT) :: timing
      !!! Current CPU/wall clock time
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
#ifdef __MPI
      INTEGER   :: info
#endif
#ifdef __ETIME
      REAL(MK) :: array(2)
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------
      REAL(MK), EXTERNAL :: etime
#endif

#ifdef __MPI
      !-------------------------------------------------------------------------
      !  Call the MPI function MPI_Wtime
      !-------------------------------------------------------------------------
      timing = MPI_Wtime()
#else
      !-------------------------------------------------------------------------
      !  Call the C routine: etime to get the cpu time
      !-------------------------------------------------------------------------
#ifdef __ETIME
      timing = etime(array)
#else
      CALL CPU_TIME(timing)
#endif
#endif

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_time_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_time_d
#endif
