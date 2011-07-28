      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_print_defines
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Prints what cpp DEFINES were used to compile the
      !                 ppm library.
      !
      !  Input        : 
      !
      !  Input/output :
      !
      !  Output       : info      (I) return status. 0 upon success.
      !                 
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_print_defines.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.9  2004/11/12 15:30:14  ivos
      !  Added MATHKEISAN and XLF.
      !
      !  Revision 1.8  2004/10/01 16:33:39  ivos
      !  cosmetics.
      !
      !  Revision 1.7  2004/10/01 16:09:12  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.6  2004/07/26 11:49:56  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.5  2004/07/26 07:45:29  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.4  2004/07/14 10:43:07  ivos
      !  Added SXF90.
      !
      !  Revision 1.3  2004/06/10 16:20:04  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.2  2004/06/02 13:39:12  ivos
      !  Changed such that output to stdout only happens for debug levels
      !  .GT.0. Log output is done always.
      !
      !  Revision 1.1  2004/03/01 12:03:16  ivos
      !  Initial implementation. Tested.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_print_defines(info)

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data, ONLY: ppm_kind_double, ppm_rank, ppm_debug
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_write
      USE ppm_module_log
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! timer
      REAL(ppm_kind_double)                  :: t0
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_print_defines',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Print the defines and log them
      !-------------------------------------------------------------------------
#ifdef __MPI
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_print_defines','__MPI defined',info)
      ENDIF
      CALL ppm_log('ppm_print_defines','__MPI defined',info)
#endif
#ifdef __Linux
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_print_defines','__Linux defined',info)
      ENDIF
      CALL ppm_log('ppm_print_defines','__Linux defined',info)
#endif
#ifdef __VECTOR
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_print_defines','__VECTOR defined',info)
      ENDIF
      CALL ppm_log('ppm_print_defines','__VECTOR defined',info)
#endif
#ifdef __SXF90
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_print_defines','__SXF90 defined',info)
      ENDIF
      CALL ppm_log('ppm_print_defines','__SXF90 defined',info)
#endif
#ifdef __ETIME
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_print_defines','__ETIME defined',info)
      ENDIF
      CALL ppm_log('ppm_print_defines','__ETIME defined',info)
#endif
#ifdef __METIS
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_print_defines','__METIS defined',info)
      ENDIF
      CALL ppm_log('ppm_print_defines','__METIS defined',info)
#endif
#ifdef __FFTW
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_print_defines','__FFTW defined',info)
      ENDIF
      CALL ppm_log('ppm_print_defines','__FFTW defined',info)
#endif
#ifdef __HYPRE
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_print_defines','__HYPRE defined',info)
      ENDIF
      CALL ppm_log('ppm_print_defines','__HYPRE defined',info)
#endif
#ifdef __XLF
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_print_defines','__XLF defined',info)
      ENDIF
      CALL ppm_log('ppm_print_defines','__XLF defined',info)
#endif
#ifdef __MATHKEISAN
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_print_defines',    &
     &        '__MATHKEISAN defined',info)
      ENDIF
      CALL ppm_log('ppm_print_defines','__MATHKEISAN defined',info)
#endif
#ifdef __CRAYFISHPACK
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_print_defines',    &
     &        '__CRAYFISHPACK defined',info)
      ENDIF
      CALL ppm_log('ppm_print_defines','__CRAYFISHPACK defined',info)
#endif
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_print_defines',    &
     &        'See ppm_define.h for descriptions.',info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_print_defines',t0,info)
      RETURN
      END SUBROUTINE ppm_print_defines
