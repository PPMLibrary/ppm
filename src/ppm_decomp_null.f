      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_decomp_null
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine performs a trivial decomposition
      !                 of the computational domain into one sub. This
      !                 one sub has the same size as the computational
      !                 domain.
      !
      !  Input        : min_phys(:)  (F) the minimum coordinate of the 
      !                                  physical/computational domain 
      !                 max_phys(:)  (F) the maximum coordinate of the 
      !                                  physical/computational domain
      !
      !  Input/output :                                            
      !
      !  Output       : min_sub(:,:) (F) the min. extent of the subdomain
      !                 max_sub(:,:) (F) the max. extent of the subdomain
      !                 nsubs        (I) the total number of subdomains
      !                 info         (I) return status
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_decomp_null.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.6  2006/09/04 18:34:41  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.5  2004/07/26 07:42:38  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.4  2004/06/10 16:19:59  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.3  2004/04/15 16:09:43  ivos
      !  Removed subs2proc assignment from decomp_null and moved it into
      !  its case in topo_mkpart. Mainly because decomp routines are not
      !  expected to do assignments (program logic).
      !
      !  Revision 1.2  2004/04/14 15:06:17  ivos
      !  Added argument checks and fixed wrong bounds in the alloc of min_sub
      !  and max_sub (leading and trailing dimension were mixed up).
      !
      !  Revision 1.1  2004/04/13 15:00:15  oingo
      !  Initial release
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_decomp_null_s(min_phys,max_phys,min_sub,max_sub,nsubs, &
     &    info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_decomp_null_d(min_phys,max_phys,min_sub,max_sub,nsubs, &
     &    info)
#endif

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"

      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: min_phys,max_phys
      REAL(MK), DIMENSION(:,:), POINTER       :: min_sub,max_sub
      INTEGER                 , INTENT(  OUT) :: nsubs
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                    :: t0
      INTEGER, DIMENSION(ppm_dim) :: ldc
      INTEGER                     :: i,iopt

      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_decomp_null',t0,info)

      !-------------------------------------------------------------------------
      !  check input arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         DO i=1,ppm_dim
            IF (max_phys(i) .LE. min_phys(i)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_decomp_null',  &
     &              'max_phys must be > min_phys',__LINE__,info)
               GOTO 9999
            ENDIF
         ENDDO
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Only one sub needs to be defined
      !-------------------------------------------------------------------------
      nsubs = 1
      
      !-------------------------------------------------------------------------
      !  This one sub has the same extent as computational box
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldc(1) = ppm_dim
      ldc(2) = nsubs
      CALL ppm_alloc(min_sub,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_decomp_null',  &
     &        'alloc of min_sub failed',__LINE__,info)
         GOTO 9999
      ENDIF
      CALL ppm_alloc(max_sub,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_decomp_null',  &
     &        'alloc of max_sub failed',__LINE__,info)
         GOTO 9999
      ENDIF

      DO i=1,ppm_dim
         min_sub(i,1) = min_phys(i)
         max_sub(i,1) = max_phys(i)
      ENDDO

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_decomp_null',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_decomp_null_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_decomp_null_d
#endif
