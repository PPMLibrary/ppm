      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_decomp_boxsplit
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine splits a (parent) box in its 4 or 8 
      !                 children for two and three-dimensional problems. 
      !                 The particles contained within the parent box is sorted 
      !                 into the respective child boxes.
      !
      !  Input        : kbox         (I) : the id of the parent box
      !
      !  Input/output : xp(:,:)      (F) : the particle coordinates
      !                 min_box(:,:) (F) : the smallest extremum of the 
      !                                    sub-domains
      !                 max_box(:,:) (F) : largest extremum of the 
      !                                    sub-domains
      !                 ppb(:)       (I) : ppb(ibox) returns the first index of 
      !                                    the particle in the box of index ibox
      !                 npbx(:)      (I) : npbx(ibox) returns the number of 
      !                                    particles in the box of index ibox
      !                 nbox         (I) : current number of boxes
      !                
      !  Output       : info         (I) : return status (zero on success)
      !                
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_decomp_boxsplit.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.10  2006/09/05 08:01:26  pchatela
      !  Proper scaling for REAL comparisons
      !  Added module_alloc to ppm_decomp_boxsplit
      !
      !  Revision 1.9  2006/09/04 18:34:41  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.8  2004/07/26 07:42:36  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.7  2004/02/05 08:59:45  walther
      !  Moved #include to column 1.
      !
      !  Revision 1.6  2004/01/26 12:11:06  walther
      !  Updated the header and are now checking INTENT(IN) variables.
      !
      !  Revision 1.5  2004/01/23 17:24:14  ivos
      !  Now includes ppm_define.h for the cpp defines and no longer sets them
      !  in the Makefile.
      !
      !  Revision 1.4  2004/01/06 13:40:51  ivos
      !  Bug fix: call to ppm_util_sort is now correct.
      !
      !  Revision 1.3  2004/01/06 12:42:39  ivos
      !  Changed npbx_temp from DIMENSION(8) to POINTER since this is now 
      !  allocated within the sort subroutines.
      !
      !  Revision 1.2  2003/12/09 12:27:27  walther
      !  bug fix in 2D call; jdx was not give a value.
      !
      !  Revision 1.1.1.1  2003/11/17 15:13:45  walther
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_decomp_boxsplit_s(xp,ppb,npbx,kbox,nbox, &
     &                                 min_box,max_box,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_decomp_boxsplit_d(xp,ppb,npbx,kbox,nbox, &
     &                                 min_box,max_box,info)
#endif

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_util_sort
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
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: kbox
      REAL(MK), DIMENSION(:,:), INTENT(INOUT) :: xp
      REAL(MK), DIMENSION(:,:), INTENT(INOUT) :: min_box,max_box
      INTEGER , DIMENSION(:)  , INTENT(INOUT) :: ppb
      INTEGER , DIMENSION(:)  , INTENT(INOUT) :: npbx
      INTEGER                 , INTENT(INOUT) :: nbox
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                        :: t0
      REAL(MK), DIMENSION(ppm_dim)    :: cen_box
      INTEGER , DIMENSION(:), POINTER :: npbx_temp
      INTEGER , DIMENSION(ppm_dim)    :: Nm
      INTEGER                         :: idx,jdx,k,iopt
      INTEGER , DIMENSION(1)          :: lda ! dummy for ppm_alloc
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_decomp_boxsplit',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         IF (kbox.LT.1.OR.kbox.GT.nbox) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_decomp_boxsplit', &
     &       'kbox must satisfy: 0 < kbox <= nbox',__LINE__,info)
            GOTO 9999
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Compute the centre of the current kbox
      !-------------------------------------------------------------------------
      DO k=1,ppm_dim
         cen_box(k) = 0.5_MK*(min_box(k,kbox) + max_box(k,kbox))
      ENDDO

      !-------------------------------------------------------------------------
      !  Define the size of the cell index (2 x 2 x 2)
      !-------------------------------------------------------------------------
      Nm  = 2
      idx = ppb(kbox)

      !-------------------------------------------------------------------------
      !  Sort the particle in two dimensions 
      !-------------------------------------------------------------------------
!      IF (ppm_dim.EQ.2) THEN
         !----------------------------------------------------------------------
         !  in two dimensions
         !----------------------------------------------------------------------
         jdx = idx + npbx(kbox) - 1
         CALL ppm_util_sort2d(xp(1:2,idx:jdx),npbx(kbox),            &
     &                        min_box(1:2,kbox),max_box(1:2,kbox), &
     &                        Nm,npbx_temp,info)
         IF (info.NE.0) GOTO 9999

         !----------------------------------------------------------------------
         !  store the pointers to the particles and the number of particles in
         !  the new boxes
         !----------------------------------------------------------------------
         ppb(nbox+1)  = ppb(kbox)
         npbx(nbox+1) = npbx_temp(1)
         DO k=2,4
            ppb(k+nbox)  = ppb(k-1+nbox) + npbx_temp(k-1)
            npbx(nbox+k) = npbx_temp(k)
         ENDDO

         !----------------------------------------------------------------------
         !  Store the corners of the new boxes
         !----------------------------------------------------------------------
         min_box(1,nbox+1) = min_box(1,kbox)
         min_box(2,nbox+1) = min_box(2,kbox)
         min_box(3,nbox+1) = min_box(3,kbox)
         max_box(1,nbox+1) = cen_box(1)
         max_box(2,nbox+1) = cen_box(2)
         max_box(3,nbox+1) = max_box(3,kbox)

         min_box(1,nbox+2) = cen_box(1)
         min_box(2,nbox+2) = min_box(2,kbox)
         min_box(3,nbox+1) = min_box(3,kbox)
         max_box(1,nbox+2) = max_box(1,kbox)
         max_box(2,nbox+2) = cen_box(2)
         max_box(3,nbox+1) = max_box(3,kbox)

         min_box(1,nbox+3) = min_box(1,kbox)
         min_box(2,nbox+3) = cen_box(2)
         min_box(3,nbox+1) = min_box(3,kbox)
         max_box(1,nbox+3) = cen_box(1)
         max_box(2,nbox+3) = max_box(2,kbox)
         max_box(3,nbox+1) = max_box(3,kbox)

         min_box(1,nbox+4) = cen_box(1)
         min_box(2,nbox+4) = cen_box(2)
         min_box(3,nbox+1) = min_box(3,kbox)
         max_box(1,nbox+4) = max_box(1,kbox)
         max_box(2,nbox+4) = max_box(2,kbox)
         max_box(3,nbox+1) = max_box(3,kbox)

         !----------------------------------------------------------------------
         !  Update the box count
         !----------------------------------------------------------------------
         nbox = nbox + 4
!      ELSE
      IF (.FALSE.)
         !----------------------------------------------------------------------
         !  in three dimensions
         !----------------------------------------------------------------------
         jdx = idx + npbx(kbox) - 1
         CALL ppm_util_sort3d(xp(1:3,idx:jdx),npbx(kbox),            &
     &                        min_box(1:3,kbox),max_box(1:3,kbox), &
     &                        Nm,npbx_temp,info)
         IF (info.NE.0) GOTO 9999

         ppb(nbox+1)  = ppb(kbox)
         npbx(nbox+1) = npbx_temp(1)
         DO k=2,8
            ppb(k+nbox)  = ppb(k-1+nbox) + npbx_temp(k-1) 
            npbx(nbox+k) = npbx_temp(k)
         ENDDO

         !----------------------------------------------------------------------
         !  Store the corners of the new boxes
         !----------------------------------------------------------------------
         min_box(1,nbox+1) = min_box(1,kbox)
         min_box(2,nbox+1) = min_box(2,kbox)
         min_box(3,nbox+1) = min_box(3,kbox)
         max_box(1,nbox+1) = cen_box(1)
         max_box(2,nbox+1) = cen_box(2)
         max_box(3,nbox+1) = cen_box(3)

         min_box(1,nbox+2) = cen_box(1)
         min_box(2,nbox+2) = min_box(2,kbox)
         min_box(3,nbox+2) = min_box(3,kbox)
         max_box(1,nbox+2) = max_box(1,kbox)
         max_box(2,nbox+2) = cen_box(2)
         max_box(3,nbox+2) = cen_box(3)

         min_box(1,nbox+3) = min_box(1,kbox)
         min_box(2,nbox+3) = cen_box(2)
         min_box(3,nbox+3) = min_box(3,kbox)
         max_box(1,nbox+3) = cen_box(1)
         max_box(2,nbox+3) = max_box(2,kbox)
         max_box(3,nbox+3) = cen_box(3)

         min_box(1,nbox+4) = cen_box(1)
         min_box(2,nbox+4) = cen_box(2)
         min_box(3,nbox+4) = min_box(3,kbox)
         max_box(1,nbox+4) = max_box(1,kbox)
         max_box(2,nbox+4) = max_box(2,kbox)
         max_box(3,nbox+4) = cen_box(3)

         min_box(1,nbox+5) = min_box(1,kbox)
         min_box(2,nbox+5) = min_box(2,kbox)
         min_box(3,nbox+5) = cen_box(3)
         max_box(1,nbox+5) = cen_box(1)
         max_box(2,nbox+5) = cen_box(2)
         max_box(3,nbox+5) = max_box(3,kbox)

         min_box(1,nbox+6) = cen_box(1)
         min_box(2,nbox+6) = min_box(2,kbox)
         min_box(3,nbox+6) = cen_box(3)
         max_box(1,nbox+6) = max_box(1,kbox)
         max_box(2,nbox+6) = cen_box(2)
         max_box(3,nbox+6) = max_box(3,kbox)

         min_box(1,nbox+7) = min_box(1,kbox)
         min_box(2,nbox+7) = cen_box(2)
         min_box(3,nbox+7) = cen_box(3)
         max_box(1,nbox+7) = cen_box(1)
         max_box(2,nbox+7) = max_box(2,kbox)
         max_box(3,nbox+7) = max_box(3,kbox)

         min_box(1,nbox+8) = cen_box(1)
         min_box(2,nbox+8) = cen_box(2)
         min_box(3,nbox+8) = cen_box(3)
         max_box(1,nbox+8) = max_box(1,kbox)
         max_box(2,nbox+8) = max_box(2,kbox)
         max_box(3,nbox+8) = max_box(3,kbox)

         !----------------------------------------------------------------------
         !  Update the box count
         !----------------------------------------------------------------------
         nbox = nbox + 8
      ENDIF 

      !-------------------------------------------------------------------------
      !  Deallocate local arrays
      !-------------------------------------------------------------------------
      iopt   = ppm_param_dealloc
      lda(1) = 0 ! dummy
      CALL ppm_alloc(npbx_temp,lda,iopt,info)
 
      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_decomp_boxsplit',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_decomp_boxsplit_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_decomp_boxsplit_d
#endif
