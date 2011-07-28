      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_neighlist_MkNeighIdx
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Creates the index offset list of cell interactions.
      !                 The interaction for the cell with itself is always
      !                 included as the first entry.
      !                  
      !  Input        : lsymm    (L) : T for using symmetry, F for full
      !                                list
      !
      !  Output       : ind(:,:) (I) : First interaction partner (box which
      !                                INTERACTS). First index: 1...3
      !                                (x,y,[z]) index shift. Second index:
      !                                interaction number 1...nnd.
      !                 jnd(:,:) (I) : Second interaction partner (box which
      !                                IS INTERACTED WITH). 
      !                                First index: 1...3
      !                                (x,y,[z]) index shift. Second index:
      !                                interaction number 1...nnd.
      !                 nnd      (I) : number of box-box interactions to be
      !                                performed.
      !                 info     (I) : return status (zero on success)
      !
      !  Remarks      : If the loops do not vectorize, maybe we need to
      !                 duplicate them and put IF(ppm_dim...) statements
      !                 around!
      !
      !                 ind and jnd are allocated inside this routine! The
      !                 lists are always (3 x nnd) since then no IF
      !                 statements in the inner loop are needed to
      !                 distinguish between 2d and 3d case. Since nnd is
      !                 <28 anyway, this should not be too much of a memory
      !                 waste.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_neighlist_MkNeighIdx.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.11  2004/10/01 16:09:11  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.10  2004/07/26 07:42:49  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.9  2004/07/16 14:46:29  ivos
      !  Added check for ppm_initialized.
      !
      !  Revision 1.8  2004/06/15 08:45:15  ivos
      !  Corrected a comment.
      !
      !  Revision 1.7  2004/06/10 16:20:02  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.6  2004/02/24 11:36:38  ivos
      !  Changed imode (INTEGER symmetry flag) argument to lsymm (LOGICAL) in
      !  order to have the same interface as the ghost routines.
      !
      !  Revision 1.5  2004/01/23 17:24:17  ivos
      !  Now includes ppm_define.h for the cpp defines and no longer sets them
      !  in the Makefile.
      !
      !  Revision 1.4  2004/01/22 14:55:44  ivos
      !  Added checks after allocs.
      !
      !  Revision 1.3  2004/01/22 13:27:56  ivos
      !  Did (1) update of the header, (2) replaced pwrite with ppm_write or
      !  ppm_error calls, (3) inserted validity check of INTENT(IN) arguments
      !  where needed.
      !
      !  Revision 1.2  2004/01/09 16:30:59  ivos
      !  Now returned interaction PAIRS for new symmetry handling.
      !
      !  Revision 1.1  2004/01/08 17:51:27  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_neighlist_MkNeighIdx(lsymm,ind,jnd,nnd,info)

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
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
      !  Arguments     
      !-------------------------------------------------------------------------
      LOGICAL                , INTENT(IN   ) :: lsymm
      INTEGER, DIMENSION(:,:), POINTER       :: ind,jnd
      INTEGER                , INTENT(  OUT) :: nnd
      INTEGER                , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER                              :: i,j,k,l,ibox,jbox,nz,idm
      REAL(MK)                             :: t0
      ! alloc
      INTEGER, DIMENSION(2)                :: lda
      INTEGER                              :: iopt
    
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_neighlist_MkNeighIdx',t0,info)
    
      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_neighlist_MkNeighIdx',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Determine number of box-box interactions needed
      !-------------------------------------------------------------------------
      nnd = 0
      ! 2D using symmetry
      IF (ppm_dim .EQ. 2 .AND. lsymm) nnd = 5
      ! 2D NOT using symmetry
      IF (ppm_dim .EQ. 2 .AND. (.NOT. lsymm)) nnd = 9
      ! 3D using symmetry
      IF (ppm_dim .EQ. 3 .AND. lsymm) nnd = 14
      ! 3D NOT using symmetry
      IF (ppm_dim .EQ. 3 .AND. (.NOT. lsymm)) nnd = 27

      !-------------------------------------------------------------------------
      !  Allocate memory for interaction lists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      lda(1) = 3
      lda(2) = nnd
      CALL ppm_alloc(ind,lda,iopt,i)
      IF (i .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_neighlist_MkNeighIdx',    &
     &        'Interaction list IND',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(jnd,lda,iopt,i)
      IF (i .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_neighlist_MkNeighIdx',    &
     &        'Interaction list JND',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Initialize ind and jnd
      !-------------------------------------------------------------------------
      DO j=1,nnd
          DO i=1,3
              ind(i,j) = 0
              jnd(i,j) = 0
          ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Set z direction according to dimensionality
      !-------------------------------------------------------------------------
      nz = 0
      IF (ppm_dim .EQ. 3) THEN
          nz = 1
      ENDIF
    
      !---------------------------------------------------------------------
      !  Compute neighbour indices
      !---------------------------------------------------------------------
      IF (lsymm) THEN
         !------------------------------------------------------------------
         !  Using symmetry
         !------------------------------------------------------------------
                          ! interaction   0 -- 0 (as initialized)

         jnd(1,2) = 1     ! interaction   0 -- 1

         ind(2,3) = 1     ! interaction   0 -- 3

         jnd(1,4) = 1     ! interaction   0 -- 4
         jnd(2,4) = 1

         ind(1,5) = 1     ! interaction   1 -- 3
         jnd(2,5) = 1
          
         IF (ppm_dim .GT. 2) THEN
             jnd(3,5)   = 1   ! interaction   0 -- 9
             ind(1,5)   = 0   ! reset to zero (1-3 will be further down)
             jnd(2,5)   = 0

             jnd(1,6)   = 1   ! interaction   0 -- 10
             jnd(3,6)   = 1

             jnd(2,7)   = 1   ! interaction   0 -- 12
             jnd(3,7)   = 1

             jnd(1,8)   = 1   ! interaction   0 -- 13
             jnd(2,8)   = 1
             jnd(3,8)   = 1

             ind(1,9)   = 1   ! interaction   1 -- 3
             jnd(2,9)   = 1

             ind(1,10)  = 1   ! interaction   1 -- 9
             jnd(3,10)  = 1

             ind(1,11)  = 1   ! interaction   1 -- 12
             jnd(2,11)  = 1
             jnd(3,11)  = 1

             ind(2,12)  = 1   ! interaction   3 -- 9
             jnd(3,12)  = 1

             ind(2,13)  = 1   ! interaction   3 -- 10
             jnd(1,13)  = 1
             jnd(3,13)  = 1

             ind(1,14)  = 1   ! interaction   4 -- 9
             ind(2,14)  = 1
             jnd(3,14)  = 1
         ENDIF
      ELSE
         !------------------------------------------------------------------
         !  Full list
         !------------------------------------------------------------------
         l    = 0
         ibox = 0
         DO k=-nz,nz
            DO j=-1,1
               DO i=-1,1
                  !---------------------------------------------------------
                  !  Add to the interaction lists
                  !---------------------------------------------------------
                  l = l + 1
                  ! center box interacts WITH all boxes around it
                  jnd(1,l) = i
                  jnd(2,l) = j
                  jnd(3,l) = k     ! will always be 0 for 2d case
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_neighlist_MkNeighIdx',t0,info)
      RETURN
      END SUBROUTINE ppm_neighlist_MkNeighIdx
