      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_util_eigen_2sym
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine computes the Eigenvalues and
      !                 orthonormal Eigenvectors of a symmetric, REAL 2x2
      !                 matrix.
      !
      !  Input        : Am(2,2)      (F) The 2x2 matrix stored in row-major
      !                                  order (i.e. first index is row
      !                                  number, 2nd is column number).
      !
      !  Input/output :                                            
      !
      !  Output       : Eval(2)      (F) Eigenvalues. Sorted with the
      !                                  largest value first.
      !                 Evec(2,2)    (F) Eigenvectors. Evec(:,1) is the
      !                                  first vector (to Eval(1)),
      !                                  Evec(:,2) the second. Eigenvectors
      !                                  are orthogonal and normalized 
      !                                  (i.e. of length one) when returned.
      !                 info         (I) return status
      !
      !  Remarks      : This routine uses explicit formulas for the
      !                 coefficients of the characteristic polynomial and
      !                 solves the resulting quadratic equation using the
      !                 routine ppm_util_quadeq_real.
      !
      !                 This routine only handles symmetric matrices, since
      !                 in this case, all Eigenvalues are REAL and all
      !                 Eigenvectors are pairwise orthogonal.
      !
      !                 For ppm_debug .GT. 1 the routine checks its own
      !                 result. Maybe - after some time - this can be
      !                 removed.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_util_eigen_2sym.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.7  2006/09/04 18:34:58  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.6  2004/12/03 17:18:12  ivos
      !  Fixed typo in debug message.
      !
      !  Revision 1.5  2004/11/04 14:49:15  ivos
      !  Some more numerics fixes. Tolerance is now returned from the
      !  equation solvers so it can be used to estimate eigenvevtor tol
      !  in the eigen routines. Error messages now include figures.
      !
      !  Revision 1.4  2004/11/03 16:29:52  ivos
      !  fixed some tolerances in an attempt to make the result checks more
      !  robust against numerical errors.
      !
      !  Revision 1.3  2004/09/22 17:23:46  ivos
      !  Changed result check to use sqrt(eps) instead of eps since the
      !  computed test functions can have an error larger than the result
      !  itself.
      !
      !  Revision 1.2  2004/09/17 14:09:03  ivos
      !  Added result checks at the end if ppm_debug.GT.1.
      !
      !  Revision 1.1  2004/09/17 11:58:48  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_eigen_2sym_s(Am,Eval,Evec,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_eigen_2sym_d(Am,Eval,Evec,info)
#endif

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_write
      USE ppm_module_util_quadeq_real
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
      REAL(MK), DIMENSION(2,2), INTENT(IN   ) :: Am
      REAL(MK), DIMENSION(2,2), INTENT(  OUT) :: Evec
      REAL(MK), DIMENSION(2  ), INTENT(  OUT) :: Eval
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                                :: t0,Etmp,a0,lmyeps,sqeps
      REAL(MK), DIMENSION(2  )                :: row
      REAL(MK), DIMENSION(3  )                :: chp
      INTEGER                                 :: i
      LOGICAL                                 :: correct
      CHARACTER(LEN=ppm_char)                 :: mesg
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_util_eigen_2sym',t0,info)
#if   __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif
      Eval = 0.0_MK
      Evec = 0.0_MK

      !-------------------------------------------------------------------------
      !  Check input arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         IF (ABS(Am(1,2)-Am(2,1)) .GT. lmyeps) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_util_eigen_2sym',     &
     &           'Matrix Am must be symmetric !',__LINE__,info)
            GOTO 9999
         ENDIF 
      ENDIF 

      !-------------------------------------------------------------------------
      !  Coefficients of the characteristic polynomial chp(l)
      !-------------------------------------------------------------------------
      chp(1) =  (Am(1,1)*Am(2,2)) - (Am(2,1)*Am(1,2))
      chp(2) = -(Am(1,1) + Am(2,2))
      chp(3) = 1.0_MK

      !-------------------------------------------------------------------------
      !  Solve the quadratic equation chp = 0 for the Eigenvalues
      !-------------------------------------------------------------------------
      CALL ppm_util_quadeq_real(chp,Eval,sqeps,info) 
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_argument,'ppm_util_eigen_2sym',     &
     &        'Matrix is not symmetric and has complex Eigenvalues. Exiting.',&
     &        __LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Compute Eigenvectors
      !-------------------------------------------------------------------------
      IF (ABS(Eval(1) - Eval(2)) .LT. lmyeps) THEN
         !---------------------------------------------------------------------
         !  Solution space is of dimension 2 => Rank of matrix is 0.
         !  ==> Any two orthogonal vectors of R2 are Evec
         !---------------------------------------------------------------------
         Evec(1,1) = 1.0_MK
         Evec(2,2) = 1.0_MK
      ELSE
         !---------------------------------------------------------------------
         !  Sort Eigenvalues in descending order
         !---------------------------------------------------------------------
         IF (Eval(2) .GT. Eval(1)) THEN
            Etmp    = Eval(1)
            Eval(1) = Eval(2)
            Eval(2) = Etmp
         ENDIF

         !---------------------------------------------------------------------
         !  Compute the two Eigenvectors. They will be orthogonal.
         !---------------------------------------------------------------------
         DO i=1,2
            !-----------------------------------------------------------------
            !  If elm(1,1) is a pivot, the second row can be eliminated
            !  to all 0 (rows are linearly dependent, otherwise Eval
            !  would not be an Eigenvalue) by a Gauss step with 
            !  elm(2,1)/elm(1,1). All we need in this case is the first 
            !  row of A-lambda*I.
            !  If elm(1,1) is not a pivot, we would first flip the rows
            !  and then proceed as above. All we need in this case is
            !  row2 of A-lambda*I.
            !-----------------------------------------------------------------
            a0 = Am(1,1) - Eval(i)
            IF (ABS(a0) .LT. lmyeps) THEN
               row(1)    = Am(2,1)
               row(2)    = Am(2,2) - Eval(i)
               ! freedom of choice in the non-pivot row
               Evec(1,i) = 1.0_MK
               ! solving the pivot row
               IF (ABS(row(1)) .LT. lmyeps) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_util_eigen_2sym', &
     &                 'No pivot found. Exiting.',__LINE__,info)
                  GOTO 9999
               ENDIF
               Evec(2,i) = -(row(2)/row(1))
               Etmp      = (Evec(2,i)*Evec(2,i))+1.0_MK
            ELSE
               row(1)    = a0
               row(2)    = Am(1,2)
               ! solving the pivot row
               Evec(1,i) = -(row(2)/row(1))
               ! freedom of choice in the non-pivot row
               Evec(2,i) = 1.0_MK
               Etmp      = (Evec(1,i)*Evec(1,i))+1.0_MK
            ENDIF

            !-----------------------------------------------------------------
            !  Normalize the vector
            !-----------------------------------------------------------------
            Etmp      = 1.0_MK/SQRT(Etmp)
            Evec(1,i) = Evec(1,i)*Etmp
            Evec(2,i) = Evec(2,i)*Etmp
         ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Check result if debug is enabled.
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
         !---------------------------------------------------------------------
         !  Check that Eigenvalues are sorted
         !---------------------------------------------------------------------
         correct = .TRUE.
         IF (Eval(2) .GT. Eval(1)) correct = .FALSE.
         IF (.NOT. correct) THEN
            info = ppm_error_warning
            CALL ppm_error(ppm_err_test_fail,'ppm_util_eigen_2sym',     &
     &           'Eigenvalues are not sorted correctly!',__LINE__,info)
         ELSE
            CALL ppm_write(ppm_rank,'ppm_util_eigen_2sym',     &
     &           'Eigenvalues are properly sorted.',info)
         ENDIF

         !---------------------------------------------------------------------
         !  Check that Eigenvalues and Eigenvectors are correct
         !---------------------------------------------------------------------
         correct = .TRUE.
         DO i=1,2
            row = MATMUL(Am,Evec(:,i)) - Eval(i)*Evec(:,i)
            IF (ABS(row(1)) .GT. sqeps) THEN
               correct = .FALSE.
               WRITE(mesg,'(A,2(A,E12.4))') 'Eigensystem is not correct.', &
     &              ' Row 1 has error: ',ABS(row(1)),' Tolerance: ',sqeps
               info = ppm_error_warning
               CALL ppm_error(ppm_err_test_fail,'ppm_util_eigen_2sym',     &
     &              mesg,__LINE__,info)
            ENDIF
            IF (ABS(row(2)) .GT. sqeps) THEN
               correct = .FALSE.
               WRITE(mesg,'(A,2(A,E12.4))') 'Eigensystem is not correct.', &
     &              ' Row 2 has error: ',ABS(row(2)),' Tolerance: ',sqeps
               info = ppm_error_warning
               CALL ppm_error(ppm_err_test_fail,'ppm_util_eigen_2sym',     &
     &              mesg,__LINE__,info)
            ENDIF
         ENDDO
         IF (correct .AND. ppm_debug .GT. 0) THEN
            WRITE(mesg,'(2A,E12.4)') 'Eigendecomposition is correct to ',   &
     &           'tolerance ',sqeps
            CALL ppm_write(ppm_rank,'ppm_util_eigen_2sym',mesg,info)
         ENDIF
         
         !---------------------------------------------------------------------
         !  Check that Eigenvectors are normalized
         !---------------------------------------------------------------------
         correct = .TRUE.
         DO i=1,2
            row(1) = Evec(1,i)*Evec(1,i) + Evec(2,i)*Evec(2,i)
            IF (ABS(row(1)-1.0_MK) .GT. lmyeps) correct = .FALSE.
         ENDDO
         IF (.NOT. correct) THEN
            info = ppm_error_warning
            CALL ppm_error(ppm_err_test_fail,'ppm_util_eigen_2sym',     &
     &           'Eigenvectors are not normalized!',__LINE__,info)
         ELSEIF (ppm_debug .GT. 0) THEN
            WRITE(mesg,'(2A,E12.4)') 'Eigenvectors are normalized to ',   &
     &           'tolerance ',lmyeps
            CALL ppm_write(ppm_rank,'ppm_util_eigen_2sym',mesg,info)
         ENDIF

         !---------------------------------------------------------------------
         !  Check that Eigenvectors are orthogonal
         !---------------------------------------------------------------------
         correct = .TRUE.
         row(1) = Evec(1,1)*Evec(1,2) + Evec(2,1)*Evec(2,2)
         IF (ABS(row(1)) .GT. sqeps) THEN
            correct = .FALSE.
            WRITE(mesg,'(A,2(A,E12.4))') 'Eigenvectors not orthogonal.',  &
     &           ' Error: ',ABS(row(1)),' Tolerance: ',sqeps
            info = ppm_error_warning
            CALL ppm_error(ppm_err_test_fail,'ppm_util_eigen_2sym',     &
     &           mesg,__LINE__,info)
         ENDIF
         IF (correct .AND. ppm_debug .GT. 0) THEN
            WRITE(mesg,'(2A,E12.4)') 'Eigenvectors are orthogonal to ',   &
     &           'tolerance ',sqeps
            CALL ppm_write(ppm_rank,'ppm_util_eigen_2sym',mesg,info)
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_eigen_2sym',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_eigen_2sym_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_eigen_2sym_d
#endif
