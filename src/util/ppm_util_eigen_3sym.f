      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_util_eigen_3sym
      !-------------------------------------------------------------------------
      ! Copyright (c) 2012 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich), 
      !                    Center for Fluid Dynamics (DTU)
      !
      !
      ! This file is part of the Parallel Particle Mesh Library (PPM).
      !
      ! PPM is free software: you can redistribute it and/or modify
      ! it under the terms of the GNU Lesser General Public License 
      ! as published by the Free Software Foundation, either 
      ! version 3 of the License, or (at your option) any later 
      ! version.
      !
      ! PPM is distributed in the hope that it will be useful,
      ! but WITHOUT ANY WARRANTY; without even the implied warranty of
      ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      ! GNU General Public License for more details.
      !
      ! You should have received a copy of the GNU General Public License
      ! and the GNU Lesser General Public License along with PPM. If not,
      ! see <http://www.gnu.org/licenses/>.
      !
      ! Parallel Particle Mesh Library (PPM)
      ! ETH Zurich
      ! CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_eigen_3sym_s(Am,Eval,Evec,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_eigen_3sym_d(Am,Eval,Evec,info)
#endif
      !!! This routine computes the Eigenvalues and
      !!! orthonormal Eigenvectors of a symmetric, REAL 3x3 matrix.
      !!!
      !!! [NOTE]
      !!! This routine uses explicit formulas for the
      !!! coefficients of the characteristic polynomial and
      !!! solves the resulting cubic equation using the
      !!! routine ppm_util_cubeq_real.
      !!!
      !!! [NOTE]
      !!! This routine only handles symmetric matrices, since
      !!! in this case, all Eigenvalues are REAL and all
      !!! Eigenvectors are pairwise orthogonal.
      !!!
      !!! [TIP]
      !!! For `ppm_debug > 1` the routine checks its own
      !!! result. Maybe - after some time - this can be removed.
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_write
      USE ppm_module_util_cubeq_real
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(3,3), INTENT(IN   ) :: Am
      !!! The 3x3 matrix stored in row-major order (i.e. first index is row
      !!! number, 2nd is column number).
      REAL(MK), DIMENSION(3,3), INTENT(  OUT) :: Evec
      !!! Eigenvalues. Sorted with the largest value first.
      REAL(MK), DIMENSION(3  ), INTENT(  OUT) :: Eval
      !!! Eigenvectors. Evec(:,1) is the first vector (to Eval(1)),
      !!! Evec(:,2) the second and Evec(:,3) the third. Eigenvectors
      !!! are orthonormal when returned.
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                                :: t0,Etmp,lmyeps,inv,t,sqeps
      REAL(MK)                                :: root_lmyeps
      REAL(MK), DIMENSION(3,3)                :: Udet,As
      REAL(MK), DIMENSION(3  )                :: row
      REAL(MK), DIMENSION(4  )                :: chp
      INTEGER , DIMENSION(3  )                :: isort
      INTEGER                                 :: i,j,ip1,np1,np2
      LOGICAL                                 :: correct
      CHARACTER(LEN=ppm_char)                 :: mesg
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_util_eigen_3sym',t0,info)
#if   __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif
      root_lmyeps = SQRT(lmyeps)
      Eval = 0.0_MK
      Evec = 0.0_MK

      !-------------------------------------------------------------------------
      !  Check input arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Precompute sub-determinants. Indices are supressed elements of Am.
      !-------------------------------------------------------------------------
      Udet      = 0.0_MK
      Udet(1,1) = (Am(2,2)*Am(3,3)) - (Am(3,2)*Am(2,3))
      Udet(1,2) = (Am(2,1)*Am(3,3)) - (Am(3,1)*Am(2,3))
      Udet(1,3) = (Am(2,1)*Am(3,2)) - (Am(3,1)*Am(2,2))

      Udet(2,2) = (Am(1,1)*Am(3,3)) - (Am(3,1)*Am(1,3))

      Udet(3,3) = (Am(1,1)*Am(2,2)) - (Am(2,1)*Am(1,2))

      !-------------------------------------------------------------------------
      !  Compute coefficients of characteristic polynomial chp(l)
      !-------------------------------------------------------------------------
      chp(1) =            Am(1,1)*Udet(1,1)
      chp(1) =  chp(1) - (Am(1,2)*Udet(1,2))
      chp(1) =  chp(1) + (Am(1,3)*Udet(1,3))
      chp(1) = -chp(1)

      chp(2) =   Udet(1,1)+ Udet(2,2)+ Udet(3,3)
      chp(3) = -(Am(1,1)  + Am(2,2)  + Am(3,3))
      chp(4) = 1.0_MK
              
      !-------------------------------------------------------------------------
      !  Solve the cubic equation chp=0 for the Eigenvalues
      !------------------------------------------------------------------------
      CALL ppm_util_cubeq_real(chp,Eval,sqeps,info)
      IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,'ppm_util_eigen_3sym',     &
     &       'Matrix is not symmetric and has complex Eigenvalues. Exiting.',&
     &       __LINE__,info)
        GOTO 9999
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Sort Eigenvalues 
      !-------------------------------------------------------------------------
      isort(1:3) = (/1,2,3/)
      IF (Eval(2) .GT. Eval(1)) isort(1:3) = (/2,1,3/)
      IF (Eval(3).GT.Eval(1).AND.Eval(3).GT.Eval(2)) THEN
        isort(3) = isort(2)
        isort(2) = isort(1)
        isort(1) = 3
      ELSEIF (Eval(3).GT.Eval(isort(2))) THEN
        isort(3) = isort(2)
        isort(2) = 3
      ENDIF
      row(1) = Eval(isort(1))
      row(2) = Eval(isort(2))
      row(3) = Eval(isort(3))
      Eval   = row

      !-------------------------------------------------------------------------
      !  Compute the Eigenvectors.
      !  Cases: 3-fold EW: Evec = (e1,e2,e3)
      !         2-fold EW: solve system of rank 1
      !         1-fold EW: solve rank 2 systems
      !-------------------------------------------------------------------------
      IF ((Eval(1).EQ.Eval(2)).AND.(Eval(1).EQ.Eval(3))) THEN
        Evec(1,1) = 1.0_MK
        Evec(2,2) = 1.0_MK
        Evec(3,3) = 1.0_MK
      ELSE
        i = 1
        DO WHILE (i .LT. 4)
            !-----------------------------------------------------------------
            !  Compute the system matrix
            !-----------------------------------------------------------------
            As      = Am
            As(1,1) = Am(1,1) - Eval(i)
            As(2,2) = Am(2,2) - Eval(i)
            As(3,3) = Am(3,3) - Eval(i)

            ip1 = i + 1
            !-----------------------------------------------------------------
            !  Check if this is a double Eigenvalue. Since they are
            !  sorted, we only need to compare i to i+1.
            !-----------------------------------------------------------------
            IF (i .LT. 3) THEN
                IF (Eval(i) .EQ. Eval(ip1)) THEN
                    !---------------------------------------------------------
                    !  Double Eigenvalue: All rows of As are linearly
                    !  dependent. Find the one with the pivot element. All
                    !  others can be made to vanish by Gauss elimination.
                    !---------------------------------------------------------
                    j = 1
                    IF (ABS(As(j,1)) .LT. lmyeps) THEN
                        j = 2
                        IF (ABS(As(j,2)) .LT. lmyeps) THEN
                            j = 3
                            IF (ABS(As(j,3)) .LT. lmyeps) THEN
                                info = ppm_error_error
                                CALL ppm_error(ppm_err_argument,               &
     &                              'ppm_util_eigen_3sym',                    &
     &                              'No pivot found. Exiting.',__LINE__,info)
                                GOTO 9999
                            ENDIF
                        ENDIF
                    ENDIF

                    !---------------------------------------------------------
                    !  Row j contains the pivot. Compute indices of
                    !  non-pivotal columns.
                    !---------------------------------------------------------
                    np1 = 2
                    np2 = 3
                    IF (j .EQ. 2) THEN
                        np1 = 3
                        np2 = 1
                    ELSEIF (j .EQ. 3) THEN
                        np1 = 1
                        np2 = 2
                    ENDIF

                    !---------------------------------------------------------
                    !  The non-pivotal elements of the solution vector can be
                    !  choosen arbitrarily, but linearly independent.
                    !---------------------------------------------------------
                    Evec(np1,i)   = 1.0_MK
                    Evec(np2,i)   = 0.0_MK
                    Evec(np1,ip1) = 0.0_MK
                    Evec(np2,ip1) = 1.0_MK

                    !---------------------------------------------------------
                    !  Compute the pivotal element from the system matrix.
                    !---------------------------------------------------------
                    inv           =  1.0_MK/As(j,j)
                    Evec(j,i)     = -As(j,np1)*inv
                    Evec(j,ip1)   = -As(j,np2)*inv

                    !---------------------------------------------------------
                    !  Normalize the first Eigenvector, taking advantage of
                    !  the special choice of non-pivotal elements.
                    !---------------------------------------------------------
                    inv = Evec(j,i)*Evec(j,i) + 1.0_MK
                    inv = 1.0_MK/SQRT(inv)
                    Evec(np1,i) =             inv
                    Evec(j  ,i) = Evec(j  ,i)*inv

                    !---------------------------------------------------------
                    !  Compute the orthogonal projection (Gram-Schmidt),
                    !  again using the special choice.
                    !---------------------------------------------------------
                    t             = Evec(j,ip1) * Evec(j,i)
                    Evec(np1,ip1) =             - t*Evec(np1,i)
                    Evec(np2,ip1) = 1.0_MK
                    Evec(j  ,ip1) = Evec(j,ip1) - t*Evec(j,i)

                    !---------------------------------------------------------
                    !  Normalize the second Eigenvector.
                    !---------------------------------------------------------
                    inv = Evec(np1,ip1)*Evec(np1,ip1) + 1.0_MK +     &
     &                    Evec(j  ,ip1)*Evec(j  ,ip1)
                    inv = 1.0_MK/SQRT(inv)
                    Evec(np1,ip1) = Evec(np1,ip1)*inv
                    Evec(np2,ip1) =               inv
                    Evec(j  ,ip1) = Evec(j  ,ip1)*inv

                    !---------------------------------------------------------
                    !  Increment the counter of Eigenvectors
                    !---------------------------------------------------------
                    i = i + 2
                ELSE
                    !---------------------------------------------------------
                    !  Single Eigenvalue: There exist exactly two linearly
                    !  independent rows.
                    !  Try if rows 1 and 2 are linearly independent by
                    !  computing the solution using the cross-produce rule
                    !---------------------------------------------------------
                    row(1)  = (As(1,2)*As(2,3)) - (As(1,3)*As(2,2))
                    row(2)  = (As(1,3)*As(2,1)) - (As(1,1)*As(2,3))
                    row(3)  = (As(1,1)*As(2,2)) - (As(1,2)*As(2,1))
                    IF ((ABS(row(1)).LT.root_lmyeps) .AND.         &
     &                  (ABS(row(2)).LT.root_lmyeps) .AND.         &
     &                  (ABS(row(3)).LT.root_lmyeps)) THEN
                        !-----------------------------------------------------
                        !  If rows 1 and 2 are linearly dependent, use
                        !  rows 1 and 3.
                        !  These are now guaranteed to be indepentent.
                        !  Otherwise we would have a double Eigenvalue.
                        !-----------------------------------------------------
                        row(1)  = (As(1,2)*As(3,3)) - (As(1,3)*As(3,2))
                        row(2)  = (As(1,3)*As(3,1)) - (As(1,1)*As(3,3))
                        row(3)  = (As(1,1)*As(3,2)) - (As(1,2)*As(3,1))
                        IF ((ABS(row(1)).LT.root_lmyeps) .AND.     &
     &                      (ABS(row(2)).LT.root_lmyeps) .AND.     &
     &                      (ABS(row(3)).LT.root_lmyeps)) THEN
                            !-------------------------------------------------
                            !  Third possibility: NOT because another
                            !  linearly independent set of rows is
                            !  available, but to get better numerics.
                            !-------------------------------------------------
                            row(1) = (As(2,2)*As(3,3)) - (As(2,3)*As(3,2))
                            row(2) = (As(2,3)*As(3,1)) - (As(2,1)*As(3,3))
                            row(3) = (As(2,1)*As(3,2)) - (As(2,2)*As(3,1))
                        ENDIF
                    ENDIF

                    !---------------------------------------------------------
                    !  Normalize the Eigenvector. They are automatically
                    !  orthogonal by the symmetry of the matrix Am.
                    !---------------------------------------------------------
                    Etmp   = (row(1)*row(1))+(row(2)*row(2))+(row(3)*row(3))
                    IF (ABS(Etmp) .LE. lmyeps) THEN
                        row(1) = 0.0_MK
                        row(2) = 0.0_MK
                        row(3) = 0.0_MK

                    ELSE
                        Etmp   = 1.0_MK/SQRT(Etmp)
                        row(1) = row(1)*Etmp
                        row(2) = row(2)*Etmp
                        row(3) = row(3)*Etmp
                    ENDIF

                    !---------------------------------------------------------
                    !  Store the Eigenvector
                    !---------------------------------------------------------
                    Evec(1,i) = row(1)
                    Evec(2,i) = row(2)
                    Evec(3,i) = row(3)
                    i = i + 1
                ENDIF            ! single or double Eval
            ELSE                ! i.LT.3
                !-------------------------------------------------------------
                !  Single Eigenvalue: There exist exactly two linearly
                !  independent rows.
                !  Try if rows 1 and 2 are linearly independent by computing
                !  the solution using the cross-produce rule
                !-------------------------------------------------------------
                row(1)  = (As(1,2)*As(2,3)) - (As(1,3)*As(2,2))
                row(2)  = (As(1,3)*As(2,1)) - (As(1,1)*As(2,3))
                row(3)  = (As(1,1)*As(2,2)) - (As(1,2)*As(2,1))
                IF ((ABS(row(1)).LT.root_lmyeps) .AND.        &
     &              (ABS(row(2)).LT.root_lmyeps) .AND.        &
     &              (ABS(row(3)).LT.root_lmyeps)) THEN
                    !---------------------------------------------------------
                    !  If rows 1 and 2 are linearly dependent, use rows 1 and
                    !  3. These are now guaranteed to be indepentent.
                    !  Otherwise we would have a double Eigenvalue.
                    !---------------------------------------------------------
                    row(1)  = (As(1,2)*As(3,3)) - (As(1,3)*As(3,2))
                    row(2)  = (As(1,3)*As(3,1)) - (As(1,1)*As(3,3))
                    row(3)  = (As(1,1)*As(3,2)) - (As(1,2)*As(3,1))
                    IF ((ABS(row(1)).LT.root_lmyeps) .AND.     &
     &                  (ABS(row(2)).LT.root_lmyeps) .AND.     &
     &                  (ABS(row(3)).LT.root_lmyeps)) THEN
                        !-----------------------------------------------------
                        !  Third possibility: NOT because another
                        !  linearly independent set of rows is
                        !  available, but to get better numerics.
                        !-----------------------------------------------------
                        row(1) = (As(2,2)*As(3,3)) - (As(2,3)*As(3,2))
                        row(2) = (As(2,3)*As(3,1)) - (As(2,1)*As(3,3))
                        row(3) = (As(2,1)*As(3,2)) - (As(2,2)*As(3,1))
                    ENDIF
                ENDIF

                !-------------------------------------------------------------
                !  Normalize the Eigenvector. They are automatically
                !  orthogonal by the symmetry of the matrix Am.
                !-------------------------------------------------------------
                Etmp   = (row(1)*row(1))+(row(2)*row(2))+(row(3)*row(3))
                IF (ABS(Etmp) .LT. lmyeps) THEN
                    row(1) = 0.0_MK
                    row(2) = 0.0_MK
                    row(3) = 0.0_MK
                ELSE
                    Etmp   = 1.0_MK/SQRT(Etmp)
                    row(1) = row(1)*Etmp
                    row(2) = row(2)*Etmp
                    row(3) = row(3)*Etmp
                ENDIF

                !-------------------------------------------------------------
                !  Store the Eigenvector
                !-------------------------------------------------------------
                Evec(1,i) = row(1)
                Evec(2,i) = row(2)
                Evec(3,i) = row(3)
                i = i + 1
            ENDIF               ! single or double Eval
        ENDDO                  ! WHILE(i.LT.4)
      ENDIF                     ! triple Eval
      
      !-------------------------------------------------------------------------
      !  Check result if debug is enabled.
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
        !---------------------------------------------------------------------
        !  Check that Eigenvalues are sorted
        !---------------------------------------------------------------------
        correct = .TRUE.
        IF (Eval(2) .GT. Eval(1)) correct = .FALSE.
        IF (Eval(3) .GT. Eval(1)) correct = .FALSE.
        IF (Eval(3) .GT. Eval(2)) correct = .FALSE.
        IF (.NOT. correct) THEN
            info = ppm_error_warning
            CALL ppm_error(ppm_err_test_fail,'ppm_util_eigen_3sym',     &
     &           'Eigenvalues are not sorted correctly!',__LINE__,info)
        ELSEIF (ppm_debug .GT. 0) THEN
            CALL ppm_write(ppm_rank,'ppm_util_eigen_3sym',     &
     &           'Eigenvalues are properly sorted.',info)
        ENDIF

        !---------------------------------------------------------------------
        !  Check that Eigenvalues and Eigenvectors are correct
        !---------------------------------------------------------------------
        correct = .TRUE.
        DO i=1,3
            row = MATMUL(Am,Evec(:,i)) - Eval(i)*Evec(:,i)
            IF (ABS(row(1)) .GT. sqeps) THEN
                correct = .FALSE.
                WRITE(mesg,'(A,2(A,E12.4))') 'Eigensystem is not correct.', &
     &              ' Row 1 has error: ',ABS(row(1)),' Tolerance: ',sqeps
                info = ppm_error_warning
                CALL ppm_error(ppm_err_test_fail,'ppm_util_eigen_3sym',     &
     &              mesg,__LINE__,info)
            ENDIF
            IF (ABS(row(2)) .GT. sqeps) THEN
                correct = .FALSE.
                WRITE(mesg,'(A,2(A,E12.4))') 'Eigensystem is not correct.', &
     &              ' Row 2 has error: ',ABS(row(2)),' Tolerance: ',sqeps
                info = ppm_error_warning
                CALL ppm_error(ppm_err_test_fail,'ppm_util_eigen_3sym',     &
     &              mesg,__LINE__,info)
            ENDIF
            IF (ABS(row(3)) .GT. sqeps) THEN
                correct = .FALSE.
                WRITE(mesg,'(A,2(A,E12.4))') 'Eigensystem is not correct.', &
     &              ' Row 3 has error: ',ABS(row(3)),' Tolerance: ',sqeps
                info = ppm_error_warning
                CALL ppm_error(ppm_err_test_fail,'ppm_util_eigen_3sym',     &
     &              mesg,__LINE__,info)
            ENDIF
        ENDDO
        IF (correct .AND. ppm_debug .GT. 0) THEN
            WRITE(mesg,'(2A,E12.4)') 'Eigendecomposition is correct to ',   &
     &           'tolerance ',sqeps
            CALL ppm_write(ppm_rank,'ppm_util_eigen_3sym',mesg,info)
        ENDIF

        !---------------------------------------------------------------------
        !  Check that Eigenvectors are normalized
        !---------------------------------------------------------------------
        correct = .TRUE.
        DO i=1,3
            row(1) = Evec(1,i)*Evec(1,i) + Evec(2,i)*Evec(2,i) +    &
     &           Evec(3,i)*Evec(3,i)
            IF (ABS(row(1)-1.0_MK) .GT. lmyeps) correct = .FALSE.
        ENDDO
        IF (.NOT. correct) THEN
            info = ppm_error_warning
            CALL ppm_error(ppm_err_test_fail,'ppm_util_eigen_3sym',     &
     &           'Eigenvectors are not normalized!',__LINE__,info)
        ELSEIF (ppm_debug .GT. 0) THEN
            WRITE(mesg,'(2A,E12.4)') 'Eigenvectors are normalized to ',   &
     &           'tolerance ',lmyeps
            CALL ppm_write(ppm_rank,'ppm_util_eigen_3sym',mesg,info)
        ENDIF

        !---------------------------------------------------------------------
        !  Check that Eigenvectors are orthogonal
        !---------------------------------------------------------------------
        correct = .TRUE.
        ! 1 -- 2
        row(1) = Evec(1,1)*Evec(1,2)+Evec(2,1)*Evec(2,2)+Evec(3,1)*Evec(3,2)
        ! 1 -- 3
        row(2) = Evec(1,1)*Evec(1,3)+Evec(2,1)*Evec(2,3)+Evec(3,1)*Evec(3,3)
        ! 2 -- 3
        row(3) = Evec(1,2)*Evec(1,3)+Evec(2,2)*Evec(2,3)+Evec(3,2)*Evec(3,3)
        IF (ABS(row(1)) .GT. sqeps) THEN
            correct = .FALSE.
            WRITE(mesg,'(A,2(A,E12.4))') 'Eigenvectors not orthogonal.',  &
     &           ' 1--2 has error: ',ABS(row(1)),' Tolerance: ',sqeps
            info = ppm_error_warning
            CALL ppm_error(ppm_err_test_fail,'ppm_util_eigen_3sym',     &
     &           mesg,__LINE__,info)
        ENDIF
        IF (ABS(row(2)) .GT. sqeps) THEN
            correct = .FALSE.
            WRITE(mesg,'(A,2(A,E12.4))') 'Eigenvectors not orthogonal.',  &
     &           ' 1--3 has error: ',ABS(row(2)),' Tolerance: ',sqeps
            info = ppm_error_warning
            CALL ppm_error(ppm_err_test_fail,'ppm_util_eigen_3sym',     &
     &           mesg,__LINE__,info)
        ENDIF
        IF (ABS(row(3)) .GT. sqeps) THEN
            correct = .FALSE.
            WRITE(mesg,'(A,2(A,E12.4))') 'Eigenvectors not orthogonal.',  &
     &           ' 2--3 has error: ',ABS(row(3)),' Tolerance: ',sqeps
            info = ppm_error_warning
            CALL ppm_error(ppm_err_test_fail,'ppm_util_eigen_3sym',     &
     &           mesg,__LINE__,info)
        ENDIF
        IF (correct .AND. ppm_debug .GT. 0) THEN
            WRITE(mesg,'(2A,E12.4)') 'Eigenvectors are orthogonal to ',   &
     &           'tolerance ',sqeps
            CALL ppm_write(ppm_rank,'ppm_util_eigen_3sym',mesg,info)
        ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_eigen_3sym',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF ((ABS(Am(1,2)-Am(2,1)) .GT. lmyeps) .OR.    &
     &        (ABS(Am(1,3)-Am(3,1)) .GT. lmyeps) .OR.    &
     &        (ABS(Am(2,3)-Am(3,2)) .GT. lmyeps)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_util_eigen_3sym',     &
     &           'Matrix Am must be symmetric !',__LINE__,info)
            GOTO 8888
         ENDIF
 8888    CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_eigen_3sym_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_eigen_3sym_d
#endif
