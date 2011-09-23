      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_cnl_MkNeighIdx
      !-------------------------------------------------------------------------
      !
      ! Parallel Particle Mesh Library (PPM)
      ! ETH Zurich
      ! CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE cnl_MkNeighIdx(lsymm,ind,jnd,nnd,ndim,info)
      !!! Creates the index offset list of cell interactions.
      !!! The interaction for the cell with itself is always included as
      !!! the first entry.
      !!!
      !!! [NOTE]
      !!! If the loops do not vectorize, maybe we need to
      !!! duplicate them and put IF(ndim...) statements
      !!! around!
      !!!
      !!! [WARNING]
      !!! `ind` and `jnd` are allocated inside this routine! The
      !!! lists are always (3 x nnd) since then no IF
      !!! statements in the inner loop are needed to
      !!! distinguish between 2D and 3D case. Since `nnd` is
      !!! <28 anyway, this should not be too much of a memory
      !!! waste.

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_error
      USE ppm_module_alloc
      IMPLICIT NONE
      
      INTEGER, PARAMETER :: MK = ppm_kind_double
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      LOGICAL                , INTENT(IN   ) :: lsymm
      !!! T for using symmetry, F for full list
      INTEGER, DIMENSION(:,:), POINTER       :: ind
      !!! First interaction partner (box which *interacts*).
      !!!
      !!! 1st index: 1...3 (x,y,[z]) index shift.                              +
      !!! 2nd index: interaction number 1...nnd.
      INTEGER, DIMENSION(:,:), POINTER       :: jnd
      !!! Second interaction partner (box which *is interacted with*).
      !!! 
      !!! 1st index: 1...3 (x,y,[z]) index shift.                              +
      !!! 2nd index: interaction number 1...nnd.
      INTEGER                , INTENT(  OUT) :: nnd
      !!! Number of box-box interactions to be performed.
      INTEGER                , INTENT(  OUT) :: info
      !!! Returns status, 0 on success
      INTEGER                , INTENT(IN   ) :: ndim
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER                              :: i,j,k,l,ibox,jbox,nz,idm
      REAL(MK)                             :: t0
      ! alloc
      INTEGER, DIMENSION(2)                :: lda
      INTEGER                              :: iopt



      !-------------------------------------------------------------------------
      !  Determine number of box-box interactions needed
      !-------------------------------------------------------------------------
      nnd = 0
      ! 2D using symmetry
      IF (ndim .EQ. 2 .AND. lsymm) nnd = 5
      ! 2D NOT using symmetry
      IF (ndim .EQ. 2 .AND. (.NOT. lsymm)) nnd = 9
      ! 3D using symmetry
      IF (ndim .EQ. 3 .AND. lsymm) nnd = 14
      ! 3D NOT using symmetry
      IF (ndim .EQ. 3 .AND. (.NOT. lsymm)) nnd = 27

      !-------------------------------------------------------------------------
      !  Allocate memory for interaction lists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      lda(1) = 3
      lda(2) = nnd
      CALL ppm_alloc(ind,lda,iopt,i)
      IF (i .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_cnl_MkNeighIdx',    &
     &        'Interaction list IND',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(jnd,lda,iopt,i)
      IF (i .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_cnl_MkNeighIdx',    &
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
      IF (ndim .EQ. 3) THEN
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

         IF (ndim .GT. 2) THEN
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
      RETURN
      END SUBROUTINE cnl_MkNeighIdx
