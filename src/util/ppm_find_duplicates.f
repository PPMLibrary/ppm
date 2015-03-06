      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_find_duplicates
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
      SUBROUTINE ppm_find_duplicates_s(adata,lda,Ndata,nident,ident,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_find_duplicates_d(adata,lda,Ndata,nident,ident,info)
      !!! This routine can be used to find duplicate entries
      !!! in a 2D array of leading dimension 2 or 3.
      !!! The comparison is done up to the precision specified
      !!! to `ppm_init`. Fast O(N) search using cell lists is used.
      !!!
      !!! [NOTE]
      !!! In the case of nident = 0, ident(:,:) will not be allocated
      !!! upon return!
#endif
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_write
      USE ppm_module_alloc
      USE ppm_module_util_rank
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
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: adata
      !!! Data to be checked for duplicates. Both single and double precision
      !!! REAL are supported. Leading dimension must be either 2 or 3.
      INTEGER                 , INTENT(IN   ) :: lda
      !!! Leading dimension of data. Must be either 2 or 3.
      INTEGER                 , INTENT(IN   ) :: Ndata
      !!! Number of data items (second dimension of the array).
      INTEGER                 , INTENT(  OUT) :: nident
      !!! Number of identities found
      INTEGER , DIMENSION(:,:), POINTER       :: ident
      !!! Indices of identical data entries. First index: 1,2 [first
      !!! and second index of the identity], second index:
      !!! 1...nident [identity ID]. This will only be allocated IF at
      !!! least one duplicate is found!
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK)                                   :: t0,lmyeps,bsum,split
      REAL(MK), DIMENSION(3)                     :: xmin,xmax,bsize,aspct
      INTEGER , DIMENSION(3)                     :: Nm
      INTEGER , DIMENSION(6)                     :: Ngl
      INTEGER , DIMENSION(2)                     :: ldc
      INTEGER , DIMENSION(:  ), POINTER          :: lpdx => NULL()
      INTEGER , DIMENSION(:  ), POINTER          :: lhbx => NULL()
      INTEGER                                    :: i,nbox,j,jdata,k,kdata
      INTEGER                                    :: iend,isize,iopt
      CHARACTER(LEN=ppm_char)                    :: mesg
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_find_duplicates',t0,info)
      nident = 0
#if   __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  If there are less than 2 points, we are done
      !-------------------------------------------------------------------------
      IF (Ndata .LT. 2) THEN
          IF (ppm_debug .GT. 1) THEN
              CALL ppm_write(ppm_rank,'ppm_find_duplicates',   &
     &            'Less than 2 data points present. Done.',info)
          ENDIF
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Determine min and max extent of cell mesh
      !-------------------------------------------------------------------------
      xmin(1:3) = HUGE(xmin(1))
      xmax(1:3) = -HUGE(xmax(1))

      IF (lda .GT. 2) THEN
          DO i=1,Ndata
              IF (adata(1,i) .GT. xmax(1)) xmax(1) = adata(1,i)
              IF (adata(1,i) .LT. xmin(1)) xmin(1) = adata(1,i)
              IF (adata(2,i) .GT. xmax(2)) xmax(2) = adata(2,i)
              IF (adata(2,i) .LT. xmin(2)) xmin(2) = adata(2,i)
              IF (adata(3,i) .GT. xmax(3)) xmax(3) = adata(3,i)
              IF (adata(3,i) .LT. xmin(3)) xmin(3) = adata(3,i)
          ENDDO
          !---------------------------------------------------------------------
          !  Add 1/100 of MAX margin in each direction
          !---------------------------------------------------------------------
          bsize(1) = MAX(ABS(xmin(1)),ABS(xmax(1)))
          bsize(2) = MAX(ABS(xmin(2)),ABS(xmax(2)))
          bsize(3) = MAX(ABS(xmin(3)),ABS(xmax(3)))
          bsize(1) = 0.01_MK*bsize(1)
          bsize(2) = 0.01_MK*bsize(2)
          bsize(3) = 0.01_MK*bsize(3)
          !---------------------------------------------------------------------
          !  In the case that all points are identical to 0.0
          !  it happens that bsize is 0.0. Set it to 100eps then.
          !---------------------------------------------------------------------
          IF (bsize(1) .LE. lmyeps) bsize(1) = 100.0_MK*lmyeps
          IF (bsize(2) .LE. lmyeps) bsize(2) = 100.0_MK*lmyeps
          IF (bsize(3) .LE. lmyeps) bsize(3) = 100.0_MK*lmyeps
          xmax(1) = xmax(1) + bsize(1)
          xmax(2) = xmax(2) + bsize(2)
          xmax(3) = xmax(3) + bsize(3)
          xmin(1) = xmin(1) - bsize(1)
          xmin(2) = xmin(2) - bsize(2)
          xmin(3) = xmin(3) - bsize(3)
      ELSE
          DO i=1,Ndata
              IF (adata(1,i) .GT. xmax(1)) xmax(1) = adata(1,i)
              IF (adata(1,i) .LT. xmin(1)) xmin(1) = adata(1,i)
              IF (adata(2,i) .GT. xmax(2)) xmax(2) = adata(2,i)
              IF (adata(2,i) .LT. xmin(2)) xmin(2) = adata(2,i)
          ENDDO
          !---------------------------------------------------------------------
          !  Add 1/100 of MAX margin in each direction
          !---------------------------------------------------------------------
          bsize(1) = MAX(ABS(xmin(1)),ABS(xmax(1)))
          bsize(2) = MAX(ABS(xmin(2)),ABS(xmax(2)))
          bsize(1) = 0.01_MK*bsize(1)
          bsize(2) = 0.01_MK*bsize(2)
          !---------------------------------------------------------------------
          !  In the case that all points are identical to 0.0
          !  it happens that bsize is 0.0. Set it to 100eps then.
          !---------------------------------------------------------------------
          IF (bsize(1) .LE. lmyeps) bsize(1) = 100.0_MK*lmyeps
          IF (bsize(2) .LE. lmyeps) bsize(2) = 100.0_MK*lmyeps
          xmax(1) = xmax(1) + bsize(1)
          xmax(2) = xmax(2) + bsize(2)
          xmin(1) = xmin(1) - bsize(1)
          xmin(2) = xmin(2) - bsize(2)
      ENDIF

      !-------------------------------------------------------------------------
      !  Debug output
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
          WRITE(mesg,'(A,2E18.4)') 'Extent of data in 1st dimension: ',  &
     &        xmin(1),xmax(1)
          CALL ppm_write(ppm_rank,'ppm_find_duplicates',mesg,info)
          WRITE(mesg,'(A,2E18.4)') 'Extent of data in 2nd dimension: ',  &
     &        xmin(2),xmax(2)
          CALL ppm_write(ppm_rank,'ppm_find_duplicates',mesg,info)
          IF (lda .GT. 2) THEN
              WRITE(mesg,'(A,2E18.4)') 'Extent of data in 3rd dimension: ',&
     &            xmin(3),xmax(3)
              CALL ppm_write(ppm_rank,'ppm_find_duplicates',mesg,info)
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Determine number of cells in each direction
      !-------------------------------------------------------------------------
      IF (lda .EQ. 2) THEN
          bsize(1) = xmax(1) - xmin(1)
          bsize(2) = xmax(2) - xmin(2)
          bsum     = bsize(1)+bsize(2)
          bsum     = 1.0_MK/bsum
          aspct(1) = bsize(1)*bsum
          aspct(2) = bsize(2)*bsum
          split    = SQRT(REAL(Ndata,MK))
          Nm(1)    = INT(aspct(1)*split)
          Nm(2)    = INT(aspct(2)*split)
          IF (Nm(1) .LT. 1) Nm(1) = 1
          IF (Nm(2) .LT. 1) Nm(2) = 1
      ELSEIF (lda .EQ. 3) THEN
          bsize(1) = xmax(1) - xmin(1)
          bsize(2) = xmax(2) - xmin(2)
          bsize(3) = xmax(3) - xmin(3)
          bsum     = bsize(1)+bsize(2)+bsize(3)
          bsum     = 1.0_MK/bsum
          aspct(1) = bsize(1)*bsum
          aspct(2) = bsize(2)*bsum
          aspct(3) = bsize(3)*bsum
          split    = REAL(Ndata,MK)**(1.0_MK/3.0_MK)
          Nm(1)    = INT(aspct(1)*split)
          Nm(2)    = INT(aspct(2)*split)
          Nm(3)    = INT(aspct(3)*split)
          IF (Nm(1) .LT. 1) Nm(1) = 1
          IF (Nm(2) .LT. 1) Nm(2) = 1
          IF (Nm(3) .LT. 1) Nm(3) = 1
      ENDIF

      !-------------------------------------------------------------------------
      !  Debug output
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
          WRITE(mesg,'(A,I8)') 'Number of cells in 1st dimension: ',Nm(1)
          CALL ppm_write(ppm_rank,'ppm_find_duplicates',mesg,info)
          WRITE(mesg,'(A,I8)') 'Number of cells in 2nd dimension: ',Nm(2)
          CALL ppm_write(ppm_rank,'ppm_find_duplicates',mesg,info)
          IF (lda .GT. 2) THEN
              WRITE(mesg,'(A,I8)') 'Number of cells in 3rd dimension: ',Nm(3)
              CALL ppm_write(ppm_rank,'ppm_find_duplicates',mesg,info)
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Total number of cells
      !-------------------------------------------------------------------------
      IF (lda .GT. 2) THEN
          nbox = Nm(1)*Nm(2)*Nm(3)
      ELSE
          nbox = Nm(1)*Nm(2)
      ENDIF

      !-------------------------------------------------------------------------
      !  Rank data entries
      !-------------------------------------------------------------------------
      Ngl(1:6) = 0
      IF (lda .EQ. 2) THEN
          CALL ppm_util_rank2d(adata(1:lda,1:Ndata),Ndata,xmin(1:2),    &
     &        xmax(1:2),Nm(1:2),Ngl(1:4),lpdx,lhbx,info)
      ELSE
          CALL ppm_util_rank3d(adata(1:lda,1:Ndata),Ndata,xmin(1:3),    &
     &        xmax(1:3),Nm(1:3),Ngl(1:6),lpdx,lhbx,info)
      ENDIF
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_find_duplicates',   &
     &          'Ranking entries failed.',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Initialize counters
      !-------------------------------------------------------------------------
      isize  = 0

      !-------------------------------------------------------------------------
      !  Find identical data points
      !-------------------------------------------------------------------------
      IF (lda .EQ. 3) THEN
          !---------------------------------------------------------------------
          !  Loop over all cells
          !---------------------------------------------------------------------
          DO i=1,nbox
              !-----------------------------------------------------------------
              !  Skip empty cells
              !-----------------------------------------------------------------
              IF ((lhbx(i+1)-lhbx(i)) .LT. 1) CYCLE
              iend = lhbx(i+1)-1
              !-----------------------------------------------------------------
              !  Loop over all entries in the cell
              !-----------------------------------------------------------------
              DO j=lhbx(i),iend
                  jdata = lpdx(j)
                  !-------------------------------------------------------------
                  !  Compare to all other entries of the same cell using
                  !  symmetry
                  !-------------------------------------------------------------
                  DO k=j+1,iend
                      kdata = lpdx(k)
                      IF ((ABS(adata(1,jdata)-adata(1,kdata)).LT.lmyeps).AND.&
     &                    (ABS(adata(2,jdata)-adata(2,kdata)).LT.lmyeps).AND.&
     &                    (ABS(adata(3,jdata)-adata(3,kdata)).LT.lmyeps)) THEN
                          !-----------------------------------------------------
                          !  jdata and kdata are identical
                          !-----------------------------------------------------
                          nident = nident + 1
                          IF (nident .GT. isize) THEN
                              !-------------------------------------------------
                              !  Grow list of identities
                              !-------------------------------------------------
                              iopt   = ppm_param_alloc_grow_preserve
                              ldc(1) = 2
                              ldc(2) = nident
                              CALL ppm_alloc(ident,ldc,iopt,info)
                              IF (info .NE. 0) THEN
                                  info = ppm_error_fatal
                                  CALL ppm_error(ppm_err_alloc,        &
     &                                'ppm_find_duplicates',           &
     &                                'duplicate list IDENT',__LINE__,info)
                                  GOTO 9999
                              ENDIF
                              isize = nident
                          ENDIF
                          ident(1,nident) = jdata
                          ident(2,nident) = kdata
                      ENDIF
                  ENDDO
              ENDDO
          ENDDO
      ELSEIF (lda .EQ. 2) THEN
          !---------------------------------------------------------------------
          !  Loop over all cells
          !---------------------------------------------------------------------
          DO i=1,nbox
              !-----------------------------------------------------------------
              !  Skip empty cells
              !-----------------------------------------------------------------
              IF ((lhbx(i+1)-lhbx(i)) .LT. 1) CYCLE
              iend = lhbx(i+1)-1
              !-----------------------------------------------------------------
              !  Loop over all entries in the cell
              !-----------------------------------------------------------------
              DO j=lhbx(i),iend
                  jdata = lpdx(j)
                  !-------------------------------------------------------------
                  !  Compare to all other entries of the same cell using
                  !  symmetry
                  !-------------------------------------------------------------
                  DO k=j+1,iend
                      kdata = lpdx(k)
                      IF ((ABS(adata(1,jdata)-adata(1,kdata)).LT.lmyeps).AND.&
     &                    (ABS(adata(2,jdata)-adata(2,kdata)).LT.lmyeps)) THEN
                          !-----------------------------------------------------
                          !  jdata and kdata are identical
                          !-----------------------------------------------------
                          nident = nident + 1
                          IF (nident .GT. isize) THEN
                              !-------------------------------------------------
                              !  Grow list of identities
                              !-------------------------------------------------
                              iopt   = ppm_param_alloc_grow_preserve
                              ldc(1) = 2
                              ldc(2) = nident
                              CALL ppm_alloc(ident,ldc,iopt,info)
                              IF (info .NE. 0) THEN
                                  info = ppm_error_fatal
                                  CALL ppm_error(ppm_err_alloc,        &
     &                                'ppm_find_duplicates',           &
     &                                'duplicate list IDENT',__LINE__,info)
                                  GOTO 9999
                              ENDIF
                              isize = nident
                          ENDIF
                          ident(1,nident) = jdata
                          ident(2,nident) = kdata
                      ENDIF
                  ENDDO
              ENDDO
          ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Deallocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(lpdx,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_alloc,'ppm_find_duplicates',    &
     &        'cell particle index list LPDX',__LINE__,info)
      ENDIF
      CALL ppm_alloc(lhbx,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_alloc,'ppm_find_duplicates',    &
     &        'cell head list LHBX',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_find_duplicates',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
        IF (.NOT. ppm_initialized) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_ppm_noinit,'ppm_find_duplicates',  &
     &          'Please call ppm_init first!',__LINE__,info)
           GOTO 8888
        ENDIF
        IF (lda .NE. 2 .AND. lda .NE. 3) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_find_duplicates', &
     &          'lda must be 2 or 3',__LINE__,info)
           GOTO 8888
        ENDIF
        IF (SIZE(adata,1) .LT. lda) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_find_duplicates', &
     &          'leading dimension of adata cannot be smaller than lda', &
     &          __LINE__,info)
           GOTO 8888
        ENDIF
        IF (Ndata .LT. 0) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_find_duplicates', &
     &          'Ndata must be => 0',__LINE__,info)
           GOTO 8888
        ENDIF
        IF (SIZE(adata,2) .LT. Ndata) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_find_duplicates', &
     &          'second dimension of adata cannot be smaller than Ndata', &
     &          __LINE__,info)
           GOTO 8888
        ENDIF
 8888   CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_find_duplicates_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_find_duplicates_d
#endif
