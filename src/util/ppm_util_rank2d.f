      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_util_rank2d
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
      SUBROUTINE ppm_util_rank2d_s(xp,np,xmin,xmax,nm,ngl,lpdx,lhbx,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_rank2d_d(xp,np,xmin,xmax,nm,ngl,lpdx,lhbx,info)
#endif
      !!! Sort particles in cells. Create index table to
      !!! particles and pointer to first particle in each cell.
      !!!
      !!! NOTE: Two do loops do not vectorize.
      !!!
      !!! NOTE: The routine uses no (0) automatic arrays.
      !!!
      !!! [NOTE]
      !!! The particles in cell ibox are:                                      +
      !!! `lpdx(lhbx(ibox):lhbx(ibox+1)-1)`                                    +
      !!! We are not using linked lists! as they do not vectorize !
      !!!
      !!! [NOTE]
      !!! We are not using automatic arrays as they do not
      !!! tell you if the resources are exhausted!
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xp
      !!! Particle coordinates
      INTEGER                 , INTENT(IN   ) :: np
      !!! Number of particles
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: xmin
      !!! Minimum extent of mesh (not including any ghost layer)
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: xmax
      !!! Maximum extent of mesh (not including any ghost layer)
      INTEGER , DIMENSION(:)  , INTENT(IN   ) :: nm
      !!! Number of mesh points (cells)
      INTEGER , DIMENSION(:)  , INTENT(IN   ) :: ngl
      !!! Number of ghost layers to add. The ghost layers will be added to the
      !!! domain passed in xmin(:) and xmax()
      !!!
      !!! * nm(1) : added layers below xmin(1)
      !!! * nm(2) : added layers below xmin(2)
      !!! * nm(3) : added layers above xmax(1)
      !!! * nm(4) : added layers above xmax(2)
      INTEGER , DIMENSION(:)  , POINTER       :: lpdx
      !!! Index of particles in cells (lpdx(icnt))
      INTEGER , DIMENSION(:)  , POINTER       :: lhbx
      !!! Pointer to first particle (in lpdx) in each cell (lhbx(nbx+1))
      INTEGER                 , INTENT(INOUT) :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      ! cell mesh spacing
      REAL(MK)                                :: rdx,rdy
      ! non-dimensional extent of mesh
      REAL(MK)                                :: x0,y0
      ! mean number of particles per cell
      REAL(MK)                                :: mean
      ! timer
      REAL(MK)                                :: t0
      ! local info level
      INTEGER                                 :: info2
      ! counters
      INTEGER                                 :: i,j,icount,ipart,icorr
      INTEGER                                 :: nbox,ibox
      ! work arrays: box idx of each particle, write pointer, number of
      ! particles per box
      INTEGER, DIMENSION(:), POINTER          :: pbox  => NULL()
      INTEGER, DIMENSION(:), POINTER          :: cbox  => NULL()
      INTEGER, DIMENSION(:), POINTER          :: npbx  => NULL()
      ! total number of cells in each direction (including ghost layers)
      INTEGER, DIMENSION(2)                   :: nmtot
      CHARACTER(LEN=ppm_char)                 :: msg
      ! dimensions for allocate
      INTEGER, DIMENSION(1)                   :: ldc
      INTEGER                                 :: iopt
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      ! store the input info (substart will reset info to 0)
      info2 = info
      CALL substart('ppm_util_rank2d',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Compute total number of mesh cells (global and in each direction)
      !-------------------------------------------------------------------------
      nmtot(1) = nm(1) + ngl(1) + ngl(3)
      nmtot(2) = nm(2) + ngl(2) + ngl(4)
      nbox  = nmtot(1)*nmtot(2)

      !-------------------------------------------------------------------------
      !  Allocate memory
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldc(1) = np
      CALL ppm_alloc(pbox,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_util_rank2d',     &
     &        'box index of particles PBOX',__LINE__,info)
         GOTO 9999
      ENDIF
      ldc(1) = nbox
      CALL ppm_alloc(cbox,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_util_rank2d',     &
     &        'work array CBOX',__LINE__,info)
         GOTO 9999
      ENDIF
      CALL ppm_alloc(npbx,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_util_rank2d',     &
     &        'number of particles per box NPBX',__LINE__,info)
         GOTO 9999
      ENDIF
      ldc(1) = nbox + 1
      CALL ppm_alloc(lhbx,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_util_rank2d',     &
     &        'first particle in each cell LHBX',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Compute mesh spacing
      !-------------------------------------------------------------------------
      rdx   = REAL(nm(1),MK)/(xmax(1) - xmin(1))
      rdy   = REAL(nm(2),MK)/(xmax(2) - xmin(2))

      !-------------------------------------------------------------------------
      !  Compute non-dimensional cell co-ordinates (min. extent of mesh)
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Due to round-off errors the lower boundary cannot be pre computed.
      !  (dach)
      !-------------------------------------------------------------------------
      x0    = xmin(1)
      y0    = xmin(2)

      !-------------------------------------------------------------------------
      !  Initialise particle counter (Number of Particle in a box)
      !-------------------------------------------------------------------------
      DO ibox=1,nbox
         npbx(ibox) = 0
      ENDDO

      !-------------------------------------------------------------------------
      !  Check extent of particles
      !-------------------------------------------------------------------------
      IF (info2 .EQ. 1) THEN
         !----------------------------------------------------------------------
         !  Count the number of particle outside the mesh. This vectorizes.
         !----------------------------------------------------------------------
         !----------------------------------------------------------------------
         !  The domain is defined not to include the upper boundary. (dach)
         !----------------------------------------------------------------------
         icount  = 0
         DO ipart=1,np
            IF (xp(1,ipart).LT.xmin(1).OR.xp(1,ipart).GE.xmax(1).OR. &
     &          xp(2,ipart).LT.xmin(2).OR.xp(2,ipart).GE.xmax(2)) THEN
               icount  = icount + 1
            ENDIF
         ENDDO

         !----------------------------------------------------------------------
         !  If non-zero, print a warning. This does not have to be bad since
         !  normally, not all the particles of a CPU reside on a single sub
         !  and cell lists are built per sub.
         !----------------------------------------------------------------------
         IF (icount .GT. 0) THEN
            WRITE(msg,'(I8,A)')icount,' particles'
            info = ppm_error_warning
            CALL ppm_error(ppm_err_part_range,'ppm_util_rank2d',msg,  &
     &          __LINE__,info)
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Find the location of the particles in the boxes. This vectorizes.
      !-------------------------------------------------------------------------
      icount = 0
      info   = 0
      icorr = 0
      DO ipart=1,np
         !----------------------------------------------------------------------
         !  This has to be a FLOOR and not an INT. The latter would give
         !  wrong results with negative box indices!
         !----------------------------------------------------------------------
         !----------------------------------------------------------------------
         !  The subtraction has to come before the multiplication due to
         !  Nummerical errors. (dach)
         !----------------------------------------------------------------------
         i = FLOOR((xp(1,ipart) - x0) * rdx) + ngl(1)
         j = FLOOR((xp(2,ipart) - y0) * rdy) + ngl(2)

         ! The calculated indices are only correct on the lower boundary.
         ! On the upper boundary it may happen that particles inside the
         ! physical subdomain get an index of a "ghost" cell
         ! We therefore have to test for that and correct this error...

         ! if particle is outside the physical domain but index belongs to a
         ! real cell -> move particle to ghost cell
         if (xp(1,ipart) .GE. xmax(1) .AND. i .LT. nm(1)+ngl(1)) THEN
            i = nm(1) + ngl(1)
            icorr = icorr + 1
         ENDIF
         if (xp(2,ipart) .GE. xmax(2) .AND. j .LT. nm(2)+ngl(2)) THEN
            j = nm(2) + ngl(2)
            icorr = icorr + 1
         ENDIF

         ! if particle is inside the physical domain but index belongs to a
         ! ghost cell -> move particle in real cell
         if (xp(1,ipart) .LT. xmax(1) .AND. i .GE. nm(1)+ngl(1)) THEN
            i = nm(1) + ngl(1) - 1
            icorr = icorr + 1
         ENDIF
         if (xp(2,ipart) .LT. xmax(2) .AND. j .GE. nm(2)+ngl(2)) THEN
            j = nm(2) + ngl(2) - 1
            icorr = icorr + 1
         ENDIF

         ! ignore particles outside the mesh (numbering is from 0...nmtot(:)-1
         ! since we are using FLOOR !!!
         IF ((i .GE. 0 .AND. i .LT. nmtot(1)) .AND.  &
     &       (j .GE. 0 .AND. j .LT. nmtot(2))) THEN
            icount       = icount + 1
            ibox         = i + 1 + j*nmtot(1)
            pbox(ipart)  = ibox
         ELSE
            ! particle is in no box
            pbox(ipart) = -1
            info        = info + 1
         ENDIF
      ENDDO


      IF (icorr.GT.0) THEN
         WRITE(msg,'(I8,A)')icorr,' particle indices corrected'
         info = ppm_error_notice
         CALL ppm_error(ppm_err_index_corr,'ppm_util_rank2d',msg,  &
     &                  __LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate the index array of proper size which is number of
      !  particles in the mesh region (does not need to be the full processor
      !  domain!)
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldc(1) = icount
      CALL ppm_alloc(lpdx,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_util_rank2d',     &
     &        'particle index list LPDX',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Count the number of particles per box (moved out of above loop
      !  since this count does not vectorize)
      !-------------------------------------------------------------------------
      DO ipart=1,np
         ibox = pbox(ipart)
         IF (ibox .GT. 0) npbx(ibox) = npbx(ibox) + 1
      ENDDO

      !-------------------------------------------------------------------------
      !  Initialize the particle box pointer and the local counter (this
      !  vectorizes)
      !-------------------------------------------------------------------------
      cbox(1) = 1
      lhbx(1) = 1
      DO ibox=2,nbox
         cbox(ibox) = cbox(ibox-1) + npbx(ibox-1)
         lhbx(ibox) = cbox(ibox)
      ENDDO
      lhbx(nbox+1)  = lhbx(nbox) + npbx(nbox)

      !-------------------------------------------------------------------------
      !  Map the particles in the correct order. This does not vectorize.
      !-------------------------------------------------------------------------
      DO ipart=1,np
         ibox = pbox(ipart)
         IF (ibox .GT. 0) THEN
            ! if particle is in any box, add it to index list
            lpdx(cbox(ibox)) = ipart
            cbox(ibox)       = cbox(ibox) + 1
         ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  Compute some statistics
      !-------------------------------------------------------------------------
      IF (info2.EQ.1) THEN
         mean = REAL(icount,MK)/REAL(nbox,MK)
         WRITE(msg,'(A,F8.2)') 'Mean number of particles per cell: ',mean
         CALL ppm_write(ppm_rank,'ppm_util_rank2d',msg,j)
      ENDIF

      !-------------------------------------------------------------------------
      !  Check that npbx adds to icount
      !-------------------------------------------------------------------------
      IF (info2.EQ.1) THEN
         !----------------------------------------------------------------------
         !  Count number of particles in each box. This vectorizes.
         !----------------------------------------------------------------------
         i = 0
         DO ibox=1,nbox
            i = i + npbx(ibox)
         ENDDO

         !----------------------------------------------------------------------
         !  Should add up to icount (number of particles in the mesh area)
         !----------------------------------------------------------------------
         IF (i.NE.icount) THEN
            WRITE(msg,'(2(A,I10))') 'icount=',icount,' Sum(npbx)=',i
            info = ppm_error_error
            CALL ppm_error(ppm_err_part_unass,'ppm_util_rank2d',msg,__LINE__,j)
            GOTO 9999
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Free work memory
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(npbx,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc,'ppm_util_rank2d',     &
     &        'number of particles per box NPBX',__LINE__,info)
      ENDIF
      CALL ppm_alloc(cbox,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc,'ppm_util_rank2d',     &
     &        'work array CBOX',__LINE__,info)
      ENDIF
      CALL ppm_alloc(pbox,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc,'ppm_util_rank2d',     &
     &        'box index of particle PBOX',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_rank2d',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
         IF (np .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_util_rank2d',  &
     &           'np must be >0',__LINE__,info)
            GOTO 8888
         ENDIF
         IF (ngl(1) .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_util_rank2d',  &
     &           'ngl(1) must be >= 0',__LINE__,info)
            GOTO 8888
         ENDIF
         IF (ngl(2) .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_util_rank2d',  &
     &           'ngl(2) must be >= 0',__LINE__,info)
            GOTO 8888
         ENDIF
         IF (ngl(3) .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_util_rank2d',  &
     &           'ngl(3) must be >= 0',__LINE__,info)
            GOTO 8888
         ENDIF
         IF (ngl(4) .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_util_rank2d',  &
     &           'ngl(4) must be >= 0',__LINE__,info)
            GOTO 8888
         ENDIF
         IF (nm(1) .LE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_util_rank2d',  &
     &           'nm(1) must be >0',__LINE__,info)
            GOTO 8888
         ENDIF
         IF (nm(2) .LE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_util_rank2d',  &
     &           'nm(2) must be >0',__LINE__,info)
            GOTO 8888
         ENDIF
         IF (xmax(1) .LE. xmin(1)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_util_rank2d',  &
     &           'xmax(1) must be > xmin(1)',__LINE__,info)
            GOTO 8888
         ENDIF
         IF (xmax(2) .LE. xmin(2)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_util_rank2d',  &
     &           'xmax(2) must be > xmin(2)',__LINE__,info)
            GOTO 8888
         ENDIF
 8888    CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_rank2d_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_rank2d_d
#endif
