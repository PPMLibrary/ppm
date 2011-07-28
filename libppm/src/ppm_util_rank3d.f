      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_util_rank3d
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Sort particles in cells. Create index table to 
      !                 particles and pointer to first particle in 
      !                 each cell.
      !
      !  Input        : xp(:,:)    (F) particle co-ordinates
      !                 Np         (I) number of particles
      !                 xmin(:)    (F) minimum extent of mesh (not including
      !                                any ghost layer)
      !                 xmax(:)    (F) maximum extent of mesh (not including
      !                                any ghost layer)
      !                 Nm(:)      (I) number of mesh points (cells)
      !                 Ngl(:)     (I) number of ghost layers to add. The
      !                                ghost layers will be added to the domain
      !                                passed in xmin(:) and xmax()
      !                                Nm(1) : added layers below xmin(1)
      !                                Nm(2) : added layers below xmin(2)
      !                                Nm(3) : added layers below xmin(3)
      !                                Nm(4) : added layers above xmax(1)
      !                                Nm(5) : added layers above xmax(2)
      !                                Nm(6) : added layers above xmax(3)
      !
      !  Input/output : info       (I) return status. If info = 1
      !                                the extent of the particles
      !                                will be checked.
      !                                If particles outside the mesh are
      !                                found, info contains their number
      !                                upon return.
      !
      !  Output       : lpdx(icnt) (I) index of particles in cells
      !                 lhbx(nbx+1)(I) pointer to first particle (in lpdx)
      !                                in each cell
      !
      !  Remarks      : Two do loops do not vectorize. 
      !
      !                 The routine uses no (0) automatic arrays.
      !
      !                 The particles in cell ibox are: 
      !                     lpdx(lhbx(ibox):lhbx(ibox+1)-1)
      !                 We are not using linked lists! as they do not
      !                 vectorize !
      !       
      !                 We are not using automatic arrays as they do not
      !                 tell you if the resources are exhausted!
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_util_rank3d.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.21  2006/09/04 18:35:00  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.19  2006/03/28 19:32:40  ivos
      !  Updated comment header.
      !
      !  Revision 1.18  2005/02/01 16:36:50  ivos
      !  Changed notice to warning if info.EQ.1 and particles are outside
      !  domain.
      !
      !  Revision 1.17  2004/12/03 17:17:47  ivos
      !  Cosmetics in Log header.
      !
      !  Revision 1.16  2004/10/28 12:38:18  davidch
      !  Fixed numerical bug in cell lists that resulted in real particles 
      !  being treated as ghosts and vice versa. The new ranking and cell 
      !  list routines are supposed to be exact. All epsilons that were added 
      !  to the domains in order to prevent the mentioned problems were
      !  removed since they are no longer needed.
      !  Modified Files:
      !      ppm_util_rank2d.f ppm_util_rank3d.f ppm_util_sort2d.f
      !      ppm_util_sort3d.f ppm_find_duplicates.f ppm_neighlist_clist.f
      !      ppm_error.h
      !
      !  Revision 1.15  2004/10/24 09:00:36  ivos
      !  Cosmetics.
      !
      !  Revision 1.14  2004/10/19 13:36:06  davidch
      !  Fixed some roundoff errors. This now solves the test cases correcty but
      !  some problems with the calling routines remain.
      !
      !  Revision 1.13  2004/10/01 16:09:14  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.12  2004/07/26 07:42:34  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.11  2004/06/15 08:39:06  ivos
      !  bugfix: replaced INT with FLOOR when computing box indices. INT causes
      !  wrong results when argument is negative.
      !
      !  Revision 1.10  2004/06/11 12:04:47  ivos
      !  Added comments about vectorization of all loops.
      !
      !  Revision 1.9  2004/06/10 16:20:05  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.8  2004/03/22 09:53:03  walther
      !  Bug fix: Np can be zero - the building of the tree relies on this.
      !
      !  Revision 1.7  2004/01/23 17:24:19  ivos
      !  Now includes ppm_define.h for the cpp defines and no longer sets them
      !  in the Makefile.
      !
      !  Revision 1.6  2004/01/22 15:02:14  ivos
      !  replaced i with info in all calls to ppm_alloc.
      !
      !  Revision 1.5  2004/01/22 14:48:08  ivos
      !  Added sussess checks to every alloc and dealloc.
      !
      !  Revision 1.4  2004/01/22 13:27:57  ivos
      !  Did (1) update of the header, (2) replaced pwrite with ppm_write or
      !  ppm_error calls, (3) inserted validity check of INTENT(IN) arguments
      !  where needed.
      !
      !  Revision 1.3  2004/01/09 09:40:22  ivos
      !  Resolved bug ID 0000010. Count of particles outside the mesh is now
      !  correct. See PPM mantis bug tracker for bug details.
      !
      !  Revision 1.2  2004/01/08 12:56:05  ivos
      !  (1) Fixed the info handling. (2) particles outside of the mesh extent are
      !  now simply ignored, but the routine does not fail. This is needed for
      !  cell lists since they are done per sub and not per processor. (3) memory
      !  for the index array is now allocated to the number of particles IN THE
      !  MESH (instead of Np).
      !
      !  Revision 1.1  2004/01/06 13:43:40  ivos
      !  Initial implementation of particle ranking routines (2d and 3d) for cell
      !  lists and sorting routines.
      !
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_rank3d_s(xp,Np,xmin,xmax,Nm,Ngl,lpdx,lhbx,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_rank3d_d(xp,Np,xmin,xmax,Nm,Ngl,lpdx,lhbx,info)
#endif
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
      INTEGER                 , INTENT(IN   ) :: Np
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: xmin,xmax
      INTEGER , DIMENSION(:)  , INTENT(IN   ) :: Nm     
      INTEGER , DIMENSION(:)  , INTENT(IN   ) :: Ngl     
      INTEGER , DIMENSION(:)  , POINTER       :: lpdx
      INTEGER , DIMENSION(:)  , POINTER       :: lhbx
      INTEGER                 , INTENT(INOUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! cell mesh spacing
      REAL(MK)                                :: rdx,rdy,rdz
      ! non-dimensional extent of mesh
      REAL(MK)                                :: x0,y0,z0
      ! mean number of particles per cell
      REAL(MK)                                :: mean
      ! timer
      REAL(MK)                                :: t0
      ! local info level
      INTEGER                                 :: info2
      ! counters
      INTEGER                                 :: i,j,k,icount,ipart,icorr
      INTEGER                                 :: n2,nbox,ibox
      ! work arrays: box idx of each particle, write pointer, number of
      ! particles per box
      INTEGER, DIMENSION(:), POINTER          :: pbox 
      INTEGER, DIMENSION(:), POINTER          :: cbox
      INTEGER, DIMENSION(:)   , POINTER       :: npbx
      ! total number of cells in each direction (including ghost layers)
      INTEGER, DIMENSION(3)                   :: Nmtot
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
      CALL substart('ppm_util_rank3d',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         IF (Np .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_util_rank3d',  &
     &           'Np must be >0',__LINE__,info)
            GOTO 9999
         ENDIF
         IF (Ngl(1) .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_util_rank2d',  &
     &           'Ngl(1) must be >= 0',__LINE__,info)
            GOTO 9999
         ENDIF
         IF (Ngl(2) .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_util_rank2d',  &
     &           'Ngl(2) must be >= 0',__LINE__,info)
            GOTO 9999
         ENDIF
         IF (Ngl(3) .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_util_rank2d',  &
     &           'Ngl(3) must be >= 0',__LINE__,info)
            GOTO 9999
         ENDIF
         IF (Ngl(4) .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_util_rank2d',  &
     &           'Ngl(4) must be >= 0',__LINE__,info)
            GOTO 9999
         ENDIF
         IF (Ngl(5) .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_util_rank2d',  &
     &           'Ngl(5) must be >= 0',__LINE__,info)
            GOTO 9999
         ENDIF
         IF (Ngl(6) .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_util_rank2d',  &
     &           'Ngl(6) must be >= 0',__LINE__,info)
            GOTO 9999
         ENDIF
         IF (Nm(1) .LE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_util_rank3d',  &
     &           'Nm(1) must be >0',__LINE__,info)
            GOTO 9999
         ENDIF
         IF (Nm(2) .LE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_util_rank3d',  &
     &           'Nm(2) must be >0',__LINE__,info)
            GOTO 9999
         ENDIF
         IF (Nm(3) .LE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_util_rank3d',  &
     &           'Nm(3) must be >0',__LINE__,info)
            GOTO 9999
         ENDIF
         IF (xmax(1) .LE. xmin(1)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_util_rank3d',  &
     &           'xmax(1) must be > xmin(1)',__LINE__,info)
            GOTO 9999
         ENDIF
         IF (xmax(2) .LE. xmin(2)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_util_rank3d',  &
     &           'xmax(2) must be > xmin(2)',__LINE__,info)
            GOTO 9999
         ENDIF
         IF (xmax(3) .LE. xmin(3)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_util_rank3d',  &
     &           'xmax(3) must be > xmin(3)',__LINE__,info)
            GOTO 9999
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Compute total number of mesh cells (global and in each direction)
      !-------------------------------------------------------------------------
      Nmtot(1) = Nm(1) + Ngl(1) + Ngl(4)
      Nmtot(2) = Nm(2) + Ngl(2) + Ngl(5)
      Nmtot(3) = Nm(3) + Ngl(3) + Ngl(6)
      nbox  = Nmtot(1) * Nmtot(2) * Nmtot(3)
      
      !-------------------------------------------------------------------------
      !  Allocate memory
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldc(1) = Np
      CALL ppm_alloc(pbox,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_util_rank3d',     &
     &        'box index of particles PBOX',__LINE__,info)
         GOTO 9999
      ENDIF
      ldc(1) = nbox
      CALL ppm_alloc(cbox,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_util_rank3d',     &
     &        'work array CBOX',__LINE__,info)
         GOTO 9999
      ENDIF
      CALL ppm_alloc(npbx,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_util_rank3d',     &
     &        'number of particles per box NPBX',__LINE__,info)
         GOTO 9999
      ENDIF
      ldc(1) = nbox + 1
      CALL ppm_alloc(lhbx,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_util_rank3d',     &
     &        'first particle in each cell LHBX',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Compute mesh spacing
      !-------------------------------------------------------------------------
      rdx   = REAL(Nm(1),MK)/(xmax(1) - xmin(1))
      rdy   = REAL(Nm(2),MK)/(xmax(2) - xmin(2))
      rdz   = REAL(Nm(3),MK)/(xmax(3) - xmin(3))

      !-------------------------------------------------------------------------
      !  Compute non-dimensional cell co-ordinates (min. extent of mesh)
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Due to round-off errors the lower boundary cannot be pre computed.
      !  (dach)
      !-------------------------------------------------------------------------
      x0    = xmin(1)
      y0    = xmin(2)
      z0    = xmin(3)

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
         DO ipart=1,Np
            IF (xp(1,ipart).LT.xmin(1).OR.xp(1,ipart).GE.xmax(1).OR. &
     &          xp(2,ipart).LT.xmin(2).OR.xp(2,ipart).GE.xmax(2).OR. &
     &          xp(3,ipart).LT.xmin(3).OR.xp(3,ipart).GE.xmax(3)) THEN
               icount  = icount + 1
            ENDIF
         ENDDO

         !----------------------------------------------------------------------
         !  If non-zero, print a warning. This does not have to be bad since
         !  normally, not all the particles of a CPU reside on a single sub
         !  and cell lists are built per sub.
         !----------------------------------------------------------------------
         IF (icount .GT. 0) THEN
            WRITE(msg,'(I8,A)') icount,' particles'
            info = ppm_error_warning
            CALL ppm_error(ppm_err_part_range,'ppm_util_rank3d',msg,   &
     &          __LINE__,info)
         ENDIF
      ENDIF
     
      !-------------------------------------------------------------------------
      !  Find the location of the particles in the boxes. This vectorizes.
      !-------------------------------------------------------------------------
      icount = 0
      info   = 0
      icorr  = 0
      n2     = Nmtot(1) * Nmtot(2)
      DO ipart=1,Np
         !----------------------------------------------------------------------
         !  This has to be a FLOOR and not an INT. The latter would give
         !  wrong results with negative box indices!
         !----------------------------------------------------------------------
         !----------------------------------------------------------------------
         !  The subtraction has to come before the multiplication due to 
         !  Nummerical errors. (dach)
         !----------------------------------------------------------------------
         i = FLOOR((xp(1,ipart) - x0) * rdx) + Ngl(1)
         j = FLOOR((xp(2,ipart) - y0) * rdy) + Ngl(2)
         k = FLOOR((xp(3,ipart) - z0) * rdz) + Ngl(3)
         
         ! The calculated indices are only correct on the lower boundary.
         ! On the upper boundary it may happen that particles inside the 
         ! physical subdomain get an index of a "ghost" cell
         ! We therefore have to test for that and correct this error...
         
         ! if particle is outside the physical domain but index belongs to a
         ! real cell -> move particle to ghost cell 
         if (xp(1,ipart) .GE. xmax(1) .AND. i .LT. Nm(1)+Ngl(1)) THEN
            i = Nm(1) + Ngl(1)
            icorr = icorr + 1
         ENDIF
         if (xp(2,ipart) .GE. xmax(2) .AND. j .LT. Nm(2)+Ngl(2)) THEN
            j = Nm(2) + Ngl(2)
            icorr = icorr + 1
         ENDIF
         if (xp(3,ipart) .GE. xmax(3) .AND. k .LT. Nm(3)+Ngl(3)) THEN
            k = Nm(3) + Ngl(3)
            icorr = icorr + 1
         ENDIF
         
         ! if particle is inside the physical domain but index belongs to a
         ! ghost cell -> move particle in real cell
         if (xp(1,ipart) .LT. xmax(1) .AND. i .GE. Nm(1)+Ngl(1)) THEN
            i = Nm(1) + Ngl(1) - 1
            icorr = icorr + 1
         ENDIF
         if (xp(2,ipart) .LT. xmax(2) .AND. j .GE. Nm(2)+Ngl(2)) THEN
            j = Nm(2) + Ngl(2) - 1
            icorr = icorr + 1
         ENDIF
         if (xp(3,ipart) .LT. xmax(3) .AND. k .GE. Nm(3)+Ngl(3)) THEN
            k = Nm(3) + Ngl(3) - 1
            icorr = icorr + 1
         ENDIF

         ! ignore particles outside the mesh (numbering is from 0...n-1
         ! since we are using INT !!!
         IF ((i .GE. 0 .AND. i .LT. Nmtot(1)) .AND.  &
     &       (j .GE. 0 .AND. j .LT. Nmtot(2)) .AND.  &
     &       (k .GE. 0 .AND. k .LT. Nmtot(3))) THEN
            icount      = icount + 1
            ibox        = i + 1 + j*Nmtot(1) + k*n2 
            pbox(ipart) = ibox
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
     &        __LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate the index array of proper size which is number of
      !  particles in the mesh region (does not need to be the full
      !  processor 
      !  domain!)
      !-------------------------------------------------------------------------
      i = 0
      iopt = ppm_param_alloc_fit
      ldc(1) = icount
      CALL ppm_alloc(lpdx,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_util_rank3d',     &
     &        'particle index list LPDX',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Count the number of particles per box (moved out of above loop
      !  since this count does not vectorize)
      !-------------------------------------------------------------------------
      DO ipart=1,Np
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
      !  Map the particles in the correct order, This does not vectorize.
      !-------------------------------------------------------------------------
      DO ipart=1,Np
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
         CALL ppm_write(ppm_rank,'ppm_util_rank3d',msg,j)
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
         !  Should add up to icount
         !----------------------------------------------------------------------
         IF (i.NE.icount) THEN
            WRITE(msg,'(2(A,I10))') 'icount=',icount,' Sum(npbx)=',i
            info = ppm_error_error
            CALL ppm_error(ppm_err_part_unass,'ppm_util_rank3d',msg,__LINE__,j)
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
         CALL ppm_error(ppm_err_dealloc,'ppm_util_rank3d',     &
     &        'number of particles per box NPBX',__LINE__,info)
      ENDIF
      CALL ppm_alloc(cbox,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc,'ppm_util_rank3d',     &
     &        'work array CBOX',__LINE__,info)
      ENDIF
      CALL ppm_alloc(pbox,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc,'ppm_util_rank3d',     &
     &        'box index of particle PBOX',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_rank3d',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_rank3d_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_rank3d_d
#endif
