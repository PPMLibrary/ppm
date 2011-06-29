      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_decomp_boxsplit_[s,d]
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich), 
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
      SUBROUTINE decomp_bsplit_s(xp,ppb,npbx,kbox,nbox, &
     &                                 min_box,max_box,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE decomp_bsplit_d(xp,ppb,npbx,kbox,nbox, &
     &                                 min_box,max_box,info)
#endif
      !!! This routine splits a (parent) box in its 4 or 8 children for
      !!! 2D and 3D problems. The particles contained within
      !!! the parent box are sorted into the respective child boxes.
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
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: kbox
      !!! ID of the parent box
      REAL(MK), DIMENSION(:,:), INTENT(INOUT) :: xp
      !!! Particle coordinates
      REAL(MK), DIMENSION(:,:), INTENT(INOUT) :: min_box
      !!! Smallest extremum of the sub-domains
      REAL(MK), DIMENSION(:,:), INTENT(INOUT) :: max_box
      !!! Largest extremum of the sub-domains
      INTEGER , DIMENSION(:)  , INTENT(INOUT) :: ppb
      !!! ppb(i_box) returns the first index of the particle in the box of
      !!! index ibox
      INTEGER , DIMENSION(:)  , INTENT(INOUT) :: npbx
      !!! npbx(i_box) returns the number of particles in the box of index ibox
      INTEGER                 , INTENT(INOUT) :: nbox
      !!! Current number of boxes
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                        :: t0
      REAL(MK), DIMENSION(ppm_dim)    :: cen_box
      INTEGER , DIMENSION(:), POINTER :: npbx_temp => NULL()
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
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
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
      IF (ppm_dim.EQ.2) THEN
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
         max_box(1,nbox+1) = cen_box(1)
         max_box(2,nbox+1) = cen_box(2)

         min_box(1,nbox+2) = cen_box(1)
         min_box(2,nbox+2) = min_box(2,kbox)
         max_box(1,nbox+2) = max_box(1,kbox)
         max_box(2,nbox+2) = cen_box(2)

         min_box(1,nbox+3) = min_box(1,kbox)
         min_box(2,nbox+3) = cen_box(2)
         max_box(1,nbox+3) = cen_box(1)
         max_box(2,nbox+3) = max_box(2,kbox)

         min_box(1,nbox+4) = cen_box(1)
         min_box(2,nbox+4) = cen_box(2)
         max_box(1,nbox+4) = max_box(1,kbox)
         max_box(2,nbox+4) = max_box(2,kbox)

         !----------------------------------------------------------------------
         !  Update the box count
         !----------------------------------------------------------------------
         nbox = nbox + 4
      ELSE
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
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_decomp_boxsplit',     &
     &        'deallocation of decomp_boxsplit failed',__LINE__,info)
          GOTO 9999
      ENDIF
      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_decomp_boxsplit',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
        IF (kbox.LT.1.OR.kbox.GT.nbox) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_decomp_boxsplit', &
     &       'kbox must satisfy: 0 < kbox <= nbox',__LINE__,info)
             GOTO 8888
        ENDIF
 8888   CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE decomp_bsplit_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE decomp_bsplit_d
#endif

