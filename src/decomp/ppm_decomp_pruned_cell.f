      !-------------------------------------------------------------------------
      !  Subroutine   :              ppm_decomp_pruned_cell
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
      SUBROUTINE decomp_pcell_s(xp,Npart,min_phys,max_phys, &
     &   ghostsize,min_sub,max_sub,nsubs,info,pcost)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE decomp_pcell_d(xp,Npart,min_phys,max_phys, &
     &   ghostsize,min_sub,max_sub,nsubs,info,pcost)
#endif
      !!! This routine performs a domain decomposition using a
      !!! pruned (incomplete) cell index list.
      !!!
      !!! [CAUTION]
      !!! The way we compute the min_sub, max_sub may have
      !!! problems with roundoff errors!

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
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xp
      !!! Position of the particles
      INTEGER                 , INTENT(IN   ) :: Npart
      !!! Number of particles
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: min_phys
      !!! Minimum coordinate of the physical/computational domain
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: max_phys
      !!! Maximum coordinate of the physical/computational domain
      REAL(MK), DIMENSION(:)  , OPTIONAL, INTENT(IN) :: pcost
      !!! Argument of length Npart, specifying the
      !!! computational cost of each particle
      REAL(MK)                , INTENT(IN   ) :: ghostsize
      !!! Size (width) of the ghost layer
      REAL(MK), DIMENSION(:,:), POINTER       :: min_sub
      !!! Minimum extent of the subdomain
      REAL(MK), DIMENSION(:,:), POINTER       :: max_sub
      !!! Maximum extent of the subdomain
      INTEGER                 , INTENT(  OUT) :: nsubs
      !!! Total number of subdomains
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(ppm_dim) :: dmx
      INTEGER , DIMENSION(ppm_dim) :: Nm
      INTEGER , DIMENSION(3)       :: ldc
      INTEGER :: i,j,k,iopt,Mm,ibox,idx,jdx,kdx,ipart
      INTEGER :: istat,n1,n2
      REAL(MK):: rdx,rdy,rdz,x0,y0,z0,rmean_npbx
      REAL(MK):: t0
      CHARACTER(ppm_char) :: mesg
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_decomp_pruned_cell',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Compute the size of the boxes and the required memory
      !-------------------------------------------------------------------------
      Mm = 1
      DO k=1,ppm_dim
         Nm(k)  = INT((max_phys(k) - min_phys(k))/ghostsize)
         dmx(k) = (max_phys(k) - min_phys(k))/REAL(Nm(k),MK)
         Mm     = Mm*Nm(k)
      ENDDO

      !-------------------------------------------------------------------------
      !  Allocate memory for counting the number of particle in the boxes (npbx)
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldc = Mm
      CALL ppm_alloc(npbx,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_decomp_pruned_cell',     &
     &                  'allocation of npbx failed',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Initialize the array
      !-------------------------------------------------------------------------
      DO k=1,Mm
         npbx(k) = 0
      ENDDO

#ifdef __MPI
      !-------------------------------------------------------------------------
      !  In parallel we need memory for the global count
      !-------------------------------------------------------------------------
      CALL ppm_alloc(npbxg,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_decomp_pruned_cell',     &
     &                  'allocation of npbxg failed',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Initialize the array
      !-------------------------------------------------------------------------
      DO k=1,Mm
         npbxg(k) = 0
      ENDDO
 
#endif

      !-------------------------------------------------------------------------
      !  Count the (local) number of particles in the boxes 
      !  (this will and cannot vectorize; of course you can sort them e.g.,
      !  using a cell lists, but then the cell list does not vectorize)
      !-------------------------------------------------------------------------
      IF     (ppm_dim.EQ.2) THEN
         rdx = 1.0_MK/dmx(1)
         rdy = 1.0_MK/dmx(2)
         x0  = min_phys(1)*rdx
         y0  = min_phys(2)*rdy
         n1  = Nm(1)
         DO ipart=1,Npart
            idx  = FLOOR(xp(1,ipart)*rdx - x0)
            jdx  = FLOOR(xp(2,ipart)*rdy - y0)
            ibox = idx + 1 + jdx*n1 
            npbx(ibox) = npbx(ibox) + 1
         ENDDO
      ELSEIF (ppm_dim.EQ.3) THEN
         rdx = 1.0_MK/dmx(1)
         rdy = 1.0_MK/dmx(2)
         rdz = 1.0_MK/dmx(3)
         x0  = min_phys(1)*rdx
         y0  = min_phys(2)*rdy
         z0  = min_phys(3)*rdz
         n1  = Nm(1)
         n2  = Nm(1)*Nm(2)
         DO ipart=1,Npart
            idx  = FLOOR(xp(1,ipart)*rdx - x0)
            jdx  = FLOOR(xp(2,ipart)*rdy - y0)
            kdx  = FLOOR(xp(3,ipart)*rdz - z0)
            ibox = idx + 1 + jdx*n1 + kdx*n2
            npbx(ibox) = npbx(ibox) + 1
         ENDDO
      ENDIF 

#ifdef __MPI
      !-------------------------------------------------------------------------
      !  Add up the number of particles in each box (the the total number of
      !  particles)
      !-------------------------------------------------------------------------
      CALL MPI_AllReduce(npbx,npbxg,Mm,MPI_INTEGER,MPI_SUM,ppm_comm,info)
      DO k=1,Mm
         npbx(k) = npbxg(k)
      ENDDO
#endif 

      !-------------------------------------------------------------------------
      !  Allocate memory for the min_sub and max_sub 
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = ppm_dim
      ldc(2) = Mm
      CALL ppm_alloc(min_sub,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_decomp_pruned_cell',     &
     &                  'allocation of min_sub failed',__LINE__,info)
         GOTO 9999
      ENDIF
      CALL ppm_alloc(max_sub,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_decomp_pruned_cell',     &
     &                  'allocation of max_sub failed',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  loop over the boxes and collect non-empty ones
      !-------------------------------------------------------------------------
      nsubs = 0
      IF     (ppm_dim.EQ.2) THEN
         !----------------------------------------------------------------------
         !  In two dimensions
         !----------------------------------------------------------------------
         DO j=1,Nm(2)
            !-------------------------------------------------------------------
            !  loop over the boxes and collect non-empty ones
            !-------------------------------------------------------------------
            DO i=1,Nm(1)
               ibox = i + (j - 1)*Nm(1)
               IF (npbx(ibox).GT.0) THEN
                  nsubs            = nsubs + 1
                  min_sub(1,nsubs) = min_phys(1) + REAL(i-1,MK)*dmx(1)
                  min_sub(2,nsubs) = min_phys(2) + REAL(j-1,MK)*dmx(2)
                  max_sub(1,nsubs) = min_phys(1) + REAL(i  ,MK)*dmx(1)
                  max_sub(2,nsubs) = min_phys(2) + REAL(j  ,MK)*dmx(2)
               ENDIF 
            ENDDO
         ENDDO
      ELSEIF (ppm_dim.EQ.3) THEN
         !----------------------------------------------------------------------
         !  In three dimensions
         !----------------------------------------------------------------------
         DO k=1,Nm(3)
            DO j=1,Nm(2)
               DO i=1,Nm(1)
                  ibox = i + (j - 1)*Nm(1) + (k - 1)*Nm(1)*Nm(2)
                  IF (npbx(ibox).GT.0) THEN
                     nsubs            = nsubs + 1
                     min_sub(1,nsubs) = min_phys(1) + REAL(i-1,MK)*dmx(1)
                     min_sub(2,nsubs) = min_phys(2) + REAL(j-1,MK)*dmx(2)
                     min_sub(3,nsubs) = min_phys(3) + REAL(k-1,MK)*dmx(3)
                     max_sub(1,nsubs) = min_phys(1) + REAL(i  ,MK)*dmx(1)
                     max_sub(2,nsubs) = min_phys(2) + REAL(j  ,MK)*dmx(2)
                     max_sub(3,nsubs) = min_phys(3) + REAL(k  ,MK)*dmx(3)
                  ENDIF 
               ENDDO
            ENDDO
         ENDDO
      ENDIF 

      !-------------------------------------------------------------------------
      !  Let us shrink the memory to fix exactly the nsubs found
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit_preserve
      ldc(1) = ppm_dim
      ldc(2) = nsubs
      CALL ppm_alloc(min_sub,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_decomp_pruned_cell',     &
     &                  'allocation of min_sub failed',__LINE__,info)
         GOTO 9999
      ENDIF
      CALL ppm_alloc(max_sub,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_decomp_pruned_cell',     &
     &                  'allocation of max_sub failed',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Free the work space again
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(npbx ,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_decomp_pruned_cell',     &
     &       'deallocation of npbx failed',__LINE__,info)
         GOTO 9999
      ENDIF
#ifdef __MPI
      CALL ppm_alloc(npbxg,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_decomp_pruned_cell',     &
     &       'deallocation of npbxg failed',__LINE__,info)
         GOTO 9999
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_decomp_pruned_cell',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
        IF (ghostsize.LE.0.0_MK) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_decomp_pruned_cell',  &
     &                     'the fifth argument must be > 0',__LINE__,info)
            GOTO 8888
        ENDIF
        DO k=1,ppm_dim
            IF (max_phys(k).LE. min_phys(k)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_decomp_pruned_cell',  &
     &                       'min_phys must be < max_phys',__LINE__,info)
               GOTO 8888
            ENDIF
        ENDDO
 8888   CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE decomp_pcell_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE decomp_pcell_d
#endif

