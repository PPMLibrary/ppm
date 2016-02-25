      !-------------------------------------------------------------------------
      !  Subroutine   :              ppm_decomp_pruned_cell
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
      USE ppm_module_mpi
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
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xp
      !!! Position of the particles
      INTEGER,                  INTENT(IN   ) :: Npart
      !!! Number of particles
      REAL(MK), DIMENSION(:),   POINTER       :: min_phys
      !!! Minimum coordinate of the physical/computational domain
      REAL(MK), DIMENSION(:),   POINTER       :: max_phys
      !!! Maximum coordinate of the physical/computational domain
      REAL(MK),                 INTENT(IN   ) :: ghostsize
      !!! Size (width) of the ghost layer
      REAL(MK), DIMENSION(:,:), POINTER       :: min_sub
      !!! Minimum extent of the subdomain
      REAL(MK), DIMENSION(:,:), POINTER       :: max_sub
      !!! Maximum extent of the subdomain
      INTEGER                 , INTENT(  OUT) :: nsubs
      !!! Total number of subdomains
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      REAL(MK), DIMENSION(:),   INTENT(IN   ), OPTIONAL :: pcost
      !!! Argument of length Npart, specifying the
      !!! computational cost of each particle
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double), DIMENSION(ppm_dim) :: dmx
      REAL(MK)                                  :: rdx,rdy,rdz
      REAL(MK)                                  :: x0,y0,z0,rmean_npbx
      REAL(ppm_kind_double)                     :: t0

      INTEGER, DIMENSION(ppm_dim) :: Nm
      INTEGER, DIMENSION(3)       :: ldc
      INTEGER                     :: i,j,k,iopt,Mm,ibox,idx,jdx,kdx,ipart
      INTEGER                     :: istat,n1,n2
#ifdef __MPI
      INTEGER                     :: request
#endif

      CHARACTER(ppm_char) :: caller = 'ppm_decomp_pruned_cell'
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

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
         Nm(k) =INT(REAL(max_phys(k)-min_phys(k),ppm_kind_double)/REAL(ghostsize,ppm_kind_double))
         dmx(k)=REAL(max_phys(k) - min_phys(k),ppm_kind_double)/REAL(Nm(k),ppm_kind_double)

         !check for round-off problems and fix them if necessary
         DO WHILE (REAL(min_phys(k),ppm_kind_double)+REAL(Nm(k)-1,ppm_kind_double)*dmx(k).LT.REAL(max_phys(k),ppm_kind_double))
            dmx(k)=dmx(k)+EPSILON(dmx(k))
         ENDDO

         check_true(<#(REAL(min_phys(k),ppm_kind_double)+REAL(Nm(k)-1,ppm_kind_double)*dmx(k).GE.REAL(max_phys(k),ppm_kind_double))#>,"round-off problem")

         Mm=Mm*Nm(k)
      ENDDO

      !-------------------------------------------------------------------------
      !  Allocate memory for counting the number of particle in the boxes (npbx)
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldc = Mm
      CALL ppm_alloc(npbx,ldc,iopt,info)
      or_fail_alloc('allocation of npbx failed',ppm_error=ppm_error_fatal)

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
      or_fail_alloc('allocation of npbxg failed')

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
      SELECT CASE (ppm_dim)
      CASE (2)
         rdx = REAL(1.0_ppm_kind_double/dmx(1),MK)
         rdy = REAL(1.0_ppm_kind_double/dmx(2),MK)
         x0  = min_phys(1)*rdx
         y0  = min_phys(2)*rdy
         n1  = Nm(1)
         DO ipart=1,Npart
            idx  = FLOOR(xp(1,ipart)*rdx - x0)
            jdx  = FLOOR(xp(2,ipart)*rdy - y0)
            ibox = idx + 1 + jdx*n1
            npbx(ibox) = npbx(ibox) + 1
         ENDDO

      CASE (3)
         rdx = REAL(1.0_ppm_kind_double/dmx(1),MK)
         rdy = REAL(1.0_ppm_kind_double/dmx(2),MK)
         rdz = REAL(1.0_ppm_kind_double/dmx(3),MK)
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

      END SELECT
#ifdef __MPI
#ifdef __MPI3
      !-------------------------------------------------------------------------
      !  Add up the number of particles in each box (the the total number of
      !  particles)
      !-------------------------------------------------------------------------
      CALL MPI_Iallreduce(npbx,npbxg,Mm,MPI_INTEGER,MPI_SUM,ppm_comm,request,info)
      or_fail_MPI("MPI_Iallreduce")
#else
      !-------------------------------------------------------------------------
      !  Add up the number of particles in each box (the the total number of
      !  particles)
      !-------------------------------------------------------------------------
      CALL MPI_AllReduce(npbx,npbxg,Mm,MPI_INTEGER,MPI_SUM,ppm_comm,info)
      or_fail_MPI("MPI_Allreduce")
      DO k=1,Mm
         npbx(k) = npbxg(k)
      ENDDO
#endif
#endif

      !-------------------------------------------------------------------------
      !  Allocate memory for the min_sub and max_sub
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = ppm_dim
      ldc(2) = Mm
      CALL ppm_alloc(min_sub,ldc,iopt,info)
      or_fail_alloc('allocation of min_sub failed')

      CALL ppm_alloc(max_sub,ldc,iopt,info)
      or_fail_alloc('allocation of max_sub failed')

#ifdef __MPI3
      CALL MPI_Wait(request,MPI_STATUS_IGNORE,info)
      or_fail_MPI("MPI_Wait")

      DO k=1,Mm
         npbx(k)=npbxg(k)
      ENDDO
#endif

      !-------------------------------------------------------------------------
      !  loop over the boxes and collect non-empty ones
      !-------------------------------------------------------------------------
      nsubs = 0
      SELECT CASE (ppm_dim)
      CASE (2)
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
                  min_sub(1,nsubs) = min_phys(1) + REAL(REAL(i-1,ppm_kind_double)*dmx(1),MK)
                  min_sub(2,nsubs) = min_phys(2) + REAL(REAL(j-1,ppm_kind_double)*dmx(2),MK)
                  max_sub(1,nsubs) = min_phys(1) + REAL(REAL(i  ,ppm_kind_double)*dmx(1),MK)
                  max_sub(2,nsubs) = min_phys(2) + REAL(REAL(j  ,ppm_kind_double)*dmx(2),MK)
               ENDIF
            ENDDO
         ENDDO

      CASE (3)
         !----------------------------------------------------------------------
         !  In three dimensions
         !----------------------------------------------------------------------
         DO k=1,Nm(3)
            DO j=1,Nm(2)
               DO i=1,Nm(1)
                  ibox = i + (j - 1)*Nm(1) + (k - 1)*Nm(1)*Nm(2)
                  IF (npbx(ibox).GT.0) THEN
                     nsubs            = nsubs + 1
                     min_sub(1,nsubs) = min_phys(1) + REAL(REAL(i-1,ppm_kind_double)*dmx(1),MK)
                     min_sub(2,nsubs) = min_phys(2) + REAL(REAL(j-1,ppm_kind_double)*dmx(2),MK)
                     min_sub(3,nsubs) = min_phys(3) + REAL(REAL(k-1,ppm_kind_double)*dmx(3),MK)
                     max_sub(1,nsubs) = min_phys(1) + REAL(REAL(i  ,ppm_kind_double)*dmx(1),MK)
                     max_sub(2,nsubs) = min_phys(2) + REAL(REAL(j  ,ppm_kind_double)*dmx(2),MK)
                     max_sub(3,nsubs) = min_phys(3) + REAL(REAL(k  ,ppm_kind_double)*dmx(3),MK)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

      END SELECT

      !-------------------------------------------------------------------------
      !  Let us shrink the memory to fix exactly the nsubs found
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit_preserve
      ldc(1) = ppm_dim
      ldc(2) = nsubs
      CALL ppm_alloc(min_sub,ldc,iopt,info)
      or_fail_alloc('allocation of min_sub failed')

      CALL ppm_alloc(max_sub,ldc,iopt,info)
      or_fail_alloc('allocation of max_sub failed')

      !-------------------------------------------------------------------------
      !  Free the work space again
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(npbx ,ldc,iopt,info)
      or_fail_dealloc('deallocation of npbx failed')

#ifdef __MPI
      CALL ppm_alloc(npbxg,ldc,iopt,info)
      or_fail_dealloc('deallocation of npbxg failed')
#endif

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
        IF (ghostsize.LE.0.0_MK) THEN
           fail('the fifth argument must be > 0',exit_point=8888)
        ENDIF
        DO k=1,ppm_dim
            IF (max_phys(k).LE. min_phys(k)) THEN
               fail('min_phys must be < max_phys',exit_point=8888)
            ENDIF
        ENDDO
      8888 CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE decomp_pcell_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE decomp_pcell_d
#endif

