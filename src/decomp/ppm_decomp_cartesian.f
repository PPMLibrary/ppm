      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_decomp_cartesian
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
      SUBROUTINE decomp_cart_s(Nm,min_phys,max_phys, &
      &          decomp,min_sub,max_sub,nsubs,info,ndom)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE decomp_cart_d(Nm,min_phys,max_phys, &
      &          decomp,min_sub,max_sub,nsubs,info,ndom)
#endif
      !!! Decomposes space using a Cartesian domain decomposition topology.
      !!!
      !!! [TIP]
      !!! In the user_defined case, the routine must be called with a non-zero
      !!! number of subdomains in `min_sub(:,:)`, `max_sub(:,:)` and this
      !!! user-defined decomposition is used unchanged.                        +
      !!!
      !!! [NOTE]
      !!! This routine is unaware of particles and simply decomposes space into
      !!! equal grid volumes. When used for particles, this should be followed
      !!! by a check of the cost balancing and possibly iterated with increasing
      !!! ndom.
      !!!
      !!! [NOTE]
      !!! The computational domain is subdivided in exactly as many subdomains
      !!! as there are processors. This is done to minimize memory usage of
      !!! ghost layers (grids need ghost layers also at intra-processor
      !!! subdomain boundaries!) and to get good convergence of the grid solver
      !!! (artificial internal boundary conditionds at subdomain boundaries
      !!! deteriorate convergence). A quad-/oct-tree decomposition of the
      !!! computational domain is made until the remaining number of subs
      !!! does not contain anymore a power of 2 factor. The remaining factor is
      !!! used to do a parallel slab decomposition, such that the
      !!! surface-to-volume ratio is minimized (in terms of numbers of grid
      !!! points). This way, nproc equally sized subdomains with minimal
      !!! collective surface nodes are created.                                +
      !!! The pencil decomposition works exactly the same
      !!! except that cuts in one direction are prohibited.
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
      INTEGER,  DIMENSION(:),   INTENT(IN   ) :: Nm
      !!! Number of mesh points (not cells) in each direction of the
      !!! global comput. domain. (including those ON the boundaries)
      REAL(MK), DIMENSION(:),   INTENT(IN   ) :: min_phys
      !!! Minimum coordinate of the physical/computational domain
      REAL(MK), DIMENSION(:),   INTENT(IN   ) :: max_phys
      !!! Maximum coordinate of the physical/computational domain
      INTEGER,                  INTENT(IN   ) :: decomp
      !!! Domain decomposition can be one of:
      !!!
      !!! * ppm_param_decomp_xpencil
      !!! * ppm_param_decomp_ypencil
      !!! * ppm_param_decomp_zpencil
      !!! * ppm_param_decomp_cuboid
      !!! * ppm_param_decomp_user_defined
      REAL(MK), DIMENSION(:,:), POINTER       :: min_sub
      !!! Minimum extend of the subs
      REAL(MK), DIMENSION(:,:), POINTER       :: max_sub
      !!! Maximum extend of the subs
      INTEGER,                  INTENT(INOUT) :: nsubs
      !!! Total number of subdomains
      INTEGER,                  INTENT(  OUT) :: info
      !!! Retuns status, 0 upon success
      INTEGER,  OPTIONAL,       INTENT(IN   ) :: ndom
      !!! Number of subdomains to be created. If not specified,
      !!! the number of subdomains will be equal to the number of processors.
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK)                     :: t0,minsv,lc
      REAL(MK)                     :: lx,ux,ly,uy,lz,uz
      REAL(MK), DIMENSION(ppm_dim) :: dx
      REAL(MK), DIMENSION(ppm_dim) :: len_phys
      ! physical extent of comput. domain
      REAL(MK)                                  :: gs

      INTEGER, DIMENSION(:,:), POINTER :: Npx    => NULL()
      ! number of grid points in each sub. index: (1:ppm_dim,1:nsub)
      INTEGER, DIMENSION(:,:), POINTER :: Npxnew => NULL()
      INTEGER                          :: i,j,k
      INTEGER                          :: tx,ty,tz
      INTEGER                          :: icut,jcut,kcut
      INTEGER                          :: nsrem,ncut
      INTEGER                          :: cutdim,constdim
      INTEGER                          :: ix,iy,iz
      INTEGER                          :: lduold,nup,ndn,iup,idn
      ! power-of-two part of the number of blocks (i.e. part for ROB)
      INTEGER                          :: rc
      ! mesh spacing
      INTEGER, DIMENSION(ppm_dim)      :: nblocks
      ! number of blocks into which to subdivide the space
      INTEGER, DIMENSION(ppm_dim)      :: surface,volume,Nc
      INTEGER, DIMENSION(2)            :: ldu
      INTEGER                          :: iopt

      CHARACTER(LEN=ppm_char) :: mesg
      CHARACTER(LEN=ppm_char) :: caller = 'ppm_decomp_cartesian'
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
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
      !  Number of subdomains
      !-------------------------------------------------------------------------
      IF (PRESENT(ndom)) THEN
         nsubs = ndom
      ELSE
         nsubs = ppm_nproc
      ENDIF

      !-------------------------------------------------------------------------
      !  Mesh spacing
      !-------------------------------------------------------------------------
      Nc(1:ppm_dim)       = Nm(1:ppm_dim)-1
      len_phys(1:ppm_dim) = max_phys(1:ppm_dim) - min_phys(1:ppm_dim)
      dx(1:ppm_dim)       = len_phys(1:ppm_dim)/REAL(Nc(1:ppm_dim),MK)

      !check for round-off problems and fix them if necessary
      DO k=1,ppm_dim
         DO WHILE (min_phys(k)+Nc(k)*dx(k).LT.max_phys(k))
            dx(k)=dx(k)+EPSILON(dx(k))
         ENDDO
      ENDDO
      check_true(<#ALL(min_phys(1:ppm_dim)+Nc(1:ppm_dim)*dx(1:ppm_dim).GE.max_phys(1:ppm_dim))#>,"round-off problem")

      !-------------------------------------------------------------------------
      !  Allocate memory for the subs if not already done so by the user
      !-------------------------------------------------------------------------
      IF (decomp .NE. ppm_param_decomp_user_defined) THEN
         iopt = ppm_param_alloc_fit
         ldu(1) = ppm_dim
         ldu(2) = nsubs
         CALL ppm_alloc(min_sub,ldu,iopt,info)
         or_fail_alloc( 'minima of subdomains MIN_SUB')

         CALL ppm_alloc(max_sub,ldu,iopt,info)
         or_fail_alloc('maxima of subdomains MAX_SUB')
      ENDIF

      !-------------------------------------------------------------------------
      !  pencil and cuboid decompositions
      !-------------------------------------------------------------------------
      SELECT CASE (decomp)
      CASE (ppm_param_decomp_xpencil, &
      &     ppm_param_decomp_ypencil, &
      &     ppm_param_decomp_zpencil, &
      &     ppm_param_decomp_cuboid)
         ! set the constant dimension (i.e. the one that will not be cut)
         constdim = 0      ! 0 means: none -> cuboids
         IF (decomp .EQ. ppm_param_decomp_xpencil) constdim = 1
         IF (decomp .EQ. ppm_param_decomp_ypencil) constdim = 2
         IF (decomp .EQ. ppm_param_decomp_zpencil) constdim = 3
         !---------------------------------------------------------------------
         !  Decomposition based on mesh
         !---------------------------------------------------------------------
         ! create as few subs as possible, but at least as many as
         ! there as processors. reason: sub borders only cause extra
         ! needs for memory (ghosts) and communication (ghost update)
         nsrem = nsubs

         !---------------------------------------------------------------------
         !  Factor into powers of two and remainder
         !---------------------------------------------------------------------
         ncut  = 0
         DO WHILE (MOD(nsrem,2) .EQ. 0)
            nsrem = nsrem/2
            ncut  = ncut + 1
         ENDDO
         ! number of remaining cuts
         nsrem = nsrem - 1
         IF (ppm_debug .GT. 0) THEN
            WRITE(mesg,'(A,I4)') 'pow-2 half-cuts: ',ncut
            CALL ppm_write(ppm_rank,caller,mesg,info)
            WRITE(mesg,'(A,I4)') 'remaining cuts : ',nsrem
            CALL ppm_write(ppm_rank,caller,mesg,info)
         ENDIF

         !---------------------------------------------------------------------
         !  Allocate memory for cut sizes (in units of the grid spacing)
         !---------------------------------------------------------------------
         iopt = ppm_param_alloc_fit
         ldu(1) = ppm_dim
         ldu(2) = 1
         CALL ppm_alloc(Npx,ldu,iopt,info)
         or_fail_alloc('Number of grid units NPX')

         !---------------------------------------------------------------------
         !  Cut all pieces in half, perpendicular to the longest
         !  direction, until powers of 2 are exhausted
         !---------------------------------------------------------------------
         Npx(1:ppm_dim,1)   = Nc(1:ppm_dim)
         nblocks(1:ppm_dim) = 0
         lduold             = 1
         DO icut=1,ncut
            !-----------------------------------------------------------------
            !  Determine direction of cut
            !-----------------------------------------------------------------
            IF (ppm_dim .GT. 2) THEN
               IF (constdim .EQ. 1) THEN
                  cutdim = 3
                  IF (Npx(2,1) .GT. Npx(3,1)) cutdim = 2
               ELSEIF (constdim .EQ. 2) THEN
                  cutdim = 3
                  IF (Npx(1,1) .GT. Npx(3,1)) cutdim = 1
               ELSEIF (constdim .EQ. 3) THEN
                  cutdim = 2
                  IF (Npx(1,1) .GT. Npx(2,1)) cutdim = 1
               ELSE
                  cutdim = 3
                  IF (Npx(2,1) .GT. Npx(3,1)) cutdim = 2
                  IF (Npx(1,1) .GT. Npx(cutdim,1)) cutdim = 1
               ENDIF
            ELSE
               IF (constdim .EQ. 1) THEN
                  cutdim = 2
               ELSEIF (constdim .EQ. 2) THEN
                  cutdim = 1
               ELSE
                  cutdim = 2
                  IF (Npx(1,1) .GT. Npx(2,1)) cutdim = 1
               ENDIF
            ENDIF
            !-----------------------------------------------------------------
            !  Cut in two pieces of integer multiple length of grid
            !  spacing
            !-----------------------------------------------------------------
            rc = 2**nblocks(cutdim)
            IF (MINVAL(Npx(cutdim,1:rc)) .LE. 2) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_bad_mesh,caller,&
               &    'Too little grid points for this number of subs.',&
               &    __LINE__,info)
               GOTO 9999
            ENDIF

            !-----------------------------------------------------------------
            !  Reallocate memory if too small
            !-----------------------------------------------------------------
            iopt = ppm_param_alloc_grow
            ldu(1) = ppm_dim
            IF (2*rc.GT.lduold) THEN
               ldu(2) = 2*rc
            ELSE
               ldu(2) = lduold
            ENDIF

            CALL ppm_alloc(Npxnew,ldu,iopt,info)
            or_fail_alloc('New number of grid units NPXNEW')

            !-----------------------------------------------------------------
            !  Subdivide the box in direction cutdim only
            !-----------------------------------------------------------------
            ! preserve all other dimensions
            DO k=1,lduold
               Npxnew(1:ppm_dim,k) = Npx(1:ppm_dim,k)
            ENDDO

            DO i=1,rc
               ! cut in direction cutdim
               Npxnew(cutdim,2*i-1) = CEILING(REAL(Npx(cutdim,i),MK)*0.5_MK)
               Npxnew(cutdim,2*i)   = Npx(cutdim,i)-Npxnew(cutdim,2*i-1)
            ENDDO

            !-----------------------------------------------------------------
            !  Allocate the new sub array
            !-----------------------------------------------------------------
            CALL ppm_alloc(Npx,ldu,iopt,info)
            or_fail_alloc('Number of grid units NPX')

            !-----------------------------------------------------------------
            !  Update the data at the end of the iteration
            !-----------------------------------------------------------------
            lduold = ldu(2)
            nblocks(cutdim) = nblocks(cutdim) + 1
            Npx = Npxnew
         ENDDO

         !---------------------------------------------------------------------
         !  Compute and store the number of subs in each direction
         !---------------------------------------------------------------------
         DO i=1,ppm_dim
            nblocks(i) = 2**nblocks(i)
            ! initialize surface and volume such that the
            ! surface-to-volume ratio is max.
            surface(i) = HUGE(surface(i))
            volume(i)  = 1
         ENDDO

         !---------------------------------------------------------------------
         !  Subdivide all pieces in nsrem equal parts for the
         !  remainder. Cut such as to yield mininum surface-to-volume
         !  ratio (in terms of number of grid points)
         !---------------------------------------------------------------------
         IF (nsrem .GT. 0) THEN
            ! number of cuts -> number of subdomains
            nsrem = nsrem+1
            !-----------------------------------------------------------------
            !  Determine aggregate surface-to-volume ratios for
            !  tested cutting in all three directions
            !-----------------------------------------------------------------
            IF (constdim .NE. 1) THEN
               surface(1) = 0
               volume(1)  = 0
               DO icut=1,nblocks(1)
                  rc = NINT(Npx(1,icut)/REAL(nsrem,MK))+1
                  IF (ppm_dim .GT. 2) THEN
                     ! total grid points in y and z
                     ty = nblocks(2)
                     DO k = 1,nblocks(2)
                        ty = ty + Npx(2,k)
                     ENDDO
                     tz = nblocks(3)
                     DO k = 1,nblocks(3)
                        tz = tz + Npx(3,k)
                     ENDDO
                     surface(1)=surface(1) + &
                     & 2*(nsrem*tz*ty+nblocks(3)*rc*ty+nblocks(2)*rc*tz)
                     volume(1) =volume(1)  + (rc*ty*tz)
                  ELSE
                     ty = nblocks(2)
                     DO k = 1,nblocks(2)
                        ty = ty + Npx(2,k)
                     ENDDO
                     surface(1)=surface(1) + 2*(nblocks(2)*rc+ty)
                     volume(1) =volume(1)  + (rc*ty)
                  ENDIF
               ENDDO
            ENDIF
            IF (constdim .NE. 2) THEN
               surface(2) = 0
               volume(2)  = 0
               DO jcut=1,nblocks(2)
                  rc = NINT(Npx(2,jcut)/REAL(nsrem,MK))+1
                  IF (ppm_dim .GT. 2) THEN
                     ! total grid points in x and z
                     tx = SUM(Npx(1,1:nblocks(1)))+nblocks(1)
                     tz = SUM(Npx(3,1:nblocks(3)))+nblocks(3)
                     surface(2)=surface(2) + &
                     & 2*(nsrem*tx*tz + nblocks(1)*rc*tz+nblocks(3)*rc*tx)
                     volume(2) =volume(2)  + (rc*tz*tx)
                  ELSE
                     tx = SUM(Npx(1,1:nblocks(1)))+nblocks(1)
                     surface(2)=surface(2) + 2*(nblocks(1)*rc*tx)
                     volume(2) =volume(2)  + (rc*tx)
                  ENDIF
               ENDDO
            ENDIF
            IF ((constdim .NE. 3) .AND. (ppm_dim .GT. 2)) THEN
               surface(3) = 0
               volume(3)  = 0
               DO kcut=1,nblocks(3)
                  rc = NINT(Npx(3,kcut)/REAL(nsrem,MK))+1
                  ! total grid points in x and y
                  tx = nblocks(1)
                  DO k = 1,nblocks(1)
                     tx = tx + Npx(1,k)
                  ENDDO
                  ty = nblocks(2)
                  DO k = 1,nblocks(2)
                     ty = ty + Npx(2,k)
                  ENDDO
                  surface(3)=surface(3) + &
                  & 2* (nsrem*ty*tx + nblocks(2)*rc*tx+nblocks(1)*rc*ty)
                  volume(3) =volume(3)  + (rc*tx*ty)
               ENDDO
            ENDIF
            ! NOTE: surface(i) and volume(i) are the total surface of
            ! all subs and their total volume (in number of grid
            ! POINTS) if the cutting would be done perpendicular to
            ! the i-th axis. To get the mean surface and volume of a
            ! sub, divide by the number of subs. Here this is not
            ! done because the surface-to-volume ratio is invariant
            ! to this division.

            !-----------------------------------------------------------------
            !  Find minimal mean surface-to-volume ratio
            !-----------------------------------------------------------------
            minsv = HUGE(minsv)
            DO i=1,ppm_dim
               IF (i .EQ. constdim) CYCLE
               lc = REAL(surface(i),MK)/REAL(volume(i),MK)
               IF (lc .LT. minsv) THEN
                  minsv  = lc
                  cutdim = i
               ENDIF
            ENDDO

            !-----------------------------------------------------------------
            !  Cut in nsrem pieces of integer multiple length of grid
            !  spacing
            !-----------------------------------------------------------------
            IF (MINVAL(Npx(cutdim,1:nblocks(cutdim))) .LE. nsrem) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_bad_mesh,caller,&
               &    'Too little grid points for this number of subs.',&
               &    __LINE__,info)
               GOTO 9999
            ENDIF
            rc = nsrem*nblocks(cutdim)

            !-----------------------------------------------------------------
            !  Grow the new grid point array if needed
            !-----------------------------------------------------------------
            iopt = ppm_param_alloc_grow
            ldu(1) = ppm_dim
            IF (rc.GT.lduold) THEN
               ldu(2) = rc
            ELSE
               ldu(2) = lduold
            ENDIF

            CALL ppm_alloc(Npxnew,ldu,iopt,info)
            or_fail_alloc('New number of grid points NPXNEW')

            lc = 1.0_MK/REAL(nsrem,MK)
            ! preserve all other dimensions
            Npxnew(1:ppm_dim,1:lduold) = Npx(1:ppm_dim,1:lduold)
            DO i=1,nblocks(cutdim)
               ! number of ceiling-size subs
               nup = MOD(Npx(cutdim,i),nsrem)
               ! number of floor-size subs
               ndn = nsrem-nup
               ! the two sizes
               iup = CEILING(REAL(Npx(cutdim,i),MK)*lc)
               idn = FLOOR(REAL(Npx(cutdim,i),MK)*lc)
               ! cut ceiling-size subs
               DO j=1,nup
                  Npxnew(cutdim,nsrem*(i-1)+j) = iup
               ENDDO
               ! cut floor-size subs
               DO j=(nup+1),(nup+ndn)
                  Npxnew(cutdim,nsrem*(i-1)+j) = idn
               ENDDO
               ! check if the decomposed grid cells sum up to the
               ! complete former subdomain
               IF ((nup*iup+ndn*idn) .NE. Npx(cutdim,i)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_mesh_miss,caller,          &
                  &    'Decomposed domains do not sum up to whole', &
                  &    __LINE__,info)
                  GOTO 9999
               ENDIF
            ENDDO

            CALL ppm_alloc(Npx,ldu,iopt,info)
            or_fail_alloc('Number of grid points NPX')

            Npx = Npxnew
            nblocks(cutdim) = nsrem*nblocks(cutdim)
          ENDIF

          !---------------------------------------------------------------------
          !  Deallocate temporary storage
          !---------------------------------------------------------------------
          iopt = ppm_param_dealloc
          CALL ppm_alloc(Npxnew,ldu,iopt,info)
          or_fail_alloc('New number of grid points NPXNEW')

          !---------------------------------------------------------------------
          !  Determine and store the subdomain boundaries in physical
          !  coordinates. Each subdomain will measure an integer
          !  multiple of the grid spacing in each direction.
          !---------------------------------------------------------------------
          ! The same value is only computed once and assigned to both
          ! min(i+1) and max(i) in order to allow later float
          ! comparisons.
          SELECT CASE (ppm_dim)
          CASE (2)
             nsubs = 0
             uy = min_phys(2)
             iy = 1
             ty = 0
             DO jcut=1,nblocks(2)
                ty = ty + Npx(2,jcut)
                ly = uy
                uy = min_phys(2) + ty*dx(2)
                ux = min_phys(1)
                ix = 1
                tx = 0
                DO icut=1,nblocks(1)
                   tx = tx + Npx(1,icut)
                   lx = ux
                   ux = min_phys(1) + tx*dx(1)
                   nsubs = nsubs + 1
                   ! Sub boundaries
                   min_sub(1,nsubs) = lx
                   max_sub(1,nsubs) = ux
                   min_sub(2,nsubs) = ly
                   max_sub(2,nsubs) = uy
                   ix = ix + Npx(1,icut)
                ENDDO
                iy = iy + Npx(2,jcut)
             ENDDO

          CASE (3)
             nsubs = 0
             uz = min_phys(3)
             iz = 1           ! start index in global mesh
             tz = 0           ! grid-point index in z direction
             DO kcut=1,nblocks(3)
                tz = tz + Npx(3,kcut)
                lz = uz      ! lower sub boundary
                uz = min_phys(3) + tz*dx(3)  ! upper sub boundary
                uy = min_phys(2)    ! the same in y ...
                iy = 1
                ty = 0
                DO jcut=1,nblocks(2)
                   ty = ty + Npx(2,jcut)
                   ly = uy
                   uy = min_phys(2) + ty*dx(2)
                   ux = min_phys(1)    ! and in x ...
                   ix = 1
                   tx = 0
                   DO icut=1,nblocks(1)
                      tx = tx + Npx(1,icut)
                      lx = ux
                      ux = min_phys(1) + tx*dx(1)
                      nsubs = nsubs + 1
                      ! Sub boundaries
                      min_sub(1,nsubs) = lx
                      max_sub(1,nsubs) = ux
                      min_sub(2,nsubs) = ly
                      max_sub(2,nsubs) = uy
                      min_sub(3,nsubs) = lz
                      max_sub(3,nsubs) = uz
                      ix = ix + Npx(1,icut)
                   ENDDO
                   iy = iy + Npx(2,jcut)
                ENDDO
                iz = iz + Npx(3,kcut)
             ENDDO

          END SELECT

          !---------------------------------------------------------------------
          !  Some diagnostics
          !---------------------------------------------------------------------
          IF (ppm_debug .GT. 0) THEN
             WRITE(mesg,'(A,I5)') 'number of subs created: ',nsubs
             CALL ppm_write(ppm_rank,caller,mesg,info)
          ENDIF
          IF (ppm_debug .GT. 1) THEN
             IF (ppm_dim .LT. 3) THEN
                DO i=1,nsubs
                   WRITE(mesg,'(I4,A,2F12.6)') i,' min: ',min_sub(1:2,i)
                   CALL ppm_write(ppm_rank,caller,mesg,info)
                   WRITE(mesg,'(I4,A,2F12.6)') i,' max: ',max_sub(1:2,i)
                   CALL ppm_write(ppm_rank,caller,mesg,info)
                   CALL ppm_write(ppm_rank,caller,  &
                   &    '------------------------------------',info)
                ENDDO
             ELSE
               DO i=1,nsubs
                  WRITE(mesg,'(I4,A,3F12.6)') i,' min: ',min_sub(1:3,i)
                  CALL ppm_write(ppm_rank,caller,mesg,info)
                  WRITE(mesg,'(I4,A,3F12.6)') i,' max: ',max_sub(1:3,i)
                  CALL ppm_write(ppm_rank,caller,mesg,info)
                  CALL ppm_write(ppm_rank,caller,  &
                  &    '------------------------------------------------',info)
               ENDDO
             ENDIF
          ENDIF

      !-------------------------------------------------------------------------
      !  user-provided decomposition
      !-------------------------------------------------------------------------
      CASE (ppm_param_decomp_user_defined)
          ! do not do anything as the user has already provided min_sub,
          ! max_sub. Trust the user that all arrays are properly
          ! allocated.
          !
          ! MAYBE CHECK HERE IF THE USER SUBS REALLY MATCH UP WITH THE GRID
          ! SPACING
      END SELECT !(decomp)

      !---------------------------------------------------------------------
      !  Deallocate memory
      !---------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(Npx,ldu,iopt,info)
      or_fail_dealloc('Number of grid points per sub NPX')

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
        DO i=1,ppm_dim
          IF (Nm(i) .LT. 2) THEN
             fail('Nm must be >1 in all space dimensions',exit_point=8888)
          ENDIF
          IF (max_phys(i) .LE. min_phys(i)) THEN
             fail('max_phys must be > min_phys',exit_point=8888)
          ENDIF
        ENDDO
        IF ((decomp .NE. ppm_param_decomp_xpencil) .AND.   &
        &   (decomp .NE. ppm_param_decomp_ypencil) .AND.   &
        &   (decomp .NE. ppm_param_decomp_zpencil) .AND.   &
        &   (decomp .NE. ppm_param_decomp_cuboid)  .AND.   &
        &   (decomp .NE. ppm_param_decomp_user_defined)) THEN
            fail('Invalid decomposition type specified',exit_point=8888)
        ENDIF
        IF ((decomp .EQ. ppm_param_decomp_zpencil) .AND.   &
        &   (ppm_dim .LT. 3)) THEN
           fail('zpencil decomposition is only possible for 3D problems',exit_point=8888)
        ENDIF
        IF ((decomp .EQ. ppm_param_decomp_user_defined) .AND.   &
        &   (nsubs .LE. 0)) THEN
           fail('at least one subdom. must be specified for user-defined decomposition',exit_point=8888)
        ENDIF
        IF (PRESENT(ndom)) THEN
           IF (ndom .LT. ppm_nproc) THEN
              fail('number of subs must be >= nproc',exit_point=8888)
           ENDIF
        ENDIF
      8888 CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE decomp_cart_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE decomp_cart_d
#endif

