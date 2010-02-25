      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_decomp_cartesian
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine decomposes space using a Cartesian domain
      !                 decomposition topology. The decomposition type is given 
      !                 by decomp and can be one of:
      !                          ppm_param_decomp_xpencil
      !                          ppm_param_decomp_ypencil
      !                          ppm_param_decomp_zpencil
      !                          ppm_param_decomp_cuboid
      !                          ppm_param_decomp_user_defined
      !
      !                 In the user_defined case, the routine must be
      !                 called with a non-zero number of subdomains in
      !                 min_sub(:,:), max_sub(:,:) and this user-defined
      !                 decomposition is used unchanged. This routine is
      !                 unaware of particles and simply decomposes space
      !                 into equal grid volumes. When used for particles,
      !                 this should be followed by a check of the cost
      !                 balancing and possibly iterated with increasing
      !                 ndom.
      !                 
      !
      !  Input        : Nm(:)        (I) the number of mesh points (not
      !                                  cells) in each direction of the
      !                                  global comput. domain. (including
      !                                  those ON the boundaries)
      !                 min_phys(:)  (F) the minimum coordinate of the 
      !                                  physical/computational domain 
      !                 max_phys(:)  (F) the maximum coordinate of the 
      !                                  physical/computational domain 
      !                 decomp       (I) Type of decomposition (see above)
      !                 ndom         (I) OPTIONAL. The number of subdomains
      !                                  to be created. If not specified,
      !                                  the number of subdomains will be
      !                                  equal to the number of processors.
      !
      !  Input/output : min_sub(:,:) (F) the minimum extend of the subs
      !                 max_sub(:,:) (F) the maximum extend of the subs
      !                 nsubs        (I) the total number of subdomains
      !
      !  Output       : info         (I) return status
      !
      !  Remarks      : The computational domain is subdivided in exactly
      !                 as many subdomains as there are processors. This is
      !                 done to minimize memory usage of ghost layers
      !                 (grids need ghost layers also at intra-processor
      !                 sub boundaries!) and to get good convergence of the
      !                 grid solver (artificial internal boundary
      !                 conditionds at sub boundaries deteriorate
      !                 convergence). A quad-/oct-tree decomposition of the
      !                 computational domain is made until the remaining
      !                 number of subs does not contain any more factor
      !                 which is a power of 2. The remaining factor is used
      !                 to do a parallel slab decomposition, such that the
      !                 surface-to-volume ratio is minimized (in terms of
      !                 numbers of grid points). This way, nproc equally
      !                 sized subdomains with minimal collective surface
      !                 nodes are created.
      !
      !                 The pencil decomposition works exactly the same
      !                 except that cuts in one direction are prohibited.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_decomp_cartesian.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.16  2006/10/12 12:55:49  ivos
      !  Bugfix: in line 343, lduold was initialized to 0 instead of 1.
      !
      !  Revision 1.15  2006/07/18 09:40:42  ivos
      !  Fixed typo
      !
      !  Revision 1.14  2006/07/18 09:22:18  ivos
      !  Unrolled some loops, added comments to variables, and cleaned code.
      !
      !  Revision 1.13  2004/09/17 12:01:26  ivos
      !  fixed argument check. Nm must be .GE.2 in all directions and not
      !  only .GE.0.
      !
      !  Revision 1.12  2004/08/05 14:31:03  hiebers
      !  Fixed bug associated with the allocation of Npxnew in loops
      !
      !  Revision 1.11  2004/07/26 07:42:36  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.10  2004/03/03 09:07:52  ivos
      !  removed xp and Npart from the argument list. This routine is not aware
      !  of particles and decomposes space based on a mesh. 
      !  If one wants particles,
      !  either decomp_bisection, decomp_tree or decomp_cartesian with an
      !  iteration can be used. topo_mkfield will support all of them to get
      !  compatible particle/mesh topologies.
      !
      !  Revision 1.9  2004/02/24 15:52:54  ivos
      !  bugfix: surface-to-volume ratio was calculated wrongly for cutdim 2.
      !
      !  Revision 1.8  2004/02/11 17:13:08  ivos
      !  Bugfix: tx was initialized to 1 instead of 0 in the 2D case.
      !
      !  Revision 1.7  2004/02/09 17:40:21  ivos
      !  Implemented (and tested) cuboid decomposition.
      !
      !  Revision 1.6  2004/02/06 16:37:58  ivos
      !  Some bugs fixed. pencil decompositions are now tested to work OK.
      !
      !  Revision 1.5  2004/02/04 17:17:19  ivos
      !  Implemented pencil decompositions. Cuboids are still missing.
      !
      !  Revision 1.4  2004/01/26 13:00:13  walther
      !  Very preliminary checkin; the routine is basically missing.
      !
      !  Revision 1.3  2004/01/22 16:17:58  walther
      !  Cosmetic update.
      !
      !  Revision 1.2  2003/12/12 15:51:56  ivos
      !  Removed topoid from the argument list as it is not needed.
      !
      !  Revision 1.1.1.1  2003/11/17 15:13:45  walther
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_decomp_cartesian_s(Nm,min_phys,max_phys,decomp,    &
     &               min_sub,max_sub,nsubs,info,ndom)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_decomp_cartesian_d(Nm,min_phys,max_phys,decomp,    &
     &               min_sub,max_sub,nsubs,info,ndom)
#endif
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
      REAL(MK), DIMENSION(:,:), POINTER       :: min_sub,max_sub
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: min_phys,max_phys
      INTEGER , DIMENSION(:)  , INTENT(IN   ) :: Nm
      INTEGER                 , INTENT(INOUT) :: nsubs
      INTEGER                 , INTENT(IN   ) :: decomp
      INTEGER , OPTIONAL      , INTENT(IN   ) :: ndom
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                          :: t0,minsv,lc,lx,ux,ly,uy,lz,uz
      INTEGER                           :: i,icut,jcut,kcut,nsrem,ncut,j,k
      INTEGER                           :: cutdim,constdim,tx,ty,tz
      INTEGER                           :: ix,iy,iz,lduold,nup,ndn,iup,idn
      ! power-of-two part of the number of blocks (i.e. part for ROB)
      INTEGER                           :: rc
      ! mesh spacing
      REAL(MK), DIMENSION(ppm_dim)      :: dx
      ! physical extent of comput. domain
      REAL(MK), DIMENSION(ppm_dim)      :: len_phys
      ! number of blocks into which to subdivide the space
      INTEGER , DIMENSION(ppm_dim)      :: nblocks
      INTEGER , DIMENSION(ppm_dim)      :: surface,volume,Nc
      ! number of grid points in each sub. index: (1:ppm_dim,1:nsub)
      INTEGER , DIMENSION(:,:), POINTER :: Npx,Npxnew
      REAL(MK)                          :: gs
      CHARACTER(LEN=ppm_char)           :: mesg
      INTEGER, DIMENSION(2)             :: ldu
      INTEGER                           :: iopt
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_decomp_cartesian',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          DO i=1,ppm_dim
              IF (Nm(i) .LT. 2) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_decomp_cartesian',  &
     &                'Nm must be >1 in all space dimensions',__LINE__,info)
                  GOTO 9999
              ENDIF
              IF (max_phys(i) .LE. min_phys(i)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_decomp_cartesian',  &
     &                'max_phys must be > min_phys',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDDO
          IF ((decomp .NE. ppm_param_decomp_xpencil) .AND.   &
     &        (decomp .NE. ppm_param_decomp_ypencil) .AND.   &
     &        (decomp .NE. ppm_param_decomp_zpencil) .AND.   &
     &        (decomp .NE. ppm_param_decomp_cuboid)  .AND.   &
     &        (decomp .NE. ppm_param_decomp_user_defined)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_decomp_cartesian',  &
     &            'Invalid decomposition type specified',__LINE__,info)
              GOTO 9999
          ENDIF
          IF ((decomp .EQ. ppm_param_decomp_zpencil) .AND.   &
     &        (ppm_dim .LT. 3)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_decomp_cartesian',  &
     &            'zpencil decomposition is only possible for 3D problems', &
     &            __LINE__,info)
              GOTO 9999
          ENDIF
          IF ((decomp .EQ. ppm_param_decomp_user_defined) .AND.   &
     &        (nsubs .LE. 0)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_decomp_cartesian',  &
     &'at least one subdom. must be specified for user-defined decomposition',&
     &            __LINE__,info)
              GOTO 9999
          ENDIF
          IF (PRESENT(ndom)) THEN
              IF (ndom .LT. ppm_nproc) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_decomp_cartesian',  &
     &                'number of subs must be >= nproc',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDIF
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
      DO i=1,ppm_dim
          Nc(i)       = Nm(i)-1
          len_phys(i) = max_phys(i) - min_phys(i)
          dx(i)       = len_phys(i)/REAL(Nc(i),MK)
      ENDDO

      !-------------------------------------------------------------------------
      !  Allocate memory for the subs if not already done so by the user
      !-------------------------------------------------------------------------
      IF (decomp .NE. ppm_param_decomp_user_defined) THEN
          iopt = ppm_param_alloc_fit
          ldu(1) = ppm_dim
          ldu(2) = nsubs
          CALL ppm_alloc(min_sub,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_decomp_cartesian',    &
     &            'minima of subdomains MIN_SUB',__LINE__,info)
              GOTO 9999
          ENDIF
          CALL ppm_alloc(max_sub,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_decomp_cartesian',    &
     &            'maxima of subdomains MAX_SUB',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  pencil and cuboid decompositions
      !-------------------------------------------------------------------------
      IF     ((decomp .EQ. ppm_param_decomp_xpencil) .OR.    &
     &        (decomp .EQ. ppm_param_decomp_ypencil) .OR.    &
     &        (decomp .EQ. ppm_param_decomp_zpencil) .OR.    &
     &        (decomp .EQ. ppm_param_decomp_cuboid)) THEN
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
              CALL ppm_write(ppm_rank,'ppm_decomp_cartesian',mesg,info)
              WRITE(mesg,'(A,I4)') 'remaining cuts : ',nsrem
              CALL ppm_write(ppm_rank,'ppm_decomp_cartesian',mesg,info)
          ENDIF

          !---------------------------------------------------------------------
          !  Allocate memory for cut sizes (in units of the grid spacing)
          !---------------------------------------------------------------------
          iopt = ppm_param_alloc_fit
          ldu(1) = ppm_dim
          ldu(2) = 1
          CALL ppm_alloc(Npx,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_decomp_cartesian',&
     &            'Number of grid units NPX',&
     &            __LINE__,info)
              GOTO 9999
          ENDIF

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
                      IF (Npx(1,1) .GT. Npx(2,1)) cutdim = 1
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
                  CALL ppm_error(ppm_err_bad_mesh,'ppm_decomp_cartesian',&
     &                'Too little grid points for this number of subs.',&
     &                __LINE__,info)
                  GOTO 9999
              ENDIF

              !-----------------------------------------------------------------
              !  Reallocate memory if too small
              !-----------------------------------------------------------------
              iopt = ppm_param_alloc_grow
              ldu(1) = ppm_dim
              IF(2*rc.GT.lduold) THEN
                 ldu(2) = 2*rc
              ELSE
                 ldu(2) = lduold
              ENDIF
              CALL ppm_alloc(Npxnew,ldu,iopt,info)
              IF (info .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'ppm_decomp_cartesian',&
     &                'New number of grid units NPXNEW',&
     &                __LINE__,info)
                  GOTO 9999
              ENDIF

              !-----------------------------------------------------------------
              !  Subdivide the box in direction cutdim only
              !-----------------------------------------------------------------
              ! preserve all other dimensions
              IF (ppm_dim .EQ. 3) THEN
                  DO k=1,lduold
                      Npxnew(1,k) = Npx(1,k)
                      Npxnew(2,k) = Npx(2,k)
                      Npxnew(3,k) = Npx(3,k)
                  ENDDO
              ELSE
                  DO k=1,lduold
                      Npxnew(1,k) = Npx(1,k)
                      Npxnew(2,k) = Npx(2,k)
                  ENDDO
              ENDIF

              DO i=1,rc
                  ! cut in direction cutdim
                  Npxnew(cutdim,2*i-1) =               &
     &                CEILING(REAL(Npx(cutdim,i),MK)*0.5_MK)
                  Npxnew(cutdim,2*i) = Npx(cutdim,i)-Npxnew(cutdim,2*i-1)
              ENDDO

              !-----------------------------------------------------------------
              !  Allocate the new sub array
              !-----------------------------------------------------------------
              CALL ppm_alloc(Npx,ldu,iopt,info)
              IF (info .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'ppm_decomp_cartesian',&
     &                'Number of grid units NPX',&
     &                __LINE__,info)
                  GOTO 9999
              ENDIF

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
                          surface(1)=surface(1) + 2*(nsrem*tz*ty+  &
     &                        nblocks(3)*rc*ty+nblocks(2)*rc*tz)
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
                          surface(2)=surface(2) + 2*(nsrem*tx*tz+  &
     &                        nblocks(1)*rc*tz+nblocks(3)*rc*tx)
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
                      surface(3)=surface(3) + 2* (nsrem*ty*tx+  &
     &                    nblocks(2)*rc*tx+nblocks(1)*rc*ty)
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
                  CALL ppm_error(ppm_err_bad_mesh,'ppm_decomp_cartesian',&
     &                'Too little grid points for this number of subs.',&
     &                __LINE__,info)
                  GOTO 9999
              ENDIF
              rc = nsrem*nblocks(cutdim)
              
              !-----------------------------------------------------------------
              !  Grow the new grid point array if needed
              !-----------------------------------------------------------------
              iopt = ppm_param_alloc_grow
              ldu(1) = ppm_dim
              IF(rc.GT.lduold) THEN
                 ldu(2) = rc
              ELSE
                 ldu(2) = lduold
              ENDIF

              CALL ppm_alloc(Npxnew,ldu,iopt,info)
              IF (info .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'ppm_decomp_cartesian',&
     &                'New number of grid points NPXNEW',&
     &                __LINE__,info)
                  GOTO 9999
              ENDIF
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
                      CALL ppm_error(ppm_err_mesh_miss,   &
     &                    'ppm_decomp_cartesian',         &
     &                    'Decomposed domains do not sum up to whole',  &
     &                    __LINE__,info)
                      GOTO 9999
                  ENDIF
              ENDDO
              CALL ppm_alloc(Npx,ldu,iopt,info)
              IF (info .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'ppm_decomp_cartesian',&
     &                'Number of grid points NPX',&
     &                __LINE__,info)
                  GOTO 9999
              ENDIF
              Npx = Npxnew
              nblocks(cutdim) = nsrem*nblocks(cutdim)
          ENDIF

          !---------------------------------------------------------------------
          !  Deallocate temporary storage
          !---------------------------------------------------------------------
          iopt = ppm_param_dealloc
          CALL ppm_alloc(Npxnew,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_decomp_cartesian',&
     &            'New number of grid points NPXNEW',&
     &            __LINE__,info)
          ENDIF

          !---------------------------------------------------------------------
          !  Determine and store the subdomain boundaries in physical
          !  coordinates. Each subdomain will measure an integer
          !  multiple of the grid spacing in each direction.
          !---------------------------------------------------------------------
          ! The same value is only computed once and assigned to both
          ! min(i+1) and max(i) in order to allow later float
          ! comparisons.
          IF (ppm_dim .EQ. 3) THEN
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
          ELSE
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
          ENDIF

          !---------------------------------------------------------------------
          !  Some diagnostics
          !---------------------------------------------------------------------
          IF (ppm_debug .GT. 0) THEN
              WRITE(mesg,'(A,I5)') 'number of subs created: ',nsubs
              CALL ppm_write(ppm_rank,'ppm_decomp_cartesian',mesg,info)
          ENDIF
          IF (ppm_debug .GT. 1) THEN
              IF (ppm_dim .LT. 3) THEN
                  DO i=1,nsubs
                      WRITE(mesg,'(I4,A,2F12.6)') i,' min: ',min_sub(1:2,i)
                      CALL ppm_write(ppm_rank,'ppm_decomp_cartesian',  &
     &                    mesg,info)
                      WRITE(mesg,'(I4,A,2F12.6)') i,' max: ',max_sub(1:2,i)
                      CALL ppm_write(ppm_rank,'ppm_decomp_cartesian',  &
     &                    mesg,info)
                      CALL ppm_write(ppm_rank,'ppm_decomp_cartesian',  &
     &                    '------------------------------------',info)
                  ENDDO
              ELSE
                  DO i=1,nsubs
                      WRITE(mesg,'(I4,A,3F12.6)') i,' min: ',min_sub(1:3,i)
                      CALL ppm_write(ppm_rank,'ppm_decomp_cartesian',  &
     &                    mesg,info)
                      WRITE(mesg,'(I4,A,3F12.6)') i,' max: ',max_sub(1:3,i)
                      CALL ppm_write(ppm_rank,'ppm_decomp_cartesian',  &
     &                    mesg,info)
                      CALL ppm_write(ppm_rank,'ppm_decomp_cartesian',  &
     &             '------------------------------------------------',info)
                  ENDDO
              ENDIF
          ENDIF

      !-------------------------------------------------------------------------
      !  user-provided decomposition
      !-------------------------------------------------------------------------
      ELSEIF (decomp .EQ. ppm_param_decomp_user_defined) THEN
          ! do not do anything as the user has already provided min_sub,
          ! max_sub. Trust the user that all arrays are properly
          ! allocated.
          !
          ! MAYBE CHECK HERE IF THE USER SUBS REALLY MATCH UP WITH THE GRID
          ! SPACING
      ENDIF

      !---------------------------------------------------------------------
      !  Deallocate memory
      !---------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(Npx,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_decomp_cartesian',&
     &        'Number of grid points per sub NPX',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_decomp_cartesian',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_decomp_cartesian_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_decomp_cartesian_d
#endif
