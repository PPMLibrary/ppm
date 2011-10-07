      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_cnl_clist
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE cnl_clist_s(xp,Np,xmin,xmax,cutoff,lsymm,clist, &
     &                                 ndim,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE cnl_clist_d(xp,Np,xmin,xmax,cutoff,lsymm,clist, &
     &                                 ndim,info)
      !!! Create cell lists for all subs of this processor.
      !!!
      !!! [NOTE]
      !!! ====================================================
      !!! Symmetry is used as follows:
      !!!
      !!! ----------------------------------------------
      !!!                                   - 2 1
      !!!                                   - 0 2
      !!!                                   - - -
      !!! ----------------------------------------------
      !!! cell 0 interacts with all cells >0
      !!! cells 2 interact with each other.
      !!!
      !!! In 3D, the TOP LAYER (larger z) above the numbered
      !!! 2 x 2 block is also included. (Plus the appropriate
      !!! diagonal interactions. see MkNeighIdx).
      !!!
      !!! The cell list is created for the current topology
      !!! topoid, which is passed by the user
      !!!
      !!! Particles in cell icell of sub isub are:                             +
      !!!   `LET a = clist(isub)%lhbx(icell)`                                  +
      !!!   `LET b = clist(isub)%lhbx(icell+1)                                 +
      !!!   `clist(isub)%lpdx(a:b-1)`
      !!! ====================================================
      !!!
      !!! [WARNING]
      !!! the sub index isub is NOT the global sub
      !!! number, but just linear from 1 to ppm_nsublist!
#endif

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_typedef
      USE ppm_module_error
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN   )    :: xp
      !!! Particle co-ordinates
      INTEGER                 , INTENT(IN   )    :: Np
      !!! Number of particles
      REAL(MK), DIMENSION(ndim), INTENT(IN  ) :: xmin,xmax
      !!! actual cell size
      REAL(MK), DIMENSION(:)  , INTENT(IN   )    :: cutoff
      !!! Cutoff in all (2,3) space directions. Actual cell size may differ,
      !!! due to round-off, but it always >= cutoff.
      LOGICAL                 , INTENT(IN   )    :: lsymm
      !!! Use symmetry?
      TYPE(t_clist), DIMENSION(1)  :: clist
      !!! Number of cells in each space direction
      !!! clist(isub)%lhbx(ibox)
      INTEGER                 , INTENT(  OUT)    :: info
      !!! Returns 0 upon success
      INTEGER                 , INTENT(IN   )    :: ndim
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! timer
      REAL(MK)                                :: t0
      ! counters
      INTEGER                                 :: i,Npdx
      ! extent of cell mesh
      REAL(MK), DIMENSION(3)                  :: cellsize
      ! number of ghostlayers
      INTEGER, DIMENSION(6)                   :: Ngl
      ! parameter for alloc
      INTEGER                                 :: lda,iopt
      INTEGER, DIMENSION(2)                   :: ldc
      LOGICAL                                 :: valid


      NULLIFY(clist(1)%lpdx)
      NULLIFY(clist(1)%lhbx)
      NULLIFY(clist(1)%Nm)
      allocate(clist(1)%Nm(3))
      clist(1)%Nm(:) = 0


      !---------------------------------------------------------------------
      !  Determine number of cell boxes and effective cell size.
      !---------------------------------------------------------------------
      DO i=1,ndim
          ! number of cells based on a cellsize = cutoff 
          clist(1)%Nm(i) = INT((xmax(i) - xmin(i))/cutoff(i))
          ! make at least one box
          IF (clist(1)%Nm(i) .LT. 1) clist(1)%Nm(i) = 1
          cellsize(i) = (xmax(i) - xmin(i))/REAL(clist(1)%Nm(i),MK)
      ENDDO
      !---------------------------------------------------------------------
      !  Find out how many ghost layers are needed:
      !  Do no longer add ghost layers to the domain since this would result
      !  in errors in the ranking due to nummerical errors
      !---------------------------------------------------------------------
      Ngl(1:6) = 0
      IF (lsymm) THEN                ! EXPLOIT SYMMETRY => only need ghost 
          DO i=1,ndim             ! layers on one side
              Ngl(ndim + i) = 1
          ENDDO
      ELSE                       ! DO NOT EXPLOIT SYMMETRY => ghost layers 
          DO i=1,ndim             ! all around
             Ngl(i) = 1
             Ngl(ndim + i) = 1
          ENDDO
      ENDIF

      !---------------------------------------------------------------------
      !  Rank the particles in this extended sub
      !---------------------------------------------------------------------
      IF (ndim .EQ. 2) THEN
          CALL cnl_rank2d(xp,Np,xmin(1:2),xmax(1:2),clist(1)%Nm(1:2),&
 &                Ngl(1:4),clist(1)%lpdx,clist(1)%lhbx,info)
          !-----------------------------------------------------------------
          !  We have to increase Nm by the ghost layers to provide the same
          !  behaviour as before the change of interface of ppm_cnl_rank
          !-----------------------------------------------------------------
          clist(1)%Nm(1) = clist(1)%Nm(1) + Ngl(1) + Ngl(3)
          clist(1)%Nm(2) = clist(1)%Nm(2) + Ngl(2) + Ngl(4)
      ELSEIF (ndim .EQ. 3) THEN
          CALL cnl_rank3d(xp,Np,xmin(1:3),xmax(1:3),clist(1)%Nm(1:3),&
 &                Ngl(1:6),clist(1)%lpdx,clist(1)%lhbx,info)
          !-----------------------------------------------------------------
          !  We have to increase Nm by the ghost layers to provide the same
          !  behaviour as before the change of interface of ppm_cnl_rank
          !-----------------------------------------------------------------
          clist(1)%Nm(1) = clist(1)%Nm(1) + Ngl(1) + Ngl(4)
          clist(1)%Nm(2) = clist(1)%Nm(2) + Ngl(2) + Ngl(5)
          clist(1)%Nm(3) = clist(1)%Nm(3) + Ngl(3) + Ngl(6)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE cnl_clist_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE cnl_clist_d
#endif
