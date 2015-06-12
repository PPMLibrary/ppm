      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_find_neigh
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
      SUBROUTINE ppm_find_neigh_s(min_phys,max_phys,bcdef, &
      &          min_sub,max_sub,nsubs,nneigh,ineigh,ghostsize,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_find_neigh_d(min_phys,max_phys,bcdef, &
      &          min_sub,max_sub,nsubs,nneigh,ineigh,ghostsize,info)
      !!! This routine find the neighbours of a sub domain.
      !!! It does that by comparing the coordinates of the
      !!! bounding box of the sub domains using an N-square
      !!! algorithm! Moreover, since the routine is called before
      !!! any mapping has been performed, *all* the subs (nsubs)
      !!! must be searched!
      !!!
      !!! NOTE: No sub is listed as a neighbor of itself.
      !!!
      !!! NOTE: [SR, 02.12.2011]
      !!!     This routine has been modified to account for ghostlayer sizes
      !!!     so that a subdomain that overlaps the ghost layer of another
      !!!     subdomain is considered as a neighbor even if it does not
      !!!     touch the subdomain itself (this is what happens when 2 subs
      !!!     that are touching at a corner are slightly shifted, leading
      !!!     to some particle interactions being missed...).
      !!!
      !!! [NOTE]
      !!! ======================================================================
      !!! The side effect of this routine is that the lists
      !!! `min_sub(:)` and `max_sub(:)` are extended to include
      !!! ghost sub domains. This is not used beyond this routine.
      !!!
      !!! The lists that are built contain unique entries.
      !!! This means that a sub can occur at most once in the
      !!! neighbor list of a particular other sub. This needs
      !!! to be taken into account when doing periodic
      !!! ghosts.
      !!!
      !!! This routine uses cell lists of the sub centers
      !!! for fast search. The size of the cells is the
      !!! extent(1:3) of the largest sub that exists.
      !!! Therefore: if we have one large sub and tons of
      !!! small ones, this routine will perform poorly.
      !!! (N**2 in the worst case).
      !!! ======================================================================

#endif
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_mpi
      USE ppm_module_copy_image_subs
      USE ppm_module_util_rank
      USE ppm_module_neighlist
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
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: min_phys
      !!! Min. extent of the physical domain
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: max_phys
      !!! Max. extent of the physical domain
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: bcdef
      !!! Boundary condition definition
      REAL(MK), DIMENSION(:,:), POINTER       :: min_sub
      !!! Min. extent of the sub domain
      REAL(MK), DIMENSION(:,:), POINTER       :: max_sub
      !!! Max. extent of the sub domain
      INTEGER , DIMENSION(:,:), POINTER       :: ineigh
      !!! Points to the `nneigh(:)` subs of isub
      INTEGER , DIMENSION(:  ), POINTER       :: nneigh
      !!! `nneigh(isub)` returns the total number
      !!! of neighbouring subs of isub
      INTEGER                 , INTENT(IN   ) :: nsubs
      !!! Total number of (real) sub domains
      REAL(MK), DIMENSION(ppm_dim),INTENT(IN) :: ghostsize
      !!! The ghostlayer size in each dimension
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(ppm_dim)            :: bsize,len_sub
      REAL(MK), DIMENSION(:,:), ALLOCATABLE :: ctrs
      REAL(MK)                                :: mx1,mx2,mx3
      REAL(MK)                                :: mn1,mn2,mn3
      INTEGER , DIMENSION(:  ), POINTER       :: subid => NULL()
      INTEGER , DIMENSION(:  ), POINTER       :: lhbx  => NULL()
      INTEGER , DIMENSION(:  ), POINTER       :: lpdx  => NULL()
      INTEGER , DIMENSION(:,:), POINTER       :: inp   => NULL()
      INTEGER , DIMENSION(:,:), POINTER       :: jnp   => NULL()
      INTEGER , DIMENSION(ppm_dim)            :: ldc,Nm,Nmtot
      INTEGER , DIMENSION(2*ppm_dim)          :: Ngl
      INTEGER                                 :: nsubsplus,nnp
      INTEGER                                 :: i,j,ii,jj,k,kk,iopt,iend
      INTEGER                                 :: ipart,jpart,n1,n2,iinter
      INTEGER                                 :: isize,ip,jp,cbox,istart
      INTEGER                                 :: jstart,jend,ibox,jbox
      INTEGER                                 :: imx,jmx,kmx
      REAL(ppm_kind_double) :: t0
      LOGICAL                                 :: isin,pbdrx,pbdry,pbdrz

      CHARACTER(LEN=ppm_char) :: caller='ppm_find_neigh'
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      !-------------------------------------------------------------------------
      !  Allocate memory for subdomain IDs
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldc(1) = nsubs
      CALL ppm_alloc(subid,ldc,iopt,info)
      or_fail_alloc('Sub IDs SUBID',ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Initialize the ID
      !-------------------------------------------------------------------------
      FORALL (i=1:nsubs) subid(i) = i

      !-------------------------------------------------------------------------
      !  Add ghost domains for the periodic system
      !-------------------------------------------------------------------------
      CALL ppm_copy_image_subs(min_phys,max_phys,bcdef, &
      &    min_sub,max_sub,nsubs,subid,nsubsplus,info)
      IF (info.NE.ppm_param_success) GOTO 9999

      !-------------------------------------------------------------------------
      !  Allocate memory for the neighbours of the subdomains
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow
      ldc(1) = nsubsplus
      CALL ppm_alloc(nneigh,ldc,iopt,info)
      or_fail_alloc('Number of neighbors NNEIGH',ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Initialize the number of neighbours to zero
      !-------------------------------------------------------------------------
      FORALL (i=1:nsubsplus) nneigh(i) = 0

      !-------------------------------------------------------------------------
      !  And allocate the pointers to the neighbours
      !-------------------------------------------------------------------------
      ldc(1) = 26 ! guessing
      ldc(2) = nsubsplus
      CALL ppm_alloc(ineigh,ldc,iopt,info)
      or_fail_alloc('List of neighbors INEIGH',ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Initialize the neighbor sub indices to undefined
      !-------------------------------------------------------------------------
      FORALL (i=1:ldc(1),j=1:nsubsplus) ineigh(i,j) = ppm_param_undefined

      !-------------------------------------------------------------------------
      !  If there are less than 2 subs there can be no neighbors
      !-------------------------------------------------------------------------
      IF (nsubs .LT. 2) GOTO 8888

      !-------------------------------------------------------------------------
      !  Allocate the center point positions
      !-------------------------------------------------------------------------
      ldc(1) = ppm_dim
      ldc(2) = nsubsplus
      ALLOCATE(ctrs(ldc(1),ldc(2)),STAT=info)
      or_fail_alloc('Sub center points CTRS',ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Store the center points and the max extent of any sub
      !-------------------------------------------------------------------------
      bsize = 0.0_MK
      SELECT CASE (ppm_dim)
      CASE (3)
          DO i=1,nsubsplus
             ctrs(1,i)  = 0.5_MK*(min_sub(1,i) + max_sub(1,i))
             ctrs(2,i)  = 0.5_MK*(min_sub(2,i) + max_sub(2,i))
             ctrs(3,i)  = 0.5_MK*(min_sub(3,i) + max_sub(3,i))
             len_sub(1) = max_sub(1,i) - min_sub(1,i)
             len_sub(2) = max_sub(2,i) - min_sub(2,i)
             len_sub(3) = max_sub(3,i) - min_sub(3,i)
             IF (len_sub(1).GT.bsize(1)) bsize(1) = len_sub(1)
             IF (len_sub(2).GT.bsize(2)) bsize(2) = len_sub(2)
             IF (len_sub(3).GT.bsize(3)) bsize(3) = len_sub(3)
          ENDDO
      CASE (2)
          DO i=1,nsubsplus
             ctrs(1,i) = 0.5_MK*(min_sub(1,i) + max_sub(1,i))
             ctrs(2,i) = 0.5_MK*(min_sub(2,i) + max_sub(2,i))
             len_sub(1) = max_sub(1,i) - min_sub(1,i)
             len_sub(2) = max_sub(2,i) - min_sub(2,i)
             IF (len_sub(1).GT.bsize(1)) bsize(1) = len_sub(1)
             IF (len_sub(2).GT.bsize(2)) bsize(2) = len_sub(2)
          ENDDO
      END SELECT

      !-------------------------------------------------------------------------
      !  Determine number of cells
      !-------------------------------------------------------------------------
      DO i=1,ppm_dim
          ! number of cells based on a cellsize = cutoff
          Nm(i) = INT((max_phys(i) - min_phys(i))/bsize(i))
          ! make at least one box
          IF (Nm(i) .LT. 1) Nm(i) = 1
          bsize(i) = (max_phys(i) - min_phys(i))/REAL(Nm(i),MK)
      ENDDO

      !-------------------------------------------------------------------------
      !  Need one layer of ghost cells. Even for non-peroidic systems as
      !  otherwise the neighbor count in the top-left cells is wrong.
      !-------------------------------------------------------------------------
      IF (ppm_dim .EQ. 3) THEN
          Ngl(1:3) = 0
          Ngl(4:6) = 1
          Nmtot(1) = Nm(1) + Ngl(1) + Ngl(4)
          Nmtot(2) = Nm(2) + Ngl(2) + Ngl(5)
          Nmtot(3) = Nm(3) + Ngl(3) + Ngl(6)
      ELSE
          Ngl(1:2) = 0
          Ngl(3:4) = 1
          Nmtot(1) = Nm(1) + Ngl(1) + Ngl(3)
          Nmtot(2) = Nm(2) + Ngl(2) + Ngl(4)
      ENDIF

      !-------------------------------------------------------------------------
      !  Build the cell list
      !-------------------------------------------------------------------------
      SELECT CASE (ppm_dim)
      CASE (3)
         CALL ppm_util_rank3d(ctrs,nsubsplus, &
         &    min_phys,max_phys,Nm,Ngl,lpdx,lhbx,info)
      CASE (2)
         CALL ppm_util_rank2d(ctrs,nsubsplus, &
         &    min_phys,max_phys,Nm,Ngl,lpdx,lhbx,info)
      END SELECT
      or_fail('Failed building the cell lists.',ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Build cell-cell index lists
      !-------------------------------------------------------------------------
      CALL ppm_neighlist_MkNeighIdx(.TRUE.,inp,jnp,nnp,info)
      or_fail('Failed building the cell-cell lists.')

      !-------------------------------------------------------------------------
      !  Save the current leading dimension of ineigh
      !-------------------------------------------------------------------------
      isize = SIZE(ineigh,1)

      !-------------------------------------------------------------------------
      !  Searching the neighbours using the cell list
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  TWO DIMENSIONS
      !-------------------------------------------------------------------------
      IF (ppm_dim.EQ.2) THEN
          n1  = Nmtot(1)
          ! loop over all REAL cells (numbering is 0...n-1)
          imx = Nmtot(1)-2
          jmx = Nmtot(2)-2
          DO j=0,jmx
              !-----------------------------------------------------------------
              !  Check if we are in a boundary cell of a periodic face.
              !  This means that there could be periodic images of subs
              !  that we already have in the list and we need to check
              !  for duplicates when adding these subs.
              !-----------------------------------------------------------------
              IF ((j.EQ.jmx).AND.(bcdef(4).EQ.ppm_param_bcdef_periodic))THEN
                  pbdry = .TRUE.
              ELSE
                  pbdry = .FALSE.
              ENDIF
              DO i=0,imx
                  IF((i.EQ.imx).AND.(bcdef(2).EQ.ppm_param_bcdef_periodic))THEN
                      pbdrx = .TRUE.
                  ELSE
                      pbdrx = .FALSE.
                  ENDIF
                  ! index of the center box
                  cbox = i + 1 + n1*j
                  ! loop over all box-box interactions
                  DO iinter=1,nnp
                      ! determine box indices for this interaction
                      ibox = cbox+(inp(1,iinter)+n1*inp(2,iinter))
                      jbox = cbox+(jnp(1,iinter)+n1*jnp(2,iinter))
                      !---------------------------------------------------------
                      !  Read indices and check if box is empty
                      !---------------------------------------------------------
                      istart = lhbx(ibox)
                      iend   = lhbx(ibox+1)-1
                      IF (iend .LT. istart) CYCLE
                      !---------------------------------------------------------
                      !  Within the box itself use symmetry and avoid
                      !  adding the particle itself to its own list
                      !---------------------------------------------------------
                      IF (ibox .EQ. jbox) THEN
                          DO ipart=istart,iend
                              ip = lpdx(ipart)
                              ii = subid(ip)
#ifdef __SXF90
!CDIR NODEP
#endif
                              DO jpart=(ipart+1),iend
                                  jp = lpdx(jpart)
                                  jj = subid(jp)
#include "topo/ppm_find_neigh_subs_2d.inc"
                              ENDDO
                          ENDDO
                      !---------------------------------------------------------
                      !  For the other boxes check all particles
                      !---------------------------------------------------------
                      ELSE
                          ! get pointers to first and last particle
                          jstart = lhbx(jbox)
                          jend   = lhbx(jbox+1)-1
                          ! skip this iinter if other box is empty
                          IF (jend .LT. jstart) CYCLE
                          ! loop over all particles inside this cell
                          DO ipart=istart,iend
                              ip = lpdx(ipart)
                              ii = subid(ip)
#ifdef __SXF90
!CDIR NODEP
#endif
                              DO jpart=jstart,jend
                                  jp = lpdx(jpart)
                                  jj = subid(jp)
#include "topo/ppm_find_neigh_subs_2d.inc"
                              ENDDO
                          ENDDO
                      ENDIF       ! ibox .EQ. jbox
                  ENDDO           ! iinter
              ENDDO               ! i
          ENDDO                   ! j

      !-------------------------------------------------------------------------
      !  THREE DIMENSIONS
      !-------------------------------------------------------------------------
      ELSE
          n1  = Nmtot(1)
          n2  = Nmtot(1)*Nmtot(2)
          ! loop over all REAL cells (numbering is 0....n-1)
          imx = Nmtot(1)-2
          jmx = Nmtot(2)-2
          kmx = Nmtot(3)-2
          DO k=0,kmx
              !-----------------------------------------------------------------
              !  Check if we are in a boundary cell of a periodic face.
              !  This means that there could be periodic images of subs
              !  that we already have in the list and we need to check
              !  for duplicates when adding these subs.
              !-----------------------------------------------------------------
              IF((k.EQ.kmx).AND.(bcdef(6).EQ.ppm_param_bcdef_periodic))THEN
                  pbdrz = .TRUE.
              ELSE
                  pbdrz = .FALSE.
              ENDIF
              DO j=0,jmx
                  IF((j.EQ.jmx).AND.(bcdef(4).EQ.ppm_param_bcdef_periodic))THEN
                      pbdry = .TRUE.
                  ELSE
                      pbdry = .FALSE.
                  ENDIF
                  DO i=0,imx
                      IF((i.EQ.imx).AND.(bcdef(2).EQ.    &
     &                    ppm_param_bcdef_periodic))THEN
                          pbdrx = .TRUE.
                      ELSE
                          pbdrx = .FALSE.
                      ENDIF
                      ! index of the center box
                      cbox = i + 1 + n1*j + n2*k
                      ! loop over all box-box interactions
                      DO iinter=1,nnp
                          ! determine box indices for this interaction
                          ibox = cbox+(inp(1,iinter)+n1*inp(2,iinter)+ &
     &                           n2*inp(3,iinter))
                          jbox = cbox+(jnp(1,iinter)+n1*jnp(2,iinter)+ &
     &                           n2*jnp(3,iinter))
                          !-----------------------------------------------------
                          !  Read indices and check if box is empty
                          !-----------------------------------------------------
                          istart = lhbx(ibox)
                          iend   = lhbx(ibox+1)-1
                          IF (iend .LT. istart) CYCLE
                          !-----------------------------------------------------
                          !  Within the box itself use symmetry and avoid
                          !  adding the particle itself to its own list
                          !-----------------------------------------------------
                          IF (ibox .EQ. jbox) THEN
                              DO ipart=istart,iend
                                  ip = lpdx(ipart)
                                  ii = subid(ip)
#ifdef __SXF90
!CDIR NODEP
#endif
                                  DO jpart=(ipart+1),iend
                                      jp = lpdx(jpart)
                                      jj = subid(jp)
#include "topo/ppm_find_neigh_subs_3d.inc"
                                  ENDDO
                              ENDDO
                          !-----------------------------------------------------
                          !  For the other boxes check all particles
                          !-----------------------------------------------------
                          ELSE
                              ! get pointers to first and last particle
                              jstart = lhbx(jbox)
                              jend   = lhbx(jbox+1)-1
                              ! skip this iinter if other box is empty
                              IF (jend .LT. jstart) CYCLE
                              ! loop over all particles inside this cell
                              DO ipart=istart,iend
                                  ip = lpdx(ipart)
                                  ii = subid(ip)
#ifdef __SXF90
!CDIR NODEP
#endif
                                  DO jpart=jstart,jend
                                      jp = lpdx(jpart)
                                      jj = subid(jp)
#include "topo/ppm_find_neigh_subs_3d.inc"
                                  ENDDO
                              ENDDO
                          ENDIF       ! ibox .EQ. jbox
                      ENDDO           ! iinter
                  ENDDO               ! i
              ENDDO                   ! j
          ENDDO                       ! k
      ENDIF                           ! ppm_dim

      !-------------------------------------------------------------------------
      !  Free the memory again
      !-------------------------------------------------------------------------
      DEALLOCATE(ctrs,STAT=info)
      or_fail_dealloc('Sub center points CTRS')

      iopt = ppm_param_dealloc
      CALL ppm_alloc(lpdx,ldc,iopt,info)
      or_fail_dealloc('cell list pointers LPDX')

      CALL ppm_alloc(lhbx,ldc,iopt,info)
      or_fail_dealloc('cell list headers LHBX')

      CALL ppm_alloc(inp,ldc,iopt,info)
      or_fail_dealloc('cell-cell interactaion index INP')

      CALL ppm_alloc(jnp,ldc,iopt,info)
      or_fail_dealloc('cell-cell interaction index JNP')

      8888 iopt = ppm_param_dealloc
      CALL ppm_alloc(subid,ldc,iopt,info)
      or_fail_dealloc('Sub IDs SUBID')

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_find_neigh_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_find_neigh_d
#endif
