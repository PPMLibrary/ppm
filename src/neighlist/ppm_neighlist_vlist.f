      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_neighlist_vlist
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
      SUBROUTINE ppm_neighlist_vlist_s(topoid,xp,np,cutoff,skin,lsymm,vlist, &
     &               nvlist,info,pidx,clist,lstore)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_neighlist_vlist_d(topoid,xp,np,cutoff,skin,lsymm,vlist, &
     &               nvlist,info,pidx,clist,lstore)
      !!! Create Verlet lists for all particles of this processor.
      !!!
      !!! TIP: Ghostparticles must be included when passing the positions 
      !!! `xp` and the array size `np` to generate Verlet lists for real/ghost
      !!! interactions.
      !!!
      !!! [NOTE]
      !!! ====================================================
      !!! The list needs to be rebuilt as soon as a particle
      !!! has moved a distance larger than `0.5*skin`. It is the
      !!! *users* responsibility to detect when this is the
      !!! case and *call* this routine again.
      !!!
      !!! vlist and nvlist are allocated in this routine. The
      !!! user just needs to pass pointers.
      !!!
      !!! the two cases for lsymm have their own duplicated
      !!! loops since the lsymm=F case does not vectorize.
      !!! lsymm=T (using symmetry) however does.
      !!! ====================================================
      !!!
      !!! NOTE: The VECTOR case was tested and found to vectorize
      !!! on the NEC SX-5 even without compiler directives.
      !!! Requires (almost) two repetitions of the main loops.
#endif

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
      USE ppm_module_check_id
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
      !!! particle co-ordinates
      INTEGER                 , INTENT(IN   ) :: np
      !!! number of particles.
      !!! The number of ghostparticles should be included to include
      !!! interactions between real and ghost-particles.
      INTEGER                 , INTENT(IN   ) :: topoid
      !!! ID of current topology
      REAL(MK)                , INTENT(IN   ) :: cutoff
      !!! cutoff radius for PP interactions
      REAL(MK)                , INTENT(IN   ) :: skin
      !!! Verlet list skin layer thickness.
      LOGICAL                 , INTENT(IN   ) :: lsymm
      !!! Use symmetry
      INTEGER, DIMENSION(:,:) , POINTER       :: vlist
      !!! Verlet list. First index: particles with which particle ip interacts.
      !!! Second index: ip. The second index only runs up to the
      !!! largest ip with non-zero nvlist. This is to save memory since the
      !!! last particles are the ghosts and they do not have a verlet list.
      !!! This is only allocated and returned if lstore is .TRUE.
      INTEGER, DIMENSION(  :) , POINTER       :: nvlist
      !!! Number of particles with which ip has to interact. Index: ip.
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      INTEGER, DIMENSION(  :) , OPTIONAL               :: pidx
      !!! OPTIONAL indices of those particles that are to be included in the
      !!! list. By default all particles are taken. If given, particles
      !!! indices in Verlet lists are relative to xp(:,pidx) and not xp(:,:)
      TYPE(ppm_t_clist), DIMENSION(:),POINTER,OPTIONAL :: clist
      !!! Cell list data structure. Pass this argument as null to force
      !!! this routine to recreate a cell list and store it in clist. Otherwise,
      !!! the cell list in clist is (re)used for the vlist being created.
      !!! PPM will use internal data structures to store the clist if this
      !!! argument is not passed.
      !!!
      !!! NOTE: use ppm_destroy_clist to deallocate the cell list.
      LOGICAL, INTENT(IN)     , OPTIONAL               :: lstore
      !!! OPTIONAL Set this to .TRUE. to store (and return) the Verlet lists in
      !!! vlist. If this is false, only nvlist is determined and returned.
      !!! Default is .TRUE.
      
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! timer
      REAL(MK)                                   :: t0
      ! effective number of particles
      REAL(MK), DIMENSION(ppm_dim)               :: min_phys,max_phys
      ! domain extents
      REAL(MK), DIMENSION(ppm_dim)               :: xmin,xmax
      ! subdomain extents 
      INTEGER                                    :: npdx
      ! counters
      INTEGER                                    :: i,idom,ibox,jbox,nbox
      INTEGER                                    :: ipart,jpart,ip,jp,maxvlen
      INTEGER                                    :: cbox,iinter,j,k
      INTEGER                                    :: ii,jj,kk
      ! coordinate difference
      REAL(MK)                                   :: dx,dy,dz
      ! inter particle distance
      REAL(MK)                                   :: dij
      ! cutoff squared
      REAL(MK)                                   :: cut2
      ! start and end particle in a box
      INTEGER                                    :: istart,iend,jstart,jend
      ! box size for helper cell list
      REAL(MK), DIMENSION(3)                     :: bsize
      ! cell neighbor lists
      INTEGER, DIMENSION(:,:), POINTER           :: inp => NULL()
      INTEGER, DIMENSION(:,:), POINTER           :: jnp => NULL()
      ! number of interactions for each cell
      INTEGER                                    :: nnp
      ! for allocate
      INTEGER, DIMENSION(2)                      :: lda
      INTEGER                                    :: iopt
      ! number of cells in all directions
      INTEGER                                    :: n1,n2,nz
      INTEGER, DIMENSION(3)                      :: lb
      INTEGER                                    :: nsbc
      LOGICAL, DIMENSION(2*ppm_dim)              :: isbc
      CHARACTER(LEN=ppm_char)                    :: mesg
      TYPE(ppm_t_clist), DIMENSION(:),POINTER    :: cl => NULL()
      ! store vlist?
      LOGICAL                                    :: lst
      LOGICAL                                    :: valid
      TYPE(ppm_t_topo)       , POINTER           :: topo => NULL()
      REAL(MK)                                   :: eps
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_neighlist_vlist',t0,info)
      
      !-------------------------------------------------------------------------
      !  If the user gave an explicit list of particles to be included, use
      !  the size of this list as the effective number of particles. Use
      !  np otherwise.
      !-------------------------------------------------------------------------
      IF (PRESENT(pidx)) THEN
          IF (np .GT. SIZE(pidx,1)) npdx = SIZE(pidx,1)
      ELSE
          npdx = np
      ENDIF
      !-------------------------------------------------------------------------
      !  Do we need to store the Verlet lists or just determine their lengths?
      !-------------------------------------------------------------------------
      IF (PRESENT(lstore)) THEN
          lst = lstore
      ELSE
          lst = .TRUE.
      ENDIF

      topo => ppm_topo(topoid)%t
      
      !-------------------------------------------------------------------------
      !  Check Arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

#if   __KIND == __DOUBLE_PRECISION
      eps = ppm_myepsd
#elif __KIND == __SINGLE_PRECISION
      eps = ppm_myepss
#endif

      
      !-------------------------------------------------------------------------
      ! Determine if there are any (non-)symmetric boundary conditions
      !-------------------------------------------------------------------------
      nsbc = 0
      isbc(:) = .FALSE. 
      DO i=1,2*ppm_dim
          SELECT CASE (topo%bcdef(i))
          CASE (ppm_param_bcdef_symmetry)
              nsbc = nsbc + 1
              isbc(i) = .TRUE.
          CASE (ppm_param_bcdef_antisymmetry)
              nsbc = nsbc + 1
              isbc(i) = .TRUE.
          CASE (ppm_param_bcdef_neumann)
              nsbc = nsbc + 1
              isbc(i) = .TRUE.
          CASE (ppm_param_bcdef_dirichlet)
              nsbc = nsbc + 1
              isbc(i) = .TRUE.
          END SELECT
      ENDDO

#if   __KIND == __DOUBLE_PRECISION
      min_phys(:) = topo%min_physd
      max_phys(:) = topo%max_physd
#elif __KIND == __SINGLE_PRECISION
      min_phys(:) = topo%min_physs
      max_phys(:) = topo%max_physs
#endif

      !-------------------------------------------------------------------------
      !  Boxes need to be cutoff+skin in all directions !
      !-------------------------------------------------------------------------
      DO i=1,ppm_dim
          bsize(i) = cutoff + skin
      ENDDO
      cut2 = bsize(1)*bsize(1)
      
      !-------------------------------------------------------------------------
      !  Generate cell lists 
      !  Check if, user is providing a cell list, to skip this step
      !-------------------------------------------------------------------------
      IF (PRESENT(clist)) THEN
          cl => clist
      ELSE
          cl => ppm_clist
      ENDIF
      IF (.NOT.(PRESENT(clist).AND.ASSOCIATED(clist))) THEN
          IF (PRESENT(pidx)) THEN
              CALL ppm_neighlist_clist(topoid,xp(:,pidx),npdx,bsize, &
     &                                 lsymm,cl,info)
          ELSE
              CALL ppm_neighlist_clist(topoid,xp,npdx,bsize,lsymm,cl,info)
          ENDIF
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_sub_failed,'ppm_neighlist_vlist',   &
     &             'Building cell lists failed.',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF
      !-------------------------------------------------------------------------
      !  Generate cell neighbor lists 
      !-------------------------------------------------------------------------
      CALL ppm_neighlist_MkNeighIdx(lsymm,inp,jnp,nnp,info)
      IF (info.NE.0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_neighlist_vlist',  &
     &        'neighlist_MkNeighIdx failed',__LINE__,info)
          GOTO 9999
      ENDIF
      !-------------------------------------------------------------------------
      !  Allocate memory for Verlet list lengths
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_grow
      lda(1) = npdx
      CALL ppm_alloc(nvlist,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_neighlist_vlist',  &
     &         'Verlet list sizes NVLIST',__LINE__,info)
          GOTO 9999
      ENDIF
      nvlist(1:npdx) = 0


      !-------------------------------------------------------------------------
      !  Determine size of Verlet lists
      !-------------------------------------------------------------------------
      DO idom=1,topo%nsublist
              !-----------------------------------------------------------------
              !  Copy subdomain extent into precision agnostic variables 
              !-----------------------------------------------------------------
#if   __KIND == __DOUBLE_PRECISION
              xmin(1) = topo%min_subd(1,idom)
              xmax(1) = topo%max_subd(1,idom)
              xmin(2) = topo%min_subd(2,idom)
              xmax(2) = topo%max_subd(2,idom)
          IF (ppm_dim .EQ. 3) THEN
              xmin(3) = topo%min_subd(3,idom)
              xmax(3) = topo%max_subd(3,idom)
          ENDIF
#elif __KIND == __SINGLE_PRECISION
              xmin(1) = topo%min_subs(1,idom)
              xmax(1) = topo%max_subs(1,idom)
              xmin(2) = topo%min_subs(2,idom)
              xmax(2) = topo%max_subs(2,idom)
          IF (ppm_dim .EQ. 3) THEN
              xmin(3) = topo%min_subs(3,idom)
              xmax(3) = topo%max_subs(3,idom)
          ENDIF
#endif
          !---------------------------------------------------------------------
          !  Lower box bound depends on symmetry and boundary condition
          !---------------------------------------------------------------------
          IF (.NOT.lsymm) THEN
              DO i=1,ppm_dim
                  IF (ABS(xmin(i)-min_phys(i)).LT.eps .AND.isbc(i)) THEN 
                      lb(i) = 0
                  ELSE
                      lb(i) = 1
                  ENDIF
              ENDDO
          ELSE
              lb(:) = 0
          ENDIF
          n1  = cl(idom)%nm(1)
          n2  = cl(idom)%nm(1)*cl(idom)%nm(2)
          IF (ppm_dim.EQ.3) THEN
              nz  = cl(idom)%nm(3)
          ELSE IF (ppm_dim .EQ. 2) THEN
              n2 = 0
              nz = lb(3)+2
          ENDIF 
          ! loop over all REAL cells (the -2 at the end does this)
          DO k=lb(3),nz-2
              DO j=lb(2),cl(idom)%nm(2)-2
                  DO i=lb(1),cl(idom)%nm(1)-2
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
                          !  Read indices and check if empty
                          !-----------------------------------------------------
                          istart = cl(idom)%lhbx(ibox)
                          iend   = cl(idom)%lhbx(ibox+1)-1
                          IF (iend .LT. istart) CYCLE
                          !-----------------------------------------------------
                          !  Within the box itself use symmetry and avoid 
                          !  adding the particle itself to its own list
                          !-----------------------------------------------------
                          IF (ibox .EQ. jbox) THEN
                              DO ipart=istart,iend
                                  ip = cl(idom)%lpdx(ipart)
                                  IF (lsymm) THEN
                                      DO jpart=(ipart+1),iend
                                          jp = cl(idom)%lpdx(jpart)
                                          ! translate to real particle
                                          ! index if needed
                                          IF (PRESENT(pidx)) THEN
                                              ii = pidx(ip)
                                              jj = pidx(jp)
                                          ELSE
                                              ii = ip
                                              jj = jp
                                          ENDIF
                                          dx = xp(1,ii) - xp(1,jj)
                                          dy = xp(2,ii) - xp(2,jj)
                                          IF (ppm_dim .GT. 2) THEN
                                              dz = xp(3,ii) - xp(3,jj)
                                              dij= (dx*dx)+(dy*dy)+(dz*dz)
                                          ELSE
                                              dz = 0.0_MK
                                              dij= (dx*dx)+(dy*dy)
                                          ENDIF
                                          IF (dij .GT. cut2) CYCLE
                                          ! add particle jp to list of
                                          ! particle ip
                                          nvlist(ip) = nvlist(ip) + 1
                                      ENDDO
                                  ELSE
#ifdef __VECTOR
                                      DO jpart=istart,iend
                                          jp = cl(idom)%lpdx(jpart)
                                          IF (jp .EQ. ip) CYCLE
#else
                                      DO jpart=(ipart+1),iend
                                          jp = cl(idom)%lpdx(jpart)
#endif
                                          ! translate to real particle
                                          ! index if needed
                                          IF (PRESENT(pidx)) THEN
                                              ii = pidx(ip)
                                              jj = pidx(jp)
                                          ELSE
                                              ii = ip
                                              jj = jp
                                          ENDIF
                                          dx = xp(1,ii) - xp(1,jj)
                                          dy = xp(2,ii) - xp(2,jj)
                                          IF (ppm_dim .GT. 2) THEN
                                              dz = xp(3,ii) - xp(3,jj)
                                              dij= (dx*dx)+(dy*dy)+(dz*dz)
                                          ELSE
                                              dz = 0.0_MK
                                              dij= (dx*dx)+(dy*dy)
                                          ENDIF
                                          IF (dij .GT. cut2) CYCLE
                                          ! add particle jp to list of
                                          ! particle ip
                                          nvlist(ip) = nvlist(ip) + 1
#ifndef __VECTOR
                                          nvlist(jp) = nvlist(jp) + 1
#endif
                                      ENDDO
                                  ENDIF
                              ENDDO
                          !-----------------------------------------------------
                          !  For the other boxes check all particles
                          !-----------------------------------------------------
                          ELSE
                              ! get pointers to first and last particle 
                              jstart = cl(idom)%lhbx(jbox)
                              jend   = cl(idom)%lhbx(jbox+1)-1
                              ! skip this iinter if empty
                              If (jend .LT. jstart) CYCLE
                              ! loop over all particles inside this cell
                              DO ipart=istart,iend
                                  ip = cl(idom)%lpdx(ipart)
                                  ! check against all particles 
                                  ! in the other cell
                                  DO jpart=jstart,jend
                                      jp = cl(idom)%lpdx(jpart)
                                      ! translate to real particle
                                      ! index if needed
                                      IF (PRESENT(pidx)) THEN
                                          ii = pidx(ip)
                                          jj = pidx(jp)
                                      ELSE
                                          ii = ip
                                          jj = jp
                                      ENDIF
                                      dx = xp(1,ii) - xp(1,jj)
                                      dy = xp(2,ii) - xp(2,jj)
                                      IF (ppm_dim .GT. 2) THEN
                                          dz = xp(3,ii) - xp(3,jj)
                                          dij= (dx*dx)+(dy*dy)+(dz*dz)
                                      ELSE
                                          dz = 0.0_MK
                                          dij= (dx*dx)+(dy*dy)
                                      ENDIF
                                      IF (dij .GT. cut2) CYCLE
                                      ! add particle jp to Verlet 
                                      ! list of particle ip
                                      nvlist(ip) = nvlist(ip) + 1
                                  ENDDO
                              ENDDO
                          ENDIF       ! ibox .EQ. jbox
                      ENDDO           ! iinter
                  ENDDO               ! i
              ENDDO                   ! j
          ENDDO                       ! k
      ENDDO                           ! idom

      !-------------------------------------------------------------------------
      !  Maximum Verlet list length
      !-------------------------------------------------------------------------
      maxvlen = MAXVAL(nvlist)
      IF (ppm_debug .GT. 0) THEN
          WRITE(mesg,'(A,I8)') 'Maximum length of Verlet lists: ',maxvlen
          CALL ppm_write(ppm_rank,'ppm_neighlist_vlist',mesg,info)
      ENDIF
      !-------------------------------------------------------------------------
      !  Only build and store the lists if needed
      !-------------------------------------------------------------------------
      IF (lst) THEN
          !---------------------------------------------------------------------
          !  Highest particle index with non-zero nvlist
          !---------------------------------------------------------------------
          ii = npdx
          DO ipart=npdx,1,-1
              IF (nvlist(ipart) .GT. 0) THEN
                  ii = ipart
                  EXIT
              ENDIF
          ENDDO
          !---------------------------------------------------------------------
          !  Allocate memory for Verlet lists
          !---------------------------------------------------------------------
          iopt = ppm_param_alloc_grow
          lda(1) = maxvlen
          lda(2) = ii
          CALL ppm_alloc(vlist,lda,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_neighlist_vlist',  &
         &         'Verlet list VLIST',__LINE__,info)
              GOTO 9999
          ENDIF
          nvlist(1:npdx) = 0
          !---------------------------------------------------------------------
          !  BUILD VERLET LISTS 
          !---------------------------------------------------------------------
          DO idom=1,topo%nsublist
              n1  = cl(idom)%nm(1)
              n2  = cl(idom)%nm(1)*cl(idom)%nm(2)
              IF (ppm_dim.EQ.3) THEN
                  nz  = cl(idom)%nm(3)
              ELSE IF (ppm_dim .EQ. 2) THEN
                  n2 = 0
                  nz = lb(3)+2
              ENDIF 
              ! get number of cells in this subdomain
              nbox = SIZE(cl(idom)%lhbx,1)-1
              ! loop over all REAL cells (the -2 at the end does this)
              DO k=lb(3),nz-2
                  DO j=lb(2),cl(idom)%nm(2)-2
                      DO i=lb(1),cl(idom)%nm(1)-2
                          ! index of the center box
                          cbox = i + 1 + n1*j + n2*k
                          ! loop over all box-box interactions
                          DO iinter=1,nnp
                              ! determine box indices for this interaction
                              ibox = cbox+(inp(1,iinter)+n1*inp(2,iinter)+ &
     &                           n2*inp(3,iinter))
                              jbox = cbox+(jnp(1,iinter)+n1*jnp(2,iinter)+ &
     &                           n2*jnp(3,iinter))
                              !-------------------------------------------------
                              !  Read indices and check if empty
                              !-------------------------------------------------
                              istart = cl(idom)%lhbx(ibox)
                              iend   = cl(idom)%lhbx(ibox+1)-1
                              IF (iend .LT. istart) CYCLE
                              !-------------------------------------------------
                              !  Within the box itself use symmetry and avoid 
                              !  adding the particle itself to its own list
                              !-------------------------------------------------
                              IF (ibox .EQ. jbox) THEN
                                  DO ipart=istart,iend
                                      ip = cl(idom)%lpdx(ipart)
                                      kk = nvlist(ip)
                                      IF (lsymm) THEN
                                          DO jpart=(ipart+1),iend
                                              jp = cl(idom)%lpdx(jpart)
                                              ! translate to real particle
                                              ! index if needed
                                              IF (PRESENT(pidx)) THEN
                                                  ii = pidx(ip)
                                                  jj = pidx(jp)
                                              ELSE
                                                  ii = ip
                                                  jj = jp
                                              ENDIF
                                              dx = xp(1,ii) - xp(1,jj)
                                              dy = xp(2,ii) - xp(2,jj)
                                              IF (ppm_dim .GT. 2) THEN
                                                  dz = xp(3,ii) - xp(3,jj)
                                                  dij= (dx*dx)+(dy*dy)+(dz*dz)
                                              ELSE
                                                  dz = 0.0_MK
                                                  dij= (dx*dx)+(dy*dy)
                                              ENDIF
                                              IF (dij .GT. cut2) CYCLE
                                              ! add particle jp to list of
                                              ! particle ip
                                              kk = kk + 1
                                              vlist(kk,ip) = jp
                                          ENDDO
                                      ELSE
#ifdef __VECTOR
                                          DO jpart=istart,iend
                                              jp = cl(idom)%lpdx(jpart)
                                              IF (jp .EQ. ip) CYCLE
#else
                                          DO jpart=(ipart+1),iend
                                              jp = cl(idom)%lpdx(jpart)
#endif
                                              ! translate to real particle
                                              ! index if needed
                                              IF (PRESENT(pidx)) THEN
                                                  ii = pidx(ip)
                                                  jj = pidx(jp)
                                              ELSE
                                                  ii = ip
                                                  jj = jp
                                              ENDIF
                                              dx = xp(1,ii) - xp(1,jj)
                                              dy = xp(2,ii) - xp(2,jj)
                                              IF (ppm_dim .GT. 2) THEN
                                                  dz = xp(3,ii) - xp(3,jj)
                                                  dij= (dx*dx)+(dy*dy)+(dz*dz)
                                              ELSE
                                                  dz = 0.0_MK
                                                  dij= (dx*dx)+(dy*dy)
                                              ENDIF
                                              IF (dij .GT. cut2) CYCLE
                                              kk = kk + 1
                                              ! add particle jp to list of
                                              ! particle ip
                                              vlist(kk,ip) = jp
#ifndef __VECTOR
                                              ! add particle ip to list of
                                              ! particle jp
                                              nvlist(jp) = nvlist(jp) + 1
                                              vlist(nvlist(jp),jp) = ip
#endif
                                          ENDDO
                                      ENDIF
                                      nvlist(ip) = kk
                                  ENDDO
                              !-------------------------------------------------
                              !  For the other boxes check all particles
                              !-------------------------------------------------
                              ELSE
                                  ! get pointers to first and last particle 
                                  jstart = cl(idom)%lhbx(jbox)
                                  jend   = cl(idom)%lhbx(jbox+1)-1
                                  ! skip this iinter if empty
                                  IF (jend .LT. jstart) CYCLE
                                  ! loop over all particles inside this cell
                                  DO ipart=istart,iend
                                      ip = cl(idom)%lpdx(ipart)
                                      kk = nvlist(ip)
                                      ! check against all particles 
                                      ! in the other cell
                                      DO jpart=jstart,jend
                                          jp = cl(idom)%lpdx(jpart)
                                          ! translate to real particle
                                          ! index if needed
                                          IF (PRESENT(pidx)) THEN
                                              ii = pidx(ip)
                                              jj = pidx(jp)
                                          ELSE
                                              ii = ip
                                              jj = jp
                                          ENDIF
                                          dx = xp(1,ii) - xp(1,jj)
                                          dy = xp(2,ii) - xp(2,jj)
                                          IF (ppm_dim .GT. 2) THEN
                                              dz = xp(3,ii) - xp(3,jj)
                                              dij= (dx*dx)+(dy*dy)+(dz*dz)
                                          ELSE
                                              dz = 0.0_MK
                                              dij= (dx*dx)+(dy*dy)
                                          ENDIF
                                          IF (dij .GT. cut2) CYCLE
                                          ! add particle jp to Verlet 
                                          ! list of particle ip
                                          kk = kk + 1
                                          vlist(kk,ip) = jp
                                      ENDDO
                                      nvlist(ip) = kk
                                  ENDDO
                              ENDIF       ! ibox .EQ. jbox
                          ENDDO           ! iinter
                      ENDDO               ! i
                  ENDDO                   ! j
              ENDDO                       ! k
          ENDDO                           ! idom
      ENDIF       ! lstore

      !-------------------------------------------------------------------------
      !  Free work cell lists
      !-------------------------------------------------------------------------
      ! JHW testing performance 20061109
      ! CALL ppm_clist_destroy(clist,info)
      ! IF (info .NE. 0) THEN
      !     info = ppm_error_error
      !     CALL ppm_error(ppm_err_dealloc,'ppm_neighlist_vlist',  &
      !&         'Failed to destroy cell lists',__LINE__,info)
      ! ENDIF
      iopt = ppm_param_dealloc
      CALL ppm_alloc(inp,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_neighlist_vlist',  &
     &         'Box interaction index INP',__LINE__,info)
      ENDIF
      CALL ppm_alloc(jnp,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_neighlist_vlist',  &
     &         'Box interaction index JNP',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_neighlist_vlist',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_neighlist_vlist',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (cutoff .LE. 0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_neighlist_vlist',  &
     &            'cutoff must be >0',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (skin .LT. 0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_neighlist_vlist',  &
     &            'skin must be >= 0',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (npdx .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_neighlist_vlist',  &
     &            'np must be >0',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (topoid .EQ. ppm_param_topo_undefined) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_neighlist_vlist',  &
     &            'Geometric topology required',__LINE__,info)
                  GOTO 8888
          ENDIF
          IF (topoid .NE. ppm_param_topo_undefined) THEN
              CALL ppm_check_topoid(topoid,valid,info)
              IF (.NOT. valid) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_neighlist_vlist',  &
     &                 'topoid out of range',__LINE__,info)
                  GOTO 8888
              ENDIF
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_neighlist_vlist_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_neighlist_vlist_d
#endif
