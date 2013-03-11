      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_topo_store
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
      SUBROUTINE ppm_topo_store_s(topoid,min_phys,max_phys,min_sub,max_sub, &
     &                            subs_bc,sub2proc,nsubs,bcdef,ghostsize,&
     &                            isublist,nsublist,nneigh,ineigh,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_topo_store_d(topoid,min_phys,max_phys,min_sub,max_sub, &
     &                            subs_bc,sub2proc,nsubs,bcdef,ghostsize,&
     &                            isublist,nsublist,nneigh,ineigh,info)
#endif
      !!! This routine stores all relevant information about
      !!! the new topology (created by ppm_decomp routines) in
      !!! a topology structure. If the given topology ID points to
      !!! an existing topology, it is replaced, else if ID == 0, a new
      !!! topology structure is allocated and appended into ppm_topo. The
      !!! new ID is returned for the user.
      !!!
      !!! [NOTE]
      !!! This routine dellocates the array: subs_bc passed to it.
      !!! The choice of placing it here is somewhat arbitrary, but
      !!! since this routine is call BOTH in ppm_topo_mkpart and
      !!! ppm_topo_mkfield - why not do it in one place - here.

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_topo_alloc
      USE ppm_module_topo_typedef
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:  ), INTENT(IN   ) :: min_phys
      !!! Min. extent of the comput. domain
      REAL(MK), DIMENSION(:  ), INTENT(IN   ) :: max_phys
      !!! Max. extent of the comput. domain
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: min_sub
      !!! Min. extent of the subdomains
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: max_sub
      !!! Max. extent of the subdomains
      INTEGER , DIMENSION(:,:), POINTER       :: subs_bc
      !!! Boundary conditions for subs on ALL processors
      INTEGER , DIMENSION(  :), INTENT(IN   ) :: sub2proc
      !!! Assignment of subs to procs
      INTEGER , DIMENSION(  :), INTENT(IN   ) :: isublist
      !!! List of subs handled by this processor
      INTEGER , DIMENSION(  :), INTENT(IN   ) :: bcdef
      !!! Boundary conditions on the computational box
      REAL(MK),                 INTENT(IN   ) :: ghostsize
      !!! Size of the ghostlayers
      INTEGER                 , INTENT(IN   ) :: nsubs
      !!! Total number of subs on all procs
      INTEGER                 , INTENT(IN   ) :: nsublist
      !!! Number of subs on current processor
      INTEGER                 , INTENT(INOUT) :: topoid
      !!! Topology ID. If 0, then create a new topology structure and return
      !!! ID. Else reallocate the topology with ID == topoid
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status. 0 upon success.
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: nneigh
      !!! Number of neighbors of each sub
      INTEGER , DIMENSION(:,:), INTENT(IN   ) :: ineigh
      !!! Neighbors of each sub.
      !!!
      !!! 1st index: 1...nneigh (neighbor index of sub)                        +
      !!! 2nd index: 1..nsubs (sub of which one wants to get the neighbors)

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER , DIMENSION(3)    :: ldc, ldl
      INTEGER                   :: i,j,k,kk,iopt,isize,iproc,isin
      INTEGER                   :: maxneigh,minbound,nsubmax,nsublistmax
      TYPE(ppm_t_topo), POINTER :: topo => NULL()
      REAL(MK)                  :: t0
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_topo_store',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Determine the maximum number of neighbours of all local subdomains
      !-------------------------------------------------------------------------
      maxneigh = 0
      DO i=1,nsublist
         j = nneigh(isublist(i))
         IF (j.GT.maxneigh) maxneigh = j
      ENDDO

      !-------------------------------------------------------------------------
      !  (Re)allocate memory for the topology structure
      !-------------------------------------------------------------------------
      CALL ppm_topo_alloc(topoid,nsubs,nsublist,maxneigh,MK,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_store', &
     &         'Could not allocate memory for Topology',__LINE__,info)
          GOTO 9999
      ENDIF
      topo => ppm_topo(topoid)%t

      !-------------------------------------------------------------------------
      !  store min_phys and max_phys
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      !-------------------------------------------------------------------------
      !  Single precision
      !-------------------------------------------------------------------------
      DO k=1,ppm_dim
         topo%min_physs(k) = min_phys(k)
         topo%max_physs(k) = max_phys(k)
      ENDDO
      topo%ghostsizes = ghostsize
#else
      !-------------------------------------------------------------------------
      !  Double precision
      !-------------------------------------------------------------------------
      DO k=1,ppm_dim
         topo%min_physd(k) = min_phys(k)
         topo%max_physd(k) = max_phys(k)
      ENDDO
      topo%ghostsized = ghostsize
#endif


      DO i=1,nsubs
          topo%sub2proc(i) = sub2proc(i)
      ENDDO

      !-------------------------------------------------------------------------
      !  initialize optmization status flags
      !-------------------------------------------------------------------------
      topo%isoptimized = .FALSE.

      !-------------------------------------------------------------------------
      !  Store the number of subdomains in the current topology
      !-------------------------------------------------------------------------
      topo%nsubs = nsubs

      !-------------------------------------------------------------------------
      !  Store the extent of the subdomains
      !-------------------------------------------------------------------------
      IF (ppm_dim.EQ.2) THEN
         !----------------------------------------------------------------------
         !  In two dimensions
         !----------------------------------------------------------------------
         DO i=1,nsubs
#if   __KIND == __SINGLE_PRECISION
            topo%min_subs(1,i) = min_sub(1,i)
            topo%min_subs(2,i) = min_sub(2,i)
            topo%max_subs(1,i) = max_sub(1,i)
            topo%max_subs(2,i) = max_sub(2,i)
#else
            topo%min_subd(1,i) = min_sub(1,i)
            topo%min_subd(2,i) = min_sub(2,i)
            topo%max_subd(1,i) = max_sub(1,i)
            topo%max_subd(2,i) = max_sub(2,i)
#endif
         ENDDO
      ELSE
         !----------------------------------------------------------------------
         !  In three dimensions
         !----------------------------------------------------------------------
         DO i=1,nsubs
#if   __KIND == __SINGLE_PRECISION
            topo%min_subs(1,i) = min_sub(1,i)
            topo%min_subs(2,i) = min_sub(2,i)
            topo%min_subs(3,i) = min_sub(3,i)
            topo%max_subs(1,i) = max_sub(1,i)
            topo%max_subs(2,i) = max_sub(2,i)
            topo%max_subs(3,i) = max_sub(3,i)
#else
            topo%min_subd(1,i) = min_sub(1,i)
            topo%min_subd(2,i) = min_sub(2,i)
            topo%min_subd(3,i) = min_sub(3,i)
            topo%max_subd(1,i) = max_sub(1,i)
            topo%max_subd(2,i) = max_sub(2,i)
            topo%max_subd(3,i) = max_sub(3,i)
#endif
         ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Store the external boundary conditions for this topology
      !-------------------------------------------------------------------------
      DO i=1,2*ppm_dim
         topo%bcdef(i) = bcdef(i)
      ENDDO

      !-------------------------------------------------------------------------
      !  The MAX of nsublist and 1 is needed to avoid allocation failures
      !  if a processor has 0 subs. Same for nsubs if there are no subs at
      !  all.
      !-------------------------------------------------------------------------
      nsubmax = MAX(nsubs,1)
      nsublistmax = MAX(nsublist,1)


      !-------------------------------------------------------------------------
      !  store the data by looping over the subs handled by the current
      !  processor: isublist(1:nsublist)
      !-------------------------------------------------------------------------
      IF (nsublist .LT. 1) THEN
          ! dummy entry if we have no subs at all
          topo%nneighsubs(1) = ppm_param_undefined
      ENDIF
      DO i=1,nsublist
         j = nneigh(isublist(i))
         topo%nneighsubs(i) = j
      ENDDO



      !-------------------------------------------------------------------------
      !  store neighboring subs by looping over the subs handled by the current
      !  processor: isublist(1:nsublist) - the neighbours are the accessed
      !  through the ineigh(1:nsublist,)
      !-------------------------------------------------------------------------
      DO i=1,nsublist
         DO j=1,topo%nneighsubs(i)
            topo%ineighsubs(j,i) = ineigh(j,isublist(i))
         ENDDO
      ENDDO


      !-------------------------------------------------------------------------
      !  And store it BC
      !-------------------------------------------------------------------------
      IF (nsubs .LT. 1) THEN
         ! dummy entries if we have no subs at all
         topo%subs_bc(1,1) = ppm_param_undefined
         topo%subs_bc(2,1) = ppm_param_undefined
         topo%subs_bc(3,1) = ppm_param_undefined
         topo%subs_bc(4,1) = ppm_param_undefined
         IF (ppm_dim.EQ.3) THEN
            topo%subs_bc(5,1) = ppm_param_undefined
            topo%subs_bc(6,1) = ppm_param_undefined
         ENDIF
      ENDIF
      DO i=1,nsubs
         topo%subs_bc(1,i) = subs_bc(1,i)
         topo%subs_bc(2,i) = subs_bc(2,i)
         topo%subs_bc(3,i) = subs_bc(3,i)
         topo%subs_bc(4,i) = subs_bc(4,i)
         IF (ppm_dim.EQ.3) THEN
            topo%subs_bc(5,i) = subs_bc(5,i)
            topo%subs_bc(6,i) = subs_bc(6,i)
         ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  And now deallocate the local array: subs_bc WHY?
      !-------------------------------------------------------------------------
!      iopt = ppm_param_dealloc
!      CALL ppm_alloc(subs_bc,ldc,iopt,info)
!      IF (info .NE. 0) THEN
!          info = ppm_error_fatal
!          CALL ppm_error(ppm_err_dealloc,'ppm_topo_store',     &
!     &        'BCs for subs on local processor failed',__LINE__,info)
!          GOTO 9999
!      ENDIF

      !-------------------------------------------------------------------------
      !  store number of subdomains handled by
      !  the current processor
      !-------------------------------------------------------------------------

      topo%nsublist = nsublist

      IF (nsublist .LT. 1) THEN
         topo%isublist(1) = ppm_param_undefined
      ENDIF
      DO i=1,nsublist
         topo%isublist(i) = isublist(i)
      ENDDO




      !-------------------------------------------------------------------------
      !  Determine and store the neighbors of this processor
      !-------------------------------------------------------------------------
      topo%nneighproc = 0

      isize = SIZE(topo%ineighproc,1)

      DO k=1,topo%nsublist
          i = topo%isublist(k)
          DO j=1,nneigh(i)
              iproc = topo%sub2proc(ineigh(j,i))
              ! if it is not myself
              IF (iproc .NE. ppm_rank) THEN
                  isin = 0
                  DO kk=1,topo%nneighproc
                      IF (topo%ineighproc(kk) .EQ. iproc) isin = 1
                  END DO
                  ! if not already in list
                  IF (isin .EQ. 0) THEN
                      topo%nneighproc = topo%nneighproc + 1
                      IF (topo%nneighproc .GT. isize) THEN
                          ! kindly ask for more memory
                          isize = isize + 1
                          ldc(1) = isize
                          iopt = ppm_param_alloc_grow_preserve
                          CALL ppm_alloc(topo%ineighproc,ldc,iopt,info)
                          IF (info .NE. 0) THEN
                              info = ppm_error_fatal
                              CALL ppm_error(ppm_err_alloc,'ppm_topo_store',&
     &                            'sub neighbor list PPM_INEIGHLIST',  &
     &                            __LINE__,info)
                              GOTO 9999
                          ENDIF
                      END IF
                      ! add iproc to the list of neighbors
                      topo%ineighproc(topo%nneighproc) = iproc
                  END IF
              END IF
          END DO
      END DO

      topo%isdefined = .TRUE.

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_topo_store',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF ((topoid .LT. 0) .OR. (topoid .GT. SIZE(ppm_topo))) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_store', &
     &             'topoid must be >= 0 and <= SIZE(ppm_topo)',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (nsubs .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_store',  &
     &            'nsubs must be >0',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (nsublist .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_store',  &
     &            'nsublist must be >0',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (nsubs .LT. nsublist) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_store',  &
     &            'total number of subs is smaller than local one',__LINE__,info)
              GOTO 8888
          ENDIF
          DO i=1,nsubs
              DO j=1,ppm_dim
                  IF(max_sub(j,i) .LE. min_sub(j,i)) THEN
                      info = ppm_error_error
                      CALL ppm_error(ppm_err_argument,'ppm_topo_store',  &
     &                    'min_sub must be < max_sub',__LINE__,info)
                      GOTO 8888
                  ENDIF
              ENDDO
          ENDDO
          DO i=1,nsubs
              IF(sub2proc(i) .LT. 0) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_topo_store',  &
     &                'found sub not assigned to any processor',__LINE__,info)
                  GOTO 8888
              ENDIF
          ENDDO
 8888     CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_topo_store_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_topo_store_d
#endif
