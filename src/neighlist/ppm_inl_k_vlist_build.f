     !-------------------------------------------------------------------------
     !  Module   :                  ppm_neighlist
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
#if   __ACTION == __COUNT
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE count_neigh_s(p_idx, clist, domain, xp, cutoff, knn, nvlist)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE count_neigh_d(p_idx, clist, domain, xp, cutoff, knn, nvlist)
#endif
#elif __ACTION == __GET
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE get_neigh_s(p_idx, clist, domain, xp, cutoff, knn, vlist, nvlist)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE get_neigh_d(p_idx, clist, domain, xp, cutoff, knn, vlist, nvlist)
#endif
#endif
      !!! Given the particle index, this subroutine locates the cell that this
      !!! cell is located in, gathers all particles in these cells and updates
      !!! verlet lists of all these particles.
          IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
          INTEGER,  INTENT(IN)                 :: p_idx
          !!! Particle index
          TYPE(ppm_clist)  , INTENT(IN)        :: clist
          !!! cell list
          REAL(MK), DIMENSION(2*ppm_dim)       :: domain
          !!! Pysical extent of whole domain including ghost layers
          REAL(MK), INTENT(IN), DIMENSION(:,:) :: xp
          !!! Particle coordinates
          REAL(MK), INTENT(IN), DIMENSION(:)   :: cutoff
          !!! Particle cutoff radii
          INTEGER,  INTENT(IN)                 :: knn
          !!! Minimum number of neighbours
#if __ACTION == __GET
          INTEGER,  DIMENSION(:,:)             :: vlist
          !!! Verlet list, where vlist(j, i) contains the jth neighbor of
          !!! ith particle
#endif
          INTEGER,  DIMENSION(:)               :: nvlist
          !!! Number of neighbors of particles. nvlist(i) contains number of
          !!! neighbors particle i has.
          INTEGER                              :: info

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
          INTEGER                              :: p_ref
          ! Reference particle
          INTEGER                              :: p_neigh
          ! Particle in the neighbor cell
          REAL(MK),    DIMENSION(ppm_dim)      :: offset_coor
          ! Offset from midpoint coordinates of reference array to get to
          ! neighbor cells.
          REAL(MK),    DIMENSION(ppm_dim)      :: c_coor
          ! Coordinates
          REAL(MK),    DIMENSION(ppm_dim)      :: n_coor
          ! Coordinates of neighbor cell
          INTEGER(ppm_kind_int64)              :: c_idx
          ! Cell index
          INTEGER(ppm_kind_int64)              :: n_idx
          ! Cell index
          INTEGER(ppm_kind_int64)              :: parentIdx
          ! Index of parent cell
          REAL(MK),    DIMENSION(ppm_dim)      :: p_coor
          ! Particle coordinates
          INTEGER                              :: p_depth
          ! Depth of the particle
          INTEGER                              :: c_depth
          ! Depth of cell
          REAL(MK), PARAMETER                  :: skin = 0._MK
          ! Skin parameter
          REAL(MK)                             :: cell_width, cell_width2
          ! Width of cell
          REAL(MK)                             :: dist
          ! distance between 2 particles
          LOGICAL                              :: isNeighbour
          REAL(MK), DIMENSION(:),POINTER        :: dist_array => NULL()

      !-------------------------------------------------------------------------
      !  Counters
      !-------------------------------------------------------------------------
          INTEGER                              :: i
          INTEGER                              :: j
          INTEGER                              :: m
          INTEGER                              :: n

          !---------------------------------------------------------------------
          !  If this particle was visited before, skip it.
          !---------------------------------------------------------------------
          IF(used(p_idx))   RETURN

#if   __KIND == __SINGLE_PRECISION
          dist_array => dist_array_s
#elif __KIND == __DOUBLE_PRECISION
          dist_array => dist_array_d
#endif
          !---------------------------------------------------------------------
          !  Initialize position on empty list to 0.
          !---------------------------------------------------------------------
          empty_pos = 0

          !---------------------------------------------------------------------
          !  Get coordinates and depth of the particle
          !---------------------------------------------------------------------
          CALL getParticleCoorDepth(p_idx, domain, p_coor, p_depth, xp, cutoff, skin)

          !---------------------------------------------------------------------
          !  Get index of the cell that this particle is located in.
          !---------------------------------------------------------------------
          c_idx   = getCellIdx(p_coor, p_depth, domain)

          !---------------------------------------------------------------------
          !  Get coordinates and depth of the particle
          !---------------------------------------------------------------------
          CALL getCellCoor_Depth(c_idx, domain, c_coor, c_depth, & 
 &                               clist%max_depth, info)

          !---------------------------------------------------------------------
          !  Compute offset coordinates, which will be used to find neighbor cells
          !---------------------------------------------------------------------
          offset_coor(1) = (domain(2) - domain(1))/(2**(c_depth-1))
          offset_coor(2) = (domain(4) - domain(3))/(2**(c_depth-1))
          DO i = 3,ppm_dim
              offset_coor(i) = (domain(2*i) - domain(2*i-1))/(2**(c_depth-1))
          END DO

          !---------------------------------------------------------------------
          !  From first cell until center cell, check whether they are empty
          !  or not and store those in empty list that contain no particles in
          !  their region.
          !---------------------------------------------------------------------
          ! From first cell to (center-1)th cell
          DO i = 1, ((3**ppm_dim - 1)/2)
              ! Compute midpoints of neighbor cell
              n_coor(1) = c_coor(1) + (ncells(i,1)*offset_coor(1))
              n_coor(2) = c_coor(2) + (ncells(i,2)*offset_coor(2))
              DO j = 3, ppm_dim
                  n_coor(j) = c_coor(j) + (ncells(i,j)*offset_coor(j))
              END DO

              ! If found coordinates are not inside domain, skip!
              IF(.NOT. inDomain(n_coor, domain))   cycle

              ! Get index of the neighbor cell
              n_idx = getCellIdx(n_coor, c_depth, domain)

              ! Store in empty list if empty
              IF(isEmpty(n_idx,clist%lookup))     CALL putInEmptyList(n_idx)
          END DO

          !---------------------------------------------------------------------
          !  Get particles that are in the same cell with input particle.
          !---------------------------------------------------------------------
          CALL getParticlesInCell(c_idx, xp, clist, own_plist, own_nplist)

          ! For each reference particle
          pref_loop: DO m = 1, own_nplist
              ! Pick a reference particle
              p_ref = own_plist(m)

              ! if its a ghost, cycle
              IF(p_ref .GT. clist%n_real_p) CYCLE

              ! Get coordinates and depth of this particle
              CALL getParticleCoorDepth(p_ref, domain, p_coor, p_depth, &
 &                                      xp, cutoff, skin)

              ! Starts at the very bottom + 1
              p_depth = clist%max_depth

              depth_loop: DO WHILE(p_depth .GE. 1)
                  ! Get the child cell this particle is located in
                  c_idx = getCellIdx(p_coor, p_depth, domain)

                  ! Get cell coordinates and depth
                  CALL getCellCoor_Depth(c_idx, domain, c_coor, c_depth, &
 &                                       clist%max_depth, info)
                  cell_width = getMinimumSideLength(domain)/(2**c_depth)
                  IF (c_depth.NE.clist%max_depth) THEN
                      cell_width2 = cell_width/2
                  ELSE
                      cell_width2 = 0._MK
                  ENDIF

                  ! Compute offset coordinates to compute midpoint coordinates
                  ! of neighbor cells at this depth.
                  offset_coor(1) = (domain(2) - domain(1))/(2**(c_depth-1))
                  offset_coor(2) = (domain(4) - domain(3))/(2**(c_depth-1))
                  DO i = 3, ppm_dim
                      offset_coor(i) = (domain(2*i) - domain(2*i-1))/&
 &                                     (2**(c_depth-1))
                  ENDDO

                  ! For each neighbor
                  DO i = 1, 3**ppm_dim
                      ! Compute midpoint coordinates of the neighbor
                      n_coor(1) = c_coor(1) + (ncells(i,1)*offset_coor(1))
                      n_coor(2) = c_coor(2) + (ncells(i,2)*offset_coor(2))
                      DO j = 3,ppm_dim
                          n_coor(j) = c_coor(j) + (ncells(i,j)*offset_coor(j))
                      ENDDO

                      ! If coordinates are not inside the domain, skip!
                      IF(.NOT. inDomain(n_coor, domain))   cycle

                      ! Get index of the neighbor cell
                      n_idx = getCellIdx(n_coor, c_depth, domain)

                      ! Get particles in the neighbor cell
                      !-----------------------------------------------------
                      ! REMARK: If this cell has no particles, it will be stored
                      !         in empty cell inside getParticlesInCell
                      !         subroutine. If it has particles only at this
                      !         depth, but no deeper particles, it will again
                      !         be saved in the empty list just  before returning.
                      !-----------------------------------------------------
                      neigh_nplist = 0
                      CALL getParticlesInCellDomain(n_idx, xp, clist, &
                          &         neigh_plist, neigh_nplist)
                      ! For each neighbor element
                      DO n = 1, neigh_nplist
                          ! Pick a candidate for neighbor particle
                          p_neigh = neigh_plist(n)
                          ! If they are neighbors 
                          CALL is_kNeighbor(p_ref,p_neigh,xp, &
                              cell_width2,cell_width,isNeighbour,dist)
                          IF (isNeighbour)  THEN
                          ! Store neighbor particle in verlet list
                          ! of reference particle
                              nvlist(p_ref) = nvlist(p_ref) + 1
#if __ACTION == __GET
                              vlist(nvlist(p_ref), p_ref) = p_neigh
                              dist_array(nvlist(p_ref)) = dist
                              dist_rank(nvlist(p_ref)) = nvlist(p_ref)
#endif
                          ENDIF
                      ENDDO
                  ENDDO
                  ! IF we have enough neighbour, exit the loop
                  IF (nvlist(p_ref).GE.knn) EXIT depth_loop
                  ! else, search further away by decrementing the depth by 1.
                  p_depth = p_depth - 1
              ENDDO depth_loop
              IF (nvlist(p_ref).LT.knn) THEN
                  WRITE(*,*) 'Could not find enough neighbours' 
                  WRITE(*,'(3(A,I0))') 'nvlist(',p_ref,') = ',nvlist(p_ref),'<',knn
                  STOP
              ENDIF
#if __ACTION == __GET
              CALL sortByDist(p_ref, dist_array(1:nvlist(p_ref)), &
                  dist_rank(1:nvlist(p_ref)), vlist)
              nvlist(p_ref) = MIN(knn,nvlist(p_ref))
#endif
          ENDDO pref_loop

          ! Label all reference particles as USED so that they wont be visited again.
          DO i = 1, own_nplist
              used(own_plist(i)) = .TRUE.
          END DO
#if   __ACTION == __COUNT
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE count_neigh_s
#elif   __KIND == __DOUBLE_PRECISION
      END SUBROUTINE count_neigh_d
#endif
#elif __ACTION == __GET
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE get_neigh_s
#elif   __KIND == __DOUBLE_PRECISION
      END SUBROUTINE get_neigh_d
#endif
#endif

