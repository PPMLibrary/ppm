     !-------------------------------------------------------------------------
     !  Subroutines :               ppm_inl_k_helpers
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
     SUBROUTINE is_kNeighbor_s(p_idx, p_neigh, xp, rc1, rc2, isNeigh, dist)
#elif __KIND == __DOUBLE_PRECISION
     SUBROUTINE is_kNeighbor_d(p_idx, p_neigh, xp, rc1, rc2, isNeigh, dist) 
#endif
      !!! Given indices of two particles, checks whether the euclidian distance
      !!! is smaller than the sum of minimum cutoff radius of these particles
      !!! and the skin, then RETURNs TRUE if so. Works for nD.
          IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_double
#endif
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
          INTEGER,  INTENT(IN)                 :: p_idx
          !!! Index of first particle
          INTEGER,  INTENT(IN)                 :: p_neigh
          !!! Index of second particle
          REAL(MK), INTENT(IN), DIMENSION(:,:) :: xp
          !!! Coordinate array of particles
          REAL(MK), INTENT(IN)                 :: rc1
          !!! Cutoff radii 1
          REAL(MK), INTENT(IN)                 :: rc2
          !!! Cutoff radii 2
          LOGICAL,  INTENT(OUT)                :: isNeigh
          !!! Return value. TRUE if particles are neighbors
          REAL(MK), INTENT(OUT)                :: dist
          !!! Return value. Distance between these 2 particles

      !---------------------------------------------------------------------
      !  Local variables and counters
      !---------------------------------------------------------------------
          INTEGER                              :: i
          ! Counter

          ! Return FALSE if they are the same particle
          IF(p_idx .EQ. p_neigh) THEN
              isNeigh = .FALSE.
              RETURN
          ENDIF

          ! Add squares of distances on each axis, then take square root of it
          dist = sqrt(sum((xp(1:ppm_dim, p_idx) - xp(1:ppm_dim, p_neigh))**2))

          ! Return TRUE if they are neighbors
          IF(dist .GT. rc1 .AND. dist .LE. rc2) THEN
              isNeigh = .TRUE.
          ELSE
              isNeigh = .FALSE.
          ENDIF

#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE is_kNeighbor_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE is_kNeighbor_d
#endif

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE is_xset_kNeighbor_s(red_idx, blue_idx, red, blue, &
 &                  rc1, rc2, isNeigh, dist)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE is_xset_kNeighbor_d(red_idx, blue_idx, red, blue, &
 &                  rc1, rc2, isNeigh, dist)
#endif
      !!! Given indices of two particles, returns whether the euclidian distance
      !!! between the 2 particles is between the 2 cutoffs rc1 and rc2
      !!!   rc1 < dist(red_idx,blue_idx) <= rc2
          IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_double
#endif
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
          INTEGER,  INTENT(IN)                 :: red_idx
          !!! Index of first particle
          INTEGER,  INTENT(IN)                 :: blue_idx
          !!! Index of second particle
          REAL(MK), INTENT(IN), DIMENSION(:,:) :: red
          !!! Coordinate array of particles (red)
          REAL(MK), INTENT(IN), DIMENSION(:,:) :: blue
          !!! Coordinate array of particles (blue)
          REAL(MK), INTENT(IN)                 :: rc1
          !!! Cutoff radii 1
          REAL(MK), INTENT(IN)                 :: rc2
          !!! Cutoff radii 2
          LOGICAL , INTENT(OUT)                :: isNeigh
          !!! Return value. TRUE if particles are neighbors
          REAL(MK), INTENT(OUT)                :: dist
          !!! Return value. Distance between particles

      !---------------------------------------------------------------------
      !  Local variables and counters
      !---------------------------------------------------------------------
          INTEGER                              :: i
          ! Counter

          dist = sqrt(sum((red(1:ppm_dim, red_idx) - blue(1:ppm_dim, blue_idx))**2))

          ! Return TRUE if they are neighbors
          IF(dist .GT. rc1 .AND. dist .LE. rc2) THEN
              isNeigh = .TRUE.
          ELSE
              isNeigh = .FALSE.
          ENDIF

#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE is_xset_kNeighbor_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE is_xset_kNeighbor_d
#endif

#if   __KIND == __SINGLE_PRECISION
      RECURSIVE SUBROUTINE getParticlesInCellDomain_s(cell_idx, xp, clist, list, nlist)
#elif   __KIND == __DOUBLE_PRECISION
      RECURSIVE SUBROUTINE getParticlesInCellDomain_d(cell_idx, xp, clist, list, nlist)
#endif
      !!! Given the cell index, this subroutine modifies the list array such
      !!! that it contains the part IDs of this cell domain and sets nlist to
      !!! number of particles in this cell domain.
      !!! NOTE: SLOW - use a k-d tree algorithm instead?
          IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
          INTEGER(ppm_kind_int64),  INTENT(IN)                   :: cell_idx
          REAL(MK),  DIMENSION(:,:),INTENT(IN)                   :: xp
          !!! this is basically a dummy argument to force fortran to generate
          !!! two versions of this routine
          TYPE(ppm_clist),          INTENT(IN)                   :: clist
          INTEGER,   DIMENSION(:),  INTENT(INOUT)                :: list
          INTEGER,                  INTENT(INOUT)                :: nlist

      !-------------------------------------------------------------------------
      !  Local variables and counters
      !-------------------------------------------------------------------------
          INTEGER(ppm_kind_int64)                      :: parentIdx
          INTEGER(ppm_kind_int64)                      :: childIdx
          INTEGER                                      :: left_end
          INTEGER                                      :: right_end
          INTEGER                                      :: border_idx
          INTEGER                                      :: i

          ! Get index of parent of this cell
          parentIdx  = parent(cell_idx)

          ! Get position on borders array that the parent is located in.
          border_idx = hash_search(clist%lookup,parentIdx)


          ! If this cell is not found in hash table, then put the cell index in
          ! empty list and return.
          IF(border_idx .EQ. htable_null)  THEN
              CALL putInEmptyList(cell_idx)
              RETURN
          END IF

          ! For 2D case
          IF(ppm_dim .EQ. 2)    THEN
              ! If the cell does not contain any particles that are in deeper
              ! levels in its region ...
              IF((clist%borders(5, border_idx) - clist%borders(1, border_idx)) &
                  .EQ. clist%borders(6, border_idx))  THEN
                  ! Put it in empty list
                  CALL putInEmptyList(cell_idx)
              ELSE
                  childIdx  = child(cell_idx)
                  DO i=0,2**ppm_dim-1
                      CALL getParticlesInCellDomain(childIdx+i, xp, clist,&
                          list, nlist)
                  ENDDO
              END IF

              ! Get index of first column on borders array.
              left_end  = 1 + cell_idx - (4*parentIdx-2)

              ! Get index of last column on borders array.
              right_end = left_end + 1

              ! If this is the top level, then get all particles
              IF(parentIdx .EQ. 0)  then
                   left_end  = 1
                   right_end = 5
              END IF
          ! For 3D case
          ELSEIF(ppm_dim .EQ. 3)   THEN
              ! If the cell does not contain any particles that are in deeper
              ! levels in its region ...
              IF((clist%borders(9, border_idx) - clist%borders(1, border_idx)) &
                  .EQ. clist%borders(10, border_idx))  THEN
                  ! Put it in empty list
                  CALL putInEmptyList(cell_idx)
              ELSE
                  childIdx  = child(cell_idx)
                  DO i=0,2**ppm_dim-1
                      CALL getParticlesInCellDomain(childIdx+i, xp, clist,&
                          list, nlist)
                  ENDDO
              END IF

              ! Get index of first column on borders array.
              left_end  = 1 + cell_idx - (8*parentIdx-6)

              ! Get index of last column on borders array.
              right_end = left_end + 1

              ! If this is the top level, then get all particles
              IF(parentIdx .EQ. 0)  then
                   left_end  = 1
                   right_end = 9
              END IF
          END IF

          ! From first column to last, get all particles and put them in the list
          DO i = (clist%borders(left_end, border_idx) + 1), &
 &                clist%borders(right_end, border_idx)
              nlist = nlist + 1
              list(nlist) = clist%rank(i)
          END DO


#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE getParticlesInCellDomain_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE getParticlesInCellDomain_d
#endif

#if   __KIND == __SINGLE_PRECISION
      RECURSIVE SUBROUTINE sortByDist_s(idx, dist, rank, vlist)
#elif __KIND == __DOUBLE_PRECISION
      RECURSIVE SUBROUTINE sortByDist_d(idx, dist, rank, vlist)
#endif
      !!! Recursive algorithm which sorts neighbours by their distances
      !!! to a center particle, in ascending order.
      !!! (swapping is actually done on rank array and vlist).
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      INTEGER,  INTENT(IN)                  :: idx
      !!! Id of the particle of interest
      REAL(MK), INTENT(IN),    DIMENSION(:) :: dist
      !!! Input array for neighbours distances
      INTEGER,  INTENT(INOUT), DIMENSION(:) :: rank
      !!! ranks of neigbours, from nearest to furtherst
      INTEGER,  INTENT(INOUT), DIMENSION(:,:) :: vlist
      !!! sorted neighbour list for the particle of interest
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER                               :: marker

      ! Keep on partitioning and sorting, recursively.
      IF(size(rank) .GT. 1) THEN
         CALL partitionByDist(idx, dist, rank, vlist, marker)
         CALL sortByDist(idx, dist, rank(:marker-1),vlist(:marker-1,:))
         CALL sortByDist(idx, dist, rank(marker:),vlist(marker:,:))
      END IF
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE sortByDist_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE sortByDist_d
#endif

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE partitionByDist_s(idx, dist, rank, vlist, marker)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE partitionByDist_d(idx, dist, rank, vlist, marker)
#endif
      !!! This subroutine partitions neighbour list by their distances and
      !!! updates marker
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      INTEGER,  INTENT(IN)                    :: idx
      !!! Id of the particle of interest
      REAL(MK), INTENT(IN),    DIMENSION(:)   :: dist
      !!! Input array for neighbours distances
      INTEGER,  INTENT(INOUT), DIMENSION(:)   :: rank
      !!! ranks of neigbours, from nearest to furtherst
      INTEGER,  INTENT(INOUT), DIMENSION(:,:) :: vlist
      !!! sorted neighbour list for the particle of interest
      INTEGER,  INTENT(OUT)                   :: marker
      !!! marker to be updated

      !---------------------------------------------------------------------
      !  Local variables and counters
      !---------------------------------------------------------------------
      INTEGER                               :: i
      INTEGER                               :: j
      INTEGER                               :: temp
      REAL(MK)                              :: pivot
      
      ! This pivoting strategy broke the code in some specific cases
      !pivot = (cutoff(rank(1)) + cutoff(rank(size(rank))))/2
      pivot = dist(rank(size(rank)/2))
      i= 0
      j= size(rank) + 1

      DO WHILE(i .LT. j)
         j = j - 1
         DO WHILE(dist(rank(j)) .GT. pivot)
            j = j-1
         END DO
         i = i + 1
         DO WHILE(dist(rank(i)) .LT. pivot)
            i = i + 1
         END DO
         IF (i .LT. j) THEN
            temp = rank(i)
            rank(i) = rank(j)
            rank(j) = temp
            temp = vlist(i,idx)
            vlist(i,idx) = vlist(j,idx)
            vlist(j,idx) = temp
         END IF
      ENDDO


      IF (i .EQ. j) THEN
         marker = i + 1
      ELSE
         marker = i
      END IF
      
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE partitionByDist_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE partitionByDist_d
#endif

