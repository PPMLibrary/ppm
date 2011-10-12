     !-------------------------------------------------------------------------
     !  Module   :                  ppm_inl_clist
     !-------------------------------------------------------------------------
     ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Lab (ETH Zurich),
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
      SUBROUTINE create_inl_clist_s(xp, Np, Mp, cutoff, skin, actual_domain, &
     & ghost_extend, lsymm, clist, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE create_inl_clist_d(xp, Np, Mp, cutoff, skin, actual_domain, &
     & ghost_extend, lsymm, clist, info)
#endif
      !!! Given particle coordinates(xp), number of all particles including
      !!! ghost particles(Mp), cutoff radii of particles(cutoff), skin parameter
      !!! and the physical extent of whole domain including ghost layers;
      !!! this subroutine creates cell lists that are to be stored in the hash
      !!! table where borders on the particle arrays are to be stored.
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      REAL(MK), INTENT(IN), DIMENSION(:,:)         :: xp
      !!! Particle coordinates array. F.e., xp(1, i) is the x-coor of particle i.
      INTEGER , INTENT(IN)                         :: Np
      !!! Number of real particles
      INTEGER , INTENT(IN)                         :: Mp
      !!! Number of all particles including ghost particles
      REAL(MK), INTENT(IN), DIMENSION(:)           :: cutoff
      !!! Particles cutoff radii
      REAL(MK), INTENT(IN)                         :: skin
      !!! Skin parameter
      REAL(MK),      DIMENSION(2*ppm_dim)          :: actual_domain
      ! Physical extent of actual domain without ghost layers.
      REAL(MK), INTENT(IN), DIMENSION(ppm_dim)     :: ghost_extend
      !!! Extra area/volume over the actual domain introduced by
      !!! ghost layers.
      LOGICAL,  INTENT(IN)                         :: lsymm
      !!! If lsymm = TRUE, verlet lists are symmetric and we have ghost
      !!! layers only in (+) directions in all axes. Else, we have ghost
      !!! layers in all directions.
      INTEGER , INTENT(OUT)                        :: info
      !!! Info to be returned. 0 if SUCCESSFUL.
      TYPE(ppm_clist), INTENT(INOUT)               :: clist

      !---------------------------------------------------------------------
      !  Local variables and parameters
      !---------------------------------------------------------------------
      REAL(MK),                DIMENSION(2*ppm_dim) :: whole_domain
      ! Physical extent of whole domain including ghost layers.
      INTEGER(ppm_kind_int64), PARAMETER            :: idx = 1
      ! Parameter used for recursive functions.
      INTEGER(ppm_kind_int64), PARAMETER            :: idx0 = 0
      ! Parameter used for recursive function for first level only.
      INTEGER                                       :: level
      ! Depth level.

      REAL(MK)                                      :: max_size
      REAL(MK)                                      :: size_diff

      !---------------------------------------------------------------------
      !  Parameters for ppm_alloc
      !---------------------------------------------------------------------
      INTEGER                                       :: iopt
      INTEGER, DIMENSION(2)                         :: lda
      INTEGER, DIMENSION(1)                         :: ldl
      INTEGER, DIMENSION(1)                         :: ldu

      !---------------------------------------------------------------------
      !  Counters
      !---------------------------------------------------------------------
      INTEGER                                       :: i
      REAL(MK)                                      :: t0

      !<<<<<<<<<<<<<<<<<<<<<<<<< Start of the code >>>>>>>>>>>>>>>>>>>>>>>>>!

      CALL substart('ppm_create_inl_clist',t0,info)

      max_size = 0.0_mk
      IF(lsymm) THEN
          DO i = 1, ppm_dim
              whole_domain(2*i-1) = actual_domain(2*i-1)
              whole_domain(2*i)   = actual_domain(2*i) + ghost_extend(i)
              max_size = MAX(max_size, (whole_domain(2*i) - whole_domain(2*i-1)))
          END DO
      ELSE
          DO i = 1, ppm_dim
              whole_domain(2*i-1) = actual_domain(2*i-1) - ghost_extend(i)
              whole_domain(2*i)   = actual_domain(2*i)   + ghost_extend(i)
              max_size = MAX(max_size, (whole_domain(2*i) - whole_domain(2*i-1)))
          END DO
      END IF 
      DO i = 1, ppm_dim
          IF ((whole_domain(2*i) - whole_domain(2*i-1)) .LE. max_size/2.0_MK) THEN
              whole_domain(2*i)   = whole_domain(2*i-1) + max_size
          END IF
      END DO

      clist%n_real_p = Np
      clist%n_all_p = Mp

      !-------------------------------------------------------------------------
      !  Allocate rank array, which will be used to sort particles, instead
      !  of modifying input arrays which are xp and cutoff.
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      lda(1) = clist%n_all_p
      CALL ppm_alloc(clist%rank, lda, iopt, info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_create_inl_clist',     &
 &                       'rank',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Create hash table and allocate borders array, then sort particles
      !  by their position. If hash table is not large enough, destroy it and
      !  deallocate borders array. Then double the size and retry, until
      !  it is successful.
      !-------------------------------------------------------------------------



      clist%grow_htable = .TRUE.
      clist%ncell = CEILING(clist%n_all_p/1.0) !Hardcoded estimation of number of cells
      DO WHILE(clist%grow_htable)
          clist%grow_htable = .FALSE.     ! set insufficient_hash_table flag to 0
          CALL create_htable(clist%lookup,clist%ncell,info) ! create hash table
          lda(2) = clist%lookup%nrow
          iopt = ppm_param_alloc_fit
          ! Number of rows of "borders" array depends on dimensionality.
          IF(ppm_dim .EQ. 2)       THEN
              lda(1) = 6
              ! Allocate "borders" array for 2D case
              CALL ppm_alloc(clist%borders, lda, iopt, info)
              IF (info.NE.0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'ppm_create_inl_clist',     &
 &                               'borders',__LINE__,info)
              END IF
          ELSEIF(ppm_dim .EQ. 3)   THEN
              lda(1) = 10
              ! Allocate "borders" array for 3D case
              CALL ppm_alloc(clist%borders, lda, iopt, info)
              IF (info.NE.0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'ppm_create_inl_clist',     &
 &                               'borders',__LINE__,info)
              END IF
          END IF

          ! Initialize array and variables
          clist%borders         = 0
          clist%borders_pos     = 0
          clist%borders_pos_max = 0

          clist%borders_pos = clist%borders_pos + 1
          CALL hash_insert(clist%lookup,idx0, clist%borders_pos, info)
          IF(ppm_dim .EQ. 2)    THEN
              clist%borders(6, clist%borders_pos)  = clist%n_all_p
          ELSE
              clist%borders(10, clist%borders_pos) = clist%n_all_p
          ENDIF

          DO i = 1, clist%n_all_p
              clist%rank(i) = i
          END DO

          ! Sort particles by their position
          CALL SortByPosition(xp, cutoff, skin, clist%rank, clist,whole_domain, idx, 0)
          IF(clist%grow_htable) THEN ! If hash table is not sufficient
              CALL destroy_htable(clist%lookup,info)        ! Destroy hash table
              iopt = ppm_param_dealloc
              lda = 0
              ! Deallocate "borders" array
              CALL ppm_alloc(clist%borders, lda, iopt, info)
              IF (info.NE.0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_dealloc,'ppm_create_cell_list',   &
 &                               'borders',__LINE__,info)
              END IF
              clist%ncell = clist%ncell*2             ! Double the number of cells
          END IF
      END DO

      !-------------------------------------------------------------------------
      !  Re-initialize rank array.
      !-------------------------------------------------------------------------
      DO i = 1, clist%n_all_p
          clist%rank(i) = i
      ENDDO

      !-------------------------------------------------------------------------
      !  Sort particles by their cutoff radii in descending order.
      !-------------------------------------------------------------------------
      CALL sortByRC(cutoff, skin, clist%rank)

      !-------------------------------------------------------------------------
      !  Get maximum depth in the cell list.
      !-------------------------------------------------------------------------
      clist%max_depth = MAX(getMaxDepth(cutoff, clist, whole_domain),1)

      !-------------------------------------------------------------------------
      !  Allocate rc_borders array, in order to store borders on rank
      !  array, that will contain particles of that level.
      !-------------------------------------------------------------------------
      ldl(1) = 0
      ldu(1) = clist%max_depth
      CALL ppm_alloc(clist%rc_borders, ldl, ldu, iopt, info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_create_inl_clist',     &
 &                       'rc_borders',__LINE__,info)
          GOTO 9999
      END IF
      clist%rc_borders(0) = 1

      !-------------------------------------------------------------------------
      !  Get borders on particle rank array for different levels, such that
      !  consecutive indices in rc_borders array will be start and end
      !  indices on the rank array, which will contain particles of that
      !  level.
      !-------------------------------------------------------------------------
      CALL getRC_Borders(cutoff, skin, clist, whole_domain,info)

      !-------------------------------------------------------------------------
      !  For every level of depth, sort particles by their cutoff radii and
      !  positions. So in the end, we have particles that are sorted by
      !  their position in their own chunk of depth.
      !-------------------------------------------------------------------------
      CALL SortByRC_Pos(xp, cutoff, skin, clist%rank(clist%rc_borders(0):&
 &                      (clist%rc_borders(1) - 1)), clist, whole_domain, &
 &                      idx0, 1, clist%rc_borders(0)-1)

      DO level = 2, clist%max_depth
          CALL SortByRC_Pos(xp, cutoff, skin, clist%rank(clist%rc_borders(level-1):&
 &                          (clist%rc_borders(level) - 1)), clist, whole_domain, &
 &                          idx, level, clist%rc_borders(level-1)-1)
      END DO


      !-------------------------------------------------------------------------
      !  Deallocate rc_borders array as its of no use anymore.
      !-------------------------------------------------------------------------
      CALL ppm_alloc(clist%rc_borders, lda, iopt, info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_dealloc,'ppm_create_inl_clist',     &
  &                      'rc_borders',__LINE__,info)
          GOTO 9999
      END IF

      !-------------------------------------------------------------------------
      !  Set last rows of borders arrays to -1 if that cell contains
      !  particles in deeper levels and otherwise set to 1.
      !-------------------------------------------------------------------------
      IF(ppm_dim .EQ. 2)  THEN
          DO i = 1, clist%borders_pos_max
              IF((clist%borders(5, i) - clist%borders(1, i)) .NE. clist%borders(6, i))  THEN
                  clist%borders(6, i) = -1
              ELSE
                  clist%borders(6, i) = 1
              END IF
          END DO
      ELSEIF(ppm_dim .EQ. 3)  THEN
          DO i = 1, clist%borders_pos_max
              IF((clist%borders(9, i) - clist%borders(1, i)) .NE. clist%borders(10, i))  THEN
                  clist%borders(10, i) = -1
              ELSE
                  clist%borders(10, i) = 1
              END IF
          END DO
      END IF

9999  CONTINUE
      CALL substop('ppm_create_inl_clist',t0,info)
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE create_inl_clist_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE create_inl_clist_d
#endif

#if __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_destroy_inl_clist(clist,info)
      !!! deallocates all arrays in clist, sets variables back to
      !!! default values and calls the destructor for the hash table
      
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      TYPE(ppm_clist), INTENT(INOUT)               :: clist
      INTEGER , INTENT(OUT)                        :: info
      !!! Info to be returned. 0 if SUCCESSFUL.
      !---------------------------------------------------------------------
      !  Local Variables
      !---------------------------------------------------------------------
      INTEGER                                       :: iopt
      INTEGER, DIMENSION(2)                         :: lda
      REAL(ppm_kind_single)                         :: t0

      CALL substart('ppm_destroy_inl_clist',t0,info)
          
      iopt = ppm_param_dealloc
      lda = 0
      !-------------------------------------------------------------------------
      !  Deallocate borders array which contains cell lists.
      !-------------------------------------------------------------------------
      CALL ppm_alloc(clist%borders, lda, iopt, info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_dealloc,'destroy_inl_vlist',   &
          &                       'borders',__LINE__,info)
          GOTO 9999
      END IF

      !-------------------------------------------------------------------------
      !  Deallocate rank array.
      !-------------------------------------------------------------------------
      CALL ppm_alloc(clist%rank,lda, iopt, info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_dealloc,'destry_inl_vlist',   &
          &                       'rank',__LINE__,info)
          GOTO 9999
      END IF

      !-------------------------------------------------------------------------
      !  Deallocate rankByPos array.
      !-------------------------------------------------------------------------
      CALL ppm_alloc(clist%rankByPos,  lda, iopt, info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_dealloc,'destroy_inl_vlist',   &
          &                       'rankByPos',__LINE__,info)
          GOTO 9999
      END IF

      !-------------------------------------------------------------------------
      !  Deallocate rc_borders
      !-------------------------------------------------------------------------
      CALL ppm_alloc(clist%rc_borders, lda, iopt, info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_dealloc,'destroy_inl_vlist',   &
          &                       'rc_borders',__LINE__,info)
          GOTO 9999
      END IF

      !-------------------------------------------------------------------------
      !  Destroy hash table.
      !-------------------------------------------------------------------------
      CALL destroy_htable(clist%lookup,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_dealloc,'create_inl_vlist',   &
          &            'Could not destroy htable',__LINE__,info)
          GOTO 9999
      END IF

      NULLIFY(clist%rc_borders,clist%rankByPos,clist%rank,clist%borders)

      clist%borders_pos      = 0
      clist%borders_pos_max  = 0
      clist%max_depth        = 0
      clist%ncell            = 0
      clist%n_real_p         = 0
      clist%n_all_p          = 0
      clist%grow_htable      = .TRUE.
      
9999  CONTINUE
      CALL substop('ppm_destroy_inl_clist',t0,info)

      END SUBROUTINE ppm_destroy_inl_clist
#endif

#if __KIND == __SINGLE_PRECISION
      PURE FUNCTION parent(idx) RESULT(parent_idx)
      !!! Given the index of child cell, returns the index of its parent.
      !!! Works for nD.
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64), INTENT(IN) :: idx
      !!! Input index
      INTEGER(ppm_kind_int64)             :: parent_idx
      !!! Index of parent cell to be returned

      parent_idx = ISHFT(idx + (2**ppm_dim-2),-ppm_dim) !/(2**ppm_dim)
      END FUNCTION parent
#endif

#if __KIND == __SINGLE_PRECISION
      PURE FUNCTION child(idx) RESULT(child_idx)
      !!! Given the index of child cell, returns the index of its first child.
      !!! Works for nD.
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64), INTENT(IN) :: idx
      !!! Input index
      INTEGER(ppm_kind_int64)             :: child_idx
      !!! Index of first child cell to be returned

      child_idx = ISHFT(idx,ppm_dim) - (2**ppm_dim-2)
      END FUNCTION child
#endif

#if __KIND == __SINGLE_PRECISION
#ifdef __DEBUG
      FUNCTION isEmpty(c_idx,lookup) RESULT(empty)
#else
      PURE FUNCTION isEmpty(c_idx,lookup) RESULT(empty)
#endif
      !!! Given the index of the cell, returns whether the cell is
      !!! empty or not.
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64), INTENT(IN) :: c_idx
      !!! Input index
      TYPE(ppm_htable), INTENT(IN)        :: lookup
      !!! hash table
      logical                             :: empty
      !!! Logical result
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64)             :: parentIdx
      ! Index of parent of input cell
      INTEGER                             :: borders_pos
      !position on borders array

      empty = .TRUE.                           ! Set empty to TRUE
      parentIdx   = parent(c_idx)              ! Get index of parent cell
      borders_pos = hash_search(lookup,parentIdx) ! Search parent in hash table
      IF(borders_pos .EQ. htable_null)  RETURN ! Return FALSE if not found
      empty = .FALSE.                          ! Set empty to FALSE and return
      END FUNCTION isEmpty
#endif

#if __KIND == __SINGLE_PRECISION
      PURE FUNCTION getCellDepth(cell_idx)  RESULT(cell_depth)
      !!! Given the cell index, return the depth of the cell.
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64), INTENT(IN) :: cell_idx
      !!! Index of cell
      INTEGER                             :: cell_depth
      !!! Depth of cell to be returned

      !---------------------------------------------------------------------
      !  Local variables and counters
      !---------------------------------------------------------------------
      INTEGER                             :: depth      ! Local variable
      INTEGER(ppm_kind_int64)             :: idx        ! Local variable

      depth = 1
      idx = cell_idx
      ! Ferit: I believe this is faster than two log operations.
      !        As intrinsic log ops are defined for real parameters.

      ! Keep on incrementing depth, until you reach the father.
      DO WHILE(idx .GT. 1)
          idx   = parent(idx)
          depth = depth + 1
      END DO

      cell_depth = depth  ! Set result to computed depth
      END FUNCTION getCellDepth
#endif

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE getCellCoor_Depth_s(cell_idx, domain, coor, cell_depth, &
 &                     max_depth,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE getCellCoor_Depth_d(cell_idx, domain, coor, cell_depth, &
 &                     max_depth,info)
#endif
      !!! Given the cell index and the domain, modifies coor and
      !!! cell_depth variables such that coor contains midpoint
      !!! coordinates of the cell and cell_depth has the depth of
      !!! the cell. Works for nD.
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64), INTENT(IN)             :: cell_idx
      !!! Index of the cell whose midpoint coordinates and depth will be
      !!! returned.
      REAL(MK), DIMENSION(2*ppm_dim), INTENT(IN)      :: domain
      !!! Physical extent of whole domain including ghost layers.
      REAL(MK), DIMENSION(ppm_dim), INTENT(INOUT)     :: coor
      !!! Midpoint coordinates of the cell to be returned.
      INTEGER,                      INTENT(INOUT)     :: cell_depth
      !!! Depth of the cell to be returned.
      INTEGER,                      INTENT(IN)        :: max_depth
      !!! maximum cell depth in cell tree
      INTEGER                                         :: info
      !!! 0 on success

      !---------------------------------------------------------------------
      !  Local variables and counters
      !---------------------------------------------------------------------
      INTEGER                                         :: i
      ! Counter
      INTEGER                                         :: j
      ! Counter
      INTEGER(ppm_kind_int64)                         :: idx
      ! Cell index
      INTEGER(ppm_kind_int64)                         :: parent_idx
      ! Index of parent
      INTEGER(ppm_kind_int64)                         :: child_idx
      ! Index of child
      INTEGER                                         :: depth
      ! Depth of cell
      INTEGER(ppm_kind_int64), DIMENSION(1:max_depth) :: levelList
      ! Array of ancestors
      REAL(MK),      DIMENSION(ppm_dim)               :: mid_coor
      ! Midpoint coordinates
      REAL(MK),      DIMENSION(ppm_dim)               :: start_coor
      ! Min. extent of cell
      REAL(MK),      DIMENSION(ppm_dim)               :: end_coor
      ! Max. extent of cell
      INTEGER                                         :: childNum
      ! Rank of the child
      REAL(MK)                                        :: t0

      IF (ppm_debug .GE. 3) THEN
          CALL substart('GetCellCoor_Depth',t0,info)
      ENDIF


      ! Set minimum and maximum physical extent of first cell,
      ! which is the domain itself. Later, as we go deeper, we restrict
      ! these coordinates to get midpoint coordinates of the input cell.
      
      start_coor(1) = domain(1)
      end_coor(1)   = domain(2)
      start_coor(2) = domain(3)
      end_coor(2)   = domain(4)
      DO i = 3, ppm_dim
          start_coor(i) = domain(2*i - 1)
          end_coor(i)   = domain(2*i)
      END DO

      ! Get depth of the cell
      cell_depth = getCellDepth(cell_idx)
      ! Set local variable idx to cell index
      idx   = cell_idx
      ! Set local variable depth to cell depth.
      depth = cell_depth

      ! Fill levelList array with indices of ancestors. So, if the
      ! cell index is 16 which is a child of 4, which is a child of 1,
      ! levelList(1:3) will contain (1, 4, 16) respectively.
      DO WHILE(idx .GT. 0)
          levelList(depth) = idx
          idx   = parent(idx)
          depth = depth - 1
      END DO

      ! Start from the oldest ancestor down to the youngest child,
      ! restricting the domain accordingly, in order to get midpoint
      ! coordinates in the end.
      DO i = 1, cell_depth - 1

          mid_coor(1) = (start_coor(1) + end_coor(1))/2
          mid_coor(2) = (start_coor(2) + end_coor(2))/2
          DO j = 3, ppm_dim
              mid_coor(j) = (start_coor(j) + end_coor(j))/2
          END DO

          parent_idx = levelList(i)
          child_idx  = levelList(i+1)
          childnum = child_idx - ((2**ppm_dim)*(parent_idx-1) + 2)

          DO j = 1, ppm_dim
              IF(MOD(childnum/(2**(ppm_dim-j)),2) .EQ. 0) THEN
                  end_coor(j)   = mid_coor(j)
              ELSE
                  start_coor(j) = mid_coor(j)
              END IF
          END DO

          coor(1) = (start_coor(1) + end_coor(1))/2
          coor(2) = (start_coor(2) + end_coor(2))/2
          DO j = 3, ppm_dim
              coor(j) = (start_coor(j) + end_coor(j))/2
          END DO
      END DO

      IF (ppm_debug .GE. 3) THEN
          CALL substop('GetCellCoor_Depth',t0,info)
      ENDIF
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE getCellCoor_Depth_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE getCellCoor_Depth_d
#endif

#if   __KIND == __SINGLE_PRECISION
      PURE FUNCTION getCellIdx_s(coor, cell_depth, domain) RESULT(cell_idx)
#elif __KIND == __DOUBLE_PRECISION
      PURE FUNCTION getCellIdx_d(coor, cell_depth, domain) RESULT(cell_idx)
#endif
      !!! Given coordinates and depth, this subroutine returns the index
      !!! of the cell.Input coordinates do not have to be midpoint
      !!! coordinates, but any coordinates that are within the physical
      !!! region of the cell or on left and/or bottom boundaries of the cell.
      !!! Works for nD.
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      REAL(MK),     INTENT(IN), DIMENSION(ppm_dim)   :: coor
      !!! Input coordinates of the cell whose index is asked for.
      INTEGER,      INTENT(IN)                       :: cell_depth
      !!! Depth of the cell whose index is asked for.
      REAL(MK),     INTENT(IN), DIMENSION(2*ppm_dim) :: domain
      !!! Physical extent of whole domain including ghost layers.

      !---------------------------------------------------------------------
      !  Counters and local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64)                        :: cell_idx
      REAL(MK),    DIMENSION(ppm_dim)                :: mid_coor
      REAL(MK),    DIMENSION(ppm_dim)                :: start_coor
      REAL(MK),    DIMENSION(ppm_dim)                :: end_coor
      INTEGER                                        :: increment
      INTEGER                                        :: i
      INTEGER                                        :: j
      INTEGER,     DIMENSION(0:ppm_dim)              :: powdim

      ! Initialize physical extent of the cell whose index is asked
      ! for. Later in the loop, this physical extent will be
      ! restricted such that it will contain input coordinates,
      ! until the desired input depth is reached, so that we have
      ! the cell index in the end.
      DO i = 1, ppm_dim
          start_coor(i) = domain(2*i - 1)
          end_coor(i)   = domain(2*i)
      END DO
      ! precompute this and save the computation later
      DO i = 0, ppm_dim
          powdim(i)=2**i
      ENDDO

      cell_idx = 1     ! Set cell_idx to greatest ancestor.
      DO i = 1, cell_depth - 1
          ! Compute current midpoint coordinates
          mid_coor(1) = (start_coor(1) + end_coor(1))/2
          mid_coor(2) = (start_coor(2) + end_coor(2))/2
          DO j = 3, ppm_dim
              mid_coor(j) = (start_coor(j) + end_coor(j))/2
          END DO

          ! Set cell_idx to first child (bottom-left)
          cell_idx  = ((powdim(ppm_dim))*(cell_idx-1) + 2)

          ! Set rank of the child to 0
          increment = 0

          ! Compute rank of the child
          increment = increment + (powdim(ppm_dim-1))*INT&
 &                     ((coor(1)-start_coor(1))/(mid_coor(1)-start_coor(1)))
          increment = increment + (powdim(ppm_dim-2))*INT&
 &                     ((coor(2)-start_coor(2))/(mid_coor(2)-start_coor(2)))
          DO j = 3, ppm_dim
              increment = increment + (powdim(ppm_dim-j))*INT&
 &                         ((coor(j)-start_coor(j))/(mid_coor(j)-start_coor(j)))
          END DO

          ! Add the rank of the child on the first child to get
          ! correct index of the child cell.
          cell_idx = cell_idx + increment
          DO j = 1, ppm_dim
              IF(MOD(increment/(powdim(ppm_dim-j)),2) .EQ. 0) THEN
                  end_coor(j)   = mid_coor(j)
              ELSE
                  start_coor(j) = mid_coor(j)
              END IF
          END DO
      END DO
#if   __KIND == __SINGLE_PRECISION
      END FUNCTION getCellIdx_s
#elif __KIND == __DOUBLE_PRECISION
      END FUNCTION getCellIdx_d
#endif

#if   __KIND == __SINGLE_PRECISION
      PURE FUNCTION getMinimumRC_s(cutoff, skin, rank) RESULT(minRC)
#elif __KIND == __DOUBLE_PRECISION
      PURE FUNCTION getMinimumRC_d(cutoff, skin, rank) RESULT(minRC)
#endif
      !!! Given the cutoff radii and rank arrays and the skin
      !!! parameter, returns the minimum cutoff radius within these arrays.
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      REAL(MK),  INTENT(IN), DIMENSION(:) :: cutoff
      !!! Cutoff radii array
      REAL(MK),  INTENT(IN)               :: skin
      !!! Skin parameter
      INTEGER,   INTENT(IN), DIMENSION(:) :: rank
      !!! Rank array, containing particles ranks
      REAL(MK)                            :: minRC
      !!! Minimum cutoff radius to be returned

      !---------------------------------------------------------------------
      !  Counters
      !---------------------------------------------------------------------
      INTEGER                             :: i ! Counter

      ! Set minimum cutoff radius to a great value.
      minRC = HUGE(1)
      ! Search for a smaller cutoff radius through whole array of particles.
      IF(size(rank).GT.0) THEN
        minRC = cutoff(rank(1))
        DO i = 2, size(rank)
            IF(cutoff(rank(i)) .LT. minRC) THEN
                minRC = cutoff(rank(i))
            END IF
        END DO
      END IF

      ! Add skin parameter on the found minimum cutoff radius.
      minRC = minRC + skin
#if   __KIND == __SINGLE_PRECISION
      END FUNCTION getMinimumRC_s
#elif __KIND == __DOUBLE_PRECISION
      END FUNCTION getMinimumRC_d
#endif

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE setSubregions_s(ownregion, subregions,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE setSubregions_d(ownregion, subregions,info)
#endif
      !!! Given the input region, subdivides it and distributes the portions
      !!! to subregions. Works for nD.
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      REAL(MK), INTENT(IN),    DIMENSION(:)            :: ownregion
      !!! Coordinates of input region containing minimum and maximum
      !!! physical extents of it.
      REAL(MK), INTENT(INOUT), DIMENSION(:,:) :: subregions
      !!! Array of subregions, which will contain minimum and maximum
      !!! physical extents of each subregion. For example, subregions(5,3)
      !!! contains the x-coordinate of maximum physical extent of 5th
      !!! subregion.
      INTEGER                                                   :: info

      !---------------------------------------------------------------------
      !  Local variables and counters
      !---------------------------------------------------------------------
      REAL(MK), DIMENSION(ppm_dim)    :: mid_coor
      REAL(MK), DIMENSION(ppm_dim)    :: start_coor
      REAL(MK), DIMENSION(ppm_dim)    :: end_coor
      INTEGER                         :: i
      INTEGER                         :: j
      REAL(MK)                        :: t0

      IF(ppm_debug.GE.3)THEN
          CALL substart('setSubregions',t0,info)
      ENDIF
      
      ! Initialize minimum and maximum physical extents to be
      ! assigned to subregions.
      DO i = 1, ppm_dim
          start_coor(i) = ownregion(2*i - 1)
          end_coor(i)   = ownregion(2*i)
          mid_coor(i)   = (start_coor(i) + end_coor(i))/2
      END DO

      ! Distribute shares accordingly.
      DO i = 1, 2**ppm_dim ! over all subregions
          DO j = 1, ppm_dim ! over all dimensions
              IF(MOD((i-1)/(2**(ppm_dim - j)) ,2) .EQ. 0) THEN
                  subregions(i, 2*j - 1) = start_coor(j)
                  subregions(i, 2*j)     = mid_coor(j)
              ELSE
                  subregions(i, 2*j - 1) = mid_coor(j)
                  subregions(i, 2*j)     = end_coor(j)
              END IF
          END DO
      END DO

      IF(ppm_debug.GE.3)THEN
          CALL substop('setSubregions',t0,info)
      ENDIF
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE setSubregions_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE setSubregions_d
#endif

#if   __KIND == __SINGLE_PRECISION
      PURE FUNCTION getMinimumSideLength_s(domain) RESULT(minLength)
#elif __KIND == __DOUBLE_PRECISION
      PURE FUNCTION getMinimumSideLength_d(domain) RESULT(minLength)
#endif
      !!! Given a physical domain, returns the longest side in it.
      !!! Works for nD.
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      REAL(MK), INTENT(IN), DIMENSION(:)         :: domain
      !!! Physical extent of the domain
      REAL(MK)                                   :: minLength
      !!! Maximum side length

      !---------------------------------------------------------------------
      !  Counters
      !---------------------------------------------------------------------
      INTEGER                                    :: i ! Counter

      minLength = domain(2) - domain(1)
      DO i = 2, ppm_dim
          minLength = min(minLength, (domain(2*i) - domain(2*i-1)))
      END DO
#if   __KIND == __SINGLE_PRECISION
      END FUNCTION getMinimumSideLength_s
#elif __KIND == __DOUBLE_PRECISION
      END FUNCTION getMinimumSideLength_d
#endif

#if   __KIND == __SINGLE_PRECISION
      RECURSIVE SUBROUTINE SortByPosition_s(xp, cutoff, skin, rank,clist,ownregion, &
 &                         idx, increment)
#elif __KIND == __DOUBLE_PRECISION
      RECURSIVE SUBROUTINE SortByPosition_d(xp, cutoff, skin, rank,clist,ownregion, &
 &                         idx, increment)
#endif
      !!! The recursive subroutine which sorts the particles by their position;
      !!! given particles coordinates, their cutoff radii, skin parameter,
      !!! rank array of particles, the domain that these particles take place,
      !!! index of current cell and increment parameter.
      !!! REMARK: increment parameter is essential for recursion, since it is
      !!! used to keep track of which chunk of particle array is being used.
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      REAL(MK),  INTENT(IN),    DIMENSION(:,:)        :: xp
      !!! Input array for particles coordinates
      REAL(MK),  INTENT(IN),    DIMENSION(:)          :: cutoff
      !!! Input array for particles cutoff radii
      REAL(MK),  INTENT(IN)                           :: skin
      !!! Skin parameter
      INTEGER, INTENT(INOUT), DIMENSION(:)            :: rank
      !!! rank array
      TYPE(ppm_clist), INTENT(INOUT)                  :: clist
      !!! cell list
      REAL(MK),  INTENT(IN),    DIMENSION(:)          :: ownregion
      !!! Region that will be used to sort particles within
      INTEGER(ppm_kind_int64),  INTENT(IN)            :: idx
      !!! Index of cell to be processed
      INTEGER,   INTENT(IN)                           :: increment
      !!! Start index in particles arrays.


      !-------------------------------------------------------------------------
      !  Local variables and counters
      !-------------------------------------------------------------------------
      INTEGER,   DIMENSION(1:2**ppm_dim)              :: incArray
      REAL(MK),  DIMENSION(1:2**ppm_dim, 1:2*ppm_dim) :: subregions
      REAL(MK),  DIMENSION(ppm_dim)                   :: mid_coor
      INTEGER                                         :: bx
      INTEGER                                         :: byLeft
      INTEGER                                         :: byRight
      INTEGER                                         :: bzBottomLeft
      INTEGER                                         :: bzBottomRight
      INTEGER                                         :: bzTopLeft
      INTEGER                                         :: bzTopRight
      REAL(MK)                                        :: minRC
      REAL(MK)                                        :: minSideLength
      INTEGER                                         :: i
      INTEGER                                         :: info

      ! If hash table is not sufficiently large, RETURN
      ! Use of grow_htable parameter globally ensures that
      ! every recursive subroutine will stop.

      IF(clist%grow_htable)   RETURN

      ! If no particles are assigned, return
      IF(size(rank) .LT. 1) RETURN

      ! Get maximum side length and the minimum cutoff radius
      minSideLength = getMinimumSideLength(ownregion)
      minRC = getMinimumRC(cutoff, skin, rank)

      ! If the cell is small enough, stop recursion for this cell.
      IF(minRC .GE. minSideLength) THEN
          RETURN
      END IF

      ! Set midpoint coordinates
      DO i = 1, ppm_dim
          mid_coor(i) = (ownregion(2*i-1) + ownregion(2*i))/2
      END DO

      ! Divide particles to 2, that are on the left of midpoint
      ! x-coordinate and on the right of it.
      CALL partition(xp, rank, mid_coor(1), bx, 1)

      ! Divide particles that are on the left of midpoint x-coordinate
      ! to 2, as on bottom of midpoint y-coordinate and on top of it.
      CALL partition(xp, rank(1:bx-1), mid_coor(2), byLeft,  2)

      ! Divide particles that are on the right of midpoint x-coordinate
      ! to 2, as on bottom of midpoint y-coordinate and on top of it.
      CALL partition(xp, rank(bx:),    mid_coor(2), byRight, 2)
      ! Increment border for upper-right portion
      byRight       = byRight       + bx      - 1

      ! If in 3D, keep on subdividing
      IF(ppm_dim .EQ. 3)  THEN
          CALL partition(xp, rank(1:byLeft-1),   mid_coor(3), bzBottomLeft,  3)
          CALL partition(xp, rank(byLeft:bx-1),  mid_coor(3), bzBottomRight, 3)
          CALL partition(xp, rank(bx:byRight-1), mid_coor(3), bzTopLeft,     3)
          CALL partition(xp, rank(byRight:),     mid_coor(3), bzTopRight,    3)

          ! Increment border indices accordingly.
          bzBottomRight = bzBottomRight + byLeft  - 1
          bzTopLeft     = bzTopLeft     + bx      - 1
          bzTopRight    = bzTopRight    + byRight - 1
      END IF

      ! Subdivide the current region, later to be assigned to recursive calls.
      CALL setSubregions(ownregion, subregions,info)

      ! Insert current cell in hash table.
      clist%borders_pos = clist%borders_pos + 1
      CALL hash_insert(clist%lookup,idx, clist%borders_pos, info)
      IF(info .NE. 0)   THEN
          clist%grow_htable = .TRUE.
          RETURN
      END IF

      ! Keep track of maximum number of cells, to be used later.
      IF(clist%borders_pos .GT. clist%borders_pos_max)   THEN
          clist%borders_pos_max = clist%borders_pos
      END IF

      ! Set last row as number of particles within the physical region
      ! of the cell, later to be used to understand whether the cell has
      ! particles in deeper levels or not.
      IF(ppm_dim .EQ. 2)       THEN
          clist%borders(6, clist%borders_pos) = size(rank)
      ELSEIF(ppm_dim .EQ. 3)   THEN
          clist%borders(10, clist%borders_pos) = size(rank)
      END IF

      ! Set increment parameters for recursive calls. This one is done
      ! genericly for 2D and 3D
      incArray(1) = increment
      incArray(  (ppm_dim-1) + 1) = increment + byLeft  - 1
      incArray(2*(ppm_dim-1) + 1) = increment + bx      - 1
      incArray(3*(ppm_dim-1) + 1) = increment + byRight - 1

      IF(ppm_dim .EQ. 3)   THEN
          incArray(2) = increment + bzBottomLeft  - 1
          incArray(4) = increment + bzBottomRight - 1
          incArray(6) = increment + bzTopLeft     - 1
          incArray(8) = increment + bzTopRight    - 1
      END IF

      ! Call recursive calls
      IF(ppm_dim .EQ. 2)       THEN
          CALL SortByPosition(xp, cutoff, skin, rank(1:byLeft-1), clist,      &
 &                    subregions(1,:), 4*idx-2, incArray(1))
          CALL SortByPosition(xp, cutoff, skin, rank(byLeft:bx-1), clist, &
 &                    subregions(2,:), 4*idx-1, incArray(2))
          CALL SortByPosition(xp, cutoff, skin, rank(bx:byRight-1), clist, &
 &                    subregions(3,:), 4*idx,   incArray(3))
          CALL SortByPosition(xp, cutoff, skin, rank(byRight:), clist, &
 &                    subregions(4,:), 4*idx+1, incArray(4))
      ELSEIF(ppm_dim .EQ. 3)   THEN
           CALL SortByPosition(xp, cutoff, skin, rank(1:bzBottomLeft-1),clist,  &
 &                    subregions(1,:), 8*idx-6, incArray(1))
           CALL SortByPosition(xp, cutoff, skin, rank(bzBottomLeft:byLeft-1),clist, &
 &                    subregions(2,:), 8*idx-5, incArray(2))
           CALL SortByPosition(xp, cutoff, skin, rank(byLeft:bzBottomRight-1),clist,&
 &                    subregions(3,:), 8*idx-4, incArray(3))
           CALL SortByPosition(xp, cutoff, skin, rank(bzBottomRight:bx-1),clist,&
 &                    subregions(4,:), 8*idx-3, incArray(4))
           CALL SortByPosition(xp, cutoff, skin, rank(bx:bzTopLeft-1),clist,  &
 &                    subregions(5,:), 8*idx-2, incArray(5))
           CALL SortByPosition(xp, cutoff, skin, rank(bzTopLeft:byRight-1),clist, &
 &                    subregions(6,:), 8*idx-1, incArray(6))
           CALL SortByPosition(xp, cutoff, skin, rank(byRight:bzTopRight-1),clist,&
 &                    subregions(7,:), 8*idx  , incArray(7))
           CALL SortByPosition(xp, cutoff, skin, rank(bzTopRight:),clist, &
 &                    subregions(8,:), 8*idx+1, incArray(8))
      END IF
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE SortByPosition_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE SortByPosition_d
#endif

#if   __KIND == __SINGLE_PRECISION
      RECURSIVE SUBROUTINE SortByRC_Pos_s(xp, cutoff, skin, rank,clist, ownregion, &
 &                         idx, level, increment)
#elif __KIND == __DOUBLE_PRECISION
      RECURSIVE SUBROUTINE SortByRC_Pos_d(xp, cutoff, skin, rank,clist, ownregion, &
 &                         idx, level, increment)
#endif
      !!! The recursive subroutine which sorts the particles by their position
      !!! and their cutoff radii; given particles coordinates, their cutoff
      !!! radii, skin parameter, rank array of particles, the domain that these
      !!! particles take place, index of current cell and increment parameter.
      !!! REMARK: increment parameter is essential for recursion, since it is
      !!! used to keep track of which chunk of particle array is being used.
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      REAL(MK),  INTENT(IN),    DIMENSION(:,:)        :: xp
      !!! Input array for particles coordinates
      REAL(MK),  INTENT(IN),    DIMENSION(:)          :: cutoff
      !!! Input array for particles cutoff radii
      REAL(MK),  INTENT(IN)                           :: skin
      !!! Skin parameter
      INTEGER, INTENT(INOUT), DIMENSION(:)            :: rank
      !!! rank array
      TYPE(ppm_clist), INTENT(INOUT)                  :: clist
      !!! cell list
      REAL(MK),  INTENT(IN),    DIMENSION(:)          :: ownregion
      !!! Region that will be used to sort particles within
      INTEGER(ppm_kind_int64),  INTENT(IN)            :: idx
      !!! Index of cell to be processed
      INTEGER,   INTENT(IN)                           :: level
      !!! Destination level that particles of this depth level
      !!! will be sorted.
      INTEGER,   INTENT(IN)                           :: increment
      !!! Start index in particles arrays.
      !---------------------------------------------------------------------
      !  Local variables and counters
      !---------------------------------------------------------------------
      INTEGER,   DIMENSION(2**ppm_dim)                   :: incArray
      REAL(MK),  DIMENSION(2**ppm_dim, 2*ppm_dim)        :: subregions
      REAL(MK),  DIMENSION(ppm_dim)                      :: mid_coor
      INTEGER                                            :: bx
      INTEGER                                            :: byLeft
      INTEGER                                            :: byRight
      INTEGER                                            :: bzBottomLeft
      INTEGER                                            :: bzBottomRight
      INTEGER                                            :: bzTopLeft
      INTEGER                                            :: bzTopRight
      INTEGER                                            :: i
      INTEGER                                            :: info

      ! If no particles are assigned, return
      IF(size(rank) .LT. 1) RETURN

      ! Set midpoint coordinates
      DO i = 1, ppm_dim
          mid_coor(i) = (ownregion(2*i-1) + ownregion(2*i))/2
      END DO

      ! Divide particles to 2, that are on the left of midpoint
      ! x-coordinate and on the right of it.
      CALL partition(xp, rank, mid_coor(1), bx, 1)

      ! Divide particles that are on the left of midpoint x-coordinate
      ! to 2, as on bottom of midpoint y-coordinate and on top of it.
      CALL partition(xp, rank(1:bx-1), mid_coor(2), byLeft,  2)

      ! Divide particles that are on the right of midpoint x-coordinate
      ! to 2, as on bottom of midpoint y-coordinate and on top of it.
      CALL partition(xp, rank(bx:),    mid_coor(2), byRight, 2)
      byRight       = byRight       + bx      - 1

      ! If in 3D, keep on subdividing
      IF(ppm_dim .EQ. 3)  THEN
          CALL partition(xp, rank(1:byLeft-1),   mid_coor(3), bzBottomLeft, 3)
          CALL partition(xp, rank(byLeft:bx-1),  mid_coor(3), bzBottomRight,3)
          CALL partition(xp, rank(bx:byRight-1), mid_coor(3), bzTopLeft,    3)
          CALL partition(xp, rank(byRight:),     mid_coor(3), bzTopRight,   3)

          bzBottomRight = bzBottomRight + byLeft  - 1
          bzTopLeft     = bzTopLeft     + bx      - 1
          bzTopRight    = bzTopRight    + byRight - 1
      END IF

      ! If the desired level is reached, update the entry on borders
      ! array, such that it contains border indices for children cells.
      IF(getCellDepth(idx) + 1 .GE. level)  THEN
          clist%borders_pos = hash_search(clist%lookup,idx)

          IF(ppm_dim .EQ. 2)       THEN
              clist%borders(1, clist%borders_pos) = increment
              clist%borders(2, clist%borders_pos) = byLeft        + increment - 1
              clist%borders(3, clist%borders_pos) = bx            + increment - 1
              clist%borders(4, clist%borders_pos) = byRight       + increment - 1
              clist%borders(5, clist%borders_pos) = size(rank)    + increment
          ELSEIF(ppm_dim .EQ. 3)   THEN
              clist%borders(1, clist%borders_pos) = increment
              clist%borders(2, clist%borders_pos) = bzBottomLeft  + increment - 1
              clist%borders(3, clist%borders_pos) = byLeft        + increment - 1
              clist%borders(4, clist%borders_pos) = bzBottomRight + increment - 1
              clist%borders(5, clist%borders_pos) = bx            + increment - 1
              clist%borders(6, clist%borders_pos) = bzTopLeft     + increment - 1
              clist%borders(7, clist%borders_pos) = byRight       + increment - 1
              clist%borders(8, clist%borders_pos) = bzTopRight    + increment - 1
              clist%borders(9, clist%borders_pos) = size(rank)    + increment
          END IF
          RETURN
      ! Else, keep on subdividing and calling recursive calls for
      ! children cells.
      ELSE
          CALL setSubregions(ownregion, subregions,info)
          incArray(1) = increment
          incArray(  (ppm_dim-1) + 1) = increment + byLeft  - 1
          incArray(2*(ppm_dim-1) + 1) = increment + bx      - 1
          incArray(3*(ppm_dim-1) + 1) = increment + byRight - 1

          IF(ppm_dim .EQ. 3)   THEN
              incArray(2) = increment + bzBottomLeft  - 1
              incArray(4) = increment + bzBottomRight - 1
              incArray(6) = increment + bzTopLeft     - 1
              incArray(8) = increment + bzTopRight    - 1
          END IF

          IF(ppm_dim .EQ. 2)  THEN
              CALL SortByRC_Pos(xp, cutoff, skin, rank(1:byLeft-1),clist, &
&                    subregions(1,:), 4*idx-2, level, incArray(1))
              CALL SortByRC_Pos(xp, cutoff, skin, rank(byLeft:bx-1),clist, &
&                    subregions(2,:), 4*idx-1, level, incArray(2))
              CALL SortByRC_Pos(xp, cutoff, skin, rank(bx:byRight-1),clist,&
&                    subregions(3,:), 4*idx, level, incArray(3))
              CALL SortByRC_Pos(xp, cutoff, skin, rank(byRight:),clist,&
&                    subregions(4,:), 4*idx+1, level, incArray(4))
          ELSEIF(ppm_dim .EQ. 3)   THEN
              CALL SortByRC_Pos(xp, cutoff, skin, rank(1:bzBottomLeft-1),clist, &
&                    subregions(1,:), 8*idx-6, level, incArray(1))
              CALL SortByRC_Pos(xp, cutoff, skin,rank(bzBottomLeft:byLeft-1),clist, &
&                    subregions(2,:), 8*idx-5, level, incArray(2))
              CALL SortByRC_Pos(xp, cutoff, skin,rank(byLeft:bzBottomRight-1),clist,&
&                    subregions(3,:), 8*idx-4, level, incArray(3))
              CALL SortByRC_Pos(xp, cutoff, skin, rank(bzBottomRight:bx-1),clist,&
&                    subregions(4,:), 8*idx-3, level, incArray(4))
              CALL SortByRC_Pos(xp, cutoff, skin, rank(bx:bzTopLeft-1),clist,  &
&                    subregions(5,:), 8*idx-2, level, incArray(5))
              CALL SortByRC_Pos(xp, cutoff, skin,rank(bzTopLeft:byRight-1),clist,&
&                    subregions(6,:), 8*idx-1, level, incArray(6))
              CALL SortByRC_Pos(xp, cutoff, skin,rank(byRight:bzTopRight-1),clist,  &
&                    subregions(7,:), 8*idx, level, incArray(7))
              CALL SortByRC_Pos(xp, cutoff, skin, rank(bzTopRight:),clist, &
&                    subregions(8,:), 8*idx+1, level, incArray(8))
          END IF
      END IF
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE SortByRC_Pos_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE SortByRC_Pos_d
#endif

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE partition_s(xp, rank, midPoint, border, axis)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE partition_d(xp, rank, midPoint, border, axis)
#endif
      !!! Given particles coordinates, rank array of particles, the midpoint
      !!! that particles will be partitioned by, axis to define on which axis
      !!! partitioning will take place and border to be updated such that
      !!! particles up to this border will be those before the midPoint and
      !!! after the border will be those after the midPoint.
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      REAL(MK), INTENT(IN),    DIMENSION(:,:) :: xp
      !!! Input array for particles coordinates
      INTEGER,  DIMENSION(:), INTENT(INOUT)   :: rank
      !!! ranks of particles
      REAL(MK), INTENT(IN)                    :: midPoint
      !!! Midpoint that particles will be partitioned by
      INTEGER,  INTENT(INOUT)                 :: border
      !!! Border index to be returned
      INTEGER,  INTENT(IN)                    :: axis
      !!! Axis on which partitioning will take place

      !---------------------------------------------------------------------
      !  Local variables and counters
      !---------------------------------------------------------------------
      INTEGER                                 :: i
      INTEGER                                 :: j
      INTEGER                                 :: temp

      IF(size(rank) .LT. 1) THEN
          border = 1
          RETURN
      END IF

      i= 0
      j= size(rank) + 1
      DO
          j = j - 1
          DO WHILE(j .GT. 0)
              IF(xp(axis, rank(j)) .LT. midPoint) exit
              j = j-1
          END DO

          i = i+1
          DO WHILE(i .LE. size(rank))
              IF(xp(axis, rank(i)) .GE. midPoint) exit
              i = i+1
          END DO

          IF (i .LT. j) THEN
              temp    = rank(i)
              rank(i) = rank(j)
              rank(j) = temp
          ELSEIF (i .EQ. j) THEN
              border = i
              RETURN
          ELSE
              border = i
              RETURN
          END IF
      END DO
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE partition_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE partition_d
#endif

#if   __KIND == __SINGLE_PRECISION
      RECURSIVE SUBROUTINE sortByRC_s(cutoff, skin, rank)
#elif __KIND == __DOUBLE_PRECISION
      RECURSIVE SUBROUTINE sortByRC_d(cutoff, skin, rank)
#endif
      !!! Recursive algorithm which sorts particles by their cutoff radii
      !!! in descending order; given cutoff radii, skin parameter and the
      !!! rank array (swapping is actually done on rank array).
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      REAL(MK), INTENT(IN),    DIMENSION(:) :: cutoff
      !!! Input array for particles cutoff radii
      REAL(MK), INTENT(IN)                  :: skin
      !!! Skin parameter
      INTEGER,  INTENT(INOUT), DIMENSION(:) :: rank
      !!! ranks of particles
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER                               :: marker

      ! Keep on partitioning and sorting, recursively.
      IF(size(rank) .GT. 1) THEN
         CALL partitionByRC(cutoff, skin, rank, marker)
         CALL sortByRC(cutoff, skin, rank(:marker-1))
         CALL sortByRC(cutoff, skin, rank(marker:))
      END IF
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE sortByRC_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE sortByRC_d
#endif

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE partitionByRC_s(cutoff, skin, rank, marker)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE partitionByRC_d(cutoff, skin, rank, marker)
#endif
      !!! This subroutine partitions particles by their cutoff radii and
      !!! updates marker; given cutoff radii, skin parameter and rank array.
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      REAL(MK), INTENT(IN),    DIMENSION(:) :: cutoff
      !!! Input array for particles cutoff radii
      REAL(MK), INTENT(IN)                  :: skin
      !!! Skin parameter
      INTEGER,  INTENT(INOUT), DIMENSION(:) :: rank
      !!! ranks of particles
      INTEGER,  INTENT(OUT)                 :: marker
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
      pivot = cutoff(rank(size(rank)/2)) + skin
      i= 0
      j= size(rank) + 1

      DO WHILE(i .LT. j)
         j = j - 1
         DO WHILE((cutoff(rank(j)) + skin) .LT. pivot)
            j = j-1
         END DO
         i = i + 1
         DO WHILE((cutoff(rank(i)) + skin) .GT. pivot)
            i = i + 1
         END DO
         IF (i .LT. j) THEN
            temp = rank(i)
            rank(i) = rank(j)
            rank(j) = temp
         END IF
      ENDDO


      IF (i .EQ. j) THEN
         marker = i + 1
      ELSE
         marker = i
      END IF
      
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE partitionByRC_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE partitionByRC_d
#endif

#if   __KIND == __SINGLE_PRECISION
      PURE FUNCTION lastIdxForRC_s(cutoff, skin, clist, inputRC) RESULT(idx)
#elif __KIND == __DOUBLE_PRECISION
      PURE FUNCTION lastIdxForRC_d(cutoff, skin, clist, inputRC) RESULT(idx)
#endif
      !!! For a given cutoff radius, returns the last index in the rank
      !!! array that is greater than or equal to this cutoff radius; given
      !!! cutoff radii of particles, skin parameter and input cutoff radius.
      !!! The result is necessary for sorting particles according to their
      !!! position for a specific level in SortByRC_Pos subroutine.
      !!! In the end, resulting index is returned by idx.
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      REAL(MK), INTENT(IN), DIMENSION(:) :: cutoff
      !!! Input array for particles cutoff radii
      REAL(MK), INTENT(IN)               :: skin
      !!! Skin parameter
      TYPE(ppm_clist), INTENT(IN)        :: clist
      !!! rank array
      REAL(MK), INTENT(IN)               :: inputRC
      !!! Input cutoff radius
      INTEGER                            :: idx
      !!! Result index

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER                            :: low
      INTEGER                            :: high
      INTEGER                            :: mid

      low  = 1
      high = size(clist%rank)
      mid = FLOOR(REAL((low + high)/2))

      DO WHILE((mid .NE. low) .OR. (mid .NE. high))
          IF((cutoff(clist%rank(mid)) + skin) .GE. inputRC)   THEN
              low = mid + 1
          ELSE
              high = mid
          END IF

          mid = FLOOR(REAL((low + high)/2))
      END DO

      idx = mid

      IF((mid .EQ. high) .AND. ((cutoff(clist%rank(mid))+skin) .GE. inputRC)) THEN
          idx = mid + 1
      END IF
#if   __KIND == __SINGLE_PRECISION
      END FUNCTION lastIdxForRC_s
#elif __KIND == __DOUBLE_PRECISION
      END FUNCTION lastIdxForRC_d
#endif

#if   __KIND == __SINGLE_PRECISION
       FUNCTION getMaxDepth_s(cutoff, clist, domain)    RESULT(depthMax)
#elif __KIND == __DOUBLE_PRECISION
       FUNCTION getMaxDepth_d(cutoff, clist, domain)    RESULT(depthMax)
#endif
      !!! This function returns the maximum depth within a domain
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      REAL(MK), INTENT(IN), DIMENSION(:)         :: cutoff
      TYPE(ppm_clist), INTENT(IN)                :: clist
      REAL(MK), INTENT(IN), DIMENSION(2*ppm_dim) :: domain
      INTEGER                                    :: depthMax

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      REAL(MK)                                   :: rc_min
      REAL(MK)                                   :: minSideLength

      rc_min = cutoff(clist%rank(size(clist%rank)))
      minSideLength = getMinimumSideLength(domain)
      depthMax = CEILING(LOG(minSideLength/rc_min)/LOG(2._MK))
#if   __KIND == __SINGLE_PRECISION
      END FUNCTION getMaxDepth_s
#elif __KIND == __DOUBLE_PRECISION
      END FUNCTION getMaxDepth_d
#endif

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE getRC_Borders_s(cutoff, skin, clist, domain,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE getRC_Borders_d(cutoff, skin, clist, domain,info)
#endif
      !!! Computes the borders on rc_borders array such that two
      !!! consecutive indices will contain borders for particle in that level.
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      REAL(MK), INTENT(IN), DIMENSION(:) :: cutoff
      !!! Input array for particles cutoff radii
      REAL(MK), INTENT(IN)               :: skin
      !!! Skin parameter
      TYPE(ppm_clist), INTENT(INOUT)     :: clist
      !!! the cell list
      REAL(MK), DIMENSION(2*ppm_dim)     :: domain
      !!! Physical extent of whole domain including ghost layers
      INTEGER                            :: info

      !---------------------------------------------------------------------
      !  Local variables and counters
      !---------------------------------------------------------------------
      INTEGER                            :: i
      REAL(MK)                           :: rc_limit
      REAL(MK)                           :: minSideLength
      INTEGER                            :: rc_border
      
      REAL(MK)                           :: t0
      
      CALL substart('getRC_Borders',t0,info)

      minSideLength = getMinimumSideLength(domain)

      DO i = 1, clist%max_depth
          ! minimum cutoff radius has to be greater than
          ! half of the maximum side length
          rc_limit  = minSideLength/2
          rc_border = lastIdxForRC(cutoff, skin, clist, rc_limit)
          clist%rc_borders(i) = rc_border
          ! maximum side length after a subdivision
          minSideLength = minSideLength/2
      END DO
      
      CALL substop('getRC_Borders',t0,info)
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE getRC_Borders_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE getRC_Borders_d
#endif
