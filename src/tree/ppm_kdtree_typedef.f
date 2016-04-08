      !----------------------------------------------------------------------
      ! K-D tree routines in Fortran 90 by:
      ! Matthew B. Kennel, Institute For Nonlinear Science,
      ! reference: http://arxiv.org/abs/physics/0408067
      ! Licensed under the Academic Free License version 1.1.
      !
      ! It has been adapted and amended for PPM library by Yaser Afshar.
      !
      ![NOTE]
      ! In this adaptation the maximum query vector size has been set to three
      ! in case of future development and higher dimensionality this routine
      ! needs to be modified
      !----------------------------------------------------------------------
      !
      ! maintain a priority queue (PQ) of data, pairs of 'priority/payload',
      ! implemented with a binary heap.  This is the type, and the 'dis' field
      ! is the priority.
      !
      TYPE :: DTYPE(kdtree_result)
        REAL(MK) :: dis
        INTEGER  :: idx
      END TYPE DTYPE(kdtree_result)

#if __KIND == __SINGLE_PRECISION
      !
      ! A heap-based priority queue lets one efficiently implement the following
      ! operations, each in log(N) time, as opposed to linear time.
      !
      ! 1)  add a datum (push a datum onto the queue, increasing its length)
      ! 2)  return the priority value of the maximum priority element
      ! 3)  pop-off (and delete) the element with the maximum priority, decreasing
      !     the size of the queue.
      ! 4)  replace the datum with the maximum priority with a supplied datum
      !     (of either higher or lower priority), maintaining the size of the
      !     queue.
      !
      !
      ! In the k-d tree case, the 'priority' is the square distance of a point in
      ! the data set to a reference point.   The goal is to keep the smallest M
      ! distances to a reference point.  The tree algorithm searches terminal
      ! nodes to decide whether to add points under consideration.
      !
      ! A priority queue is useful here because it lets one quickly return the
      ! largest distance currently existing in the list.  If a new candidate
      ! distance is smaller than this, then the new candidate ought to replace
      ! the old candidate.  In priority queue terms, this means removing the
      ! highest priority element, and inserting the new one.
      !
      ! Algorithms based on Cormen, Leiserson, Rivest, _Introduction
      ! to Algorithms_, 1990, with further optimization by the Kennel.
      !
      ! This module is not written in the most clear way, but is implemented such
      ! for speed, as it its operations will be called many times during searches
      ! of large numbers of neighbors.
      !
      TYPE :: pq
        INTEGER                                      :: heap_size = 0
        ! The priority queue consists of elements
        ! priority(1:heap_size), with associated payload(:).
        !
        ! There are heap_size active elements.
        ! Assumes the allocation is always sufficient.
        ! Will NOT increase it to match.

        TYPE(kdtree_result_s), DIMENSION(:), POINTER :: elems_s => NULL()
        TYPE(kdtree_result_d), DIMENSION(:), POINTER :: elems_d => NULL()
      CONTAINS
        PROCEDURE :: pq_create_s
        PROCEDURE :: pq_create_d
        PROCEDURE :: pq_insert_s
        PROCEDURE :: pq_insert_d
        PROCEDURE :: pq_delete
        PROCEDURE :: pq_extract_max_s
        PROCEDURE :: pq_extract_max_d
        PROCEDURE :: pq_max_s
        PROCEDURE :: pq_max_d
        PROCEDURE :: pq_replace_max_s
        PROCEDURE :: pq_replace_max_d

        GENERIC   :: create      => pq_create_s,pq_create_d
        GENERIC   :: insert      => pq_insert_s,pq_insert_d
        GENERIC   :: delete      => pq_delete
        GENERIC   :: extract_max => pq_extract_max_s,pq_extract_max_d
        GENERIC   :: max         => pq_max_s,pq_max_d
        GENERIC   :: replace_max => pq_replace_max_s,pq_replace_max_d

      END TYPE pq
#endif

      TYPE :: DTYPE(interval)
        REAL(MK) :: lower
        REAL(MK) :: upper
      END TYPE DTYPE(interval)

      TYPE :: DTYPE(tree_node)
        ! an internal tree node
        INTEGER                                       :: cut_dim
        ! the dimension to cut
        REAL(MK)                                      :: cut_val
        ! where to cut the dimension
        REAL(MK)                                      :: cut_val_left
        REAL(MK)                                      :: cut_val_right
        ! improved cutoffs knowing the spread in child boxes.

        INTEGER                                       :: l
        INTEGER                                       :: u

        TYPE(DTYPE(tree_node)),               POINTER :: left  => NULL()
        TYPE(DTYPE(tree_node)),               POINTER :: right => NULL()
        TYPE(DTYPE(interval)),  DIMENSION(:), POINTER :: box   => NULL()
        ! child pointers
        ! Points included in this node are indexes[k] with k \in [l,u]

      END TYPE DTYPE(tree_node)

      TYPE :: DTYPE(kdtree)
        ! Global information about the tree, one per tree
        INTEGER                           :: dimen=0
        INTEGER                           :: n=0
        ! dimensionality and total # of points

        REAL(MK), DIMENSION(:,:), POINTER :: the_data => NULL()
        ! pointer to the actual data array
        !
        !  IMPORTANT NOTE:  IT IS DIMENSIONED   the_data(1:d,1:N)
        !  The search time is dominated by the
        !  evaluation of distances in the terminal nodes.  Putting all
        !  vector components in consecutive memory location improves
        !  memory cache locality, and hence search speed, and may enable
        !  vectorization on some processors and compilers.

        INTEGER,  DIMENSION(:),   POINTER :: ind => NULL()
        !  permuted index into the data, so that indexes[l..u] of some
        ! bucket represent the indexes of the actual points in that
        ! bucket.

        LOGICAL                           :: sort      = .FALSE.
        ! do we always sort output results?
        LOGICAL                           :: rearrange = .FALSE.

        REAL(MK), DIMENSION(:,:), POINTER :: rearranged_data => NULL()
        ! if (rearrange .eqv. .TRUE.) THEN rearranged_data has been
        ! created so that rearranged_data(:,i) = the_data(:,ind(i)),
        ! permitting search to use more cache-friendly rearranged_data, at
        ! some initial computation and storage cost.

        TYPE(DTYPE(tree_node)),   POINTER :: root => NULL()
        ! root pointer of the tree

      CONTAINS
        PROCEDURE :: create                     => DTYPE(kdtree_create)
        PROCEDURE :: destroy                    => DTYPE(kdtree_destroy)
        PROCEDURE :: destroy_node               => DTYPE(kdtree_destroy_node)
        PROCEDURE :: build                      => DTYPE(build_tree)
        PROCEDURE :: build_range                => DTYPE(build_tree_for_range)
        PROCEDURE :: spread_in_coordinate       => DTYPE(spread_in_coordinate)
        PROCEDURE :: select_on_coordinate_value => DTYPE(select_on_coordinate_value)

      END TYPE DTYPE(kdtree)

      TYPE :: DTYPE(tree_search_record)
        !
        ! One of these is created for each search.
        !
        !
        ! Many fields are copied from the tree structure, in order to
        ! speed up the search.
        !
        INTEGER                                             :: dimen
        INTEGER                                             :: nn
        INTEGER                                             :: nfound
        INTEGER                                             :: centeridx=999
        INTEGER                                             :: correltime=9999
        ! exclude points within 'correltime' of 'centeridx', if centeridx >= 0
        INTEGER                                             :: nalloc
        ! how much allocated for results(:)?
        INTEGER,                    DIMENSION(:),   POINTER :: ind => NULL()
        ! temp pointer to indexes

        REAL(MK)                                            :: ballsize
        REAL(MK),                   DIMENSION(3)            :: qv
        ! query vector (at PPM we do not have dimension more than 3)
        REAL(MK),                   DIMENSION(:,:), POINTER :: data => NULL()
        ! temp pointer to data

        TYPE(DTYPE(kdtree_result)), DIMENSION(:),   POINTER :: results => NULL()
        ! results
        TYPE(pq)                                            :: pq

        LOGICAL                                             :: rearrange
        ! are the data rearranged or original?
        LOGICAL                                             :: overflow
        ! did the # of points found overflow the storage provided?
      END TYPE DTYPE(tree_search_record)

#undef  __KIND
#undef    MK
#undef   _MK
#undef    DTYPE
