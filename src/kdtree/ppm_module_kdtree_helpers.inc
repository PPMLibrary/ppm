  type DTYPE(kdtree2_result)
      ! a pair of distances, indexes
      real(kdkind)    :: dis!=0.0
      integer :: idx!=-1   Initializers cause some bugs in compilers.
  end type DTYPE(kdtree2_result)
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
  ! to Algorithms_, 1990, with further optimization by the author.
  !
  ! Originally informed by a C implementation by Sriranga Veeraraghavan.
  !
  ! This module is not written in the most clear way, but is implemented such
  ! for speed, as it its operations will be called many times during searches
  ! of large numbers of neighbors.
  !
  type DTYPE(pq)
      !
      ! The priority queue consists of elements
      ! priority(1:heap_size), with associated payload(:).
      !
      ! There are heap_size active elements. 
      ! Assumes the allocation is always sufficient.  Will NOT increase it
      ! to match.
      integer :: heap_size = 0
      type(DTYPE(kdtree2_result)), pointer :: elems(:) 
  end type DTYPE(pq)

  public :: DTYPE(kdtree2_result)

  public :: DTYPE(pq)

