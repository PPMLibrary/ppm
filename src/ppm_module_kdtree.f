module ppm_module_kdtree
  USE ppm_module_typedef
  USE ppm_module_kdtree_helpers
  INTEGER, PARAMETER :: kdkind = ppm_kind_double
  ! K-D tree routines in Fortran 90 by Matt Kennel.
  ! Original program was written in Sather by Steve Omohundro and
  ! Matt Kennel.  Only the Euclidean metric is supported. 
  !
  !
  ! This module is identical to 'kd_tree', except that the order
  ! of subscripts is reversed in the data file.
  ! In otherwords for an embedding of N D-dimensional vectors, the
  ! data file is here, in natural Fortran order  data(1:D, 1:N)
  ! because Fortran lays out columns first,
  !
  ! whereas conventionally (C-style) it is data(1:N,1:D)
  ! as in the original kd_tree module. 
  !
  !-------------DATA TYPE, CREATION, DELETION---------------------
  public :: kdtree2, kdtree2_result, tree_node, kdtree2_create, kdtree2_destroy
  !---------------------------------------------------------------
  !-------------------SEARCH ROUTINES-----------------------------
  public :: kdtree2_n_nearest,kdtree2_n_nearest_around_point
  ! Return fixed number of nearest neighbors around arbitrary vector,
  ! or extant point in dataset, with decorrelation window. 
  !
  public :: kdtree2_r_nearest, kdtree2_r_nearest_around_point
  ! Return points within a fixed ball of arb vector/extant point 
  !
  public :: kdtree2_sort_results
  ! Sort, in order of increasing distance, rseults from above.
  !
  public :: kdtree2_r_count, kdtree2_r_count_around_point 
  ! Count points within a fixed ball of arb vector/extant point 
  !
  public :: kdtree2_n_nearest_brute_force, kdtree2_r_nearest_brute_force
  ! brute force of kdtree2_[n|r]_nearest
  !----------------------------------------------------------------


  integer, parameter :: bucket_size = 12
  ! The maximum number of points to keep in a terminal node.

  type interval
      real(kdkind) :: lower,upper
  end type interval

  type :: tree_node
      ! an internal tree node
      private
      integer :: cut_dim
      ! the dimension to cut
      real(kdkind) :: cut_val
      ! where to cut the dimension
      real(kdkind) :: cut_val_left, cut_val_right  
      ! improved cutoffs knowing the spread in child boxes.
      integer :: l, u
      type (tree_node), pointer :: left, right
      type(interval), pointer :: box(:) => null()
      ! child pointers
      ! Points included in this node are indexes[k] with k \in [l,u] 


  end type tree_node

  type :: kdtree2
      ! Global information about the tree, one per tree
      integer :: dimen=0, n=0
      ! dimensionality and total # of points
      real(kdkind), pointer :: the_data(:,:) => null()
      ! pointer to the actual data array 
      ! 
      !  IMPORTANT NOTE:  IT IS DIMENSIONED   the_data(1:d,1:N)
      !  which may be opposite of what may be conventional.
      !  This is, because in Fortran, the memory layout is such that
      !  the first dimension is in sequential order.  Hence, with
      !  (1:d,1:N), all components of the vector will be in consecutive
      !  memory locations.  The search time is dominated by the
      !  evaluation of distances in the terminal nodes.  Putting all
      !  vector components in consecutive memory location improves
      !  memory cache locality, and hence search speed, and may enable 
      !  vectorization on some processors and compilers. 

      integer, pointer :: ind(:) => null()
      ! permuted index into the data, so that indexes[l..u] of some
      ! bucket represent the indexes of the actual points in that
      ! bucket.
      logical       :: sort = .false.
      ! do we always sort output results?
      logical       :: rearrange = .false. 
      real(kdkind), pointer :: rearranged_data(:,:) => null()
      ! if (rearrange .eqv. .true.) then rearranged_data has been
      ! created so that rearranged_data(:,i) = the_data(:,ind(i)),
      ! permitting search to use more cache-friendly rearranged_data, at
      ! some initial computation and storage cost.
      type (tree_node), pointer :: root => null()
      ! root pointer of the tree
  end type kdtree2


  type :: tree_search_record
      !
      ! One of these is created for each search.
      !
      private
      ! 
      ! Many fields are copied from the tree structure, in order to
      ! speed up the search.
      !
      integer           :: dimen   
      integer           :: nn, nfound
      real(kdkind)      :: ballsize
      integer           :: centeridx=999, correltime=9999
      ! exclude points within 'correltime' of 'centeridx', iff centeridx >= 0
      integer           :: nalloc  ! how much allocated for results(:)?
      logical           :: rearrange  ! are the data rearranged or original? 
      ! did the # of points found overflow the storage provided?
      logical           :: overflow
      real(kdkind), pointer :: qv(:)  ! query vector
      type(kdtree2_result), pointer :: results(:) ! results
      type(pq) :: pq
      real(kdkind), pointer :: data(:,:)  ! temp pointer to data
      integer, pointer      :: ind(:)     ! temp pointer to indexes
  end type tree_search_record

  private
  ! everything else is private.

  type(tree_search_record), save, target :: sr   ! A GLOBAL VARIABLE for search

contains

include 'kdtree/ppm_kdtree.f'

end module ppm_module_kdtree

