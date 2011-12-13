  type DTYPE(interval)
      real(kdkind) :: lower,upper
  end type DTYPE(interval)

  type :: DTYPE(tree_node)
      ! an internal tree node
      private
      integer :: cut_dim
      ! the dimension to cut
      real(kdkind) :: cut_val
      ! where to cut the dimension
      real(kdkind) :: cut_val_left, cut_val_right  
      ! improved cutoffs knowing the spread in child boxes.
      integer :: l, u
      type (DTYPE(tree_node)), pointer :: left, right
      type(DTYPE(interval)), pointer :: box(:) => null()
      ! child pointers
      ! Points included in this node are indexes[k] with k \in [l,u] 
  end type DTYPE(tree_node)

  type :: DTYPE(kdtree2)
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
      type (DTYPE(tree_node)), pointer :: root => null()
      ! root pointer of the tree
  end type DTYPE(kdtree2)


  type :: DTYPE(tree_search_record)
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
      type(DTYPE(kdtree2_result)), pointer :: results(:) ! results
      type(DTYPE(pq)) :: pq
      real(kdkind), pointer :: data(:,:)  ! temp pointer to data
      integer, pointer      :: ind(:)     ! temp pointer to indexes
  end type DTYPE(tree_search_record)
