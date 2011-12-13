  function DTYPE(kdtree2_create)(input_data,dim,sort,rearrange) result (mr)
    !
    ! create the actual tree structure, given an input array of data.
    !
    ! Note, input data is input_data(1:d,1:N), NOT the other way around.
    ! THIS IS THE REVERSE OF THE PREVIOUS VERSION OF THIS MODULE.
    ! The reason for it is cache friendliness, improving performance.
    !
    ! Optional arguments:  If 'dim' is specified, then the tree
    !                      will only search the first 'dim' components
    !                      of input_data, otherwise, dim is inferred
    !                      from SIZE(input_data,1).
    !
    !                      if sort .eqv. .true. then output results
    !                      will be sorted by increasing distance.
    !                      default=.false., as it is faster to not sort.
    !                      
    !                      if rearrange .eqv. .true. then an internal
    !                      copy of the data, rearranged by terminal node,
    !                      will be made for cache friendliness. 
    !                      default=.true., as it speeds searches, but
    !                      building takes longer, and extra memory is used.
    !
    ! .. Function Return Cut_value ..
    type (DTYPE(kdtree2)), pointer :: mr
    integer, intent(in), optional      :: dim
    logical, intent(in), optional      :: sort
    logical, intent(in), optional      :: rearrange
    ! ..
    ! .. Array Arguments ..
    real(kdkind), target :: input_data(:,:)
    !
    integer :: i
    ! ..
    allocate (mr)
    mr%the_data => input_data
    ! pointer assignment

    if (present(dim)) then
       mr%dimen = dim
    else
       mr%dimen = size(input_data,1)
    end if
    mr%n = size(input_data,2)

    if (mr%dimen > mr%n) then
       !  unlikely to be correct
       write (*,*) 'KD_TREE_TRANS: likely user error.'
       write (*,*) 'KD_TREE_TRANS: You passed in matrix with D=',mr%dimen
       write (*,*) 'KD_TREE_TRANS: and N=',mr%n
       write (*,*) 'KD_TREE_TRANS: note, that new format is data(1:D,1:N)'
       write (*,*) 'KD_TREE_TRANS: with usually N >> D.'
       write (*,*) '               If N =approx= D, then a k-d tree'
       write (*,*) 'KD_TREE_TRANS: is not an appropriate data structure.'
       stop
    end if

    call DTYPE(build_tree)(mr)

    if (present(sort)) then
       mr%sort = sort
    else
       mr%sort = .false.
    endif

    if (present(rearrange)) then
       mr%rearrange = rearrange
    else
       mr%rearrange = .true.
    endif

    if (mr%rearrange) then
       allocate(mr%rearranged_data(mr%dimen,mr%n))
       do i=1,mr%n
          mr%rearranged_data(:,i) = mr%the_data(:, &
           mr%ind(i))
       enddo
    else
       nullify(mr%rearranged_data)
    endif

  end function DTYPE(kdtree2_create)

    subroutine DTYPE(build_tree)(tp)
      type (DTYPE(kdtree2)), pointer :: tp
      ! ..
      integer :: j
      type(DTYPE(tree_node)), pointer :: dummy => null()
      ! ..
      allocate (tp%ind(tp%n))
      forall (j=1:tp%n)
         tp%ind(j) = j
      end forall
      tp%root => DTYPE(build_tree_for_range)(tp,1,tp%n, dummy)
    end subroutine DTYPE(build_tree)

    recursive function DTYPE(build_tree_for_range)(tp,l,u,parent) result (res)
      ! .. Function Return Cut_value ..
      type (DTYPE(tree_node)), pointer :: res
      ! ..
      ! .. Structure Arguments ..
      type (DTYPE(kdtree2)), pointer :: tp
      type (DTYPE(tree_node)),pointer           :: parent
      ! ..
      ! .. Scalar Arguments ..
      integer, intent (In) :: l, u
      ! ..
      ! .. Local Scalars ..
      integer :: i, c, m, dimen
      logical :: recompute
      real(kdkind)    :: average

!!$      If (.False.) Then 
!!$         If ((l .Lt. 1) .Or. (l .Gt. tp%n)) Then
!!$            Stop 'illegal L value in build_tree_for_range'
!!$         End If
!!$         If ((u .Lt. 1) .Or. (u .Gt. tp%n)) Then
!!$            Stop 'illegal u value in build_tree_for_range'
!!$         End If
!!$         If (u .Lt. l) Then
!!$            Stop 'U is less than L, thats illegal.'
!!$         End If
!!$      Endif
!!$      
      ! first compute min and max
      dimen = tp%dimen
      allocate (res)
      allocate(res%box(dimen))

      ! First, compute an APPROXIMATE bounding box of all points associated with this node.
      if ( u < l ) then
         ! no points in this box
         nullify(res)
         return
      end if

      if ((u-l)<=bucket_size) then
         !
         ! always compute true bounding box for terminal nodes.
         !
         do i=1,dimen
            call DTYPE(spread_in_coordinate)(tp,i,l,u,res%box(i))
         end do
         res%cut_dim = 0
         res%cut_val = 0.0
         res%l = l
         res%u = u
         res%left =>null()
         res%right => null() 
      else
         ! 
         ! modify approximate bounding box.  This will be an
         ! overestimate of the true bounding box, as we are only recomputing 
         ! the bounding box for the dimension that the parent split on.
         !
         ! Going to a true bounding box computation would significantly
         ! increase the time necessary to build the tree, and usually
         ! has only a very small difference.  This box is not used
         ! for searching but only for deciding which coordinate to split on.
         !
         do i=1,dimen
            recompute=.true.
            if (associated(parent)) then
               if (i .ne. parent%cut_dim) then
                  recompute=.false.
               end if
            endif
            if (recompute) then
               call DTYPE(spread_in_coordinate)(tp,i,l,u,res%box(i))
            else
               res%box(i) = parent%box(i)
            endif
         end do
         

         c = maxloc(res%box(1:dimen)%upper-res%box(1:dimen)%lower,1)
         !
         ! c is the identity of which coordinate has the greatest spread.
         !
         
         if (.false.) then
            ! select exact median to have fully balanced tree.
            m = (l+u)/2
            call DTYPE(select_on_coordinate)(tp%the_data,tp%ind,c,m,l,u)
         else
            !
            ! select point halfway between min and max, as per A. Moore,
            ! who says this helps in some degenerate cases, or 
            ! actual arithmetic average. 
            !
            if (.true.) then
               ! actually compute average
               average = sum(tp%the_data(c,tp%ind(l:u))) / real(u-l+1,kdkind)
            else
               average = (res%box(c)%upper + res%box(c)%lower)/2.0
            endif
               
            res%cut_val = average
            m = DTYPE(select_on_coordinate_value)(tp%the_data,tp%ind,c,average,l,u)
         endif
            
         ! moves indexes around
         res%cut_dim = c
         res%l = l
         res%u = u
!         res%cut_val = tp%the_data(c,tp%ind(m))

         res%left => DTYPE(build_tree_for_range)(tp,l,m,res)
         res%right => DTYPE(build_tree_for_range)(tp,m+1,u,res)

         if (associated(res%right) .eqv. .false.) then
            res%box = res%left%box
            res%cut_val_left = res%left%box(c)%upper
            res%cut_val = res%cut_val_left
         elseif (associated(res%left) .eqv. .false.) then
            res%box = res%right%box
            res%cut_val_right = res%right%box(c)%lower
            res%cut_val = res%cut_val_right
         else
            res%cut_val_right = res%right%box(c)%lower
            res%cut_val_left = res%left%box(c)%upper
            res%cut_val = (res%cut_val_left + res%cut_val_right)/2


            ! now remake the true bounding box for self.  
            ! Since we are taking unions (in effect) of a tree structure,
            ! this is much faster than doing an exhaustive
            ! search over all points
            res%box%upper = max(res%left%box%upper,res%right%box%upper)
            res%box%lower = min(res%left%box%lower,res%right%box%lower) 
         endif
      end if
    end function DTYPE(build_tree_for_range)

    integer function DTYPE(select_on_coordinate_value)(v,ind,c,alpha,li,ui) &
     result(res)
      ! Move elts of ind around between l and u, so that all points
      ! <= than alpha (in c cooordinate) are first, and then
      ! all points > alpha are second. 

      !
      ! Algorithm (matt kennel). 
      !
      ! Consider the list as having three parts: on the left,
      ! the points known to be <= alpha.  On the right, the points
      ! known to be > alpha, and in the middle, the currently unknown
      ! points.   The algorithm is to scan the unknown points, starting
      ! from the left, and swapping them so that they are added to
      ! the left stack or the right stack, as appropriate.
      ! 
      ! The algorithm finishes when the unknown stack is empty. 
      !
      ! .. Scalar Arguments ..
      integer, intent (In) :: c, li, ui
      real(kdkind), intent(in) :: alpha
      ! ..
      real(kdkind) :: v(1:,1:)
      integer :: ind(1:)
      integer :: tmp  
      ! ..
      integer :: lb, rb
      !
      ! The points known to be <= alpha are in
      ! [l,lb-1]
      !
      ! The points known to be > alpha are in
      ! [rb+1,u].  
      !
      ! Therefore we add new points into lb or
      ! rb as appropriate.  When lb=rb
      ! we are done.  We return the location of the last point <= alpha.
      !
      ! 
      lb = li; rb = ui

      do while (lb < rb)
         if ( v(c,ind(lb)) <= alpha ) then
            ! it is good where it is.
            lb = lb+1
         else
            ! swap it with rb.
            tmp = ind(lb); ind(lb) = ind(rb); ind(rb) = tmp
            rb = rb-1
         endif
      end do
      
      ! now lb .eq. ub 
      if (v(c,ind(lb)) <= alpha) then
         res = lb
      else
         res = lb-1
      endif
      
    end function DTYPE(select_on_coordinate_value)

    subroutine DTYPE(select_on_coordinate)(v,ind,c,k,li,ui)
      ! Move elts of ind around between l and u, so that the kth
      ! element
      ! is >= those below, <= those above, in the coordinate c.
      ! .. Scalar Arguments ..
      integer, intent (In) :: c, k, li, ui
      ! ..
      integer :: i, l, m, s, t, u
      ! ..
      real(kdkind) :: v(:,:)
      integer :: ind(:)
      ! ..
      l = li
      u = ui
      do while (l<u)
         t = ind(l)
         m = l
         do i = l + 1, u
            if (v(c,ind(i))<v(c,t)) then
               m = m + 1
               s = ind(m)
               ind(m) = ind(i)
               ind(i) = s
            end if
         end do
         s = ind(l)
         ind(l) = ind(m)
         ind(m) = s
         if (m<=k) l = m + 1
         if (m>=k) u = m - 1
      end do
    end subroutine DTYPE(select_on_coordinate)

   subroutine DTYPE(spread_in_coordinate)(tp,c,l,u,interv) 
      ! the spread in coordinate 'c', between l and u. 
      !
      ! Return lower bound in 'smin', and upper in 'smax', 
      ! ..
      ! .. Structure Arguments ..
      type (DTYPE(kdtree2)), pointer :: tp
      type(DTYPE(interval)), intent(out) :: interv
      ! ..
      ! .. Scalar Arguments ..
      integer, intent (In) :: c, l, u
      ! ..
      ! .. Local Scalars ..
      real(kdkind) :: last, lmax, lmin, t, smin,smax
      integer :: i, ulocal
      ! ..
      ! .. Local Arrays ..
      real(kdkind), pointer :: v(:,:)
      integer, pointer :: ind(:)
      ! ..
      v => tp%the_data(1:,1:)
      ind => tp%ind(1:)
      smin = v(c,ind(l))
      smax = smin

      ulocal = u

      do i = l + 2, ulocal, 2
         lmin = v(c,ind(i-1))
         lmax = v(c,ind(i))
         if (lmin>lmax) then
            t = lmin
            lmin = lmax
            lmax = t
         end if
         if (smin>lmin) smin = lmin
         if (smax<lmax) smax = lmax
      end do
      if (i==ulocal+1) then
         last = v(c,ind(ulocal))
         if (smin>last) smin = last
         if (smax<last) smax = last
      end if

      interv%lower = smin
      interv%upper = smax

    end subroutine DTYPE(spread_in_coordinate)


  subroutine DTYPE(kdtree2_destroy)(tp)
    ! Deallocates all memory for the tree, except input data matrix
    ! .. Structure Arguments ..
    type (DTYPE(kdtree2)), pointer :: tp
    ! ..
    call DTYPE(destroy_node)(tp%root)

    deallocate (tp%ind)
    nullify (tp%ind)

    if (tp%rearrange) then
       deallocate(tp%rearranged_data)
       nullify(tp%rearranged_data)
    endif

    deallocate(tp)
    return

  contains
    recursive subroutine DTYPE(destroy_node)(np)
      ! .. Structure Arguments ..
      type (DTYPE(tree_node)), pointer :: np
      ! ..
      ! .. Intrinsic Functions ..
      intrinsic ASSOCIATED
      ! ..
      if (associated(np%left)) then
         call DTYPE(destroy_node)(np%left)
         nullify (np%left)
      end if
      if (associated(np%right)) then
         call DTYPE(destroy_node)(np%right)
         nullify (np%right)
      end if
      if (associated(np%box)) deallocate(np%box)
      deallocate(np)
      return
      
    end subroutine DTYPE(destroy_node)

  end subroutine DTYPE(kdtree2_destroy)

  subroutine DTYPE(kdtree2_n_nearest)(tp,qv,nn,results)
    ! Find the 'nn' vectors in the tree nearest to 'qv' in euclidean norm
    ! returning their indexes and distances in 'indexes' and 'distances'
    ! arrays already allocated passed to this subroutine.
    type (DTYPE(kdtree2)), pointer      :: tp
    real(kdkind), target, intent (In)    :: qv(:)
    integer, intent (In)         :: nn
    type(DTYPE(kdtree2_result)), target :: results(:)


    sr%ballsize = huge(1.0)
    sr%qv => qv
    sr%nn = nn
    sr%nfound = 0
    sr%centeridx = -1
    sr%correltime = 0
    sr%overflow = .false. 

    sr%results => results

    sr%nalloc = nn   ! will be checked

    sr%ind => tp%ind
    sr%rearrange = tp%rearrange
    if (tp%rearrange) then
       sr%Data => tp%rearranged_data
    else
       sr%Data => tp%the_data
    endif
    sr%dimen = tp%dimen

    call DTYPE(validate_query_storage)(nn) 
    sr%pq = pq_create(results)

    call DTYPE(search)(tp%root)

    if (tp%sort) then
       call kdtree2_sort_results(nn, results)
    endif
!    deallocate(sr%pqp)
    return
  end subroutine DTYPE(kdtree2_n_nearest)

  subroutine DTYPE(kdtree2_n_nearest_around_point)(tp,idxin,correltime,nn,results)
    ! Find the 'nn' vectors in the tree nearest to point 'idxin',
    ! with correlation window 'correltime', returing results in
    ! results(:), which must be pre-allocated upon entry.
    type (DTYPE(kdtree2)), pointer        :: tp
    integer, intent (In)           :: idxin, correltime, nn
    type(DTYPE(kdtree2_result)), target   :: results(:)

    allocate (sr%qv(tp%dimen))
    sr%qv = tp%the_data(:,idxin) ! copy the vector
    sr%ballsize = huge(1.0)       ! the largest real(kdkind) number
    sr%centeridx = idxin
    sr%correltime = correltime

    sr%nn = nn
    sr%nfound = 0

    sr%dimen = tp%dimen
    sr%nalloc = nn

    sr%results => results

    sr%ind => tp%ind
    sr%rearrange = tp%rearrange

    if (sr%rearrange) then
       sr%Data => tp%rearranged_data
    else
       sr%Data => tp%the_data
    endif

    call DTYPE(validate_query_storage)(nn)
    sr%pq = pq_create(results)

    call DTYPE(search)(tp%root)

    if (tp%sort) then
       call kdtree2_sort_results(nn, results)
    endif
    deallocate (sr%qv)
    return
  end subroutine DTYPE(kdtree2_n_nearest_around_point)

  subroutine DTYPE(kdtree2_r_nearest)(tp,qv,r2,nfound,nalloc,results) 
    ! find the nearest neighbors to point 'idxin', within SQUARED
    ! Euclidean distance 'r2'.   Upon ENTRY, nalloc must be the
    ! size of memory allocated for results(1:nalloc).  Upon
    ! EXIT, nfound is the number actually found within the ball. 
    !
    !  Note that if nfound .gt. nalloc then more neighbors were found
    !  than there were storage to store.  The resulting list is NOT
    !  the smallest ball inside norm r^2 
    !
    ! Results are NOT sorted unless tree was created with sort option.
    type (DTYPE(kdtree2)), pointer      :: tp
    real(kdkind), target, intent (In)    :: qv(:)
    real(kdkind), intent(in)             :: r2
    integer, intent(out)         :: nfound
    integer, intent (In)         :: nalloc
    type(DTYPE(kdtree2_result)), target :: results(:)

    !
    sr%qv => qv
    sr%ballsize = r2
    sr%nn = 0      ! flag for fixed ball search
    sr%nfound = 0
    sr%centeridx = -1
    sr%correltime = 0

    sr%results => results

    call DTYPE(validate_query_storage)(nalloc)
    sr%nalloc = nalloc
    sr%overflow = .false. 
    sr%ind => tp%ind
    sr%rearrange= tp%rearrange

    if (tp%rearrange) then
       sr%Data => tp%rearranged_data
    else
       sr%Data => tp%the_data
    endif
    sr%dimen = tp%dimen

    !
    !sr%dsl = Huge(sr%dsl)    ! set to huge positive values
    !sr%il = -1               ! set to invalid indexes
    !

    call DTYPE(search)(tp%root)
    nfound = sr%nfound
    if (tp%sort) then
       call kdtree2_sort_results(nfound, results)
    endif

    if (sr%overflow) then
       write (*,*) 'KD_TREE_TRANS: warning! return from kdtree2_r_nearest found more neighbors'
       write (*,*) 'KD_TREE_TRANS: than storage was provided for.  Answer is NOT smallest ball'
       write (*,*) 'KD_TREE_TRANS: with that number of neighbors!  I.e. it is wrong.'
    endif

    return
  end subroutine DTYPE(kdtree2_r_nearest)

  subroutine DTYPE(kdtree2_r_nearest_around_point)(tp,idxin,correltime,r2,&
   nfound,nalloc,results)
    !
    ! Like kdtree2_r_nearest, but around a point 'idxin' already existing
    ! in the data set. 
    ! 
    ! Results are NOT sorted unless tree was created with sort option.
    !
    type (DTYPE(kdtree2)), pointer      :: tp
    integer, intent (In)         :: idxin, correltime, nalloc
    real(kdkind), intent(in)             :: r2
    integer, intent(out)         :: nfound
    type(DTYPE(kdtree2_result)), target :: results(:)
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic HUGE
    ! ..
    allocate (sr%qv(tp%dimen))
    sr%qv = tp%the_data(:,idxin) ! copy the vector
    sr%ballsize = r2
    sr%nn = 0    ! flag for fixed r search
    sr%nfound = 0
    sr%centeridx = idxin
    sr%correltime = correltime

    sr%results => results

    sr%nalloc = nalloc
    sr%overflow = .false.

    call DTYPE(validate_query_storage)(nalloc)

    !    sr%dsl = HUGE(sr%dsl)    ! set to huge positive values
    !    sr%il = -1               ! set to invalid indexes

    sr%ind => tp%ind
    sr%rearrange = tp%rearrange

    if (tp%rearrange) then
       sr%Data => tp%rearranged_data
    else
       sr%Data => tp%the_data
    endif
    sr%rearrange = tp%rearrange
    sr%dimen = tp%dimen

    !
    !sr%dsl = Huge(sr%dsl)    ! set to huge positive values
    !sr%il = -1               ! set to invalid indexes
    !

    call DTYPE(search)(tp%root)
    nfound = sr%nfound
    if (tp%sort) then
       call kdtree2_sort_results(nfound,results)
    endif

    if (sr%overflow) then
       write (*,*) 'KD_TREE_TRANS: warning! return from kdtree2_r_nearest found more neighbors'
       write (*,*) 'KD_TREE_TRANS: than storage was provided for.  Answer is NOT smallest ball'
       write (*,*) 'KD_TREE_TRANS: with that number of neighbors!  I.e. it is wrong.'
    endif

    deallocate (sr%qv)
    return
  end subroutine DTYPE(kdtree2_r_nearest_around_point)

  function DTYPE(kdtree2_r_count)(tp,qv,r2) result(nfound)
    ! Count the number of neighbors within square distance 'r2'. 
    type (DTYPE(kdtree2)), pointer   :: tp
    real(kdkind), target, intent (In) :: qv(:)
    real(kdkind), intent(in)          :: r2
    integer                   :: nfound
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic HUGE
    ! ..
    sr%qv => qv
    sr%ballsize = r2

    sr%nn = 0       ! flag for fixed r search
    sr%nfound = 0
    sr%centeridx = -1
    sr%correltime = 0
    
    nullify(sr%results) ! for some reason, FTN 95 chokes on '=> null()'

    sr%nalloc = 0            ! we do not allocate any storage but that's OK
                             ! for counting.
    sr%ind => tp%ind
    sr%rearrange = tp%rearrange
    if (tp%rearrange) then
       sr%Data => tp%rearranged_data
    else
       sr%Data => tp%the_data
    endif
    sr%dimen = tp%dimen

    !
    !sr%dsl = Huge(sr%dsl)    ! set to huge positive values
    !sr%il = -1               ! set to invalid indexes
    !
    sr%overflow = .false.

    call DTYPE(search)(tp%root)

    nfound = sr%nfound

    return
  end function DTYPE(kdtree2_r_count)

  function DTYPE(kdtree2_r_count_around_point)(tp,idxin,correltime,r2) &
   result(nfound)
    ! Count the number of neighbors within square distance 'r2' around
    ! point 'idxin' with decorrelation time 'correltime'.
    !
    type (DTYPE(kdtree2)), pointer :: tp
    integer, intent (In)    :: correltime, idxin
    real(kdkind), intent(in)        :: r2
    integer                 :: nfound
    ! ..
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic HUGE
    ! ..
    allocate (sr%qv(tp%dimen))
    sr%qv = tp%the_data(:,idxin)
    sr%ballsize = r2

    sr%nn = 0       ! flag for fixed r search
    sr%nfound = 0
    sr%centeridx = idxin
    sr%correltime = correltime
    nullify(sr%results)

    sr%nalloc = 0            ! we do not allocate any storage but that's OK
                             ! for counting.

    sr%ind => tp%ind
    sr%rearrange = tp%rearrange

    if (sr%rearrange) then
       sr%Data => tp%rearranged_data
    else
       sr%Data => tp%the_data
    endif
    sr%dimen = tp%dimen

    !
    !sr%dsl = Huge(sr%dsl)    ! set to huge positive values
    !sr%il = -1               ! set to invalid indexes
    !
    sr%overflow = .false.

    call DTYPE(search)(tp%root)

    nfound = sr%nfound

    return
  end function DTYPE(kdtree2_r_count_around_point)



  function DTYPE(square_distance)(d, iv,qv) result (res)
    ! distance between iv[1:n] and qv[1:n] 
    ! .. Function Return Value ..
    ! re-implemented to improve vectorization.
    real(kdkind) :: res
    ! ..
    ! ..
    ! .. Scalar Arguments ..
    integer :: d
    ! ..
    ! .. Array Arguments ..
    real(kdkind) :: iv(:),qv(:)
    ! ..
    ! ..
    res = sum( (iv(1:d)-qv(1:d))**2 )
  end function DTYPE(square_distance)
  
  recursive subroutine DTYPE(search)(node)
    !
    ! This is the innermost core routine of the kd-tree search.  Along
    ! with "process_terminal_node", it is the performance bottleneck. 
    !
    ! This version uses a logically complete secondary search of
    ! "box in bounds", whether the sear
    !
    type (DTYPE(Tree_node)), pointer          :: node
    ! ..
    type(DTYPE(tree_node)),pointer            :: ncloser, nfarther
    !
    integer                            :: cut_dim, i
    ! ..
    real(kdkind)                               :: qval, dis
    real(kdkind)                               :: ballsize
    real(kdkind), pointer           :: qv(:)
    type(DTYPE(interval)), pointer :: box(:) 

    if ((associated(node%left) .and. associated(node%right)) .eqv. .false.) then
       ! we are on a terminal node
       if (sr%nn .eq. 0) then
          call DTYPE(process_terminal_node_fixedball)(node)
       else
          call DTYPE(process_terminal_node)(node)
       endif
    else
       ! we are not on a terminal node
       qv => sr%qv(1:)
       cut_dim = node%cut_dim
       qval = qv(cut_dim)

       if (qval < node%cut_val) then
          ncloser => node%left
          nfarther => node%right
          dis = (node%cut_val_right - qval)**2
!          extra = node%cut_val - qval
       else
          ncloser => node%right
          nfarther => node%left
          dis = (node%cut_val_left - qval)**2
!          extra = qval- node%cut_val_left
       endif

       if (associated(ncloser)) call DTYPE(search)(ncloser)

       ! we may need to search the second node. 
       if (associated(nfarther)) then
          ballsize = sr%ballsize
!          dis=extra**2
          if (dis <= ballsize) then
             !
             ! we do this separately as going on the first cut dimen is often
             ! a good idea.
             ! note that if extra**2 < sr%ballsize, then the next
             ! check will also be false. 
             !
             box => node%box(1:)
             do i=1,sr%dimen
                if (i .ne. cut_dim) then
                   dis = dis + DTYPE(dis2_from_bnd)(qv(i),box(i)%lower,box(i)%upper)
                   if (dis > ballsize) then
                      return
                   endif
                endif
             end do
             
             !
             ! if we are still here then we need to search mroe.
             !
             call DTYPE(search)(nfarther)
          endif
       endif
    end if
  end subroutine DTYPE(search)


  real(kdkind) function DTYPE(dis2_from_bnd)(x,amin,amax) result (res)
    real(kdkind), intent(in) :: x, amin,amax

    if (x > amax) then
       res = (x-amax)**2;
       return
    else
       if (x < amin) then
          res = (amin-x)**2;
          return
       else
          res = 0.0
          return
       endif
    endif
    return
  end function DTYPE(dis2_from_bnd)

  logical function DTYPE(box_in_search_range)(node, sr) result(res)
    !
    ! Return the distance from 'qv' to the CLOSEST corner of node's
    ! bounding box
    ! for all coordinates outside the box.   Coordinates inside the box
    ! contribute nothing to the distance.
    !
    type (DTYPE(tree_node)), pointer :: node
    type (DTYPE(tree_search_record)), pointer :: sr

    integer :: dimen, i
    real(kdkind)    :: dis, ballsize
    real(kdkind)    :: l, u

    dimen = sr%dimen
    ballsize = sr%ballsize
    dis = 0.0
    res = .true.
    do i=1,dimen
       l = node%box(i)%lower
       u = node%box(i)%upper
       dis = dis + (DTYPE(dis2_from_bnd)(sr%qv(i),l,u))
       if (dis > ballsize) then
          res = .false.
          return
       endif
    end do
    res = .true.
    return
  end function DTYPE(box_in_search_range)


  subroutine DTYPE(process_terminal_node)(node)
    !
    ! Look for actual near neighbors in 'node', and update
    ! the search results on the sr data structure.
    !
    type (DTYPE(tree_node)), pointer          :: node
    !
    real(kdkind), pointer          :: qv(:)
    integer, pointer       :: ind(:)
    real(kdkind), pointer          :: data(:,:)
    !
    integer                :: dimen, i, indexofi, k, centeridx, correltime
    real(kdkind)                   :: ballsize, sd, newpri
    logical                :: rearrange
    type(DTYPE(pq)), pointer      :: pqp 
    !
    ! copy values from sr to local variables
    !
    !
    ! Notice, making local pointers with an EXPLICIT lower bound
    ! seems to generate faster code.
    ! why?  I don't know.
    qv => sr%qv(1:) 
    pqp => sr%pq
    dimen = sr%dimen
    ballsize = sr%ballsize 
    rearrange = sr%rearrange
    ind => sr%ind(1:)
    data => sr%Data(1:,1:)     
    centeridx = sr%centeridx
    correltime = sr%correltime

    !    doing_correl = (centeridx >= 0)  ! Do we have a decorrelation window? 
    !    include_point = .true.    ! by default include all points
    ! search through terminal bucket.

    mainloop: do i = node%l, node%u
       if (rearrange) then
          sd = 0.0
          do k = 1,dimen
             sd = sd + (data(k,i) - qv(k))**2
             if (sd>ballsize) cycle mainloop
          end do
          indexofi = ind(i)  ! only read it if we have not broken out
       else
          indexofi = ind(i)
          sd = 0.0
          do k = 1,dimen
             sd = sd + (data(k,indexofi) - qv(k))**2
             if (sd>ballsize) cycle mainloop
          end do
       endif

       if (centeridx > 0) then ! doing correlation interval?
          if (abs(indexofi-centeridx) < correltime) cycle mainloop
       endif


       ! 
       ! two choices for any point.  The list so far is either undersized,
       ! or it is not.
       !
       ! If it is undersized, then add the point and its distance
       ! unconditionally.  If the point added fills up the working
       ! list then set the sr%ballsize, maximum distance bound (largest distance on
       ! list) to be that distance, instead of the initialized +infinity. 
       !
       ! If the running list is full size, then compute the
       ! distance but break out immediately if it is larger
       ! than sr%ballsize, "best squared distance" (of the largest element),
       ! as it cannot be a good neighbor. 
       !
       ! Once computed, compare to best_square distance.
       ! if it is smaller, then delete the previous largest
       ! element and add the new one. 

       if (sr%nfound .lt. sr%nn) then
          !
          ! add this point unconditionally to fill list.
          !
          sr%nfound = sr%nfound +1 
          newpri = pq_insert(pqp,sd,indexofi)
          if (sr%nfound .eq. sr%nn) ballsize = newpri
          ! we have just filled the working list.
          ! put the best square distance to the maximum value
          ! on the list, which is extractable from the PQ. 

       else
          !
          ! now, if we get here,
          ! we know that the current node has a squared
          ! distance smaller than the largest one on the list, and
          ! belongs on the list. 
          ! Hence we replace that with the current one.
          !
          ballsize = pq_replace_max(pqp,sd,indexofi)
       endif
    end do mainloop
    !
    ! Reset sr variables which may have changed during loop
    !
    sr%ballsize = ballsize 

  end subroutine DTYPE(process_terminal_node)

  subroutine DTYPE(validate_query_storage)(n)
    !
    ! make sure we have enough storage for n
    !
    integer, intent(in) :: n

    if (size(sr%results,1) .lt. n) then
       write (*,*) 'KD_TREE_TRANS:  not provide enough storage for results(1:n)'
       stop
       return
    endif

    return
  end subroutine DTYPE(validate_query_storage)

  subroutine DTYPE(process_terminal_node_fixedball)(node)
    !
    ! Look for actual near neighbors in 'node', and update
    ! the search results on the sr data structure, i.e.
    ! save all within a fixed ball.
    !
    type (DTYPE(tree_node)), pointer          :: node
    !
    real(kdkind), pointer          :: qv(:)
    integer, pointer       :: ind(:)
    real(kdkind), pointer          :: data(:,:)
    !
    integer                :: nfound
    integer                :: dimen, i, indexofi, k
    integer                :: centeridx, correltime, nn
    real(kdkind)                   :: ballsize, sd
    logical                :: rearrange

    !
    ! copy values from sr to local variables
    !
    qv => sr%qv(1:)
    dimen = sr%dimen
    ballsize = sr%ballsize 
    rearrange = sr%rearrange
    ind => sr%ind(1:)
    data => sr%Data(1:,1:)
    centeridx = sr%centeridx
    correltime = sr%correltime
    nn = sr%nn ! number to search for
    nfound = sr%nfound

    ! search through terminal bucket.
    mainloop: do i = node%l, node%u

       ! 
       ! two choices for any point.  The list so far is either undersized,
       ! or it is not.
       !
       ! If it is undersized, then add the point and its distance
       ! unconditionally.  If the point added fills up the working
       ! list then set the sr%ballsize, maximum distance bound (largest distance on
       ! list) to be that distance, instead of the initialized +infinity. 
       !
       ! If the running list is full size, then compute the
       ! distance but break out immediately if it is larger
       ! than sr%ballsize, "best squared distance" (of the largest element),
       ! as it cannot be a good neighbor. 
       !
       ! Once computed, compare to best_square distance.
       ! if it is smaller, then delete the previous largest
       ! element and add the new one. 

       ! which index to the point do we use? 

       if (rearrange) then
          sd = 0.0
          do k = 1,dimen
             sd = sd + (data(k,i) - qv(k))**2
             if (sd>ballsize) cycle mainloop
          end do
          indexofi = ind(i)  ! only read it if we have not broken out
       else
          indexofi = ind(i)
          sd = 0.0
          do k = 1,dimen
             sd = sd + (data(k,indexofi) - qv(k))**2
             if (sd>ballsize) cycle mainloop
          end do
       endif

       if (centeridx > 0) then ! doing correlation interval?
          if (abs(indexofi-centeridx)<correltime) cycle mainloop
       endif

       nfound = nfound+1
       if (nfound .gt. sr%nalloc) then
          ! oh nuts, we have to add another one to the tree but
          ! there isn't enough room.
          sr%overflow = .true.
       else
          sr%results(nfound)%dis = sd
          sr%results(nfound)%idx = indexofi
       endif
    end do mainloop
    !
    ! Reset sr variables which may have changed during loop
    !
    sr%nfound = nfound
  end subroutine DTYPE(process_terminal_node_fixedball)

  subroutine DTYPE(kdtree2_n_nearest_brute_force)(tp,qv,nn,results) 
    ! find the 'n' nearest neighbors to 'qv' by exhaustive search.
    ! only use this subroutine for testing, as it is SLOW!  The
    ! whole point of a k-d tree is to avoid doing what this subroutine
    ! does.
    type (DTYPE(kdtree2)), pointer :: tp
    real(kdkind), intent (In)       :: qv(:)
    integer, intent (In)    :: nn
    type(DTYPE(kdtree2_result))    :: results(:) 

    integer :: i, j, k
    real(kdkind), allocatable :: all_distances(:)
    ! ..
    allocate (all_distances(tp%n))
    do i = 1, tp%n
       all_distances(i) = DTYPE(square_distance)(tp%dimen,qv,tp%the_data(:,i))
    end do
    ! now find 'n' smallest distances
    do i = 1, nn
       results(i)%dis =  huge(1.0)
       results(i)%idx = -1
    end do
    do i = 1, tp%n
       if (all_distances(i)<results(nn)%dis) then
          ! insert it somewhere on the list
          do j = 1, nn
             if (all_distances(i)<results(j)%dis) exit
          end do
          ! now we know 'j'
          do k = nn - 1, j, -1
             results(k+1) = results(k)
          end do
          results(j)%dis = all_distances(i)
          results(j)%idx = i
       end if
    end do
    deallocate (all_distances)
  end subroutine DTYPE(kdtree2_n_nearest_brute_force)
  

  subroutine DTYPE(kdtree2_r_nearest_brute_force)(tp,qv,r2,nfound,results) 
    ! find the nearest neighbors to 'qv' with distance**2 <= r2 by exhaustive search.
    ! only use this subroutine for testing, as it is SLOW!  The
    ! whole point of a k-d tree is to avoid doing what this subroutine
    ! does.
    type (DTYPE(kdtree2)), pointer :: tp
    real(kdkind), intent (In)       :: qv(:)
    real(kdkind), intent (In)       :: r2
    integer, intent(out)    :: nfound
    type(DTYPE(kdtree2_result))    :: results(:) 

    integer :: i, nalloc
    real(kdkind), allocatable :: all_distances(:)
    ! ..
    allocate (all_distances(tp%n))
    do i = 1, tp%n
       all_distances(i) = DTYPE(square_distance)(tp%dimen,qv,tp%the_data(:,i))
    end do
    
    nfound = 0
    nalloc = size(results,1)

    do i = 1, tp%n
       if (all_distances(i)< r2) then
          ! insert it somewhere on the list
          if (nfound .lt. nalloc) then
             nfound = nfound+1
             results(nfound)%dis = all_distances(i)
             results(nfound)%idx = i
          endif
       end if
    enddo
    deallocate (all_distances)

    call kdtree2_sort_results(nfound,results)


  end subroutine DTYPE(kdtree2_r_nearest_brute_force)

  subroutine DTYPE(kdtree2_sort_results)(nfound,results)
    !  Use after search to sort results(1:nfound) in order of increasing 
    !  distance.
    integer, intent(in)          :: nfound
    type(DTYPE(kdtree2_result)), target :: results(:) 
    !
    !

    !THIS IS BUGGY WITH INTEL FORTRAN
    !    If (nfound .Gt. 1) Call heapsort(results(1:nfound)%dis,results(1:nfound)%ind,nfound)
    !
    if (nfound .gt. 1) call DTYPE(heapsort_struct)(results,nfound)

    return
  end subroutine DTYPE(kdtree2_sort_results)

  subroutine DTYPE(heapsort)(a,ind,n)
    !
    ! Sort a(1:n) in ascending order, permuting ind(1:n) similarly.
    ! 
    ! If ind(k) = k upon input, then it will give a sort index upon output.
    !
    integer,intent(in)          :: n
    real(kdkind), intent(inout)         :: a(:) 
    integer, intent(inout)      :: ind(:)

    !
    !
    real(kdkind)        :: value   ! temporary for a value from a()
    integer     :: ivalue  ! temporary for a value from ind()

    integer     :: i,j
    integer     :: ileft,iright

    ileft=n/2+1
    iright=n

    !    do i=1,n
    !       ind(i)=i
    ! Generate initial idum array
    !    end do

    if(n.eq.1) return                  

    do 
       if(ileft > 1)then
          ileft=ileft-1
          value=a(ileft); ivalue=ind(ileft)
       else
          value=a(iright); ivalue=ind(iright)
          a(iright)=a(1); ind(iright)=ind(1)
          iright=iright-1
          if (iright == 1) then
             a(1)=value;ind(1)=ivalue
             return
          endif
       endif
       i=ileft
       j=2*ileft
       do while (j <= iright) 
          if(j < iright) then
             if(a(j) < a(j+1)) j=j+1
          endif
          if(value < a(j)) then
             a(i)=a(j); ind(i)=ind(j)
             i=j
             j=j+j
          else
             j=iright+1
          endif
       end do
       a(i)=value; ind(i)=ivalue
    end do
  end subroutine DTYPE(heapsort)

  subroutine DTYPE(heapsort_struct)(a,n)
    !
    ! Sort a(1:n) in ascending order
    ! 
    !
    integer,intent(in)                 :: n
    type(DTYPE(kdtree2_result)),intent(inout) :: a(:)

    !
    !
    type(DTYPE(kdtree2_result)) :: value ! temporary value

    integer     :: i,j
    integer     :: ileft,iright

    ileft=n/2+1
    iright=n

    !    do i=1,n
    !       ind(i)=i
    ! Generate initial idum array
    !    end do

    if(n.eq.1) return                  

    do 
       if(ileft > 1)then
          ileft=ileft-1
          value=a(ileft)
       else
          value=a(iright)
          a(iright)=a(1)
          iright=iright-1
          if (iright == 1) then
             a(1) = value
             return
          endif
       endif
       i=ileft
       j=2*ileft
       do while (j <= iright) 
          if(j < iright) then
             if(a(j)%dis < a(j+1)%dis) j=j+1
          endif
          if(value%dis < a(j)%dis) then
             a(i)=a(j); 
             i=j
             j=j+j
          else
             j=iright+1
          endif
       end do
       a(i)=value
    end do
  end subroutine DTYPE(heapsort_struct)
