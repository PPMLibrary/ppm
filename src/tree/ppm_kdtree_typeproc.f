
      SUBROUTINE DTYPE(pq_create)(this,results_in)
        !!!
        !!! Create a priority queue from ALREADY allocated
        !!! array pointers for storage.  NOTE! It will NOT
        !!! add any alements to the heap, i.e. any existing
        !!! data in the input arrays will NOT be used and may
        !!! be overwritten.
        !!!
        !!! usage:
        !!!    CALL  pq%pq_create(x,k)

        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(pq)                                         :: this

        TYPE(DTYPE(kdtree2_result)), DIMENSION(:), TARGET :: results_in
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        INTEGER :: nalloc,info

        CHARACTER(LEN=ppm_char) :: caller='pq_create'

        nalloc = SIZE(results_in,1)
        IF (nalloc.LT.1) THEN
           fail("PQ_CREATE: error, input arrays must be allocated.", &
           &   ppm_err_wrong_dim,exit_point=no,ppm_error=ppm_error_fatal)
        ENDIF

        this%DTYPE(elems) => results_in
        this%heap_size = 0

      END SUBROUTINE DTYPE(pq_create)

      FUNCTION DTYPE(pq_insert)(this,dis_in,idx) RESULT(dis_out)
        !
        ! Insert a new element and return the new maximum priority,
        ! which may or may not be the same as the old maximum priority.
        !
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(pq)               :: this

        REAL(MK), INTENT(IN   ) :: dis_in

        INTEGER,  INTENT(IN   ) :: idx

        REAL(MK)                :: dis_out
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        REAL(MK) :: parentdis

        INTEGER :: i, parent

        this%heap_size = this%heap_size + 1
        i = this%heap_size

        DO WHILE (i.GT.1)
           parent = INT(i/2)
           parentdis = this%DTYPE(elems)(parent)%dis
           IF (dis_in.GT.parentdis) THEN
              ! move what was in i's parent into i.
              this%DTYPE(elems)(i)=DTYPE(kdtree2_result)(parentdis,this%DTYPE(elems)(parent)%idx)
              i = parent
           ELSE
              EXIT
           ENDIF
        ENDDO

        ! insert the element at the determined position
        this%DTYPE(elems)(i) = DTYPE(kdtree2_result)(dis_in,idx)

        dis_out = this%DTYPE(elems)(1)%dis
        RETURN
      END FUNCTION DTYPE(pq_insert)


      SUBROUTINE DTYPE(heapify)(a,i_in)
        !
        ! take a heap rooted at 'i' and force it to be in the
        ! heap canonical form.   This is performance critical
        ! and has been tweaked a little to reflect this.
        !
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(pq), INTENT(INOUT) :: a

        INTEGER,   INTENT(IN   ) :: i_in
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        TYPE(DTYPE(kdtree2_result)) :: temp

        REAL(MK) :: pri_i, pri_l, pri_r, pri_largest

        INTEGER :: i, l, r, largest

        i = i_in

        bigloop:  DO
           l = 2*i ! left(i)
           r = l+1 ! right(i)
           !
           ! set 'largest' to the index of either i, l, r
           ! depending on whose priority is largest.
           !
           ! note that l or r can be larger than the heap size
           ! in which case they do not count.

           ! does left child have higher priority?
           IF (l.GT.a%heap_size) THEN
              ! we know that i is the largest as both l and r are invalid.
              EXIT bigloop
           ELSE
              pri_i = a%DTYPE(elems)(i)%dis
              pri_l = a%DTYPE(elems)(l)%dis
              IF (pri_l.GT.pri_i) THEN
                 largest = l
                 pri_largest = pri_l
              ELSE
                 largest = i
                 pri_largest = pri_i
              ENDIF

              !
              ! between i and l we have a winner
              ! now choose between that and r.
              !
              IF (r.LE.a%heap_size) THEN
                 pri_r = a%DTYPE(elems)(r)%dis
                 IF (pri_r.GT.pri_largest) THEN
                    largest = r
                 ENDIF
              ENDIF
           ENDIF

           IF (largest.NE.i) THEN
              ! swap data in nodes largest and i, THEN heapify

              temp = a%DTYPE(elems)(i)
              a%DTYPE(elems)(i) = a%DTYPE(elems)(largest)
              a%DTYPE(elems)(largest) = temp
              !
              ! Canonical heapify() algorithm has tail-ecursive call:
              !
              !        call heapify(a,largest)
              ! we will simulate with cycle
              !
              i = largest
              CYCLE bigloop ! continue the loop
           ELSE
              RETURN   ! break from the loop
           ENDIF
        ENDDO bigloop
        RETURN
      END SUBROUTINE DTYPE(heapify)

#if  __KIND == __SINGLE_PRECISION
      SUBROUTINE pq_delete(this,i,info)
        !
        ! delete item with index 'i'
        !
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(pq)                :: this

        INTEGER,   INTENT(IN   ) :: i

        INTEGER,   INTENT(  OUT) :: info
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------

        start_subroutine("pq_delete")

        IF (i.LT.1.OR.i.GT.this%heap_size) THEN
           fail("PQ_DELETE: error, attempt to remove out of bounds element.", &
           & ppm_error=ppm_error_fatal)
        ENDIF

        IF (ASSOCIATED(this%elems_s)) THEN
           ! swap the item to be deleted with the last element
           ! and shorten heap by one.
           this%elems_s(i) = this%elems_s(this%heap_size)
           this%heap_size  = this%heap_size - 1

           CALL heapify_s(this,i)
        ELSE
           ! swap the item to be deleted with the last element
           ! and shorten heap by one.
           this%elems_d(i) = this%elems_d(this%heap_size)
           this%heap_size  = this%heap_size - 1

           CALL heapify_d(this,i)
        ENDIF

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        end_subroutine()
        RETURN
      END SUBROUTINE pq_delete
#endif

      SUBROUTINE DTYPE(pq_extract_max)(this,e,info)
        !
        ! return the priority and payload of maximum priority
        ! element, and remove it from the queue.
        ! (equivalent to 'pop()' on a stack)
        !
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(pq)                                  :: this

        TYPE(DTYPE(kdtree2_result)), INTENT(  OUT) :: e

        INTEGER,                     INTENT(  OUT) :: info
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------

        start_subroutine("pq_extract_max")

        IF (this%heap_size.GE.1) THEN
           !
           ! return max as first element
           !
           e = this%DTYPE(elems)(1)

           !
           ! move last element to first
           !
           this%DTYPE(elems)(1) = this%DTYPE(elems)(this%heap_size)
           this%heap_size = this%heap_size-1

           CALL DTYPE(heapify)(this,1)

           GOTO 9999
         ELSE
           fail('PQ_EXTRACT_MAX: error, attempted to pop non-positive PQ', &
           &    ppm_error=ppm_error_fatal)
         ENDIF

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        end_subroutine()
        RETURN
      END SUBROUTINE DTYPE(pq_extract_max)

      SUBROUTINE DTYPE(pq_max)(this,e,info)
        !
        ! return the priority and its payload of the maximum priority element
        ! on the queue, which should be the first one, if it is
        ! in heapified form.
        !
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(pq)                                  :: this

        TYPE(DTYPE(kdtree2_result)), INTENT(  OUT) :: e

        INTEGER,                     INTENT(  OUT) :: info
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------

        start_subroutine("pq_max")

        IF (this%heap_size.GT.0) THEN
           !
           ! return max as first element
           !
           e = this%DTYPE(elems)(1)
           GOTO 9999
         ELSE
           fail('PQ_MAX: ERROR, heap_size < 1',ppm_error=ppm_error_fatal)
         ENDIF

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        end_subroutine()
        RETURN
      END SUBROUTINE DTYPE(pq_max)

      FUNCTION DTYPE(pq_replace_max)(this,dis,idx) RESULT(pqrmax)
        !
        ! Replace the extant maximum priority element
        ! in the PQ with (dis,idx).
        ! Return the new maximum priority, which may be larger
        ! or smaller than the old one.
        !
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(pq)               :: this

        REAL(MK), INTENT(IN   ) :: dis

        INTEGER,  INTENT(IN   ) :: idx

        REAL(MK)                :: pqrmax

        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        REAL(MK) :: prichild,prichildp1

        INTEGER :: parent,child,N

        N=this%heap_size

        IF (N.GE.1) THEN
           parent =1
           child=2

           loop: DO WHILE (child.LE.N)
              prichild = this%DTYPE(elems)(child)%dis
              !
              ! posibly child+1 has higher priority, and if
              ! so, get it, and increment child.
              !

              IF (child.LT.N) THEN
                 prichildp1 = this%DTYPE(elems)(child+1)%dis
                 IF (prichild.LT.prichildp1) THEN
                    child = child+1
                    prichild = prichildp1
                 ENDIF
              ENDIF
              IF (dis.GE.prichild) THEN
                 EXIT loop
                 ! we have a proper place for our new element,
                 ! bigger than either children's priority.
              ELSE
                 ! move child into parent.
                 this%DTYPE(elems)(parent) = this%DTYPE(elems)(child)
                 parent = child
                 child = 2*parent
              ENDIF
           ENDDO loop

           this%DTYPE(elems)(parent) = DTYPE(kdtree2_result)(dis,idx)

           pqrmax = this%DTYPE(elems)(1)%dis
        ELSE
           this%DTYPE(elems)(1) = DTYPE(kdtree2_result)(dis,idx)

           pqrmax = dis
        ENDIF

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END FUNCTION DTYPE(pq_replace_max)

      SUBROUTINE DTYPE(kdtree2_create)(this,input_data,info,dim,sort,rearrange)
        !
        ! create the actual tree structure, given an input array of data.
        !
        ! Note, input data is input_data(1:d,1:N).
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
        !
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(DTYPE(kdtree2))                   :: this

        REAL(MK), DIMENSION(:,:), TARGET        :: input_data

        INTEGER,                  INTENT(  OUT) :: info
        ! .. Array Arguments ..
        INTEGER,  OPTIONAL,       INTENT(IN   ) :: dim

        LOGICAL,  OPTIONAL,       INTENT(IN   ) :: sort
        LOGICAL,  OPTIONAL,       INTENT(IN   ) :: rearrange

        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------

        INTEGER :: i,j,ldu(2),iopt

        start_subroutine("kdtree2_create")

        this%the_data => input_data
        ! pointer assignment

#ifdef __DEBUG
        check_true(<#ppm_dim.EQ.SIZE(input_data,1)#>,"ppm_dim is different from input data size")
#endif
        this%dimen=MERGE(dim,ppm_dim,PRESENT(dim))
        this%n    =SIZE(input_data,2)

        CALL this%build(info)
        or_fail("building the tree has failed.")

        this%sort     =MERGE(sort,     .FALSE.,PRESENT(sort))
        this%rearrange=MERGE(rearrange,.TRUE., PRESENT(rearrange))

        IF (this%rearrange) THEN
           iopt=ppm_param_alloc_fit
           ldu(1)=this%dimen
           ldu(2)=this%n
           CALL ppm_alloc(this%rearranged_data,ldu,iopt,info)
           or_fail_alloc("rearranged_data allocation has failed.", &
           & exit_point=no,ppm_error=ppm_error_fatal)

           FORALL (i=1:this%dimen,j=1:this%n) this%rearranged_data(i,j) = this%the_data(i,this%ind(j))
        ENDIF

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        end_subroutine()
        RETURN
      END SUBROUTINE DTYPE(kdtree2_create)

      SUBROUTINE DTYPE(kdtree2_destroy)(this,info)
        !
        ! Deallocates all memory for the tree, except input data matrix
        !
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(DTYPE(kdtree2))  :: this

        INTEGER, INTENT(  OUT) :: info

        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------

        INTEGER :: i,iopt,ldu(2)

        start_subroutine("kdtree2_destroy")

        CALL this%destroy_node(this%root)

        iopt=ppm_param_dealloc
        ldu=0
        CALL ppm_alloc(this%ind,ldu,iopt,info)
        or_fail_dealloc("this%ind",ppm_error=ppm_error_fatal)

        IF (this%rearrange) THEN
           CALL ppm_alloc(this%rearranged_data,ldu,iopt,info)
           or_fail_dealloc("this%rearranged_data",ppm_error=ppm_error_fatal)
        ENDIF

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        end_subroutine()
        RETURN
      END SUBROUTINE DTYPE(kdtree2_destroy)

      RECURSIVE SUBROUTINE DTYPE(kdtree2_destroy_node)(this,tn)
        !
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(DTYPE(kdtree2))           :: this

        TYPE(DTYPE(tree_node)), POINTER :: tn

        IF (ASSOCIATED(tn%left)) THEN
           CALL this%destroy_node(tn%left)
           NULLIFY (tn%left)
        ENDIF
        IF (ASSOCIATED(tn%right)) THEN
           CALL this%destroy_node(tn%right)
           NULLIFY (tn%right)
        ENDIF
        IF (ASSOCIATED(tn%box)) DEALLOCATE(tn%box)
        DEALLOCATE(tn)
        NULLIFY(tn)
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        RETURN
      END SUBROUTINE DTYPE(kdtree2_destroy_node)

      SUBROUTINE DTYPE(build_tree)(this,info)
        !
        !
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(DTYPE(kdtree2))  :: this

        INTEGER, INTENT(  OUT) :: info
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        TYPE(DTYPE(tree_node)), POINTER :: dummy

        INTEGER :: i,iopt,ldu(1)

        start_subroutine("build_tree")

        iopt=ppm_param_alloc_fit
        ldu=this%n
        CALL ppm_alloc(this%ind,ldu,iopt,info)
        or_fail_alloc("kdtree2 ind allocation failed.",ppm_error=ppm_error_fatal)

        FORALL (i=1:this%n) this%ind(i)=i

        NULLIFY(dummy)
        CALL this%build_range(1,this%n,dummy,this%root)

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        end_subroutine()
        RETURN
      END SUBROUTINE DTYPE(build_tree)

      RECURSIVE SUBROUTINE DTYPE(build_tree_for_range)(this,l,u,parent,res)
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(DTYPE(kdtree2))                 :: this
        ! .. Structure Arguments ..

        INTEGER,                INTENT(IN   ) :: l
        INTEGER,                INTENT(IN   ) :: u

        TYPE(DTYPE(tree_node)), POINTER       :: parent
        TYPE(DTYPE(tree_node)), POINTER       :: res
        ! .. Function Return Cut_value ..
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        REAL(MK) :: average

        INTEGER :: i,c,m,dimen
        INTEGER :: info

        CHARACTER(LEN=ppm_char) :: caller="build_tree_for_range"

        LOGICAL :: recompute

        ! First, compute an APPROXIMATE bounding box of all points associated with this node.
        IF (u.LT.l) THEN
           ! no points in this box
           NULLIFY(res)
           RETURN
        ENDIF

        ! first compute min and max
        dimen = this%dimen

        ALLOCATE(res,STAT=info)
        or_fail_alloc("Allocation of the result pointer has failed.", &
        & exit_point=no,ppm_error=ppm_error_fatal)

        ALLOCATE(res%box(dimen),STAT=info)
        or_fail_alloc("res%box allocation has failed.", &
        & exit_point=no,ppm_error=ppm_error_fatal)

        IF ((u-l).LE.bucket_size) THEN
           !
           ! always compute true bounding box for terminal nodes.
           !
           DO i=1,dimen
              CALL this%spread_in_coordinate(i,l,u,res%box(i))
           ENDDO
           res%cut_dim = 0
           res%cut_val = 0.0_MK
           res%l = l
           res%u = u
        ELSE
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
           DO i=1,dimen
              recompute=.TRUE.
              IF (ASSOCIATED(parent)) THEN
                 IF (i.NE.parent%cut_dim) THEN
                    recompute=.FALSE.
                 ENDIF
              ENDIF
              IF (recompute) THEN
                 CALL this%spread_in_coordinate(i,l,u,res%box(i))
              ELSE
                 res%box(i) = parent%box(i)
              ENDIF
           ENDDO


           c = MAXLOC(res%box(1:dimen)%upper-res%box(1:dimen)%lower,1)
           !
           ! c is the identity of which coordinate has the greatest spread.
           !

           !
           ! select point halfway between min and max, as per A.
           ! Moore, who says this helps in some degenerate cases,
           ! or actual arithmetic average.
           !
           ! actually compute average
           average = SUM(this%the_data(c,this%ind(l:u)))/REAL(u-l+1,MK)

           ! select exact median to have fully balanced tree.
           ! m = (l+u)/2
           ! CALL this%select_on_coordinate(c,m,l,u)
           ! or
           ! average = (res%box(c)%upper + res%box(c)%lower)/2.0_MK

           res%cut_val = average
           m = this%select_on_coordinate_value(c,average,l,u)

           ! moves indexes around
           res%cut_dim = c
           res%l = l
           res%u = u

           CALL this%build_range(l  ,m,res,res%left)
           CALL this%build_range(m+1,u,res,res%right)

           IF (.NOT.ASSOCIATED(res%right)) THEN
              res%box          = res%left%box
              res%cut_val_left = res%left%box(c)%upper
              res%cut_val      = res%cut_val_left
           ELSE IF (.NOT.ASSOCIATED(res%left)) THEN
              res%box           = res%right%box
              res%cut_val_right = res%right%box(c)%lower
              res%cut_val       = res%cut_val_right
           ELSE
              res%cut_val_right = res%right%box(c)%lower
              res%cut_val_left  = res%left%box(c)%upper
              res%cut_val       = (res%cut_val_left + res%cut_val_right)/2.0_MK

              ! now remake the true bounding box for self.
              ! Since we are taking unions (in effect) of a tree structure,
              ! this is much faster than doing an exhaustive
              ! search over all points
              res%box%upper = MAX(res%left%box%upper,res%right%box%upper)
              res%box%lower = MIN(res%left%box%lower,res%right%box%lower)
           ENDIF
        ENDIF !(u-l).LE.bucket_size

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
        RETURN
      END SUBROUTINE DTYPE(build_tree_for_range)

      SUBROUTINE DTYPE(spread_in_coordinate)(this,c,l,u,interv)
        !
        ! the spread in coordinate 'c', between l and u.
        !
        ! Return lower bound in 'smin', and upper in 'smax',
        !
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(DTYPE(kdtree2))                :: this
        ! .. Structure Arguments ..
        INTEGER,               INTENT(IN   ) :: c
        INTEGER,               INTENT(IN   ) :: l
        INTEGER,               INTENT(IN   ) :: u
        ! .. Scalar Arguments ..
        TYPE(DTYPE(interval)), INTENT(  OUT) :: interv
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        REAL(MK), DIMENSION(:,:), POINTER :: tdata
        REAL(MK)                          :: last,lmax,lmin,tmp,smin,smax

        INTEGER, DIMENSION(:), POINTER :: ind
        INTEGER                        :: i,ulocal

        tdata => this%the_data
        ind   => this%ind

        smin = tdata(c,ind(l))
        smax = smin

        ulocal = u

        DO i=l+2,ulocal,2
           lmin = tdata(c,ind(i-1))
           lmax = tdata(c,ind(i))
           IF (lmin.GT.lmax) THEN
              tmp  = lmin
              lmin = lmax
              lmax = tmp
           ENDIF
           IF (smin.GT.lmin) smin = lmin
           IF (smax.LT.lmax) smax = lmax
        ENDDO
        IF (i.EQ.ulocal+1) THEN
           last = tdata(c,ind(ulocal))
           IF (smin.GT.last) smin = last
           IF (smax.LT.last) smax = last
        ENDIF

        interv = DTYPE(interval)(smin,smax)

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        RETURN
      END SUBROUTINE DTYPE(spread_in_coordinate)


      FUNCTION DTYPE(select_on_coordinate_value)(this,c,alpha,li,ui) RESULT(res)
        !
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
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(DTYPE(kdtree2))   :: this

        INTEGER,  INTENT(IN   ) :: c

        REAL(MK), INTENT(IN   ) :: alpha

        INTEGER,  INTENT(IN   ) :: li
        INTEGER,  INTENT(IN   ) :: ui

        INTEGER                 :: res

        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        REAL(MK), DIMENSION(:,:), POINTER :: tdata

        INTEGER, DIMENSION(:), POINTER :: ind
        INTEGER                        :: tmp
        INTEGER                        :: lb, rb

        tdata => this%the_data
        ind   => this%ind

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
        lb = li
        rb = ui

        DO WHILE (lb.LT.rb)
           IF (tdata(c,ind(lb)).LE.alpha) THEN
              ! it is good where it is.
              lb=lb+1
           ELSE
              ! swap it with rb.
              tmp    =ind(lb)
              ind(lb)=ind(rb)
              ind(rb)=tmp
              rb     =rb-1
           ENDIF
        ENDDO
        ! now lb .eq. ub
        IF (tdata(c,ind(lb)).LE.alpha) THEN
           res = lb
        ELSE
           res = lb-1
        ENDIF

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        RETURN
      END FUNCTION DTYPE(select_on_coordinate_value)

#undef  __KIND
#undef    MK
#undef   _MK
#undef    DTYPE
