    function DTYPE(pq_create)(results_in) result(res)
        !
        ! Create a priority queue from ALREADY allocated
        ! array pointers for storage.  NOTE! It will NOT
        ! add any alements to the heap, i.e. any existing
        ! data in the input arrays will NOT be used and may
        ! be overwritten.
        ! 
        ! usage:
        !    real(kdkind), pointer :: x(:)
        !    integer, pointer :: k(:)
        !    allocate(x(1000),k(1000))
        !    pq => pq_create(x,k)
        !
        type(DTYPE(kdtree2_result)), target:: results_in(:) 
        type(DTYPE(pq)) :: res
        !
        !
        integer :: nalloc

        nalloc = size(results_in,1)
        if (nalloc .lt. 1) then
            write (*,*) 'PQ_CREATE: error, input arrays must be allocated.'
        end if
        res%elems => results_in
        res%heap_size = 0
        return
    end function DTYPE(pq_create)

    !
    ! operations for getting parents and left + right children
    ! of elements in a binary heap.
    !


    subroutine DTYPE(heapify)(a,i_in)
        !
        ! take a heap rooted at 'i' and force it to be in the
        ! heap canonical form.   This is performance critical 
        ! and has been tweaked a little to reflect this.
        !
        type(DTYPE(pq)),pointer   :: a
        integer, intent(in) :: i_in
        !
        integer :: i, l, r, largest

        real(kdkind)    :: pri_i, pri_l, pri_r, pri_largest


        type(DTYPE(kdtree2_result)) :: temp

        i = i_in

        bigloop:  do
            l = 2*i ! left(i)
            r = l+1 ! right(i)
            ! 
            ! set 'largest' to the index of either i, l, r
            ! depending on whose priority is largest.
            !
            ! note that l or r can be larger than the heap size
            ! in which case they do not count.


            ! does left child have higher priority? 
            if (l .gt. a%heap_size) then
                ! we know that i is the largest as both l and r are invalid.
                exit 
            else
                pri_i = a%elems(i)%dis
                pri_l = a%elems(l)%dis 
                if (pri_l .gt. pri_i) then
                    largest = l
                    pri_largest = pri_l
                else
                    largest = i
                    pri_largest = pri_i
                endif

                !
                ! between i and l we have a winner
                ! now choose between that and r.
                !
                if (r .le. a%heap_size) then
                    pri_r = a%elems(r)%dis
                    if (pri_r .gt. pri_largest) then
                        largest = r
                    endif
                endif
            endif

            if (largest .ne. i) then
                ! swap data in nodes largest and i, then heapify

                temp = a%elems(i)
                a%elems(i) = a%elems(largest)
                a%elems(largest) = temp 
                ! 
                ! Canonical heapify() algorithm has tail-ecursive call: 
                !
                !        call heapify(a,largest)   
                ! we will simulate with cycle
                !
                i = largest
                cycle bigloop ! continue the loop 
            else
                return   ! break from the loop
            end if
        enddo bigloop
        return
    end subroutine DTYPE(heapify)

    subroutine DTYPE(pq_max)(a,e) 
        !
        ! return the priority and its payload of the maximum priority element
        ! on the queue, which should be the first one, if it is 
        ! in heapified form.
        !
        type(DTYPE(pq)),pointer :: a
        type(DTYPE(kdtree2_result)),intent(out)  :: e

        if (a%heap_size .gt. 0) then
            e = a%elems(1) 
        else
            write (*,*) 'PQ_MAX: ERROR, heap_size < 1'
            stop
        endif
        return
    end subroutine DTYPE(pq_max)

    real(kdkind) function DTYPE(pq_maxpri)(a)
        type(DTYPE(pq)), pointer :: a

        if (a%heap_size .gt. 0) then
            DTYPE(pq_maxpri) = a%elems(1)%dis
        else
            write (*,*) 'PQ_MAX_PRI: ERROR, heapsize < 1'
            stop
        endif
        return
    end function DTYPE(pq_maxpri)

    subroutine DTYPE(pq_extract_max)(a,e)
        !
        ! return the priority and payload of maximum priority
        ! element, and remove it from the queue.
        ! (equivalent to 'pop()' on a stack)
        !
        type(DTYPE(pq)),pointer :: a
        type(DTYPE(kdtree2_result)), intent(out) :: e

        if (a%heap_size .ge. 1) then
            !
            ! return max as first element
            !
            e = a%elems(1) 

            !
            ! move last element to first
            !
            a%elems(1) = a%elems(a%heap_size) 
            a%heap_size = a%heap_size-1
            call DTYPE(heapify)(a,1)
            return
        else
            write (*,*) 'PQ_EXTRACT_MAX: error, attempted to pop non-positive PQ'
            stop
        end if

    end subroutine DTYPE(pq_extract_max)


    real(kdkind) function DTYPE(pq_insert)(a,dis,idx) 
        !
        ! Insert a new element and return the new maximum priority,
        ! which may or may not be the same as the old maximum priority.
        !
        type(DTYPE(pq)),pointer  :: a
        real(kdkind), intent(in) :: dis
        integer, intent(in) :: idx
        !    type(kdtree2_result), intent(in) :: e
        !
        integer :: i, isparent
        real(kdkind)    :: parentdis
        !

        a%heap_size = a%heap_size + 1
        i = a%heap_size

        do while (i .gt. 1)
            isparent = int(i/2)
            parentdis = a%elems(isparent)%dis
            if (dis .gt. parentdis) then
                ! move what was in i's parent into i.
                a%elems(i)%dis = parentdis
                a%elems(i)%idx = a%elems(isparent)%idx
                i = isparent
            else
                exit
            endif
        end do

        ! insert the element at the determined position
        a%elems(i)%dis = dis
        a%elems(i)%idx = idx

        DTYPE(pq_insert) = a%elems(1)%dis 
        return
        !    end if

    end function DTYPE(pq_insert)

    subroutine DTYPE(pq_adjust_heap)(a,i)
        type(DTYPE(pq)),pointer  :: a
        integer, intent(in) :: i
        !
        ! nominally arguments (a,i), but specialize for a=1
        !
        ! This routine assumes that the trees with roots 2 and 3 are already heaps,
        ! i.e. the children of '1' are heaps.  When the procedure is completed, the
        ! tree rooted at 1 is a heap.
        real(kdkind) :: prichild
        integer :: parent, child, N

        type(DTYPE(kdtree2_result)) :: e

        e = a%elems(i) 

        parent = i
        child = 2*i
        N = a%heap_size

        do while (child .le. N)
            if (child .lt. N) then
                if (a%elems(child)%dis .lt. a%elems(child+1)%dis) then
                    child = child+1
                endif
            endif
            prichild = a%elems(child)%dis
            if (e%dis .ge. prichild) then
                exit 
            else
                ! move child into parent.
                a%elems(parent) = a%elems(child) 
                parent = child
                child = 2*parent
            end if
        end do
        a%elems(parent) = e
        return
    end subroutine DTYPE(pq_adjust_heap)


    real(kdkind) function DTYPE(pq_replace_max)(a,dis,idx) 
        !
        ! Replace the extant maximum priority element
        ! in the PQ with (dis,idx).  Return
        ! the new maximum priority, which may be larger
        ! or smaller than the old one.
        !
        type(DTYPE(pq)),pointer         :: a
        real(kdkind), intent(in) :: dis
        integer, intent(in) :: idx
        !    type(kdtree2_result), intent(in) :: e
        ! not tested as well!

        integer :: parent, child, N
        real(kdkind)    :: prichild, prichildp1

        type(DTYPE(kdtree2_result)) :: etmp

        if (.true.) then
            N=a%heap_size
            if (N .ge. 1) then
                parent =1
                child=2

                loop: do while (child .le. N)
                    prichild = a%elems(child)%dis

                    !
                    ! posibly child+1 has higher priority, and if
                    ! so, get it, and increment child.
                    !

                    if (child .lt. N) then
                        prichildp1 = a%elems(child+1)%dis
                        if (prichild .lt. prichildp1) then
                            child = child+1
                            prichild = prichildp1
                        endif
                    endif

                    if (dis .ge. prichild) then
                        exit loop  
                        ! we have a proper place for our new element, 
                        ! bigger than either children's priority.
                    else
                        ! move child into parent.
                        a%elems(parent) = a%elems(child) 
                        parent = child
                        child = 2*parent
                    end if
                end do loop
                a%elems(parent)%dis = dis
                a%elems(parent)%idx = idx
                DTYPE(pq_replace_max) = a%elems(1)%dis
            else
                a%elems(1)%dis = dis
                a%elems(1)%idx = idx
                DTYPE(pq_replace_max) = dis
            endif
        else
            !
            ! slower version using elementary pop and push operations.
            !
            call pq_extract_max(a,etmp) 
            etmp%dis = dis
            etmp%idx = idx
            DTYPE(pq_replace_max) = pq_insert(a,dis,idx)
        endif
        return
    end function DTYPE(pq_replace_max)

    subroutine DTYPE(pq_delete)(a,i)
        ! 
        ! delete item with index 'i'
        !
        type(DTYPE(pq)),pointer :: a
        integer           :: i

        if ((i .lt. 1) .or. (i .gt. a%heap_size)) then
            write (*,*) &
                'PQ_DELETE: error, attempt to remove out of bounds element.'
            stop
        endif

        ! swap the item to be deleted with the last element
        ! and shorten heap by one.
        a%elems(i) = a%elems(a%heap_size) 
        a%heap_size = a%heap_size - 1

        call DTYPE(heapify)(a,i)

    end subroutine DTYPE(pq_delete)
