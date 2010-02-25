      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_util_collapse_list
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Generates a list whose entries are unique
      !                 from a list of integer vectors.
      !
      !  Input        : inlist(:,:)  (I) original list of n-plets. 1st
      !                                  index: 1...n; 2nd: 1...nin.
      !                 nin          (I) length of inlist.
      !
      !  Output       : outlist(:,:) (I) collapsed list of n-plets. The 
      !                                  list is sorted and has the same
      !                                  LBOUND as inlist.
      !                 nout         (I) length of outlist, i.e. UBOUND
      !                                  = LBOUND + nout - 1.
      !                 idx(:)       (I) Index array such that 
      !                                   outlist(i) = inlist(idx(i))
      !                                  OPTIONAL. Only returned if
      !                                  present.
      !                 info         (I) return status. 
      !
      !  Remarks      : 
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_util_collapse_list.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.3  2005/03/12 01:28:24  ivos
      !  Added nin,nout, and idx.
      !
      !  Revision 1.2  2005/02/10 16:56:25  ivos
      !  Removed TAB characters in source text.
      !
      !  Revision 1.1  2005/02/09 16:29:53  pchatela
      !  Initial insertion
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_util_collapse_list(inlist,nin,outlist,nout,info,idx)

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_util_qsort
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      ! non-negative integer list to be inverted
      INTEGER, DIMENSION(:,:) , POINTER       :: inlist
      ! inverse list
      INTEGER, DIMENSION(:,:) , POINTER       :: outlist
      INTEGER                 , INTENT(IN   ) :: nin
      INTEGER, DIMENSION(:  ) , POINTER, OPTIONAL :: idx
      ! return status
      INTEGER                 , INTENT(  OUT) :: nout,info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)                 :: t0
      INTEGER, DIMENSION(2)                 :: ldl,ldu
      INTEGER                               :: iopt,l1,u1
      INTEGER                               :: i,j,tmpmin,inmin,inmax
      INTEGER                               :: tmpmax,tmpfact,prevdist
      INTEGER, DIMENSION(:), POINTER        :: ids, sort
      INTEGER, DIMENSION(:), POINTER        :: spans,mins
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialization
      !-------------------------------------------------------------------------
      CALL substart('ppm_util_collapse_list',t0,info)
      
      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT.ASSOCIATED(inlist)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_collapse_list',  &
     &            'inlist must be associated/allocated',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (SIZE(inlist,2) .LT. nin) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_collapse_list',  &
     &            'inlist must be at least of length nin',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Determine geometry of inlist
      !-------------------------------------------------------------------------
      inmin = LBOUND(inlist,2)
      inmax = inmin + nin - 1
      l1    = LBOUND(inlist,1)
      u1    = UBOUND(inlist,1)
      
      !-------------------------------------------------------------------------
      !  If inlist contains only one element, we are done
      !-------------------------------------------------------------------------
      IF (nin .LT. 2) THEN
          nout = nin
          iopt = ppm_param_alloc_fit
          ldl(1) = l1
          ldu(1) = u1
          ldl(2) = inmin
          ldu(2) = inmin + nout - 1
          CALL ppm_alloc(outlist,ldl,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_util_collapse_list',        &
     &            'collapsed list OUTLIST',__LINE__,info)
              GOTO 9999
          ENDIF
          outlist = inlist
          IF (PRESENT(idx)) THEN
              iopt = ppm_param_alloc_fit
              ldl(1) = inmin
              ldu(1) = inmin + nout - 1
              CALL ppm_alloc(idx,ldl,ldu,iopt,info)
              IF (info .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'ppm_util_collapse_list',    &
     &                'Permutation index IDX',__LINE__,info)
                  GOTO 9999
              ENDIF
              idx(inmin) = inmin
          ENDIF
          GOTO 9999 
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate work memory
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldl(1) = inmin
      ldu(1) = inmax
      CALL ppm_alloc(ids,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_collapse_list',        &
     &        'indices list IDS',__LINE__,info)
          GOTO 9999
      ENDIF
      ldl(1) = l1
      ldu(1) = u1
      CALL ppm_alloc(spans,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_collapse_list',        &
     &        'spans list SPANS',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(mins,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_collapse_list',        &
     &        'mins list MINS',__LINE__,info)
          GOTO 9999
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Determine the min and max values in inlist
      !  and compute the spans of the values, needed to generate unique keys
      !-------------------------------------------------------------------------
      DO j=l1,u1
          tmpmax = MAXVAL(inlist(j,:))
          tmpmin = MINVAL(inlist(j,:))
          mins(j) = tmpmin
          spans(j) = tmpmax-tmpmin+1
      ENDDO
      
      !-------------------------------------------------------------------------
      !  Set up the list to carry out the sorting
      !-------------------------------------------------------------------------
      ids(:) = inlist(l1,:)-mins(l1)
      tmpfact = spans(l1)
      DO j=l1+1,u1
           ids(:) = ids(:) + tmpfact * (inlist(j,:)-mins(j))
           tmpfact = tmpfact * spans(j)
      ENDDO
      
      !-------------------------------------------------------------------------
      !  Quicksort the IDs
      !-------------------------------------------------------------------------
      CALL ppm_util_qsort(ids,sort,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_util_collapse_list',        &
     &        'qsort of IDS failed',__LINE__,info)
          GOTO 9999
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Compute new size
      !-------------------------------------------------------------------------
      prevdist = ids(sort(LBOUND(ids,1)))
      nout = 1
      DO i = LBOUND(ids,1)+1,UBOUND(ids,1)
          IF (ids(sort(i)).NE.prevdist) THEN
              nout = nout + 1
              prevdist = ids(sort(i))
          ENDIF
      ENDDO
      
      !-------------------------------------------------------------------------
      !  Allocate new collapsed list 
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldl(1) = l1
      ldu(1) = u1
      ldl(2) = inmin
      ldu(2) = inmin + nout - 1
      CALL ppm_alloc(outlist,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_collapse_list',        &
     &        'collapsed list OUTLIST',__LINE__,info)
          GOTO 9999
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Allocate index list if needed
      !-------------------------------------------------------------------------
      IF (PRESENT(idx)) THEN
          iopt = ppm_param_alloc_fit
          ldl(1) = inmin
          ldu(1) = inmin + nout - 1
          CALL ppm_alloc(idx,ldl,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_util_collapse_list',        &
     &            'Permutation index IDX',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Fill new collapsed list
      !-------------------------------------------------------------------------
      i = inmin
      j = LBOUND(sort,1)
      outlist(:,i) = inlist(:,sort(j))
      prevdist = ids(sort(j))
      DO 
         IF ((i-inmin+1).EQ.nout) EXIT
         j = j+1
         IF (ids(sort(j)).NE.prevdist) THEN
             i = i+1
             outlist(:,i) = inlist(:,sort(j))
             IF (PRESENT(idx)) THEN
                 idx(i)   = sort(j)
             ENDIF
             prevdist     = ids(sort(j))
         ENDIF
      ENDDO
      
      !-------------------------------------------------------------------------
      !  Deallocate stuff
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(ids,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_util_collapse_list',        &
     &        'indices list IDS',__LINE__,info)
      ENDIF
      CALL ppm_alloc(spans,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_util_collapse_list',        &
     &        'spans list SPANS',__LINE__,info)
      ENDIF
      CALL ppm_alloc(mins,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_util_collapse_list',        &
     &        'mins list MINS',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_collapse_list',t0,info)
      RETURN
      END SUBROUTINE ppm_util_collapse_list

