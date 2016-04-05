      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_util_collapse_list
      !-------------------------------------------------------------------------
      ! Copyright (c) 2012 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
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

      SUBROUTINE ppm_util_collapse_list(inlist,nin,outlist,nout,info,idx)
      !!! Generates a list whose entries are unique from a list of
      !!! integer vectors.
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
      INTEGER, DIMENSION(:,:) , POINTER       :: inlist
      !!! Original list of n-plets.
      !!!
      !!! 1st index: 1...n                                                     +
      !!! 2nd index: 1...nin.
      INTEGER, DIMENSION(:,:) , POINTER       :: outlist
      !!! Collapsed list of n-plets. The list is sorted and has the same
      !!! LBOUND as inlist.
      INTEGER                 , INTENT(IN   ) :: nin
      !!! Length of inlist.
      INTEGER, DIMENSION(:  ) , POINTER, OPTIONAL :: idx
      !!! Index array such that outlist(i) = inlist(idx(i))
      INTEGER                 , INTENT(  OUT) :: nout
      !!! Length of outlist, i.e. UBOUND = LBOUND + nout - 1.
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)                 :: t0
      INTEGER, DIMENSION(2)                 :: ldl,ldu
      INTEGER                               :: iopt,l1,u1
      INTEGER                               :: i,j,tmpmin,inmin,inmax
      INTEGER                               :: tmpmax,tmpfact,prevdist
      INTEGER, DIMENSION(:), POINTER        :: ids   => NULL()
      INTEGER, DIMENSION(:), POINTER        :: sort  => NULL()
      INTEGER, DIMENSION(:), POINTER        :: spans => NULL()
      INTEGER, DIMENSION(:), POINTER        :: mins  => NULL()
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
        CALL check
        IF (info .NE. 0) GOTO 9999
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
      CONTAINS
      SUBROUTINE check
          IF (.NOT.ASSOCIATED(inlist)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_collapse_list',  &
     &            'inlist must be associated/allocated',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (SIZE(inlist,2) .LT. nin) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_collapse_list',  &
     &            'inlist must be at least of length nin',__LINE__,info)
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_util_collapse_list

