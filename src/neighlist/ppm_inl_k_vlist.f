      !-------------------------------------------------------------------------
      !  Subroutines :               ppm_inl_k_vlist
      !-------------------------------------------------------------------------
      ! Copyright (c) 2014 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
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
      !-------------------------------------------------------------------------
      !  MOSAIC Group
      !  Max Planck Institute of Molecular Cell Biology and Genetics
      !  Pfotenhauerstr. 108, 01307 Dresden, Germany
      !-------------------------------------------------------------------------

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
      SUBROUTINE DTYPE(process_terminal_node)(node)
      ! Look for actual near neighbors in 'node', and update
      ! the search results on the sr data structure.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------

      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      TYPE(DTYPE(tree_node)), POINTER :: node

      !-------------------------------------------------------------------------
      !  Local variables, arrays and counters
      !-------------------------------------------------------------------------
      TYPE(DTYPE(tree_search_record)), POINTER :: sr

      TYPE(pq), POINTER :: pqp

      REAL(MK), DIMENSION(3)            :: qv
      REAL(MK), DIMENSION(:,:), POINTER :: data
      REAL(MK)                          :: ballsize,sd,newpri

      INTEGER, DIMENSION(:), POINTER :: ind
      INTEGER                        :: dimen,i,indexofi,k
      INTEGER                        :: centeridx,correltime,nn

      LOGICAL :: rearrange

      sr => DTYPE(sr)

      !
      ! copy values from sr to local variables
      !
      qv         = sr%qv
      pqp        => sr%pq
      dimen      = sr%dimen
      ballsize   = sr%ballsize
      rearrange  = sr%rearrange
      ind        => sr%ind(1:)
      data       => sr%data(1:,1:)
      centeridx  = sr%centeridx
      correltime = sr%correltime

      ! search through terminal bucket.
      mainloop: DO i = node%l,node%u
         IF (rearrange) THEN
            sd =    (data(1,i) - qv(1))**2
            IF (sd.GT.ballsize) CYCLE mainloop

            sd = sd+(data(2,i) - qv(2))**2
            IF (sd.GT.ballsize) CYCLE mainloop

            DO k = 3,dimen
               sd = sd + (data(k,i) - qv(k))**2
               IF (sd.GT.ballsize) CYCLE mainloop
            ENDDO

            indexofi = ind(i)  ! only read it if we have not broken out
         ELSE
            indexofi = ind(i)

            sd =    (data(1,indexofi) - qv(1))**2
            IF (sd.GT.ballsize) CYCLE mainloop

            sd = sd+(data(2,indexofi) - qv(2))**2
            IF (sd.GT.ballsize) CYCLE mainloop

            DO k = 3,dimen
               sd = sd + (data(k,indexofi) - qv(k))**2
               IF (sd.GT.ballsize) CYCLE mainloop
            ENDDO
         ENDIF !rearrange

         IF (centeridx.GT.0) THEN ! doing correlation interval?
            IF (ABS(indexofi-centeridx).LT.correltime) CYCLE mainloop
         ENDIF

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

         IF (sr%nfound.LT.sr%nn) THEN
            ! add this point unconditionally to fill list.
            sr%nfound = sr%nfound +1
            newpri = pqp%insert(sd,indexofi)
            IF (sr%nfound.EQ.sr%nn) ballsize = newpri
            ! we have just filled the working list.
            ! put the best square distance to the maximum value
            ! on the list, which is extractable from the PQ.
         ELSE
            !
            ! now, if we get here,
            ! we know that the current node has a squared
            ! distance smaller than the largest one on the list, and
            ! belongs on the list.
            ! Hence we replace that with the current one.
            !
            ballsize = pqp%replace_max(sd,indexofi)
         ENDIF
      ENDDO mainloop
      !
      ! Reset sr variables which may have changed during loop
      !
      sr%ballsize = ballsize

      END SUBROUTINE DTYPE(process_terminal_node)

      SUBROUTINE DTYPE(process_terminal_node_fixedball)(node)
      ! Look for actual near neighbors in 'node', and update
      ! the search results on the sr data structure, i.e.
      ! save all within a fixed ball.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------

      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      TYPE(DTYPE(tree_node)), POINTER :: node

      !-------------------------------------------------------------------------
      !  Local variables, arrays and counters
      !-------------------------------------------------------------------------
      TYPE(DTYPE(tree_search_record)), POINTER :: sr

      REAL(MK), DIMENSION(3)            :: qv
      REAL(MK), DIMENSION(:,:), POINTER :: data
      REAL(MK)                          :: ballsize,sd

      INTEGER, DIMENSION(:), POINTER :: ind
      INTEGER                        :: nfound
      INTEGER                        :: dimen,i,indexofi,k
      INTEGER                        :: centeridx, correltime, nn

      LOGICAL :: rearrange

      sr        => DTYPE(sr)
      !
      ! copy values from sr to local variables
      !
      qv        = sr%qv
      dimen     = sr%dimen
      ballsize  = sr%ballsize
      rearrange = sr%rearrange
      ind       => sr%ind(1:)
      data      => sr%data(1:,1:)
      centeridx = sr%centeridx
      correltime= sr%correltime
      nn        = sr%nn ! number to search for
      nfound    = sr%nfound

      ! search through terminal bucket.
      mainloop: DO i = node%l,node%u
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

         IF (rearrange) THEN
            sd =    (data(1,i) - qv(1))**2
            IF (sd.GT.ballsize) CYCLE mainloop

            sd = sd+(data(2,i) - qv(2))**2
            IF (sd.GT.ballsize) CYCLE mainloop

            DO k = 3,dimen
               sd = sd + (data(k,i) - qv(k))**2
               IF (sd.GT.ballsize) CYCLE mainloop
            ENDDO

            indexofi = ind(i)  ! only read it if we have not broken out
         ELSE
            indexofi = ind(i)

            sd =    (data(1,indexofi) - qv(1))**2
            IF (sd.GT.ballsize) CYCLE mainloop

            sd = sd+(data(2,indexofi) - qv(2))**2
            IF (sd.GT.ballsize) CYCLE mainloop

            DO k = 3,dimen
               sd = sd + (data(k,indexofi) - qv(k))**2
               IF (sd.GT.ballsize) CYCLE mainloop
            ENDDO
         ENDIF !rearrange

         IF (centeridx.GT.0) THEN ! doing correlation interval?
            IF (ABS(indexofi-centeridx).LT.correltime) CYCLE mainloop
         ENDIF

         nfound = nfound+1
         IF (nfound.GT.sr%nalloc) THEN
            ! oh nuts, we have to add another one to the tree but
            ! there isn't enough room.
            sr%overflow = .TRUE.
         ELSE
            sr%results(nfound)=DTYPE(kdtree_result)(sd,indexofi)
         ENDIF
      ENDDO mainloop
      !
      ! Reset sr variables which may have changed during loop
      !
      sr%nfound = nfound
      END SUBROUTINE DTYPE(process_terminal_node_fixedball)

      RECURSIVE SUBROUTINE DTYPE(search)(node)
      ! This is the innermost core routine of the kd-tree search.  Along
      ! with "process_terminal_node", it is the performance bottleneck.
      !
      ! This version uses a logically complete secondary search of
      ! "box in bounds", whether the sear
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------

      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      TYPE(DTYPE(tree_node)), POINTER :: node

      !-------------------------------------------------------------------------
      !  Local variables, arrays and counters
      !-------------------------------------------------------------------------
      TYPE(DTYPE(interval)), DIMENSION(:), POINTER :: box

      TYPE(DTYPE(tree_search_record)), POINTER :: sr

      TYPE(DTYPE(tree_node)), POINTER :: ncloser
      TYPE(DTYPE(tree_node)), POINTER :: nfarther

      REAL(MK), DIMENSION(3) :: qv
      REAL(MK)               :: qval,dis,tdis,ballsize

      INTEGER :: cut_dim, i

      sr => DTYPE(sr)

      IF (.NOT.(ASSOCIATED(node%left).AND.ASSOCIATED(node%right))) THEN
         ! we are on a terminal node
         IF (sr%nn.EQ.0) THEN
            CALL DTYPE(process_terminal_node_fixedball)(node)
         ELSE
            CALL DTYPE(process_terminal_node)(node)
         ENDIF
      ELSE
         ! we are not on a terminal node
         qv      =  sr%qv
         cut_dim =  node%cut_dim
         qval    =  qv(cut_dim)

         IF (qval.LT.node%cut_val) THEN
            ncloser  => node%left
            nfarther => node%right
            dis = (node%cut_val_right - qval)**2
         ELSE
            ncloser  => node%right
            nfarther => node%left
            dis = (node%cut_val_left - qval)**2
         ENDIF

         IF (ASSOCIATED(ncloser)) CALL DTYPE(search)(ncloser)

         ! we may need to search the second node.
         IF (ASSOCIATED(nfarther)) THEN
            ballsize = sr%ballsize
            IF (dis.LE.ballsize) THEN
               !
               ! we do this separately as going on the first cut dimen is often
               ! a good idea.
               !
               box => node%box(1:)

               SELECT CASE (sr%dimen)
               CASE (2)
                  SELECT CASE (cut_dim)
                  CASE (1)
                     IF      (qv(2).GT.box(2)%upper) THEN
                        tdis=qv(2)-box(2)%upper
                        dis =dis + tdis*tdis
                     ELSE IF (qv(2).LT.box(2)%lower) THEN
                        tdis=box(2)%lower-qv(2)
                        dis =dis + tdis*tdis
                     ENDIF
                     IF (dis.GT.ballsize) THEN
                        RETURN
                     ENDIF

                  CASE (2)
                     IF      (qv(1).GT.box(1)%upper) THEN
                        tdis=qv(1)-box(1)%upper
                        dis =dis + tdis*tdis
                     ELSE IF (qv(1).LT.box(1)%lower) THEN
                        tdis=box(1)%lower-qv(1)
                        dis =dis + tdis*tdis
                     ENDIF
                     IF (dis.GT.ballsize) THEN
                        RETURN
                     ENDIF

                  END SELECT

               CASE (3)
                  SELECT CASE (cut_dim)
                  CASE (1)
                     IF      (qv(2).GT.box(2)%upper) THEN
                        tdis=qv(2)-box(2)%upper
                        dis =dis + tdis*tdis
                     ELSE IF (qv(2).LT.box(2)%lower) THEN
                        tdis=box(2)%lower-qv(2)
                        dis =dis + tdis*tdis
                     ENDIF
                     IF (dis.GT.ballsize) THEN
                        RETURN
                     ENDIF

                     IF      (qv(3).GT.box(3)%upper) THEN
                        tdis=qv(3)-box(3)%upper
                        dis =dis + tdis*tdis
                     ELSE IF (qv(3).LT.box(3)%lower) THEN
                        tdis=box(3)%lower-qv(3)
                        dis =dis + tdis*tdis
                     ENDIF
                     IF (dis.GT.ballsize) THEN
                        RETURN
                     ENDIF

                  CASE (2)
                     IF      (qv(1).GT.box(1)%upper) THEN
                        tdis=qv(1)-box(1)%upper
                        dis =dis + tdis*tdis
                     ELSE IF (qv(1).LT.box(1)%lower) THEN
                        tdis=box(1)%lower-qv(1)
                        dis =dis + tdis*tdis
                     ENDIF
                     IF (dis.GT.ballsize) THEN
                        RETURN
                     ENDIF

                     IF      (qv(3).GT.box(3)%upper) THEN
                        tdis=qv(3)-box(3)%upper
                        dis =dis + tdis*tdis
                     ELSE IF (qv(3).LT.box(3)%lower) THEN
                        tdis=box(3)%lower-qv(3)
                        dis =dis + tdis*tdis
                     ENDIF
                     IF (dis.GT.ballsize) THEN
                        RETURN
                     ENDIF

                  CASE (3)
                     IF      (qv(1).GT.box(1)%upper) THEN
                        tdis=qv(1)-box(1)%upper
                        dis =dis + tdis*tdis
                     ELSE IF (qv(1).LT.box(1)%lower) THEN
                        tdis=box(1)%lower-qv(1)
                        dis =dis + tdis*tdis
                     ENDIF
                     IF (dis.GT.ballsize) THEN
                        RETURN
                     ENDIF

                     IF      (qv(2).GT.box(2)%upper) THEN
                        tdis=qv(2)-box(2)%upper
                        dis =dis + tdis*tdis
                     ELSE IF (qv(2).LT.box(2)%lower) THEN
                        tdis=box(2)%lower-qv(2)
                        dis =dis + tdis*tdis
                     ENDIF
                     IF (dis.GT.ballsize) THEN
                        RETURN
                     ENDIF

                  END SELECT !(cut_dim)

               END SELECT !(sr%dimen)
               !
               ! if we are still here then we need to search more.
               !
               CALL DTYPE(search)(nfarther)

            ENDIF !dis.LE.ballsize
         ENDIF !ASSOCIATED(nfarther)
      ENDIF !.NOT.(ASSOCIATED(node%left).AND.ASSOCIATED(node%right))
      END SUBROUTINE DTYPE(search)

      SUBROUTINE DTYPE(kdtree_n_nearest)(tp,qv,nn,results,info)
      ! Find the 'nn' vectors in the tree nearest to 'qv' in euclidean norm
      ! returning their indexes and distances in 'indexes' and 'distances'
      ! arrays already allocated passed to this subroutine.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------

      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      TYPE(DTYPE(kdtree)),                      POINTER       :: tp

      REAL(MK),                   DIMENSION(:), INTENT(IN   ) :: qv

      INTEGER,                                  INTENT(IN   ) :: nn

      TYPE(DTYPE(kdtree_result)), DIMENSION(:), TARGET        :: results

      INTEGER,                                  INTENT(  OUT) :: info

      !-------------------------------------------------------------------------
      !  Local variables, arrays and counters
      !-------------------------------------------------------------------------
      TYPE(DTYPE(tree_search_record)), POINTER :: sr

      REAL(ppm_kind_double) :: t0
      REAL(MK),   PARAMETER :: big=REAL(HUGE(1.0),MK)

      CHARACTER(LEN=ppm_char) :: caller='kdtree_n_nearest'

      CALL substart(caller,t0,info)

      sr            => DTYPE(sr)

      sr%dimen      = tp%dimen
      sr%nn         = nn
      sr%nfound     = 0
      sr%centeridx  =-1
      sr%correltime = 0
      sr%nalloc     = nn   ! will be checked
      sr%ind        => tp%ind
      sr%ballsize   = big
      sr%qv(1:tp%dimen)=qv(1:tp%dimen)
      IF (tp%rearrange) THEN
         sr%data    => tp%rearranged_data
      ELSE
         sr%data    => tp%the_data
      ENDIF

      sr%results    => results
      sr%rearrange  = tp%rearrange
      sr%overflow   = .FALSE.

      check_false(<#SIZE(results,1).LT.nn#>, &
      & "KD_TREE_TRANS:  you did not provide enough storage for results(1:nn)")

      !This part is done only for intel compiler
      ASSOCIATE(pq_ => sr%pq)
        CALL pq_%create(results)
      END ASSOCIATE

      CALL DTYPE(search)(tp%root)

      IF (tp%sort) THEN
         CALL qsort(results(1:nn),info)
         or_fail("kdtree_util_qsort")
      ENDIF

      ! DEALLOCATE(sr%pq)
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE DTYPE(kdtree_n_nearest)

      SUBROUTINE DTYPE(kdtree_n_nearest_around_point)(tp,idxin,correltime,nn,results,info)
      ! Find the 'nn' vectors in the tree nearest to point 'idxin',
      ! with correlation window 'correltime', returning results in
      ! results(:), which must be pre-allocated upon entry.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------

      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      TYPE(DTYPE(kdtree)),                      POINTER       :: tp

      INTEGER,                                  INTENT(IN   ) :: idxin
      INTEGER,                                  INTENT(IN   ) :: correltime
      INTEGER,                                  INTENT(IN   ) :: nn

      TYPE(DTYPE(kdtree_result)), DIMENSION(:), TARGET        :: results

      INTEGER,                                  INTENT(  OUT) :: info

      !-------------------------------------------------------------------------
      !  Local variables, arrays and counters
      !-------------------------------------------------------------------------
      TYPE(DTYPE(tree_search_record)), POINTER :: sr

      REAL(ppm_kind_double) :: t0
      REAL(MK), PARAMETER   :: big=REAL(HUGE(1.0),MK)

      INTEGER :: iopt,ldu(1)

      CHARACTER(LEN=ppm_char) :: caller='kdtree_n_nearest_around_point'

      CALL substart(caller,t0,info)

      sr            => DTYPE(sr)
      sr%dimen      = tp%dimen
      sr%nn         = nn
      sr%nfound     = 0
      sr%centeridx  =idxin
      sr%correltime =correltime
      sr%nalloc     = nn   ! will be checked
      sr%ind        => tp%ind
      sr%ballsize   = big
      sr%qv(1:tp%dimen)=tp%the_data(:,idxin) ! copy the vector
      IF (tp%rearrange) THEN
         sr%data    => tp%rearranged_data
      ELSE
         sr%data    => tp%the_data
      ENDIF

      sr%results    => results
      sr%rearrange  = tp%rearrange

      check_false(<#SIZE(results,1).LT.nn#>, &
      & "KD_TREE_TRANS:  you did not provide enough storage for results(1:nn)")

      !This part is done only for intel compiler
      ASSOCIATE(pq_ => sr%pq)
        CALL pq_%create(results)
      END ASSOCIATE

      CALL DTYPE(search)(tp%root)

      IF (tp%sort) THEN
         CALL qsort(results(1:nn),info)
         or_fail("kdtree_util_qsort")
      ENDIF

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE DTYPE(kdtree_n_nearest_around_point)

      SUBROUTINE DTYPE(kdtree_r_nearest)(tp,qv,r2,nalloc,nfound,results,info)
      ! Find the nearest neighbors to point 'idxin', within SQUARED
      ! Euclidean distance 'r2'.   Upon ENTRY, nalloc must be the
      ! size of memory allocated for results(1:nalloc).  Upon
      ! EXIT, nfound is the number actually found within the ball.
      !
      !  Note that if nfound .GT. nalloc then more neighbors were found
      !  than there were storage to store.  The resulting list is NOT
      !  the smallest ball inside norm r^2
      !
      ! Results are NOT sorted unless tree was created with sort option.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------

      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      TYPE(DTYPE(kdtree)),                      POINTER       :: tp

      REAL(MK),                   DIMENSION(:), INTENT(IN   ) :: qv
      REAL(MK),                                 INTENT(IN   ) :: r2

      INTEGER,                                  INTENT(IN   ) :: nalloc
      INTEGER,                                  INTENT(  OUT) :: nfound

      TYPE(DTYPE(kdtree_result)), DIMENSION(:), TARGET        :: results

      INTEGER,                                  INTENT(  OUT) :: info

      !-------------------------------------------------------------------------
      !  Local variables, arrays and counters
      !-------------------------------------------------------------------------
      TYPE(DTYPE(tree_search_record)), POINTER :: sr

      REAL(ppm_kind_double) :: t0

      INTEGER :: iopt,ldu(1)

      CHARACTER(LEN=ppm_char) :: caller='kdtree_r_nearest'

      CALL substart(caller,t0,info)

      sr            => DTYPE(sr)

      sr%dimen      = tp%dimen
      sr%nn         = 0 ! flag for fixed ball search
      sr%nfound     = 0
      sr%centeridx  =-1
      sr%correltime = 0

      check_false(<#SIZE(results,1).LT.nalloc#>, &
      & "KD_TREE_TRANS:  you did not provide enough storage for results(1:nalloc)")

      sr%nalloc     = nalloc
      sr%ind        => tp%ind
      sr%ballsize   = r2
      sr%qv(1:tp%dimen)=qv(1:tp%dimen)
      IF (tp%rearrange) THEN
         sr%data    => tp%rearranged_data
      ELSE
         sr%data    => tp%the_data
      ENDIF

      sr%results    => results
      sr%rearrange  = tp%rearrange
      sr%overflow   = .FALSE.

      CALL DTYPE(search)(tp%root)

      nfound = sr%nfound

      IF (tp%sort) THEN
         CALL qsort(results(1:nfound),info)
         or_fail("kdtree_util_qsort")
      ENDIF

      IF (sr%overflow) THEN
         stdout("KD_TREE_TRANS: warning! return from kdtree_r_nearest found more neighbors")
         stdout("KD_TREE_TRANS: than storage was provided for.  Answer is NOT smallest ball")
         stdout("KD_TREE_TRANS: with that number of neighbors!  I.e. it is wrong.")
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE DTYPE(kdtree_r_nearest)

      SUBROUTINE DTYPE(kdtree_r_nearest_around_point)(tp, &
      &          idxin,correltime,nalloc,r2,nfound,results,info)
      ! Like kdtree_r_nearest, but around a point 'idxin' already existing
      ! in the data set.
      !
      ! Results are NOT sorted unless tree was created with sort option.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------

      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      TYPE(DTYPE(kdtree)),                      POINTER       :: tp

      INTEGER,                                  INTENT(IN   ) :: idxin
      INTEGER,                                  INTENT(IN   ) :: correltime
      INTEGER,                                  INTENT(IN   ) :: nalloc

      REAL(MK),                                 INTENT(IN   ) :: r2

      INTEGER,                                  INTENT(  OUT) :: nfound

      TYPE(DTYPE(kdtree_result)), DIMENSION(:), TARGET        :: results

      INTEGER,                                  INTENT(  OUT) :: info

      !-------------------------------------------------------------------------
      !  Local variables, arrays and counters
      !-------------------------------------------------------------------------
      TYPE(DTYPE(tree_search_record)), POINTER :: sr

      REAL(ppm_kind_double) :: t0

      INTEGER :: iopt,ldu(1)

      CHARACTER(LEN=ppm_char) :: caller='kdtree_r_nearest_around_point'

      CALL substart(caller,t0,info)

      sr            => DTYPE(sr)
      sr%dimen      = tp%dimen
      sr%nn         = 0 ! flag for fixed r search
      sr%nfound     = 0
      sr%centeridx  =idxin
      sr%correltime =correltime

      check_false(<#SIZE(results,1).LT.nalloc#>, &
      & "KD_TREE_TRANS:  you did not provide enough storage for results(1:nalloc)")

      sr%nalloc     = nalloc
      sr%ind        => tp%ind
      sr%ballsize   = r2
      sr%qv(1:tp%dimen)=tp%the_data(:,idxin) ! copy the vector
      IF (tp%rearrange) THEN
         sr%data    => tp%rearranged_data
      ELSE
         sr%data    => tp%the_data
      ENDIF

      sr%results    => results
      sr%rearrange  = tp%rearrange
      sr%overflow   = .FALSE.

      CALL DTYPE(search)(tp%root)

      nfound = sr%nfound

      IF (tp%sort) THEN
         CALL qsort(results(1:nfound),info)
         or_fail("kdtree_util_qsort")
      ENDIF

      IF (sr%overflow) THEN
         stdout("KD_TREE_TRANS: warning! return from kdtree_r_nearest found more neighbors")
         stdout("KD_TREE_TRANS: than storage was provided for.  Answer is NOT smallest ball")
         stdout("KD_TREE_TRANS: with that number of neighbors!  I.e. it is wrong.")
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE DTYPE(kdtree_r_nearest_around_point)

      SUBROUTINE DTYPE(kdtree_r_count)(tp,qv,r2,nfound,info)
      ! Count the number of neighbors within square distance 'r2'.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------

      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      TYPE(DTYPE(kdtree)),    POINTER       :: tp

      REAL(MK), DIMENSION(:), INTENT(IN   ) :: qv
      REAL(MK),               INTENT(IN   ) :: r2

      INTEGER,                INTENT(  OUT) :: nfound
      INTEGER,                INTENT(  OUT) :: info

      !-------------------------------------------------------------------------
      !  Local variables, arrays and counters
      !-------------------------------------------------------------------------
      TYPE(DTYPE(tree_search_record)), POINTER :: sr

      REAL(ppm_kind_double) :: t0

      INTEGER :: iopt,ldu(1)

      CHARACTER(LEN=ppm_char) :: caller='kdtree_r_count'

      CALL substart(caller,t0,info)

      sr            => DTYPE(sr)
      sr%dimen      = tp%dimen
      sr%nn         = 0 ! flag for fixed r search
      sr%nfound     = 0
      sr%centeridx  =-1
      sr%correltime = 0
      sr%nalloc     = 0 ! we do not allocate any storage but that's OK
                        ! for counting.
      sr%ind        => tp%ind
      sr%ballsize   = r2
      sr%qv(1:tp%dimen)=qv(1:tp%dimen)
      IF (tp%rearrange) THEN
         sr%data    => tp%rearranged_data
      ELSE
         sr%data    => tp%the_data
      ENDIF

      NULLIFY(sr%results)

      sr%rearrange  = tp%rearrange
      sr%overflow   = .FALSE.

      CALL DTYPE(search)(tp%root)

      nfound = sr%nfound

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE DTYPE(kdtree_r_count)

      SUBROUTINE DTYPE(kdtree_r_count_around_point)(tp,idxin,correltime,r2,nfound,info)
      ! Count the number of neighbors within square distance 'r2' around
      ! point 'idxin' with decorrelation time 'correltime'.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------

      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      TYPE(DTYPE(kdtree)), POINTER       :: tp

      INTEGER,             INTENT(IN   ) :: idxin
      INTEGER,             INTENT(IN   ) :: correltime

      REAL(MK),            INTENT(IN   ) :: r2

      INTEGER,             INTENT(  OUT) :: nfound
      INTEGER,             INTENT(  OUT) :: info

      !-------------------------------------------------------------------------
      !  Local variables, arrays and counters
      !-------------------------------------------------------------------------
      TYPE(DTYPE(tree_search_record)), POINTER :: sr

      REAL(ppm_kind_double) :: t0

      INTEGER :: iopt,ldu(1)

      CHARACTER(LEN=ppm_char) :: caller='kdtree_r_count_around_point'

      CALL substart(caller,t0,info)

      sr            => DTYPE(sr)
      sr%dimen      = tp%dimen
      sr%nn         = 0 ! flag for fixed r search
      sr%nfound     = 0
      sr%centeridx  =idxin
      sr%correltime =correltime
      sr%nalloc     = 0
      ! we do not allocate any storage but that's OK for counting.
      sr%ind        => tp%ind
      sr%ballsize   = r2
      sr%qv(1:tp%dimen)=tp%the_data(:,idxin) ! copy the vector
      IF (tp%rearrange) THEN
         sr%data    => tp%rearranged_data
      ELSE
         sr%data    => tp%the_data
      ENDIF

      NULLIFY(sr%results)
      sr%rearrange  = tp%rearrange
      sr%overflow   = .FALSE.

      CALL DTYPE(search)(tp%root)

      nfound = sr%nfound

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE DTYPE(kdtree_r_count_around_point)

      SUBROUTINE DTYPE(kdtree_n_nearest_brute_force)(tp,qv,nn,results,info)
      ! Find the 'n' nearest neighbors to 'qv' by exhaustive search.
      ! only use this subroutine for testing, as it is SLOW!  The
      ! whole point of a k-d tree is to avoid doing what this subroutine
      ! does.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------

      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      TYPE(DTYPE(kdtree)),                      POINTER       :: tp

      REAL(MK),                   DIMENSION(:), INTENT(IN   ) :: qv

      INTEGER,                                  INTENT(IN   ) :: nn

      TYPE(DTYPE(kdtree_result)), DIMENSION(:), TARGET        :: results

      INTEGER,                                  INTENT(  OUT) :: info

      !-------------------------------------------------------------------------
      !  Local variables, arrays and counters
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)               :: t0
      REAL(MK), PARAMETER                 :: big=REAL(HUGE(1.0),MK)
      REAL(MK), DIMENSION(:), ALLOCATABLE :: all_distances

      INTEGER :: i,j,k,d

      CHARACTER(LEN=ppm_char) :: caller='kdtree_n_nearest_brute_force'

      CALL substart(caller,t0,info)

      ALLOCATE(all_distances(tp%n),STAT=info)
      or_fail_alloc("all_distances")

      d=tp%dimen
      FORALL (i=1:tp%n) all_distances(i) = SUM((qv(1:d)-tp%the_data(1:d,i))**2)

      ! now find 'n' smallest distances
      FORALL (i=1:nn) results(i)=DTYPE(kdtree_result)(big,-1)

      DO i =1,tp%n
         IF (all_distances(i).LT.results(nn)%dis) THEN
            ! insert it somewhere on the list
            DO j=1,nn
               IF (all_distances(i).LT.results(j)%dis) EXIT
            ENDDO
            ! now we know 'j'
            DO k = nn - 1, j, -1
               results(k+1) = results(k)
            ENDDO
            results(j)=DTYPE(kdtree_result)(all_distances(i),i)
         ENDIF
      ENDDO

      DEALLOCATE (all_distances,STAT=info)
      or_fail_dealloc("all_distances")
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE DTYPE(kdtree_n_nearest_brute_force)

      SUBROUTINE DTYPE(kdtree_r_nearest_brute_force)(tp,qv,r2,nfound,results,info)
      ! Find the nearest neighbors to 'qv' with distance**2 <= r2 by exhaustive search.
      ! only use this subroutine for testing, as it is SLOW!  The
      ! whole point of a k-d tree is to avoid doing what this subroutine
      ! does.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_util_sort
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      TYPE(DTYPE(kdtree)),                      POINTER       :: tp

      REAL(MK),                   DIMENSION(:), INTENT(IN   ) :: qv
      REAL(MK),                                 INTENT(IN   ) :: r2

      INTEGER,                                  INTENT(  OUT) :: nfound

      TYPE(DTYPE(kdtree_result)), DIMENSION(:), TARGET        :: results

      INTEGER,                                  INTENT(  OUT) :: info

      !-------------------------------------------------------------------------
      !  Local variables, arrays and counters
      !-------------------------------------------------------------------------
      TYPE(DTYPE(tree_search_record)), POINTER :: sr

      REAL(ppm_kind_double)               :: t0
      REAL(MK), DIMENSION(:), ALLOCATABLE :: all_distances

      INTEGER :: iopt,ldu(1),i,d,nalloc

      CHARACTER(LEN=ppm_char) :: caller='kdtree_r_nearest_brute_force'

      CALL substart(caller,t0,info)

      ALLOCATE(all_distances(tp%n),STAT=info)
      or_fail_alloc("all_distances")

      d=tp%dimen
      FORALL (i=1:tp%n) all_distances(i) = SUM((qv(1:d)-tp%the_data(1:d,i))**2)

      nfound = 0
      nalloc = SIZE(results,1)

      DO i = 1, tp%n
         IF (all_distances(i).LT.r2) THEN
            ! insert it somewhere on the list
            IF (nfound.LT.nalloc) THEN
               nfound = nfound+1
               results(nfound)=DTYPE(kdtree_result)(all_distances(i),i)
            ENDIF
         ENDIF
      ENDDO

      DEALLOCATE (all_distances,STAT=info)
      or_fail_dealloc("all_distances")

      CALL qsort(results(1:nfound),info)
      or_fail("kdtree_util_qsort")

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE DTYPE(kdtree_r_nearest_brute_force)

      SUBROUTINE DTYPE(kdtree_util_qsort)(inlist,info)
      !!! From a list of values sorts them into ascending order.
      !!! [NOTE]
      !! Yaser
      !!! this is only the modification of the hybrid QuickSort algorithm
      !!! which will sort the input array in an ascending order.
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------

      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      TYPE(DTYPE(kdtree_result)), DIMENSION(:), INTENT(INOUT) :: inlist
      !!! List to be sorted
      INTEGER,                                   INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      TYPE(DTYPE(kdtree_result))                            :: datap
      TYPE(DTYPE(kdtree_result)), DIMENSION(:), ALLOCATABLE :: inlistarray

      REAL(ppm_kind_double) :: t0

      INTEGER, PARAMETER                 :: m = 9
      !   QuickSort Cutoff
      !   Quit QuickSort-ing when a subsequence contains M or fewer
      !   elements and finish off at end with straight insertion sort.
      !   According to Knuth, V.3, the optimum value of M is around 9.
      INTEGER                            :: stklength, dn
      INTEGER, DIMENSION(:), ALLOCATABLE :: lstk
      INTEGER, DIMENSION(:), ALLOCATABLE :: rstk
      INTEGER, DIMENSION(:), ALLOCATABLE :: inlistidx
      ! Permutation list
      INTEGER                            :: indexp,indext
      INTEGER                            :: i,j,n,l,r,p
      INTEGER                            :: istk
      INTEGER                            :: inlistu,inlistl

      CHARACTER(LEN=ppm_char) :: caller='kdtree_util_qsort'
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialization
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      inlistl = LBOUND(inlist,1)
      inlistu = UBOUND(inlist,1)
      n = inlistu-inlistl+1

      ALLOCATE(inlistidx(inlistl:inlistu),STAT=info)
      or_fail_alloc('indices list inlistidx',ppm_error=ppm_error_fatal)

      FORALL (i=inlistl:inlistu) inlistidx(i)=i

      !-------------------------------------------------------------------------
      ! Compute the log of the list length, it is in fact a bound for the length
      ! of the stack of intervals
      !-------------------------------------------------------------------------
      stklength = 1
      dn = n
      DO
         IF (dn.EQ.0) EXIT
         dn = ISHFT(dn,-1)
         stklength = stklength + 1
      ENDDO
      stklength = 2*stklength

      !-------------------------------------------------------------------------
      ! Allocate the stack of intervals
      !-------------------------------------------------------------------------
      ALLOCATE(lstk(stklength),rstk(stklength),STAT=info)
      or_fail_alloc('Failed to allocate left indices list LSTK & right indices list RSTK!', &
      & ppm_error=ppm_error_fatal)

      ! If array is short, skip QuickSort and go directly to
      ! the straight insertion sort.
      IF (n.LE.m) GOTO 900

      !-------------------------------------------------------------------------
      !     QuickSort
      !
      !     The Qn:s correspond roughly to steps in Algorithm Q,
      !     Knuth, V.3, PP.116-117, modified to select the median
      !     of the first, last, and middle elements as the pivot
      !     key (in Knuths notation, K).  Also modified to leave
      !     data in place and produce an outlist array (here:OUTLIST).  To
      !     simplify comments, let INLIST[I]=INLIST(outlist(I)).

      ! Q1 : Initialize
      istk = 0
      l = inlistl
      r = inlistu

      200 CONTINUE

      ! Q2: Sort the subsequence INLIST[L]..INLIST[R].
      !
      !     At this point, INLIST[l] <= INLIST[m] <= INLIST[r] for all l < L,
      !     r > R, and L <= m <= R.  (First time through, there is no
      !     DATA for l < L or r > R.)

      i=l
      j=r

      ! Q2.5: Select pivot key
      !
      !     Let the pivot, P, be the midpoint of this subsequence,
      !     P=(L+R)/2; then rearrange outlist(L), outlist(P), and outlist(R)
      !     so the corresponding INLIST values are in increasing order.
      !     The pivot key, DATAP, is then DATA[P].

      p=(l+r)/2

      indexp=inlistidx(p)
      datap=inlist(indexp)

      IF (inlist(inlistidx(l))%dis.GT.datap%dis) THEN
         inlistidx(p)=inlistidx(l)
         inlistidx(l)=indexp
         indexp=inlistidx(p)
         datap=inlist(indexp)
      ENDIF

      IF (datap%dis.GT.inlist(inlistidx(r))%dis) THEN
         IF (inlist(inlistidx(l))%dis.GT.inlist(inlistidx(r))%dis) THEN
            inlistidx(p)=inlistidx(l)
            inlistidx(l)=inlistidx(r)
         ELSE
            inlistidx(p)=inlistidx(r)
         ENDIF
         inlistidx(r)=indexp
         indexp=inlistidx(p)
         datap=inlist(indexp)
      ENDIF

      !     Now we swap values between the right and left sides and/or
      !     move DATAP until all smaller values are on the left and all
      !     larger values are on the right.  Neither the left or right
      !     side will be internally ordered yet; however, DATAP will be
      !     in its final position.
      300 CONTINUE

      ! Q3: Search for datum on left >= DATAP
      !
      !     At this point, DATA[L] <= DATAP.  We can therefore start scanning
      !     up from L, looking for a value >= DATAP (this scan is guaranteed
      !     to terminate since we initially placed DATAP near the middle of
      !     the subsequence).

      i=i+1
      IF (inlist(inlistidx(i))%dis.LT.datap%dis) GOTO 300

      400 CONTINUE

      ! Q4: Search for datum on right <= DATAP
      !
      !     At this point, DATA[R] >= DATAP.  We can therefore start scanning
      !     down from R, looking for a value <= DATAP (this scan is guaranteed
      !     to terminate since we initially placed DATAP near the middle of
      !     the subsequence).
      j=j-1
      IF (inlist(inlistidx(j))%dis.GT.datap%dis) GOTO 400

      ! Q5: Have the two scans collided?
      IF (i.LT.j) THEN
         ! Q6: No, interchange DATA[I] <--> DATA[J] and continue
         !PRINT *, ' Interchange '
         indext=inlistidx(i)
         inlistidx(i)=inlistidx(j)
         inlistidx(j)=indext
         GOTO 300
      ELSE

      ! Q7: Yes, select next subsequence to sort
      !
      !     At this point, I >= J and DATA[l] <= DATA[I] == DATAP <= DATA[r],
      !     for all L <= l < I and J < r <= R.  If both subsequences are
      !     more than M elements long, push the longer one on the stack and
      !     go back to QuickSort the shorter; if only one is more than M
      !     elements long, go back and QuickSort it; otherwise, pop a
      !     subsequence off the stack and QuickSort it.
         !PRINT *, ' Collision '
         IF (r-j .GE. i-l .AND. i-l .GT. m) THEN
            istk=istk+1
            lstk(istk)=j+1
            rstk(istk)=r
            r=i-1
         ELSE IF (i-l .GT. r-j .AND. r-j .GT. m) THEN
            istk=istk+1
            lstk(istk)=l
            rstk(istk)=i-1
            l=j+1
         ELSE IF (r-j .GT. m) THEN
            l=j+1
         ELSE IF (i-l .GT. m) THEN
            r=i-1
         ELSE
            ! Q8: Pop the stack, or terminate QuickSort if empty
            IF (istk.LT.1) GOTO 900
            l=lstk(istk)
            r=rstk(istk)
            istk=istk-1
         ENDIF
         GOTO 200
      ENDIF

      900 CONTINUE

      !-------------------------------------------------------------------------
      ! Q9: Straight Insertion sort
      DO i=inlistl+1,inlistu
         IF (inlist(inlistidx(i-1))%dis.GT.inlist(inlistidx(i))%dis) THEN
            indexp=inlistidx(i)
            datap=inlist(indexp)
            p=i-1

      920   CONTINUE

            inlistidx(p+1) = inlistidx(p)
            p=p-1
            IF (p.GT.(inlistl-1)) THEN
               IF (inlist(inlistidx(p))%dis.GT.datap%dis) GOTO 920
            ENDIF
            inlistidx(p+1) = indexp
         ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      ! Deallocate the stack of intervals
      !-------------------------------------------------------------------------
      DEALLOCATE(lstk,rstk,STAT=info)
      or_fail_dealloc('left indices list LSTK & right indices list RSTK')

      ALLOCATE(inlistarray(inlistl:inlistu),STAT=info)
      or_fail_alloc("inlistarray")

      FORALL(i=inlistl:inlistu)
         inlistarray(i)=inlist(inlistidx(i))
      END FORALL

      FORALL(i=inlistl:inlistu) inlist(i)=inlistarray(i)

      !-------------------------------------------------------------------------
      ! Deallocate the stack of intervals
      !-------------------------------------------------------------------------
      DEALLOCATE(inlistarray,inlistidx,STAT=info)
      or_fail_dealloc('Failed to deallocate inlistarray & inlistidx!')

      !===================================================================
      !
      !     All done

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN

      END SUBROUTINE DTYPE(kdtree_util_qsort)

#undef  __KIND
#undef    MK
#undef   _MK
#undef    DTYPE
