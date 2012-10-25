      minclude ppm_header(ppm_module_name="loadbal_choose_sub")

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE loadbal_choose_sub_s(topoid,ideal_load,candidate_sublist,     &
     &                                send_sublist,csublist_size,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE loadbal_choose_sub_d(topoid,ideal_load,candidate_sublist,     &
     &                                send_sublist,csublist_size,info)
#endif
      !!! This routine chooses the best possible sequence of subs to be
      !!! sent to the underloaded processor.
      !!!
      !!! [NOTE]
      !!! This routine is called by ppm_loadbal_bc
      !!!

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_loadbal
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_write
      USE ppm_module_util_commopt
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   )     :: topoid
      !!! The topology where the DLB will take place
      REAL(MK)                                    :: ideal_load
      !!! Ideal load that would balance the workload between overloaded and
      !!! underloaded processors
      INTEGER,DIMENSION(:),POINTER,INTENT(IN   )  :: candidate_sublist
      !!! List of the candidate subdomains (local IDs)
      INTEGER,DIMENSION(:),POINTER,INTENT(  OUT)  :: send_sublist
      !!! Final sequence of the subs to be sent (global IDs)
      INTEGER                    , INTENT(IN   )  :: csublist_size
      !!! Size of the candidate subs list
      INTEGER                    , INTENT(  OUT)  :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER                                 :: isub,i,j,max_loc,counter
      INTEGER                                 :: iopt,first2go
      INTEGER,DIMENSION(1)                    :: temp,ldu
      REAL(MK)                                :: uber_load,old_ideal_load
      REAL(MK)                                :: new_ideal_load,temp_load
      REAL(MK),DIMENSION(:),POINTER           :: candidate_costlist
      REAL(MK),DIMENSION(:),POINTER            :: modified_candidates
      TYPE(ppm_t_topo),POINTER                :: topo => NULL()
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      start_subroutine("loadbal_choose_sub")
      !-------------------------------------------------------------------------
      !  Get the topology first
      !-------------------------------------------------------------------------
      topo    => ppm_topo(topoid)%t
      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF
      !-------------------------------------------------------------------------
      !  Conditions when to send a sub
      !  1. DO NOT send anything, leave the send_list unallocated
      !-------------------------------------------------------------------------
      uber_load = 2._MK * ideal_load
!      stdout("csublist_size",csublist_size)
!!      stdout("candidate costs",'ppm_loadbal_subcostd(candidate_sublist(1:csublist_size))')
!      stdout("ideal load",ideal_load)

      !-------------------------------------------------------------------------
      !  Create the candidates' cost list
      !-------------------------------------------------------------------------

      iopt = ppm_param_alloc_fit
      ldu(1) = csublist_size
      CALL ppm_alloc(candidate_costlist,ldu,iopt,info)
      or_fail_alloc("candidate_costlist")
      DO i=1,csublist_size
#if    __KIND == __SINGLE_PRECISION
         candidate_costlist(i) = ppm_loadbal_subcosts(candidate_sublist(i))
#else
         candidate_costlist(i) = ppm_loadbal_subcostd(candidate_sublist(i))
#endif
      ENDDO
!      stdout("candidate_costlist:",'candidate_costlist')
!      stdout("MAXVAL(candidate_costlist)",'MAXVAL(candidate_costlist)')
!      stdout("MINVAL(candidate_costlist)",'MINVAL(candidate_costlist)')
      !-------------------------------------------------------------------------
      !  Check different cases..
      !  1. If the lightest sub is more than twice costly than the ideal load
      !     then don't send anything. Skip DLB this time
      !-------------------------------------------------------------------------
#if    __KIND == __SINGLE_PRECISION
      IF (MINVAL(candidate_costlist).GT.uber_load .OR. &
     &           csublist_size.EQ.0) THEN
#else
      IF (MINVAL(candidate_costlist).GT.uber_load .OR. &
     &           csublist_size.EQ.0) THEN
#endif
        stdout("Min balancing load > 2*ideal_load.. No DLB this time!")
        iopt = ppm_param_alloc_fit
        ldu(1) = 1
        CALL ppm_alloc(send_sublist,ldu,iopt,info)
        or_fail_alloc("send_sublist")
        send_sublist(1) = -1
        GOTO 9999
      !-------------------------------------------------------------------------
      !  2. Send only ONE sub if even the minimum subcost is more than the
      !     ideal_load (it's intrinsicly less than uber_load since I check this
      !     before)
      !     The sub should be the one with minimum cost!
      !-------------------------------------------------------------------------
      ELSE IF (MINVAL(candidate_costlist).GT.ideal_load) THEN
        stdout("Balancing load > ideal_load.. I send a SINGLE sub")
        iopt = ppm_param_alloc_fit
        ldu(1) = 1
        CALL ppm_alloc(send_sublist,ldu,iopt,info)
        or_fail_alloc("send_sublist")
        temp            = MINLOC(candidate_costlist)
        ! send the global ID
        send_sublist(1) = topo%isublist(candidate_sublist(temp(1)))

      ELSE
        stdout("ELSE case...")
        iopt = ppm_param_alloc_fit
        ldu(1) = csublist_size
        CALL ppm_alloc(send_sublist,ldu,iopt,info)
        or_fail_alloc("send_sublist")
        old_ideal_load = ideal_load
        new_ideal_load = ideal_load - ppm_myepsd
        !-----------------------------------------------------------------------
        !  Create a candidates list that will hold the absolute distances
        !  away from the ideal load
        !-----------------------------------------------------------------------
        iopt = ppm_param_alloc_fit
        ldu(1) = csublist_size
        CALL ppm_alloc(modified_candidates,ldu,iopt,info)
        or_fail_alloc("modified_candidates")
        counter = 0
        send_sublist = -1
        modified_candidates = 0._MK
!        stdout("old_ideal",old_ideal_load,"new_ideal",new_ideal_load)
        DO WHILE(new_ideal_load.LT.old_ideal_load .AND. &
     &           MINVAL(modified_candidates).LE.new_ideal_load )

            !-------------------------------------------------------------------
            !  Initialize arrays
            !
            !-------------------------------------------------------------------
            old_ideal_load = new_ideal_load
            modified_candidates =  ABS(candidate_costlist - old_ideal_load)
!            print*,'modified candidates:',modified_candidates
            IF (MINVAL(modified_candidates).LT.new_ideal_load) THEN
                temp = MINLOC(modified_candidates)
                first2go = temp(1)
                ! Add this sub to the list
                counter = counter + 1
                send_sublist(counter) = candidate_sublist(first2go)
                new_ideal_load = modified_candidates(first2go)
!                stdout("old_ideal",old_ideal_load,"new_ideal",new_ideal_load)
                candidate_costlist(first2go) = 999._MK*old_ideal_load

           ELSE
                GOTO 6666
           ENDIF
        ENDDO
 6666   CONTINUE
        iopt = ppm_param_alloc_fit_preserve
        ldu(1) = counter
        CALL ppm_alloc(send_sublist,ldu,iopt,info)
        or_fail_alloc("send_sublist")
        !----------------------------------------------------------------------
        !  Return global IDs of the subs in the send_sublist
        !----------------------------------------------------------------------
        DO i=1,SIZE(send_sublist,1)
            send_sublist(i)  = topo%isublist(send_sublist(i))
        ENDDO
!      ELSE
!        stdout("An unknown case...oh wait..have I forgotten to implement this??")
      ENDIF


      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      end_subroutine()
      RETURN
      CONTAINS
      SUBROUTINE check
        IF (.NOT. ppm_initialized) THEN
            fail("Please call ppm_init first!", exit_point=8888)
        ENDIF

        IF (ideal_load .LT. 0.0_MK) THEN
            fail("ideal_load must be >= 0.0",exit_point=8888)
        ENDIF

8888   CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE loadbal_choose_sub_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE loadbal_choose_sub_d
#endif
