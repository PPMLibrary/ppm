      minclude ppm_header(ppm_module_name="loadbal_choose_sub")

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE loadbal_choose_sub_s(ideal_load,candidate_sublist,     &
     &                                send_sublist,csublist_size,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE loadbal_choose_sub_d(ideal_load,candidate_sublist,     &
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
      REAL(MK)                                    :: ideal_load
      !!! Ideal load that would balance the workload between overloaded and
      !!! underloaded processors
      INTEGER,DIMENSION(:),POINTER,INTENT(IN   )  :: candidate_sublist
      !!! List of the candidate subdomains
      INTEGER,DIMENSION(:),POINTER,INTENT(  OUT)  :: send_sublist
      !!! Final sequence of the subs to be sent
      INTEGER                    , INTENT(IN   )  :: csublist_size
      INTEGER                    , INTENT(  OUT)  :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER                                 :: isub,i,j,max_loc,counter
      INTEGER,DIMENSION(1)                    :: temp
      REAL(MK)                                :: uber_load
      REAL(MK),DIMENSION(:),POINTER           :: candidate_costlist
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      start_subroutine("loadbal_choose_sub")

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
#if    __KIND == __SINGLE_PRECISION
      IF (MINVAL(ppm_loadbal_subcosts).GT.uber_load) THEN
#else
      IF (MINVAL(ppm_loadbal_subcostd).GT.uber_load) THEN
#endif
        ALLOCATE(send_sublist(1),STAT=info)
            or_fail_alloc("send_sublist")
        send_sublist(1) = -1
        GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Create the candidates' cost list
      !-------------------------------------------------------------------------
      counter = 0
      ALLOCATE(candidate_costlist(1:csublist_size),STAT=info)
            or_fail_alloc("candidate_costlist")
      DO i=1,csublist_size
#if    __KIND == __SINGLE_PRECISION
         candidate_costlist(i) = ppm_loadbal_subcosts(candidate_sublist(i))
#else
         candidate_costlist(i) = ppm_loadbal_subcostd(candidate_sublist(i))
#endif
      ENDDO
      !-------------------------------------------------------------------------
      !  Conditions when to send a sub
      !  2. Find the most costly subdomain that is nearly as equal as the
      !     ideal_load
      !-------------------------------------------------------------------------
      IF (MAXVAL(candidate_costlist).LT.ideal_load) THEN

        ALLOCATE(send_sublist(1:csublist_size),STAT=info)
            or_fail_alloc("send_sublist")
        send_sublist = -1
        !-----------------------------------------------------------------------
        !  Find the sub with max cost and
        !-----------------------------------------------------------------------
        DO WHILE (MAXVAL(candidate_costlist).LT.ideal_load)
            counter = counter + 1
            !-------------------------------------------------------------------
            !  SubID with the max cost is added to the sequence
            !-------------------------------------------------------------------
            temp                  = MAXLOC(candidate_costlist)
            send_sublist(counter) = temp(1)
            !-------------------------------------------------------------------
            !  Update the "remaining" of the ideal load to be sent
            !-------------------------------------------------------------------
            ideal_load = ideal_load - MAXVAL(candidate_costlist)
            !-------------------------------------------------------------------
            !  Updates candidates such that the max becomes min.
            !-------------------------------------------------------------------
            candidate_costlist(MAXLOC(candidate_costlist)) = 0._MK
        ENDDO
      !-------------------------------------------------------------------------
      !  Conditions when to send a sub
      !  3. Send only ONE sub if even the minimum subcost is more than the
      !     ideal_load (it's intrinsicly less than uber_load since I check this
      !     before)
      !     The sub should be the one with minimum cost!
      !-------------------------------------------------------------------------
      ELSE IF (MINVAL(candidate_costlist).GT.ideal_load) THEN
        ALLOCATE(send_sublist(1),STAT=info)
            or_fail_alloc("send_sublist")
        temp            = MINLOC(candidate_costlist)
        send_sublist(1) = temp(1)
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
