      minclude ppm_header(ppm_module_name="loadbal_recvsub")

      SUBROUTINE loadbal_recvsub(topoid,prec,info)
      !!! This routine is called by the underloaded processor to receive a sub
      !!! from an overlaoded neighbor
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


      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: topoid
      !!! Topology ID
      INTEGER                 , INTENT(IN   ) :: prec
      !!! Precision. One of:
      !!! * ppm_kind_single
      !!! * ppm_kind_double
      INTEGER                 , INTENT(INOUT) :: info
      !!! return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER                            :: i,isub,jsub,iproc,nprocs_of_this_sub
      INTEGER                            :: tag1,nneighsubs,iopt
      INTEGER  , DIMENSION(:  ), POINTER :: ineighsubs => NULL()
      INTEGER                            :: new_nsub
      INTEGER,DIMENSION(2)               :: ldu2
      INTEGER,DIMENSION(1)               :: ldu1
      REAL(ppm_kind_single)              :: subcosts
      REAL(ppm_kind_double)              :: subcostd
      INTEGER                            :: numelm2=2
      REAL(ppm_kind_single),DIMENSION(2) :: min_subs2,max_subs2
      REAL(ppm_kind_double),DIMENSION(2) :: min_subd2,max_subd2
      INTEGER                            :: numelm3=3
      REAL(ppm_kind_single),DIMENSION(3) :: min_subs3,max_subs3
      REAL(ppm_kind_double),DIMENSION(3) :: min_subd3,max_subd3
      TYPE(ppm_t_topo),POINTER           :: topo => NULL()
#ifdef __MPI
      INTEGER                            :: MPTYPE
      INTEGER, DIMENSION(MPI_STATUS_SIZE):: status
#endif
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      start_subroutine("loadbal_recvsub")


      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

#ifdef __MPI
      !-------------------------------------------------------------------------
      !  Define MPI data type
      !-------------------------------------------------------------------------
      IF (prec .EQ. ppm_kind_single) THEN
        MPTYPE = MPI_REAL
      ELSE IF(prec .EQ. ppm_kind_double) THEN
        MPTYPE = MPI_DOUBLE_PRECISION
      ENDIF
#endif
      !-------------------------------------------------------------------------
      !  Get the topology first
      !-------------------------------------------------------------------------
      topo   => ppm_topo(topoid)%t

      !-------------------------------------------------------------------------
      !  We need to have more than 1 processor so that DLB makes sense
      !-------------------------------------------------------------------------
      IF (ppm_nproc.EQ.1) THEN
        stdout("WARNING: you are using only one process(or). DLB will quit!")
        GOTO 9999
      ENDIF
      !-------------------------------------------------------------------------
      !  Receive subdomain info from the most underloaded processor
      !-------------------------------------------------------------------------
      IF (ppm_dim .EQ. 2) THEN
          IF (prec .EQ. ppm_kind_single) THEN
            tag1 = 400
            CALL MPI_Recv(min_subs2,numelm2,MPTYPE,ppm_loadbal_sendrank,tag1,&
     &              ppm_comm,status,info)
            or_fail("min_sub recv failed!")
            tag1 = 500
            CALL MPI_Recv(max_subs2,numelm2,MPTYPE,ppm_loadbal_sendrank,tag1,&
     &              ppm_comm,status,info)
            or_fail("max_sub recv failed!")
            tag1 = 600
            CALL MPI_Recv(subcosts,1,MPTYPE,ppm_loadbal_sendrank,tag1,&
     &              ppm_comm,status,info)
            or_fail("subcost recv failed!")

          ELSE IF(prec .EQ. ppm_kind_double) THEN
            tag1 = 400
            CALL MPI_Recv(min_subd2,numelm2,MPTYPE,ppm_loadbal_sendrank,tag1,&
     &              ppm_comm,status,info)
             or_fail("min_sub recv failed!")
            tag1 = 500
            CALL MPI_Recv(max_subd2,numelm2,MPTYPE,ppm_loadbal_sendrank,tag1,&
     &              ppm_comm,status,info)
             or_fail("max_sub recv failed!")
            tag1 = 600
            CALL MPI_Recv(subcostd,1,MPTYPE,ppm_loadbal_sendrank,tag1,&
     &              ppm_comm,status,info)
             or_fail("subcost recv failed!")

          ENDIF
      ELSE IF (ppm_dim .EQ. 3) THEN
          IF (prec .EQ. ppm_kind_single) THEN
            tag1 = 400
            CALL MPI_Recv(min_subs3,numelm3,MPTYPE,ppm_loadbal_sendrank,tag1,&
     &              ppm_comm,status,info)
            or_fail("min_sub recv failed!")
            tag1 = 500
            CALL MPI_Recv(max_subs3,numelm3,MPTYPE,ppm_loadbal_sendrank,tag1,&
     &              ppm_comm,status,info)
            or_fail("max_sub recv failed!")
            tag1 = 600
            CALL MPI_Recv(subcosts,1,MPTYPE,ppm_loadbal_sendrank,tag1,&
     &              ppm_comm,status,info)
            or_fail("subcost recv failed!")

          ELSE IF(prec .EQ. ppm_kind_double) THEN
            tag1 = 400
            CALL MPI_Recv(min_subd3,numelm3,MPTYPE,ppm_loadbal_sendrank,tag1,&
     &              ppm_comm,status,info)
             or_fail("min_sub recv failed!")
            tag1 = 500
            CALL MPI_Recv(max_subd3,numelm3,MPTYPE,ppm_loadbal_sendrank,tag1,&
     &              ppm_comm,status,info)
             or_fail("max_sub recv failed!")
            tag1 = 600
            CALL MPI_Recv(subcostd,1,MPTYPE,ppm_loadbal_sendrank,tag1,&
     &              ppm_comm,status,info)
             or_fail("subcost recv failed!")
          ENDIF

      ENDIF
      tag1 = 700
      CALL MPI_Recv(nneighsubs,1,MPI_INTEGER,ppm_loadbal_sendrank,tag1,ppm_comm,&
     &              status,info)
        or_fail("nneighsubs recv failed!")
      !-------------------------------------------------------------------------
      !  Prepare the local array for receiving ineighsubs
      !-------------------------------------------------------------------------
      ALLOCATE(ineighsubs(1:nneighsubs),STAT=info)
        or_fail_alloc("ineighsubs")
      tag1 = 800
      CALL MPI_Recv(ineighsubs,nneighsubs,MPI_INTEGER,ppm_loadbal_sendrank,  &
     &              tag1,ppm_comm,status,info)
        or_fail("ineighsubs recv failed!")
      stdout("subdomain info received")
      !-------------------------------------------------------------------------
      !  Add the received subdomain by expanding the array by one element
      !-------------------------------------------------------------------------
      new_nsub= topo%nsublist + 1
      iopt    = ppm_param_alloc_grow_preserve
      ldu2(1) = ppm_dim
      ldu2(2) = new_nsub
      ldu1(1) = new_nsub
      IF(prec .EQ. ppm_kind_single) THEN
          CALL ppm_alloc(topo%min_subs,ldu2,iopt,info)
            or_fail_alloc("min_subs re-alloc failed!")
          CALL ppm_alloc(topo%max_subs,ldu2,iopt,info)
            or_fail_alloc("max_subs re-alloc failed!")
          CALL ppm_alloc(topo%sub_costs,ldu1,iopt,info)
            or_fail_alloc("sub_costs re-alloc failed!")
      ELSE IF(prec .EQ. ppm_kind_double) THEN
          CALL ppm_alloc(topo%min_subd,ldu2,iopt,info)
            or_fail_alloc("min_subd re-alloc failed!")
          CALL ppm_alloc(topo%max_subd,ldu2,iopt,info)
            or_fail_alloc("max_subd re-alloc failed!")
          CALL ppm_alloc(topo%sub_costd,ldu1,iopt,info)
            or_fail_alloc("sub_costd re-alloc failed!")
      ENDIF
      CALL ppm_alloc(topo%nneighsubs,ldu1,iopt,info)
        or_fail_alloc("nneighsubs re-alloc failed!")
      CALL ppm_alloc(topo%isublist,ldu1,iopt,info)
        or_fail_alloc("isublist re-alloc failed!")
      CALL ppm_alloc(topo%ineighsubs,ldu2,iopt,info)
        or_fail_alloc("ineighsubs re-alloc failed!")

      !-------------------------------------------------------------------------
      !  Copy subdomain-related information
      !-------------------------------------------------------------------------
      IF (ppm_dim .EQ. 2) THEN
          IF(prec .EQ. ppm_kind_single) THEN
              topo%min_subs(:,new_nsub) = min_subs2
              topo%max_subs(:,new_nsub) = max_subs2
              topo%sub_costs( new_nsub) = subcosts
          ELSE IF(prec .EQ. ppm_kind_double) THEN
              topo%min_subd(:,new_nsub) = min_subd2
              topo%max_subd(:,new_nsub) = max_subd2
              topo%sub_costd( new_nsub) = subcostd
          ENDIF
      ELSE IF (ppm_dim .EQ. 3) THEN
          IF(prec .EQ. ppm_kind_single) THEN
              topo%min_subs(:,new_nsub) = min_subs3
              topo%max_subs(:,new_nsub) = max_subs3
              topo%sub_costs( new_nsub) = subcosts
          ELSE IF(prec .EQ. ppm_kind_double) THEN
              topo%min_subd(:,new_nsub) = min_subd3
              topo%max_subd(:,new_nsub) = max_subd3
              topo%sub_costd( new_nsub) = subcostd
          ENDIF
      ENDIF
      topo%nneighsubs(new_nsub) = nneighsubs
      topo%ineighsubs(:,new_nsub) = ineighsubs
      topo%nsublist = new_nsub
      topo%isublist(new_nsub) = new_nsub ! that's the new ID

      !-----------------------------------------------------------------------------
      !  Receive the particle buffer, first the size then the data
      !-----------------------------------------------------------------------------
      tag1 = 900
      CALL MPI_Recv(ppm_loadbal_subnpart,1,MPI_INTEGER,ppm_loadbal_sendrank,tag1, &
     &              ppm_comm,status,info)
        or_fail("# of particle coordinates recv failed!")

      !-----------------------------------------------------------------------------
      !  Before we send the buffer, we should also trim it so that the size will be
      !  shorter (we are allowed to send only an 1D array)
      !-----------------------------------------------------------------------------
      ldu1(1) = ppm_dim*ppm_loadbal_subnpart
      IF (prec .EQ. ppm_kind_single) THEN
        CALL ppm_alloc(ppm_loadbal_xps,ldu1,iopt,info)
            or_fail_alloc("ppm_loadbal_xps re-alloc failed!")
      ELSE IF (prec .EQ. ppm_kind_double) THEN
        CALL ppm_alloc(ppm_loadbal_xpd,ldu1,iopt,info)
            or_fail_alloc("ppm_loadbal_xpd re-alloc failed!")
      ENDIF
      tag1 = 1000
      IF (prec .EQ. ppm_kind_single) THEN
        CALL MPI_Recv(ppm_loadbal_xps,ldu1(1),MPTYPE,ppm_loadbal_sendrank,&
     &              tag1,ppm_comm,status,info)
            or_fail("ppm_loadbal_xps recv failed!")
      ELSE IF (prec .EQ. ppm_kind_double) THEN
        CALL MPI_Recv(ppm_loadbal_xpd,ldu1(1),MPTYPE,ppm_loadbal_sendrank,&
     &              tag1,ppm_comm,status,info)
            or_fail("ppm_loadbal_xpd recv failed!")
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



 8888   CONTINUE
      END SUBROUTINE check
      END SUBROUTINE loadbal_recvsub
