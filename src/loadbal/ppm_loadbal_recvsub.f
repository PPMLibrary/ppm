      minclude ppm_header(ppm_module_name="loadbal_recvsub")

      SUBROUTINE loadbal_recvsub(topoid,isub,sender,prec,info)
      !!! This routine is called by the underloaded processor to receive a sub
      !!! from an overloaded neighbor
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
      USE ppm_module_particles_typedef
      USE ppm_module_interfaces
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
      INTEGER                 , INTENT(IN   ) :: isub
      !!! GLOBAL ID of the sub to be received
      INTEGER                 , INTENT(IN   ) :: sender
      !!! MPI rank of the sending overloaded process
      INTEGER                 , INTENT(IN   ) :: prec
      !!! Precision. One of:
      !!! * ppm_kind_single
      !!! * ppm_kind_double
      INTEGER                 , INTENT(INOUT) :: info
      !!! return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER                            :: i,jsub,iproc,nprocs_of_this_sub
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
      INTEGER,DIMENSION(:),POINTER       :: list_add_parts => NULL()
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

      tag1 = 700
!      stdout("sender",sender)
      CALL MPI_Recv(nneighsubs,1,MPI_INTEGER,sender,tag1,ppm_comm,&
     &              status,info)
        or_fail("nneighsubs recv failed!")
      !-------------------------------------------------------------------------
      !  Prepare the local array for receiving ineighsubs
      !-------------------------------------------------------------------------
      ALLOCATE(ineighsubs(1:nneighsubs),STAT=info)
        or_fail_alloc("ineighsubs")
      tag1 = 800
      CALL MPI_Recv(ineighsubs,nneighsubs,MPI_INTEGER,sender,  &
     &              tag1,ppm_comm,status,info)
        or_fail("ineighsubs recv failed!")
!      stdout("subdomain info received")
      !-------------------------------------------------------------------------
      !  Add the received subdomain by expanding the array by one element
      !-------------------------------------------------------------------------
      new_nsub= topo%nsublist + 1
      iopt    = ppm_param_alloc_grow_preserve
      ldu2(1) = ppm_dim
      ldu2(2) = new_nsub
      ldu1(1) = new_nsub

      CALL ppm_alloc(topo%isublist,ldu1,iopt,info)
        or_fail_alloc("isublist re-alloc failed!")
      CALL ppm_alloc(topo%nneighsubs,ldu1,iopt,info)
        or_fail_alloc("nneighsubs re-alloc failed!")
      !-------------------------------------------------------------------------
      !  It may be the case where I have to increase the size of
      !  topo%nneighsubs
      !-------------------------------------------------------------------------
      ldu2(1) = nneighsubs !,MAXVAL(topo%nneighsubs))
      ldu2(2) = new_nsub
      CALL ppm_alloc(topo%ineighsubs,ldu2,iopt,info)
        or_fail_alloc("ineighsubs re-alloc failed!")


      topo%nneighsubs(new_nsub) = nneighsubs
      topo%ineighsubs(1:nneighsubs,new_nsub) = ineighsubs
      topo%nsublist = new_nsub
      topo%isublist(new_nsub) = isub ! that's the new ID

      stdout("Sub",'topo%isublist(new_nsub)'," received: total subs:",'topo%nsublist')
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
