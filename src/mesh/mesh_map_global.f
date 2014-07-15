      SUBROUTINE equi_mesh_map_global(this,target_mesh,info)
      !!! This routine maps field data between two topologies using a global
      !!! mapping (i.e. every processor communicates with every other).
      !!! Source mesh must be on the current field topology. Global lists
      !!! with all mesh blocks that have to be sent and/or received are
      !!! built in this routine. `mesh_map_push`, `mesh_map_pop` and
      !!! `mesh_map_send` will use these lists.
      !!!
      !!! [NOTE]
      !!! The first part of the send/recv lists contains on-processor data.
      !!!
      !!! [CAUTION]
      !!! *Side effect*: this routine uses the *same* global send/recv buffers,
      !!! pointers and lists as the particle mapping routines (reason:
      !!! these buffers can be large => memory issues). Field and particle
      !!! mappings can therefore never overlap, but one must be completed
      !!! before the other starts.
      !!!
      !!! [NOTE]
      !!! A map_field_partial could be constructed as follows: every processor
      !!! known both `isendlist` and `irecvlist`. They could now use
      !!! `ppm_util_commopt` to find a good communication sequence and compress
      !!! the send/recv loop. Also, the loop when building the local send/recv
      !!! lists below would just need to go over the neighbors of a sub
      !!! instead of all subs.  Difficulty: the order of mesh blocks in the send
      !!! and recv buffer could not match up any more and push/pop could fail.
      !!! *check this!*
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data_mesh
      USE ppm_module_topo_typedef
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_t_equi_mesh)                 :: this
      CLASS(ppm_t_equi_mesh_), INTENT(IN   ) :: target_mesh
      !!! target Mesh
      INTEGER,                 INTENT(  OUT) :: info
      !!! Returns status, 0 upon success

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      TYPE(ppm_t_topo), POINTER :: topo
      TYPE(ppm_t_topo), POINTER :: target_topo

      REAL(ppm_kind_double) :: t0

      INTEGER, DIMENSION(ppm_dim) :: Offset
      INTEGER, DIMENSION(2)       :: ldu
      INTEGER                     :: i,j,isub,sendrank,recvrank
      INTEGER                     :: iopt,iset,ibuffer,pdim
      INTEGER                     :: nsendlist
      INTEGER                     :: nrecvlist

      CHARACTER(ppm_char) :: caller='mesh_map_global'

      LOGICAL :: valid

      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)
      pdim = ppm_dim

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         CALL check
         IF (info.NE.0) GOTO 9999
      ENDIF

      topo        => ppm_topo(this%topoid)%t
      target_topo => ppm_topo(target_mesh%topoid)%t

      IF (ppm_buffer_set.GT.0) THEN
         fail('Buffer was not empty. Possible loss of data!',ppm_err_map_incomp, &
         & exit_point=no,ppm_error=ppm_error_warning)
      ENDIF

      !-------------------------------------------------------------------------
      !  Save the map type for the subsequent calls (used for checks!)
      !-------------------------------------------------------------------------
      ppm_map_type = ppm_param_map_global

      !-------------------------------------------------------------------------
      !  Check if origin and target meshes are compatible (i.e. have the
      !  same number of grid points in the whole comput. domain)
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         IF (ANY(this%Nm.NE.target_mesh%Nm)) THEN
            stdout("Source and destination meshes are of different size",ppm_error_notice)
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for the local temporary sendlists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = topo%nsublist
      CALL ppm_alloc(isendfromsub,ldu,iopt,info)
      or_fail_alloc('local send source sub list ISENDFROMSUB')

      CALL ppm_alloc(isendtosub,ldu,iopt,info)
      or_fail_alloc('local send destination sub list ISENDTOSUB')

      ldu(1) = pdim
      ldu(2) = topo%nsublist
      CALL ppm_alloc(isendpatchid,ldu,iopt,info)
      or_fail_alloc("isendpatchid,ldu")

      CALL ppm_alloc(isendblkstart,ldu,iopt,info)
      or_fail_alloc('local send block start list ISENDBLKSTART')

      CALL ppm_alloc(isendblksize,ldu,iopt,info)
      or_fail_alloc('local send block size list ISENDBLKSIZE')

      CALL ppm_alloc(ioffset,ldu,iopt,info)
      or_fail_alloc('local send block offset list IOFFSET')

      !-------------------------------------------------------------------------
      !  Find intersecting mesh domains to be sent
      !-------------------------------------------------------------------------
      nsendlist = 0
      !TOCHECK
      Offset    = 0
      DO i=1,topo%nsublist
         isub = topo%isublist(i)
         DO j=1,target_topo%nsubs
            CALL this%block_intersect(target_mesh,isub,j,Offset, &
            &    nsendlist,isendfromsub,isendtosub,isendpatchid, &
            &    isendblkstart,isendblksize,ioffset,info)
            or_fail("block_intersect failed")
          ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Allocate memory for the local temporary receive lists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = target_topo%nsublist
      CALL ppm_alloc(irecvfromsub,ldu,iopt,info)
      or_fail_alloc('local recv source sub list IRECVFROMSUB')

      CALL ppm_alloc(irecvtosub,ldu,iopt,info)
      or_fail_alloc('local recv destination sub list IRECVTOSUB')

      ldu(1) = pdim
      ldu(2) = target_topo%nsublist
      CALL ppm_alloc(irecvpatchid,ldu,iopt,info)
      or_fail_alloc("isendpatchid,ldu")

      CALL ppm_alloc(irecvblkstart,ldu,iopt,info)
      or_fail_alloc('local recv block start list IRECVBLKSTART')

      CALL ppm_alloc(irecvblksize,ldu,iopt,info)
      or_fail_alloc('local recv block size list IRECVBLKSIZE')

      !-------------------------------------------------------------------------
      !  Find intersecting mesh domains to be received
      !-------------------------------------------------------------------------
      nrecvlist = 0
      ! loop over fromtopo first and THEN over totopo in order to get the
      ! same ordering of mesh blocks as in the sendlist. This is crucial
      ! for the push and the pop to work properly
      DO j=1,topo%nsubs
         DO i=1,target_topo%nsublist
            isub = target_topo%isublist(i)
            CALL this%block_intersect(target_mesh,j,isub,offset, &
            &    nrecvlist,irecvfromsub,irecvtosub,irecvpatchid, &
            &    irecvblkstart,irecvblksize,ioffset,info)
            or_fail("block_intersect failed")
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Allocate memory for the global mesh sendlists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = nsendlist
      CALL ppm_alloc(ppm_mesh_isendfromsub,ldu,iopt,info)
      or_fail_alloc('source send sub list PPM_MESH_ISENDFROMSUB')

      ldu(1) = pdim
      ldu(2) = nsendlist
      CALL ppm_alloc(ppm_mesh_isendblkstart,ldu,iopt,info)
      or_fail_alloc('send block start list PPM_MESH_ISENDBLKSTART')

      CALL ppm_alloc(ppm_mesh_isendblksize,ldu,iopt,info)
      or_fail_alloc('send block size list PPM_MESH_ISENDBLKSIZE')

      CALL ppm_alloc(ppm_mesh_isendpatchid,ldu,iopt,info)
      or_fail_alloc("send patch id list ppm_mesh_isendpatchid")

      !-------------------------------------------------------------------------
      !  Allocate memory for the global mesh receive lists
      !-------------------------------------------------------------------------
      ldu(1) = nrecvlist
      CALL ppm_alloc(ppm_mesh_irecvtosub,ldu,iopt,info)
      or_fail_alloc('destination recv sub list PPM_MESH_IRECVTOSUB')

      ldu(1) = pdim
      ldu(2) = nrecvlist
      CALL ppm_alloc(ppm_mesh_irecvblkstart,ldu,iopt,info)
      or_fail_alloc('recv block start list PPM_MESH_IRECVBLKSTART')

      CALL ppm_alloc(ppm_mesh_irecvblksize,ldu,iopt,info)
      or_fail_alloc('recv block size list PPM_MESH_IRECVBLKSIZE')

      CALL ppm_alloc(ppm_mesh_irecvpatchid,ldu,iopt,info)
      or_fail_alloc("recv patch id list ppm_mesh_irecvpatchid")

      !-------------------------------------------------------------------------
      !  Allocate memory for the global send/recv lists
      !-------------------------------------------------------------------------
      ldu(1) = ppm_nproc
      CALL ppm_alloc(ppm_isendlist,ldu,iopt,info)
      or_fail_alloc('global send rank list PPM_ISENDLIST')

      CALL ppm_alloc(ppm_irecvlist,ldu,iopt,info)
      or_fail_alloc('global recv rank list PPM_IRECVLIST')

      ldu(1) = ppm_nproc + 1
      CALL ppm_alloc(ppm_psendbuffer,ldu,iopt,info)
      or_fail_alloc('global send buffer pointer PPM_PSENDBUFFER')

      CALL ppm_alloc(ppm_precvbuffer,ldu,iopt,info)
      or_fail_alloc('global recv buffer pointer PPM_PRECVBUFFER')

      !-------------------------------------------------------------------------
      !  Reset the number of buffer entries
      !-------------------------------------------------------------------------
      ppm_buffer_set = 0

      !-------------------------------------------------------------------------
      !  loop over all processors; First is the processor itself
      !-------------------------------------------------------------------------
      sendrank           = ppm_rank - 1
      recvrank           = ppm_rank + 1
      ppm_psendbuffer(1) = 1
      ppm_precvbuffer(1) = 1
      ppm_nsendlist      = 0
      ppm_nrecvlist      = 0
      ppm_nsendbuffer    = 0    ! data length in buffer
      ppm_nrecvbuffer    = 0    ! data length in buffer
      iset               = 0
      ibuffer            = 0

      DO i=1,ppm_nproc

         !----------------------------------------------------------------------
         !  compute the next processor
         !----------------------------------------------------------------------
         sendrank = sendrank + 1
         IF (sendrank.GE.ppm_nproc) sendrank = sendrank - ppm_nproc
         recvrank = recvrank - 1
         IF (recvrank.LT.        0) recvrank = recvrank + ppm_nproc

         !----------------------------------------------------------------------
         !  Store the processor to which we will send to
         !----------------------------------------------------------------------
         ppm_nsendlist                = ppm_nsendlist + 1
         ppm_isendlist(ppm_nsendlist) = sendrank

         !----------------------------------------------------------------------
         !  Store the processor to which we will recv from
         !----------------------------------------------------------------------
         ppm_nrecvlist                = ppm_nrecvlist + 1
         ppm_irecvlist(ppm_nrecvlist) = recvrank

         !----------------------------------------------------------------------
         !  Find all mesh blocks that are sent and store them.
         !  To get all mesh blocks that are sent in round
         !  i=1,ppm_nsendlist:
         !     dest_rank = ppm_isendlist(i)
         !     ilo = ppm_psendbuffer(i)
         !     ihi = ppm_psendbuffer(i+1)-1
         !     DO j=ilo,ihi
         !       SEND(ppm_mesh_isendblk(:,j) of sub ppm_mesh_isendfromsub(j) to
         !         processor dest_rank to sub ppm_mesh_isendtosub(j))
         !     ENDDO
         !----------------------------------------------------------------------
         ibuffer = ppm_nsendlist + 1
         ppm_psendbuffer(ibuffer) = ppm_psendbuffer(ppm_nsendlist)
         DO j=1,nsendlist
            IF (target_topo%sub2proc(isendtosub(j)) .EQ. sendrank) THEN
               ppm_psendbuffer(ibuffer) = ppm_psendbuffer(ibuffer) + 1
               iset = ppm_psendbuffer(ibuffer) - 1
               ppm_mesh_isendfromsub(iset)         = isendfromsub(j)
               ppm_mesh_isendblkstart(1:pdim,iset) = isendblkstart(1:pdim,j)
               ppm_mesh_isendblksize(1:pdim,iset)  = isendblksize(1:pdim,j)
               ppm_mesh_isendpatchid(1:pdim,iset)  = isendpatchid(1:pdim,j)

               IF (ppm_debug .GT. 1) THEN
                  SELECT CASE (ppm_dim)
                  CASE (2)
                     stdout_f('(2(A,2I4),A,I3)'," sending ",'isendblkstart(1:2,j)', &
                     & " of size ",'isendblksize(1:2,j)'," to ",sendrank)
                  CASE (3)
                     stdout_f('(2(A,3I4),A,I3)'," sending ",'isendblkstart(1:3,j)', &
                     & " of size ",'isendblksize(1:3,j)'," to ",sendrank)
                  END SELECT
               ENDIF
            ENDIF
         ENDDO

         !----------------------------------------------------------------------
         !  Find all mesh blocks that are received and store them.
         !  To get all mesh blocks that are received in round
         !  i=1,ppm_nrecvlist:
         !     source_rank = ppm_irecvlist(i)
         !     ilo = ppm_precvbuffer(i)
         !     ihi = ppm_precvbuffer(i+1)-1
         !     DO j=ilo,ihi
         !       RECV(ppm_mesh_irecvblk(:,j) of sub ppm_mesh_irecvfromsub(j)
         !         from processor source_rank to sub ppm_mesh_irecvtosub(j))
         !     ENDDO
         !----------------------------------------------------------------------
         ibuffer = ppm_nrecvlist + 1
         ppm_precvbuffer(ibuffer) = ppm_precvbuffer(ppm_nrecvlist)
         DO j=1,nrecvlist
            IF (topo%sub2proc(irecvfromsub(j)) .EQ. recvrank) THEN
               ppm_precvbuffer(ibuffer) = ppm_precvbuffer(ibuffer) + 1
               iset = ppm_precvbuffer(ibuffer) - 1
               ppm_mesh_irecvtosub(iset)           = irecvtosub(j)
               ppm_mesh_irecvblkstart(1:pdim,iset) = irecvblkstart(1:pdim,j)
               ppm_mesh_irecvblksize(1:pdim,iset)  = irecvblksize(1:pdim,j)
               ppm_mesh_irecvpatchid(1:pdim,iset)  = irecvpatchid(1:pdim,j)

               IF (ppm_debug .GT. 1) THEN
                  SELECT CASE (ppm_dim)
                  CASE (2)
                     stdout_f('(2(A,2I4),A,I3)'," recving ",'irecvblksize(1:2,j)', &
                     & " of size ",'irecvblkstart(1:2,j)'," from ",recvrank)
                  CASE (3)
                     stdout_f('(2(A,3I4),A,I3)'," recving ",'irecvblksize(1:3,j)', &
                     & " of size ",'irecvblkstart(1:3,j)'," from ",recvrank)
                  END SELECT
               ENDIF
            ENDIF
         ENDDO

      ENDDO ! i=1,ppm_nproc

      !-------------------------------------------------------------------------
      !  Deallocate memory of the local lists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_dealloc
      CALL ppm_alloc(isendfromsub,ldu,iopt,info)
      or_fail_dealloc('local send source sub list ISENDFROMSUB')

      CALL ppm_alloc(isendtosub,ldu,iopt,info)
      or_fail_dealloc('local send destination sub list ISENDTOSUB')

      CALL ppm_alloc(isendblkstart,ldu,iopt,info)
      or_fail_dealloc('local send block start list ISENDBLKSTART')

      CALL ppm_alloc(isendblksize,ldu,iopt,info)
      or_fail_dealloc('local send block size list ISENDBLKSIZE')

      CALL ppm_alloc(isendpatchid,ldu,iopt,info)
      or_fail_dealloc('local send patch id list ISENDPATCHID')

      CALL ppm_alloc(irecvfromsub,ldu,iopt,info)
      or_fail_dealloc('local recv source sub list IRECVFROMSUB')

      CALL ppm_alloc(irecvtosub,ldu,iopt,info)
      or_fail_dealloc('local recv destination sub list IRECVTOSUB')

      CALL ppm_alloc(irecvblkstart,ldu,iopt,info)
      or_fail_dealloc('local recv block start list IRECVBLKSTART')

      CALL ppm_alloc(irecvblksize,ldu,iopt,info)
      or_fail_dealloc('local recv block size list IRECVBLKSIZE')

      CALL ppm_alloc(irecvpatchid,ldu,iopt,info)
      or_fail_dealloc('local recv patch id list IRECVPATCHID')

      CALL ppm_alloc(ioffset,ldu,iopt,info)
      or_fail_dealloc('local recv block offset list IOFFSET')

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
        USE ppm_module_check_id
        IMPLICIT NONE
        CALL ppm_check_topoid(this%topoid,valid,info)
        IF (.NOT. valid) THEN
           fail("topoid not valid",ppm_err_argument,exit_point=8888)
        ENDIF
        CALL ppm_check_topoid(target_mesh%topoid,valid,info)
        IF (.NOT. valid) THEN
           fail("target topoid not valid",ppm_err_argument,exit_point=8888)
        ENDIF
      8888 CONTINUE
      END SUBROUTINE check
      END SUBROUTINE equi_mesh_map_global
