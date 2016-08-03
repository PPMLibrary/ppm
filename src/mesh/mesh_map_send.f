      SUBROUTINE equi_mesh_map_send(this,info)
      !!! This routine performs the actual send/recv of the
      !!! mesh blocks and all pushed data.
      !!!
      !!! [NOTE]
      !!! The first part of the buffer contains the on processor data.
      !!!
      !!! [NOTE]
      !!! Using non-blocking communication (ISend,Irecv)
      !!! could possibly speed-up this routine since only a
      !!! Send/Recv/SendRecv is done where needed (we know a
      !!! priori how much data we are going to receive).
      !!! Skipped actions could then already loop to the next
      !!! item and would not have to wait for others.
      !!!
      !!! [TIP]
      !!! If you want to only send certain meshnodes you have to first call
      !!! ppm_map_field_push for a logical mask, then call ppm_map_field_pop
      !!! and recover the new mask. First 3 (2 for 2D meshes)
      !!! indices are mesh (i,j[,k]), last one is subid for all subs on
      !!! the local processor.

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_mpi
      USE ppm_module_mapping_typedef, ONLY : ppm_mesh_isendpatchid, &
      &   ppm_mesh_isendblksize,ppm_mesh_irecvblksize
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_t_equi_mesh) :: this
      !!!
      INTEGER, INTENT(  OUT) :: info
      !!! Return status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER(ppm_kind_int64), DIMENSION(1) :: ld
      INTEGER,                 DIMENSION(3) :: ldu
      INTEGER                               :: i,j,k,bdim
      INTEGER                               :: iopt,tag1,Ndata,msend,mrecv
      INTEGER(ppm_kind_int64)               :: ibuffer,jbuffer
      INTEGER(ppm_kind_int64)               :: offs
      INTEGER(ppm_kind_int64)               :: allsend,allrecv

      start_subroutine("mesh_map_send")

      ! warn if buffer is empty
      !TODO
      IF (ppm_buffer_set.LT.1) THEN
         IF (ppm_debug.GT.1) THEN
            fail('Buffer is empty: skipping send!',ppm_err_buffer_empt,exit_point=no,ppm_error=ppm_error_notice)
            info = 0
         ENDIF
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(nsend,ldu,iopt,info)
      or_fail_alloc("nsend")

      CALL ppm_alloc(psend,ldu,iopt,info)
      or_fail_alloc("psend")

      ldu(1) = ppm_nrecvlist
      CALL ppm_alloc(nrecv,ldu,iopt,info)
      or_fail_alloc("nrecv")

      CALL ppm_alloc(precv,ldu,iopt,info)
      or_fail_alloc("precv")

      ldu(1) = ppm_nrecvlist
      ldu(2) = ppm_buffer_set
      CALL ppm_alloc(pp,ldu,iopt,info)
      or_fail_alloc("pp")

      ldu(1) = ppm_nsendlist
      ldu(2) = ppm_buffer_set
      CALL ppm_alloc(qq,ldu,iopt,info)
      or_fail_alloc("qq")

      !-------------------------------------------------------------------------
      !  Count the size of the buffer that will not be send
      !-------------------------------------------------------------------------
      Ndata = 0
      !---------------------------------------------------------------------
      !  access mesh blocks belonging to the 1st processor
      !---------------------------------------------------------------------
      DO j=ppm_psendbuffer(1),ppm_psendbuffer(2)-1
         !------------------------------------------------------------------
         !  Get the number of mesh points in this block
         !------------------------------------------------------------------
         Ndata = Ndata + PRODUCT(ppm_mesh_isendblksize(1:ppm_dim,j))
      ENDDO

      k = SUM(ppm_buffer_dim(1:ppm_buffer_set))*Ndata

      !-------------------------------------------------------------------------
      !  Initialize the buffer counters
      !-------------------------------------------------------------------------
      ppm_nrecvbuffer = INT(k,ppm_kind_int64)
      nsend(1)        = k
      nrecv(1)        = k
      psend(1)        = Ndata
      precv(1)        = Ndata
      mrecv           = -1
      msend           = -1

      !-------------------------------------------------------------------------
      !  Count the size of the buffers that will be sent
      !-------------------------------------------------------------------------
      DO k=2,ppm_nsendlist
          nsend(k) = 0
          psend(k) = 0
          Ndata    = 0

          IF (ppm_lsendlist(k)) THEN
             !---------------------------------------------------------------------
             !  Number of mesh points to be sent off to the k-th processor in
             !  the sendlist
             !---------------------------------------------------------------------
             DO i=ppm_psendbuffer(k),ppm_psendbuffer(k+1)-1
                Ndata = Ndata + PRODUCT(ppm_mesh_isendblksize(1:ppm_dim,i))
             ENDDO

             !---------------------------------------------------------------------
             !  Store the number of mesh points in psend
             !---------------------------------------------------------------------
             psend(k) = Ndata

             !---------------------------------------------------------------------
             !  Store the size of the data to be sent
             !---------------------------------------------------------------------
             nsend(k) = SUM(ppm_buffer_dim(1:ppm_buffer_set))*Ndata
          ENDIF

          !---------------------------------------------------------------------
          !  Find the maximum buffer length (for the allocate)
          !---------------------------------------------------------------------
          msend = MAX(msend,nsend(k))
          IF (ppm_debug.GT.1) THEN
             stdout_f('(A,I9)',"msend = ",msend)
          ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  Count the size of the buffers that will be received
      !  (For meshes we know this a priori. This is different from
      !  particles)
      !-------------------------------------------------------------------------
      DO k=2,ppm_nrecvlist
         nrecv(k) = 0
         precv(k) = 0
         Ndata    = 0

         IF (ppm_lrecvlist(k)) THEN
            !---------------------------------------------------------------------
            !  Number of mesh points to be received from the k-th processor in
            !  the recvlist
            !---------------------------------------------------------------------
            DO i=ppm_precvbuffer(k),ppm_precvbuffer(k+1)-1
               Ndata = Ndata + PRODUCT(ppm_mesh_irecvblksize(1:ppm_dim,i))
            ENDDO

            !---------------------------------------------------------------------
            !  Store the number of mesh points in precv
            !---------------------------------------------------------------------
            precv(k) = Ndata

            !---------------------------------------------------------------------
            !  Store the size of the data to be received
            !---------------------------------------------------------------------
            nrecv(k) = SUM(ppm_buffer_dim(1:ppm_buffer_set))*Ndata
         ENDIF

         !---------------------------------------------------------------------
         !  Find the maximum buffer length (for the allocate)
         !---------------------------------------------------------------------
         mrecv = MAX(mrecv,nrecv(k))
         IF (ppm_debug.GT.1) THEN
            stdout_f('(A,I9)',"mrecv = ",mrecv)
         ENDIF

         !---------------------------------------------------------------------
         !  Increment the total receive buffer count
         !---------------------------------------------------------------------
         ppm_nrecvbuffer = ppm_nrecvbuffer + INT(nrecv(k),ppm_kind_int64)

         IF (ppm_debug.GT.1) THEN
            stdout_f('(A,I9)',"ppm_nrecvbuffer = ",ppm_nrecvbuffer)
         ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  Allocate the memory for the copy of the particle buffer
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow
      ld(1) = ppm_nrecvbuffer
      IF (ppm_kind.EQ.ppm_kind_double) THEN
         CALL ppm_alloc(ppm_recvbufferd,ld,iopt,info)
      ELSE
         CALL ppm_alloc(ppm_recvbuffers,ld,iopt,info)
      ENDIF
      or_fail_alloc("global receive buffer PPM_RECVBUFFER")

      !-------------------------------------------------------------------------
      !  Allocate memory for the smaller send and receive buffer
      !-------------------------------------------------------------------------
      ! only allocate if there actually is anything to be sent/recvd with
      ! other processors. Otherwise mrecv and msend would both still be -1
      ! (as initialized above) and the alloc would throw a FATAL. This was
      ! Bug ID 000012.
      IF (ppm_nrecvlist.GT.1) THEN
         iopt   = ppm_param_alloc_grow
         ldu(1) = MAX(mrecv,1)
         IF (ppm_kind.EQ.ppm_kind_double) THEN
            CALL ppm_alloc(recvd,ldu,iopt,info)
         ELSE
            CALL ppm_alloc(recvs,ldu,iopt,info)
         ENDIF
         or_fail_alloc("local receive buffer recv")
      ENDIF

      IF (ppm_nsendlist.GT.1) THEN
         iopt   = ppm_param_alloc_grow
         ldu(1) = MAX(msend,1)
         IF (ppm_kind.EQ.ppm_kind_double) THEN
            CALL ppm_alloc(sendd,ldu,iopt,info)
         ELSE
            CALL ppm_alloc(sends,ldu,iopt,info)
         ENDIF
         or_fail_alloc("local send buffer send")
      ENDIF

      !-------------------------------------------------------------------------
      !  Sum of all mesh points that will be sent and received
      !-------------------------------------------------------------------------
      allsend = 0_ppm_kind_int64
      allrecv = 0_ppm_kind_int64
      DO i=1,ppm_nsendlist
         allsend = allsend + INT(psend(i),ppm_kind_int64)
      ENDDO
      DO i=1,ppm_nrecvlist
         allrecv = allrecv + INT(precv(i),ppm_kind_int64)
      ENDDO

      !-------------------------------------------------------------------------
      !  Compute the pointer to the position of the data in the main send
      !  buffer
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.1) THEN
         stdout_f('(A,I9)',"ppm_buffer_set = ",ppm_buffer_set)
      ENDIF
      bdim = 0
      offs = 0_ppm_kind_int64
      DO k=1,ppm_buffer_set
         offs    = offs + allsend*INT(bdim,ppm_kind_int64)
         qq(1,k) = offs + 1_ppm_kind_int64
         bdim    = ppm_buffer_dim(k)
         DO j=2,ppm_nsendlist
            qq(j,k) = qq(j-1,k) + INT(psend(j-1),ppm_kind_int64)*INT(bdim,ppm_kind_int64)
            IF (ppm_debug.GT.1) THEN
               stdout_f('(A,I9)',"qq(j,k) = ",'qq(j,k)')
            ENDIF
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Compute the pointer to the position of the data in the main receive
      !  buffer
      !-------------------------------------------------------------------------
      bdim = 0
      offs = 0_ppm_kind_int64
      DO k=1,ppm_buffer_set
         offs    = offs + allrecv*INT(bdim,ppm_kind_int64)
         pp(1,k) = offs + 1_ppm_kind_int64
         bdim    = ppm_buffer_dim(k)
         DO j=2,ppm_nrecvlist
            pp(j,k) = pp(j-1,k) + INT(precv(j-1),ppm_kind_int64)*INT(bdim,ppm_kind_int64)
            IF (ppm_debug.GT.1) THEN
               stdout_f('(A,I9)',"pp(j,k) = ",'pp(j,k)')
               stdout_f('(A,I9,A,I4)',"precv(j-1) = ",'precv(j-1)',", bdim = ",bdim)
            ENDIF
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  First copy the on processor data - which is in the first buffer
      !-------------------------------------------------------------------------
      IF (ppm_kind.EQ.ppm_kind_double) THEN
         DO k=1,ppm_buffer_set
            ibuffer = pp(1,k) - 1_ppm_kind_int64
            jbuffer = qq(1,k) - 1_ppm_kind_int64
            DO j=1,psend(1)*ppm_buffer_dim(k)
               ibuffer                  = ibuffer + 1_ppm_kind_int64
               jbuffer                  = jbuffer + 1_ppm_kind_int64
               ppm_recvbufferd(ibuffer) = ppm_sendbufferd(jbuffer)
            ENDDO
         ENDDO
      ELSE
         DO k=1,ppm_buffer_set
            ibuffer = pp(1,k) - 1_ppm_kind_int64
            jbuffer = qq(1,k) - 1_ppm_kind_int64
            DO j=1,psend(1)*ppm_buffer_dim(k)
               ibuffer                  = ibuffer + 1_ppm_kind_int64
               jbuffer                  = jbuffer + 1_ppm_kind_int64
               ppm_recvbuffers(ibuffer) = ppm_sendbuffers(jbuffer)
            ENDDO
         ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  loop over the processors in the ppm_isendlist(); skip the first entry
      !  which is the local processor
      !-------------------------------------------------------------------------
      IF (ppm_kind.EQ.ppm_kind_double) THEN
         !----------------------------------------------------------------------
         !  For each destination processor in the sendlist
         !----------------------------------------------------------------------
         nsendd_loop: DO k=2,ppm_nsendlist
            !-------------------------------------------------------------------
            !  Collect each buffer set
            !-------------------------------------------------------------------
            ibuffer = 0_ppm_kind_int64
            DO j=1,ppm_buffer_set
               jbuffer = qq(k,j) - 1_ppm_kind_int64
               !----------------------------------------------------------------
               !  Collect the data into the send buffer
               !----------------------------------------------------------------
               DO i=1,psend(k)*ppm_buffer_dim(j)
                  ibuffer        = ibuffer + 1_ppm_kind_int64
                  jbuffer        = jbuffer + 1_ppm_kind_int64
                  sendd(ibuffer) = ppm_sendbufferd(jbuffer)
               ENDDO
            ENDDO

            !-------------------------------------------------------------------
            !  Perform the actual send/recv as needed
            !-------------------------------------------------------------------
#ifdef __MPI
            ! The following IF is needed in order to skip "dummy"
            ! communication rounds where the current processor has to wait
            ! (only needed in the partial mapping).
            IF (ppm_isendlist(k).GE.0.AND.ppm_irecvlist(k).GE.0) THEN
               tag1 = k
               IF (psend(k).GT.0.AND.precv(k).GT.0) THEN
                  IF (ppm_debug.GT.1) THEN
                     stdout_f('(2(A,I3))',"SendRecv to ",'ppm_isendlist(k)'," from ",'ppm_irecvlist(k)')
                  ENDIF
                  CALL MPI_SendRecv(sendd,nsend(k),ppm_mpi_kind,  &
                  &    ppm_isendlist(k),tag1,recvd,nrecv(k),ppm_mpi_kind,  &
                  &    ppm_irecvlist(k),tag1,ppm_comm,MPI_STATUS_IGNORE,info)
                  or_fail_MPI("MPI_SendRecv")
               ELSEIF (psend(k).GT.0.AND.precv(k).EQ.0) THEN
                  IF (ppm_debug.GT.1) THEN
                     stdout_f('(A,I3)',"Send to ",'ppm_isendlist(k)')
                  ENDIF
                  CALL MPI_Send(sendd,nsend(k),ppm_mpi_kind,      &
                  &    ppm_isendlist(k),tag1,ppm_comm,info)
                  or_fail_MPI("MPI_Send")
               ELSEIF (psend(k).EQ.0.AND.precv(k).GT.0) THEN
                  IF (ppm_debug.GT.1) THEN
                     stdout_f('(A,I3)',"Recv from ",'ppm_irecvlist(k)')
                  ENDIF
                  CALL MPI_Recv(recvd,nrecv(k),ppm_mpi_kind,      &
                  &    ppm_irecvlist(k),tag1,ppm_comm,MPI_STATUS_IGNORE,info)
                  or_fail_MPI("MPI_Recv")
               ELSE
                  ! do nothing
               ENDIF
            ENDIF
#else
            recvd = sendd
#endif

            !-------------------------------------------------------------------
            !  Store the data back in the recv buffer
            !-------------------------------------------------------------------
            ibuffer = 0_ppm_kind_int64
            DO j=1,ppm_buffer_set
               jbuffer = pp(k,j) - 1_ppm_kind_int64
               DO i=1,precv(k)*ppm_buffer_dim(j)
                  ibuffer                  = ibuffer + 1_ppm_kind_int64
                  jbuffer                  = jbuffer + 1_ppm_kind_int64
                  ppm_recvbufferd(jbuffer) = recvd(ibuffer)
               ENDDO
            ENDDO
         ENDDO nsendd_loop
      ELSE
         !----------------------------------------------------------------------
         !  For each destination processor in the sendlist
         !----------------------------------------------------------------------
         nsends_loop: DO k=2,ppm_nsendlist
            !-------------------------------------------------------------------
            !  Collect each buffer set
            !-------------------------------------------------------------------
            ibuffer = 0_ppm_kind_int64
            DO j=1,ppm_buffer_set
               jbuffer = qq(k,j) - 1_ppm_kind_int64
               !----------------------------------------------------------------
               !  Collect the data into the send buffer
               !----------------------------------------------------------------
               DO i=1,psend(k)*ppm_buffer_dim(j)
                  ibuffer        = ibuffer + 1_ppm_kind_int64
                  jbuffer        = jbuffer + 1_ppm_kind_int64
                  sends(ibuffer) = ppm_sendbuffers(jbuffer)
               ENDDO
            ENDDO

            !-------------------------------------------------------------------
            !  Perform the actual send/recv as needed
            !-------------------------------------------------------------------
#ifdef __MPI
            ! The following IF is needed in order to skip "dummy"
            ! communication rounds where the current processor has to wait
            ! (only needed in the partial mapping).
            IF (ppm_isendlist(k).GE.0.AND.ppm_irecvlist(k).GE.0) THEN
               tag1 = k
               IF (psend(k).GT.0.AND.precv(k).GT.0) THEN
                  IF (ppm_debug.GT.1) THEN
                     stdout_f('(2(A,I3))',"SendRecv to ",'ppm_isendlist(k)'," from ",'ppm_irecvlist(k)')
                  ENDIF
                  CALL MPI_SendRecv(sends,nsend(k),ppm_mpi_kind,  &
                  &    ppm_isendlist(k),tag1,recvs,nrecv(k),ppm_mpi_kind,  &
                  &    ppm_irecvlist(k),tag1,ppm_comm,MPI_STATUS_IGNORE,info)
                  or_fail_MPI("MPI_SendRecv")
               ELSEIF (psend(k).GT.0.AND.precv(k).EQ.0) THEN
                  IF (ppm_debug.GT.1) THEN
                     stdout_f('(A,I3)',"Send to ",'ppm_isendlist(k)')
                  ENDIF
                  CALL MPI_Send(sends,nsend(k),ppm_mpi_kind,      &
                  &    ppm_isendlist(k),tag1,ppm_comm,info)
                  or_fail_MPI("MPI_Send")
               ELSEIF (psend(k).EQ.0.AND.precv(k).GT.0) THEN
                  IF (ppm_debug.GT.1) THEN
                     stdout_f('(A,I3)',"Recv from ",'ppm_irecvlist(k)')
                  ENDIF
                  CALL MPI_Recv(recvs,nrecv(k),ppm_mpi_kind,      &
                  &    ppm_irecvlist(k),tag1,ppm_comm,MPI_STATUS_IGNORE,info)
                  or_fail_MPI("MPI_Recv")
               ELSE
                  ! do nothing
               ENDIF
            ENDIF
#else
            recvs = sends
#endif

            !-------------------------------------------------------------------
            !  Store the data back in the recv buffer
            !-------------------------------------------------------------------
            ibuffer = 0_ppm_kind_int64
            DO j=1,ppm_buffer_set
               jbuffer = pp(k,j) - 1_ppm_kind_int64
               DO i=1,precv(k)*ppm_buffer_dim(j)
                  ibuffer                  = ibuffer + 1_ppm_kind_int64
                  jbuffer                  = jbuffer + 1_ppm_kind_int64
                  ppm_recvbuffers(jbuffer) = recvs(ibuffer)
               ENDDO
            ENDDO
         ENDDO nsends_loop
      ENDIF

      !-------------------------------------------------------------------------
      !  Deallocate the send buffer to save memory
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      IF (ppm_kind.EQ.ppm_kind_single) THEN
         CALL ppm_alloc(ppm_sendbuffers,ldu,iopt,info)
         or_fail_dealloc("ppm_sendbuffer")

         IF (ppm_nrecvlist.GT.1) THEN
            CALL ppm_alloc(recvs,ldu,iopt,info)
            or_fail_dealloc("recvs")
         ENDIF

         IF (ppm_nsendlist.GT.1) THEN
            CALL ppm_alloc(sends,ldu,iopt,info)
            or_fail_dealloc("sends")
         ENDIF
      ELSE
         CALL ppm_alloc(ppm_sendbufferd,ldu,iopt,info)
         or_fail_dealloc("ppm_sendbuffer")

         IF (ppm_nrecvlist.GT.1) THEN
            CALL ppm_alloc(recvd,ldu,iopt,info)
            or_fail_dealloc("recvd")
         ENDIF

         IF (ppm_nsendlist.GT.1) THEN
            CALL ppm_alloc(sendd,ldu,iopt,info)
            or_fail_dealloc("sendd")
         ENDIF
      ENDIF


      !-------------------------------------------------------------------------
      !  Deallocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(nsend,ldu,iopt,info)
      or_fail_dealloc("nsend")
      CALL ppm_alloc(psend,ldu,iopt,info)
      or_fail_dealloc("psend")
      CALL ppm_alloc(nrecv,ldu,iopt,info)
      or_fail_dealloc("nrecv")
      CALL ppm_alloc(precv,ldu,iopt,info)
      or_fail_dealloc("precv")
      CALL ppm_alloc(   pp,ldu,iopt,info)
      or_fail_dealloc("pp")
      CALL ppm_alloc(   qq,ldu,iopt,info)
      or_fail_dealloc("qq")

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      end_subroutine()
      RETURN
      END SUBROUTINE equi_mesh_map_send
