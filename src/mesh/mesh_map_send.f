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
      USE ppm_module_data_mesh
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
      CLASS(ppm_t_equi_mesh)               :: this
      !!!
      INTEGER              , INTENT(  OUT) :: info
      !!! Return status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3) :: ldu
      INTEGER               :: i,j,k,ibuffer,jbuffer,bdim,offs
      INTEGER               :: iopt,tag1,Ndata,msend,mrecv,allsend,allrecv
      CHARACTER(ppm_char)   :: mesg
#ifdef __MPI
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: commstat
      INTEGER                             :: commreq
#endif
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------
      start_subroutine("mesh_map_send")

      ! warn if buffer is empty
      !TODO
      IF (ppm_buffer_set .LT. 1) THEN
        IF (ppm_debug .GT. 1) THEN
            info = ppm_error_notice
            CALL ppm_error(ppm_err_buffer_empt,caller,    &
     &          'Buffer is empty: skipping send!',__LINE__,info)
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
      IF (ppm_dim .LT. 3) THEN
          !---------------------------------------------------------------------
          !  access mesh blocks belonging to the 1st processor
          !---------------------------------------------------------------------
          DO j=ppm_psendbuffer(1),ppm_psendbuffer(2)-1
             !------------------------------------------------------------------
             !  Get the number of mesh points in this block
             !------------------------------------------------------------------
             Ndata = Ndata + (ppm_mesh_isendblksize(1,j)*    &
     &                ppm_mesh_isendblksize(2,j))
          ENDDO
      ELSE
          !---------------------------------------------------------------------
          !  access mesh blocks belonging to the 1st processor
          !---------------------------------------------------------------------
          DO j=ppm_psendbuffer(1),ppm_psendbuffer(2)-1
             !------------------------------------------------------------------
             !  Get the number of mesh points in this block
             !------------------------------------------------------------------
             Ndata = Ndata + (ppm_mesh_isendblksize(1,j)*    &
     &                ppm_mesh_isendblksize(2,j)*ppm_mesh_isendblksize(3,j))
          ENDDO
      ENDIF
      ibuffer = 0
      DO j=1,ppm_buffer_set
         bdim     = ppm_buffer_dim(j)
         ibuffer  = ibuffer  + bdim*Ndata
      ENDDO

      !-------------------------------------------------------------------------
      !  Initialize the buffer counters
      !-------------------------------------------------------------------------
      ppm_nrecvbuffer = ibuffer
      nsend(1)        = ibuffer
      nrecv(1)        = ibuffer
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

          !---------------------------------------------------------------------
          !  Number of mesh points to be sent off to the k-th processor in
          !  the sendlist
          !---------------------------------------------------------------------
          IF (ppm_dim .LT. 3) THEN
              DO i=ppm_psendbuffer(k),ppm_psendbuffer(k+1)-1
                  Ndata = Ndata + (ppm_mesh_isendblksize(1,i)*    &
     &                ppm_mesh_isendblksize(2,i))
              ENDDO
          ELSE
              DO i=ppm_psendbuffer(k),ppm_psendbuffer(k+1)-1
                  Ndata = Ndata + (ppm_mesh_isendblksize(1,i)*    &
     &                ppm_mesh_isendblksize(2,i)*ppm_mesh_isendblksize(3,i))
              ENDDO
          ENDIF

          !---------------------------------------------------------------------
          !  Store the number of mesh points in psend
          !---------------------------------------------------------------------
          psend(k) = Ndata

          !---------------------------------------------------------------------
          !  Store the size of the data to be sent
          !---------------------------------------------------------------------
          DO j=1,ppm_buffer_set
              nsend(k) = nsend(k) + (ppm_buffer_dim(j)*Ndata)
          ENDDO

          !---------------------------------------------------------------------
          !  Find the maximum buffer length (for the allocate)
          !---------------------------------------------------------------------
          msend = MAX(msend,nsend(k))
          IF (ppm_debug .GT. 1) THEN
              WRITE(mesg,'(A,I9)') 'msend = ',msend
              CALL ppm_write(ppm_rank,caller,mesg,info)
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

          !---------------------------------------------------------------------
          !  Number of mesh points to be received from the k-th processor in
          !  the recvlist
          !---------------------------------------------------------------------
          IF (ppm_dim .LT. 3) THEN
              DO i=ppm_precvbuffer(k),ppm_precvbuffer(k+1)-1
                  Ndata = Ndata + (ppm_mesh_irecvblksize(1,i)*    &
     &                ppm_mesh_irecvblksize(2,i))
              ENDDO
          ELSE
              DO i=ppm_precvbuffer(k),ppm_precvbuffer(k+1)-1
                  Ndata = Ndata + (ppm_mesh_irecvblksize(1,i)*    &
     &                ppm_mesh_irecvblksize(2,i)*ppm_mesh_irecvblksize(3,i))
              ENDDO
          ENDIF

          !---------------------------------------------------------------------
          !  Store the number of mesh points in precv
          !---------------------------------------------------------------------
          precv(k) = Ndata

          !---------------------------------------------------------------------
          !  Store the size of the data to be received
          !---------------------------------------------------------------------
          DO j=1,ppm_buffer_set
              nrecv(k) = nrecv(k) + (ppm_buffer_dim(j)*Ndata)
          ENDDO

          !---------------------------------------------------------------------
          !  Find the maximum buffer length (for the allocate)
          !---------------------------------------------------------------------
          mrecv = MAX(mrecv,nrecv(k))
          IF (ppm_debug .GT. 1) THEN
              WRITE(mesg,'(A,I9)') 'mrecv = ',mrecv
              CALL ppm_write(ppm_rank,caller,mesg,info)
          ENDIF

          !---------------------------------------------------------------------
          !  Increment the total receive buffer count
          !---------------------------------------------------------------------
          ppm_nrecvbuffer = ppm_nrecvbuffer + nrecv(k)

          IF (ppm_debug .GT. 1) THEN
              WRITE(mesg,'(A,I9)') 'ppm_nrecvbuffer = ',ppm_nrecvbuffer
              CALL ppm_write(ppm_rank,caller,mesg,info)
          ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  Allocate the memory for the copy of the particle buffer
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow
      ldu(1) = ppm_nrecvbuffer
      IF (ppm_kind.EQ.ppm_kind_double) THEN
         CALL ppm_alloc(ppm_recvbufferd,ldu,iopt,info)
      ELSE
         CALL ppm_alloc(ppm_recvbuffers,ldu,iopt,info)
      ENDIF
          or_fail_alloc("global receive buffer PPM_RECVBUFFER")

      !-------------------------------------------------------------------------
      !  Allocate memory for the smaller send and receive buffer
      !-------------------------------------------------------------------------
      ! only allocate if there actually is anything to be sent/recvd with
      ! other processors. Otherwise mrecv and msend would both still be -1
      ! (as initialized above) and the alloc would throw a FATAL. This was
      ! Bug ID 000012.
      IF (ppm_nrecvlist .GT. 1) THEN
          iopt   = ppm_param_alloc_grow
          ldu(1) = MAX(mrecv,1)
          IF (ppm_kind.EQ.ppm_kind_double) THEN
             CALL ppm_alloc(recvd,ldu,iopt,info)
          ELSE
             CALL ppm_alloc(recvs,ldu,iopt,info)
          ENDIF
              or_fail_alloc("local receive buffer recv")
      ENDIF

      IF (ppm_nsendlist .GT. 1) THEN
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
      allsend = SUM(psend(1:ppm_nsendlist))
      allrecv = SUM(precv(1:ppm_nrecvlist))

      !-------------------------------------------------------------------------
      !  Compute the pointer to the position of the data in the main send
      !  buffer
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
          WRITE(mesg,'(A,I9)') 'ppm_buffer_set=',ppm_buffer_set
          CALL ppm_write(ppm_rank,caller,mesg,info)
      ENDIF
      bdim = 0
      offs = 0
      DO k=1,ppm_buffer_set
         offs    = offs + allsend*bdim
         qq(1,k) = offs + 1
         bdim    = ppm_buffer_dim(k)
         DO j=2,ppm_nsendlist
            qq(j,k) = qq(j-1,k) + psend(j-1)*bdim
            IF (ppm_debug .GT. 1) THEN
                WRITE(mesg,'(A,I9)') 'qq(j,k)=',qq(j,k)
                CALL ppm_write(ppm_rank,caller,mesg,info)
            ENDIF
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Compute the pointer to the position of the data in the main receive
      !  buffer
      !-------------------------------------------------------------------------
      bdim = 0
      offs = 0
      DO k=1,ppm_buffer_set
         offs    = offs + allrecv*bdim
         pp(1,k) = offs + 1
         bdim    = ppm_buffer_dim(k)
         DO j=2,ppm_nrecvlist
            pp(j,k) = pp(j-1,k) + precv(j-1)*bdim
            IF (ppm_debug .GT. 1) THEN
                WRITE(mesg,'(A,I9)') 'pp(j,k)=',pp(j,k)
                CALL ppm_write(ppm_rank,caller,mesg,info)
                WRITE(mesg,'(A,I9,A,I4)') 'precv(j-1)=',precv(j-1),', bdim=',bdim
                CALL ppm_write(ppm_rank,caller,mesg,info)
            ENDIF
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  First copy the on processor data - which is in the first buffer
      !-------------------------------------------------------------------------
      IF (ppm_kind.EQ.ppm_kind_double) THEN
         DO k=1,ppm_buffer_set
            ibuffer = pp(1,k) - 1
            jbuffer = qq(1,k) - 1
            DO j=1,psend(1)*ppm_buffer_dim(k)
               ibuffer                  = ibuffer + 1
               jbuffer                  = jbuffer + 1
               ppm_recvbufferd(ibuffer) = ppm_sendbufferd(jbuffer)
            ENDDO
         ENDDO
      ELSE
         DO k=1,ppm_buffer_set
            ibuffer = pp(1,k) - 1
            jbuffer = qq(1,k) - 1
            DO j=1,psend(1)*ppm_buffer_dim(k)
               ibuffer                  = ibuffer + 1
               jbuffer                  = jbuffer + 1
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
         DO k=2,ppm_nsendlist
            !-------------------------------------------------------------------
            !  Collect each buffer set
            !-------------------------------------------------------------------
            ibuffer = 0
            DO j=1,ppm_buffer_set
               jbuffer = qq(k,j) - 1
               !----------------------------------------------------------------
               !  Collect the data into the send buffer
               !----------------------------------------------------------------
               DO i=1,psend(k)*ppm_buffer_dim(j)
                  ibuffer        = ibuffer + 1
                  jbuffer        = jbuffer + 1
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
            IF (ppm_isendlist(k) .GE. 0 .AND. ppm_irecvlist(k) .GE. 0) THEN
                tag1 = k
                IF (psend(k) .GT. 0 .AND. precv(k) .GT. 0) THEN
                    IF (ppm_debug .GT. 1) THEN
                        WRITE(mesg,'(2(A,I3))') 'SendRecv to ',ppm_isendlist(k),&
     &                      ' from ',ppm_irecvlist(k)
                        CALL ppm_write(ppm_rank,caller,mesg,info)
                    ENDIF
                    CALL MPI_Irecv(recvd,nrecv(k),ppm_mpi_kind,
                    CALL MPI_SendRecv(sendd,nsend(k),ppm_mpi_kind,  &
     &                 ppm_isendlist(k),tag1,recvd,nrecv(k),ppm_mpi_kind,  &
     &                 ppm_irecvlist(k),tag1,ppm_comm,commstat,info)
                       or_fail_MPI("MPI_SendRecv")
                ELSEIF (psend(k) .GT. 0 .AND. precv(k) .EQ. 0) THEN
                    IF (ppm_debug .GT. 1) THEN
                        WRITE(mesg,'(A,I3)') 'Send to ',ppm_isendlist(k)
                        CALL ppm_write(ppm_rank,caller,mesg,info)
                    ENDIF
                    CALL MPI_Send(sendd,nsend(k),ppm_mpi_kind,      &
     &                 ppm_isendlist(k),tag1,ppm_comm,info)
                       or_fail_MPI("MPI_Send")
                ELSEIF (psend(k) .EQ. 0 .AND. precv(k) .GT. 0) THEN
                    IF (ppm_debug .GT. 1) THEN
                        WRITE(mesg,'(A,I3)') 'Recv from ',ppm_irecvlist(k)
                        CALL ppm_write(ppm_rank,caller,mesg,info)
                    ENDIF
                    CALL MPI_Recv(recvd,nrecv(k),ppm_mpi_kind,      &
     &                 ppm_irecvlist(k),tag1,ppm_comm,commstat,info)
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
            ibuffer = 0
            DO j=1,ppm_buffer_set
               jbuffer = pp(k,j) - 1
               DO i=1,precv(k)*ppm_buffer_dim(j)
                  ibuffer                  = ibuffer + 1
                  jbuffer                  = jbuffer + 1
                  ppm_recvbufferd(jbuffer) = recvd(ibuffer)
               ENDDO
            ENDDO
         ENDDO
      ELSE
         !----------------------------------------------------------------------
         !  For each destination processor in the sendlist
         !----------------------------------------------------------------------
         DO k=2,ppm_nsendlist
            !-------------------------------------------------------------------
            !  Collect each buffer set
            !-------------------------------------------------------------------
            ibuffer = 0
            DO j=1,ppm_buffer_set
               jbuffer = qq(k,j) - 1
               !----------------------------------------------------------------
               !  Collect the data into the send buffer
               !----------------------------------------------------------------
               DO i=1,psend(k)*ppm_buffer_dim(j)
                  ibuffer        = ibuffer + 1
                  jbuffer        = jbuffer + 1
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
            IF (ppm_isendlist(k) .GE. 0 .AND. ppm_irecvlist(k) .GE. 0) THEN
                tag1 = k
                IF (psend(k) .GT. 0 .AND. precv(k) .GT. 0) THEN
                    CALL MPI_SendRecv(sends,nsend(k),ppm_mpi_kind,  &
     &                 ppm_isendlist(k),tag1,recvs,nrecv(k),ppm_mpi_kind,  &
     &                 ppm_irecvlist(k),tag1,ppm_comm,commstat,info)
                       or_fail_MPI("MPI_SendRecv")
                ELSEIF (psend(k) .GT. 0 .AND. precv(k) .EQ. 0) THEN
                    CALL MPI_Send(sends,nsend(k),ppm_mpi_kind,      &
     &                 ppm_isendlist(k),tag1,ppm_comm,info)
                       or_fail_MPI("MPI_Send")
                ELSEIF (psend(k) .EQ. 0 .AND. precv(k) .GT. 0) THEN
                    CALL MPI_Recv(recvs,nrecv(k),ppm_mpi_kind,      &
     &                 ppm_irecvlist(k),tag1,ppm_comm,commstat,info)
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
            ibuffer = 0
            DO j=1,ppm_buffer_set
               jbuffer = pp(k,j) - 1
               DO i=1,precv(k)*ppm_buffer_dim(j)
                  ibuffer                  = ibuffer + 1
                  jbuffer                  = jbuffer + 1
                  ppm_recvbuffers(jbuffer) = recvs(ibuffer)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Deallocate the send buffer to save memory
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      IF (ppm_kind .EQ. ppm_kind_single) THEN
          CALL ppm_alloc(ppm_sendbuffers,ldu,iopt,info)
      ELSE
          CALL ppm_alloc(ppm_sendbufferd,ldu,iopt,info)
      ENDIF
          or_fail_dealloc("ppm_sendbuffer")

      !-------------------------------------------------------------------------
      !  Deallocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(nsend,ldu,iopt,info)
      or_fail_dealloc("nsend")
      CALL ppm_alloc(nrecv,ldu,iopt,info)
      or_fail_dealloc("nrecv")
      CALL ppm_alloc(psend,ldu,iopt,info)
      or_fail_dealloc("psend")
      CALL ppm_alloc(precv,ldu,iopt,info)
      or_fail_dealloc("precv")
      CALL ppm_alloc(   pp,ldu,iopt,info)
      or_fail_dealloc("pp")
      CALL ppm_alloc(   qq,ldu,iopt,info)
      or_fail_dealloc("qq")
      CALL ppm_alloc(recvd,ldu,iopt,info)
      or_fail_dealloc("recvd")
      CALL ppm_alloc(recvs,ldu,iopt,info)
      or_fail_dealloc("recvs")
      CALL ppm_alloc(sendd,ldu,iopt,info)
      or_fail_dealloc("sendd")
      CALL ppm_alloc(sends,ldu,iopt,info)
      or_fail_dealloc("sends")

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      end_subroutine()
      RETURN
      END SUBROUTINE equi_mesh_map_send
