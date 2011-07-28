      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_map_part_send
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine performs the actual send/recv of the 
      !                 particles and all pushed data.
      !
      !  Input        : Npart   (I) the old number of particles on the processor
      !
      !  Input/output : 
      ! 
      !  output       : Mpart   (I) the new number of particles on
      !                             processor after the send/recv
      !                 info    (I) return status. 0 upon success.
      !
      !  Remarks      : The first part of the buffer contains the on processor
      !                 data. The packing could be performed more efficiently.
      !
      !                 The two send/recv of the integers should be combined
      !                 into one send receive
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_map_part_send.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:56  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.26  2006/09/28 15:23:29  ivos
      !  The bug fix of Jens that gives correct Mpart upon exit of
      !  map_part_ghost_get.
      !
      !  Revision 1.25  2006/02/03 09:41:27  ivos
      !  Added the PRELIMINARY ghost_put functionality. Still needs clean-up,
      !  but should work.
      !
      !  Revision 1.24  2004/11/11 15:26:18  ivos
      !  Moved allocatable work arrays to module.
      !
      !  Revision 1.23  2004/10/14 09:26:27  davidch
      !  send and receive buffer are now never smaller than 1
      !
      !  Revision 1.22  2004/10/01 16:33:37  ivos
      !  cosmetics.
      !
      !  Revision 1.21  2004/10/01 16:09:08  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.20  2004/09/21 15:20:34  hiebers
      !  erased quotation mark in revision comment
      !
      !  Revision 1.19  2004/09/21 14:38:27  ivos
      !  Added deallocate statements for ppm_sendbuffer and ppm_recvbuffer
      !  in order to save memory (Simone suggestion).
      !
      !  Revision 1.18  2004/07/26 07:42:47  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.17  2004/05/28 11:45:02  walther
      !  Added Npart to the argument list and support for particle ghosts.
      !
      !  Revision 1.16  2004/02/24 09:14:04  gonnetp
      !  send and recieve arrays must be at least length 1 (serial).
      !
      !  Revision 1.15  2004/02/20 15:05:39  walther
      !  Major bug fix in the pointer arrays: qq() and pp() - before this fix, the
      !  push and pop would only work for the xp.
      !
      !  Revision 1.14  2004/02/20 10:26:44  walther
      !  Bug fix: the sends, recvs and sendd, recvd buffers were declared with the
      !  wrong kind. Removed old commented code and aligned some MPI calls.
      !
      !  Revision 1.13  2004/02/11 17:01:41  ivos
      !  Bugfix: psend and precv replaced by psend(k), precv(k) in the debug
      !  write statements. mesg string used to overflow otherwise...
      !
      !  Revision 1.12  2004/01/28 12:26:02  ivos
      !  Bugfix: resolved bug ID 0000012 (attempt to allocate array of negative
      !  size if there was no communication needed).
      !
      !  Revision 1.11  2004/01/27 14:52:40  ivos
      !  Fixed some misplaced continuation characters.
      !
      !  Revision 1.10  2004/01/23 17:24:17  ivos
      !  Now includes ppm_define.h for the cpp defines and no longer sets them
      !  in the Makefile.
      !
      !  Revision 1.9  2004/01/23 11:31:23  ivos
      !  Cleanup: (1) updated headers, (2) inserted ppm_error and ppm_write,
      !  (3) added argument checking, (4) added checks after every alloc.
      !
      !  Revision 1.8  2003/12/17 14:00:21  ivos
      !  bug fix: psend and precv explicitly need to be set to 0 if the round 
      !  is skipped.
      !
      !  Revision 1.7  2003/12/17 10:22:11  ivos
      !  Removed error of uninitialized MPI tag in SENDRECV (caused newer 
      !  versions of MPI to crash with a broken pipe exception). 
      !  Removed tag3 and tag4.
      !
      !  Revision 1.6  2003/12/16 15:22:56  ivos
      !  bug fix: all lists are now associated even if we only have 1 processor.
      !  This is needed for pop not work if only one CPU is used.
      !
      !  Revision 1.5  2003/12/16 12:21:17  ivos
      !  replaced monster IF with GOTO.
      !
      !  Revision 1.4  2003/12/12 17:17:51  hiebers
      !  Bugfix: assigned ppm_nrevcbuffer
      !
      !  Revision 1.3  2003/12/11 17:32:16  hiebers
      !  enabled routine for serial processing
      !
      !  Revision 1.2  2003/12/09 09:54:23  ivos
      !  Added IF statements to implement ppm_icommseq() type of protocols with
      !  ranks <0 marking "wait rounds".
      !
      !  Revision 1.1.1.1  2003/11/17 15:13:45  walther
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_map_part_send(Npart,Mpart,info)

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif 
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: Npart
      INTEGER                 , INTENT(  OUT) :: Mpart,info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3) :: ldl,ldu
      INTEGER               :: i,j,k,idom,nbuffer,ibuffer,jbuffer,bdim
      INTEGER               :: iopt,count,tag1,tag2,qpart,msend,mrecv
      INTEGER               :: npart_send,npart_recv
      CHARACTER(ppm_char)   :: mesg
      REAL(ppm_kind_double) :: t0
#ifdef __MPI
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
#endif
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_part_send',t0,info)

      !-------------------------------------------------------------------------
      !  Allocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit 
      ldu(1) = MAX(ppm_nsendlist,1)
      CALL ppm_alloc(nsend,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_send',     &
     &        'send counter NSEND',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(nrecv,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_send',     &
     &        'receive counter NRECV',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(psend,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_send',     &
     &        'particle send counter PSEND',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(precv,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_send',     &
     &        'particle receive counter PRECV',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = MAX(ppm_nsendlist,1)
      ldu(2) = ppm_buffer_set 
      CALL ppm_alloc(pp,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_send',     &
     &        'work buffer PP',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(qq,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_send',     &
     &        'work buffer QQ',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Count the size of the buffer that will NOT be send
      !-------------------------------------------------------------------------
      ibuffer = 0
      qpart   = (ppm_psendbuffer(2) - ppm_psendbuffer(1))
      DO j=1,ppm_buffer_set
         bdim      = ppm_buffer_dim(j)
         ibuffer  = ibuffer  + bdim*qpart
      ENDDO 

      !-------------------------------------------------------------------------
      !  Initialize the counter for the total set of new particles 
      !JHW 20060928     IF (ppm_map_type.EQ.ppm_param_map_ghost_get) THEN
      !-------------------------------------------------------------------------
      IF (ppm_map_type.EQ.ppm_param_map_ghost_get.OR. &
          ppm_map_type.EQ.ppm_param_map_ghost_put) THEN
         Mpart = qpart + Npart
      ELSE
         Mpart = qpart
      ENDIF

      !-------------------------------------------------------------------------
      !  Count the size of the buffer that WILL be sent
      !-------------------------------------------------------------------------
      ppm_nrecvbuffer = ibuffer
      nsend(1)        = ibuffer
      nrecv(1)        = ibuffer
      psend(1)        = qpart
      precv(1)        = qpart
      mrecv           = -1
      msend           = -1

      DO k=2,ppm_nsendlist
         nsend(k) = 0
         psend(k) = 0

         !----------------------------------------------------------------------
         !  The number of particles send off to the k-th processor in the 
         !  sendlist
         !----------------------------------------------------------------------
         qpart    = (ppm_psendbuffer(k+1) - ppm_psendbuffer(k))

         !----------------------------------------------------------------------
         !  Store the number of particles to send
         !----------------------------------------------------------------------
         psend(k) = qpart ! store the number of particles to send 

         !----------------------------------------------------------------------
         !  Store the size of the data to be send
         !----------------------------------------------------------------------
         DO j=1,ppm_buffer_set
            nsend(k) = nsend(k) + ppm_buffer_dim(j)*qpart
         ENDDO 

#ifdef __MPI
         !----------------------------------------------------------------------
         !  Make a send/recv of the number of particles and data size to that 
         !  has to be send/recv
         !----------------------------------------------------------------------
         ! The following IF is needed in order to skip "dummy"
         ! communication rounds where the current processor has to wait
         ! (only needed in the partial mapping).
         IF (ppm_isendlist(k) .GE. 0 .AND. ppm_irecvlist(k) .GE. 0) THEN
             tag1 = 100
             tag2 = 200
             IF (ppm_debug .GT. 1) THEN
                 WRITE(mesg,'(A,I5,2(A,I9))') 'sending to ',   &
     &               ppm_isendlist(k),', nsend=',nsend(k),', psend=',psend(k)
                 CALL ppm_write(ppm_rank,'ppm_map_part_send',mesg,info)
             ENDIF
             CALL MPI_SendRecv(nsend(k),1,MPI_INTEGER,ppm_isendlist(k),tag1, &
     &                         nrecv(k),1,MPI_INTEGER,ppm_irecvlist(k),tag1, &
     &                         ppm_comm,status,info)
             CALL MPI_SendRecv(psend(k),1,MPI_INTEGER,ppm_isendlist(k),tag2, &
     &                         precv(k),1,MPI_INTEGER,ppm_irecvlist(k),tag2, &
     &                         ppm_comm,status,info)

             IF (ppm_debug .GT. 1) THEN
                 WRITE(mesg,'(A,I5,2(A,I9))') 'received from ',   &
     &               ppm_irecvlist(k),', nrecv=',nrecv(k),', precv=',precv(k)
                 CALL ppm_write(ppm_rank,'ppm_map_part_send',mesg,info)
             ENDIF
         ELSE
             ! skip this round, i.e. neither send nor receive any
             ! particles.
             nrecv(k) = 0
             precv(k) = 0
         ENDIF
#endif
         !----------------------------------------------------------------------
         !  Find the required (maximum) size of the send/recv buffers
         !----------------------------------------------------------------------
         mrecv = MAX(mrecv,nrecv(k))
         msend = MAX(msend,nsend(k))
         IF (ppm_debug .GT. 1) THEN
             WRITE(mesg,'(2(A,I9))') 'mrecv=',mrecv,', msend=',msend
             CALL ppm_write(ppm_rank,'ppm_map_part_send',mesg,info)
         ENDIF

         !----------------------------------------------------------------------
         !  Increment the total receive buffer
         !----------------------------------------------------------------------
         ppm_nrecvbuffer = ppm_nrecvbuffer + nrecv(k)

         !----------------------------------------------------------------------
         !  Increment the total number of particle to receive
         !----------------------------------------------------------------------
         Mpart           = Mpart           + precv(k)

         IF (ppm_debug .GT. 1) THEN
             WRITE(mesg,'(A,I9)') 'ppm_nrecvbuffer=',ppm_nrecvbuffer
             CALL ppm_write(ppm_rank,'ppm_map_part_send',mesg,info)
             WRITE(mesg,'(A,I9)') 'Mpart=',Mpart
             CALL ppm_write(ppm_rank,'ppm_map_part_send',mesg,info)
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
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_send',     &
     &        'global receive buffer PPM_RECVBUFFER',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for the smaller send and receive buffer
      !-------------------------------------------------------------------------
      ! only allocate if there actually is anything to be sent/recvd with
      ! other processors. Otherwise mrecv and msend would both still be -1
      ! (as initialized above) and the alloc would throw a FATAL. This was
      ! Bug ID 000012.
      IF (ppm_nsendlist .GT. 1) THEN
          iopt   = ppm_param_alloc_grow
          ldu(1) = MAX(mrecv,1)
          IF (ppm_kind.EQ.ppm_kind_double) THEN
             CALL ppm_alloc(recvd,ldu,iopt,info)
          ELSE
             CALL ppm_alloc(recvs,ldu,iopt,info)
          ENDIF 
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_map_part_send',     &
     &            'local receive buffer RECV',__LINE__,info)
              GOTO 9999
          ENDIF

          ldu(1) = MAX(msend,1)
          IF (ppm_kind.EQ.ppm_kind_double) THEN
             CALL ppm_alloc(sendd,ldu,iopt,info)
          ELSE
             CALL ppm_alloc(sends,ldu,iopt,info)
          ENDIF 
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_map_part_send',     &
     &            'local send buffer SEND',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Debugging print of the number of sets in the buffer
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
          WRITE(mesg,'(A,I9)') 'ppm_buffer_set=',ppm_buffer_set
          CALL ppm_write(ppm_rank,'ppm_map_part_send',mesg,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Compute the total number of particles to send 
      !  As it is now, this should be equal to Npart
      !-------------------------------------------------------------------------
      npart_send = 0
      DO j=1,ppm_nsendlist
         npart_send = npart_send + psend(j)
      ENDDO

      !-------------------------------------------------------------------------
      !  Compute the pointer to the position of the data in the main send 
      !  buffer 
      !-------------------------------------------------------------------------
      DO k=1,ppm_buffer_set
         IF (k.EQ.1) THEN
            qq(1,k) = 1
         ELSE
            qq(1,k) = qq(1,k-1) + npart_send*ppm_buffer_dim(k-1)
         ENDIF
         bdim    = ppm_buffer_dim(k)
         DO j=2,ppm_nsendlist
            qq(j,k) = qq(j-1,k) + psend(j-1)*bdim
            IF (ppm_debug .GT. 1) THEN
                WRITE(mesg,'(A,I9)') 'qq(j,k)=',qq(j,k)
                CALL ppm_write(ppm_rank,'ppm_map_part_send',mesg,info)
            ENDIF
         ENDDO
      ENDDO


      !-------------------------------------------------------------------------
      !  Compute the total number of particles to receive
      !  As it is now, this should be equal to Mpart
      !-------------------------------------------------------------------------
      npart_recv = 0
      DO j=1,ppm_nsendlist
         npart_recv = npart_recv + precv(j)
      ENDDO

      !-------------------------------------------------------------------------
      !  Compute the pointer to the position of the data in the main receive 
      !  buffer 
      !-------------------------------------------------------------------------
      DO k=1,ppm_buffer_set
         IF (k.EQ.1) THEN
            pp(1,k) = 1 
         ELSE
            pp(1,k) = pp(1,k-1) + npart_recv*ppm_buffer_dim(k-1)
         ENDIF
         bdim    = ppm_buffer_dim(k)
         DO j=2,ppm_nsendlist
            pp(j,k) = pp(j-1,k) + precv(j-1)*bdim
            IF (ppm_debug .GT. 1) THEN
                WRITE(mesg,'(A,I9)') 'pp(j,k)=',pp(j,k)
                CALL ppm_write(ppm_rank,'ppm_map_part_send',mesg,info)
                WRITE(mesg,'(A,I9,A,I4)') 'precv(j-1)=',precv(j-1),', bdim=',bdim
                CALL ppm_write(ppm_rank,'ppm_map_part_send',mesg,info)
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
         !  For each send/recv 
         !----------------------------------------------------------------------
         DO k=2,ppm_nsendlist
            !-------------------------------------------------------------------
            !  Collect each data type (xp, vp, etc)
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
            !  Perform the actual send/recv
            !-------------------------------------------------------------------
#ifdef __MPI
            ! The following IF is needed in order to skip "dummy"
            ! communication rounds where the current processor has to wait
            ! (only needed in the partial mapping).
            IF (ppm_isendlist(k) .GE. 0 .AND. ppm_irecvlist(k) .GE. 0) THEN
                tag1 = 300
                CALL MPI_SendRecv( &
     &             sendd,nsend(k),ppm_mpi_kind,ppm_isendlist(k),tag1, &
     &             recvd,nrecv(k),ppm_mpi_kind,ppm_irecvlist(k),tag1, &
     &             ppm_comm,status,info)
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
         !  For each send/recv 
         !----------------------------------------------------------------------
         DO k=2,ppm_nsendlist
            !-------------------------------------------------------------------
            !  Collect each data type (xp, vp, etc)
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
            !  Perform the actual send/recv
            !-------------------------------------------------------------------
#ifdef __MPI
            ! The following IF is needed in order to skip "dummy"
            ! communication rounds where the current processor has to wait
            ! (only needed in the partial mapping).
            IF (ppm_isendlist(k) .GE. 0 .AND. ppm_irecvlist(k) .GE. 0) THEN
                tag1 = 400
                CALL MPI_SendRecv(                                    &
     &             sends,nsend(k),ppm_mpi_kind,ppm_isendlist(k),tag1, & 
     &             recvs,nrecv(k),ppm_mpi_kind,ppm_irecvlist(k),tag1, &
     &             ppm_comm,status,info)

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
      !  before we through away the precv() data let us store it for later use:
      !  when sending ghosts back (ppm_map_part_ghost_put())
      !-------------------------------------------------------------------------
      IF (ppm_map_type.EQ.ppm_param_map_ghost_get) THEN
         ldu(1) = ppm_nsendlist + 1
         iopt   = ppm_param_alloc_grow
         CALL ppm_alloc(ppm_precvbuffer,ldu,iopt,info)
         IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
     &           'global recv buffer pointer PPM_PRECVBUFFER',__LINE__,info)
             GOTO 9999
         ENDIF
!print*,'in ppm_map_part_send'
         ppm_precvbuffer(1) = Npart + 1
         DO k=1,ppm_nsendlist
            ppm_precvbuffer(k+1) = ppm_precvbuffer(k) + precv(k)
         ENDDO

         DO k=1,ppm_nsendlist+1
!print*,'ppm_precvbuffer(1): ',k,ppm_precvbuffer(k)
         enddo
      ENDIF

      !-------------------------------------------------------------------------
      !  low level debugging
      !-------------------------------------------------------------------------
!write(mesg,'(a,i4.4)') 'recvbuf',ppm_rank
!open(10,file=mesg)
!do i=1,ppm_nrecvbuffer,2
!   write(10,*) ppm_recvbuffers(i),ppm_recvbuffers(i+1)
!enddo
!close(10)

      !-------------------------------------------------------------------------
      !  Deallocate the send buffer to save memory
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      IF (ppm_kind .EQ. ppm_kind_single) THEN
          CALL ppm_alloc(ppm_sendbuffers,ldu,iopt,info)
      ELSE
          CALL ppm_alloc(ppm_sendbufferd,ldu,iopt,info)
      ENDIF
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_send',     &
     &        'send buffer PPM_SENDBUFFER',__LINE__,info)
      ENDIF
    
      !-------------------------------------------------------------------------
      !  Deallocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(nsend,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_part_send',     &
     &        'send counter NSEND',__LINE__,info)
      ENDIF
      CALL ppm_alloc(nrecv,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_part_send',     &
     &        'receive counter NRECV',__LINE__,info)
      ENDIF
      CALL ppm_alloc(psend,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_part_send',     &
     &        'particle send counter PSEND',__LINE__,info)
      ENDIF
      CALL ppm_alloc(precv,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_part_send',     &
     &        'particle receive counter PRECV',__LINE__,info)
      ENDIF
      CALL ppm_alloc(   pp,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_part_send',     &
     &        'work array PP',__LINE__,info)
      ENDIF
      CALL ppm_alloc(   qq,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_part_send',     &
     &        'work array QQ',__LINE__,info)
      ENDIF
      CALL ppm_alloc(recvd,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_part_send',     &
     &        'local receive buffer RECVD',__LINE__,info)
      ENDIF
      CALL ppm_alloc(recvs,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_part_send',     &
     &        'local receive buffer RECVS',__LINE__,info)
      ENDIF
      CALL ppm_alloc(sendd,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_part_send',     &
     &        'local send buffer SENDD',__LINE__,info)
      ENDIF
      CALL ppm_alloc(sends,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_part_send',     &
     &        'local send buffer SENDS',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_part_send',t0,info)
      RETURN
      END SUBROUTINE ppm_map_part_send
