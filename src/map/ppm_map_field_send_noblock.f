      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_map_field_send_noblock
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

      SUBROUTINE ppm_map_field_send_noblock(info)
      !!! This routine performs the actual send/recv of the
      !!! mesh blocks and all pushed data. This is the non-blocking version of
      !!! `ppm_map_field_send`.
      !!!
      !!! [WARNING]
      !!! This routine has not been tested, reviewed, or checked. Comments
      !!! and documentation are wrong. This routine might kill your cat or
      !!! worse.
      !!!
      !!! [NOTE]
      !!! The first part of the buffer contains the on processor data.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      integer, parameter :: mk = ppm_kind_double
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER              , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3) :: ldu
      INTEGER               :: i,j,k,ibuffer,jbuffer,bdim,offs
      INTEGER               :: iopt,tag1,Ndata,msend,mrecv,allsend,allrecv
      CHARACTER(ppm_char)   :: mesg
      REAL(ppm_kind_double) :: t0
#ifdef __MPI
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: commstat
#endif
      !-----------------------------------------------------
      !  Non-blocking
      !-----------------------------------------------------
      INTEGER               :: kbuffer, kreq, rindex
      INTEGER, DIMENSION(:), POINTER :: klist    => NULL()
      INTEGER, DIMENSION(:), POINTER :: isendoff => NULL()
      INTEGER, DIMENSION(:), POINTER :: irecvoff => NULL()
      INTEGER, DIMENSION(:), POINTER :: request  => NULL()
      ! 3 lines for new separate request arrays
      INTEGER               :: kreqrecv, kreqsend
      INTEGER, DIMENSION(:), POINTER :: sendreqs => NULL()
      INTEGER, DIMENSION(:), POINTER :: recvreqs => NULL()
      INTEGER, DIMENSION(:), POINTER :: klistsend => NULL()
      INTEGER, DIMENSION(:), POINTER :: klistrecv => NULL()
      REAL(ppm_kind_double), DIMENSION(:), POINTER :: isendd => NULL()
      REAL(ppm_kind_double), DIMENSION(:), POINTER :: irecvd => NULL()
      REAL(ppm_kind_single), DIMENSION(:), POINTER :: isends => NULL()
      REAL(ppm_kind_single), DIMENSION(:), POINTER :: irecvs => NULL()

      LOGICAL               :: rflag
      REAL(ppm_kind_double) :: nonblock_time0, nonblock_timer, nonblock_timeout

      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_field_send_noblock',t0,info)

      ! skip if the buffer is empty
      IF (ppm_buffer_set .LT. 1) THEN
        IF (ppm_debug .GT. 1) THEN
        CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',  &
     &      'Buffer is empty: skipping send!',info)
        ENDIF
        GOTO 9999
      ENDIF


      !-------------------------------------------------------------------------
      !  Allocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(nsend,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
     &        'send counter NSEND',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(psend,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
     &        'particle send counter PSEND',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = ppm_nrecvlist
      CALL ppm_alloc(nrecv,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
     &        'receive counter NRECV',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(precv,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
     &        'particle receive counter PRECV',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = ppm_nrecvlist
      ldu(2) = ppm_buffer_set
      CALL ppm_alloc(pp,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
     &        'work buffer PP',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = ppm_nsendlist
      ldu(2) = ppm_buffer_set
      CALL ppm_alloc(qq,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
     &        'work buffer QQ',__LINE__,info)
          GOTO 9999
      ENDIF

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
              CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
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
              CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
          ENDIF

          !---------------------------------------------------------------------
          !  Increment the total receive buffer count
          !---------------------------------------------------------------------
          ppm_nrecvbuffer = ppm_nrecvbuffer + nrecv(k)

          IF (ppm_debug .GT. 1) THEN
              WRITE(mesg,'(A,I9)') 'ppm_nrecvbuffer = ',ppm_nrecvbuffer
              CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
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
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
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
      IF (ppm_nrecvlist .GT. 1) THEN
          iopt   = ppm_param_alloc_grow
          ldu(1) = MAX(mrecv,1)
          IF (ppm_kind.EQ.ppm_kind_double) THEN
             CALL ppm_alloc(recvd,ldu,iopt,info)
          ELSE
             CALL ppm_alloc(recvs,ldu,iopt,info)
          ENDIF
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
     &            'local receive buffer RECV',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      IF (ppm_nsendlist .GT. 1) THEN
          ldu(1) = MAX(msend,1)
          IF (ppm_kind.EQ.ppm_kind_double) THEN
             CALL ppm_alloc(sendd,ldu,iopt,info)
          ELSE
             CALL ppm_alloc(sends,ldu,iopt,info)
          ENDIF
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
     &            'local send buffer SEND',__LINE__,info)
              GOTO 9999
          ENDIF
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
          CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
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
                CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
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
                CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
                WRITE(mesg,'(A,I9,A,I4)') 'precv(j-1)=',precv(j-1),', bdim=',bdim
                CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
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

      !-----------------------------------------------------
      !  Allocate buffers for non-blocking communication
      !-----------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = 2*ppm_nsendlist
      CALL ppm_alloc(klist,ldu,iopt,info)
      IF(info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
              &            'request array',__LINE__,info)
         GOTO 9999
      ENDIF
      iopt   = ppm_param_alloc_fit
      ldu(1) = 2*ppm_nsendlist
      CALL ppm_alloc(request,ldu,iopt,info)
      IF(info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
              &            'request array',__LINE__,info)
         GOTO 9999
      ENDIF

      ! 4 new arrays for separate request lists
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(klistsend,ldu,iopt,info)
      IF(info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
              &            'k list for send requests',__LINE__,info)
         GOTO 9999
      ENDIF
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(klistrecv,ldu,iopt,info)
      IF(info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
              &            'k list for recv requests',__LINE__,info)
         GOTO 9999
      ENDIF
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(sendreqs,ldu,iopt,info)
      IF(info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
              &            'send requests array',__LINE__,info)
         GOTO 9999
      ENDIF
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(recvreqs,ldu,iopt,info)
      IF(info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
              &            'recv requests array',__LINE__,info)
         GOTO 9999
      ENDIF

      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(isendoff,ldu,iopt,info)
      IF(info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
              &            'isendoff array',__LINE__,info)
         GOTO 9999
      ENDIF
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(irecvoff,ldu,iopt,info)
      IF(info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
              &            'irecvoff array',__LINE__,info)
         GOTO 9999
      ENDIF

      IF(ppm_kind.EQ.ppm_kind_double) THEN
         iopt   = ppm_param_alloc_fit
         ldu(1) = SUM(nsend(2:ppm_nsendlist))
         CALL ppm_alloc(isendd,ldu,iopt,info)
         IF(info.NE.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
                 &            'isendd array',__LINE__,info)
            GOTO 9999
         ENDIF
         iopt   = ppm_param_alloc_fit
         ldu(1) = SUM(nrecv(2:ppm_nsendlist))
         CALL ppm_alloc(irecvd,ldu,iopt,info)
         IF(info.NE.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
                 &            'irecvd array',__LINE__,info)
            GOTO 9999
         ENDIF
      ELSE
         iopt   = ppm_param_alloc_fit
         ldu(1) = SUM(nsend(2:ppm_nsendlist))
         CALL ppm_alloc(isends,ldu,iopt,info)
         IF(info.NE.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
                 &            'isendd array',__LINE__,info)
            GOTO 9999
         ENDIF
         iopt   = ppm_param_alloc_fit
         ldu(1) = SUM(nrecv(2:ppm_nsendlist))
         CALL ppm_alloc(irecvs,ldu,iopt,info)
         IF(info.NE.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
                 &            'irecvd array',__LINE__,info)
            GOTO 9999
         ENDIF
      END IF


      !-----------------------------------------------------
      !  loop over cpus in isendlist, skip the first guy which is the local
      !  processor, put the stuff in one (bigish) send buffer and store the
      !  index
      !-----------------------------------------------------
      IF (ppm_kind.EQ.ppm_kind_double) THEN
         kbuffer = 0
         DO k=2,ppm_nsendlist
            IF(psend(k).GT.0) THEN
                isendoff(k) = kbuffer+1
            ELSE
                isendoff(k) = -1
            END IF
            !-----------------------------------------------------
            !  collect each buffer set
            !-----------------------------------------------------
            DO j=1,ppm_buffer_set
               jbuffer = qq(k,j) - 1
               DO i=1,psend(k)*ppm_buffer_dim(j)
                  kbuffer = kbuffer + 1
                  jbuffer = jbuffer + 1
                  isendd(kbuffer) = ppm_sendbufferd(jbuffer)
               END DO
            END DO
         END DO
         kbuffer = 1
         DO k=2,ppm_nsendlist
            !-----------------------------------------------------
            !  construct irecvoff
            !-----------------------------------------------------
            IF(precv(k).GT.0) THEN
               irecvoff(k) = kbuffer
               kbuffer = kbuffer + nrecv(k)
            ELSE
               irecvoff(k) = -1
            END IF
         END DO

#ifdef __MPI
         kreqsend = 0
         kreqrecv = 0
         DO k=2,ppm_nsendlist
            !-----------------------------------------------------
            !  perform the actual send/recv as needed,
            !  skip zero-communication rounds (see ppm_map_field_send.f)
            !-----------------------------------------------------
            IF(ppm_isendlist(k).GE.0.AND.ppm_irecvlist(k).GE.0) THEN
               tag1 = k
               IF(psend(k).GT.0) THEN
                  kreqsend = kreqsend + 1
                  CALL MPI_ISend(isendd(isendoff(k)),nsend(k),ppm_mpi_kind, &
                       & ppm_isendlist(k),tag1,ppm_comm,sendreqs(kreqsend),info)
                  klistsend(kreqsend) = k
               END IF
               IF(precv(k).GT.0) THEN
                  kreqrecv = kreqrecv + 1
                  CALL MPI_IRecv(irecvd(irecvoff(k)),nrecv(k),ppm_mpi_kind, &
                       & ppm_irecvlist(k),tag1,ppm_comm,recvreqs(kreqrecv),info)
                  klistrecv(kreqrecv) = k
               END IF
            ELSE
               !-----------------------------------------------------
               ! do nothing
               !-----------------------------------------------------
            END IF
         END DO
         IF (ppm_debug .GT. 1) THEN
            WRITE(mesg,'(A,I9,I9)') 'send/recv requests= ', kreqsend, kreqrecv
            CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
         ENDIF
#else
#error not implemented for usage without MPI. take ppm_map_field_send.f instead
#endif

      ELSE
         !-----------------------------------------------------
         !  single precision case
         !-----------------------------------------------------
         kbuffer = 0
         DO k=2,ppm_nsendlist
            IF(psend(k).GT.0) THEN
                isendoff(k) = kbuffer+1
            ELSE
                isendoff(k) = -1
            END IF
            !-----------------------------------------------------
            !  collect each buffer set
            !-----------------------------------------------------
            DO j=1,ppm_buffer_set
               jbuffer = qq(k,j) - 1
               DO i=1,psend(k)*ppm_buffer_dim(j)
                  kbuffer = kbuffer + 1
                  jbuffer = jbuffer + 1
                  isends(kbuffer) = ppm_sendbuffers(jbuffer)
               END DO
               IF ((psend(k).GT.0).AND.(ppm_debug .GT. 1)) THEN
                  jbuffer = qq(k,j)
                  WRITE(mesg,'(A,I9,I9,F10.4,A,I9)') 'Initial send hash ',&
     &                  k,j,SUM(ppm_sendbuffers(jbuffer:jbuffer+psend(k)*&
     &                  ppm_buffer_dim(j)-1)),' to ',ppm_isendlist(k)
                  CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
               ENDIF
            END DO
         END DO
         kbuffer = 1
         DO k=2,ppm_nsendlist
            !-----------------------------------------------------
            !  construct irecvoff
            !-----------------------------------------------------
            IF(precv(k).GT.0) THEN
               irecvoff(k) = kbuffer
               kbuffer = kbuffer + nrecv(k)
            ELSE
               irecvoff(k) = -1
            END IF
         END DO

#ifdef __MPI
         kreqsend = 0
         kreqrecv = 0
         DO k=2,ppm_nsendlist
            !-----------------------------------------------------
            !  perform the actual send/recv as needed,
            !  skip zero-communication rounds (see ppm_map_field_send.f)
            !-----------------------------------------------------
            IF(ppm_isendlist(k).GE.0.AND.ppm_irecvlist(k).GE.0) THEN
               tag1 = k
               IF(psend(k).GT.0) THEN
                  kreqsend = kreqsend + 1
                  CALL MPI_ISend(isends(isendoff(k)),nsend(k),ppm_mpi_kind, &
                       & ppm_isendlist(k),tag1,ppm_comm,sendreqs(kreqsend),info)
                  klistsend(kreqsend) = k
               END IF
               IF(precv(k).GT.0) THEN
                  kreqrecv = kreqrecv + 1
                  CALL MPI_IRecv(irecvs(irecvoff(k)),nrecv(k),ppm_mpi_kind, &
                       & ppm_irecvlist(k),tag1,ppm_comm,recvreqs(kreqrecv),info)
                  klistrecv(kreqrecv) = k
               END IF
            ELSE
               !-----------------------------------------------------
               ! do nothing
               !-----------------------------------------------------
            END IF
            END DO
         IF (ppm_debug .GT. 1) THEN
            WRITE(mesg,'(A,I9,I9)') 'send/recv requests= ', kreqsend, kreqrecv
            CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
         ENDIF
#else
#error not implemented for usage without MPI. take ppm_map_field_send.f instead
#endif
      ENDIF

      IF(ppm_kind.EQ.ppm_kind_double) THEN
#undef __NONBLOCKDBG
#ifdef __NONBLOCKDBG
         CALL MPI_Barrier(ppm_comm,info)
         nonblock_time0 = MPI_Wtime(info)
         WRITE(*,*) ppm_rank, ': nonblocking ddebug started'
         nonblock_timeout = 20.0_mk
         DO
            CALL MPI_TestAny(kreq,request,rindex,rflag,MPI_STATUSES_IGNORE,info)
            IF (rflag.EQV..TRUE.) THEN
               IF (rindex.EQ.MPI_UNDEFINED) THEN
                  nonblock_timer = MPI_Wtime(info) - nonblock_time0
                  WRITE(*,*) ppm_rank, ': all non-blocking stuff completed in ', nonblock_timer
                  EXIT
               ELSE
                  nonblock_timer = MPI_Wtime(info) - nonblock_time0
               ENDIF
            ELSE
               nonblock_timer = MPI_Wtime(info) - nonblock_time0
               IF (nonblock_timer.GT.nonblock_timeout) THEN
                  WRITE(*,*) ppm_rank, ': non-blocking stuff timed out'
                  DO k=1,kreq
                     CALL MPI_Test(request(k),rflag,MPI_STATUS_IGNORE,info)
                     IF (.NOT.rflag) THEN
                        WRITE(*,*) ppm_rank, ': was still waiting for recv from ', ppm_irecvlist(klist(k))
                     ENDIF
                  ENDDO
                  CALL MPI_Abort(ppm_comm,info)
               ENDIF
            ENDIF
         ENDDO

#endif
         DO
            CALL MPI_WaitAny(kreqrecv,recvreqs,rindex,MPI_STATUSES_IGNORE,info)
            IF(rindex.EQ.MPI_UNDEFINED) EXIT
            k = klistrecv(rindex)
            kbuffer = irecvoff(k)-1
            IF (ppm_debug .GT. 1) THEN
                WRITE(mesg,'(A,I9,A,I9)') 'IRecv ', rindex, 'completed, k = ', k
                CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
            ENDIF
            DO j=1,ppm_buffer_set
               jbuffer = pp(k,j) - 1
               DO i=1,precv(k)*ppm_buffer_dim(j)
                  jbuffer = jbuffer + 1
                  kbuffer = kbuffer + 1
                  ppm_recvbufferd(jbuffer) = irecvd(kbuffer)
               END DO
            END DO
         END DO
         IF (ppm_debug .GT. 1) THEN
             WRITE(mesg,'(A)') 'All IRecv completed'
             CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
         ENDIF
         CALL MPI_Waitall(kreqsend,sendreqs,MPI_STATUSES_IGNORE,info)
         IF (ppm_debug .GT. 1) THEN
             WRITE(mesg,'(A)') 'All ISend completed'
             CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
         ENDIF


         !-----------------------------------------------------
         !  deallocate non-blocking things
         !-----------------------------------------------------
         iopt = ppm_param_dealloc
         CALL ppm_alloc(klist,ldu,iopt,info)
         CALL ppm_alloc(request,ldu,iopt,info)
         ! 4 new arrays
         CALL ppm_alloc(klistsend,ldu,iopt,info)
         CALL ppm_alloc(sendreqs,ldu,iopt,info)
         CALL ppm_alloc(klistrecv,ldu,iopt,info)
         CALL ppm_alloc(recvreqs,ldu,iopt,info)

         CALL ppm_alloc(isendoff,ldu,iopt,info)
         CALL ppm_alloc(irecvoff,ldu,iopt,info)
         CALL ppm_alloc(isendd,ldu,iopt,info)
         CALL ppm_alloc(irecvd,ldu,iopt,info)
      ELSE

#ifdef __NONBLOCKDBG
         CALL MPI_Barrier(ppm_comm,info)
         nonblock_time0 = MPI_Wtime(info)
         WRITE(*,*) ppm_rank, ': nonblocking ddebug started'
         nonblock_timeout = 20.0_mk
         DO
            CALL MPI_TestAny(kreq,request,rindex,rflag,MPI_STATUSES_IGNORE,info)
            IF (rflag.EQV..TRUE.) THEN
               IF (rindex.EQ.MPI_UNDEFINED) THEN
                  nonblock_timer = MPI_Wtime(info) - nonblock_time0
                  WRITE(*,*) ppm_rank, ': all non-blocking stuff completed in ', nonblock_timer
                  EXIT
               ELSE
                  nonblock_timer = MPI_Wtime(info) - nonblock_time0
               ENDIF
            ELSE
               nonblock_timer = MPI_Wtime(info) - nonblock_time0
               IF (nonblock_timer.GT.nonblock_timeout) THEN
                  WRITE(*,*) ppm_rank, ': non-blocking stuff timed out'
                  DO k=1,kreq
                     CALL MPI_Test(request(k),rflag,MPI_STATUS_IGNORE,info)
                     IF (.NOT.rflag) THEN
                        WRITE(*,*) ppm_rank, ': was still waiting for recv from ', ppm_irecvlist(klist(k))
                     ENDIF
                  ENDDO
                  CALL MPI_Abort(ppm_comm,info)
               ENDIF
            ENDIF
         ENDDO
#endif

         DO
            CALL MPI_WaitAny(kreqrecv,recvreqs,rindex,MPI_STATUSES_IGNORE,info)
            IF(rindex.EQ.MPI_UNDEFINED) EXIT
            k = klistrecv(rindex)
            kbuffer = irecvoff(k)-1
            IF (ppm_debug .GT. 1) THEN
                WRITE(mesg,'(A,I9,A,I9)') 'IRecv ', rindex, 'completed, k = ', k
                CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
            ENDIF
            DO j=1,ppm_buffer_set
               jbuffer = pp(k,j) - 1
               DO i=1,precv(k)*ppm_buffer_dim(j)
                  jbuffer = jbuffer + 1
                  kbuffer = kbuffer + 1
                  ppm_recvbuffers(jbuffer) = irecvs(kbuffer)
               END DO
               IF ((precv(k).GT.0).AND.(ppm_debug .GT. 1)) THEN
                  jbuffer = pp(k,j)
                  WRITE(mesg,'(A,I9,I9,F10.4,A,I9)') 'Final recv hash ',&
     &            k,j,SUM(ppm_recvbuffers(jbuffer:jbuffer+precv(k)* &
     &            ppm_buffer_dim(j)-1)),' from ',ppm_irecvlist(k)
                  CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
               ENDIF
            END DO
         END DO
         IF (ppm_debug .GT. 1) THEN
             WRITE(mesg,'(A)') 'All IRecv completed'
             CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
         ENDIF
         CALL MPI_Waitall(kreqsend,sendreqs,MPI_STATUSES_IGNORE,info)
         IF (ppm_debug .GT. 1) THEN
             WRITE(mesg,'(A)') 'All ISend completed'
             CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
         ENDIF
         !-----------------------------------------------------
         !  deallocate non-blocking things
         !-----------------------------------------------------
         iopt = ppm_param_dealloc
         CALL ppm_alloc(klist,ldu,iopt,info)
         CALL ppm_alloc(request,ldu,iopt,info)
         ! 4 new arrays
         CALL ppm_alloc(klistsend,ldu,iopt,info)
         CALL ppm_alloc(sendreqs,ldu,iopt,info)
         CALL ppm_alloc(klistrecv,ldu,iopt,info)
         CALL ppm_alloc(recvreqs,ldu,iopt,info)

         CALL ppm_alloc(isendoff,ldu,iopt,info)
         CALL ppm_alloc(irecvoff,ldu,iopt,info)
         CALL ppm_alloc(isends,ldu,iopt,info)
         CALL ppm_alloc(irecvs,ldu,iopt,info)
      END IF

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
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
     &        'send buffer PPM_SENDBUFFER',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Deallocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(nsend,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_send_noblock',     &
     &        'send counter NSEND',__LINE__,info)
      ENDIF
      CALL ppm_alloc(nrecv,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_send_noblock',     &
     &        'receive counter NRECV',__LINE__,info)
      ENDIF
      CALL ppm_alloc(psend,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_send_noblock',     &
     &        'particle send counter PSEND',__LINE__,info)
      ENDIF
      CALL ppm_alloc(precv,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_send_noblock',     &
     &        'particle receive counter PRECV',__LINE__,info)
      ENDIF
      CALL ppm_alloc(   pp,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_send_noblock',     &
     &        'work array PP',__LINE__,info)
      ENDIF
      CALL ppm_alloc(   qq,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_send_noblock',     &
     &        'work array QQ',__LINE__,info)
      ENDIF
      CALL ppm_alloc(recvd,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_send_noblock',     &
     &        'local receive buffer RECVD',__LINE__,info)
      ENDIF
      CALL ppm_alloc(recvs,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_send_noblock',     &
     &        'local receive buffer RECVS',__LINE__,info)
      ENDIF
      CALL ppm_alloc(sendd,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_send_noblock',     &
     &        'local send buffer SENDD',__LINE__,info)
      ENDIF
      CALL ppm_alloc(sends,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_send_noblock',     &
     &        'local send buffer SENDS',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_field_send_noblock',t0,info)
      RETURN
      END SUBROUTINE ppm_map_field_send_noblock
