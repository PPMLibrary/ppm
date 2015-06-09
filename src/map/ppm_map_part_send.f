      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_map_part_send
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

      SUBROUTINE ppm_map_part_send(Npart,Mpart,info)
      !!! This routine performs the actual send/recv of the particles and all
      !!! pushed data.
      !!!
      !!! [NOTE]
      !!! The first part of the buffer contains the on processor data. The
      !!! packing could be performed more efficiently.
      !!!
      !!! [NOTE]
      !!! The two send/recv of the integers should be combined into one
      !!! send/receive

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      USE ppm_module_mpi
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER, INTENT(IN   ) :: Npart
      !!! The old number of particles on the processor
      INTEGER, INTENT(  OUT) :: Mpart
      !!! The new number of particles on processor after the send/recv
      INTEGER, INTENT(  OUT) :: info
      !!! Return status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      INTEGER, DIMENSION(3)               :: ldl,ldu
      INTEGER                             :: i,j,k,idom,bdim,sbdim
      INTEGER                             :: nbuffer,ibuffer,jbuffer
      INTEGER                             :: iopt,count,tag1,qpart,msend,mrecv
      INTEGER                             :: npart_send,npart_recv

      CHARACTER(ppm_char) :: caller='ppm_map_part_send'

      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)
      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.3) THEN
         CALL check
         IF (info.NE.0) GOTO 9999
      ENDIF

      ! skip if the buffer is empty
      IF (ppm_buffer_set.LT.1) THEN
         fail('Buffer is empty: skipping send!', &
         & ppm_err_buffer_empt,ppm_error=ppm_error_notice)
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_grow
      ldu(1) = MAX(ppm_nsendlist,1)
      IF ((ppm_nsendlist.NE.old_nsendlist) .OR.  &
      &  (ppm_buffer_set.NE.old_buffer_set)) THEN
         old_nsendlist = ppm_nsendlist
         old_buffer_set = ppm_buffer_set
         CALL ppm_alloc(nsend,ldu,iopt,info)
         or_fail_alloc('send counter NSEND',ppm_error=ppm_error_fatal)

         CALL ppm_alloc(nrecv,ldu,iopt,info)
         or_fail_alloc('receive counter NRECV',ppm_error=ppm_error_fatal)

         CALL ppm_alloc(psend,ldu,iopt,info)
         or_fail_alloc('particle send counter PSEND',ppm_error=ppm_error_fatal)

         CALL ppm_alloc(precv,ldu,iopt,info)
         or_fail_alloc('particle receive counter PRECV',ppm_error=ppm_error_fatal)

         ldu(2) = ppm_buffer_set
         CALL ppm_alloc(pp,ldu,iopt,info)
         or_fail_alloc('work buffer PP',ppm_error=ppm_error_fatal)

         CALL ppm_alloc(qq,ldu,iopt,info)
         or_fail_alloc('work buffer QQ',ppm_error=ppm_error_fatal)
      ENDIF
      !-------------------------------------------------------------------------
      !  Count the total size of the buffer dimensions
      !  It is handy... simplifies many summations, avoids loops, etc...
      !-------------------------------------------------------------------------
      sbdim = SUM(ppm_buffer_dim(1:ppm_buffer_set))

      !-------------------------------------------------------------------------
      !  Count the size of the buffer that will NOT be sent
      !-------------------------------------------------------------------------
      qpart   = (ppm_psendbuffer(2) - ppm_psendbuffer(1))
      ibuffer = sbdim*qpart

      !-------------------------------------------------------------------------
      !  Initialize the counter for the total set of new particles
      !JHW 20060928     IF (ppm_map_type.EQ.ppm_param_map_ghost_get) THEN
      !-------------------------------------------------------------------------
      IF (ppm_map_type.EQ.ppm_param_map_ghost_get.OR. &
      &   ppm_map_type.EQ.ppm_param_map_ghost_put) THEN
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
         nsend(k) = sbdim*qpart

#ifdef __MPI
         !----------------------------------------------------------------------
         !  Make a send/recv of the number of particles and data size to that
         !  has to be send/recv
         !----------------------------------------------------------------------
         ! The following IF is needed in order to skip "dummy"
         ! communication rounds where the current processor has to wait
         ! (only needed in the partial mapping).
         IF (ppm_isendlist(k).GE.0.AND.ppm_irecvlist(k).GE.0) THEN
             tag1 = 100
             IF (ppm_debug .GT. 1) THEN
                stdout_f('(A,I5,2(A,I9))',"sending to ",'ppm_isendlist(k)', &
                & ", nsend=",'nsend(k)',", psend=",'psend(k)')
             ENDIF
             CALL MPI_SendRecv(psend(k),1,MPI_INTEGER,ppm_isendlist(k),tag1, &
             &    precv(k),1,MPI_INTEGER,ppm_irecvlist(k),tag1,ppm_comm,     &
             &    MPI_STATUS_IGNORE,info)

             ! Compute nrecv(k) from precv(k)
             nrecv(k) = sbdim*precv(k)

             IF (ppm_debug .GT. 1) THEN
                stdout_f('(A,I5,2(A,I9))',"received from ",'ppm_irecvlist(k)', &
                & ", nrecv=",'nrecv(k)',", precv=",'precv(k)')
             ENDIF
         ELSE
             ! skip this round, i.e. neither send nor receive any
             ! particles.
             nrecv(k) = 0
             precv(k) = 0
         ENDIF
#endif

      ENDDO
      !----------------------------------------------------------------------
      !  Find the required (maximum) size of the send/recv buffers
      !----------------------------------------------------------------------
      DO k=2,ppm_nsendlist
         IF (mrecv.LE.nrecv(k)) mrecv = nrecv(k)
      ENDDO
      DO k=2,ppm_nsendlist
         IF (msend.LE.nsend(k)) msend = nsend(k)
      ENDDO
      !----------------------------------------------------------------------
      !  Increment the total receive buffer
      !----------------------------------------------------------------------
      DO k=2,ppm_nsendlist
         ppm_nrecvbuffer = ppm_nrecvbuffer + nrecv(k)
      ENDDO
      !----------------------------------------------------------------------
      !  Increment the total number of particle to receive
      !----------------------------------------------------------------------
      DO k=2,ppm_nsendlist
         Mpart           = Mpart           + precv(k)
      ENDDO

      IF (ppm_debug .GT. 1) THEN
          stdout_f('(2(A,I9))',"mrecv=",mrecv,", msend=",msend)
          stdout_f('(A,I9)',"ppm_nrecvbuffer=",ppm_nrecvbuffer)
          stdout_f('(A,I9)',"Mpart=",Mpart)
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate the memory for the copy of the particle buffer
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow
      ldu(1) = ppm_nrecvbuffer
      SELECT CASE (ppm_kind)
      CASE (ppm_kind_double)
         CALL ppm_alloc(ppm_recvbufferd,ldu,iopt,info)
      CASE DEFAULT
         CALL ppm_alloc(ppm_recvbuffers,ldu,iopt,info)
      END SELECT
      or_fail_alloc('global receive buffer PPM_RECVBUFFER', &
      & ppm_error=ppm_error_fatal)

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
         SELECT CASE (ppm_kind)
         CASE (ppm_kind_double)
            CALL ppm_alloc(recvd,ldu,iopt,info)
            or_fail_alloc('local receive buffer RECV',ppm_error=ppm_error_fatal)

            ldu(1) = MAX(msend,1)
            CALL ppm_alloc(sendd,ldu,iopt,info)

         CASE DEFAULT
            CALL ppm_alloc(recvs,ldu,iopt,info)
            or_fail_alloc('local receive buffer RECV',ppm_error=ppm_error_fatal)

            ldu(1) = MAX(msend,1)
            CALL ppm_alloc(sends,ldu,iopt,info)

         END SELECT
         or_fail_alloc('local send buffer SEND',ppm_error=ppm_error_fatal)
      ENDIF

      !-------------------------------------------------------------------------
      !  Debugging print of the number of sets in the buffer
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
         stdout_f('(A,I9)',"ppm_buffer_set=",ppm_buffer_set)
      ENDIF

      !-------------------------------------------------------------------------
      !  Compute the total number of particles to send
      !  As it is now, this should be equal to Npart
      !-------------------------------------------------------------------------
      npart_send = SUM(psend(1:ppm_nsendlist))

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
               stdout_f('(A,I9)',"qq(j,k)=",'qq(j,k)')
            ENDIF
         ENDDO
      ENDDO


      !-------------------------------------------------------------------------
      !  Compute the total number of particles to receive
      !  As it is now, this should be equal to Mpart
      !-------------------------------------------------------------------------
      npart_recv=SUM(precv(1:ppm_nsendlist))

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
               stdout_f('(A,I9)',"pp(j,k)=",'pp(j,k)')
               stdout_f('(A,I9,A,I4)',"precv(j-1)=",'precv(j-1)',", bdim=",bdim)
            ENDIF
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  First copy the on processor data - which is in the first buffer
      !-------------------------------------------------------------------------
      SELECT CASE (ppm_kind)
      CASE (ppm_kind_double)
         DO k=1,ppm_buffer_set
            ibuffer = pp(1,k) - 1
            jbuffer = qq(1,k) - 1
            DO j=1,psend(1)*ppm_buffer_dim(k)
               ibuffer                  = ibuffer + 1
               jbuffer                  = jbuffer + 1
               ppm_recvbufferd(ibuffer) = ppm_sendbufferd(jbuffer)
            ENDDO
         ENDDO
      CASE DEFAULT
         DO k=1,ppm_buffer_set
            ibuffer = pp(1,k) - 1
            jbuffer = qq(1,k) - 1
            DO j=1,psend(1)*ppm_buffer_dim(k)
               ibuffer                  = ibuffer + 1
               jbuffer                  = jbuffer + 1
               ppm_recvbuffers(ibuffer) = ppm_sendbuffers(jbuffer)
            ENDDO
         ENDDO
      END SELECT

      !-------------------------------------------------------------------------
      !  loop over the processors in the ppm_isendlist(); skip the first entry
      !  which is the local processor
      !-------------------------------------------------------------------------
      SELECT CASE (ppm_kind)
      CASE (ppm_kind_double)
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
                &    sendd,nsend(k),ppm_mpi_kind,ppm_isendlist(k),tag1, &
                &    recvd,nrecv(k),ppm_mpi_kind,ppm_irecvlist(k),tag1, &
                &    ppm_comm,MPI_STATUS_IGNORE,info)
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
      CASE DEFAULT
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
                CALL MPI_SendRecv(                                      &
                &    sends,nsend(k),ppm_mpi_kind,ppm_isendlist(k),tag1, &
                &    recvs,nrecv(k),ppm_mpi_kind,ppm_irecvlist(k),tag1, &
                &    ppm_comm,MPI_STATUS_IGNORE,info)
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
      END SELECT

      !-------------------------------------------------------------------------
      !  before we through away the precv() data let us store it for later use:
      !  when sending ghosts back (ppm_map_part_ghost_put())
      !-------------------------------------------------------------------------
      IF (ppm_map_type.EQ.ppm_param_map_ghost_get) THEN
         ldu(1) = ppm_nsendlist + 1
         iopt   = ppm_param_alloc_grow
         CALL ppm_alloc(ppm_precvbuffer,ldu,iopt,info)
         or_fail_alloc('global recv buffer pointer PPM_PRECVBUFFER', &
         & ppm_error=ppm_error_fatal)

         ppm_precvbuffer(1) = Npart + 1
         DO k=1,ppm_nsendlist
            ppm_precvbuffer(k+1) = ppm_precvbuffer(k) + precv(k)
         ENDDO
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
      SELECT CASE (ppm_kind)
      CASE (ppm_kind_single)
          CALL ppm_alloc(ppm_sendbuffers,ldu,iopt,info)
      CASE DEFAULT
          CALL ppm_alloc(ppm_sendbufferd,ldu,iopt,info)
      END SELECT
      or_fail_dealloc('send buffer PPM_SENDBUFFER',exit_point=no)


      !-------------------------------------------------------------------------
      !  Deallocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
!      CALL ppm_alloc(nsend,ldu,iopt,info)
!      IF (info .NE. 0) THEN
!          info = ppm_error_error
!          CALL ppm_error(ppm_err_dealloc,caller,     &
!     &        'send counter NSEND',__LINE__,info)
!      ENDIF
!      CALL ppm_alloc(nrecv,ldu,iopt,info)
!      IF (info .NE. 0) THEN
!          info = ppm_error_error
!          CALL ppm_error(ppm_err_dealloc,caller,     &
!     &        'receive counter NRECV',__LINE__,info)
!      ENDIF
!      CALL ppm_alloc(psend,ldu,iopt,info)
!      IF (info .NE. 0) THEN
!          info = ppm_error_error
!          CALL ppm_error(ppm_err_dealloc,caller,     &
!     &        'particle send counter PSEND',__LINE__,info)
!      ENDIF
!      CALL ppm_alloc(precv,ldu,iopt,info)
!      IF (info .NE. 0) THEN
!          info = ppm_error_error
!          CALL ppm_error(ppm_err_dealloc,caller,     &
!     &        'particle receive counter PRECV',__LINE__,info)
!      ENDIF
!      CALL ppm_alloc(   pp,ldu,iopt,info)
!      IF (info .NE. 0) THEN
!          info = ppm_error_error
!          CALL ppm_error(ppm_err_dealloc,caller,     &
!     &        'work array PP',__LINE__,info)
!      ENDIF
!      CALL ppm_alloc(   qq,ldu,iopt,info)
!      IF (info .NE. 0) THEN
!          info = ppm_error_error
!          CALL ppm_error(ppm_err_dealloc,caller,     &
!     &        'work array QQ',__LINE__,info)
!      ENDIF
      CALL ppm_alloc(recvd,ldu,iopt,info)
      or_fail_dealloc('local receive buffer RECVD',exit_point=no)

      CALL ppm_alloc(recvs,ldu,iopt,info)
      or_fail_dealloc('local receive buffer RECVS',exit_point=no)

      CALL ppm_alloc(sendd,ldu,iopt,info)
      or_fail_dealloc('local send buffer SENDD',exit_point=no)

      CALL ppm_alloc(sends,ldu,iopt,info)
      or_fail_dealloc('local send buffer SENDS',exit_point=no)

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
        IF (Npart .LT. 0) THEN
           fail('Npart must be >=0',exit_point=8888)
        ENDIF
      8888 CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_map_part_send
