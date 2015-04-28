      !-------------------------------------------------------------------------
      !  Subroutine   :                 map_part_send
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
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

      SUBROUTINE DTYPE(map_part_send)(Pc,mapID,info)
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
      USE ppm_module_mpi
      IMPLICIT NONE

      DEFINE_MK()
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(DTYPE(ppm_t_particles))           :: Pc
      !!! Particle set
      INTEGER                 , INTENT(IN   ) :: mapID
      !!! mapping ID
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      TYPE(DTYPE(ppm_t_part_mapping)), POINTER :: map

      REAL(MK), DIMENSION(:), POINTER  :: send
      REAL(MK), DIMENSION(:), POINTER  :: recv
      REAL(ppm_kind_double)            :: t0

      INTEGER, DIMENSION(:,:), POINTER  :: pp
      INTEGER, DIMENSION(:,:), POINTER  :: qq
      INTEGER, DIMENSION(3)             :: ldl,ldu
      INTEGER                           :: i,j,k,idom,nbuffer,ibuffer,jbuffer
      INTEGER                           :: bdim,sbdim
      INTEGER                           :: iopt,count,tag1,qpart,msend,mrecv
      INTEGER                           :: npart_send,npart_recv,Mpart
#ifdef __MPI
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
#endif

      CHARACTER(ppm_char) :: mesg
      CHARACTER(ppm_char) :: caller ='map_part_send'

      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)
      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 3) THEN
         CALL check
         IF (info .NE. 0) GOTO 9999
      ENDIF

      map => Pc%maps%vec(mapID)%t

      ! skip if the buffer is empty
      IF (map%ppm_buffer_set .LT. 1) THEN
         IF (ppm_debug .GT. 1) THEN
            fail('Buffer is empty: skipping send!',ppm_err_buffer_empt,exit_point=no,ppm_error=ppm_error_notice)
            info=0
         ENDIF
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_grow
      ldu(1) = MAX(map%ppm_nsendlist,1)
      IF ((map%ppm_nsendlist.NE.map%old_nsendlist) .OR.  &
      &   (map%ppm_buffer_set.NE.map%old_buffer_set)) THEN
          map%old_nsendlist = map%ppm_nsendlist
          map%old_buffer_set = map%ppm_buffer_set
          CALL ppm_alloc(map%nsend,ldu,iopt,info)
          or_fail_alloc('send counter NSEND',ppm_error=ppm_error_fatal)

          CALL ppm_alloc(map%nrecv,ldu,iopt,info)
          or_fail_alloc('receive counter NRECV',ppm_error=ppm_error_fatal)

          CALL ppm_alloc(map%psend,ldu,iopt,info)
          or_fail_alloc('particle send counter PSEND',ppm_error=ppm_error_fatal)

          CALL ppm_alloc(map%precv,ldu,iopt,info)
          or_fail_alloc('particle receive counter PRECV',ppm_error=ppm_error_fatal)

          ldu(2) = map%ppm_buffer_set
          CALL ppm_alloc(map%pp,ldu,iopt,info)
          or_fail_alloc('work buffer PP',ppm_error=ppm_error_fatal)

          CALL ppm_alloc(map%qq,ldu,iopt,info)
          or_fail_alloc('work buffer QQ',ppm_error=ppm_error_fatal)
      ENDIF

      send => map%send
      recv => map%recv
      pp   => map%pp
      qq   => map%qq

      !-------------------------------------------------------------------------
      !  Count the total size of the buffer dimensions
      !  It is handy... simplifies many summations, avoids loops, etc...
      !-------------------------------------------------------------------------
      sbdim = SUM(map%ppm_buffer_dim(1:map%ppm_buffer_set))

      !-------------------------------------------------------------------------
      !  Count the size of the buffer that will NOT be sent
      !-------------------------------------------------------------------------
      qpart   = (map%ppm_psendbuffer(2) - map%ppm_psendbuffer(1))
      ibuffer = sbdim*qpart

      !-------------------------------------------------------------------------
      !  Initialize the counter for the total set of new particles
      !JHW 20060928     IF (ppm_map_type.EQ.ppm_param_map_ghost_get) THEN
      !-------------------------------------------------------------------------
      IF (map%ppm_map_type.EQ.ppm_param_map_ghost_get.OR. &
      &   map%ppm_map_type.EQ.ppm_param_map_ghost_put) THEN
         Mpart = qpart + map%oldNpart
      ELSE
         Mpart = qpart
      ENDIF

      !-------------------------------------------------------------------------
      !  Count the size of the buffer that WILL be sent
      !-------------------------------------------------------------------------
      map%ppm_nrecvbuffer = ibuffer
      map%nsend(1)        = ibuffer
      map%nrecv(1)        = ibuffer
      map%psend(1)        = qpart
      map%precv(1)        = qpart
      mrecv               = -1
      msend               = -1

      DO k=2,map%ppm_nsendlist

         !----------------------------------------------------------------------
         !  The number of particles send off to the k-th processor in the
         !  sendlist
         !----------------------------------------------------------------------
         qpart=(map%ppm_psendbuffer(k+1)-map%ppm_psendbuffer(k))

         !----------------------------------------------------------------------
         !  Store the number of particles to send
         !----------------------------------------------------------------------
         map%psend(k) = qpart ! store the number of particles to send

         !----------------------------------------------------------------------
         !  Store the size of the data to be send
         !----------------------------------------------------------------------
         map%nsend(k) = sbdim*qpart

#ifdef __MPI
         !----------------------------------------------------------------------
         !  Make a send/recv of the number of particles and data size to that
         !  has to be send/recv
         !----------------------------------------------------------------------
         ! The following IF is needed in order to skip "dummy"
         ! communication rounds where the current processor has to wait
         ! (only needed in the partial mapping).
         IF (map%ppm_isendlist(k) .GE. 0 .AND. map%ppm_irecvlist(k) .GE. 0) THEN
             tag1 = 100
             IF (ppm_debug .GT. 1) THEN
                stdout_f('(A,I5,2(A,I9))',"sending to ",'map%ppm_isendlist(k)', &
                & ", nsend=",'map%nsend(k)',", psend=",'map%psend(k)')
             ENDIF
             CALL MPI_SendRecv(map%psend(k),1,MPI_INTEGER,map%ppm_isendlist(k),tag1, &
             & map%precv(k),1,MPI_INTEGER,map%ppm_irecvlist(k),tag1, &
             & ppm_comm,status,info)
             or_fail_MPI("MPI_SendRecv")

             ! Compute nrecv(k) from precv(k)
             map%nrecv(k) = sbdim*map%precv(k)

             IF (ppm_debug .GT. 1) THEN
                stdout_f('(A,I5,2(A,I9))',"received from ",'map%ppm_irecvlist(k)', &
                & ", nrecv=",'map%nrecv(k)',", precv=",'map%precv(k)')
             ENDIF
         ELSE
             ! skip this round, i.e. neither send nor receive any
             ! particles.
             map%nrecv(k) = 0
             map%precv(k) = 0
         ENDIF
#endif

      ENDDO
      !----------------------------------------------------------------------
      !  Find the required (maximum) size of the send/recv buffers
      !----------------------------------------------------------------------
      DO k=2,map%ppm_nsendlist
         IF (mrecv.LE.map%nrecv(k)) mrecv = map%nrecv(k)
      ENDDO
      DO k=2,map%ppm_nsendlist
         IF (msend.LE.map%nsend(k)) msend = map%nsend(k)
      ENDDO
      !----------------------------------------------------------------------
      !  Increment the total receive buffer
      !----------------------------------------------------------------------
      DO k=2,map%ppm_nsendlist
         map%ppm_nrecvbuffer = map%ppm_nrecvbuffer + map%nrecv(k)
      ENDDO
      !----------------------------------------------------------------------
      !  Increment the total number of particle to receive
      !----------------------------------------------------------------------
      DO k=2,map%ppm_nsendlist
         Mpart = Mpart + map%precv(k)
      ENDDO

      IF (ppm_debug .GT. 1) THEN
         stdout_f('(2(A,I9))',"mrecv=",mrecv,", msend=",msend)
         stdout_f('(A,I9)',"ppm_nrecvbuffer=",'map%ppm_nrecvbuffer')
         stdout_f('(A,I9)',"Mpart=",Mpart)
      ENDIF

      map%newNpart = Mpart

      !-------------------------------------------------------------------------
      !  Allocate the memory for the copy of the particle buffer
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow
      ldu(1) = map%ppm_nrecvbuffer
      CALL ppm_alloc(map%ppm_recvbuffer,ldu,iopt,info)
      or_fail_alloc('global receive buffer PPM_RECVBUFFER',ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Allocate memory for the smaller send and receive buffer
      !-------------------------------------------------------------------------
      ! only allocate if there actually is anything to be sent/recvd with
      ! other processors. Otherwise mrecv and msend would both still be -1
      ! (as initialized above) and the alloc would throw a FATAL. This was
      ! Bug ID 000012.
      IF (map%ppm_nsendlist .GT. 1) THEN
         iopt   = ppm_param_alloc_grow
         ldu(1) = MAX(mrecv,1)
         CALL ppm_alloc(map%recv,ldu,iopt,info)
         or_fail_alloc('local receive buffer RECV',ppm_error=ppm_error_fatal)

         recv => map%recv

         ldu(1) = MAX(msend,1)
         CALL ppm_alloc(map%send,ldu,iopt,info)
         or_fail_alloc('local send buffer SEND',ppm_error=ppm_error_fatal)

         send => map%send
      ENDIF

      !-------------------------------------------------------------------------
      !  Debugging print of the number of sets in the buffer
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
         stdout_f('(A,I9)',"ppm_buffer_set=",'map%ppm_buffer_set')
      ENDIF

      !-------------------------------------------------------------------------
      !  Compute the total number of particles to send
      !  As it is now, this should be equal to Npart
      !-------------------------------------------------------------------------
      npart_send =SUM(map%psend(1:map%ppm_nsendlist))

      !-------------------------------------------------------------------------
      !  Compute the pointer to the position of the data in the main send
      !  buffer
      !-------------------------------------------------------------------------
      DO k=1,map%ppm_buffer_set
         IF (k.EQ.1) THEN
            qq(1,k) = 1
         ELSE
            qq(1,k) = qq(1,k-1) + npart_send*map%ppm_buffer_dim(k-1)
         ENDIF
         bdim    = map%ppm_buffer_dim(k)
         DO j=2,map%ppm_nsendlist
            qq(j,k) = qq(j-1,k) + map%psend(j-1)*bdim
            IF (ppm_debug .GT. 1) THEN
               stdout_f('(A,I9)',"qq(j,k)=",'qq(j,k)')
            ENDIF
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Compute the total number of particles to receive
      !  As it is now, this should be equal to Mpart
      !-------------------------------------------------------------------------
      npart_recv=SUM(map%precv(1:map%ppm_nsendlist))

      !-------------------------------------------------------------------------
      !  Compute the pointer to the position of the data in the main receive
      !  buffer
      !-------------------------------------------------------------------------
      DO k=1,map%ppm_buffer_set
         IF (k.EQ.1) THEN
            pp(1,k) = 1
         ELSE
            pp(1,k) = pp(1,k-1) + npart_recv*map%ppm_buffer_dim(k-1)
         ENDIF
         bdim    = map%ppm_buffer_dim(k)
         DO j=2,map%ppm_nsendlist
            pp(j,k) = pp(j-1,k) + map%precv(j-1)*bdim
            IF (ppm_debug .GT. 1) THEN
               stdout_f('(A,I9)',"pp(j,k)=",'pp(j,k)')
               stdout_f('(A,I9,A,I4)',"precv(j-1)=",'map%precv(j-1)',", bdim=",bdim)
            ENDIF
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  First copy the on processor data - which is in the first buffer
      !-------------------------------------------------------------------------
      DO k=1,map%ppm_buffer_set
         ibuffer = pp(1,k) - 1
         jbuffer = qq(1,k) - 1
         DO j=1,map%psend(1)*map%ppm_buffer_dim(k)
            ibuffer                  = ibuffer + 1
            jbuffer                  = jbuffer + 1
            map%ppm_recvbuffer(ibuffer) = map%ppm_sendbuffer(jbuffer)
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  loop over the processors in the ppm_isendlist(); skip the first entry
      !  which is the local processor
      !-------------------------------------------------------------------------
      !----------------------------------------------------------------------
      !  For each send/recv
      !----------------------------------------------------------------------
      DO k=2,map%ppm_nsendlist
         !-------------------------------------------------------------------
         !  Collect each data type (xp, vp, etc)
         !-------------------------------------------------------------------
         ibuffer = 0
         DO j=1,map%ppm_buffer_set
            jbuffer = qq(k,j) - 1
            !----------------------------------------------------------------
            !  Collect the data into the send buffer
            !----------------------------------------------------------------
            DO i=1,map%psend(k)*map%ppm_buffer_dim(j)
               ibuffer        = ibuffer + 1
               jbuffer        = jbuffer + 1
               send(ibuffer) = map%ppm_sendbuffer(jbuffer)
            ENDDO
         ENDDO

         !-------------------------------------------------------------------
         !  Perform the actual send/recv
         !-------------------------------------------------------------------
#ifdef __MPI
         ! The following IF is needed in order to skip "dummy"
         ! communication rounds where the current processor has to wait
         ! (only needed in the partial mapping).
         IF (map%ppm_isendlist(k) .GE. 0 .AND. map%ppm_irecvlist(k) .GE. 0) THEN
            tag1 = 300
            CALL MPI_SendRecv(send,map%nsend(k),ppm_mpi_kind, &
            &    map%ppm_isendlist(k),tag1,recv,map%nrecv(k), &
            &    ppm_mpi_kind,map%ppm_irecvlist(k),tag1,      &
            &    ppm_comm,status,info)
            or_fail_MPI("MPI_SendRecv")
         ENDIF
#else
         recv = send
#endif
         !-------------------------------------------------------------------
         !  Store the data back in the recv buffer
         !-------------------------------------------------------------------
         ibuffer = 0
         DO j=1,map%ppm_buffer_set
            jbuffer = pp(k,j) - 1
            DO i=1,map%precv(k)*map%ppm_buffer_dim(j)
               ibuffer                  = ibuffer + 1
               jbuffer                  = jbuffer + 1
               map%ppm_recvbuffer(jbuffer) = recv(ibuffer)
           ENDDO
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  before we through away the precv() data let us store it for later use:
      !  when sending ghosts back (ppm_map_part_ghost_put())
      !-------------------------------------------------------------------------
      IF (map%ppm_map_type.EQ.ppm_param_map_ghost_get) THEN
         ldu(1) = map%ppm_nsendlist + 1
         iopt   = ppm_param_alloc_grow
         CALL ppm_alloc(map%ppm_precvbuffer,ldu,iopt,info)
         or_fail_alloc('global recv buffer pointer PPM_PRECVBUFFER',ppm_error=ppm_error_fatal)

         map%ppm_precvbuffer(1) = map%oldNpart + 1
         DO k=1,map%ppm_nsendlist
            map%ppm_precvbuffer(k+1) = map%ppm_precvbuffer(k) + map%precv(k)
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
      CALL ppm_alloc(map%ppm_sendbuffer,ldu,iopt,info)
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
      CALL ppm_alloc(map%recv,ldu,iopt,info)
      or_fail_dealloc('local receive buffer RECV',exit_point=no)

      CALL ppm_alloc(map%send,ldu,iopt,info)
      or_fail_dealloc('local send buffer SEND',exit_point=no)

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
        IF (.NOT. Pc%maps%exists(mapID)) THEN
           fail('Invalid mapID: mapping does not exist.',ppm_err_wrong_dim,exit_point=8888)
        ENDIF
        IF (Pc%maps%vec(mapID)%t%oldNpart .LT. 0) THEN
          fail('Npart must be >=0',exit_point=8888)
        ENDIF
      8888 CONTINUE
      END SUBROUTINE check
      END SUBROUTINE DTYPE(map_part_send)

#undef __KIND
