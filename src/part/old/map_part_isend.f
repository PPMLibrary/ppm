      !-------------------------------------------------------------------------
      !  Subroutine   :                 map_part_isend
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
      !-------------------------------------------------------------------------
      !  MOSAIC Group
      !  Max Planck Institute of Molecular Cell Biology and Genetics
      !  Pfotenhauerstr. 108, 01307 Dresden, Germany
      !-------------------------------------------------------------------------

      SUBROUTINE DTYPE(map_part_isend)(Pc,mapID,info)
      !!! This routine performs the actual send/recv of the particles and all
      !!! pushed data.
      !!!
      !!! [NOTE]
      !!! The first part of the buffer contains the on processor data. The
      !!! packing could be performed more efficiently.
      !!!
      !!! [TIP] Yaser:
      !!! This is not the best implementation
      !!! I could not figure out the general derived type to use for MPI_Isend &
      !!! MPI_Irecv to avoid several mpi call, so this implementation will call
      !!! mpi for each ppm_buffer_set, but using nonblocking send & recv it will
      !!! prevent several buffer copy which is very useful on big data !
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
      CLASS(DTYPE(ppm_t_particles)) :: Pc
      !!! Particle set
      INTEGER,        INTENT(IN   ) :: mapID
      !!! mapping ID
      INTEGER,        INTENT(  OUT) :: info
      !!! Return status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      TYPE(DTYPE(ppm_t_part_mapping)), POINTER :: map

      REAL(MK) :: t0

      INTEGER, DIMENSION(3)                :: ldu
      INTEGER                              :: i,j,k,l,m
      INTEGER                              :: bdim,sbdim
      INTEGER                              :: iopt,tag,nsr
      INTEGER                              :: qpart,Mpart
      INTEGER                              :: npart_send,npart_recv
#ifdef __MPI
      INTEGER, DIMENSION(:,:), POINTER     :: pp
      INTEGER, DIMENSION(:,:), POINTER     :: qq
      INTEGER, DIMENSION(:),   ALLOCATABLE :: request
#endif

      CHARACTER(ppm_char) :: caller ='map_part_isend'

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
      IF (ppm_debug.GT.3) THEN
         CALL check
         IF (info.NE.0) GOTO 9999
      ENDIF

      map => Pc%maps%vec(mapID)%t

      ! skip if the buffer is empty
      IF (map%ppm_buffer_set.LT.1) THEN
         IF (ppm_debug.GT.1) THEN
            fail('Buffer is empty: skipping send!', &
            & ppm_err_buffer_empt,exit_point=no,ppm_error=ppm_error_notice)
            info=0
         ENDIF
         GOTO 9999
      ELSE
         !-------------------------------------------------------------------------
         !  Debugging print of the number of sets in the buffer
         !-------------------------------------------------------------------------
         IF (ppm_debug.GT.1) THEN
            stdout_f('(A,I9)',"ppm_buffer_set=",'map%ppm_buffer_set')
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate
      !-------------------------------------------------------------------------
      IF ((map%ppm_nsendlist.NE.map%old_nsendlist) .OR. &
      &   (map%ppm_buffer_set.NE.map%old_buffer_set)) THEN
         iopt = ppm_param_alloc_grow
         ldu(1) = MAX(map%ppm_nsendlist,1)

         CALL ppm_alloc(map%psend,ldu,iopt,info)
         or_fail_alloc('particle send counter PSEND',ppm_error=ppm_error_fatal)

         CALL ppm_alloc(map%precv,ldu,iopt,info)
         or_fail_alloc('particle receive counter PRECV',ppm_error=ppm_error_fatal)
      ENDIF

      !-------------------------------------------------------------------------
      !  Count the size of the buffer that will NOT be sent
      !-------------------------------------------------------------------------
      qpart=map%ppm_psendbuffer(2)-map%ppm_psendbuffer(1)

      map%psend(1)=qpart
      map%precv(1)=qpart

#ifdef __MPI
      ALLOCATE(request(ppm_nproc*map%ppm_buffer_set*2),STAT=info)
      or_fail_alloc("request")

      !counter for number of send & recv
      nsr=0

      DO k=2,map%ppm_nsendlist
         !----------------------------------------------------------------------
         !  Store the number of particles send off to the k-th processor in the
         !  sendlist
         !----------------------------------------------------------------------
         map%psend(k) = map%ppm_psendbuffer(k+1)-map%ppm_psendbuffer(k)

         !----------------------------------------------------------------------
         !  Make a send/recv of the number of particles and data size to that
         !  has to be send/recv
         !----------------------------------------------------------------------
         ! The following IF is needed in order to skip "dummy"
         ! communication rounds where the current processor has to wait
         ! (only needed in the partial mapping).
         IF (map%ppm_isendlist(k).GE.0.AND.map%ppm_irecvlist(k).GE.0) THEN
            IF (ppm_rank.GT.map%ppm_isendlist(k)) THEN
               tag=ppm_rank*(ppm_rank-1)/2+map%ppm_isendlist(k)+k
            ELSE
               tag=map%ppm_isendlist(k)*(map%ppm_isendlist(k)-1)/2+ppm_rank+k
            ENDIF

            nsr=nsr+1
            CALL MPI_Isend(map%psend(k),1,MPI_INTEGER,map%ppm_isendlist(k), &
            &    tag,ppm_comm,request(nsr),info)
            or_fail_MPI("MPI_Isend")

            IF (ppm_debug.GT.1) THEN
               stdout_f('(A,I5,A,I9)',"sending to ",'map%ppm_isendlist(k)', &
               & ", psend=",'map%psend(k)')
            ENDIF

            IF (ppm_rank.GT.map%ppm_irecvlist(k)) THEN
               tag=ppm_rank*(ppm_rank-1)/2+map%ppm_irecvlist(k)+k
            ELSE
               tag=map%ppm_irecvlist(k)*(map%ppm_irecvlist(k)-1)/2+ppm_rank+k
            ENDIF

            nsr=nsr+1
            CALL MPI_Irecv(map%precv(k),1,MPI_INTEGER,map%ppm_irecvlist(k), &
            &    tag,ppm_comm,request(nsr),info)
            or_fail_MPI("MPI_Irecv")
         ELSE
            ! skip this round, i.e. neither send nor receive any
            ! particles.
            map%precv(k) = 0
         ENDIF
      ENDDO !k=2,map%ppm_nsendlist

      !-------------------------------------------------------------------------
      !  Allocate
      !-------------------------------------------------------------------------
      IF ((map%ppm_buffer_set.NE.map%old_buffer_set).OR. &
      &   (map%ppm_nsendlist.NE.map%old_nsendlist)) THEN
         map%old_nsendlist = map%ppm_nsendlist
         map%old_buffer_set = map%ppm_buffer_set

         iopt = ppm_param_alloc_grow
         ldu(1) = MAX(map%ppm_nsendlist,1)
         ldu(2) = map%ppm_buffer_set

         CALL ppm_alloc(map%pp,ldu,iopt,info)
         or_fail_alloc('work buffer PP',ppm_error=ppm_error_fatal)

         CALL ppm_alloc(map%qq,ldu,iopt,info)
         or_fail_alloc('work buffer QQ',ppm_error=ppm_error_fatal)
      ENDIF

      pp => map%pp
      qq => map%qq
#endif
      !-------------------------------------------------------------------------
      !  Count the total size of the buffer dimensions
      !  It is handy... simplifies many summations, avoids loops, etc...
      !-------------------------------------------------------------------------
      sbdim = SUM(map%ppm_buffer_dim(1:map%ppm_buffer_set))

#ifdef __MPI
      !-------------------------------------------------------------------------
      !  Compute the total number of particles to send
      !  As it is now, this should be equal to Npart
      !-------------------------------------------------------------------------
      npart_send=SUM(map%psend(1:map%ppm_nsendlist))

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
         bdim = map%ppm_buffer_dim(k)
         DO j=2,map%ppm_nsendlist
            qq(j,k) = qq(j-1,k) + map%psend(j-1)*bdim
            IF (ppm_debug .GT. 1) THEN
               stdout_f('(A,I9)',"qq(j,k)=",'qq(j,k)')
            ENDIF
         ENDDO
      ENDDO

      IF (nsr.GT.0) THEN
         CALL MPI_Waitall(nsr,request(1:nsr),MPI_STATUSES_IGNORE,info)
         or_fail_MPI("MPI_Waitall")
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Compute the total number of particles to receive
      !  As it is now, this should be equal to Mpart
      !-------------------------------------------------------------------------
      npart_recv=SUM(map%precv(1:map%ppm_nsendlist))

      !-------------------------------------------------------------------------
      !  Count the size of the receive buffer
      !----------------------------------------------------------------------
      map%ppm_nrecvbuffer=sbdim*npart_recv

      !-------------------------------------------------------------------------
      !  Counter for the total set of new particles
      !  JHW 20060928
      !-------------------------------------------------------------------------
      SELECT CASE (map%ppm_map_type)
      CASE (ppm_param_map_ghost_get,ppm_param_map_ghost_put)
         !----------------------------------------------------------------------
         !  Count the total number of particle to receive
         !----------------------------------------------------------------------
         Mpart = map%oldNpart+npart_recv

      CASE DEFAULT
         !----------------------------------------------------------------------
         !  Count the total number of particle to receive
         !----------------------------------------------------------------------
         Mpart = npart_recv

      END SELECT

      IF (ppm_debug .GT. 1) THEN
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

#ifdef __MPI
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

      !counter for number of send & recv
      nsr=0

      DO k=1,map%ppm_buffer_set
         bdim = map%ppm_buffer_dim(k)
         !-------------------------------------------------------------------------
         !  loop over the processors in the ppm_isendlist(); skip the first entry
         !  which is the local processor
         !-------------------------------------------------------------------------
         !----------------------------------------------------------------------
         !  For each send/recv
         !----------------------------------------------------------------------
         DO i=2,map%ppm_nsendlist
            !-------------------------------------------------------------------
            !  Perform the actual send/recv
            !-------------------------------------------------------------------
            ! The following IF is needed in order to skip "dummy"
            ! communication rounds where the current processor has to wait
            ! (only needed in the partial mapping).
            IF (map%ppm_isendlist(i).LT.0.OR.map%ppm_irecvlist(i).LT.0) CYCLE
            IF (map%psend(i).GT.0) THEN
               IF (ppm_rank.GT.map%ppm_isendlist(i)) THEN
                  tag=ppm_rank*(ppm_rank-1)/2+map%ppm_isendlist(i)+k
               ELSE
                  tag=map%ppm_isendlist(i)*(map%ppm_isendlist(i)-1)/2+ppm_rank+k
               ENDIF

               nsr=nsr+1
               CALL MPI_Isend(map%ppm_sendbuffer(qq(i,k)),map%psend(i)*bdim, &
               &    ppm_mpi_kind,map%ppm_isendlist(i),tag,ppm_comm,          &
               &    request(nsr),info)
               or_fail_MPI("MPI_Isend")
            ENDIF
            IF (map%precv(i).GT.0) THEN
               IF (ppm_rank.GT.map%ppm_irecvlist(i)) THEN
                  tag=ppm_rank*(ppm_rank-1)/2+map%ppm_irecvlist(i)+k
               ELSE
                  tag=map%ppm_irecvlist(i)*(map%ppm_irecvlist(i)-1)/2+ppm_rank+k
               ENDIF

               nsr=nsr+1
               CALL MPI_Irecv(map%ppm_recvbuffer(pp(i,k)),map%precv(i)*bdim, &
               &    ppm_mpi_kind,map%ppm_irecvlist(i),tag,ppm_comm,          &
               &    request(nsr),info)
               or_fail_MPI("MPI_Irecv")
            ENDIF
         ENDDO !i=2,map%ppm_nsendlist
      ENDDO !k=1,map%ppm_buffer_set
      DO k=1,map%ppm_buffer_set
         !----------------------------------------------------------------------
         !  For each processor
         !----------------------------------------------------------------------
         bdim = map%ppm_buffer_dim(k)
         l=qq(1,k)-1
         m=pp(1,k)-1
         DO j=1,map%psend(1)*bdim
            map%ppm_recvbuffer(m+j)=map%ppm_sendbuffer(l+j)
         ENDDO
      ENDDO !k=1,map%ppm_buffer_set
#else
      map%ppm_recvbuffer=map%ppm_sendbuffer
#endif

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

#ifdef __MPI
      IF (nsr.GT.0) THEN
         CALL MPI_Waitall(nsr,request(1:nsr),MPI_STATUSES_IGNORE,info)
         or_fail_MPI("MPI_Waitall")
      ENDIF

      !-------------------------------------------------------------------------
      !  Deallocate
      !-------------------------------------------------------------------------
      DEALLOCATE(request,STAT=info)
      or_fail_dealloc("request")
#endif

      !-------------------------------------------------------------------------
      !  Deallocate the send buffer to save memory
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(map%ppm_sendbuffer,ldu,iopt,info)
      or_fail_dealloc('send buffer PPM_SENDBUFFER',exit_point=no)

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
      END SUBROUTINE DTYPE(map_part_isend)

#undef __KIND
