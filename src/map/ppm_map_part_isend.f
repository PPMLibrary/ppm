      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_map_part_isend
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
      !-------------------------------------------------------------------------
      !  MOSAIC Group
      !  Max Planck Institute of Molecular Cell Biology and Genetics
      !  Pfotenhauerstr. 108, 01307 Dresden, Germany
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_map_part_isend(Npart,Mpart,info)
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

      INTEGER, DIMENSION(2)              :: ldu
      INTEGER                            :: i,j,k,l,m
      INTEGER                            :: bdim,sbdim
      INTEGER                            :: iopt,tag,nsr
      INTEGER                            :: qpart
      INTEGER                            :: npart_send,npart_recv
#ifdef __MPI
      INTEGER, DIMENSION(:), ALLOCATABLE :: request
#endif

      CHARACTER(ppm_char) :: caller='ppm_map_part_isend'
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

      ! skip if the buffer is empty
      IF (ppm_buffer_set.LT.1) THEN
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
            stdout_f('(A,I9)',"ppm_buffer_set=",ppm_buffer_set)
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate
      !-------------------------------------------------------------------------
      IF ((ppm_nsendlist.NE.old_nsendlist) .OR.  &
      &  (ppm_buffer_set.NE.old_buffer_set)) THEN
         iopt = ppm_param_alloc_grow
         ldu(1) = MAX(ppm_nsendlist,1)

         CALL ppm_alloc(psend,ldu,iopt,info)
         or_fail_alloc('particle send counter PSEND',ppm_error=ppm_error_fatal)

         CALL ppm_alloc(precv,ldu,iopt,info)
         or_fail_alloc('particle receive counter PRECV',ppm_error=ppm_error_fatal)
      ENDIF

      !-------------------------------------------------------------------------
      !  Count the size of the buffer that will NOT be sent
      !-------------------------------------------------------------------------
      qpart   =ppm_psendbuffer(2)-ppm_psendbuffer(1)
      psend(1)=qpart
      precv(1)=qpart

#ifdef __MPI
      ALLOCATE(request(ppm_nproc*ppm_buffer_set*2),STAT=info)
      or_fail_alloc("request")

      !counter for number of send & recv
      nsr=0

      DO k=2,ppm_nsendlist
         !----------------------------------------------------------------------
         !  Store the number of particles send off to the k-th processor in the
         !  sendlist
         !----------------------------------------------------------------------
         psend(k)=ppm_psendbuffer(k+1)-ppm_psendbuffer(k)

         !----------------------------------------------------------------------
         !  Make a send/recv of the number of particles and data size to that
         !  has to be send/recv
         !----------------------------------------------------------------------
         ! The following IF is needed in order to skip "dummy"
         ! communication rounds where the current processor has to wait
         ! (only needed in the partial mapping).
         IF (ppm_isendlist(k).GE.0.AND.ppm_irecvlist(k).GE.0) THEN
           IF (ppm_rank.GT.ppm_irecvlist(k)) THEN
               tag=ppm_rank*(ppm_rank-1)/2+ppm_irecvlist(k)+k
            ELSE
               tag=ppm_irecvlist(k)*(ppm_irecvlist(k)-1)/2+ppm_rank+k
            ENDIF

            nsr=nsr+1
            CALL MPI_Irecv(precv(k),1,MPI_INTEGER,ppm_irecvlist(k), &
            &    tag,ppm_comm,request(nsr),info)
            or_fail_MPI("MPI_Irecv")

            IF (ppm_rank.GT.ppm_isendlist(k)) THEN
               tag=ppm_rank*(ppm_rank-1)/2+ppm_isendlist(k)+k
            ELSE
               tag=ppm_isendlist(k)*(ppm_isendlist(k)-1)/2+ppm_rank+k
            ENDIF

            nsr=nsr+1
            CALL MPI_Isend(psend(k),1,MPI_INTEGER,ppm_isendlist(k), &
            &    tag,ppm_comm,request(nsr),info)
            or_fail_MPI("MPI_Isend")

            IF (ppm_debug.GT.1) THEN
               stdout_f('(A,I5,A,I9)',"sending to ",'ppm_isendlist(k)', &
               & ", psend=",'psend(k)')
            ENDIF
         ELSE
            ! skip this round, i.e. neither send nor receive any
            ! particles.
            precv(k) = 0
         ENDIF
      ENDDO !k=2,ppm_nsendlist

      !-------------------------------------------------------------------------
      !  Allocate
      !-------------------------------------------------------------------------
      IF ((ppm_buffer_set.NE.old_buffer_set).OR. &
      &   (ppm_nsendlist.NE.old_nsendlist)) THEN
         old_nsendlist = ppm_nsendlist
         old_buffer_set = ppm_buffer_set

         iopt = ppm_param_alloc_grow
         ldu(1) = MAX(ppm_nsendlist,1)
         ldu(2) = ppm_buffer_set

         CALL ppm_alloc(pp,ldu,iopt,info)
         or_fail_alloc('work buffer PP',ppm_error=ppm_error_fatal)

         CALL ppm_alloc(qq,ldu,iopt,info)
         or_fail_alloc('work buffer QQ',ppm_error=ppm_error_fatal)
      ENDIF

#endif
      !-------------------------------------------------------------------------
      !  Count the total size of the buffer dimensions
      !  It is handy... simplifies many summations, avoids loops, etc...
      !-------------------------------------------------------------------------
      sbdim = SUM(ppm_buffer_dim(1:ppm_buffer_set))

#ifdef __MPI
      !-------------------------------------------------------------------------
      !  Compute the total number of particles to send
      !  As it is now, this should be equal to Npart
      !-------------------------------------------------------------------------
      npart_send=SUM(psend(1:ppm_nsendlist))

      !-------------------------------------------------------------------------
      !  Compute the pointer to the position of the data in the main send
      !  buffer
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
         DO k=1,ppm_buffer_set
            IF (k.EQ.1) THEN
               qq(1,k) = 1
            ELSE
               qq(1,k) = qq(1,k-1) + npart_send*ppm_buffer_dim(k-1)
            ENDIF
            bdim = ppm_buffer_dim(k)
            DO j=2,ppm_nsendlist
               qq(j,k) = qq(j-1,k) + psend(j-1)*bdim
               stdout_f('(A,I9)',"qq(j,k)=",'qq(j,k)')
            ENDDO
         ENDDO
      ELSE
         DO k=1,ppm_buffer_set
            IF (k.EQ.1) THEN
               qq(1,k) = 1
            ELSE
               qq(1,k) = qq(1,k-1) + npart_send*ppm_buffer_dim(k-1)
            ENDIF
            bdim = ppm_buffer_dim(k)
            DO j=2,ppm_nsendlist
               qq(j,k) = qq(j-1,k) + psend(j-1)*bdim
            ENDDO
         ENDDO
      ENDIF

      IF (nsr.GT.0) THEN
         CALL MPI_Waitall(nsr,request(1:nsr),MPI_STATUSES_IGNORE,info)
         or_fail_MPI("MPI_Waitall")
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Compute the total number of particles to receive
      !  As it is now, this should be equal to Mpart
      !-------------------------------------------------------------------------
      npart_recv=SUM(precv(1:ppm_nsendlist))

      !-------------------------------------------------------------------------
      !  Count the size of the receive buffer
      !----------------------------------------------------------------------
      ppm_nrecvbuffer=sbdim*npart_recv

      !-------------------------------------------------------------------------
      !  Counter for the total set of new particles
      !  JHW 20060928
      !-------------------------------------------------------------------------
      SELECT CASE (ppm_map_type)
      CASE (ppm_param_map_ghost_get,ppm_param_map_ghost_put)
         !----------------------------------------------------------------------
         !  Count the total number of particle to receive
         !----------------------------------------------------------------------
         Mpart = Npart+npart_recv

      CASE DEFAULT
         !----------------------------------------------------------------------
         !  Count the total number of particle to receive
         !----------------------------------------------------------------------
         Mpart = npart_recv

      END SELECT

      IF (ppm_debug .GT. 1) THEN
         stdout_f('(A,I9)',"ppm_nrecvbuffer=",'ppm_nrecvbuffer')
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
      or_fail_alloc('global receive buffer PPM_RECVBUFFER',ppm_error=ppm_error_fatal)

#ifdef __MPI
      !-------------------------------------------------------------------------
      !  Compute the pointer to the position of the data in the main receive
      !  buffer
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.1) THEN
         DO k=1,ppm_buffer_set
            IF (k.EQ.1) THEN
               pp(1,k) = 1
            ELSE
               pp(1,k) = pp(1,k-1) + npart_recv*ppm_buffer_dim(k-1)
            ENDIF
            bdim = ppm_buffer_dim(k)
            DO j=2,ppm_nsendlist
               pp(j,k) = pp(j-1,k) + precv(j-1)*bdim
               stdout_f('(A,I9)',"pp(j,k)=",'pp(j,k)')
               stdout_f('(A,I9,A,I4)',"precv(j-1)=",'precv(j-1)',", bdim=",bdim)
            ENDDO
         ENDDO
      ELSE
         DO k=1,ppm_buffer_set
            IF (k.EQ.1) THEN
               pp(1,k) = 1
            ELSE
               pp(1,k) = pp(1,k-1) + npart_recv*ppm_buffer_dim(k-1)
            ENDIF
            bdim = ppm_buffer_dim(k)
            DO j=2,ppm_nsendlist
               pp(j,k) = pp(j-1,k) + precv(j-1)*bdim
            ENDDO
         ENDDO
      ENDIF

      !counter for number of send & recv
      nsr=0
#endif
      IF (ppm_kind.EQ.ppm_kind_double) THEN
#ifdef __MPI
         DO k=1,ppm_buffer_set
            !----------------------------------------------------------------------
            !  For each processor
            !----------------------------------------------------------------------
            bdim = ppm_buffer_dim(k)
            l=qq(1,k)-1
            m=pp(1,k)-1
            DO j=1,psend(1)*bdim
               ppm_recvbufferd(m+j)=ppm_sendbufferd(l+j)
            ENDDO
         ENDDO !k=1,ppm_buffer_set

         DO k=1,ppm_buffer_set
            bdim = ppm_buffer_dim(k)
            !-------------------------------------------------------------------------
            !  loop over the processors in the ppm_isendlist(); skip the first entry
            !  which is the local processor
            !-------------------------------------------------------------------------
            !----------------------------------------------------------------------
            !  For each send/recv
            !----------------------------------------------------------------------
            DO i=2,ppm_nsendlist
               !-------------------------------------------------------------------
               !  Perform the actual send/recv
               !-------------------------------------------------------------------
               ! The following IF is needed in order to skip "dummy"
               ! communication rounds where the current processor has to wait
               ! (only needed in the partial mapping).
               IF (ppm_isendlist(i).LT.0.OR.ppm_irecvlist(i).LT.0) CYCLE
               IF (precv(i).GT.0) THEN
                  IF (ppm_rank.GT.ppm_irecvlist(i)) THEN
                     tag=ppm_rank*(ppm_rank-1)/2+ppm_irecvlist(i)+k
                  ELSE
                     tag=ppm_irecvlist(i)*(ppm_irecvlist(i)-1)/2+ppm_rank+k
                  ENDIF

                  nsr=nsr+1
                  CALL MPI_Irecv(ppm_recvbufferd(pp(i,k)),precv(i)*bdim, &
                  &    ppm_mpi_kind,ppm_irecvlist(i),tag,ppm_comm,request(nsr),info)
                  or_fail_MPI("MPI_Irecv")
               ENDIF
               IF (psend(i).GT.0) THEN
                  IF (ppm_rank.GT.ppm_isendlist(i)) THEN
                     tag=ppm_rank*(ppm_rank-1)/2+ppm_isendlist(i)+k
                  ELSE
                     tag=ppm_isendlist(i)*(ppm_isendlist(i)-1)/2+ppm_rank+k
                  ENDIF

                  nsr=nsr+1
                  CALL MPI_Isend(ppm_sendbufferd(qq(i,k)),psend(i)*bdim, &
                  &    ppm_mpi_kind,ppm_isendlist(i),tag,ppm_comm,request(nsr),info)
                  or_fail_MPI("MPI_Isend")
               ENDIF
            ENDDO !i=2,ppm_nsendlist
         ENDDO !k=1,ppm_buffer_set
#else
         ppm_recvbufferd=ppm_sendbufferd
#endif
      ELSE
#ifdef __MPI
         DO k=1,ppm_buffer_set
            !----------------------------------------------------------------------
            !  For each processor
            !----------------------------------------------------------------------
            bdim = ppm_buffer_dim(k)
            l=qq(1,k)-1
            m=pp(1,k)-1
            DO j=1,psend(1)*bdim
               ppm_recvbuffers(m+j)=ppm_sendbuffers(l+j)
            ENDDO
         ENDDO !k=1,ppm_buffer_set

         DO k=1,ppm_buffer_set
            bdim = ppm_buffer_dim(k)
            !-------------------------------------------------------------------------
            !  loop over the processors in the ppm_isendlist(); skip the first entry
            !  which is the local processor
            !-------------------------------------------------------------------------
            !----------------------------------------------------------------------
            !  For each send/recv
            !----------------------------------------------------------------------
            DO i=2,ppm_nsendlist
               !-------------------------------------------------------------------
               !  Perform the actual send/recv
               !-------------------------------------------------------------------
               ! The following IF is needed in order to skip "dummy"
               ! communication rounds where the current processor has to wait
               ! (only needed in the partial mapping).
               IF (ppm_isendlist(i).LT.0.OR.ppm_irecvlist(i).LT.0) CYCLE
               IF (precv(i).GT.0) THEN
                  IF (ppm_rank.GT.ppm_irecvlist(i)) THEN
                     tag=ppm_rank*(ppm_rank-1)/2+ppm_irecvlist(i)+k
                  ELSE
                     tag=ppm_irecvlist(i)*(ppm_irecvlist(i)-1)/2+ppm_rank+k
                  ENDIF

                  nsr=nsr+1
                  CALL MPI_Irecv(ppm_recvbuffers(pp(i,k)),precv(i)*bdim, &
                  &    ppm_mpi_kind,ppm_irecvlist(i),tag,ppm_comm,request(nsr),info)
                  or_fail_MPI("MPI_Irecv")
               ENDIF
               IF (psend(i).GT.0) THEN
                  IF (ppm_rank.GT.ppm_isendlist(i)) THEN
                     tag=ppm_rank*(ppm_rank-1)/2+ppm_isendlist(i)+k
                  ELSE
                     tag=ppm_isendlist(i)*(ppm_isendlist(i)-1)/2+ppm_rank+k
                  ENDIF

                  nsr=nsr+1
                  CALL MPI_Isend(ppm_sendbuffers(qq(i,k)),psend(i)*bdim, &
                  &    ppm_mpi_kind,ppm_isendlist(i),tag,ppm_comm,request(nsr),info)
                  or_fail_MPI("MPI_Isend")
               ENDIF
            ENDDO !i=2,ppm_nsendlist
         ENDDO !k=1,ppm_buffer_set
#else
         ppm_recvbuffers=ppm_sendbuffers
#endif
      ENDIF

      !-------------------------------------------------------------------------
      !  before we through away the precv() data let us store it for later use:
      !  when sending ghosts back (ppm_map_part_ghost_put())
      !-------------------------------------------------------------------------
      IF (ppm_map_type.EQ.ppm_param_map_ghost_get) THEN
         ldu(1) = ppm_nsendlist + 1
         iopt   = ppm_param_alloc_grow
         CALL ppm_alloc(ppm_precvbuffer,ldu,iopt,info)
         or_fail_alloc('global recv buffer pointer PPM_PRECVBUFFER',ppm_error=ppm_error_fatal)

         ppm_precvbuffer(1) = Npart + 1
         DO k=1,ppm_nsendlist
            ppm_precvbuffer(k+1) = ppm_precvbuffer(k) + precv(k)
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
      SELECT CASE (ppm_kind)
      CASE (ppm_kind_single)
         CALL ppm_alloc(ppm_sendbuffers,ldu,iopt,info)

      CASE DEFAULT
         CALL ppm_alloc(ppm_sendbufferd,ldu,iopt,info)

      END SELECT
      or_fail_dealloc('send buffer PPM_SENDBUFFER',exit_point=no)

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
        IF (Npart.LT.0) THEN
           fail('Npart must be >=0',exit_point=8888)
        ENDIF
      8888 CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_map_part_isend
