      !-------------------------------------------------------------------------
      !  Subroutine   :                 equi_mesh_map_isend
      !-------------------------------------------------------------------------
      ! Copyright (c) 2012 CSE Lab (ETH Zurich), MOSAIC Group (Max Planck
      ! Institute of Molecular Cell Biology and Genetics), Center for Fluid
      ! Dynamics (DTU)
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
      SUBROUTINE equi_mesh_map_isend(this,info)
      !!! This routine performs the actual send/recv of the
      !!! mesh blocks and all pushed data.
      !!! This is the simplest non-blocking version of `equi_mesh_map_send`
      !!! using Iscatterv.
      !!!
      !!! [TIP]
      !!! If you want to only send certain meshnodes you have to first call
      !!! ppm_map_field_push for a logical mask, then call ppm_map_field_pop
      !!! and recover the new mask. First 3 (2 for 2D meshes)
      !!! indices are mesh (i,j[,k]), last one is subid for all subs on
      !!! the local processor.
      !!!
      !!! [WARNING] Yaser:
      !!! This is not the best implementation
      !!! I could not figure out the general derive type to use for MPI_Iscatterv
      !!! to avoid several mpi call, so this implementation will call mpi for each
      !!! ppm_buffer_set, but using MPI_Iscatterv it will prevent several buffer
      !!! copy which is very useful on big data !
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
      INTEGER,               INTENT(  OUT) :: info
      !!! Return status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3)                :: ldu
      INTEGER                              :: i,j,k,bdim,offs
      INTEGER                              :: iopt,Ndata
      INTEGER                              :: allsend,allrecv
#ifdef __MPI
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: request
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: senddispls,recvdispls
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: sendcounts,recvcounts
#endif
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------
      start_subroutine("mesh_map_isend")

      ! warn if buffer is empty
      IF (ppm_buffer_set.LT.1) THEN
         IF (ppm_debug.GT.1) THEN
            fail('Buffer is empty: skipping send!',ppm_err_buffer_empt,exit_point=no,ppm_error=ppm_error_notice)
            info = 0
         ENDIF
         GOTO 9999
      ELSE
         IF (ppm_debug.GT.1) THEN
            stdout_f('(A,I9)',"ppm_buffer_set=",ppm_buffer_set)
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(psend,ldu,iopt,info)
      or_fail_alloc("psend")

      ldu(1) = ppm_nrecvlist
      CALL ppm_alloc(precv,ldu,iopt,info)
      or_fail_alloc("precv")

#ifdef __MPI
      ldu(1) = ppm_nrecvlist
      ldu(2) = ppm_buffer_set
      CALL ppm_alloc(pp,ldu,iopt,info)
      or_fail_alloc("pp")

      ldu(1) = ppm_nsendlist
      ldu(2) = ppm_buffer_set
      CALL ppm_alloc(qq,ldu,iopt,info)
      or_fail_alloc("qq")
#endif
      !-------------------------------------------------------------------------
      !  Initialize the buffer counters
      !-------------------------------------------------------------------------
      ppm_nrecvbuffer = 0

      !-------------------------------------------------------------------------
      !  Count the size of the buffers that will be sent
      !-------------------------------------------------------------------------
      DO k=1,ppm_nsendlist
         Ndata = 0

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
      ENDDO

      !-------------------------------------------------------------------------
      !  Count the size of the buffers that will be received
      !  (For meshes we know this a priori. This is different from
      !  particles)
      !-------------------------------------------------------------------------
      DO k=1,ppm_nrecvlist
         Ndata = 0

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
         !  Increment the total receive buffer count
         !---------------------------------------------------------------------
         ppm_nrecvbuffer = ppm_nrecvbuffer + SUM(ppm_buffer_dim(1:ppm_buffer_set))*Ndata
      ENDDO

      IF (ppm_debug.GT.1) THEN
         stdout_f('(A,I9)',"ppm_nrecvbuffer = ",ppm_nrecvbuffer)
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate the memory for the copy of the buffer
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow
      ldu(1) = ppm_nrecvbuffer
      IF (ppm_kind.EQ.ppm_kind_double) THEN
         CALL ppm_alloc(ppm_recvbufferd,ldu,iopt,info)
      ELSE
         CALL ppm_alloc(ppm_recvbuffers,ldu,iopt,info)
      ENDIF
      or_fail_alloc("global receive buffer PPM_RECVBUFFER")

#ifdef __MPI
      !-------------------------------------------------------------------------
      !  Sum of all mesh points that will be sent and received
      !-------------------------------------------------------------------------
      allsend = SUM(psend(1:ppm_nsendlist))
      allrecv = SUM(precv(1:ppm_nrecvlist))

      !-------------------------------------------------------------------------
      !  Compute the pointer to the position of the data in the main send
      !  buffer
      !-------------------------------------------------------------------------
      bdim = 0
      offs = 0
      DO k=1,ppm_buffer_set
         offs    = offs + allsend*bdim
         qq(1,k) = offs + 1
         bdim    = ppm_buffer_dim(k)
         DO j=2,ppm_nsendlist
            qq(j,k) = qq(j-1,k) + psend(j-1)*bdim
            IF (ppm_debug.GT.1) THEN
               stdout_f('(A,I9)',"qq(j,k)=",'qq(j,k)')
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
            IF (ppm_debug.GT.1) THEN
               stdout_f('(A,I9)',"pp(j,k)=",'pp(j,k)')
               stdout_f('(A,I9,A,I4)',"precv(j-1)=",'precv(j-1)',", bdim=",bdim)
            ENDIF
         ENDDO
      ENDDO

      ALLOCATE(senddispls(ppm_nsendlist,ppm_buffer_set), &
      &        recvdispls(ppm_nrecvlist,ppm_buffer_set), &
      &        sendcounts(ppm_nsendlist,ppm_buffer_set), &
      &        recvcounts(ppm_nrecvlist,ppm_buffer_set),STAT=info)
      or_fail_alloc("senddispls, recvdispls, sendcounts & recvcounts")

      ALLOCATE(request(ppm_nproc,ppm_buffer_set),STAT=info)
      or_fail_alloc("request")

      !-------------------------------------------------------------------------
      ! Compute real processor displacements and counts
      !-------------------------------------------------------------------------
      DO k=1,ppm_buffer_set
         bdim = ppm_buffer_dim(k)
         DO i=1,ppm_nsendlist
            j=ppm_isendlist(i)+1
            senddispls(j,k)=qq(i,k)-1
            sendcounts(j,k)=psend(i)*bdim
         ENDDO
         DO i=1,ppm_nrecvlist
            j=ppm_irecvlist(i)+1
            recvdispls(j,k)=pp(i,K)
            recvcounts(j,k)=precv(i)*bdim
         ENDDO
      ENDDO

#endif
      IF (ppm_kind.EQ.ppm_kind_double) THEN
#ifdef __MPI
         DO k=1,ppm_buffer_set
            !----------------------------------------------------------------------
            !  For each processor
            !  The outcome is as if the root executed n send operations for each
            !  ppm_buffer_set.
            !  The send buffer is ignored for all the nonroot processes.
            !  All arguments to the function are significant on process root,
            !  while on the other processes, only arguments ppm_recvbufferd,
            !  recvcounts, ppm_mpi_kind, root, ppm_comm are significant.
            !----------------------------------------------------------------------
            DO i=1,ppm_nproc
               IF (i-1.NE.ppm_rank.AND.recvcounts(i,k).EQ.0) THEN
                  request(i,k)=MPI_SUCCESS
               ELSE
                  CALL MPI_Iscatterv(ppm_sendbufferd,sendcounts(:,k),senddispls(:,k), &
                  &    ppm_mpi_kind,ppm_recvbufferd(recvdispls(i,k)),recvcounts(i,k), &
                  &    ppm_mpi_kind,i-1,ppm_comm,request(i,k),info)
                  or_fail_MPI("MPI_Iscatterv")
               ENDIF
            ENDDO
         ENDDO
#else
         ppm_recvbufferd=ppm_sendbufferd
#endif
      ELSE
#ifdef __MPI
         DO k=1,ppm_buffer_set
            !----------------------------------------------------------------------
            !  For each processor
            !  The outcome is as if the root executed n send operations for each
            !  ppm_buffer_set.
            !  The send buffer is ignored for all the nonroot processes.
            !  All arguments to the function are significant on process root,
            !  while on the other processes, only arguments ppm_recvbuffers,
            !  recvcounts, ppm_mpi_kind, root, ppm_comm are significant.
            !----------------------------------------------------------------------
            DO i=1,ppm_nproc
               IF (i-1.NE.ppm_rank.AND.recvcounts(i,k).EQ.0) THEN
                  request(i,k)=MPI_SUCCESS
               ELSE
                  CALL MPI_Iscatterv(ppm_sendbuffers,sendcounts(:,k),senddispls(:,k), &
                  &    ppm_mpi_kind,ppm_recvbuffers(recvdispls(i,k)),recvcounts(i,k), &
                  &    ppm_mpi_kind,i-1,ppm_comm,request(i,k),info)
                  or_fail_MPI("MPI_Iscatterv")
               ENDIF
            ENDDO
         ENDDO
#else
         ppm_recvbuffers=ppm_sendbuffers
#endif
      ENDIF

      !-------------------------------------------------------------------------
      !  Deallocate the send buffer to save memory
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      IF (ppm_kind.EQ.ppm_kind_single) THEN
         CALL ppm_alloc(recvs,ldu,iopt,info)
         or_fail_dealloc("recvs")

         CALL ppm_alloc(sends,ldu,iopt,info)
         or_fail_dealloc("sends")
      ELSE
         CALL ppm_alloc(recvd,ldu,iopt,info)
         or_fail_dealloc("recvd")

         CALL ppm_alloc(sendd,ldu,iopt,info)
         or_fail_dealloc("sendd")
      ENDIF

      !-------------------------------------------------------------------------
      !  Deallocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(nsend,ldu,iopt,info)
      or_fail_dealloc("nsend")
      CALL ppm_alloc(nrecv,ldu,iopt,info)
      or_fail_dealloc("nrecv")

#ifdef __MPI
      k=ppm_buffer_set*ppm_nproc
      CALL MPI_Waitall(k,RESHAPE(request,(/k/)),MPI_STATUSES_IGNORE,info)
      or_fail_MPI("MPI_Waitall")
#endif

      !-------------------------------------------------------------------------
      !  Deallocate the send buffer to save memory
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      IF (ppm_kind.EQ.ppm_kind_single) THEN
         CALL ppm_alloc(ppm_sendbuffers,ldu,iopt,info)
         or_fail_dealloc("ppm_sendbuffer")
      ELSE
         CALL ppm_alloc(ppm_sendbufferd,ldu,iopt,info)
         or_fail_dealloc("ppm_sendbuffer")
      ENDIF

      !-------------------------------------------------------------------------
      !  Deallocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(psend,ldu,iopt,info)
      or_fail_dealloc("psend")
      CALL ppm_alloc(precv,ldu,iopt,info)
      or_fail_dealloc("precv")
#ifdef __MPI
      CALL ppm_alloc(   pp,ldu,iopt,info)
      or_fail_dealloc("pp")
      CALL ppm_alloc(   qq,ldu,iopt,info)
      or_fail_dealloc("qq")

      DEALLOCATE(senddispls,recvdispls,sendcounts,recvcounts,request,STAT=info)
      or_fail_dealloc("senddispls,recvdispls,sendcounts,recvcounts & request")
#endif
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      end_subroutine()
      RETURN
      END SUBROUTINE equi_mesh_map_isend
