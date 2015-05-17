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
      !!! [TIP] Yaser:
      !!! This is not the best implementation
      !!! I could not figure out the general derive type to use for MPI_Isend &
      !!! MPI_Irecv to avoid several mpi call, so this implementation will call
      !!! mpi for each ppm_buffer_set, but using nonblocking send & recv it will
      !!! prevent several buffer copy which is very useful on big data !
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_mpi
      USE ppm_module_mapping_typedef, ONLY : ppm_mesh_isendblksize, &
      &   ppm_mesh_irecvblksize
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
      INTEGER, DIMENSION(3)              :: ldu
      INTEGER                            :: i,j,k,l,m,nsr,tag
      INTEGER                            :: bdim,offs
      INTEGER                            :: iopt,Ndata
      INTEGER                            :: allsend,allrecv
#ifdef __MPI
      INTEGER, DIMENSION(:), ALLOCATABLE :: request
#endif
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------
      start_subroutine("mesh_map_isend")

      ! warn if buffer is empty
      IF (ppm_buffer_set.LT.1) THEN
         IF (ppm_debug.GT.1) THEN
            fail('Buffer is empty: skipping send!', &
            & ppm_err_buffer_empt,exit_point=no,    &
            & ppm_error=ppm_error_notice)
            info = 0
         ENDIF
         GOTO 9999
      ELSE
         IF (ppm_debug.GT.1) THEN
            stdout_f('(A,I9)',"ppm_buffer_set=",ppm_buffer_set)
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Initialize the buffer counters
      !-------------------------------------------------------------------------
      ppm_nrecvbuffer = 0

      !-------------------------------------------------------------------------
      !  Allocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldu(1) = ppm_nrecvlist
      CALL ppm_alloc(precv,ldu,iopt,info)
      or_fail_alloc("precv")

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
      !  Sum of all mesh points that will be received
      !-------------------------------------------------------------------------
      allrecv = SUM(precv(1:ppm_nrecvlist))

      !-------------------------------------------------------------------------
      !  Allocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldu(1) = ppm_nrecvlist
      ldu(2) = ppm_buffer_set
      CALL ppm_alloc(pp,ldu,iopt,info)
      or_fail_alloc("pp")

      !-------------------------------------------------------------------------
      !  Compute the pointer to the position of the data in the main receive
      !  buffer
      !-------------------------------------------------------------------------
      bdim = 0
      offs = 0
      IF (ppm_debug.GT.1) THEN
         DO k=1,ppm_buffer_set
            offs    = offs + allrecv*bdim
            pp(1,k) = offs + 1
            bdim    = ppm_buffer_dim(k)
            DO j=2,ppm_nrecvlist
               pp(j,k) = pp(j-1,k) + precv(j-1)*bdim
               stdout_f('(A,I9)',"pp(j,k)=",'pp(j,k)')
               stdout_f('(A,I9,A,I4)',"precv(j-1)=",'precv(j-1)',", bdim=",bdim)
            ENDDO
         ENDDO
      ELSE
         DO k=1,ppm_buffer_set
            offs    = offs + allrecv*bdim
            pp(1,k) = offs + 1
            bdim    = ppm_buffer_dim(k)
            DO j=2,ppm_nrecvlist
               pp(j,k) = pp(j-1,k) + precv(j-1)*bdim
            ENDDO
         ENDDO
      ENDIF

      ALLOCATE(request(ppm_nproc*ppm_buffer_set*2),STAT=info)
      or_fail_alloc("request")

      !counter for number of send & recv
      nsr=0
#endif

      IF (ppm_kind.EQ.ppm_kind_double) THEN
#ifdef __MPI
         DO k=1,ppm_buffer_set
            bdim = ppm_buffer_dim(k)
            DO i=2,ppm_nsendlist
               !----------------------------------------------------------------------
               ! The following IF is needed in order to skip "dummy"
               ! communication rounds where the current processor has to wait
               ! (only needed in the partial mapping).
               !----------------------------------------------------------------------
               IF (ppm_isendlist(i).LT.0.OR.ppm_irecvlist(i).LT.0) CYCLE
               IF (precv(i).GT.0) THEN
                  IF (ppm_rank.GT.ppm_irecvlist(i)) THEN
                     tag=ppm_rank*(ppm_rank-1)/2+ppm_irecvlist(i)+k
                  ELSE
                     tag=ppm_irecvlist(i)*(ppm_irecvlist(i)-1)/2+ppm_rank+k
                  ENDIF

                  nsr=nsr+1
                  CALL MPI_Irecv(ppm_recvbufferd(pp(i,k)),precv(i)*bdim,ppm_mpi_kind, &
                  &    ppm_irecvlist(i),tag,ppm_comm,request(nsr),info)
                  or_fail_MPI("MPI_Irecv")
               ENDIF
            ENDDO !i=2,ppm_nsendlist
         ENDDO !k=1,ppm_buffer_set
#else
         ppm_recvbufferd=ppm_sendbufferd
#endif
      ELSE
#ifdef __MPI
         DO k=1,ppm_buffer_set
            bdim = ppm_buffer_dim(k)
            DO i=2,ppm_nsendlist
               !----------------------------------------------------------------------
               ! The following IF is needed in order to skip "dummy"
               ! communication rounds where the current processor has to wait
               ! (only needed in the partial mapping).
               !----------------------------------------------------------------------
               IF (ppm_isendlist(i).LT.0.OR.ppm_irecvlist(i).LT.0) CYCLE
               IF (precv(i).GT.0) THEN
                  IF (ppm_rank.GT.ppm_irecvlist(i)) THEN
                     tag=ppm_rank*(ppm_rank-1)/2+ppm_irecvlist(i)+k
                  ELSE
                     tag=ppm_irecvlist(i)*(ppm_irecvlist(i)-1)/2+ppm_rank+k
                  ENDIF

                  nsr=nsr+1
                  CALL MPI_Irecv(ppm_recvbuffers(pp(i,k)),precv(i)*bdim,ppm_mpi_kind, &
                  &    ppm_irecvlist(i),tag,ppm_comm,request(nsr),info)
                  or_fail_MPI("MPI_Irecv")
               ENDIF
            ENDDO
         ENDDO !k=1,ppm_buffer_set
#else
         ppm_recvbuffers=ppm_sendbuffers
#endif
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(psend,ldu,iopt,info)
      or_fail_alloc("psend")

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

#ifdef __MPI
      !-------------------------------------------------------------------------
      !  Sum of all mesh points that will be sent
      !-------------------------------------------------------------------------
      allsend = SUM(psend(1:ppm_nsendlist))

      !-------------------------------------------------------------------------
      !  Allocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      ldu(2) = ppm_buffer_set
      CALL ppm_alloc(qq,ldu,iopt,info)
      or_fail_alloc("qq")

      !-------------------------------------------------------------------------
      !  Compute the pointer to the position of the data in the main send
      !  buffer
      !-------------------------------------------------------------------------
      bdim = 0
      offs = 0
      IF (ppm_debug.GT.1) THEN
         DO k=1,ppm_buffer_set
            offs    = offs + allsend*bdim
            qq(1,k) = offs + 1
            bdim    = ppm_buffer_dim(k)
            DO j=2,ppm_nsendlist
               qq(j,k) = qq(j-1,k) + psend(j-1)*bdim
               stdout_f('(A,I9)',"qq(j,k)=",'qq(j,k)')
            ENDDO
         ENDDO
      ELSE
         DO k=1,ppm_buffer_set
            offs    = offs + allsend*bdim
            qq(1,k) = offs + 1
            bdim    = ppm_buffer_dim(k)
            DO j=2,ppm_nsendlist
               qq(j,k) = qq(j-1,k) + psend(j-1)*bdim
            ENDDO
         ENDDO
      ENDIF

      IF (ppm_kind.EQ.ppm_kind_double) THEN
         DO k=1,ppm_buffer_set
            bdim = ppm_buffer_dim(k)
            !----------------------------------------------------------------------
            !  For each processor
            !----------------------------------------------------------------------
            l=qq(1,k)-1
            m=pp(1,k)-1
            DO j=1,psend(1)*bdim
               ppm_recvbufferd(m+j)=ppm_sendbufferd(l+j)
            ENDDO
         ENDDO !k=1,ppm_buffer_set

         DO k=1,ppm_buffer_set
            bdim = ppm_buffer_dim(k)
            DO i=2,ppm_nsendlist
               !----------------------------------------------------------------------
               ! The following IF is needed in order to skip "dummy"
               ! communication rounds where the current processor has to wait
               ! (only needed in the partial mapping).
               !----------------------------------------------------------------------
               IF (ppm_isendlist(i).LT.0.OR.ppm_irecvlist(i).LT.0) CYCLE
               IF (psend(i).GT.0) THEN
                  IF (ppm_rank.GT.ppm_isendlist(i)) THEN
                     tag=ppm_rank*(ppm_rank-1)/2+ppm_isendlist(i)+k
                  ELSE
                     tag=ppm_isendlist(i)*(ppm_isendlist(i)-1)/2+ppm_rank+k
                  ENDIF

                  nsr=nsr+1
                  CALL MPI_Isend(ppm_sendbufferd(qq(i,k)),psend(i)*bdim,ppm_mpi_kind, &
                  &    ppm_isendlist(i),tag,ppm_comm,request(nsr),info)
                  or_fail_MPI("MPI_Isend")
               ENDIF
            ENDDO !i=2,ppm_nsendlist
         ENDDO !k=1,ppm_buffer_set
      ELSE
         DO k=1,ppm_buffer_set
            bdim = ppm_buffer_dim(k)
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
            DO i=2,ppm_nsendlist
               !----------------------------------------------------------------------
               ! The following IF is needed in order to skip "dummy"
               ! communication rounds where the current processor has to wait
               ! (only needed in the partial mapping).
               !----------------------------------------------------------------------
               IF (ppm_isendlist(i).LT.0.OR.ppm_irecvlist(i).LT.0) CYCLE
               IF (psend(i).GT.0) THEN
                  IF (ppm_rank.GT.ppm_isendlist(i)) THEN
                     tag=ppm_rank*(ppm_rank-1)/2+ppm_isendlist(i)+k
                  ELSE
                     tag=ppm_isendlist(i)*(ppm_isendlist(i)-1)/2+ppm_rank+k
                  ENDIF

                  nsr=nsr+1
                  CALL MPI_Isend(ppm_sendbuffers(qq(i,k)),psend(i)*bdim,ppm_mpi_kind, &
                  &    ppm_isendlist(i),tag,ppm_comm,request(nsr),info)
                  or_fail_MPI("MPI_Isend")
               ENDIF
            ENDDO
         ENDDO !k=1,ppm_buffer_set
      ENDIF
#endif

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

      CALL ppm_alloc(psend,ldu,iopt,info)
      or_fail_dealloc("psend")

      CALL ppm_alloc(precv,ldu,iopt,info)
      or_fail_dealloc("precv")

#ifdef __MPI
      CALL ppm_alloc(   pp,ldu,iopt,info)
      or_fail_dealloc("pp")

      CALL ppm_alloc(   qq,ldu,iopt,info)
      or_fail_dealloc("qq")

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
      IF (ppm_kind.EQ.ppm_kind_single) THEN
         CALL ppm_alloc(ppm_sendbuffers,ldu,iopt,info)
         or_fail_dealloc("ppm_sendbuffer")
      ELSE
         CALL ppm_alloc(ppm_sendbufferd,ldu,iopt,info)
         or_fail_dealloc("ppm_sendbuffer")
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      end_subroutine()
      RETURN
      END SUBROUTINE equi_mesh_map_isend
