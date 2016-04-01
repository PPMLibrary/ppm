      !--*- f90 -*--------------------------------------------------------------
      !  Subroutine   :                 ppm_netstat
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

      SUBROUTINE ppm_netstat(topoid,latency,bandwidth,info)
      !!! This routine provides a simple means to measure the total latency and
      !!! mean bandwidth using the current communication schedule

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------

      USE ppm_module_data
      USE ppm_module_topo
      USE ppm_module_util_commopt
      USE ppm_module_util_time
      USE ppm_module_alloc
      USE ppm_module_write
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_mpi
      IMPLICIT NONE

      INTEGER, PARAMETER :: MK = ppm_kind_double
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER,  INTENT(IN   ) :: topoid
      !!! topology ID to which we are currently mapped
      REAL(MK), INTENT(  OUT) :: latency
      !!! Network latency in seconds
      REAL(MK), INTENT(  OUT) :: bandwidth
      !!! Network bandwidth in Bytes per second
      INTEGER,  INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  local vars
      !-------------------------------------------------------------------------
      TYPE(ppm_t_topo), POINTER :: topo

      REAL(ppm_kind_double)           :: t0
      REAL(MK)                        :: tstart,tend
      REAL(MK), DIMENSION(:), POINTER :: latencies
      REAL(MK), DIMENSION(:), POINTER :: bandwidths

      INTEGER                             :: i,irank,icomm,tag
      INTEGER                             :: iopt
      INTEGER                             :: src,dest
      INTEGER, PARAMETER                  :: lncomm = 1000 !N latency comms
      INTEGER, PARAMETER                  :: bncomm = 100  !N bandwidth comms
      INTEGEr, DIMENSION(1)               :: lda
      INTEGER, DIMENSION(:), POINTER      :: buf
      INTEGER, PARAMETER                  :: small = 1
      INTEGER, PARAMETER                  :: big   = 8*(1024**2)

      CHARACTER(LEN=ppm_char) :: caller = 'ppm_netstat'

      CALL substart(caller,t0,info)

      topo => ppm_topo(topoid)%t

      !----------------------------------------------------------------------
      !  first check if the optimal communication protocol is known
      !----------------------------------------------------------------------
      IF (.NOT.topo%isoptimized) THEN
         !-------------------------------------------------------------------
         !  if not: determine it
         !-------------------------------------------------------------------
         CALL ppm_util_commopt(topoid,info)
         IF (ppm_debug.GT.1) THEN
            DO i=1,topo%nneighproc
               stdout_f('(A,I4)',"have neighbor: ",'topo%ineighproc(i)')
            ENDDO
            DO i=1,topo%ncommseq
               stdout_f('(A,I4)',"communicate: ",'topo%icommseq(i)')
            ENDDO
         ENDIF
      ENDIF
      !-------------------------------------------------------------------------
      !  Allocate result arrays
      !-------------------------------------------------------------------------
      NULLIFY(latencies,bandwidths)
      iopt   = ppm_param_alloc_fit
      lda(1) = topo%ncommseq
      CALL ppm_alloc(latencies,lda,iopt,info)
      or_fail_alloc('latencies')

      CALL ppm_alloc(bandwidths,lda,iopt,info)
      or_fail_alloc('bandwidths')

      latencies(:) = 0.0_MK
      bandwidths(:) = 0.0_MK

      IF (ppm_nproc.EQ.1) THEN
         fail('Cannot compute latency/bandwidth on 1 proc',ppm_err_test_fail, &
         & exit_point=8888,ppm_error=ppm_error_warning)
      ENDIF
#ifdef __MPI
      !-------------------------------------------------------------------------
      !  Measure Latency
      !-------------------------------------------------------------------------
      ALLOCATE(buf(small),STAT=info)
      or_fail_alloc('latency buffer')

      buf(:) = 42

      DO icomm=1,topo%ncommseq
         irank = topo%icommseq(icomm)
         tag = 100
         dest = MAX(ppm_rank,irank)
         src = MIN(ppm_rank,irank)

         CALL MPI_Barrier(ppm_comm,info)
         or_fail_MPI("MPI_Barrier")

         CALL ppm_util_time(tstart)

         DO i=1,lncomm
            IF (ppm_rank.EQ.src) THEN
               CALL MPI_Send(buf,small,MPI_INTEGER,dest,i,ppm_comm,info)
               or_fail_MPI("MPI Send failed!")
            ELSE
               CALL MPI_Recv(buf,small,MPI_INTEGER,src,i,ppm_comm,MPI_STATUS_IGNORE,info)
               or_fail_MPI("MPI Recv failed!")
            ENDIF

!           CALL MPI_SendRecv(sbuf,small,MPI_INTEGER,irank,tag, &
!           &    rbuf,small,MPI_INTEGER,irank,tag,ppm_comm,MPI_STATUS_IGNORE,info)

         ENDDO

         CALL ppm_util_time(tend)

         latencies(icomm) = (tend-tstart)/(2*lncomm)
      ENDDO


      DEALLOCATE(buf,STAT=info)
      or_fail_dealloc('dealloc latency buffer')

      !-------------------------------------------------------------------------
      !  Measure Bandwidth
      !-------------------------------------------------------------------------
      ALLOCATE(buf(big),STAT=info)
      or_fail_alloc('bandwidth buffer')

      buf(:) = 42

      DO icomm=1,topo%ncommseq
         irank = topo%icommseq(icomm)
         tag = 200

         CALL MPI_Barrier(ppm_comm,info)
         or_fail_MPI("MPI_Barrier")

         CALL ppm_util_time(tstart)

         DO i=1,bncomm
            IF (ppm_rank.EQ.src) THEN
               CALL MPI_Send(buf,big,MPI_INTEGER,dest,i,ppm_comm,info)
               or_fail_MPI("MPI Send failed!")
            ELSE
               CALL MPI_Recv(buf,big,MPI_INTEGER,src,i,ppm_comm,MPI_STATUS_IGNORE,info)
               or_fail_MPI("MPI Recv failed!")
            ENDIF
!           CALL MPI_SendRecv(sbuf,big,MPI_INTEGER,irank,tag, &
!           &    rbuf,big,MPI_INTEGER,irank,tag,ppm_comm,MPI_STATUS_IGNORE,info)

         ENDDO

         CALL ppm_util_time(tend)

         bandwidths(icomm) = big/(((tend-tstart)-latencies(icomm))/bncomm)
      ENDDO


      DEALLOCATE(buf,STAT=info)
      or_fail_dealloc("dealloc bandwidth buffer")
#else
      fail('Bandwidth and latency default to 0',ppm_err_nompi,exit_point=no, &
      & ppm_error=ppm_error_warning)
#endif
      !-------------------------------------------------------------------------
      !  Computer summary statistics
      !-------------------------------------------------------------------------
      8888 CONTINUE

      bandwidth = SUM(bandwidths(1:topo%ncommseq))
      bandwidth = 4*bandwidth/topo%ncommseq

      latency = SUM(latencies(1:topo%ncommseq))

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)

      END SUBROUTINE ppm_netstat
