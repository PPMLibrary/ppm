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

      ! arguments
      INTEGER,            INTENT(IN)    :: topoid
      !!! topology ID to which we are currently mapped
      REAL(MK),           INTENT(OUT)   :: latency
      !!! Network latency in seconds
      REAL(MK),           INTENT(OUT)   :: bandwidth
      !!! Network bandwidth in Bytes per second
      INTEGER,            INTENT(OUT)   :: info

      ! local vars
      TYPE(ppm_t_topo), POINTER         :: topo
      INTEGER                           :: i,irank,icomm,tag
      REAL(MK)                          :: t0
      REAL(MK)                          :: tstart,tend
      REAL(MK), DIMENSION(:), POINTER   :: latencies => NULL()
      REAL(MK), DIMENSION(:), POINTER   :: bandwidths => NULL()
      INTEGER                           :: iproc
#ifdef __MPI
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: stat
#endif
      INTEGER                           :: iopt
      INTEGER                           :: src,dest
      INTEGER, PARAMETER                :: lncomm = 1000 !N latency comms
      INTEGER, PARAMETER                :: bncomm = 100  !N bandwidth comms
      INTEGEr, DIMENSION(1)             :: lda
      INTEGER, DIMENSION(:), POINTER    :: buf
!      INTEGER, DIMENSION(:), POINTER    :: rbuf => NULL()
      INTEGER, PARAMETER                :: small = 1
      INTEGER, PARAMETER                :: big   = 8*(1024**2)
      CHARACTER(LEN=32), PARAMETER      :: caller = 'ppm_netstat'
      CHARACTER(LEN=32)                 :: mesg

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
        IF (ppm_debug .GT. 1) THEN
            DO i=1,topo%nneighproc
                WRITE(mesg,'(A,I4)') 'have neighbor: ', topo%ineighproc(i)
                CALL ppm_write(ppm_rank,caller,mesg,info)
            ENDDO
            DO i=1,topo%ncommseq
                WRITE(mesg,'(A,I4)') 'communicate: ', topo%icommseq(i)
                CALL ppm_write(ppm_rank,caller,mesg,info)
            ENDDO
        ENDIF
      ENDIF
      !-------------------------------------------------------------------------
      !  Allocate result arrays
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      lda(1) = topo%ncommseq
      CALL ppm_alloc(latencies,lda,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_alloc,caller,     &
     &        'latencies',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(bandwidths,lda,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_alloc,caller,     &
     &        'bandwidths',__LINE__,info)
          GOTO 9999
      ENDIF
      latencies(:) = 0.0_MK
      bandwidths(:) = 0.0_MK

      IF (ppm_nproc.EQ.1) THEN
          info = ppm_error_warning
          CALL ppm_error(ppm_err_test_fail,caller,     &
     &        'Cannot compute latency/bandwidth on 1 proc',__LINE__,info)
          GOTO 8888
      ENDIF
#ifdef __MPI
      !-------------------------------------------------------------------------
      !  Measure Latency
      !-------------------------------------------------------------------------
      ALLOCATE(buf(small),STAT=info)
      IF (info.NE.0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_alloc,caller,     &
     &        'latency buffer',__LINE__,info)
          GOTO 9999
      ENDIF
      buf(:) = 42

      DO icomm=1,topo%ncommseq
          irank = topo%icommseq(icomm)
          tag = 100
          dest = MAX(ppm_rank,irank)
          src = MIN(ppm_rank,irank)

          CALL MPI_Barrier(ppm_comm,info)
          tstart = MPI_Wtime()
          DO i=1,lncomm
              IF (ppm_rank.EQ.src) THEN
                  CALL MPI_Send(buf,small,MPI_INTEGER,dest,i,&
                  &             ppm_comm,info)
              ELSE
                  CALL MPI_Recv(buf,small,MPI_INTEGER,src,i,&
                  &             ppm_comm,stat,info)
              ENDIF
!             CALL MPI_SendRecv(sbuf,small,MPI_INTEGER,irank,tag, &
!     &                         rbuf,small,MPI_INTEGER,irank,tag, &
!     &                         ppm_comm,stat,info)

          ENDDO

          tend = MPI_Wtime()

          latencies(icomm) = (tend-tstart)/(2*lncomm)
      ENDDO


      DEALLOCATE(buf,stat=info)
      IF (info.NE.0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_alloc,caller,     &
     &        'dealloc latency buffer',__LINE__,info)
          GOTO 9999
      ENDIF
      !-------------------------------------------------------------------------
      !  Measure Bandwidth
      !-------------------------------------------------------------------------
      ALLOCATE(buf(big),stat=info)
      IF (info.NE.0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_alloc,caller,     &
     &        'bandwidth buffer',__LINE__,info)
          GOTO 9999
      ENDIF

      buf(:) = 42

      DO icomm=1,topo%ncommseq
          irank = topo%icommseq(icomm)
          tag = 200

          CALL MPI_Barrier(ppm_comm,info)
          tstart = MPI_Wtime()

          DO i=1,bncomm
              IF (ppm_rank.EQ.src) THEN
                  CALL MPI_Send(buf,big,MPI_INTEGER,dest,i,&
                  &             ppm_comm,info)
              ELSE
                  CALL MPI_Recv(buf,big,MPI_INTEGER,src,i,&
                  &             ppm_comm,stat,info)
              ENDIF

!             CALL MPI_SendRecv(sbuf,big,MPI_INTEGER,irank,tag, &
!             &                 rbuf,big,MPI_INTEGER,irank,tag, &
!             &                 ppm_comm,stat,info)

          ENDDO

          tend = MPI_Wtime()
          bandwidths(icomm) = big/( ((tend-tstart)-latencies(icomm))/bncomm)
      ENDDO


      DEALLOCATE(buf,stat=info)
      IF (info.NE.0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_alloc,caller,     &
     &        'dealloc bandwidth buffer',__LINE__,info)
          GOTO 9999
      ENDIF
#else
      info = ppm_error_warning
      CALL ppm_error(ppm_err_nompi,caller,     &
     &    'Bandwidth and latency default to 0',__LINE__,info)

#endif
      !-------------------------------------------------------------------------
      !  Computer summary statistics
      !-------------------------------------------------------------------------
8888  CONTINUE
      bandwidth = 0.0_mk
      latency = 0.0_mk
      DO i=1,topo%ncommseq
          latency = latency + latencies(i)
          bandwidth = bandwidth + bandwidths(i)
          !write(*,'(3I,2F14.2)'),i,ppm_rank,topo%icommseq(i),&
          !&                      latencies(i),bandwidths(i)
      ENDDO
      bandwidth = 4*bandwidth / topo%ncommseq


      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
9999  CONTINUE
      CALL substop(caller,t0,info)

END SUBROUTINE ppm_netstat
