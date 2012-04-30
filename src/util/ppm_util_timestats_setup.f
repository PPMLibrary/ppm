      !-------------------------------------------------------------------------
      !  Subroutine   :                     ppm_timestats_setup
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
      SUBROUTINE ppm_tstats_setup(nstats,info,nsamples)
      !!! Calls ppm_util_time, pop the last tic out of the buffer
      !!! and returns the difference.
      !!! Optionally, print the results on stdout.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_write
      USE ppm_module_error
      IMPLICIT NONE
      INTEGER, PARAMETER :: mk = ppm_kind_double
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                , INTENT(IN   ) :: nstats
      !!! Returns status, 0 upon success
      INTEGER                , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      INTEGER, OPTIONAL      , INTENT(IN   ) :: nsamples
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3)               :: ldu
      REAL(mk)               :: t0,t1
      !!! Current CPU clock time
      CHARACTER(LEN=ppm_char)                :: cbuf
      CHARACTER(LEN=ppm_char)                :: caller = 'ppm_tstats_setup'
    
      ppm_ntstats = nstats
      IF (PRESENT(nsamples)) THEN
        ppm_tstats_nsamples = nsamples
      ENDIF

      ALLOCATE(ppm_tstats(ppm_ntstats),STAT=info)
      IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
        &             'Could not allocate ppm_tstats',__LINE__,info)
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      RETURN
      END SUBROUTINE ppm_tstats_setup
      
      SUBROUTINE ppm_tstats_add(tlabel,id,info)
      !!! Calls ppm_util_time, pop the last tic out of the buffer
      !!! and returns the difference.
      !!! Optionally, print the results on stdout.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_write
      USE ppm_module_error
      IMPLICIT NONE
      INTEGER, PARAMETER :: mk = ppm_kind_double
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      CHARACTER(LEN=*)       , INTENT(IN   ) :: tlabel
      INTEGER                , INTENT(  OUT) :: id
      INTEGER                , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3)               :: ldu
      REAL(mk)               :: t0,t1
      !!! Current CPU clock time
      CHARACTER(LEN=ppm_char)                :: cbuf
      CHARACTER(LEN=ppm_char)                :: caller = 'ppm_tstats_add'
      INTEGER,                          SAVE :: istat = 1
      INTEGER                                :: iopt
     
      IF (istat.GT.ppm_ntstats) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
        &             'ppm_tstats is not large enough',__LINE__,info)
      ENDIF
      ldu = ppm_tstats_nsamples
      iopt = ppm_param_alloc_fit 
      CALL ppm_alloc(ppm_tstats(istat)%times,ldu,iopt,info)
      IF (info .NE. 0) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_alloc,caller,'ppm_tstats(istat)%times',__LINE__,info)
        GOTO 9999
      ENDIF
      ppm_tstats(istat)%times = 0.0_mk
      ppm_tstats(istat)%label = tlabel
      
      id = istat
      istat = istat + 1
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      RETURN
      END SUBROUTINE ppm_tstats_add
      
      SUBROUTINE ppm_tstats_collect(filename,info)
      !!! Calls ppm_util_time, pop the last tic out of the buffer
      !!! and returns the difference.
      !!! Optionally, print the results on stdout.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_write
      USE ppm_module_error
      
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      
      INTEGER, PARAMETER :: mk = ppm_kind_double
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      CHARACTER(LEN=*)       , INTENT(IN   ) :: filename
      INTEGER                , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3)               :: ldu
      REAL(mk)               :: t0,t1
      !!! Current CPU clock time
      CHARACTER(LEN=ppm_char)                :: cbuf
      CHARACTER(LEN=ppm_char)                :: caller = 'ppm_tstats_collect'
      REAL(mk), DIMENSION(ppm_nproc)         :: outvec
      INTEGER                                :: nsamples,istat,istep
      INTEGER                                :: iopt
      INTEGER                                :: iunit = 42
      CHARACTER(128)                         :: sfmt
      CHARACTER(1)                           :: tab = char(9)
      nsamples = ppm_tstats_nsamples

      ldu = nsamples
      iopt = ppm_param_alloc_fit 
      DO istat=1,ppm_ntstats
        CALL ppm_alloc(ppm_tstats(istat)%tmin,ldu,iopt,info)
        IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,caller,'ppm_tstats(istat)%tmin',__LINE__,info)
          GOTO 9999
        ENDIF
        ppm_tstats(istat)%tmin = 0.0_mk
        CALL ppm_alloc(ppm_tstats(istat)%tmax,ldu,iopt,info)
        IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,caller,'ppm_tstats(istat)%tmax',__LINE__,info)
          GOTO 9999
        ENDIF
        ppm_tstats(istat)%tmax = 0.0_mk
        CALL ppm_alloc(ppm_tstats(istat)%tmean,ldu,iopt,info)
        IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,caller,'ppm_tstats(istat)%tmean',__LINE__,info)
          GOTO 9999
        ENDIF
        ppm_tstats(istat)%tmean = 0.0_mk

#ifdef __MPI
        call MPI_Reduce(ppm_tstats(istat)%times,ppm_tstats(istat)%tmin,&
        &         nsamples,ppm_mpi_kind,MPI_MIN,0,ppm_comm,info)
        call MPI_Reduce(ppm_tstats(istat)%times,ppm_tstats(istat)%tmax,&
        &         nsamples,ppm_mpi_kind,MPI_MAX,0,ppm_comm,info)
        call MPI_Reduce(ppm_tstats(istat)%times,ppm_tstats(istat)%tmean,&
        &         nsamples,ppm_mpi_kind,MPI_SUM,0,ppm_comm,info)
        ppm_tstats(istat)%tmean(1:nsamples) = ppm_tstats(istat)%tmean(1:nsamples)&
        &                                   / ppm_nproc

#else
        ppm_tstats(istat)%tmin = ppm_tstats(istat)%times
        ppm_tstats(istat)%tmax = ppm_tstats(istat)%times
        ppm_tstats(istat)%tmean = ppm_tstats(istat)%times
#endif

      ENDDO
      IF (ppm_rank.EQ.0) THEN

        OPEN(iunit,file=filename)
        ! write header
        DO istat = 1,ppm_ntstats
          IF (istat.LT.ppm_ntstats) THEN
            WRITE(iunit,'(6A)',advance='no') &
            &     '"',ppm_tstats(istat)%label(1:LEN_TRIM(ppm_tstats(istat)%label)),'"'&
            &     ,tab,tab,tab
          ELSE
            WRITE(iunit,'(3A)') &
            &     '"',ppm_tstats(istat)%label(1:LEN_TRIM(ppm_tstats(istat)%label)),'"'
          ENDIF
        ENDDO
        ! write timings
        DO istep = 1,nsamples
          DO istat = 1,ppm_ntstats
            WRITE(iunit,'(E14.6,A)',advance='no') ppm_tstats(istat)%tmin(istep),tab
            WRITE(iunit,'(E14.6,A)',advance='no') ppm_tstats(istat)%tmax(istep),tab
            IF (istat.LT.ppm_ntstats) THEN
              WRITE(iunit,'(E14.6,A)',advance='no') ppm_tstats(istat)%tmean(istep),tab
            ELSE
              WRITE(iunit,'(E14.6)',advance='no') ppm_tstats(istat)%tmean(istep)
            ENDIF
          ENDDO
          WRITE(iunit,'(A)') ''
        ENDDO
        CLOSE(iunit)
      ENDIF
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      RETURN
      END SUBROUTINE ppm_tstats_collect
