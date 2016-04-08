      !-------------------------------------------------------------------------
      !  Subroutine   :                     ppm_util_printstats
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

      SUBROUTINE ppm_util_printstats(info,verbose)
      !!! Printout statistics
      !!! Optionally, print the results for all ranks
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_write
      USE ppm_module_mpi
      USE ppm_module_typedef
      IMPLICIT NONE

      INTEGER, PARAMETER :: MK = ppm_kind_double
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER                , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      LOGICAL                , INTENT(IN   ), OPTIONAL :: verbose
      !!! not implemented yet
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3)    :: ldu
      INTEGER                  :: i,isize
#ifdef __MPI
      INTEGER                  :: MPTYPE = MPI_DOUBLE_PRECISION
#endif
      REAL(MK)                 :: t0,t1,t2,t3
      !!! Current CPU clock time
      CHARACTER(LEN=ppm_char)  :: mesg
      CHARACTER(LEN=ppm_char)  :: caller = 'ppm_util_printstats'

      info = 0
      !-------------------------------------------------------------------------
      !  Call ppm_util_time
      !-------------------------------------------------------------------------

      isize = ppm_tstats_idx
#ifdef __MPI
      CALL MPI_Reduce(ppm_tstats_times(1:isize),&
          ppm_tstats_times_max(1:isize),isize,MPTYPE,MPI_MAX,0,ppm_comm,info)
      CALL MPI_Reduce(ppm_tstats_times(1:isize),&
          ppm_tstats_times_min(1:isize),isize,MPTYPE,MPI_MIN,0,ppm_comm,info)
      CALL MPI_Reduce(ppm_tstats_times(1:isize),&
          ppm_tstats_times_avg(1:isize),isize,MPTYPE,MPI_SUM,0,ppm_comm,info)

      IF (ppm_rank .EQ. 0) THEN
          stdout("----- TIMINGS - [max/avg/min] ----")
          DO i = 1,isize
              mesg = ppm_tstats_labels(i)
              t1 = ppm_tstats_times_max(i)
              t2 = ppm_tstats_times_avg(i) / REAL(ppm_nproc,8)
              t3 = ppm_tstats_times_min(i)

              stdout_f('(A,A,A)',"  ",'TRIM(ADJUSTL(mesg))',":")
              stdout_f('(A,3(E11.4,A))',"   ",t1," / ",t2," / ",t3)
          ENDDO
          stdout("--------------")
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      RETURN
      END SUBROUTINE ppm_util_printstats
