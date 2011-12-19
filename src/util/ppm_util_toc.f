      !-------------------------------------------------------------------------
      !  Subroutine   :                     ppm_util_toc
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

      SUBROUTINE ppm_util_toc(mesg,diff_t,info,verbose)
      !!! Calls ppm_util_time, pop the last tic out of the buffer
      !!! and returns the difference.
      !!! Optionally, print the results on stdout.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_write
      IMPLICIT NONE
      INTEGER, PARAMETER :: MK = ppm_kind_double
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      CHARACTER(LEN=*)       , INTENT(IN   ) :: mesg
      !!! Returns status, 0 upon success
      REAL(MK)               , INTENT(  OUT) :: diff_t
      !!! Difference in time between tic and toc
      INTEGER                , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      LOGICAL                , INTENT(IN   ), OPTIONAL :: verbose
      !!! Difference in time between tic and toc
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3)               :: ldu
      REAL(MK)               :: t0,t1
      !!! Current CPU clock time
      CHARACTER(LEN=ppm_char)                :: cbuf
      CHARACTER(LEN=ppm_char)                :: caller = 'ppm_util_toc'
      
      info = 0
      !-------------------------------------------------------------------------
      !  Call ppm_util_time
      !-------------------------------------------------------------------------
      IF (ppm_btic%idx.LE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_sub_failed,'ppm_util_toc',&
              'wrong number of tics and tocs',__LINE__,info)
          GOTO 9999
      ELSE
          t0 = ppm_btic%tic(ppm_btic%idx)
          ppm_btic%idx = ppm_btic%idx -1
      ENDIF
      CALL ppm_util_time(t1)

      ! difference between current time and last entry in the tic buffer
      diff_t = t1-t0

      ! increment counter for statistics
      ppm_tstats_idx = ppm_tstats_idx + 1

      ! check that the array in ppm_btic%tic is large enough - reallocate if not
      IF (ppm_tstats_idx.GT.max_size) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_sub_failed,caller,&
          'array ppm_t_tstats%stat overflow - max_size too small',__LINE__,info)
          GOTO 9999
      ENDIF

      ppm_tstats_times(ppm_tstats_idx)  = diff_t
      ppm_tstats_labels(ppm_tstats_idx) = mesg

      IF (PRESENT(verbose)) THEN
          IF (verbose) THEN
              WRITE(cbuf,'(A,A,E17.7,A)') mesg, ' took ',diff_t,' seconds'
              CALL ppm_write(ppm_rank,caller,cbuf,info)
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      RETURN
      END SUBROUTINE ppm_util_toc
