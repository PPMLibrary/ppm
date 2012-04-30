      !-------------------------------------------------------------------------
      !  Subroutine   :                     ppm_tstats_tic
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

      SUBROUTINE ppm_tstats_tic(id,step,info)
      !!! Puts the current cpu time in a buffer
      !!! Calls ppm_util_time,
      !!! which uses either `MPI_Wtime`, f90 `CPU_TIME` or `etime`,
      !!! based on how PPM was configured (see
      !!! ./configure --help):
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_write
      USE ppm_module_error
      USE ppm_module_util_time
      IMPLICIT NONE
      INTEGER, PARAMETER :: MK = ppm_kind_double
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                , INTENT(IN   ) :: id
      INTEGER                , INTENT(IN   ) :: step
      INTEGER                , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      CHARACTER(LEN=ppm_char)             :: caller = 'ppm_tstats_tic'
      INTEGER, DIMENSION(3)               :: ldu
      INTEGER                             :: iopt
      INTEGER                             :: istat
      
      info = 0
      !-------------------------------------------------------------------------
      !  Call ppm_util_time
      !-------------------------------------------------------------------------

      ! check that the array in ppm_btic%tic is large enough - reallocate if not
      IF (step.GT.ppm_tstats_nsamples) THEN
          ldu(1) = 2 * ppm_tstats_nsamples
          iopt = ppm_param_alloc_grow_preserve
          DO istat = 1,ppm_ntstats
            CALL ppm_alloc(ppm_tstats(istat)%times,ldu,iopt,info)
            IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,caller,'growing ppm_tstats',__LINE__,info)
              GOTO 9999
            ENDIF
            ppm_tstats(istat)% &
            &   times(ppm_tstats_nsamples+1:2*ppm_tstats_nsamples) = 0.0_mk
          ENDDO
          ppm_tstats_nsamples = 2 * ppm_tstats_nsamples
      ENDIF

      ! add a timing entry in the buffer
      CALL ppm_util_time(ppm_tstats(id)%times(step))

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      RETURN
      END SUBROUTINE ppm_tstats_tic
