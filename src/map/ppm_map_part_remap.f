      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_map_part_remap
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
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_remap_s(topoid,xp,Npart,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_remap_d(topoid,xp,Npart,info)
#endif
      !!! This routine maps the particles onto the given topology.
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_topo_check
      IMPLICIT NONE

#if    __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:)  , INTENT(IN   ) :: xp
      !!! The position of the particles
      INTEGER                   , INTENT(IN   ) :: Npart
      !!! The number of particles (on processor)
      INTEGER                   , INTENT(IN   ) :: topoid
      !!! Topology identifier of destination topology.
      INTEGER                   , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      LOGICAL :: topo_ok

      CHARACTER(LEN=ppm_char) :: caller="ppm_map_part_remap"
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         CALL check
         IF (info.NE.0) GOTO 9999
      ENDIF


      ! if there is still some data left in the buffer, warn the user
      IF (ppm_buffer_set.GT.0) THEN
         fail("Buffer was not empty. Possible loss of data!",ppm_err_map_incomp, &
         & ppm_error=ppm_error_warning,exit_point=no)
      ENDIF


      !-------------------------------------------------------------------------
      !  Check if current topology is ok
      !-------------------------------------------------------------------------
      CALL ppm_topo_check(topoid,xp,Npart,topo_ok,info)
      or_fail("ppm_topo_check")

      !-------------------------------------------------------------------------
      !  if topology is ok remap particles otherwise abort remapping
      !-------------------------------------------------------------------------
      IF (topo_ok) THEN
         CALL ppm_map_part_global(topoid,xp,Npart,info)
      ELSE
         fail("Particles are not on current topology. Mapping discarded.",ppm_err_topo_missm)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF (Npart.LE.0) THEN
             fail("Npart must be >0",exit_point=8888)
          ENDIF

          CALL ppm_check_topoid(topoid,topo_ok,info)
          IF (.NOT.topo_ok) THEN
             fail("topoid out of range",exit_point=8888)
          ENDIF
      8888 CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_remap_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_part_remap_d
#endif
