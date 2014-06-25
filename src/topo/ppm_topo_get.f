      !--*- f90 -*--------------------------------------------------------------
      !  Subroutine   :                   ppm_topo_get
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

      SUBROUTINE ppm_topo_get(topoid,topo,info)
      !!! This routine returns a pointer to the topology
      !!! pointed at by topoid.
      !!!
      !!! [TIP]
      !!! The user should never need to manually modify the topology structure.
      !!! If this routine is needed, then it most probably means that either the
      !!! library is missing a feature or the client code is not well designed!
      !!!
      !!! [WARNING]
      !!! In the current implementation you can seriously break things if you
      !!! *change* the members of the returned topology as only a pointer is
      !!! returned, rather than a copy of the internal data.

      ! TODO safer implementation
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
      USE ppm_module_check_id
      USE ppm_module_topo_typedef
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER                  , INTENT(IN   ) :: topoid
      !!! Topology ID
      TYPE(ppm_t_topo), POINTER                :: topo
      !!! Returns the topology pointed at by `topoid`
      INTEGER                  , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)     :: t0
      LOGICAL                   :: valid
!       TYPE(ppm_t_topo), POINTER :: in_topo => NULL()
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_topo_get',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Copy the topology
      !-------------------------------------------------------------------------
      !in_topo => ppm_topo(topoid)%t

      topo => ppm_topo(topoid)%t

!      topo%ID = topoid
!      topo%isdefined   = in_topo%isdefined
!      topo%min_physs   = in_topo%min_physs
!      topo%max_physs   = in_topo%max_physs
!      topo%min_physd   = in_topo%min_physd
!      topo%max_physd   = in_topo%max_physd
!      topo%bcdef       = in_topo%bcdef
!      topo%nsubs       = in_topo%nsubs
!      topo%min_subs    = in_topo%min_subs
!      topo%max_subs    = in_topo%max_subs
!      topo%min_subd    = in_topo%min_subd
!      topo%max_subd    = in_topo%max_subd
!      topo%sub_costs   = in_topo%sub_costs
!      topo%sub_costd   = in_topo%sub_costd
!      topo%sub2proc    = in_topo%sub2proc
!      topo%nsublist    = in_topo%nsublist
!      topo%isublist    = in_topo%isublist
!      topo%subs_bc     = in_topo%subs_bc
!      topo%nneighsubs  = in_topo%nneighsubs
!      topo%ineighsubs  = in_topo%ineighsubs
!      topo%nneighproc  = in_topo%nneighproc
!      topo%ineighproc  = in_topo%ineighproc
!      topo%isoptimized = in_topo%isoptimized
!      topo%ncommseq    = in_topo%ncommseq
!      topo%icommseq    = in_topo%icommseq
!      topo%max_meshid  = in_topo%max_meshid
!      topo%mesh        = in_topo%mesh ! TODO: check if this works


      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_topo_get',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_topo_get',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 8888
          ENDIF
          CALL ppm_check_topoid(topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_get',  &
     &            'Topology ID is invalid!',__LINE__,info)
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_topo_get
