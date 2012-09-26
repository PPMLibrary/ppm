      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_dynamic_loadbal
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

      SUBROUTINE dynamic_loadbal(decomp,info)

      !!! This subroutine balances the workload amongst the processes by sending
      !!! or receiving subdomains locally.
      !!! Depending on the type of the domain decomposition used the DLB
      !!! algorithm takes respective action:
      !!! - Cuboid domain decomposition:
      !!! * An overloaded processor sends some of its subdomains to its neighbors
      !!!   with lower workload.
      !!! - Tree-like (e.g. tree, pencil and slab) decompositions:
      !!! * The subdomain that will be sent for load balancing gets refined by 
      !!! going one level up in the octree representation of the domain decomp
      !!! and a new repartitioning is done on this subdomain such that we have
      !!! more subdomains than before and thus we will be able to tune the balance
      !!! finer by sending smaller subdomains to underloaded neighbors.
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_loadbal
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: decomp
      !!! Decomposition type
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_dynamic_loadbal',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Depending on the decomposition type we will choose one of the DLB 
      !  strategies
      !-------------------------------------------------------------------------
      IF (                                            &
          !---------------------------------------------------------------------
          !  If we deal with a Cuboid/cartesian/user_defined decomp,
          !  use only existing subdomains to balance the computational load
          !---------------------------------------------------------------------
          decomp .EQ. ppm_param_decomp_cuboid    .OR. &
          decomp .EQ. ppm_param_decomp_cartesian .OR. &
          decomp .EQ. ppm_param_decomp_user_defined) THEN



 
      ENDIF
      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_dynamic_loadbal',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
 8888   CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE dynamic_loadbal_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE dynamic_loadbal_d
#endif

