      !--*- f90 -*--------------------------------------------------------------
      !  Module       :            ppm_module_interp_m2p_vecd
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

#define __DOUBLE_PRECISION 2
#define __2D               3
#define __3D               4
#define __VEC              5

      MODULE ppm_module_interp_m2p_vecd
      !!! Contains the mesh to particle interpolation routines. Currently we
      !!! support 2nd order B-spline and MP4 interpolation schemes.
        !-------------------------------------------------------------------------
        !  Modules
        !-------------------------------------------------------------------------
        USE ppm_module_data
        USE ppm_module_error
        USE ppm_module_alloc
        USE ppm_module_substart
        USE ppm_module_substop
        USE ppm_module_interfaces, ONLY : ppm_t_equi_mesh_,ppm_t_field_, &
        &   ppm_t_subpatch_
        IMPLICIT NONE

        PRIVATE

        !-----------------------------------------------------------------------
        !  Interface
        !-----------------------------------------------------------------------
        !-----------------------------------------------------------------------
        !
        !-----------------------------------------------------------------------
        PUBLIC :: m2p_interp_bsp2_dv_2d
        PUBLIC :: m2p_interp_bsp2_dv_3d

        PUBLIC :: m2p_interp_mp4_dv_2d
        PUBLIC :: m2p_interp_mp4_dv_3d

      CONTAINS

#define __KIND  __DOUBLE_PRECISION
#define __DIME  __2D
#define __MODE  __VEC
        ! 2D VEC DOUBLE
#include "interpolate/m2p_interp_bsp2.f"
#include "interpolate/m2p_interp_mp4.f"
#undef  __MODE
#undef  __DIME

#define __DIME  __3D
#define __MODE  __VEC
        ! 3D VEC DOUBLE
#include "interpolate/m2p_interp_bsp2.f"
#include "interpolate/m2p_interp_mp4.f"
#undef  __MODE
#undef  __DIME
#undef  __KIND

#undef __SINGLE_PRECISION
#undef __DOUBLE_PRECISION
#undef __2D
#undef __3D
#undef __VEC

      END MODULE ppm_module_interp_m2p_vecd

