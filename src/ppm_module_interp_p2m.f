      !--*- f90 -*--------------------------------------------------------------
      !  Module       :            ppm_module_rmsh_remesh
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

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2
#define __2D               3
#define __3D               4
#define __VEC              5
#define __SCA              6

      MODULE ppm_module_interp_p2m
      !!! Contains the particle to mesh interpolation routines. Currently we
      !!! support 2nd order B-spline and MP4 interpolation schemes.
      
        !-----------------------------------------------------------------------
        !  Interface
        !-----------------------------------------------------------------------
        INTERFACE ppm_interp_p2m
           ! 2d scalar
           MODULE PROCEDURE p2m_ss_2d
           MODULE PROCEDURE p2m_ds_2d
           ! 2d vector
           MODULE PROCEDURE p2m_sv_2d
           MODULE PROCEDURE p2m_dv_2d
           ! 3d scalar
           MODULE PROCEDURE p2m_ss_3d
           MODULE PROCEDURE p2m_ds_3d
           ! 3d vector
           MODULE PROCEDURE p2m_sv_3d
           MODULE PROCEDURE p2m_dv_3d
        END INTERFACE

        INTERFACE ppm_interp_p2m_reno
           ! 2d scalar
           MODULE PROCEDURE p2m_renorm_ss_3d
           MODULE PROCEDURE p2m_renorm_ds_3d
           ! 3d vector
           MODULE PROCEDURE p2m_renorm_sv_3d
           MODULE PROCEDURE p2m_renorm_dv_3d
        END INTERFACE

        INTERFACE p2m_interp_bsp2
            MODULE PROCEDURE p2m_interp_bsp2_ss_2d
            MODULE PROCEDURE p2m_interp_bsp2_ds_2d
            MODULE PROCEDURE p2m_interp_bsp2_sv_2d
            MODULE PROCEDURE p2m_interp_bsp2_dv_2d
            MODULE PROCEDURE p2m_interp_bsp2_ss_3d
            MODULE PROCEDURE p2m_interp_bsp2_ds_3d
            MODULE PROCEDURE p2m_interp_bsp2_sv_3d
            MODULE PROCEDURE p2m_interp_bsp2_dv_3d
        END INTERFACE

        INTERFACE p2m_interp_mp4
            MODULE PROCEDURE p2m_interp_mp4_ss_2d
            MODULE PROCEDURE p2m_interp_mp4_ds_2d
            MODULE PROCEDURE p2m_interp_mp4_sv_2d
            MODULE PROCEDURE p2m_interp_mp4_dv_2d
            MODULE PROCEDURE p2m_interp_mp4_ss_3d
            MODULE PROCEDURE p2m_interp_mp4_ds_3d
            MODULE PROCEDURE p2m_interp_mp4_sv_3d
            MODULE PROCEDURE p2m_interp_mp4_dv_3d
        END INTERFACE
        
      CONTAINS

        
#define __KIND  __SINGLE_PRECISION
#define __DIME  __2D
#define __MODE  __SCA
        ! 2D SCA SINGLE
#include "interpolate/ppm_interp_p2m.f"
#include "interpolate/p2m_interp_bsp2.f"
#include "interpolate/p2m_interp_mp4.f"
#undef  __MODE
#define __MODE  __VEC
        ! 2D VEC SINGLE
#include "interpolate/ppm_interp_p2m.f"
#include "interpolate/p2m_interp_bsp2.f"
#include "interpolate/p2m_interp_mp4.f"
#undef  __MODE
#undef  __DIME
        
#define __DIME  __3D
#define __MODE  __SCA
        ! 3D SCA SINGLE
#include "interpolate/ppm_interp_p2m.f"
#include "interpolate/ppm_interp_p2m_renorm.f"
#include "interpolate/p2m_interp_bsp2.f"
#include "interpolate/p2m_interp_mp4.f"
#undef  __MODE
#define __MODE  __VEC
        ! 3D VEC SINGLE
#include "interpolate/ppm_interp_p2m.f"
#include "interpolate/ppm_interp_p2m_renorm.f"
#include "interpolate/p2m_interp_bsp2.f"
#include "interpolate/p2m_interp_mp4.f"
#undef  __MODE
#undef  __DIME
#undef  __KIND


#define __KIND  __DOUBLE_PRECISION
#define __DIME  __2D
#define __MODE  __SCA
        ! 2D SCA DOUBLE
#include "interpolate/ppm_interp_p2m.f"
#include "interpolate/p2m_interp_bsp2.f"
#include "interpolate/p2m_interp_mp4.f"
#undef  __MODE
#define __MODE  __VEC
        ! 2D VEC DOUBLE
#include "interpolate/ppm_interp_p2m.f"
#include "interpolate/p2m_interp_bsp2.f"
#include "interpolate/p2m_interp_mp4.f"
#undef  __MODE
#undef  __DIME
        
#define __DIME  __3D
#define __MODE  __SCA
        ! 3D SCA DOUBLE
#include "interpolate/ppm_interp_p2m.f"
#include "interpolate/ppm_interp_p2m_renorm.f"
#include "interpolate/p2m_interp_bsp2.f"
#include "interpolate/p2m_interp_mp4.f"
#undef  __MODE
#define __MODE  __VEC
        ! 3D VEC DOUBLE
#include "interpolate/ppm_interp_p2m.f"
#include "interpolate/ppm_interp_p2m_renorm.f"
#include "interpolate/p2m_interp_bsp2.f"
#include "interpolate/p2m_interp_mp4.f"
#undef  __MODE
#undef  __DIME
#undef  __KIND        





#undef __SINGLE_PRECISION 
#undef __DOUBLE_PRECISION 
#undef __2D               
#undef __3D               
#undef __VEC              
#undef __SCA              
        
      END MODULE ppm_module_interp_p2m

