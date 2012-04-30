      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                 ppm_module_rmsh
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

      !-------------------------------------------------------------------------
      !  Define types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2
#define __2D               3
#define __3D               4
#define __VEC              5
#define __SCA              6

      MODULE ppm_module_rmsh
#ifdef COMPILEME
      !!! This module contains all interfaces to the remeshing routines. For
      !!! convenience all interpolation modules are `USE`d by this module.
      !!! Therefore its not necessary to include `ppm_module_interp_*` in the
      !!! client application. 
         !----------------------------------------------------------------------
         !  PPM modules
         !----------------------------------------------------------------------
         USE ppm_module_interp_m2p
         USE ppm_module_interp_p2m
         USE ppm_module_topo_typedef
         USE ppm_module_mesh_typedef
         
        !-----------------------------------------------------------------------
        !  Define interface ppm_rmsh_comp_weights
        !-----------------------------------------------------------------------
        INTERFACE ppm_rmsh_comp_weights
           MODULE PROCEDURE ppm_rmsh_comp_weights_s
           MODULE PROCEDURE ppm_rmsh_comp_weights_d
        END INTERFACE

        !-----------------------------------------------------------------------
        !  Define interface ppm_rmsh_create_part
        !-----------------------------------------------------------------------
        INTERFACE ppm_rmsh_create_part
           ! 2d scalar
           MODULE PROCEDURE ppm_rmsh_create_part_sss_2d
           MODULE PROCEDURE ppm_rmsh_create_part_ssv_2d
           MODULE PROCEDURE ppm_rmsh_create_part_dss_2d
           MODULE PROCEDURE ppm_rmsh_create_part_dsv_2d
           ! 2d vector
           MODULE PROCEDURE ppm_rmsh_create_part_svs_2d
           MODULE PROCEDURE ppm_rmsh_create_part_svv_2d
           MODULE PROCEDURE ppm_rmsh_create_part_dvs_2d
           MODULE PROCEDURE ppm_rmsh_create_part_dvv_2d
           ! 3d scalar
           MODULE PROCEDURE ppm_rmsh_create_part_sss_3d
           MODULE PROCEDURE ppm_rmsh_create_part_ssv_3d
           MODULE PROCEDURE ppm_rmsh_create_part_dss_3d
           MODULE PROCEDURE ppm_rmsh_create_part_dsv_3d
           ! 3d vector
           MODULE PROCEDURE ppm_rmsh_create_part_svs_3d
           MODULE PROCEDURE ppm_rmsh_create_part_svv_3d
           MODULE PROCEDURE ppm_rmsh_create_part_dvs_3d
           MODULE PROCEDURE ppm_rmsh_create_part_dvv_3d
        END INTERFACE
       
        !-----------------------------------------------------------------------
        !  Define interface ppm_rmsh_remesh
        !-----------------------------------------------------------------------
        INTERFACE ppm_rmsh_remesh
           ! 2d scalar
           MODULE PROCEDURE ppm_rmsh_remesh_ss_2d
           MODULE PROCEDURE ppm_rmsh_remesh_ds_2d
           ! 2d vector
           MODULE PROCEDURE ppm_rmsh_remesh_sv_2d
           MODULE PROCEDURE ppm_rmsh_remesh_dv_2d
           ! 3d scalar
           MODULE PROCEDURE ppm_rmsh_remesh_ss_3d
           MODULE PROCEDURE ppm_rmsh_remesh_ds_3d
           ! 3d vector
           MODULE PROCEDURE ppm_rmsh_remesh_sv_3d
           MODULE PROCEDURE ppm_rmsh_remesh_dv_3d
        END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
        CONTAINS

#define __KIND  __SINGLE_PRECISION
#include "rmsh/ppm_rmsh_comp_weights.f"
#undef  __KIND        
#define __KIND __DOUBLE_PRECISION 
#include "rmsh/ppm_rmsh_comp_weights.f"
#undef  __KIND        

#define __KIND  __SINGLE_PRECISION
#define __DIME  __2D
#define __MODE  __SCA
        ! 2D SCA SINGLE
#define __MODE2  __SCA
#include "rmsh/ppm_rmsh_create_part.f"
#undef __MODE2
#define __MODE2  __VEC
#include "rmsh/ppm_rmsh_create_part.f"
#undef __MODE2

#undef  __MODE
#define __MODE  __VEC
        ! 2D VEC SINGLE
#define __MODE2  __SCA
#include "rmsh/ppm_rmsh_create_part.f"
#undef __MODE2
#define __MODE2  __VEC
#include "rmsh/ppm_rmsh_create_part.f"
#undef __MODE2
#undef  __MODE
#undef  __DIME
        
#define __DIME  __3D
#define __MODE  __SCA
        ! 3D SCA SINGLE
#define __MODE2 __SCA
#include "rmsh/ppm_rmsh_create_part.f"
#undef __MODE2
#define __MODE2 __VEC
#include "rmsh/ppm_rmsh_create_part.f"
#undef __MODE2
#undef  __MODE
#define __MODE  __VEC
        ! 3D VEC SINGLE
#define __MODE2 __SCA
#include "rmsh/ppm_rmsh_create_part.f"
#undef __MODE2
#define __MODE2 __VEC
#include "rmsh/ppm_rmsh_create_part.f"
#undef __MODE2
#undef  __MODE
#undef  __DIME
#undef  __KIND


#define __KIND  __DOUBLE_PRECISION
#define __DIME  __2D
#define __MODE  __SCA
        ! 2D SCA DOUBLE
#define __MODE2 __SCA
#include "rmsh/ppm_rmsh_create_part.f"
#undef __MODE2
#define __MODE2 __VEC
#include "rmsh/ppm_rmsh_create_part.f"
#undef __MODE2
#undef  __MODE
#define __MODE  __VEC
        ! 2D VEC DOUBLE
#define __MODE2 __SCA
#include "rmsh/ppm_rmsh_create_part.f"
#undef __MODE2
#define __MODE2 __VEC
#include "rmsh/ppm_rmsh_create_part.f"
#undef __MODE2
#undef  __MODE
#undef  __DIME
        
#define __DIME  __3D
#define __MODE  __SCA
        ! 3D SCA DOUBLE
#define __MODE2 __SCA
#include "rmsh/ppm_rmsh_create_part.f"
#undef __MODE2
#define __MODE2 __VEC
#include "rmsh/ppm_rmsh_create_part.f"
#undef __MODE2
#undef  __MODE
#define __MODE  __VEC
        ! 3D VEC DOUBLE
#define __MODE2 __SCA
#include "rmsh/ppm_rmsh_create_part.f"
#undef __MODE2
#define __MODE2 __VEC
#include "rmsh/ppm_rmsh_create_part.f"
#undef __MODE2
#undef  __MODE
#undef  __DIME
#undef  __KIND        

#define __KIND  __SINGLE_PRECISION
#define __DIME  __2D
#define __MODE  __SCA
        ! 2D SCA SINGLE
#include "rmsh/ppm_rmsh_remesh.f"
#undef  __MODE
#define __MODE  __VEC
        ! 2D VEC SINGLE
#include "rmsh/ppm_rmsh_remesh.f"
#undef  __MODE
#undef  __DIME
        
#define __DIME  __3D
#define __MODE  __SCA
        ! 3D SCA SINGLE
#include "rmsh/ppm_rmsh_remesh.f"
#undef  __MODE
#define __MODE  __VEC
        ! 3D VEC SINGLE
#include "rmsh/ppm_rmsh_remesh.f"
#undef  __MODE
#undef  __DIME
#undef  __KIND


#define __KIND  __DOUBLE_PRECISION
#define __DIME  __2D
#define __MODE  __SCA
        ! 2D SCA DOUBLE
#include "rmsh/ppm_rmsh_remesh.f"
#undef  __MODE
#define __MODE  __VEC
        ! 2D VEC DOUBLE
#include "rmsh/ppm_rmsh_remesh.f"
#undef  __MODE
#undef  __DIME
        
#define __DIME  __3D
#define __MODE  __SCA
        ! 3D SCA DOUBLE
#include "rmsh/ppm_rmsh_remesh.f"
#undef  __MODE
#define __MODE  __VEC
        ! 3D VEC DOUBLE
#include "rmsh/ppm_rmsh_remesh.f"
#undef  __MODE
#undef  __DIME
#undef  __KIND        


#endif
      END MODULE ppm_module_rmsh

