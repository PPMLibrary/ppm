     !-------------------------------------------------------------------------
     !  Module   :                  ppm_module_sop
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
     MODULE ppm_module_sop
     !!! This module contains routines and functions used for Self-Organizing
     !!! Particles

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

#define debug_verbosity 2

     USE ppm_module_data
     USE ppm_module_typedef
     USE ppm_module_particles
     USE ppm_module_sop_typedef
     USE ppm_module_write
     USE ppm_module_substart
     USE ppm_module_substop
#include "ppm_define.h"

     IMPLICIT NONE

     INTERFACE sop_approx_wp
         MODULE PROCEDURE sop_approx_wp_1d
         MODULE PROCEDURE sop_approx_wp_2d
     END INTERFACE
     INTERFACE sop_dump_debug
         MODULE PROCEDURE sop_dump_1d
         MODULE PROCEDURE sop_dump_1di
         MODULE PROCEDURE sop_dump_2d
         MODULE PROCEDURE sop_dump_2di
     END INTERFACE
     !INTERFACE sop_check_debug
     !MODULE PROCEDURE sop_check_1d
     !MODULE PROCEDURE sop_check_1di
     !MODULE PROCEDURE sop_check_2d
     !END INTERFACE


     !TODO: duplicate everything such that it works also
     ! for single precision

     PRIVATE

     INTEGER, PARAMETER :: prec = ppm_kind_double

     !private:: define_sop_args

     !====================================================================!
     ! Variable describing the 'state' of the system
     !====================================================================!
     LOGICAL            :: Psi_saturates ! true when the forces have
     !          to be truncated to prevent particles from flying across the domain
     !                (happens only when 2 particles are very close to each other)

     REAL(prec)         :: Psi_global ! mean interaction 
     !                   potential of the particles
     REAL(prec)         :: Psi_global_old 
     REAL(prec)         :: Psi_max  
     !====================================================================!
     ! random numbers
     !====================================================================!
     INTEGER                               :: seedsize
     INTEGER,  DIMENSION(:), ALLOCATABLE   :: seed
     REAL(prec), DIMENSION(:), ALLOCATABLE :: randnb
     INTEGER                               :: randnb_i
     !====================================================================!
     ! maths
     !====================================================================!
     REAL(prec),PARAMETER :: PI = ACOS(-1._prec)


     PUBLIC sop_adapt_particles, sop_init_opts

     CONTAINS

     !SUBROUTINE define_sop_args

     !CALL arg(fuse_radius, 'fuse_radius', &
     !ctrl_name = 'fuse_radius',       &
     !long_flag      =  '--sop-fuse-r',                 &
     !default   =  0._prec )

     !CALL arg(attractive_radius0, 'attractive_radius0', &
     !ctrl_name = 'attractive_radius',       &
     !long_flag      =  '--sop-attract-rad',                 &
     !default   =  0._prec )

     !CALL arg(adaptivity_criterion, 'adaptivity_criterion', &
     !ctrl_name = 'adaptivity_criterion',       &
     !long_flag      =  '--sop-adapt-crit',                 &
     !default   =  8._prec )

     !CALL arg(nneigh_theo, 'nneigh_theo', &
     !ctrl_name = 'nneigh_theo',       &
     !long_flag      =  '--nneigh-theo',                 &
     !default   =  24 )

     !CALL arg(rcp_over_D, 'rcp_over_D', &
     !ctrl_name = 'rcp_over_D',     &
     !long_flag      =  '--rcp-over-D',                 &
     !default   =  2._prec ,         &
     !min       =  0._prec           )

     !CALL arg(maximum_D, 'maximum_D', &
     !ctrl_name = 'maximum_D',     &
     !long_flag      =  '--sop-max-D',                 &
     !default   =  1._prec ,         &
     !min       =  0._prec           )

     !CALL arg(scale_D, 'scale_D', &
     !ctrl_name = 'scale_D',     &
     !long_flag      =  '--scale-D',                 &
     !default   =  1._prec ,         &
     !min       =  0._prec           )

     !CALL arg(minimum_D, 'minimum_D', &
     !ctrl_name = 'minimum_D',     &
     !long_flag      =  '--sop-min-D',                 &
     !default   =  0._prec ,         &
     !min       =  0._prec            )

     !CALL arg(write_pdb, 'write_pdb', &
     !ctrl_name = 'write_pdb',     &
     !long_flag      =  '--sop-write-pdb',                 &
     !default   =  .FALSE.)

     !CALL arg(write_xyz, 'write_xyz', &
     !ctrl_name = 'write_xyz',     &
     !long_flag =  '--sop-write-xyz',                 &
     !default   =  .FALSE.)

     !END SUBROUTINE define_sop_args


#define __KIND __DOUBLE_PRECISION
#include "sop/ppm_sop_helpers.f"

#define __KIND __DOUBLE_PRECISION
#include "sop/sop_adapt_particles.f"

#define __KIND __DOUBLE_PRECISION
#include "sop/sop_compute_D.f"

#define __KIND __DOUBLE_PRECISION
#include "sop/sop_dump_debug.f"

     !!#include "sop/sop_check_debug.f"

#define __KIND __DOUBLE_PRECISION
#include "sop/sop_fuse_particles.f"

#define __KIND __DOUBLE_PRECISION
#include "sop/sop_gradient_descent.f"

#define __KIND __DOUBLE_PRECISION
#include "sop/sop_gradient_psi.f"

!#define __KIND __DOUBLE_PRECISION
!#include "sop/sop_interpolate_particles.f"

#define __KIND __DOUBLE_PRECISION
#include "sop/sop_interpolate.f"

#define __KIND __DOUBLE_PRECISION
#include "sop/sop_potential_psi.f"


#define __KIND __DOUBLE_PRECISION
#include "sop/sop_spawn_particles.f"

#define __KIND __DOUBLE_PRECISION
#define __LDA 1
#include "sop/sop_approx_wp.f"
#define __KIND __DOUBLE_PRECISION
#define __LDA 2
#include "sop/sop_approx_wp.f"


#undef __KIND
#undef __DIM

     END MODULE ppm_module_sop
