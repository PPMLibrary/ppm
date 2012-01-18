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
     !!!
     !!!WARNING: this module needs DC operators. 
     !!! PPM needs to be compiled with --enable-dcops

#ifdef __DCOPS


#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

#define __MORSE 1
#define __SCHRADER 2
#define __REPULSIVE 3

#ifndef __SOP_POTENTIAL
#define __SOP_POTENTIAL __MORSE
#endif

#define debug_verbosity 2
#define __USE_RANDOMNUMBERS 1
     !method for minimisation of the interaction potential
     ! __USE_LBFGS  for L-BFGS
     ! __USE_SD     for steepest descent using a line search
#define __USE_SD 1
#undef __USE_LBFGS 
!#undef __USE_SD 
!#define __USE_LBFGS 1
#define __USE_DEL_METHOD2 1

     USE ppm_module_data
     USE ppm_module_typedef
     USE ppm_module_particles_typedef
     USE ppm_module_particles
     USE ppm_module_sop_typedef
     USE ppm_module_write
     USE ppm_module_substart
     USE ppm_module_substop

     IMPLICIT NONE

     INTERFACE sop_approx_wp
         MODULE PROCEDURE sop_approx_wp_1d_s
         MODULE PROCEDURE sop_approx_wp_1d_d
         MODULE PROCEDURE sop_approx_wp_2d_s
         MODULE PROCEDURE sop_approx_wp_2d_d
     END INTERFACE
     INTERFACE sop_dump_debug
         MODULE PROCEDURE sop_dump_1d_s
         MODULE PROCEDURE sop_dump_1d_d
         MODULE PROCEDURE sop_dump_1di
         MODULE PROCEDURE sop_dump_2d_s
         MODULE PROCEDURE sop_dump_2d_d
         MODULE PROCEDURE sop_dump_2di
     END INTERFACE

     INTERFACE sop_adapt_particles
         MODULE PROCEDURE sop_adapt_particles_s
         MODULE PROCEDURE sop_adapt_particles_d
     END INTERFACE
     INTERFACE sop_close_neighbours
         MODULE PROCEDURE sop_close_neighbours_s
         MODULE PROCEDURE sop_close_neighbours_d
     END INTERFACE
     INTERFACE sop_compute_D
         MODULE PROCEDURE sop_compute_D_s
         MODULE PROCEDURE sop_compute_D_d
     END INTERFACE
     INTERFACE sop_fuse_particles
         MODULE PROCEDURE sop_fuse_particles_s
         MODULE PROCEDURE sop_fuse_particles_d
     END INTERFACE
     INTERFACE sop_fuse2_particles
         MODULE PROCEDURE sop_fuse2_particles_s
         MODULE PROCEDURE sop_fuse2_particles_d
     END INTERFACE
     INTERFACE sop_gradient_descent
         MODULE PROCEDURE sop_gradient_descent_s
         MODULE PROCEDURE sop_gradient_descent_d
     END INTERFACE
     INTERFACE sop_gradient_descent_ls
         MODULE PROCEDURE sop_gradient_descent_ls_s
         MODULE PROCEDURE sop_gradient_descent_ls_d
     END INTERFACE
     INTERFACE sop_gradient_psi
         MODULE PROCEDURE sop_gradient_psi_s
         MODULE PROCEDURE sop_gradient_psi_d
     END INTERFACE
     INTERFACE sop_interpolate
         MODULE PROCEDURE sop_interpolate_s
         MODULE PROCEDURE sop_interpolate_d
     END INTERFACE
     INTERFACE sop_potential_psi
         MODULE PROCEDURE sop_potential_psi_s
         MODULE PROCEDURE sop_potential_psi_d
     END INTERFACE
     INTERFACE sop_spawn_particles
         MODULE PROCEDURE sop_spawn_particles_s
         MODULE PROCEDURE sop_spawn_particles_d
     END INTERFACE
     INTERFACE sop_spawn2_particles
         MODULE PROCEDURE sop_spawn2_particles_s
         MODULE PROCEDURE sop_spawn2_particles_d
     END INTERFACE
     INTERFACE sop_init_opts
         MODULE PROCEDURE sop_init_opts_s
         MODULE PROCEDURE sop_init_opts_d
     END INTERFACE
     INTERFACE sop_init_stats
         MODULE PROCEDURE sop_init_stats_s
         MODULE PROCEDURE sop_init_stats_d
     END INTERFACE
     INTERFACE sop_plot_potential
         MODULE PROCEDURE sop_plot_potential_s
         MODULE PROCEDURE sop_plot_potential_d
     END INTERFACE
     INTERFACE check_duplicates
         MODULE PROCEDURE check_duplicates_s
         MODULE PROCEDURE check_duplicates_d
     END INTERFACE
     INTERFACE sop_voronoi_MC
         MODULE PROCEDURE sop_voronoi_MC_s
         MODULE PROCEDURE sop_voronoi_MC_d
     END INTERFACE

    !----------------------------------------------------------------------
    ! Private variables for the module
    !----------------------------------------------------------------------

     PRIVATE

     INTEGER , DIMENSION(3)  :: ldc
     !!! Number of elements in all dimensions for allocation

     !====================================================================!
     ! Variable describing the 'state' of the system
     !====================================================================!
!     REAL(ppm_kind_single)         :: Psi_global_d
!     REAL(ppm_kind_double)         :: Psi_global_s
     ! total interaction potential of the particles
!     REAL(ppm_kind_single)         :: Psi_global_old_s,Psi_max_s
!     REAL(ppm_kind_double)         :: Psi_global_old_d,Psi_max_d

     INTEGER            :: adapt_wpgradid
     ! id of where the gradient of the field is stored 
     INTEGER            :: fuse_id = 0
     ! id of the property storing which particles are candidate for fusion

     INTEGER            :: voronoi_id = 0
     ! id of the property storing the voronoi volume of each particle

     INTEGER            :: nb_neigh_id = 0

     INTEGER            :: potential_before_id = 0
     ! for debugging only

     INTEGER            :: memory_used_total
     ! Evaluation of the number of bytes occupied in memory

     LOGICAL            :: adaptation_ok

     !====================================================================!
     ! maths
     !====================================================================!
     !REAL(prec),PARAMETER :: PI = ACOS(-1._prec)


     PUBLIC sop_adapt_particles, sop_init_opts, sop_init_stats,&
         &  sop_plot_potential, sop_compute_D, sop_fuse2_particles


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

#define  __KIND __SINGLE_PRECISION
#define  DTYPE(a) a/**/_s
#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_single
#include "sop/ppm_sop_helpers.f"
#include "sop/sop_adapt_particles.f"
#include "sop/sop_close_neighbours.f"
#include "sop/sop_compute_D.f"
#include "sop/sop_dump_debug.f"
#include "sop/sop_fuse_particles.f"
#include "sop/sop_fuse2_particles.f"
#include "sop/sop_gradient_descent.f"
#include "sop/sop_gradient_psi.f"
#include "sop/sop_interpolate.f"
#include "sop/sop_potential_psi.f"
#include "sop/sop_spawn_particles.f"
#include "sop/sop_spawn2_particles.f"
#include "sop/sop_voronoi_MC.f"
#define __LDA 1
#include "sop/sop_approx_wp.f"
#define __LDA 2
#include "sop/sop_approx_wp.f"
#undef   __KIND
#undef   DTYPE
#undef   DEFINE_MK

#define  __KIND __DOUBLE_PRECISION
#define  DTYPE(a) a/**/_d
#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_double
#include "sop/ppm_sop_helpers.f"
#include "sop/sop_adapt_particles.f"
#include "sop/sop_close_neighbours.f"
#include "sop/sop_compute_D.f"
#include "sop/sop_dump_debug.f"
#include "sop/sop_fuse_particles.f"
#include "sop/sop_fuse2_particles.f"
#include "sop/sop_gradient_descent.f"
#include "sop/sop_gradient_psi.f"
#include "sop/sop_interpolate.f"
#include "sop/sop_potential_psi.f"
#include "sop/sop_spawn_particles.f"
#include "sop/sop_spawn2_particles.f"
#include "sop/sop_voronoi_MC.f"
#define __LDA 1
#include "sop/sop_approx_wp.f"
#define __LDA 2
#include "sop/sop_approx_wp.f"
#undef __KIND
#undef   DTYPE
#undef   DEFINE_MK

#endif
     END MODULE ppm_module_sop
