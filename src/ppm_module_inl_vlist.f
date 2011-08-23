     !--*- f90 -*-------------------------------------------------------------- 
     !  Module   :                  ppm_module_inl_vlist
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
      MODULE ppm_module_inl_vlist
      !!! This module contains routines and functions used in creating verlet lists
      !!! by use of inhomogeneous cell lists.

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2
#define __COUNT 3
#define __GET 4
#define __YES 1
#define __NO 2
        !-------------------------------------------------------------------------
        !  Used modules
        !-------------------------------------------------------------------------
        USE ppm_module_data
        USE ppm_module_inl_clist
        USE ppm_module_alloc
        USE ppm_module_error
        USE ppm_module_substart
        USE ppm_module_substop

        !-------------------------------------------------------------------------
        !  Declaration of arrays
        !-------------------------------------------------------------------------
        INTEGER,       DIMENSION(:),   POINTER :: own_plist => NULL()
        INTEGER,       DIMENSION(:),   POINTER :: neigh_plist => NULL()
        INTEGER(ppm_kind_int64), DIMENSION(:), POINTER :: empty_list => NULL()
        INTEGER,       DIMENSION(:,:), POINTER :: ncells => NULL()
        LOGICAL,       DIMENSION(:) ,  POINTER :: used => NULL()

        !-------------------------------------------------------------------------
        !  Declaration of variables
        !-------------------------------------------------------------------------
        TYPE(ppm_clist), SAVE :: clist
        INTEGER               :: own_nplist
        INTEGER               :: neigh_nplist
        INTEGER               :: empty_pos
        INTEGER               :: max_nneigh

        !-------------------------------------------------------------------------
        !  Declaration of interfaces
        !-------------------------------------------------------------------------
        INTERFACE ppm_inl_vlist
            MODULE PROCEDURE inl_vlist_s_aniso
            MODULE PROCEDURE inl_vlist_d_aniso
            MODULE PROCEDURE inl_vlist_s
            MODULE PROCEDURE inl_vlist_d
        END INTERFACE

        INTERFACE create_inl_vlist
            MODULE PROCEDURE create_inl_vlist_s_aniso
            MODULE PROCEDURE create_inl_vlist_d_aniso
            MODULE PROCEDURE create_inl_vlist_s
            MODULE PROCEDURE create_inl_vlist_d
        END INTERFACE

        INTERFACE getVerletLists
            MODULE PROCEDURE getVerletLists_s_aniso
            MODULE PROCEDURE getVerletLists_d_aniso
            MODULE PROCEDURE getVerletLists_s
            MODULE PROCEDURE getVerletLists_d
        END INTERFACE

        INTERFACE inDomain
            MODULE PROCEDURE inDomain_s
            MODULE PROCEDURE inDomain_d
        END INTERFACE

        INTERFACE isNeighbor
            MODULE PROCEDURE isNeighbor_s_aniso
            MODULE PROCEDURE isNeighbor_d_aniso
            MODULE PROCEDURE isNeighbor_s
            MODULE PROCEDURE isNeighbor_d
        END INTERFACE
        
        INTERFACE cross_neighbor
            MODULE PROCEDURE cross_neighbor_s
            MODULE PROCEDURE cross_neighbor_d
        END INTERFACE

        INTERFACE count_neigh
            MODULE PROCEDURE count_neigh_s_aniso
            MODULE PROCEDURE count_neigh_d_aniso
            MODULE PROCEDURE count_neigh_s
            MODULE PROCEDURE count_neigh_d
        END INTERFACE

        INTERFACE count_neigh_sym 
            MODULE PROCEDURE count_neigh_sym_s_aniso
            MODULE PROCEDURE count_neigh_sym_d_aniso
            MODULE PROCEDURE count_neigh_sym_s
            MODULE PROCEDURE count_neigh_sym_d
        END INTERFACE

        INTERFACE get_neigh
            MODULE PROCEDURE get_neigh_s_aniso
            MODULE PROCEDURE get_neigh_d_aniso
            MODULE PROCEDURE get_neigh_s
            MODULE PROCEDURE get_neigh_d
        END INTERFACE

        INTERFACE get_neigh_sym
            MODULE PROCEDURE get_neigh_sym_s_aniso
            MODULE PROCEDURE get_neigh_sym_d_aniso
            MODULE PROCEDURE get_neigh_sym_s
            MODULE PROCEDURE get_neigh_sym_d
        END INTERFACE

        INTERFACE getParticleCoorDepth
            MODULE PROCEDURE getParticleCoorDepth_s_aniso
            MODULE PROCEDURE getParticleCoorDepth_d_aniso
            MODULE PROCEDURE getParticleCoorDepth_s
            MODULE PROCEDURE getParticleCoorDepth_d
        END INTERFACE

        INTERFACE getParticlesInCell
            MODULE PROCEDURE getParticlesInCell_s
            MODULE PROCEDURE getParticlesInCell_d
        END INTERFACE

        INTERFACE getSubdomainParticles
            MODULE PROCEDURE getSubdomainParticles_s_aniso
            MODULE PROCEDURE getSubdomainParticles_d_aniso
            MODULE PROCEDURE getSubdomainParticles_s
            MODULE PROCEDURE getSubdomainParticles_d
        END INTERFACE

        PRIVATE :: create_inl_vlist
        PRIVATE :: getVerletLists
        PRIVATE :: count_neigh, count_neigh_sym
        PRIVATE :: get_neigh, get_neigh_sym
        PRIVATE :: getSubdomainParticles
        PRIVATE :: getParticlesInCell
        PRIVATE :: getParticleCoorDepth
        PRIVATE :: inDomain
        PRIVATE :: isNeighbor
        !-------------------------------------------------------------------------
        !  Privatizing arrays, variables and parameters
        !-------------------------------------------------------------------------
!!! to be completed when test driver is removed!

        CONTAINS
#define __ANISO __YES
#define __KIND __SINGLE_PRECISION
#include "neighlist/ppm_inl_vlist.f"
#include "neighlist/ppm_inl_helpers.f"
#define __ACTION __COUNT
#include "neighlist/ppm_inl_vlist_build.f"
#undef __ACTION
#define __ACTION __GET
#include "neighlist/ppm_inl_vlist_build.f"
#undef __ACTION
#undef  __KIND
#define __KIND __DOUBLE_PRECISION
#include "neighlist/ppm_inl_vlist.f"
#include "neighlist/ppm_inl_helpers.f"
#define __ACTION __COUNT
#include "neighlist/ppm_inl_vlist_build.f"
#undef __ACTION
#define __ACTION __GET
#include "neighlist/ppm_inl_vlist_build.f"
#undef __ACTION
#undef  __KIND
#undef  __ANISO

#define __ANISO __NO
#define __KIND __SINGLE_PRECISION
#include "neighlist/ppm_inl_vlist.f"
#include "neighlist/ppm_inl_helpers.f"
#define __ACTION __COUNT
#include "neighlist/ppm_inl_vlist_build.f"
#undef __ACTION
#define __ACTION __GET
#include "neighlist/ppm_inl_vlist_build.f"
#undef __ACTION
#undef  __KIND
#define __KIND __DOUBLE_PRECISION
#include "neighlist/ppm_inl_vlist.f"
#include "neighlist/ppm_inl_helpers.f"
#define __ACTION __COUNT
#include "neighlist/ppm_inl_vlist_build.f"
#undef __ACTION
#define __ACTION __GET
#include "neighlist/ppm_inl_vlist_build.f"
#undef __ACTION
#undef  __KIND

#undef __SINGLE_PRECISION
#undef __DOUBLE_PRECISION
#undef __COUNT
#undef __GET
#undef  __ANISO
#undef  __YES
#undef  __NO

      END MODULE ppm_module_inl_vlist
