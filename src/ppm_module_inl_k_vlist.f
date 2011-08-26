     !--*- f90 -*-------------------------------------------------------------- 
     !  Module   :                  ppm_module_inl_k_vlist
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
      MODULE ppm_module_inl_k_vlist
      !!! This module contains routines and functions used in creating verlet lists
      !!! by use of inhomogeneous cell lists.

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2
#define __COUNT 3
#define __GET 4
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
        REAL(ppm_kind_single),DIMENSION(:),POINTER,PRIVATE :: dist_array_s => NULL()
        REAL(ppm_kind_double),DIMENSION(:),POINTER,PRIVATE :: dist_array_d => NULL()
        INTEGER,       DIMENSION(:),   POINTER :: dist_rank => NULL()

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
        INTERFACE ppm_inl_k_vlist
            MODULE PROCEDURE inl_k_vlist_s
            MODULE PROCEDURE inl_k_vlist_d
        END INTERFACE

        INTERFACE create_inl_vlist
            MODULE PROCEDURE create_inl_vlist_s
            MODULE PROCEDURE create_inl_vlist_d
        END INTERFACE

        INTERFACE getVerletLists
            MODULE PROCEDURE getVerletLists_s
            MODULE PROCEDURE getVerletLists_d
        END INTERFACE

        INTERFACE inDomain
            MODULE PROCEDURE inDomain_s
            MODULE PROCEDURE inDomain_d
        END INTERFACE

        INTERFACE is_kNeighbor
            MODULE PROCEDURE is_kNeighbor_s
            MODULE PROCEDURE is_kNeighbor_d
        END INTERFACE
        
        INTERFACE cross_neighbor
            MODULE PROCEDURE cross_neighbor_s
            MODULE PROCEDURE cross_neighbor_d
        END INTERFACE

        INTERFACE count_neigh
            MODULE PROCEDURE count_neigh_s
            MODULE PROCEDURE count_neigh_d
        END INTERFACE

        !INTERFACE count_neigh_sym 
            !MODULE PROCEDURE count_neigh_sym_s
            !MODULE PROCEDURE count_neigh_sym_d
        !END INTERFACE

        INTERFACE get_neigh
            MODULE PROCEDURE get_neigh_s
            MODULE PROCEDURE get_neigh_d
        END INTERFACE

        !INTERFACE get_neigh_sym
            !MODULE PROCEDURE get_neigh_sym_s
            !MODULE PROCEDURE get_neigh_sym_d
        !END INTERFACE

        INTERFACE getParticleCoorDepth
            MODULE PROCEDURE getParticleCoorDepth_s
            MODULE PROCEDURE getParticleCoorDepth_d
        END INTERFACE

        INTERFACE getParticlesInCell
            MODULE PROCEDURE getParticlesInCell_s
            MODULE PROCEDURE getParticlesInCell_d
        END INTERFACE

        INTERFACE getParticlesInCellDomain
            MODULE PROCEDURE getParticlesInCellDomain_s
            MODULE PROCEDURE getParticlesInCellDomain_d
        END INTERFACE

        INTERFACE getSubdomainParticles
            MODULE PROCEDURE getSubdomainParticles_s
            MODULE PROCEDURE getSubdomainParticles_d
        END INTERFACE

        INTERFACE SortByDist
            MODULE PROCEDURE SortByDist_s
            MODULE PROCEDURE SortByDist_d
        END INTERFACE

        INTERFACE partitionByDist
            MODULE PROCEDURE partitionByDist_s
            MODULE PROCEDURE partitionByDist_d
        END INTERFACE


        PRIVATE :: create_inl_vlist
        PRIVATE :: getVerletLists
        PRIVATE :: count_neigh!, count_neigh_sym
        PRIVATE :: get_neigh!, get_neigh_sym
        PRIVATE :: getSubdomainParticles
        PRIVATE :: getParticlesInCell
        PRIVATE :: getParticlesInCellDomain
        PRIVATE :: getParticleCoorDepth
        PRIVATE :: inDomain
        PRIVATE :: is_kNeighbor
        PRIVATE :: SortByDist
        PRIVATE :: partitionByDist
        PRIVATE :: cross_neighbor
        !-------------------------------------------------------------------------
        !  Privatizing arrays, variables and parameters
        !-------------------------------------------------------------------------
!!! to be completed when test driver is removed!

        CONTAINS

#define __KIND __SINGLE_PRECISION
#include "neighlist/ppm_inl_k_vlist.f"
#include "neighlist/ppm_inl_helpers.f"
#include "neighlist/ppm_inl_k_helpers.f"
#define __ACTION __COUNT
#include "neighlist/ppm_inl_k_vlist_build.f"
#undef __ACTION
#define __ACTION __GET
#include "neighlist/ppm_inl_k_vlist_build.f"
#undef __ACTION
#undef  __KIND
#define __KIND __DOUBLE_PRECISION
#include "neighlist/ppm_inl_k_vlist.f"
#include "neighlist/ppm_inl_helpers.f"
#include "neighlist/ppm_inl_k_helpers.f"
#define __ACTION __COUNT
#include "neighlist/ppm_inl_k_vlist_build.f"
#undef __ACTION
#define __ACTION __GET
#include "neighlist/ppm_inl_k_vlist_build.f"
#undef __ACTION
#undef  __KIND

#undef __SINGLE_PRECISION
#undef __DOUBLE_PRECISION
#undef __COUNT
#undef __GET
      END MODULE ppm_module_inl_k_vlist
