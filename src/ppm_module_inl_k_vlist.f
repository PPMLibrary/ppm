      !--*- f90 -*--------------------------------------------------------------
      !  Module   :                  ppm_module_inl_k_vlist
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
      !-------------------------------------------------------------------------
      !  MOSAIC Group
      !  Max Planck Institute of Molecular Cell Biology and Genetics
      !  Pfotenhauerstr. 108, 01307 Dresden, Germany
      !-------------------------------------------------------------------------
      MODULE ppm_module_inl_k_vlist
      !!! This module contains routines and functions used in creating verlet lists
      !!! by use of kdtree.

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2
        !-------------------------------------------------------------------------
        !  Used modules
        !-------------------------------------------------------------------------
        USE ppm_module_data
        USE ppm_module_error
        USE ppm_module_alloc
        USE ppm_module_write
        USE ppm_module_substart
        USE ppm_module_substop
        USE ppm_module_kdtree
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Privatizing arrays, variables and parameters
        !-------------------------------------------------------------------------
        PRIVATE

        !-------------------------------------------------------------------------
        !  Declaration of arrays
        !-------------------------------------------------------------------------

        !-------------------------------------------------------------------------
        !  Declaration of variables
        !-------------------------------------------------------------------------

        !-------------------------------------------------------------------------
        !  Declaration of interfaces
        !-------------------------------------------------------------------------
        INTERFACE process_terminal_node
            MODULE PROCEDURE process_terminal_node_s
            MODULE PROCEDURE process_terminal_node_d
        END INTERFACE

        INTERFACE process_terminal_node_fixedball
            MODULE PROCEDURE process_terminal_node_fixedball_s
            MODULE PROCEDURE process_terminal_node_fixedball_d
        END INTERFACE

        INTERFACE search
            MODULE PROCEDURE search_s
            MODULE PROCEDURE search_d
        END INTERFACE

        INTERFACE kdtree_n_nearest
            MODULE PROCEDURE kdtree_n_nearest_s
            MODULE PROCEDURE kdtree_n_nearest_d
        END INTERFACE

        INTERFACE kdtree_n_nearest_around_point
            MODULE PROCEDURE kdtree_n_nearest_around_point_s
            MODULE PROCEDURE kdtree_n_nearest_around_point_d
        END INTERFACE

        INTERFACE kdtree_r_nearest
            MODULE PROCEDURE kdtree_r_nearest_s
            MODULE PROCEDURE kdtree_r_nearest_d
        END INTERFACE

        INTERFACE kdtree_r_nearest_around_point
            MODULE PROCEDURE kdtree_r_nearest_around_point_s
            MODULE PROCEDURE kdtree_r_nearest_around_point_d
        END INTERFACE

        INTERFACE kdtree_r_count
            MODULE PROCEDURE kdtree_r_count_s
            MODULE PROCEDURE kdtree_r_count_d
        END INTERFACE

        INTERFACE kdtree_r_count_around_point
            MODULE PROCEDURE kdtree_r_count_around_point_s
            MODULE PROCEDURE kdtree_r_count_around_point_d
        END INTERFACE

        INTERFACE kdtree_n_nearest_brute_force
            MODULE PROCEDURE kdtree_n_nearest_brute_force_s
            MODULE PROCEDURE kdtree_n_nearest_brute_force_d
        END INTERFACE

        INTERFACE kdtree_r_nearest_brute_force
            MODULE PROCEDURE kdtree_r_nearest_brute_force_s
            MODULE PROCEDURE kdtree_r_nearest_brute_force_d
        END INTERFACE

        INTERFACE qsort
            MODULE PROCEDURE kdtree_util_qsort_s
            MODULE PROCEDURE kdtree_util_qsort_d
        END INTERFACE

        PUBLIC :: kdtree_n_nearest
        PUBLIC :: kdtree_n_nearest_around_point
        PUBLIC :: kdtree_r_nearest
        PUBLIC :: kdtree_r_nearest_around_point
        PUBLIC :: kdtree_r_count
        PUBLIC :: kdtree_r_count_around_point
        PUBLIC :: kdtree_n_nearest_brute_force
        PUBLIC :: kdtree_r_nearest_brute_force

        CONTAINS

#define  DTYPE(a) a/**/_s
#define __KIND __SINGLE_PRECISION
#define  MK ppm_kind_single
#define  _MK _ppm_kind_single
#include "neighlist/ppm_inl_k_vlist.f"

#define  DTYPE(a) a/**/_d
#define __KIND __DOUBLE_PRECISION
#define  MK ppm_kind_double
#define  _MK _ppm_kind_double
#include "neighlist/ppm_inl_k_vlist.f"

#undef __SINGLE_PRECISION
#undef __DOUBLE_PRECISION
      END MODULE ppm_module_inl_k_vlist
