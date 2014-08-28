      !--*- f90 -*--------------------------------------------------------------
      !  Module   :                  ppm_module_kdtree
      !-------------------------------------------------------------------------
      ! Copyright (c) 2012 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      !                    Center for Fluid Dynamics (DTU)
      !
      !
      ! This file is part of the Parallel Particle Mesh Library (PPM).
      !
      !  PPM is free software: you can redistribute it and/or modify
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

      MODULE ppm_module_kdtree
      !----------------------------------------------------------------------
      ! K-D tree routines in Fortran 90 by:
      ! Matthew B. Kennel, Institute For Nonlinear Science,
      ! reference: http://arxiv.org/abs/physics/0408067
      !
      ! It has been adapted and amended for PPM library by Yaser Afshar.
      !
      ![NOTE]
      ! In this adaptation the maximum query vector size has been set to three
      ! in case of future development and higher dimensionality this routine
      ! needs to be modified
      !----------------------------------------------------------------------
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

      USE ppm_module_data
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE

      PRIVATE

      !----------------------------------------------------------------------
      ! Global variables
      !----------------------------------------------------------------------
      INTEGER, PARAMETER :: bucket_size = 12
      ! The maximum number of points to keep in a terminal node.

#define  DTYPE(a) a/**/_s
#define __KIND __SINGLE_PRECISION
#define  MK ppm_kind_single
#define  _MK _ppm_kind_single
#include "tree/ppm_kdtree_typedef.f"

#define  DTYPE(a) a/**/_d
#define __KIND __DOUBLE_PRECISION
#define  MK ppm_kind_double
#define  _MK _ppm_kind_double
#include "tree/ppm_kdtree_typedef.f"

      TYPE(tree_search_record_s), SAVE, TARGET :: sr_s
      TYPE(tree_search_record_d), SAVE, TARGET :: sr_d
      ! A GLOBAL VARIABLE for search

      PUBLIC :: kdtree2_result_s,     kdtree2_result_d
      PUBLIC :: pq
      PUBLIC :: interval_s,           interval_d
      PUBLIC :: tree_node_s,          tree_node_d
      PUBLIC :: kdtree2_s,            kdtree2_d
      PUBLIC :: tree_search_record_s, tree_search_record_d
      PUBLIC :: sr_s,                 sr_d

      CONTAINS

#define  DTYPE(a) a/**/_s
#define __KIND __SINGLE_PRECISION
#define  MK ppm_kind_single
#define  _MK _ppm_kind_single
#include "tree/ppm_kdtree_typeproc.f"

#define  DTYPE(a) a/**/_d
#define __KIND __DOUBLE_PRECISION
#define  MK ppm_kind_double
#define  _MK _ppm_kind_double
#include "tree/ppm_kdtree_typeproc.f"

      END MODULE ppm_module_kdtree

