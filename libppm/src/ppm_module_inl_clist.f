      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                 ppm_module_inl_clist
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

      MODULE ppm_module_inl_clist
      !!! This module is used to create cell lists stored in various depths, providing
      !!! ease of USE to form inhomogeneous cell lists

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

        !-------------------------------------------------------------------------
        !  Used modules
        !-------------------------------------------------------------------------
        USE ppm_module_data
        USE ppm_module_inl_hash
        USE ppm_module_alloc
        USE ppm_module_error
        USE ppm_module_substart
        USE ppm_module_substop
        USE ppm_module_mktopo

        !-------------------------------------------------------------------------
        !  Declaration of parameters
        !-------------------------------------------------------------------------

        !-------------------------------------------------------------------------
        !  Declaration of arrays
        !-------------------------------------------------------------------------
        INTEGER, DIMENSION(:,:), POINTER :: borders  => NULL()
        !!! contains the boundaries in the particle rank array separating the
        !!! particles belonging to different cells.
        !!!
        !!! in 2D this is a 6xn_cells array
        !!! The column indeces of borders(:,k) contain indeces to the rank array 
        !!! specifying which particles belong in which subcell of cell k
        !!! -----------
        !!! |2:3 |4:5 |
        !!! |----|----|
        !!! |1:2 |3:4 |
        !!! -----------
        !!! e.g.
        !!! index 1 (+1) to index 2 are all particles belonging to lower left subcell
        !!!   --> rank(borders(1,k)+1:borders(2,k))
        !!!
        !!! index 6 is 1 if the cell contains particles in deeper levels,
        !!! otherwise it is -1
        INTEGER, DIMENSION(:),   POINTER :: rank       => NULL()
        ! rank of particles
        INTEGER, DIMENSION(:),   POINTER :: rankByPos  => NULL()
        ! rank of particles
        INTEGER, DIMENSION(:),   POINTER :: rc_borders => NULL()

        !-------------------------------------------------------------------------
        !  Declaration of variables
        !-------------------------------------------------------------------------
        INTEGER :: borders_pos
        INTEGER :: borders_pos_max
        INTEGER :: max_depth
        INTEGER :: ncell
        INTEGER :: n_real_p
        INTEGER :: n_all_p
        INTEGER :: insuf_hash_table

        INTERFACE ppm_create_inl_clist
            MODULE PROCEDURE create_inl_clist_s
            MODULE PROCEDURE create_inl_clist_d
        END INTERFACE

        INTERFACE getCellCoor_Depth
            MODULE PROCEDURE getCellCoor_Depth_s
            MODULE PROCEDURE getCellCoor_Depth_d
        END INTERFACE

        INTERFACE getCellIdx
            MODULE PROCEDURE getCellIdx_s
            MODULE PROCEDURE getCellIdx_d
        END INTERFACE

        INTERFACE getMinimumRC
            MODULE PROCEDURE getMinimumRC_s
            MODULE PROCEDURE getMinimumRC_d
        END INTERFACE

        INTERFACE setSubregions
            MODULE PROCEDURE setSubregions_s
            MODULE PROCEDURE setSubregions_d
        END INTERFACE

        INTERFACE SortByPosition
            MODULE PROCEDURE SortByPosition_s
            MODULE PROCEDURE SortByPosition_d
        END INTERFACE

        INTERFACE SortByRC_Pos
            MODULE PROCEDURE SortByRC_Pos_s
            MODULE PROCEDURE SortByRC_Pos_d
        END INTERFACE

        INTERFACE partition
            MODULE PROCEDURE partition_s
            MODULE PROCEDURE partition_d
        END INTERFACE

        INTERFACE sortByRC
            MODULE PROCEDURE sortByRC_s
            MODULE PROCEDURE sortByRC_d
        END INTERFACE

        INTERFACE partitionByRC
            MODULE PROCEDURE partitionByRC_s
            MODULE PROCEDURE partitionByRC_d
        END INTERFACE

        INTERFACE lastIdxForRC
            MODULE PROCEDURE lastIdxForRC_s
            MODULE PROCEDURE lastIdxForRC_d
        END INTERFACE

        INTERFACE getRC_Borders
            MODULE PROCEDURE getRC_Borders_s
            MODULE PROCEDURE getRC_Borders_d
        END INTERFACE

        INTERFACE getMinimumSideLength
            MODULE PROCEDURE getMinimumSideLength_s
            MODULE PROCEDURE getMinimumSideLength_d
        END INTERFACE

        INTERFACE getMaxDepth
            MODULE PROCEDURE getMaxDepth_s
            MODULE PROCEDURE getMaxDepth_d
        END INTERFACE

        !-------------------------------------------------------------------------
        !  Privatizing arrays, variables and parameters of the module
        !-------------------------------------------------------------------------

        CONTAINS

#define __KIND __SINGLE_PRECISION
#include "neighlist/ppm_inl_clist.f"
#undef  __KIND
#define __KIND __DOUBLE_PRECISION
#include "neighlist/ppm_inl_clist.f"
#undef  __KIND

#undef __SINGLE_PRECISION
#undef __DOUBLE_PRECISION

      END MODULE ppm_module_inl_clist
