      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                  ppm_module_alloc
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

      !-------------------------------------------------------------------------
      !  Define types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION           1
#define __DOUBLE_PRECISION           2
#define __INTEGER                    3
#define __LOGICAL                    4
#define __SINGLE_PRECISION_COMPLEX   5
#define __DOUBLE_PRECISION_COMPLEX   6
#define __LONGINT                    7

      MODULE ppm_module_alloc
      !!! This module provides all PPM allocation routines. These routines allow
      !!! for allocation, reallocation, growing and shrinking of one- to five-
      !!! dimensional arrays in all intrinsic data types.
      !!!
      !!! When extending the PPM core or numerical routines these routines *must*
      !!! be used. PPM clients may take advantage of these functionalities.

         USE ppm_module_data, ONLY: ppm_kind_single,ppm_kind_double,&
  &                                 ppm_kind_int64

         !----------------------------------------------------------------------
         !  Work arrays for reallocation.
         !----------------------------------------------------------------------
         REAL(ppm_kind_single)   , DIMENSION(:        ), POINTER :: work_1ds
         REAL(ppm_kind_double)   , DIMENSION(:        ), POINTER :: work_1dd
         COMPLEX(ppm_kind_single), DIMENSION(:        ), POINTER :: work_1dsc
         COMPLEX(ppm_kind_double), DIMENSION(:        ), POINTER :: work_1ddc
         INTEGER                 , DIMENSION(:        ), POINTER :: work_1di
         INTEGER(ppm_kind_int64) , DIMENSION(:        ), POINTER :: work_1dli
         LOGICAL                 , DIMENSION(:        ), POINTER :: work_1dl

         REAL(ppm_kind_single)   , DIMENSION(:,:      ), POINTER :: work_2ds
         REAL(ppm_kind_double)   , DIMENSION(:,:      ), POINTER :: work_2dd
         COMPLEX(ppm_kind_single), DIMENSION(:,:      ), POINTER :: work_2dsc
         COMPLEX(ppm_kind_double), DIMENSION(:,:      ), POINTER :: work_2ddc
         INTEGER                 , DIMENSION(:,:      ), POINTER :: work_2di
         INTEGER(ppm_kind_int64) , DIMENSION(:,:      ), POINTER :: work_2dli
         LOGICAL                 , DIMENSION(:,:      ), POINTER :: work_2dl

         REAL(ppm_kind_single)   , DIMENSION(:,:,:    ), POINTER :: work_3ds
         REAL(ppm_kind_double)   , DIMENSION(:,:,:    ), POINTER :: work_3dd
         COMPLEX(ppm_kind_single), DIMENSION(:,:,:    ), POINTER :: work_3dsc
         COMPLEX(ppm_kind_double), DIMENSION(:,:,:    ), POINTER :: work_3ddc
         INTEGER                 , DIMENSION(:,:,:    ), POINTER :: work_3di
         INTEGER(ppm_kind_int64) , DIMENSION(:,:,:    ), POINTER :: work_3dli
         LOGICAL                 , DIMENSION(:,:,:    ), POINTER :: work_3dl

         REAL(ppm_kind_single)   , DIMENSION(:,:,:,:  ), POINTER :: work_4ds
         REAL(ppm_kind_double)   , DIMENSION(:,:,:,:  ), POINTER :: work_4dd
         COMPLEX(ppm_kind_single), DIMENSION(:,:,:,:  ), POINTER :: work_4dsc
         COMPLEX(ppm_kind_double), DIMENSION(:,:,:,:  ), POINTER :: work_4ddc
         INTEGER                 , DIMENSION(:,:,:,:  ), POINTER :: work_4di
         INTEGER(ppm_kind_int64) , DIMENSION(:,:,:,:  ), POINTER :: work_4dli
         LOGICAL                 , DIMENSION(:,:,:,:  ), POINTER :: work_4dl

         REAL(ppm_kind_single)   , DIMENSION(:,:,:,:,:), POINTER :: work_5ds
         REAL(ppm_kind_double)   , DIMENSION(:,:,:,:,:), POINTER :: work_5dd
         COMPLEX(ppm_kind_single), DIMENSION(:,:,:,:,:), POINTER :: work_5dsc
         COMPLEX(ppm_kind_double), DIMENSION(:,:,:,:,:), POINTER :: work_5ddc
         INTEGER                 , DIMENSION(:,:,:,:,:), POINTER :: work_5di
         INTEGER(ppm_kind_int64) , DIMENSION(:,:,:,:,:), POINTER :: work_5dli
         LOGICAL                 , DIMENSION(:,:,:,:,:), POINTER :: work_5dl

         PRIVATE :: work_1ds,work_1dd,work_1dsc,work_1ddc,work_1di,work_1dl
         PRIVATE :: work_2ds,work_2dd,work_2dsc,work_2ddc,work_2di,work_2dl
         PRIVATE :: work_3ds,work_3dd,work_3dsc,work_3ddc,work_3di,work_3dl
         PRIVATE :: work_4ds,work_4dd,work_4dsc,work_4ddc,work_4di,work_4dl
         PRIVATE :: work_5ds,work_5dd,work_5dsc,work_5ddc,work_5di,work_5dl
         PRIVATE :: work_1dli,work_2dli,work_3dli,work_4dli,work_5dli

         !----------------------------------------------------------------------
         !  Define interfaces to allocation routines
         !----------------------------------------------------------------------
         INTERFACE ppm_alloc
            MODULE PROCEDURE alloc_1d_s
            MODULE PROCEDURE alloc_1d_d
            MODULE PROCEDURE alloc_1d_sc
            MODULE PROCEDURE alloc_1d_dc
            MODULE PROCEDURE alloc_1d_i
            MODULE PROCEDURE alloc_1d_li
            MODULE PROCEDURE alloc_1d_l
            MODULE PROCEDURE alloc_1dl_s
            MODULE PROCEDURE alloc_1dl_d
            MODULE PROCEDURE alloc_1dl_sc
            MODULE PROCEDURE alloc_1dl_dc
            MODULE PROCEDURE alloc_1dl_i
            MODULE PROCEDURE alloc_1dl_li
            MODULE PROCEDURE alloc_1dl_l

            MODULE PROCEDURE alloc_2d_s
            MODULE PROCEDURE alloc_2d_d
            MODULE PROCEDURE alloc_2d_sc
            MODULE PROCEDURE alloc_2d_dc
            MODULE PROCEDURE alloc_2d_i
            MODULE PROCEDURE alloc_2d_li
            MODULE PROCEDURE alloc_2d_l
            MODULE PROCEDURE alloc_2dl_s
            MODULE PROCEDURE alloc_2dl_d
            MODULE PROCEDURE alloc_2dl_sc
            MODULE PROCEDURE alloc_2dl_dc
            MODULE PROCEDURE alloc_2dl_i
            MODULE PROCEDURE alloc_2dl_li
            MODULE PROCEDURE alloc_2dl_l

            MODULE PROCEDURE alloc_3d_s
            MODULE PROCEDURE alloc_3d_d
            MODULE PROCEDURE alloc_3d_sc
            MODULE PROCEDURE alloc_3d_dc
            MODULE PROCEDURE alloc_3d_i
            MODULE PROCEDURE alloc_3d_li
            MODULE PROCEDURE alloc_3d_l
            MODULE PROCEDURE alloc_3dl_s
            MODULE PROCEDURE alloc_3dl_d
            MODULE PROCEDURE alloc_3dl_sc
            MODULE PROCEDURE alloc_3dl_dc
            MODULE PROCEDURE alloc_3dl_i
            MODULE PROCEDURE alloc_3dl_li
            MODULE PROCEDURE alloc_3dl_l

            MODULE PROCEDURE alloc_4d_s
            MODULE PROCEDURE alloc_4d_d
            MODULE PROCEDURE alloc_4d_sc
            MODULE PROCEDURE alloc_4d_dc
            MODULE PROCEDURE alloc_4d_i
            MODULE PROCEDURE alloc_4d_li
            MODULE PROCEDURE alloc_4d_l
            MODULE PROCEDURE alloc_4dl_s
            MODULE PROCEDURE alloc_4dl_d
            MODULE PROCEDURE alloc_4dl_sc
            MODULE PROCEDURE alloc_4dl_dc
            MODULE PROCEDURE alloc_4dl_i
            MODULE PROCEDURE alloc_4dl_li
            MODULE PROCEDURE alloc_4dl_l

            MODULE PROCEDURE alloc_5d_s
            MODULE PROCEDURE alloc_5d_d
            MODULE PROCEDURE alloc_5d_sc
            MODULE PROCEDURE alloc_5d_dc
            MODULE PROCEDURE alloc_5d_i
            MODULE PROCEDURE alloc_5d_li
            MODULE PROCEDURE alloc_5d_l
            MODULE PROCEDURE alloc_5dl_s
            MODULE PROCEDURE alloc_5dl_d
            MODULE PROCEDURE alloc_5dl_sc
            MODULE PROCEDURE alloc_5dl_dc
            MODULE PROCEDURE alloc_5dl_i
            MODULE PROCEDURE alloc_5dl_li
            MODULE PROCEDURE alloc_5dl_l

            MODULE PROCEDURE ppm_alloc_topo
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "alloc/ppm_alloc_1d.f"
#include "alloc/ppm_alloc_2d.f"
#include "alloc/ppm_alloc_3d.f"
#include "alloc/ppm_alloc_4d.f"
#include "alloc/ppm_alloc_5d.f"
#include "alloc/ppm_alloc_1dl.f"
#include "alloc/ppm_alloc_2dl.f"
#include "alloc/ppm_alloc_3dl.f"
#include "alloc/ppm_alloc_4dl.f"
#include "alloc/ppm_alloc_5dl.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "alloc/ppm_alloc_1d.f"
#include "alloc/ppm_alloc_2d.f"
#include "alloc/ppm_alloc_3d.f"
#include "alloc/ppm_alloc_4d.f"
#include "alloc/ppm_alloc_5d.f"
#include "alloc/ppm_alloc_1dl.f"
#include "alloc/ppm_alloc_2dl.f"
#include "alloc/ppm_alloc_3dl.f"
#include "alloc/ppm_alloc_4dl.f"
#include "alloc/ppm_alloc_5dl.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
#include "alloc/ppm_alloc_1d.f"
#include "alloc/ppm_alloc_2d.f"
#include "alloc/ppm_alloc_3d.f"
#include "alloc/ppm_alloc_4d.f"
#include "alloc/ppm_alloc_5d.f"
#include "alloc/ppm_alloc_1dl.f"
#include "alloc/ppm_alloc_2dl.f"
#include "alloc/ppm_alloc_3dl.f"
#include "alloc/ppm_alloc_4dl.f"
#include "alloc/ppm_alloc_5dl.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "alloc/ppm_alloc_1d.f"
#include "alloc/ppm_alloc_2d.f"
#include "alloc/ppm_alloc_3d.f"
#include "alloc/ppm_alloc_4d.f"
#include "alloc/ppm_alloc_5d.f"
#include "alloc/ppm_alloc_1dl.f"
#include "alloc/ppm_alloc_2dl.f"
#include "alloc/ppm_alloc_3dl.f"
#include "alloc/ppm_alloc_4dl.f"
#include "alloc/ppm_alloc_5dl.f"
#undef __KIND

#define __KIND __INTEGER
#include "alloc/ppm_alloc_1d.f"
#include "alloc/ppm_alloc_2d.f"
#include "alloc/ppm_alloc_3d.f"
#include "alloc/ppm_alloc_4d.f"
#include "alloc/ppm_alloc_5d.f"
#include "alloc/ppm_alloc_1dl.f"
#include "alloc/ppm_alloc_2dl.f"
#include "alloc/ppm_alloc_3dl.f"
#include "alloc/ppm_alloc_4dl.f"
#include "alloc/ppm_alloc_5dl.f"
#undef __KIND

#define __KIND __LONGINT
#include "alloc/ppm_alloc_1d.f"
#include "alloc/ppm_alloc_2d.f"
#include "alloc/ppm_alloc_3d.f"
#include "alloc/ppm_alloc_4d.f"
#include "alloc/ppm_alloc_5d.f"
#include "alloc/ppm_alloc_1dl.f"
#include "alloc/ppm_alloc_2dl.f"
#include "alloc/ppm_alloc_3dl.f"
#include "alloc/ppm_alloc_4dl.f"
#include "alloc/ppm_alloc_5dl.f"
#undef __KIND

#define __KIND __LOGICAL
#include "alloc/ppm_alloc_1d.f"
#include "alloc/ppm_alloc_2d.f"
#include "alloc/ppm_alloc_3d.f"
#include "alloc/ppm_alloc_4d.f"
#include "alloc/ppm_alloc_5d.f"
#include "alloc/ppm_alloc_1dl.f"
#include "alloc/ppm_alloc_2dl.f"
#include "alloc/ppm_alloc_3dl.f"
#include "alloc/ppm_alloc_4dl.f"
#include "alloc/ppm_alloc_5dl.f"
#undef __KIND

#include "alloc/ppm_alloc_topo.f"

#include "alloc/ppm_alloc_argcheck.f"

      END MODULE ppm_module_alloc
