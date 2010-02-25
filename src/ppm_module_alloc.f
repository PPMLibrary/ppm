      !-------------------------------------------------------------------------
      !  Module       :                  ppm_module_alloc
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the utility
      !                 routines.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_alloc.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:56  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.2  2004/11/11 15:23:01  ivos
      !  Moved allocatable work data to the module.
      !
      !  Revision 1.1  2004/07/26 07:29:24  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
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

      MODULE ppm_module_alloc

         USE ppm_module_data, ONLY: ppm_kind_single,ppm_kind_double

         !----------------------------------------------------------------------
         !  Work arrays for reallocation.
         !----------------------------------------------------------------------
         REAL(ppm_kind_single)   , DIMENSION(:        ), POINTER :: work_1ds
         REAL(ppm_kind_double)   , DIMENSION(:        ), POINTER :: work_1dd
         COMPLEX(ppm_kind_single), DIMENSION(:        ), POINTER :: work_1dsc
         COMPLEX(ppm_kind_double), DIMENSION(:        ), POINTER :: work_1ddc
         INTEGER                 , DIMENSION(:        ), POINTER :: work_1di
         LOGICAL                 , DIMENSION(:        ), POINTER :: work_1dl

         REAL(ppm_kind_single)   , DIMENSION(:,:      ), POINTER :: work_2ds
         REAL(ppm_kind_double)   , DIMENSION(:,:      ), POINTER :: work_2dd
         COMPLEX(ppm_kind_single), DIMENSION(:,:      ), POINTER :: work_2dsc
         COMPLEX(ppm_kind_double), DIMENSION(:,:      ), POINTER :: work_2ddc
         INTEGER                 , DIMENSION(:,:      ), POINTER :: work_2di
         LOGICAL                 , DIMENSION(:,:      ), POINTER :: work_2dl

         REAL(ppm_kind_single)   , DIMENSION(:,:,:    ), POINTER :: work_3ds
         REAL(ppm_kind_double)   , DIMENSION(:,:,:    ), POINTER :: work_3dd
         COMPLEX(ppm_kind_single), DIMENSION(:,:,:    ), POINTER :: work_3dsc
         COMPLEX(ppm_kind_double), DIMENSION(:,:,:    ), POINTER :: work_3ddc
         INTEGER                 , DIMENSION(:,:,:    ), POINTER :: work_3di
         LOGICAL                 , DIMENSION(:,:,:    ), POINTER :: work_3dl

         REAL(ppm_kind_single)   , DIMENSION(:,:,:,:  ), POINTER :: work_4ds
         REAL(ppm_kind_double)   , DIMENSION(:,:,:,:  ), POINTER :: work_4dd
         COMPLEX(ppm_kind_single), DIMENSION(:,:,:,:  ), POINTER :: work_4dsc
         COMPLEX(ppm_kind_double), DIMENSION(:,:,:,:  ), POINTER :: work_4ddc
         INTEGER                 , DIMENSION(:,:,:,:  ), POINTER :: work_4di
         LOGICAL                 , DIMENSION(:,:,:,:  ), POINTER :: work_4dl

         REAL(ppm_kind_single)   , DIMENSION(:,:,:,:,:), POINTER :: work_5ds
         REAL(ppm_kind_double)   , DIMENSION(:,:,:,:,:), POINTER :: work_5dd
         COMPLEX(ppm_kind_single), DIMENSION(:,:,:,:,:), POINTER :: work_5dsc
         COMPLEX(ppm_kind_double), DIMENSION(:,:,:,:,:), POINTER :: work_5ddc
         INTEGER                 , DIMENSION(:,:,:,:,:), POINTER :: work_5di
         LOGICAL                 , DIMENSION(:,:,:,:,:), POINTER :: work_5dl

         PRIVATE :: work_1ds,work_1dd,work_1dsc,work_1ddc,work_1di,work_1dl
         PRIVATE :: work_2ds,work_2dd,work_2dsc,work_2ddc,work_2di,work_2dl
         PRIVATE :: work_3ds,work_3dd,work_3dsc,work_3ddc,work_3di,work_3dl
         PRIVATE :: work_4ds,work_4dd,work_4dsc,work_4ddc,work_4di,work_4dl
         PRIVATE :: work_5ds,work_5dd,work_5dsc,work_5ddc,work_5di,work_5dl

         !----------------------------------------------------------------------
         !  Define interfaces to allocation routines
         !----------------------------------------------------------------------
         INTERFACE ppm_alloc
            MODULE PROCEDURE ppm_alloc_1ds
            MODULE PROCEDURE ppm_alloc_1dd
            MODULE PROCEDURE ppm_alloc_1dsc
            MODULE PROCEDURE ppm_alloc_1ddc
            MODULE PROCEDURE ppm_alloc_1di
            MODULE PROCEDURE ppm_alloc_1dl
            MODULE PROCEDURE ppm_alloc_1dls
            MODULE PROCEDURE ppm_alloc_1dld
            MODULE PROCEDURE ppm_alloc_1dlsc
            MODULE PROCEDURE ppm_alloc_1dldc
            MODULE PROCEDURE ppm_alloc_1dli
            MODULE PROCEDURE ppm_alloc_1dll

            MODULE PROCEDURE ppm_alloc_2ds
            MODULE PROCEDURE ppm_alloc_2dd
            MODULE PROCEDURE ppm_alloc_2dsc
            MODULE PROCEDURE ppm_alloc_2ddc
            MODULE PROCEDURE ppm_alloc_2di
            MODULE PROCEDURE ppm_alloc_2dl
            MODULE PROCEDURE ppm_alloc_2dls
            MODULE PROCEDURE ppm_alloc_2dld
            MODULE PROCEDURE ppm_alloc_2dlsc
            MODULE PROCEDURE ppm_alloc_2dldc
            MODULE PROCEDURE ppm_alloc_2dli
            MODULE PROCEDURE ppm_alloc_2dll

            MODULE PROCEDURE ppm_alloc_3ds
            MODULE PROCEDURE ppm_alloc_3dd
            MODULE PROCEDURE ppm_alloc_3dsc
            MODULE PROCEDURE ppm_alloc_3ddc
            MODULE PROCEDURE ppm_alloc_3di
            MODULE PROCEDURE ppm_alloc_3dl
            MODULE PROCEDURE ppm_alloc_3dls
            MODULE PROCEDURE ppm_alloc_3dld
            MODULE PROCEDURE ppm_alloc_3dlsc
            MODULE PROCEDURE ppm_alloc_3dldc
            MODULE PROCEDURE ppm_alloc_3dli
            MODULE PROCEDURE ppm_alloc_3dll
            
            MODULE PROCEDURE ppm_alloc_4ds
            MODULE PROCEDURE ppm_alloc_4dd
            MODULE PROCEDURE ppm_alloc_4dsc
            MODULE PROCEDURE ppm_alloc_4ddc
            MODULE PROCEDURE ppm_alloc_4di
            MODULE PROCEDURE ppm_alloc_4dl
            MODULE PROCEDURE ppm_alloc_4dls
            MODULE PROCEDURE ppm_alloc_4dld
            MODULE PROCEDURE ppm_alloc_4dlsc
            MODULE PROCEDURE ppm_alloc_4dldc
            MODULE PROCEDURE ppm_alloc_4dli
            MODULE PROCEDURE ppm_alloc_4dll

            MODULE PROCEDURE ppm_alloc_5ds
            MODULE PROCEDURE ppm_alloc_5dd
            MODULE PROCEDURE ppm_alloc_5dsc
            MODULE PROCEDURE ppm_alloc_5ddc
            MODULE PROCEDURE ppm_alloc_5di
            MODULE PROCEDURE ppm_alloc_5dl
            MODULE PROCEDURE ppm_alloc_5dls
            MODULE PROCEDURE ppm_alloc_5dld
            MODULE PROCEDURE ppm_alloc_5dlsc
            MODULE PROCEDURE ppm_alloc_5dldc
            MODULE PROCEDURE ppm_alloc_5dli
            MODULE PROCEDURE ppm_alloc_5dll
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS
 
#define __KIND __SINGLE_PRECISION
#include "ppm_alloc_1d.f"
#include "ppm_alloc_2d.f"
#include "ppm_alloc_3d.f"
#include "ppm_alloc_4d.f"
#include "ppm_alloc_5d.f"
#include "ppm_alloc_1dl.f"
#include "ppm_alloc_2dl.f"
#include "ppm_alloc_3dl.f"
#include "ppm_alloc_4dl.f"
#include "ppm_alloc_5dl.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_alloc_1d.f"
#include "ppm_alloc_2d.f"
#include "ppm_alloc_3d.f"
#include "ppm_alloc_4d.f"
#include "ppm_alloc_5d.f"
#include "ppm_alloc_1dl.f"
#include "ppm_alloc_2dl.f"
#include "ppm_alloc_3dl.f"
#include "ppm_alloc_4dl.f"
#include "ppm_alloc_5dl.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_alloc_1d.f"
#include "ppm_alloc_2d.f"
#include "ppm_alloc_3d.f"
#include "ppm_alloc_4d.f"
#include "ppm_alloc_5d.f"
#include "ppm_alloc_1dl.f"
#include "ppm_alloc_2dl.f"
#include "ppm_alloc_3dl.f"
#include "ppm_alloc_4dl.f"
#include "ppm_alloc_5dl.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_alloc_1d.f"
#include "ppm_alloc_2d.f"
#include "ppm_alloc_3d.f"
#include "ppm_alloc_4d.f"
#include "ppm_alloc_5d.f"
#include "ppm_alloc_1dl.f"
#include "ppm_alloc_2dl.f"
#include "ppm_alloc_3dl.f"
#include "ppm_alloc_4dl.f"
#include "ppm_alloc_5dl.f"
#undef __KIND

#define __KIND __INTEGER
#include "ppm_alloc_1d.f"
#include "ppm_alloc_2d.f"
#include "ppm_alloc_3d.f"
#include "ppm_alloc_4d.f"
#include "ppm_alloc_5d.f"
#include "ppm_alloc_1dl.f"
#include "ppm_alloc_2dl.f"
#include "ppm_alloc_3dl.f"
#include "ppm_alloc_4dl.f"
#include "ppm_alloc_5dl.f"
#undef __KIND

#define __KIND __LOGICAL
#include "ppm_alloc_1d.f"
#include "ppm_alloc_2d.f"
#include "ppm_alloc_3d.f"
#include "ppm_alloc_4d.f"
#include "ppm_alloc_5d.f"
#include "ppm_alloc_1dl.f"
#include "ppm_alloc_2dl.f"
#include "ppm_alloc_3dl.f"
#include "ppm_alloc_4dl.f"
#include "ppm_alloc_5dl.f"
#undef __KIND

      END MODULE ppm_module_alloc
