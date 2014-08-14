#define __SINGLE_PRECISION         1
#define __DOUBLE_PRECISION         2
#define __INTEGER                  3
#define __LOGICAL                  4
#define __SINGLE_PRECISION_COMPLEX 5
#define __DOUBLE_PRECISION_COMPLEX 6

      MODULE ppm_module_mapping_typedef
      !!! This module defines the mapping data type and provides
      !!!  the basic mapping routines for particles and meshes; namely
      !!!  the push, send and pop routine.


         !----------------------------------------------------------------------
         !  Modules
         !----------------------------------------------------------------------
      !   USE ppm_module_map_part_util
      !   USE ppm_module_map_part_ghost
      !   USE ppm_module_map_part_global
      !   USE ppm_module_map_part_partial
         USE ppm_module_interfaces
         USE ppm_module_data
         USE ppm_module_substart
         USE ppm_module_substop
         IMPLICIT NONE

         PRIVATE

         PUBLIC :: ppm_t_ptr_part_mapping_s,ppm_t_ptr_part_mapping_d
         PUBLIC :: ppm_t_part_mapping_s,ppm_t_part_mapping_d


         !----------------------------------------------------------------------
         !  Type declaration
         !----------------------------------------------------------------------

#define  DTYPE(a) a/**/_s
#define  MK ppm_kind_single
#define  _MK _ppm_kind_single
#include "map/mapping_typedef.inc"

#define  DTYPE(a) a/**/_d
#define  MK ppm_kind_double
#define  _MK _ppm_kind_double
#include "map/mapping_typedef.inc"

         !----------------------------------------------------------------------
         !  Type-bound procedures
         !----------------------------------------------------------------------
         CONTAINS

#define  DTYPE(a) a/**/_s
#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_single
#include "map/mapping_typeproc.f"

#define  DTYPE(a) a/**/_d
#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_double
#include "map/mapping_typeproc.f"

      END MODULE ppm_module_mapping_typedef
