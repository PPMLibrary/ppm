      MODULE ppm_module_cnl

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2
      USE ppm_module_data

#define interface_s_d(a) \
      INTERFACE a ;\
          MODULE PROCEDURE a/**/_s; \
          MODULE PROCEDURE a/**/_d; \
      END INTERFACE

      interface_s_d(cnl_clist)
      interface_s_d(cnl_vlist)
      interface_s_d(cnl_vlist_build_2d)
      interface_s_d(cnl_vlist_build_3d)
      interface_s_d(cnl_rank2d)
      interface_s_d(cnl_rank3d)

      TYPE t_clist
          !!! Cell list data structure
          INTEGER(ppm_kind_int64), DIMENSION(:), POINTER    :: nm  => NULL()
          !!! Number of cells in x,y,(z) direction (including the ghosts 
          !!! cells) in each subdomain. 
          INTEGER(ppm_kind_int64), DIMENSION(:), POINTER    :: lpdx => NULL()
          !!! particle index list
          INTEGER(ppm_kind_int64), DIMENSION(:), POINTER    :: lhbx => NULL()
          !!! first particle in each cell
      END TYPE

      PUBLIC cnl_vlist


      CONTAINS

#define  __KIND __SINGLE_PRECISION
#define  DTYPE(a) a/**/_s
#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_single
#include "neighlist/cnl_clist.f"
#include "neighlist/cnl_vlist.f"
#include "neighlist/cnl_rank2d.f"
#include "neighlist/cnl_rank3d.f"
#define __DIM 2
#include "neighlist/cnl_vlist_build.f"
#define __DIM 3
#include "neighlist/cnl_vlist_build.f"
#undef   DEFINE_MK
#undef   DTYPE
#undef   __KIND

#define  __KIND __DOUBLE_PRECISION
#define  DTYPE(a) a/**/_d
#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_double
#include "neighlist/cnl_clist.f"
#include "neighlist/cnl_vlist.f"
#include "neighlist/cnl_rank2d.f"
#include "neighlist/cnl_rank3d.f"
#define __DIM 2
#include "neighlist/cnl_vlist_build.f"
#define __DIM 3
#include "neighlist/cnl_vlist_build.f"
#undef   DEFINE_MK
#undef   DTYPE
#undef   __KIND

#include "neighlist/cnl_MkNeighIdx.f"


      END MODULE ppm_module_cnl
