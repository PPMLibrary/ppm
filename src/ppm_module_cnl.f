      module ppm_module_cnl

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2
      use ppm_module_data


      interface cnl_clist
          module procedure cnl_clist_s
          module procedure cnl_clist_d
      end interface cnl_clist
      
      interface cnl_rank2d
          module procedure cnl_rank2d_s
          module procedure cnl_rank2d_d
      end interface cnl_rank2d
      
      interface cnl_rank3d
          module procedure cnl_rank3d_s
          module procedure cnl_rank3d_d
      end interface cnl_rank3d

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



      contains
#define __KIND __SINGLE_PRECISION
#include "neighlist/cnl_clist.f"
#include "neighlist/cnl_rank2d.f"
#include "neighlist/cnl_rank3d.f"
#undef  __KIND
#define __KIND __DOUBLE_PRECISION
#include "neighlist/cnl_clist.f"
#include "neighlist/cnl_rank2d.f"
#include "neighlist/cnl_rank3d.f"
#undef  __KIND

#include "neighlist/cnl_MkNeighIdx.f"

#define __DIM 2
#include "neighlist/cnl_vlist_build.f"
#define __DIM 3
#include "neighlist/cnl_vlist_build.f"

#include "neighlist/cnl_vlist.f"

      end module ppm_module_cnl
