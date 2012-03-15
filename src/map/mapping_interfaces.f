#define CONTAINER DTYPE(ppm_c_pmaps)_
#define __CONTAINER(a) DTYPE(ppm_c_pmaps)__/**/a
#define VEC_TYPE DTYPE(ppm_t_part_mapping)_
#include "cont/collection_interfaces.f"


SUBROUTINE DTYPE(map_create)_(map,source_topoid,target_topoid,info)
    !!! Constructor for particle mapping data structure
    IMPORT DTYPE(ppm_t_part_mapping)_
    CLASS(DTYPE(ppm_t_part_mapping)_)   :: map
    INTEGER,                INTENT(IN) :: source_topoid
    INTEGER,                INTENT(IN) :: target_topoid
    INTEGER,               INTENT(OUT) :: info
END SUBROUTINE DTYPE(map_create)_

SUBROUTINE DTYPE(map_destroy)_(map,info)
    IMPORT DTYPE(ppm_t_part_mapping)_
    CLASS(DTYPE(ppm_t_part_mapping)_)   :: map
    INTEGER,            INTENT(  OUT)  :: info
END SUBROUTINE DTYPE(map_destroy)_


#undef DTYPE
#undef DEFINE_MK
