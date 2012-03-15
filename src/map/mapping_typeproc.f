#define CONTAINER DTYPE(ppm_c_pmaps)
#define __CONTAINER(a) DTYPE(ppm_c_pmaps)_/**/a
#define VEC_TYPE DTYPE(ppm_t_part_mapping)
#include "cont/collection_typeproc.f"


SUBROUTINE DTYPE(map_create)(map,source_topoid,target_topoid,info)
    !!! Constructor for particle mapping data structure
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_part_mapping))   :: map
    INTEGER,                INTENT(IN) :: source_topoid
    INTEGER,                INTENT(IN) :: target_topoid
    INTEGER,               INTENT(OUT) :: info

    REAL(KIND(1.D0))                   :: t0
    CHARACTER(LEN=ppm_char)            :: caller = 'map_create'



    CALL substart(caller,t0,info)

    map%source_topoid = source_topoid
    map%target_topoid = target_topoid

    CALL substop(caller,t0,info)

    9999  CONTINUE

END SUBROUTINE DTYPE(map_create)

SUBROUTINE DTYPE(map_destroy)(map,info)
    CLASS(DTYPE(ppm_t_part_mapping))   :: map
    INTEGER,            INTENT(  OUT)  :: info
    !!! Returns status, 0 upon success.
    REAL(KIND(1.D0))                   :: t0
    CHARACTER(LEN=ppm_char)            :: caller = 'map_destroy'


    CALL substart(caller,t0,info)

    map%source_topoid = -1
    map%target_topoid = -1

    dealloc_pointer(map%send)
    dealloc_pointer(map%recv)
    dealloc_pointer(map%nsend)
    dealloc_pointer(map%nrecv)
    dealloc_pointer(map%psend)
    dealloc_pointer(map%precv)
    dealloc_pointer(map%pp)
    dealloc_pointer(map%qq)

    map%oldNpart = 0
    map%newNpart = 0

    CALL substop(caller,t0,info)
    9999  CONTINUE

END SUBROUTINE DTYPE(map_destroy)


#undef DTYPE
#undef DEFINE_MK
