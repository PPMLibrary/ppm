
!----------------------------------------------------------------------
!  Generic mapping data type
!----------------------------------------------------------------------
TYPE,EXTENDS(DTYPE(ppm_t_mapping)_) :: DTYPE(ppm_t_mapping)
END TYPE

TYPE,EXTENDS(DTYPE(ppm_t_part_mapping)_) :: DTYPE(ppm_t_part_mapping)
    CONTAINS
    PROCEDURE     :: create => DTYPE(map_create)
    PROCEDURE     :: destroy => DTYPE(map_destroy)
END TYPE 
!----------------------------------------------------------------------
! Container for particle mappings
!----------------------------------------------------------------------
#define CONTAINER DTYPE(ppm_c_pmaps)
#define __CONTAINER(a) DTYPE(ppm_c_pmaps)_/**/a
#define VEC_TYPE DTYPE(ppm_t_part_mapping)
#include "cont/collection_template.inc"


TYPE,EXTENDS(DTYPE(ppm_t_mesh_mapping)_) :: DTYPE(ppm_t_mesh_mapping)

    !             
    !             PROCEDURE       :: mesh_send
    !             GENERIC, PUBLIC ::  mesh_pop => &
    !                    DTYPE(ppm_map_mesh_pop_1d),&
    !                    DTYPE(ppm_map_mesh_pop_1dc),&
    !                    ppm_map_mesh_pop_1di,&
    !                    ppm_map_mesh_pop_1dl,&
    !                    DTYPE(ppm_map_mesh_pop_2d),&
    !                    DTYPE(ppm_map_mesh_pop_2dc),&
    !                    ppm_map_mesh_pop_2di,&
    !                    ppm_map_mesh_pop_2dl
    !             GENERIC, PUBLIC ::  mesh_push => &
    !                    DTYPE(ppm_map_mesh_push_1d),&
    !                    DTYPE(ppm_map_mesh_push_1dc),&
    !                    ppm_map_mesh_push_1di,&
    !                    ppm_map_mesh_push_1dl,&
    !                    DTYPE(ppm_map_mesh_push_2d),&
    !                    DTYPE(ppm_map_mesh_push_2dc),&
    !                    ppm_map_mesh_push_2di,&
    !                    ppm_map_mesh_push_2dl
END TYPE
#undef   MK
#undef   _MK
#undef   DTYPE
