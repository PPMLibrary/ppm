minclude ppm_create_collection_interfaces(DTYPE(part_mapping)_,DTYPE(part_mapping)_)

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
