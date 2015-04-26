      minclude ppm_create_collection_interfaces(DTYPE(part_mapping)_,DTYPE(part_mapping)_)

      SUBROUTINE DTYPE(map_part_create)_(map,source_topoid,target_topoid,info)
          !!! Constructor for particle mapping data structure
          IMPORT :: DTYPE(ppm_t_part_mapping)_
          IMPLICIT NONE
          CLASS(DTYPE(ppm_t_part_mapping)_) :: map
          INTEGER,            INTENT(IN   ) :: source_topoid
          INTEGER,            INTENT(IN   ) :: target_topoid
          INTEGER,            INTENT(  OUT) :: info
      END SUBROUTINE DTYPE(map_part_create)_

      SUBROUTINE DTYPE(map_part_destroy)_(map,info)
          IMPORT :: DTYPE(ppm_t_part_mapping)_
          IMPLICIT NONE
          CLASS(DTYPE(ppm_t_part_mapping)_) :: map
          INTEGER,            INTENT(  OUT) :: info
      END SUBROUTINE DTYPE(map_part_destroy)_
#if __MYTYPE == __MAPPINGTYPE
      minclude ppm_create_collection_interfaces(mesh_mapping_,mesh_mapping_)

      SUBROUTINE map_mesh_create_(map,source_topoid,target_topoid,info)
          !!! Constructor for particle mapping data structure
          IMPORT :: ppm_t_mesh_mapping_
          IMPLICIT NONE
          CLASS(ppm_t_mesh_mapping_) :: map
          INTEGER,     INTENT(IN   ) :: source_topoid
          INTEGER,     INTENT(IN   ) :: target_topoid
          INTEGER,     INTENT(  OUT) :: info
      END SUBROUTINE map_mesh_create_

      SUBROUTINE map_mesh_destroy_(map,info)
          IMPORT :: ppm_t_mesh_mapping_
          IMPLICIT NONE
          CLASS(ppm_t_mesh_mapping_) :: map
          INTEGER,     INTENT(  OUT) :: info
      END SUBROUTINE map_mesh_destroy_
#endif

#undef __MYTYPE
#undef DTYPE
#undef DEFINE_MK
