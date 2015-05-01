#if __MYTYPE == __MAPPINGTYPE
      !----------------------------------------------------------------------
      !  Generic mapping data type
      !----------------------------------------------------------------------
      TYPE,EXTENDS(ppm_t_mapping_) :: ppm_t_mapping
      END TYPE
#endif

      TYPE,EXTENDS(DTYPE(ppm_t_part_mapping)_) :: DTYPE(ppm_t_part_mapping)
      CONTAINS
          PROCEDURE :: create  => DTYPE(map_part_create)
          PROCEDURE :: destroy => DTYPE(map_part_destroy)
      END TYPE
      !----------------------------------------------------------------------
      ! Container for particle mappings
      !----------------------------------------------------------------------
minclude ppm_create_collection(DTYPE(part_mapping),DTYPE(part_mapping),generate="extend")

#if __MYTYPE == __MAPPINGTYPE
      TYPE,EXTENDS(ppm_t_mesh_mapping_) :: ppm_t_mesh_mapping
      CONTAINS
          PROCEDURE :: create  => map_mesh_create
          PROCEDURE :: destroy => map_mesh_destroy
!          PROCEDURE :: mesh_send
!          GENERIC :: mesh_pop => &
!             &          DTYPE(ppm_map_mesh_pop_1d),  &
!             &          DTYPE(ppm_map_mesh_pop_1dc), &
!             &          ppm_map_mesh_pop_1di,        &
!             &          ppm_map_mesh_pop_1dl,        &
!             &          DTYPE(ppm_map_mesh_pop_2d),  &
!             &          DTYPE(ppm_map_mesh_pop_2dc), &
!             &          ppm_map_mesh_pop_2di,        &
!             &          ppm_map_mesh_pop_2dl
!             GENERIC :: mesh_push =>                 &
!             &          DTYPE(ppm_map_mesh_push_1d), &
!             &          DTYPE(ppm_map_mesh_push_1dc),&
!             &          ppm_map_mesh_push_1di,       &
!             &          ppm_map_mesh_push_1dl,       &
!             &          DTYPE(ppm_map_mesh_push_2d), &
!             &          DTYPE(ppm_map_mesh_push_2dc),&
!             &          ppm_map_mesh_push_2di,       &
!             &          ppm_map_mesh_push_2dl
      END TYPE
      !----------------------------------------------------------------------
      ! Container for mesh mappings
      !----------------------------------------------------------------------
minclude ppm_create_collection(mesh_mapping,mesh_mapping,generate="extend")

#endif
#undef __MYTYPE

#undef   MK
#undef  _MK
#undef  DTYPE
