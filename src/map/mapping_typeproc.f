      SUBROUTINE DTYPE(map_part_create)(map,source_topoid,target_topoid,info)
        !!! Constructor for particle mapping data structure
        IMPLICIT NONE

        DEFINE_MK()
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(DTYPE(ppm_t_part_mapping)) :: map

        INTEGER,           INTENT(IN   ) :: source_topoid
        INTEGER,           INTENT(IN   ) :: target_topoid
        INTEGER,           INTENT(  OUT) :: info

        start_subroutine("map_part_create")

        map%source_topoid = source_topoid
        map%target_topoid = target_topoid

        end_subroutine()
      END SUBROUTINE DTYPE(map_part_create)

      SUBROUTINE DTYPE(map_part_destroy)(map,info)
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(DTYPE(ppm_t_part_mapping)) :: map

        INTEGER,           INTENT(  OUT) :: info
        !!! Returns status, 0 upon success.

        start_subroutine("map_part_destroy")

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

        end_subroutine()
      END SUBROUTINE DTYPE(map_part_destroy)

minclude ppm_create_collection_procedures(DTYPE(part_mapping),DTYPE(part_mapping)_)

#if __MYTYPE == __MAPPINGTYPE
      SUBROUTINE map_mesh_create(map,source_topoid,target_topoid,info)
        !!! Constructor for particle mapping data structure
        IMPLICIT NONE

        DEFINE_MK()
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(ppm_t_mesh_mapping) :: map

        INTEGER,    INTENT(IN   ) :: source_topoid
        INTEGER,    INTENT(IN   ) :: target_topoid
        INTEGER,    INTENT(  OUT) :: info

        start_subroutine("map_mesh_create")

        map%source_topoid = source_topoid
        map%target_topoid = target_topoid

        end_subroutine()
      END SUBROUTINE map_mesh_create

      SUBROUTINE map_mesh_destroy(map,info)
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(ppm_t_mesh_mapping) :: map

        INTEGER,    INTENT(  OUT) :: info
        !!! Returns status, 0 upon success.

        start_subroutine("map_mesh_destroy")

        map%source_topoid = -1
        map%target_topoid = -1

        dealloc_pointer(map%nsend)
        dealloc_pointer(map%nrecv)
        dealloc_pointer(map%psend)
        dealloc_pointer(map%precv)
        dealloc_pointer(map%pp)
        dealloc_pointer(map%qq)

        dealloc_pointer(map%ghostsize)
        dealloc_pointer(map%ghost_fromsub)
        dealloc_pointer(map%ghost_tosub)
        dealloc_pointer(map%ghost_patchid)
        dealloc_pointer(map%ghost_blkstart)
        dealloc_pointer(map%ghost_blksize)
        dealloc_pointer(map%ghost_blk)

        map%ghost_nsend = 0
        map%ghost_nrecv = 0

        dealloc_pointer(map%ghost_recvtosub)
        dealloc_pointer(map%ghost_recvpatchid)
        dealloc_pointer(map%ghost_recvblkstart)
        dealloc_pointer(map%ghost_recvblksize)
        dealloc_pointer(map%ghost_recvblk)

        end_subroutine()
      END SUBROUTINE map_mesh_destroy

      !
minclude ppm_create_collection_procedures(mesh_mapping,mesh_mapping_)

#endif
#undef __MYTYPE

#undef DTYPE
#undef DEFINE_MK
