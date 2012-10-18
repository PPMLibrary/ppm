minclude ppm_create_collection_interfaces(DTYPE(part_prop)_,DTYPE(part_prop)_)
minclude ppm_create_collection_interfaces(DTYPE(neighlist)_,DTYPE(neighlist)_)
minclude ppm_create_collection_interfaces(DTYPE(particles)_,DTYPE(particles)_)
!minclude ppm_create_collection_interfaces(DTYPE(ppm_t_sop)_)

!CREATE ENTRY
SUBROUTINE DTYPE(prop_create)_(prop,datatype,parts,npart,lda,name,flags,info,field,zero)
    !!! Constructor for particle property data structure
    IMPORT DTYPE(ppm_t_part_prop)_,ppm_param_length_pptflags,ppm_t_field_
    IMPORT ppm_t_discr_kind 
    CLASS(DTYPE(ppm_t_part_prop)_)      :: prop
    INTEGER,                INTENT(IN) :: datatype
    CLASS(ppm_t_discr_kind), TARGET, INTENT(IN) :: parts
    INTEGER,                INTENT(IN) :: npart
    INTEGER,                INTENT(IN) :: lda
    CHARACTER(LEN=*),       INTENT(IN) :: name
    LOGICAL, DIMENSION(ppm_param_length_pptflags),INTENT(IN) :: flags
    INTEGER,               INTENT(OUT) :: info
    CLASS(ppm_t_field_),OPTIONAL,TARGET, INTENT(IN) :: field
    LOGICAL, OPTIONAL,     INTENT( IN) :: zero
END SUBROUTINE
!DESTROY ENTRY
SUBROUTINE DTYPE(prop_destroy)_(prop,info)
    IMPORT DTYPE(ppm_t_part_prop)_
    CLASS(DTYPE(ppm_t_part_prop)_)      :: prop
    INTEGER,                                INTENT(  OUT)  :: info
END SUBROUTINE

SUBROUTINE DTYPE(prop_print_info)_(prop,info,level,fileunit,propid)
    IMPORT DTYPE(ppm_t_part_prop)_
    CLASS(DTYPE(ppm_t_part_prop)_)                          :: prop
    INTEGER,                            INTENT(  OUT)      :: info
    INTEGER, OPTIONAL,                  INTENT(IN   )      :: level
    INTEGER, OPTIONAL,                  INTENT(IN   )      :: fileunit
    INTEGER, OPTIONAL,                  INTENT(IN   )      :: propid
END SUBROUTINE

SUBROUTINE DTYPE(get_xp)_(this,xp,info,with_ghosts)
    IMPORT DTYPE(ppm_t_particles)_, MK
    CLASS(DTYPE(ppm_t_particles)_)               :: this
    REAL(MK),DIMENSION(:,:),POINTER,INTENT(OUT)  :: xp
    INTEGER,                        INTENT(OUT)  :: info
    LOGICAL,OPTIONAL                             :: with_ghosts
END SUBROUTINE

SUBROUTINE DTYPE(set_xp)_(this,xp,info,read_only,ghosts_ok)
    IMPORT DTYPE(ppm_t_particles)_, MK
    CLASS(DTYPE(ppm_t_particles)_)   :: this
    LOGICAL,OPTIONAL                 :: read_only
    LOGICAL,OPTIONAL                 :: ghosts_ok
    REAL(MK),DIMENSION(:,:),POINTER  :: xp
    INTEGER,            INTENT(OUT)  :: info

END SUBROUTINE

SUBROUTINE DTYPE(part_prop_create)_(this,info,field,part_prop,discr_data,&
        dtype,name,lda,zero,with_ghosts)
    IMPORT DTYPE(ppm_t_particles)_,ppm_t_field_,&
        DTYPE(ppm_t_part_prop)_,ppm_t_discr_data
    !!! Adds a property to an existing particle set
    CLASS(DTYPE(ppm_t_particles)_)        :: this
    INTEGER,               INTENT(OUT)    :: info
    CLASS(ppm_t_field_),OPTIONAL,INTENT(IN   ) :: field
    !CLASS(DTYPE(ppm_t_part_prop)_),OPTIONAL,POINTER,INTENT(OUT):: discr_data
    CLASS(DTYPE(ppm_t_part_prop)_),OPTIONAL,POINTER,INTENT(OUT):: part_prop
    CLASS(ppm_t_discr_data),OPTIONAL,POINTER,INTENT(OUT):: discr_data
    INTEGER, OPTIONAL,      INTENT(IN   ) :: dtype
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN ) :: name
    INTEGER, OPTIONAL,      INTENT(IN   ) :: lda
    LOGICAL, OPTIONAL                     :: zero
    LOGICAL, OPTIONAL                     :: with_ghosts
END SUBROUTINE

SUBROUTINE DTYPE(part_prop_destroy)_(this,prop,info)
    IMPORT DTYPE(ppm_t_particles)_,ppm_t_discr_data
    !!! Destroy a property from an existing particle set
    CLASS(DTYPE(ppm_t_particles)_)         :: this
    CLASS(ppm_t_discr_data), INTENT(INOUT) :: prop
    INTEGER,                INTENT(OUT)    :: info
END SUBROUTINE

SUBROUTINE DTYPE(part_prop_realloc)_(Pc,prop,info,with_ghosts,datatype,lda)
    IMPORT DTYPE(ppm_t_particles)_,DTYPE(ppm_t_part_prop)_
    !!! Reallocate the property array to the correct size
    !!! (e.g. if the number of particles has changed or if the type
    !!! of the data changes)
    CLASS(DTYPE(ppm_t_particles)_)        :: Pc
    CLASS(DTYPE(ppm_t_part_prop)_)        :: prop
    INTEGER,               INTENT(OUT)    :: info
    LOGICAL, OPTIONAL                     :: with_ghosts
    INTEGER, OPTIONAL                     :: datatype
    INTEGER, OPTIONAL                     :: lda
END SUBROUTINE


SUBROUTINE DTYPE(part_neigh_create)_(this,Part_src,info,&
        name,skin,symmetry,cutoff,Nlist)
    IMPORT DTYPE(ppm_t_particles)_,MK,DTYPE(ppm_t_neighlist)_
    CLASS(DTYPE(ppm_t_particles)_)                     :: this
    CLASS(DTYPE(ppm_t_particles)_),TARGET, INTENT(IN)  :: Part_src
    INTEGER,               INTENT(OUT)                 :: info
    CHARACTER(LEN=*) , OPTIONAL                        :: name
    REAL(MK), OPTIONAL                                 :: skin
    REAL(MK), OPTIONAL                                 :: cutoff
    LOGICAL, OPTIONAL                                  :: symmetry    
    CLASS(DTYPE(ppm_t_neighlist)_),POINTER,OPTIONAL,INTENT(OUT) :: Nlist
END SUBROUTINE

SUBROUTINE DTYPE(part_neigh_destroy)_(this,Nlist,info)
    IMPORT DTYPE(ppm_t_particles)_,DTYPE(ppm_t_neighlist)_
    !!! Destroy a property from an existing particle set
    CLASS(DTYPE(ppm_t_particles)_)                       :: this
    CLASS(DTYPE(ppm_t_neighlist)_),POINTER,INTENT(INOUT) :: Nlist
    INTEGER,                               INTENT(OUT)   :: info
END SUBROUTINE


SUBROUTINE DTYPE(part_create)_(Pc,Npart,info,name)
    IMPORT DTYPE(ppm_t_particles)_
    CLASS(DTYPE(ppm_t_particles)_)                         :: Pc
    INTEGER,                                INTENT(IN   )  :: Npart
    INTEGER,                                INTENT(  OUT)  :: info
    CHARACTER(LEN=*) , OPTIONAL                            :: name
END SUBROUTINE


SUBROUTINE DTYPE(part_destroy)_(Pc,info)
    IMPORT DTYPE(ppm_t_particles)_
    CLASS(DTYPE(ppm_t_particles)_)                          :: Pc
    INTEGER,                                INTENT(  OUT)  :: info
END SUBROUTINE

SUBROUTINE DTYPE(particles_initialize2d)_(Pc,Npart_global,info,&
        distrib,topoid,minphys,maxphys,cutoff,name)
    IMPORT DTYPE(ppm_t_particles)_,MK,ppm_dim
    CLASS(DTYPE(ppm_t_particles)_)                          :: Pc
    INTEGER,                             INTENT(INOUT)     :: Npart_global
    INTEGER,                             INTENT(  OUT)     :: info
    INTEGER,OPTIONAL,                    INTENT(IN   )     :: distrib
    INTEGER,OPTIONAL,                    INTENT(IN   )     :: topoid
    REAL(MK),DIMENSION(ppm_dim),OPTIONAL,INTENT(IN   )     :: minphys
    REAL(MK),DIMENSION(ppm_dim),OPTIONAL,INTENT(IN   )     :: maxphys
    REAL(MK),                   OPTIONAL,INTENT(IN   )     :: cutoff
    CHARACTER(LEN=*),           OPTIONAL,INTENT(IN   )     :: name
END SUBROUTINE

SUBROUTINE DTYPE(particles_initialize3d)_(Pc,Npart_global,info,&
        distrib,topoid,minphys,maxphys,cutoff,name)
    IMPORT DTYPE(ppm_t_particles)_,MK,ppm_dim
    CLASS(DTYPE(ppm_t_particles)_)                          :: Pc
    INTEGER,                             INTENT(INOUT)     :: Npart_global
    INTEGER,                             INTENT(  OUT)     :: info
    INTEGER,OPTIONAL,                    INTENT(IN   )     :: distrib
    INTEGER,OPTIONAL,                    INTENT(IN   )     :: topoid
    REAL(MK),DIMENSION(ppm_dim),OPTIONAL,INTENT(IN   )     :: minphys
    REAL(MK),DIMENSION(ppm_dim),OPTIONAL,INTENT(IN   )     :: maxphys
    REAL(MK),                   OPTIONAL,INTENT(IN   )     :: cutoff
    CHARACTER(LEN=*),           OPTIONAL,INTENT(IN   )     :: name
END SUBROUTINE

!!temporary hack to deal with both 2d and 3d
SUBROUTINE DTYPE(part_initialize)_(Pc,Npart_global,info,&
        distrib,topoid,minphys,maxphys,cutoff,name)
    IMPORT DTYPE(ppm_t_particles)_,MK,ppm_dim
    CLASS(DTYPE(ppm_t_particles)_)                          :: Pc
    INTEGER,                            INTENT(INOUT)      :: Npart_global
    INTEGER,                            INTENT(  OUT)      :: info
    INTEGER,OPTIONAL,                    INTENT(IN   )     :: distrib
    INTEGER,OPTIONAL,                    INTENT(IN   )     :: topoid
    REAL(MK),DIMENSION(ppm_dim),OPTIONAL,INTENT(IN   )     :: minphys
    REAL(MK),DIMENSION(ppm_dim),OPTIONAL,INTENT(IN   )     :: maxphys
    REAL(MK),                   OPTIONAL,INTENT(IN   )     :: cutoff
    CHARACTER(LEN=*),           OPTIONAL,INTENT(IN   )     :: name
END SUBROUTINE

SUBROUTINE DTYPE(part_print_info)_(Pc,info,level,fileunit)
    IMPORT DTYPE(ppm_t_particles)_
    CLASS(DTYPE(ppm_t_particles)_)                          :: Pc
    INTEGER,                            INTENT(  OUT)      :: info
    INTEGER, OPTIONAL,                  INTENT(IN   )      :: level
    INTEGER, OPTIONAL,                  INTENT(IN   )      :: fileunit
END SUBROUTINE

SUBROUTINE DTYPE(part_del_parts)_(Pc,list_del_parts,nb_del,info)
    IMPORT DTYPE(ppm_t_particles)_
    CLASS(DTYPE(ppm_t_particles)_)                         :: Pc
    INTEGER,DIMENSION(:),POINTER,           INTENT(IN   )  :: list_del_parts
    INTEGER,                                INTENT(IN   )  :: nb_del
    INTEGER,                                INTENT(  OUT)  :: info
END SUBROUTINE

SUBROUTINE DTYPE(part_prop_push)_(Pc,prop_id,info)
    IMPORT DTYPE(ppm_t_particles)_
    CLASS(DTYPE(ppm_t_particles)_)                          :: Pc
    INTEGER,                                INTENT(IN   )  :: prop_id
    INTEGER,                                INTENT(  OUT)  :: info
END SUBROUTINE

SUBROUTINE DTYPE(part_prop_pop)_(Pc,prop_id,Npart_new,info)
    IMPORT DTYPE(ppm_t_particles)_
    CLASS(DTYPE(ppm_t_particles)_)                          :: Pc
    INTEGER,                                INTENT(IN   )  :: prop_id
    INTEGER,                                INTENT(IN   )  :: Npart_new
    INTEGER,                                INTENT(  OUT)  :: info
END SUBROUTINE

FUNCTION DTYPE(has_neighlist)_(this,Part) RESULT(res)
    IMPORT DTYPE(ppm_t_particles)_,ppm_t_discr_kind
    CLASS(DTYPE(ppm_t_particles)_),TARGET          :: this
    CLASS(ppm_t_discr_kind),OPTIONAL,TARGET        :: Part
    LOGICAL                                        :: res
END FUNCTION

FUNCTION DTYPE(get_neighlist)_(this,Part) RESULT(NList)
    IMPORT DTYPE(ppm_t_particles)_,DTYPE(ppm_t_neighlist)_,ppm_t_discr_kind
    CLASS(DTYPE(ppm_t_particles)_),TARGET          :: this
    CLASS(ppm_t_discr_kind),OPTIONAL,TARGET        :: Part
    CLASS(DTYPE(ppm_t_neighlist)_),POINTER         :: NList
END FUNCTION

SUBROUTINE DTYPE(get_vlist)_(this,nvlist,vlist,info,NList)
    IMPORT DTYPE(ppm_t_particles)_,DTYPE(ppm_t_neighlist)_
    CLASS(DTYPE(ppm_t_particles)_)               :: this
    INTEGER,DIMENSION(:),POINTER,  INTENT( OUT)  :: nvlist
    INTEGER,DIMENSION(:,:),POINTER,INTENT( OUT)  :: vlist
    INTEGER,                       INTENT(INOUT) :: info
    CLASS(DTYPE(ppm_t_neighlist)_),OPTIONAL,TARGET :: NList
END SUBROUTINE

SUBROUTINE DTYPE(get_nvlist)_(this,nvlist,info,NList)
    IMPORT DTYPE(ppm_t_particles)_,DTYPE(ppm_t_neighlist)_
    CLASS(DTYPE(ppm_t_particles)_)               :: this
    INTEGER,DIMENSION(:),POINTER,  INTENT( OUT)  :: nvlist
    INTEGER,                       INTENT(INOUT) :: info
    CLASS(DTYPE(ppm_t_neighlist)_),OPTIONAL,TARGET :: NList
END SUBROUTINE

SUBROUTINE DTYPE(part_map)_(Pc,info,global,topoid)
    IMPORT DTYPE(ppm_t_particles)_
    CLASS(DTYPE(ppm_t_particles)_)                          :: Pc
    INTEGER,                            INTENT(  OUT)      :: info
    LOGICAL, OPTIONAL                                   :: global
    INTEGER, OPTIONAL                                   :: topoid
END SUBROUTINE

SUBROUTINE DTYPE(part_map_positions)_(Pc,info,global,topoid)
    IMPORT DTYPE(ppm_t_particles)_
    CLASS(DTYPE(ppm_t_particles)_)                      :: Pc
    INTEGER,                            INTENT(  OUT)   :: info
    LOGICAL, OPTIONAL                                   :: global
    INTEGER, OPTIONAL                                   :: topoid
END SUBROUTINE

SUBROUTINE DTYPE(part_map_push)_(Pc,info,Field)
    IMPORT DTYPE(ppm_t_particles)_,ppm_t_field_
    CLASS(DTYPE(ppm_t_particles)_)                      :: Pc
    INTEGER,                            INTENT(  OUT)   :: info
    CLASS(ppm_t_field_),OPTIONAL                        :: Field
END SUBROUTINE

SUBROUTINE DTYPE(part_map_send)_(Pc,info)
    IMPORT DTYPE(ppm_t_particles)_
    CLASS(DTYPE(ppm_t_particles)_)                      :: Pc
    INTEGER,                            INTENT(  OUT)   :: info
END SUBROUTINE

SUBROUTINE DTYPE(part_map_pop)_(Pc,info,Field)
    IMPORT DTYPE(ppm_t_particles)_,ppm_t_field_
    CLASS(DTYPE(ppm_t_particles)_)                      :: Pc
    INTEGER,                            INTENT(  OUT)   :: info
    CLASS(ppm_t_field_),OPTIONAL                        :: Field
END SUBROUTINE

SUBROUTINE DTYPE(part_map_pop_positions)_(Pc,info)
    IMPORT DTYPE(ppm_t_particles)_
    CLASS(DTYPE(ppm_t_particles)_)                      :: Pc
    INTEGER,                            INTENT(  OUT)   :: info
END SUBROUTINE

SUBROUTINE DTYPE(part_map_ghosts)_(Pc,info,ghostsize)
    IMPORT DTYPE(ppm_t_particles)_,MK
    CLASS(DTYPE(ppm_t_particles)_)                       :: Pc
    INTEGER,                            INTENT(  OUT)   :: info
    REAL(MK), OPTIONAL                                  :: ghostsize
END SUBROUTINE

SUBROUTINE DTYPE(part_map_ghost_get)_(Pc,info,ghostsize)
    IMPORT DTYPE(ppm_t_particles)_,MK,ppm_t_field_
    CLASS(DTYPE(ppm_t_particles)_)                      :: Pc
    INTEGER,                            INTENT(  OUT)   :: info
    REAL(MK), OPTIONAL                                  :: ghostsize
END SUBROUTINE

SUBROUTINE DTYPE(part_map_ghost_push)_(Pc,info,Field,ghostsize)
    IMPORT DTYPE(ppm_t_particles)_,MK,ppm_t_field_
    CLASS(DTYPE(ppm_t_particles)_)                      :: Pc
    INTEGER,                            INTENT(  OUT)   :: info
    CLASS(ppm_t_field_),OPTIONAL                        :: Field
    REAL(MK), OPTIONAL                                  :: ghostsize
END SUBROUTINE

SUBROUTINE DTYPE(part_map_ghost_send)_(Pc,info)
    IMPORT DTYPE(ppm_t_particles)_,ppm_t_field_
    CLASS(DTYPE(ppm_t_particles)_)                      :: Pc
    INTEGER,                            INTENT(  OUT)   :: info
END SUBROUTINE

SUBROUTINE DTYPE(part_map_ghost_pop)_(Pc,info,Field,ghostsize)
    IMPORT DTYPE(ppm_t_particles)_,MK,ppm_t_field_
    CLASS(DTYPE(ppm_t_particles)_)                      :: Pc
    INTEGER,                            INTENT(  OUT)   :: info
    CLASS(ppm_t_field_),OPTIONAL                        :: Field
    REAL(MK), OPTIONAL                                  :: ghostsize
END SUBROUTINE

SUBROUTINE DTYPE(part_map_ghost_pop_pos)_(Pc,info,ghostsize)
    IMPORT DTYPE(ppm_t_particles)_,MK,ppm_t_field_
    CLASS(DTYPE(ppm_t_particles)_)                      :: Pc
    INTEGER,                            INTENT(  OUT)   :: info
    REAL(MK), OPTIONAL                                  :: ghostsize
END SUBROUTINE

SUBROUTINE DTYPE(part_apply_bc)_(Pc,info)
    IMPORT DTYPE(ppm_t_particles)_
    CLASS(DTYPE(ppm_t_particles)_)                         :: Pc
    INTEGER,                            INTENT(  OUT)     :: info
END SUBROUTINE

SUBROUTINE DTYPE(part_move)_(Pc,disp,info)
    IMPORT DTYPE(ppm_t_particles)_, MK
    CLASS(DTYPE(ppm_t_particles)_)                         :: Pc
    REAL(MK), DIMENSION(:,:), POINTER,  INTENT(IN   )     :: disp
    INTEGER,                            INTENT(  OUT)     :: info
END SUBROUTINE

SUBROUTINE DTYPE(part_neighlist)_(this,info,P_xset,name,skin,symmetry,cutoff,&
        lstore,incl_ghosts,knn)
    IMPORT DTYPE(ppm_t_particles)_,MK
    CLASS(DTYPE(ppm_t_particles)_),TARGET                  :: this
    INTEGER,                            INTENT(  OUT)      :: info
    CLASS(DTYPE(ppm_t_particles)_),OPTIONAL,TARGET         :: P_xset
    REAL(MK), OPTIONAL                                     :: skin
    LOGICAL, OPTIONAL                                      :: symmetry    
    REAL(MK), OPTIONAL                                     :: cutoff
    CHARACTER(LEN=*) , OPTIONAL                            :: name
    LOGICAL, OPTIONAL,                  INTENT(IN   )      :: lstore
    LOGICAL, OPTIONAL,                  INTENT(IN   )      :: incl_ghosts
    INTEGER, OPTIONAL,                  INTENT(IN   )      :: knn
END SUBROUTINE

SUBROUTINE DTYPE(part_set_cutoff)_(Pc,cutoff,info,Nlist)
    IMPORT DTYPE(ppm_t_particles)_,MK,DTYPE(ppm_t_neighlist)_
    CLASS(DTYPE(ppm_t_particles)_)            :: Pc
    REAL(MK),                 INTENT(IN   )  :: cutoff
    INTEGER,                  INTENT(   OUT) :: info
    CLASS(DTYPE(ppm_t_neighlist)_),OPTIONAL,INTENT(INOUT) :: Nlist
END SUBROUTINE

SUBROUTINE DTYPE(part_comp_global_index)_(Pc,info)
    IMPORT DTYPE(ppm_t_particles)_
    CLASS(DTYPE(ppm_t_particles)_)            :: Pc
    INTEGER,                  INTENT(   OUT) :: info
END SUBROUTINE


SUBROUTINE DTYPE(part_map_create)_(Pc,id,source_topoid,target_topoid,info)
    IMPORT DTYPE(ppm_t_particles)_
    !!! Adds a property to an existing particle set
    CLASS(DTYPE(ppm_t_particles)_)         :: Pc
    INTEGER,                INTENT(  OUT) :: id
    INTEGER,                INTENT(IN   ) :: source_topoid
    INTEGER,                INTENT(IN   ) :: target_topoid
    INTEGER,               INTENT(OUT)    :: info

END SUBROUTINE 

SUBROUTINE DTYPE(part_map_destroy)_(Pc,id,info)
    IMPORT DTYPE(ppm_t_particles)_
    CLASS(DTYPE(ppm_t_particles)_)        :: Pc
    INTEGER,                INTENT(INOUT) :: id
    INTEGER,                INTENT(  OUT) :: info

END SUBROUTINE

!DESTROY ENTRY
SUBROUTINE DTYPE(neigh_destroy)_(neigh,info)
    IMPORT DTYPE(ppm_t_neighlist)_
    CLASS(DTYPE(ppm_t_neighlist)_)      :: neigh
    INTEGER,             INTENT(  OUT)  :: info
END SUBROUTINE

FUNCTION DTYPE(has_ghosts)_(this,Field) RESULT(res)
    IMPORT DTYPE(ppm_t_particles)_,ppm_t_field_
    CLASS(DTYPE(ppm_t_particles)_)                 :: this
    CLASS(ppm_t_field_),OPTIONAL                   :: Field
    LOGICAL                                        :: res
END FUNCTION

SUBROUTINE DTYPE(part_get_discr)_(this,Field,prop,info)
    IMPORT DTYPE(ppm_t_particles)_,ppm_t_field_,DTYPE(ppm_t_part_prop)_,ppm_t_discr_data
    CLASS(DTYPE(ppm_t_particles)_)                         :: this
    CLASS(ppm_t_field_),TARGET,             INTENT(IN   )  :: Field
    !CLASS(DTYPE(ppm_t_part_prop)_),POINTER, INTENT(  OUT)  :: prop
    CLASS(ppm_t_discr_data),POINTER,        INTENT(  OUT)  :: prop
    INTEGER,                                INTENT(  OUT)  :: info
END SUBROUTINE
SUBROUTINE DTYPE(part_prop_zero)_(this,Field,info)
    IMPORT DTYPE(ppm_t_particles)_,ppm_t_field_
    CLASS(DTYPE(ppm_t_particles)_)                         :: this
    CLASS(ppm_t_field_),TARGET,             INTENT(IN   )  :: Field
    INTEGER,                                INTENT(  OUT)  :: info
END SUBROUTINE
SUBROUTINE DTYPE(part_p2m)_(this,Mesh,Field,kernel,info,p2m_bcdef)
    IMPORT DTYPE(ppm_t_particles)_,ppm_t_field_,ppm_t_equi_mesh_
    CLASS(DTYPE(ppm_t_particles)_)                  :: this
    CLASS(ppm_t_equi_mesh_)                         :: Mesh
    CLASS(ppm_t_field_)                             :: Field
    INTEGER                     ,     INTENT(IN   ) :: kernel
    INTEGER                     ,     INTENT(  OUT) :: info
    INTEGER, DIMENSION(:  )     , POINTER, OPTIONAL :: p2m_bcdef
END SUBROUTINE
SUBROUTINE DTYPE(part_remesh)_(this,Mesh,kernel,info,p2m_bcdef,cutoff_val,cutoff_field)
    IMPORT DTYPE(ppm_t_particles)_,ppm_t_field_,ppm_t_equi_mesh_, MK
    CLASS(DTYPE(ppm_t_particles)_)                  :: this
    CLASS(ppm_t_equi_mesh_)                         :: Mesh
    INTEGER                     ,     INTENT(IN   ) :: kernel
    INTEGER                     ,     INTENT(  OUT) :: info
    INTEGER, DIMENSION(:  )     , POINTER, OPTIONAL :: p2m_bcdef
    REAL(MK), DIMENSION(2),  OPTIONAL, INTENT(IN)   :: cutoff_val
    CLASS(ppm_t_field_),    OPTIONAL, INTENT(IN)    :: cutoff_field
END SUBROUTINE
