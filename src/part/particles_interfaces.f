minclude define_abstract_collection_interfaces(DTYPE(ppm_t_part_prop)_)
minclude define_abstract_collection_interfaces(DTYPE(ppm_t_neighlist)_)
minclude define_abstract_collection_interfaces(DTYPE(ppm_t_particles)_)
!minclude define_abstract_collection_interfaces(DTYPE(ppm_t_sop)_)

!CREATE ENTRY
SUBROUTINE DTYPE(prop_create)_(prop,datatype,npart,lda,name,flags,info,field,zero)
    !!! Constructor for particle property data structure
    IMPORT DTYPE(ppm_t_part_prop)_,ppm_param_length_pptflags,ppm_t_field_
    CLASS(DTYPE(ppm_t_part_prop)_)      :: prop
    INTEGER,                INTENT(IN) :: datatype
    INTEGER,                INTENT(IN) :: npart
    INTEGER,                INTENT(IN) :: lda
    CHARACTER(LEN=*),       INTENT(IN) :: name
    LOGICAL, DIMENSION(ppm_param_length_pptflags),INTENT(IN) :: flags
    INTEGER,               INTENT(OUT) :: info
    CLASS(ppm_t_field_),OPTIONAL,TARGET, INTENT(IN) :: field
    LOGICAL, OPTIONAL,     INTENT( IN) :: zero
    !!! if true, then initialize the data to zero
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
    !!! Data structure containing the particles
    INTEGER,                            INTENT(  OUT)      :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    !  Optional arguments
    !-------------------------------------------------------------------------
    INTEGER, OPTIONAL,                  INTENT(IN   )      :: level
    !!! indentation level at which to printout the info. Default = 0
    INTEGER, OPTIONAL,                  INTENT(IN   )      :: fileunit
    !!! Already open file unit for printout. Default = stdout
    INTEGER, OPTIONAL,                  INTENT(IN   )      :: propid
    !!! id of this property in the parent struct
END SUBROUTINE

!!----------------------------------------------------------------
!! Procedures for Particle Sets DS
!!----------------------------------------------------------------
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

SUBROUTINE DTYPE(part_prop_create)_(this,info,field,discr_data,&
        dtype,name,lda,zero,with_ghosts)
    IMPORT DTYPE(ppm_t_particles)_,ppm_t_field_,DTYPE(ppm_t_part_prop)_
    !!! Adds a property to an existing particle set
    CLASS(DTYPE(ppm_t_particles)_)        :: this
    INTEGER,               INTENT(OUT)    :: info
    CLASS(ppm_t_field_),OPTIONAL,INTENT(IN   ) :: field
    CLASS(DTYPE(ppm_t_part_prop)_),OPTIONAL,POINTER,INTENT(OUT):: discr_data
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

SUBROUTINE DTYPE(part_prop_realloc)_(Pc,id,info,with_ghosts,datatype,lda)
    IMPORT DTYPE(ppm_t_particles)_
    !!! Reallocate the property array to the correct size
    !!! (e.g. if the number of particles has changed or if the type
    !!! of the data changes)
    CLASS(DTYPE(ppm_t_particles)_)         :: Pc
    INTEGER,                INTENT(IN   ) :: id
    INTEGER,               INTENT(OUT)    :: info
    LOGICAL, OPTIONAL                     :: with_ghosts
    !!! if true, then allocate with Mpart instead of the default size of Npart
    INTEGER, OPTIONAL                     :: datatype
    !!! deallocate the old data array and allocate a new one,
    !!! possibly of a different data type.
    INTEGER, OPTIONAL                     :: lda
    !!! deallocate the old data array and allocate a new one,
    !!! possibly of a different dimension
END SUBROUTINE


SUBROUTINE DTYPE(part_neigh_create)_(Pc,id,info,&
        P_id,name,skin,symmetry,cutoff)
    IMPORT DTYPE(ppm_t_particles)_, MK
    !!! Create a data structure to store a neighbour list
    CLASS(DTYPE(ppm_t_particles)_)         :: Pc
    INTEGER,               INTENT(  OUT)  :: id
    INTEGER,               INTENT(OUT)    :: info
    INTEGER, OPTIONAL                     :: P_id    
    !!! Id of the set of particles that this neighbor list refers to
    !!! The default, 0, stands for "self"
    CHARACTER(LEN=*) , OPTIONAL           :: name
    !!! name of this neighbour list
    REAL(MK), OPTIONAL                    :: skin
    REAL(MK), OPTIONAL                    :: cutoff
    LOGICAL, OPTIONAL                     :: symmetry    
END SUBROUTINE

SUBROUTINE DTYPE(part_neigh_destroy)_(Pc,id,info)
    IMPORT DTYPE(ppm_t_particles)_
    !!! Destroy a property from an existing particle set
    CLASS(DTYPE(ppm_t_particles)_)         :: Pc
    INTEGER,                INTENT(INOUT) :: id
    INTEGER,               INTENT(OUT)    :: info
END SUBROUTINE


SUBROUTINE DTYPE(part_create)_(Pc,Npart,info,name)
    IMPORT DTYPE(ppm_t_particles)_
    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    CLASS(DTYPE(ppm_t_particles)_)                          :: Pc
    !!! Data structure containing the particles
    INTEGER,                                INTENT(IN   )  :: Npart
    !!! Number of particles 
    INTEGER,                                INTENT(  OUT)  :: info
    !!! Returns status, 0 upon success.
    !-------------------------------------------------------------------------
    !  Optional arguments
    !-------------------------------------------------------------------------
    CHARACTER(LEN=*) , OPTIONAL                            :: name
    !!! give a name to this Particle set
END SUBROUTINE


SUBROUTINE DTYPE(part_destroy)_(Pc,info)
    IMPORT DTYPE(ppm_t_particles)_
    CLASS(DTYPE(ppm_t_particles)_)                          :: Pc
    !!! Data structure containing the Pc
    INTEGER,                                INTENT(  OUT)  :: info
    !!! Returns status, 0 upon success.
END SUBROUTINE

SUBROUTINE DTYPE(particles_initialize2d)_(Pc,Npart_global,info,&
        distrib,topoid,minphys,maxphys,cutoff,name)
    IMPORT DTYPE(ppm_t_particles)_,MK,ppm_dim
    CLASS(DTYPE(ppm_t_particles)_)                          :: Pc
    INTEGER,                             INTENT(INOUT)     :: Npart_global
    !!! total number of particles that will be initialized
    INTEGER,                             INTENT(  OUT)     :: info
    !!! Returns status, 0 upon success.
    INTEGER,OPTIONAL,                    INTENT(IN   )     :: distrib
    !!! type of initial distribution. One of
    !!! ppm_param_part_init_cartesian (default)
    !!! ppm_param_part_init_random
    INTEGER,OPTIONAL,                    INTENT(IN   )     :: topoid
    !!! topology id (used only to get the extent of the physical domain)
    REAL(MK),DIMENSION(ppm_dim),OPTIONAL,INTENT(IN   )     :: minphys
    !!! extent of the physical domain. Only if topoid is not present.
    REAL(MK),DIMENSION(ppm_dim),OPTIONAL,INTENT(IN   )     :: maxphys
    !!! extent of the physical domain. Only if topoid is not present.
    REAL(MK),                   OPTIONAL,INTENT(IN   )     :: cutoff
    !!! cutoff of the particles
    CHARACTER(LEN=*),           OPTIONAL,INTENT(IN   )     :: name
    !!! name for this set of particles
END SUBROUTINE

SUBROUTINE DTYPE(particles_initialize3d)_(Pc,Npart_global,info,&
        distrib,topoid,minphys,maxphys,cutoff,name)
    IMPORT DTYPE(ppm_t_particles)_,MK,ppm_dim
    CLASS(DTYPE(ppm_t_particles)_)                          :: Pc
    INTEGER,                             INTENT(INOUT)     :: Npart_global
    !!! total number of particles that will be initialized
    INTEGER,                             INTENT(  OUT)     :: info
    !!! Returns status, 0 upon success.
    INTEGER,OPTIONAL,                    INTENT(IN   )     :: distrib
    !!! type of initial distribution. One of
    !!! ppm_param_part_init_cartesian (default)
    !!! ppm_param_part_init_random
    INTEGER,OPTIONAL,                    INTENT(IN   )     :: topoid
    !!! topology id (used only to get the extent of the physical domain)
    REAL(MK),DIMENSION(ppm_dim),OPTIONAL,INTENT(IN   )     :: minphys
    !!! extent of the physical domain. Only if topoid is not present.
    REAL(MK),DIMENSION(ppm_dim),OPTIONAL,INTENT(IN   )     :: maxphys
    !!! extent of the physical domain. Only if topoid is not present.
    REAL(MK),                   OPTIONAL,INTENT(IN   )     :: cutoff
    !!! cutoff of the particles
    CHARACTER(LEN=*),           OPTIONAL,INTENT(IN   )     :: name
    !!! name for this set of particles
END SUBROUTINE

!!temporary hack to deal with both 2d and 3d
SUBROUTINE DTYPE(part_initialize)_(Pc,Npart_global,info,&
        distrib,topoid,minphys,maxphys,cutoff,name)
    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    IMPORT DTYPE(ppm_t_particles)_,MK,ppm_dim
    CLASS(DTYPE(ppm_t_particles)_)                          :: Pc
    !!! Data structure containing the particles
    INTEGER,                            INTENT(INOUT)      :: Npart_global
    !!! total number of particles that will be initialized
    INTEGER,                            INTENT(  OUT)      :: info
    !!! Return status, on success 0.
    INTEGER,OPTIONAL,                    INTENT(IN   )     :: distrib
    !!! type of initial distribution. One of
    !!! ppm_param_part_init_cartesian (default)
    !!! ppm_param_part_init_random
    INTEGER,OPTIONAL,                    INTENT(IN   )     :: topoid
    !!! topology id (used only to get the extent of the physical domain)
    REAL(MK),DIMENSION(ppm_dim),OPTIONAL,INTENT(IN   )     :: minphys
    !!! extent of the physical domain. Only if topoid is not present.
    REAL(MK),DIMENSION(ppm_dim),OPTIONAL,INTENT(IN   )     :: maxphys
    !!! extent of the physical domain. Only if topoid is not present.
    REAL(MK),                   OPTIONAL,INTENT(IN   )     :: cutoff
    !!! cutoff of the particles
    CHARACTER(LEN=*),           OPTIONAL,INTENT(IN   )     :: name
    !!! name for this set of particles
END SUBROUTINE

SUBROUTINE DTYPE(part_print_info)_(Pc,info,level,fileunit)
    IMPORT DTYPE(ppm_t_particles)_
    CLASS(DTYPE(ppm_t_particles)_)                          :: Pc
    !!! Data structure containing the particles
    INTEGER,                            INTENT(  OUT)      :: info
    !!! Return status, on success 0.
    INTEGER, OPTIONAL,                  INTENT(IN   )      :: level
    !!! indentation level at which to printout the info. Default = 0
    INTEGER, OPTIONAL,                  INTENT(IN   )      :: fileunit
    !!! Already open file unit for printout. Default = stdout
END SUBROUTINE

SUBROUTINE DTYPE(part_del_parts)_(Pc,list_del_parts,nb_del,info)
    IMPORT DTYPE(ppm_t_particles)_
    CLASS(DTYPE(ppm_t_particles)_)                         :: Pc
    !!! Data structure containing the particles
    INTEGER,DIMENSION(:),POINTER,           INTENT(IN   )  :: list_del_parts
    !!! list of particles to be deleted
    INTEGER,                                INTENT(IN   )  :: nb_del
    !!! number of particles to be deleted
    INTEGER,                                INTENT(  OUT)  :: info
    !!! Returns status, 0 upon success.
END SUBROUTINE

!ESTABLISH RELATIONSHIP BETWEEN A PARTICLE SET AND A FIELD
SUBROUTINE DTYPE(part_set_rel)_(this,field,info)
    IMPORT ppm_t_field_,DTYPE(ppm_t_particles)_
    CLASS(DTYPE(ppm_t_particles)_)     :: this
    CLASS(ppm_t_field_)                :: field
    INTEGER,               INTENT(OUT) :: info
END SUBROUTINE

SUBROUTINE DTYPE(part_prop_push)_(Pc,prop_id,info)
    IMPORT DTYPE(ppm_t_particles)_
    CLASS(DTYPE(ppm_t_particles)_)                          :: Pc
    !!! Data structure containing the particles
    INTEGER,                                INTENT(IN   )  :: prop_id
    !!! id of the property to be pushed
    INTEGER,                                INTENT(  OUT)  :: info
    !!! Returns status, 0 upon success.
END SUBROUTINE

SUBROUTINE DTYPE(part_prop_pop)_(Pc,prop_id,Npart_new,info)
    IMPORT DTYPE(ppm_t_particles)_
    CLASS(DTYPE(ppm_t_particles)_)                          :: Pc
    !!! Data structure containing the particles
    INTEGER,                                INTENT(IN   )  :: prop_id
    !!! id of the property to be pushed
    INTEGER,                                INTENT(IN   )  :: Npart_new
    !!! number of particles to pop from the buffer
    INTEGER,                                INTENT(  OUT)  :: info
    !!! Returns status, 0 upon success.
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

SUBROUTINE DTYPE(part_mapping)_(Pc,info,debug,global,topoid)
    IMPORT DTYPE(ppm_t_particles)_
    CLASS(DTYPE(ppm_t_particles)_)                          :: Pc
    !!! Data structure containing the particles
    INTEGER,                            INTENT(  OUT)      :: info
    !!! Return status, on success 0.
    LOGICAL, OPTIONAL                                   :: debug
    !!! IF true, printout more
    LOGICAL, OPTIONAL                                   :: global
    !!! does a global mapping. Default is false (i.e. partial mapping)
    INTEGER, OPTIONAL                                   :: topoid
    !!! topology id
END SUBROUTINE

SUBROUTINE DTYPE(part_mapping_ghosts)_(Pc,info,ghostsize,debug)
    IMPORT DTYPE(ppm_t_particles)_,MK
    CLASS(DTYPE(ppm_t_particles)_)                       :: Pc
    !!! Data structure containing the particles
    INTEGER,                            INTENT(  OUT)   :: info
    !!! Return status, on success 0.
    REAL(MK), OPTIONAL                                  :: ghostsize
    !!! size of the ghost layers. Default is to use the particles cutoff
    LOGICAL, OPTIONAL                                   :: debug
    !!! IF true, printout more
END SUBROUTINE

SUBROUTINE DTYPE(part_apply_bc)_(Pc,info)
    IMPORT DTYPE(ppm_t_particles)_
    CLASS(DTYPE(ppm_t_particles)_)                         :: Pc
    !!! Data structure containing the particles
    INTEGER,                            INTENT(  OUT)     :: info
    !!! Return status, on success 0.
END SUBROUTINE

SUBROUTINE DTYPE(part_move)_(Pc,disp,info)
    IMPORT DTYPE(ppm_t_particles)_, MK
    CLASS(DTYPE(ppm_t_particles)_)                         :: Pc
    !!! Data structure containing the particles
    REAL(MK), DIMENSION(:,:), POINTER,  INTENT(IN   )     :: disp
    !!! Data structure containing the particles
    INTEGER,                            INTENT(  OUT)     :: info
    !!! Return status, on success 0.
END SUBROUTINE

SUBROUTINE DTYPE(part_neighlist)_(Pc,info,nlid,lstore,incl_ghosts,knn)
    IMPORT DTYPE(ppm_t_particles)_
    CLASS(DTYPE(ppm_t_particles)_)                          :: Pc
    !!! Data structure containing the particles
    INTEGER,                            INTENT(  OUT)      :: info
    !!! Return status, on success 0.
    INTEGER, OPTIONAL,                  INTENT(INOUT)      :: nlid
    !!! which neighbour list are we computing. Default is 1
    LOGICAL, OPTIONAL,                  INTENT(IN   )      :: lstore
    !!! store verlet lists
    LOGICAL, OPTIONAL,                  INTENT(IN   )      :: incl_ghosts
    !!! if true, then verlet lists are computed for all particles, incl. ghosts.
    !!! Default is false.
    INTEGER, OPTIONAL,                  INTENT(IN   )      :: knn
    !!! if present, neighbour lists are constructed such that each particle
    !!! has at least knn neighbours.
END SUBROUTINE

SUBROUTINE DTYPE(part_set_cutoff)_(Pc,cutoff,info,Nlist)
    IMPORT DTYPE(ppm_t_particles)_,MK,DTYPE(ppm_t_neighlist)_
    CLASS(DTYPE(ppm_t_particles)_)            :: Pc
    REAL(MK),                 INTENT(IN   )  :: cutoff
    !!! cutoff radius
    INTEGER,                  INTENT(   OUT) :: info
    !!! return status. On success, 0
    CLASS(DTYPE(ppm_t_neighlist)_),OPTIONAL,INTENT(INOUT) :: Nlist
    !!! ID of the neighbor list for which this cutoff radius
    !!! applies. Default is ppm_param_default_nlID
END SUBROUTINE

SUBROUTINE DTYPE(part_comp_global_index)_(Pc,info)
    IMPORT DTYPE(ppm_t_particles)_
    CLASS(DTYPE(ppm_t_particles)_)            :: Pc
    INTEGER,                  INTENT(   OUT) :: info
    !!! return status. On success, 0
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
    IMPORT DTYPE(ppm_t_particles)_,ppm_t_field_,DTYPE(ppm_t_part_prop)_
    CLASS(DTYPE(ppm_t_particles)_)                         :: this
    CLASS(ppm_t_field_),TARGET,             INTENT(IN   )  :: Field
    CLASS(DTYPE(ppm_t_part_prop)_),POINTER, INTENT(  OUT)  :: prop
    INTEGER,                                INTENT(  OUT)  :: info
END SUBROUTINE
