SUBROUTINE subpatch_get_field_2d_rd_(this,wp,Field,info)
    !!! gets pointer to data
    IMPORT ppm_t_subpatch_,ppm_t_field_,ppm_kind_double
    CLASS(ppm_t_subpatch_)              :: this
    CLASS(ppm_t_field_)                  :: Field
    REAL(ppm_kind_double),DIMENSION(:,:),POINTER :: wp
    INTEGER,                 INTENT(OUT) :: info
END SUBROUTINE

SUBROUTINE subpatch_get_field_3d_rd_(this,wp,Field,info)
    !!! gets pointer to data
    IMPORT ppm_t_subpatch_,ppm_t_field_,ppm_kind_double
    CLASS(ppm_t_subpatch_)              :: this
    CLASS(ppm_t_field_)                  :: Field
    REAL(ppm_kind_double),DIMENSION(:,:,:),POINTER :: wp
    INTEGER,                 INTENT(OUT) :: info
END SUBROUTINE

SUBROUTINE subpatch_data_create_(this,discr_data,Nmp,info)
    !!! Constructor for subdomain data 
    IMPORT ppm_t_subpatch_data_,ppm_t_mesh_discr_data_
    CLASS(ppm_t_subpatch_data_)             :: this
    CLASS(ppm_t_mesh_discr_data_),TARGET,  INTENT(IN) :: discr_data
    INTEGER,DIMENSION(:),POINTER,INTENT(IN) :: Nmp
    INTEGER,                    INTENT(OUT) :: info
END SUBROUTINE
!DESTROY
SUBROUTINE subpatch_data_destroy_(this,info)
    !!! Destructor for subdomain data
    IMPORT ppm_t_subpatch_data_
    CLASS(ppm_t_subpatch_data_)        :: this
    INTEGER,               INTENT(OUT) :: info
END SUBROUTINE

!CREATE
SUBROUTINE subpatch_create_(p,meshID,istart,iend,istart_g,iend_g,info)
    !!! Constructor for subpatch
    IMPORT ppm_t_subpatch_,ppm_kind_double
    CLASS(ppm_t_subpatch_)             :: p
    INTEGER                            :: meshID
    INTEGER,DIMENSION(:)               :: istart
    INTEGER,DIMENSION(:)               :: iend
    INTEGER,DIMENSION(:)               :: istart_g
    INTEGER,DIMENSION(:)               :: iend_g
    INTEGER,               INTENT(OUT) :: info
END SUBROUTINE

!DESTROY
SUBROUTINE subpatch_destroy_(p,info)
    !!! Destructor for subpatch
    IMPORT ppm_t_subpatch_
    CLASS(ppm_t_subpatch_)              :: p
    INTEGER,               INTENT(OUT) :: info
END SUBROUTINE

!CREATE
SUBROUTINE subpatch_A_create_(this,vecsize,info,patchid)
    !!! Destructor for subdomain data data structure
    IMPORT ppm_t_A_subpatch_
    CLASS(ppm_t_A_subpatch_)           :: this
    INTEGER                            :: vecsize
    INTEGER,               INTENT(OUT) :: info
    INTEGER,OPTIONAL,      INTENT(IN ) :: patchid
END SUBROUTINE

!DESTROY
SUBROUTINE subpatch_A_destroy_(this,info)
    !!! Destructor for subdomain data data structure
    IMPORT ppm_t_A_subpatch_
    CLASS(ppm_t_A_subpatch_)           :: this
    INTEGER,               INTENT(OUT) :: info
END SUBROUTINE

SUBROUTINE equi_mesh_create_(this,topoid,Offset,info,Nm,h,ghostsize)
    IMPORT ppm_t_equi_mesh_,ppm_kind_double
    CLASS(ppm_t_equi_mesh_)                 :: this
    !!! cartesian mesh object
    INTEGER                 , INTENT(IN   ) :: topoid
    !!! Topology ID for which mesh has been created 
    REAL(ppm_kind_double), DIMENSION(:), INTENT(IN   ) :: Offset
    !!! Offset in each dimension
    INTEGER                 , INTENT(  OUT) :: info
    !!! Returns status, 0 upon success
    INTEGER,DIMENSION(:),              OPTIONAL,INTENT(IN   ) :: Nm
    !!! Global number of mesh points in the whole comput. domain
    !!! Makes sense only if the computational domain is bounded.
    !!! Note: Exactly one of Nm and h should be specified
    REAL(ppm_kind_double),DIMENSION(:),OPTIONAL,INTENT(IN   ) :: h
    !!! Mesh spacing
    !!! Note: Exactly one of Nm and h should be specified
    INTEGER,DIMENSION(:),              OPTIONAL,INTENT(IN   ) :: ghostsize
END SUBROUTINE

SUBROUTINE equi_mesh_destroy_(this,info)
    IMPORT ppm_t_equi_mesh_
    CLASS(ppm_t_equi_mesh_)                 :: this
    INTEGER                 , INTENT(  OUT) :: info
END SUBROUTINE

SUBROUTINE equi_mesh_def_patch_(this,patch,info,patchid,infinite)
    !!! Add a patch to a mesh
    IMPORT ppm_t_equi_mesh_,ppm_kind_double
    CLASS(ppm_t_equi_mesh_)                 :: this
    REAL(ppm_kind_double),DIMENSION(:)      :: patch
    !!! Positions of the corners of the patch
    !!! (x1,y1,z1,x2,y2,z2), where 1 is the lower-left-bottom corner
    !!! and 2 is the upper-right-top corner.
    INTEGER                 , INTENT(  OUT) :: info
    !!! Returns status, 0 upon success
    INTEGER, OPTIONAL                       :: patchid
    !!! id of the patch, if we want one.
    LOGICAL, OPTIONAL                       :: infinite
    !!! true if the patch should cover the whole computational domain
END SUBROUTINE

SUBROUTINE equi_mesh_def_uniform_(this,info,patchid)
    !!! Add a uniform patch to a mesh
    IMPORT ppm_t_equi_mesh_,ppm_kind_double
    CLASS(ppm_t_equi_mesh_)                 :: this
    INTEGER                 , INTENT(  OUT) :: info
    !!! Returns status, 0 upon success
    INTEGER, OPTIONAL                       :: patchid
    !!! id of the (uniform) patch, if we want one.
END SUBROUTINE

FUNCTION equi_mesh_new_subpatch_data_ptr_(this,info) RESULT(sp)
    !!! returns a pointer to a new subpatch_data object
    IMPORT ppm_t_equi_mesh_,ppm_t_subpatch_data_
    IMPLICIT NONE
    CLASS(ppm_t_equi_mesh_)                 :: this
    !!! cartesian mesh object
    CLASS(ppm_t_subpatch_data_),POINTER     :: sp
    INTEGER                 , INTENT(  OUT) :: info
    !!! Returns status, 0 upon success
END FUNCTION 
!CREATE
SUBROUTINE field_info_create_(this,field,info)
    IMPORT ppm_t_field_info_,ppm_t_main_abstr
    CLASS(ppm_t_field_info_)              :: this
    CLASS(ppm_t_main_abstr),TARGET,INTENT(IN ) :: field
    INTEGER,                  INTENT(OUT) :: info
END SUBROUTINE
!DESTROY
SUBROUTINE field_info_destroy_(this,info)
    IMPORT ppm_t_field_info_
    CLASS(ppm_t_field_info_)           :: this
    INTEGER,               INTENT(OUT) :: info
END SUBROUTINE 
!LIST OF DISCRETIZED FIELDS
FUNCTION equi_mesh_list_of_fields_(this,info) RESULT(fids)
    IMPORT ppm_t_equi_mesh_
    CLASS(ppm_t_equi_mesh_)         :: this
    INTEGER,DIMENSION(:),POINTER    :: fids
    INTEGER                         :: info
END FUNCTION
!ESTABLISH RELATIONSHIP BETWEEN MESH AND FIELD
SUBROUTINE equi_mesh_set_rel_(this,field,info)
    IMPORT ppm_t_field_,ppm_t_equi_mesh_
    CLASS(ppm_t_equi_mesh_)            :: this
    CLASS(ppm_t_field_)                :: field
    !!! this mesh is discretized on that field
    INTEGER,               INTENT(OUT)  :: info
END SUBROUTINE
!GHOST GET
SUBROUTINE equi_mesh_map_ghost_get_(this,info,ghostsize)
    IMPORT ppm_t_equi_mesh_
    CLASS(ppm_t_equi_mesh_)                 :: this
    INTEGER                 , INTENT(  OUT) :: info
    INTEGER, OPTIONAL, DIMENSION(:)   , INTENT(IN   ) :: ghostsize
END SUBROUTINE
SUBROUTINE equi_mesh_block_intersect_(this,to_mesh,isub,jsub,offset,&
        ghostsize,nsendlist,isendfromsub,isendtosub,isendpatchid,&
        isendblkstart,isendblksize,ioffset,info,lsymm)
    IMPORT ppm_t_equi_mesh_
    CLASS(ppm_t_equi_mesh_)                 :: this
    CLASS(ppm_t_equi_mesh_)                 :: to_mesh
    INTEGER                 , INTENT(IN   ) :: isub
    INTEGER                 , INTENT(IN   ) :: jsub
    INTEGER, DIMENSION(:)   , INTENT(IN   ) :: ghostsize
    INTEGER, DIMENSION(:)   , INTENT(IN   ) :: offset
    INTEGER, DIMENSION(:)   , POINTER       :: isendfromsub
    INTEGER, DIMENSION(:)   , POINTER       :: isendtosub
    INTEGER, DIMENSION(:,:) , POINTER       :: isendpatchid
    INTEGER, DIMENSION(:,:) , POINTER       :: ioffset
    INTEGER, DIMENSION(:,:) , POINTER       :: isendblkstart
    INTEGER, DIMENSION(:,:) , POINTER       :: isendblksize
    INTEGER                 , INTENT(INOUT) :: nsendlist
    INTEGER                 , INTENT(  OUT) :: info
    LOGICAL, DIMENSION(3)   , INTENT(IN   ), OPTIONAL :: lsymm
END SUBROUTINE
SUBROUTINE equi_mesh_map_ghost_init_(this,ghostsize,info)
    IMPORT ppm_t_equi_mesh_
    CLASS(ppm_t_equi_mesh_)                 :: this
    INTEGER, DIMENSION(:)   , INTENT(IN   ) :: ghostsize
    INTEGER                 , INTENT(  OUT) :: info
END SUBROUTINE
SUBROUTINE equi_mesh_map_ghost_push_(this,field,info)
    IMPORT ppm_t_equi_mesh_,ppm_t_field_
    CLASS(ppm_t_equi_mesh_)            :: this
    CLASS(ppm_t_field_)                :: field
    INTEGER,               INTENT(OUT) :: info
END SUBROUTINE 
SUBROUTINE equi_mesh_map_send_(this,info)
    IMPORT ppm_t_equi_mesh_
    CLASS(ppm_t_equi_mesh_)            :: this
    INTEGER,               INTENT(OUT) :: info
END SUBROUTINE 
SUBROUTINE equi_mesh_map_ghost_pop_(this,field,info)
    IMPORT ppm_t_equi_mesh_,ppm_t_field_
    CLASS(ppm_t_equi_mesh_)            :: this
    CLASS(ppm_t_field_)                :: field
    INTEGER,               INTENT(OUT) :: info
END SUBROUTINE 
SUBROUTINE mesh_discr_data_create_(this,field,info)
    IMPORT ppm_t_mesh_discr_data_,ppm_t_field_
    CLASS(ppm_t_mesh_discr_data_)       :: this
    CLASS(ppm_t_field_),    INTENT(IN)  :: field
    INTEGER,                INTENT(OUT) :: info
END SUBROUTINE
SUBROUTINE mesh_discr_data_destroy_(this,info)
    IMPORT ppm_t_mesh_discr_data_
    CLASS(ppm_t_mesh_discr_data_)       :: this
    INTEGER,                INTENT(OUT) :: info
END SUBROUTINE
SUBROUTINE equi_mesh_create_prop_(this,field,info,discr_data)
    IMPORT ppm_t_equi_mesh_,ppm_t_mesh_discr_data_,ppm_t_field_
    CLASS(ppm_t_equi_mesh_)                           :: this
    CLASS(ppm_t_field_),                  INTENT(IN)  :: field
    INTEGER,                              INTENT(OUT) :: info
    CLASS(ppm_t_mesh_discr_data_),POINTER,INTENT(OUT) :: discr_data
END SUBROUTINE
