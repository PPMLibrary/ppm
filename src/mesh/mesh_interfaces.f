minclude ppm_get_field_interface_template(2,s)
minclude ppm_get_field_interface_template(3,s)
minclude ppm_get_field_interface_template(4,s)
minclude ppm_get_field_interface_template(2,d)
minclude ppm_get_field_interface_template(3,d)
minclude ppm_get_field_interface_template(4,d)

SUBROUTINE subpatch_data_create_(this,discr_data,sp,info)
    !!! Constructor for subdomain data 
    IMPORT ppm_t_subpatch_data_,ppm_t_mesh_discr_data_,ppm_t_subpatch_
    CLASS(ppm_t_subpatch_data_)             :: this
    CLASS(ppm_t_mesh_discr_data_),TARGET,  INTENT(IN) :: discr_data
    CLASS(ppm_t_subpatch_),      INTENT(IN) :: sp
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
SUBROUTINE subpatch_create_(p,mesh,isub,istart,iend,pstart,pend,&
        istart_p,iend_p,ghostsize,bcdef,info)
    !!! Constructor for subpatch
    IMPORT ppm_t_subpatch_,ppm_t_equi_mesh_,ppm_kind_double
    CLASS(ppm_t_subpatch_)             :: p
    CLASS(ppm_t_equi_mesh_),TARGET     :: mesh
    INTEGER                            :: isub
    INTEGER,DIMENSION(:)               :: istart
    INTEGER,DIMENSION(:)               :: iend
    REAL(ppm_kind_double),DIMENSION(:) :: pstart
    REAL(ppm_kind_double),DIMENSION(:) :: pend
    INTEGER,DIMENSION(:)               :: istart_p
    INTEGER,DIMENSION(:)               :: iend_p
    INTEGER,DIMENSION(:)               :: ghostsize
    INTEGER,DIMENSION(:)               :: bcdef
    INTEGER,               INTENT(OUT) :: info
END SUBROUTINE

!DESTROY
SUBROUTINE subpatch_destroy_(p,info)
    !!! Destructor for subpatch
    IMPORT ppm_t_subpatch_
    CLASS(ppm_t_subpatch_)              :: p
    INTEGER,               INTENT(OUT) :: info
END SUBROUTINE

!GET POSITION
PURE FUNCTION subpatch_get_pos2d_(p,i,j) RESULT (pos)
    IMPORT ppm_t_subpatch_,ppm_kind_double,ppm_dim
    CLASS(ppm_t_subpatch_), INTENT(IN) :: p
    INTEGER,                INTENT(IN) :: i
    INTEGER,                INTENT(IN) :: j
    REAL(ppm_kind_double),DIMENSION(ppm_dim) :: pos
END FUNCTION
PURE FUNCTION subpatch_get_pos3d_(p,i,j,k) RESULT (pos)
    IMPORT ppm_t_subpatch_,ppm_kind_double,ppm_dim
    CLASS(ppm_t_subpatch_), INTENT(IN) :: p
    INTEGER,                INTENT(IN) :: i
    INTEGER,                INTENT(IN) :: j
    INTEGER,                INTENT(IN) :: k
    REAL(ppm_kind_double),DIMENSION(ppm_dim) :: pos
END FUNCTION

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

SUBROUTINE equi_mesh_create_(this,topoid,Offset,info,Nm,h,ghostsize,name)
    IMPORT ppm_t_equi_mesh_,ppm_kind_double
    CLASS(ppm_t_equi_mesh_)                 :: this
    INTEGER                 , INTENT(IN   ) :: topoid
    REAL(ppm_kind_double), DIMENSION(:), INTENT(IN   ) :: Offset
    INTEGER                 , INTENT(  OUT) :: info
    INTEGER,DIMENSION(:),              OPTIONAL,INTENT(IN   ) :: Nm
    REAL(ppm_kind_double),DIMENSION(:),OPTIONAL,INTENT(IN   ) :: h
    INTEGER,DIMENSION(:),              OPTIONAL,INTENT(IN   ) :: ghostsize
    CHARACTER(LEN=*), OPTIONAL,                 INTENT(IN   ) :: name
    !!! name of this mesh
END SUBROUTINE

SUBROUTINE equi_mesh_destroy_(this,info)
    IMPORT ppm_t_equi_mesh_
    CLASS(ppm_t_equi_mesh_)                 :: this
    INTEGER                 , INTENT(  OUT) :: info
END SUBROUTINE

SUBROUTINE equi_mesh_def_patch_(this,patch,info,patchid,infinite,bcdef)
    !!! Add a patch to a mesh
    IMPORT ppm_t_equi_mesh_,ppm_kind_double,ppm_dim
    CLASS(ppm_t_equi_mesh_)                 :: this
    REAL(ppm_kind_double),DIMENSION(:)      :: patch
    INTEGER                 , INTENT(  OUT) :: info
    INTEGER, OPTIONAL                       :: patchid
    LOGICAL, OPTIONAL                       :: infinite
    INTEGER, OPTIONAL,    DIMENSION(2*ppm_dim)  :: bcdef
END SUBROUTINE

SUBROUTINE equi_mesh_def_uniform_(this,info,patchid)
    !!! Add a uniform patch to a mesh
    IMPORT ppm_t_equi_mesh_,ppm_kind_double
    CLASS(ppm_t_equi_mesh_)                 :: this
    INTEGER                 , INTENT(  OUT) :: info
    INTEGER, OPTIONAL                       :: patchid
END SUBROUTINE

FUNCTION equi_mesh_new_subpatch_data_ptr_(this,info) RESULT(sp)
    !!! returns a pointer to a new subpatch_data object
    IMPORT ppm_t_equi_mesh_,ppm_t_subpatch_data_
    IMPLICIT NONE
    CLASS(ppm_t_equi_mesh_)                 :: this
    CLASS(ppm_t_subpatch_data_),POINTER     :: sp
    INTEGER                 , INTENT(  OUT) :: info
END FUNCTION 
!LIST OF DISCRETIZED FIELDS
FUNCTION equi_mesh_list_of_fields_(this,info) RESULT(fids)
    IMPORT ppm_t_equi_mesh_
    CLASS(ppm_t_equi_mesh_)         :: this
    INTEGER,DIMENSION(:),POINTER    :: fids
    INTEGER                         :: info
END FUNCTION
!GHOST GET
SUBROUTINE equi_mesh_map_ghost_get_(this,info)
    IMPORT ppm_t_equi_mesh_
    CLASS(ppm_t_equi_mesh_)                 :: this
    INTEGER                 , INTENT(  OUT) :: info
END SUBROUTINE
!GHOST PUT
SUBROUTINE equi_mesh_map_ghost_put_(this,info)
    IMPORT ppm_t_equi_mesh_
    CLASS(ppm_t_equi_mesh_)                 :: this
    INTEGER                 , INTENT(  OUT) :: info
END SUBROUTINE
SUBROUTINE equi_mesh_block_intersect_(this,to_mesh,isub_loc,isub,jsub,offset,&
        nsendlist,isendfromsub,isendtosub,isendpatchid,&
        isendblkstart,isendblksize,ioffset,info,lsymm)
    IMPORT ppm_t_equi_mesh_
    CLASS(ppm_t_equi_mesh_)                 :: this
    CLASS(ppm_t_equi_mesh_)                 :: to_mesh
    INTEGER                 , INTENT(IN   ) :: isub_loc
    INTEGER                 , INTENT(IN   ) :: isub
    INTEGER                 , INTENT(IN   ) :: jsub
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
SUBROUTINE equi_mesh_map_ghost_init_(this,info)
    IMPORT ppm_t_equi_mesh_
    CLASS(ppm_t_equi_mesh_)                 :: this
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
    CLASS(ppm_t_mesh_discr_data_)          :: this
    CLASS(ppm_t_field_),TARGET,INTENT(IN)  :: field
    INTEGER,                   INTENT(OUT) :: info
END SUBROUTINE
SUBROUTINE mesh_discr_data_destroy_(this,info)
    IMPORT ppm_t_mesh_discr_data_
    CLASS(ppm_t_mesh_discr_data_)       :: this
    INTEGER,                INTENT(OUT) :: info
END SUBROUTINE
SUBROUTINE equi_mesh_create_prop_(this,field,discr_data,info,p_idx)
    IMPORT ppm_t_equi_mesh_,ppm_t_mesh_discr_data_,ppm_t_field_
    CLASS(ppm_t_equi_mesh_)                           :: this
    CLASS(ppm_t_field_),                  INTENT(IN)  :: field
    CLASS(ppm_t_mesh_discr_data_),POINTER,INTENT(OUT) :: discr_data
    INTEGER,                              INTENT(OUT) :: info
    INTEGER, OPTIONAL,                    INTENT(OUT) :: p_idx
END SUBROUTINE
SUBROUTINE equi_mesh_print_vtk_(this,filename,info)
    IMPORT ppm_t_equi_mesh_
    CLASS(ppm_t_equi_mesh_)                           :: this
    CHARACTER(LEN=*)                                  :: filename
    INTEGER,                              INTENT(OUT) :: info
END SUBROUTINE
SUBROUTINE equi_mesh_m2p_(this,Part,Field,kernel,info)
    IMPORT ppm_t_equi_mesh_,ppm_t_particles_d_,ppm_t_field_
    CLASS(ppm_t_equi_mesh_)                          :: this
    CLASS(ppm_t_particles_d_)                        :: Part
    CLASS(ppm_t_field_)                              :: Field
    INTEGER                        , INTENT(IN   ) :: kernel
    INTEGER                        , INTENT(  OUT) :: info
END SUBROUTINE
SUBROUTINE equi_mesh_prop_zero_(this,Field,info)
    IMPORT ppm_t_equi_mesh_,ppm_t_field_
    CLASS(ppm_t_equi_mesh_)                           :: this
    CLASS(ppm_t_field_)                               :: Field
    INTEGER,                              INTENT(OUT) :: info
END SUBROUTINE
