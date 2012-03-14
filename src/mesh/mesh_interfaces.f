SUBROUTINE subpatch_get_field_2d_rd_(this,wp,Field,info)
    !!! gets pointer to data
    IMPORT ppm_t_subpatch_,ppm_t_field_,ppm_kind_double
    CLASS(ppm_t_subpatch_)              :: this
    CLASS(ppm_t_field_)                  :: Field
    REAL(ppm_kind_double),DIMENSION(:,:),POINTER :: wp
END SUBROUTINE

SUBROUTINE subpatch_get_field_3d_rd_(this,wp,Field,info)
    !!! gets pointer to data
    IMPORT ppm_t_subpatch_,ppm_t_field_,ppm_kind_double
    CLASS(ppm_t_subpatch_)              :: this
    CLASS(ppm_t_field_)                  :: Field
    REAL(ppm_kind_double),DIMENSION(:,:,:),POINTER :: wp
END SUBROUTINE

SUBROUTINE subpatch_data_create_(pdata,datatype,lda,Nmp,info)
    !!! Constructor for subdomain data 
    IMPORT ppm_t_subpatch_data_
    CLASS(ppm_t_subpatch_data_)              :: pdata
    INTEGER,                     INTENT(IN) :: datatype
    INTEGER,                     INTENT(IN) :: lda
    !!! number of data components per mesh node
    INTEGER,DIMENSION(:),POINTER,INTENT(IN) :: Nmp
    !!! number of mesh nodes in each dimension on this patch
    INTEGER,                    INTENT(OUT) :: info
END SUBROUTINE
!DESTROY
SUBROUTINE subpatch_data_destroy_(pdata,info)
    !!! Destructor for subdomain data
    IMPORT ppm_t_subpatch_data_
    CLASS(ppm_t_subpatch_data_)     :: pdata
    INTEGER,               INTENT(OUT) :: info
END SUBROUTINE

!CREATE
SUBROUTINE subpatch_create_(p,meshid,istart,iend,info)
    !!! Constructor for subpatch
    IMPORT ppm_t_subpatch_,ppm_kind_double
    CLASS(ppm_t_subpatch_)              :: p
    INTEGER                            :: meshid
    INTEGER,DIMENSION(:)               :: istart
    INTEGER,DIMENSION(:)               :: iend
    INTEGER,               INTENT(OUT) :: info
END SUBROUTINE

!DESTROY
SUBROUTINE subpatch_destroy_(p,info)
    !!! Destructor for subpatch
    IMPORT ppm_t_subpatch_
    CLASS(ppm_t_subpatch_)              :: p
    INTEGER,               INTENT(OUT) :: info
END SUBROUTINE


SUBROUTINE equi_mesh_create_(this,topoid,Nm,Offset,info)
    IMPORT ppm_t_equi_mesh_,ppm_kind_double
    CLASS(ppm_t_equi_mesh_)                  :: this
    !!! cartesian mesh object
    INTEGER                 , INTENT(IN   ) :: topoid
    !!! Topology ID for which mesh has been created 
    INTEGER , DIMENSION(:),   INTENT(IN   ) :: Nm
    !!! Global number of mesh points in the whole comput. domain
    REAL(ppm_kind_double), DIMENSION(:), INTENT(IN   ) :: Offset
    !!! Offset in each dimension
    INTEGER                 , INTENT(  OUT) :: info
    !!! Returns status, 0 upon success
END SUBROUTINE

SUBROUTINE equi_mesh_destroy_(this,info)
    IMPORT ppm_t_equi_mesh_
    CLASS(ppm_t_equi_mesh_)                  :: this
    INTEGER                 , INTENT(  OUT) :: info
END SUBROUTINE

SUBROUTINE equi_mesh_add_patch_(this,patch,info,patchid)
    !!! Add a patch to a mesh
    IMPORT ppm_t_equi_mesh_,ppm_kind_double
    CLASS(ppm_t_equi_mesh_)                  :: this
    REAL(ppm_kind_double),DIMENSION(:)      :: patch
    !!! Positions of the corners of the patch
    !!! (x1,y1,z1,x2,y2,z2), where 1 is the lower-left-bottom corner
    !!! and 2 is the upper-right-top corner.
    INTEGER                 , INTENT(  OUT) :: info
    !!! Returns status, 0 upon success
    INTEGER, OPTIONAL                       :: patchid
    !!! id of the patch, if we want one.
END SUBROUTINE
