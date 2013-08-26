      !minclude define_abstract_collection_interfaces(DTYPE(ppm_t_dcop)_)
      
      !!CREATE
      !SUBROUTINE DTYPE(dcop_create)_(op,nterms,coeffs,degree,order,&
              !name,with_ghosts,vector,interp,Part,nlid,info)
          !IMPORT DTYPE(ppm_t_dcop)_,MK,ppm_t_main_abstr
          !CLASS(DTYPE(ppm_t_dcop)_)              :: op
          !INTEGER,                INTENT(IN   ) :: nterms
          !REAL(MK),DIMENSION(:),  INTENT(IN   ) :: coeffs
          !INTEGER,DIMENSION(:),   INTENT(IN   ) :: degree
          !INTEGER,DIMENSION(:),   INTENT(IN   ) :: order
          !LOGICAL,                INTENT(IN   ) :: with_ghosts
          !LOGICAL,                INTENT(IN   ) :: vector
          !LOGICAL,                INTENT(IN   ) :: interp
          !CLASS(ppm_t_main_abstr),POINTER, INTENT(IN   ) :: Part
          !INTEGER,                INTENT(IN   ) :: nlid
          !CHARACTER(LEN=*)                      :: name
          !INTEGER,                INTENT(OUT)   :: info
      !END SUBROUTINE
      !!DESTROY
      !SUBROUTINE DTYPE(dcop_destroy)_(this,info)
          !IMPORT DTYPE(ppm_t_dcop)_    
          !CLASS(DTYPE(ppm_t_dcop)_)          :: this
          !INTEGER,               INTENT(OUT) :: info
      !END SUBROUTINE
      !!INIT
      
      !!COMPUTE OPERATOR ON A DISCRETIZATION (PARTICLE SET OR MESH)
      !SUBROUTINE DTYPE(dcop_compute)_(this,Field_src,Field_to,info)
          !IMPORT DTYPE(ppm_t_dcop)_,ppm_t_main_abstr    
          !CLASS(DTYPE(ppm_t_dcop)_)                    :: this
          !CLASS(ppm_t_main_abstr),TARGET,INTENT(IN)    :: Field_src
          !CLASS(ppm_t_main_abstr),TARGET,INTENT(INOUT) :: Field_to
          !INTEGER,                       INTENT(OUT)   :: info
      !END SUBROUTINE

