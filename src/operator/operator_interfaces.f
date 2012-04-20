!--------------------------------------
!GENERIC OPERATOR CLASS
!--------------------------------------

!CREATE
SUBROUTINE operator_create_(this,nterms,coeffs,degree,info,name)
    IMPORT ppm_t_operator_,ppm_t_main_abstr,ppm_kind_double
    CLASS(ppm_t_operator_)                :: this
    INTEGER,                INTENT(IN   ) :: nterms
    REAL(ppm_kind_double),DIMENSION(:),INTENT(IN):: coeffs
    INTEGER,DIMENSION(:),   INTENT(IN   ) :: degree
    INTEGER,                INTENT(OUT)   :: info
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)  :: name
END SUBROUTINE
!DESTROY
SUBROUTINE operator_destroy_(this,info)
    IMPORT ppm_t_operator_,ppm_t_main_abstr
    CLASS(ppm_t_operator_)             :: this
    INTEGER,               INTENT(OUT) :: info
END SUBROUTINE
!INIT
SUBROUTINE operator_discretize_on_(this,Discr,op_discr,info,method,&
        iargs,rargs,largs,with_ghosts,vector,interp,order) 
    IMPORT ppm_t_operator_,ppm_t_main_abstr,ppm_t_operator_discr,&
           ppm_kind_double
    CLASS(ppm_t_operator_)                     :: this
    CLASS(ppm_t_main_abstr),  INTENT(IN)       :: Discr
    CLASS(ppm_t_operator_discr),POINTER        :: op_discr
    INTEGER,                  INTENT(OUT)      :: info
    CHARACTER(LEN=*),    OPTIONAL,INTENT(IN)   :: method
    INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN)   :: iargs
    REAL(ppm_kind_double),DIMENSION(:),OPTIONAL,INTENT(IN)   :: rargs
    LOGICAL,DIMENSION(:),OPTIONAL,INTENT(IN)   :: largs
    LOGICAL,             OPTIONAL,INTENT(IN)   :: with_ghosts
    LOGICAL,             OPTIONAL,INTENT(IN)   :: interp
    LOGICAL,             OPTIONAL,INTENT(IN)   :: vector
    INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN)   :: order
END SUBROUTINE

!--------------------------------------
!DISCRETIZED OPERATORS
!--------------------------------------

SUBROUTINE operator_discr_destroy_(this,info)
    IMPORT ppm_t_operator_discr_
    CLASS(ppm_t_operator_discr_)                 :: this
    INTEGER,                       INTENT(OUT)   :: info
END SUBROUTINE

!COMPUTE OPERATOR ON A DISCRETIZATION (PARTICLE SET OR MESH)
SUBROUTINE operator_discr_compute_(this,Field_src,Field_to,info)
    IMPORT ppm_t_operator_discr_,ppm_t_main_abstr
    CLASS(ppm_t_operator_discr_)                 :: this
    CLASS(ppm_t_main_abstr),TARGET,INTENT(IN)    :: Field_src
    CLASS(ppm_t_main_abstr),TARGET,INTENT(INOUT) :: Field_to
    INTEGER,                       INTENT(OUT)   :: info
END SUBROUTINE

