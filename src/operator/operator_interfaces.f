!--------------------------------------
!GENERIC OPERATOR CLASS
!--------------------------------------

!CREATE
SUBROUTINE operator_create_(this,nterms,coeffs,degree,info,name)
    IMPORT ppm_t_operator_,ppm_kind_double
    CLASS(ppm_t_operator_)                :: this
    INTEGER,                INTENT(IN   ) :: nterms
    REAL(ppm_kind_double),DIMENSION(:),INTENT(IN):: coeffs
    INTEGER,DIMENSION(:),   INTENT(IN   ) :: degree
    INTEGER,                INTENT(OUT)   :: info
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)  :: name
END SUBROUTINE
!DESTROY
SUBROUTINE operator_destroy_(this,info)
    IMPORT ppm_t_operator_
    CLASS(ppm_t_operator_)             :: this
    INTEGER,               INTENT(OUT) :: info
END SUBROUTINE
!INIT
SUBROUTINE operator_discretize_on_(this,Discr_to,op_discr,opts,&
        info,Discr_src,nn_sq) 
    IMPORT ppm_t_operator_,ppm_t_discr_kind,ppm_t_operator_discr,&
           ppm_kind_double,ppm_t_options,ppm_t_discr_data
    CLASS(ppm_t_operator_)                     :: this
    CLASS(ppm_t_discr_kind),TARGET,  INTENT(IN):: Discr_to
    CLASS(ppm_t_operator_discr),POINTER        :: op_discr
    CLASS(ppm_t_options)                       :: opts
    INTEGER,                  INTENT(OUT)      :: info
    CLASS(ppm_t_discr_kind),OPTIONAL,TARGET,INTENT(IN):: Discr_src
    CLASS(ppm_t_discr_data),OPTIONAL,        INTENT(IN):: nn_sq
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

