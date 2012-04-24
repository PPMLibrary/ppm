#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

#define __REAL 3 
#define __COMPLEX 4 
#define __INTEGER 5 
#define __LONGINT 6 
#define __LOGICAL 7 
#define __CHAR 8 

#define __crash_on_null_pointers  1
#undef __WITH_KDTREE

      MODULE ppm_module_operator_typedef
      !!! Declares differential operator data types
      !!!
      !!! [NOTE]
      !!! Most of the declared variables in this module should not be accessed
      !!! directly by the PPM client developer, they are used internally in the
      !!! library.

         !----------------------------------------------------------------------
         !  Modules
         !----------------------------------------------------------------------
         USE ppm_module_alloc
         USE ppm_module_data
         USE ppm_module_container_typedef
         USE ppm_module_interfaces
         USE ppm_module_error
         USE ppm_module_write
         USE ppm_module_substart
         USE ppm_module_substop

         IMPLICIT NONE
         !----------------------------------------------------------------------
         ! Global parameters
         !----------------------------------------------------------------------
         !----------------------------------------------------------------------
         ! Local variables
         !----------------------------------------------------------------------

         INTEGER, PRIVATE, DIMENSION(3)    :: ldc
         !!! Number of elements in all dimensions for allocation
         !----------------------------------------------------------------------
         ! Type declaration
         !----------------------------------------------------------------------
#define  DTYPE(a) a/**/_s
#define  MK ppm_kind_single
#define  _MK _ppm_kind_single
#include "operator/dcop_typedef.f"

#define  DTYPE(a) a/**/_d
#define  MK ppm_kind_double
#define  _MK _ppm_kind_double
#include "operator/dcop_typedef.f"


TYPE,EXTENDS(ppm_t_operator_) :: ppm_t_operator
    !!! generic operator type
    CONTAINS
    PROCEDURE :: create        => operator_create
    PROCEDURE :: destroy       => operator_destroy
    PROCEDURE :: discretize_on => operator_discretize_on
END TYPE
minclude define_collection_type(ppm_t_operator)


         !----------------------------------------------------------------------
         ! Type-bound procedures
         !----------------------------------------------------------------------


         CONTAINS
         
minclude define_collection_procedures(ppm_t_operator)

#define DTYPE(a) a/**/_s
#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_single
#include "operator/dcop_typeproc.f"

#define DTYPE(a) a/**/_d
#define  DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_double
#include "operator/dcop_typeproc.f"


!CREATE
SUBROUTINE operator_create(this,nterms,coeffs,degree,info,name)
    !!! Create a description for a differential operator
    CLASS(ppm_t_operator)                 :: this
    INTEGER,                INTENT(IN   ) :: nterms
    !!! Number of terms in the linear combination
    REAL(ppm_kind_double),DIMENSION(:),INTENT(IN):: coeffs
    !!! Multiplicative coefficients of each term in the linear combination of
    !!! differential operators
    INTEGER,DIMENSION(:),   INTENT(IN   ) :: degree
    !!! Degree of differentiation of each term
    INTEGER,                INTENT(OUT)   :: info
    !!! Returns status, 0 upon success.
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)  :: name
    !!! name for this operator

    start_subroutine("operator_create")

    !Check arguments
    IF (MINVAL(degree).LT.0) THEN
        fail("invalid degree: must be positive")
    ENDIF
    IF (SIZE(degree).NE.ppm_dim*nterms) THEN
        fail("wrong number of terms in degree argument")
    ENDIF
    IF (SIZE(coeffs).NE.nterms) THEN
        fail("wrong number of terms in coeffs argument")
    ENDIF


    !allocate operators descriptors
    ldc(1) = ppm_dim * nterms
    CALL ppm_alloc(this%degree,ldc,ppm_param_alloc_fit,info)
        or_fail_alloc("this%degree")
    ldc(1) = nterms
    CALL ppm_alloc(this%coeffs,ldc,ppm_param_alloc_fit,info)
        or_fail_alloc("this%coeffs")
    this%coeffs = coeffs 
    this%degree = degree 
    this%nterms = nterms 
    IF (PRESENT(name)) THEN
        this%name = TRIM(ADJUSTL(name))
    ELSE
        this%name = "no_name"
    ENDIF


    end_subroutine()
END SUBROUTINE
!DESTROY
SUBROUTINE operator_destroy(this,info)
    CLASS(ppm_t_operator)             :: this
    INTEGER,               INTENT(OUT) :: info
    start_subroutine("operator_destroy")

    this%nterms = 0
    dealloc_pointer(this%degree)
    dealloc_pointer(this%coeffs)
    this%name = ""

    end_subroutine()
END SUBROUTINE
!COMPUTE
SUBROUTINE operator_discretize_on(this,Discr_src,op_discr,info,method,&
        iargs,rargs,largs,with_ghosts,vector,interp,order,Discr_to) 
    CLASS(ppm_t_operator)                      :: this
    CLASS(ppm_t_discr_kind),TARGET,  INTENT(IN):: Discr_src
    CLASS(ppm_t_operator_discr),POINTER        :: op_discr
    INTEGER,                  INTENT(OUT)      :: info
    CHARACTER(LEN=*),    OPTIONAL,INTENT(IN)   :: method
    INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN)   :: iargs
    REAL(ppm_kind_double),DIMENSION(:),OPTIONAL,INTENT(IN)   :: rargs
    LOGICAL,DIMENSION(:),OPTIONAL,INTENT(IN)   :: largs
    LOGICAL,             OPTIONAL,INTENT(IN)   :: with_ghosts
    LOGICAL,             OPTIONAL,INTENT(IN)   :: vector
    LOGICAL,             OPTIONAL,INTENT(IN)   :: interp
    INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN)   :: order
    !!! Order of approximation for each term
    CLASS(ppm_t_discr_kind),OPTIONAL,TARGET, INTENT(IN):: Discr_to
    

    CLASS(ppm_t_discr_kind),POINTER :: Discr2 => NULL()

    start_subroutine("operator_discretize_on")

    IF (PRESENT(Discr_to)) THEN
        Discr2 => Discr_to
    ELSE
        Discr2 => Discr_src
    ENDIF
    
    SELECT TYPE(Discr_src)
    CLASS IS (ppm_t_equi_mesh_)
        fail("mesh methods not yet implemented")
    CLASS IS (ppm_t_particles_s_)
        SELECT CASE (method)
        CASE ("DC-PSE")
            allocate(ppm_t_dcop_s::op_discr,stat=info)
            fail("single-precision DCops not yet implemented (no big deal: just waiting for templating features...")
        CASE ("PSE")
            !allocate(ppm_t_pseop::op_discr,stat=info)
            fail("PSE method not yet implemented")
        CASE DEFAULT
            fail("Method not recognized")
        END SELECT
    CLASS IS (ppm_t_particles_d_)
        SELECT CASE (method)
        CASE ("DC-PSE")
            allocate(ppm_t_dcop_d::op_discr,stat=info)
            CALL op_discr%create(Discr_src,Discr2,info,this%nterms,&
                with_ghosts,vector,interp,order=order)

        CASE ("PSE")
            !allocate(ppm_t_pseop::op_discr,stat=info)
            fail("PSE method not yet implemented")
        CASE DEFAULT
            fail("Method not recognized")
        END SELECT
    CLASS DEFAULT
       fail("Class not recognized for discretization argument")
    END SELECT

    end_subroutine()
END SUBROUTINE


         END MODULE ppm_module_operator_typedef

