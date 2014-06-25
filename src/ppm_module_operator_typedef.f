#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

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
      USE ppm_module_data
      USE ppm_module_util_functions
      USE ppm_module_interfaces
      USE ppm_module_particles_typedef
      USE ppm_module_vbp_typedef
      IMPLICIT NONE

      !----------------------------------------------------------------------
      ! Global parameters
      !----------------------------------------------------------------------
      !----------------------------------------------------------------------
      ! Local variables
      !----------------------------------------------------------------------

      INTEGER, PRIVATE, DIMENSION(3) :: ldc
      !!! Number of elements in all dimensions for allocation
      !----------------------------------------------------------------------
      ! Type declaration
      !----------------------------------------------------------------------
      TYPE,EXTENDS(ppm_t_operator_) :: ppm_t_operator
      !!! generic operator type
      CONTAINS
         PROCEDURE :: create        => operator_create
         PROCEDURE :: destroy       => operator_destroy
         PROCEDURE :: discretize_on => operator_discretize_on
      END TYPE
minclude ppm_create_collection(operator,operator,generate="extend")


#define  DTYPE(a) a/**/_s
#define  MK ppm_kind_single
#define  _MK _ppm_kind_single
#include "operator/dcop_typedef.f"

#define  DTYPE(a) a/**/_d
#define  MK ppm_kind_double
#define  _MK _ppm_kind_double
#include "operator/dcop_typedef.f"

      !----------------------------------------------------------------------
      ! Type-bound procedures
      !----------------------------------------------------------------------

      CONTAINS

minclude ppm_create_collection_procedures(operator,operator_)

#define DTYPE(a) a/**/_s
#define DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_single
#define __KIND __SINGLE_PRECISION
#include "operator/dcop_helpers.f"
#define __DIM 2
#include "operator/dcop_comp_weights.f"
#define __DIM 3
#include "operator/dcop_comp_weights.f"
#include "operator/dcop_typeproc.f"
#undef __KIND


#define DTYPE(a) a/**/_d
#define DEFINE_MK() INTEGER, PARAMETER :: MK = ppm_kind_double
#define __KIND __DOUBLE_PRECISION
#include "operator/dcop_helpers.f"
#define __DIM 2
#include "operator/dcop_comp_weights.f"
#define __DIM 3
#include "operator/dcop_comp_weights.f"
#include "operator/dcop_typeproc.f"
#undef __KIND


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
          CLASS(ppm_t_operator)               :: this
          INTEGER,              INTENT(  OUT) :: info
          !!! return 0 on success

          start_subroutine("operator_destroy")

          this%nterms = 0
          dealloc_pointer(this%degree)
          dealloc_pointer(this%coeffs)
          this%name = ""

          end_subroutine()
      END SUBROUTINE
      !COMPUTE
      SUBROUTINE operator_discretize_on(this,Discr_to,op_discr,opts,info, &
      &          Discr_src,nn_sq)
          CLASS(ppm_t_operator)                      :: this
          CLASS(ppm_t_discr_kind),TARGET,INTENT(IN)  :: Discr_to
          !!! Discretization on which the operator is applied (where the results are
          !!! stored)
          CLASS(ppm_t_operator_discr),POINTER        :: op_discr
          !!! Resulting discretized operator (on which %compute() is then called)
          CLASS(ppm_t_options)                       :: opts
          !!! Arguments for the discretized operator
          !!!    - method (ppm_param_op_fd, ppm_param_op_pse, ppm_param_op_dcpse)
          !!!    - with_ghosts (true if this operator is also applied to
          !!!             ghost particles, default is false)
          !!!    - vector (if the result should be a vector, default is false)
          !!!    - interp (if the operator is an interpolant, default is false)
          !!!    - order  (order of approximation, default is 2)
          !!!    - c      (ratio h/epsilon in the DC operators, default is 0.5)
          INTEGER,                       INTENT(OUT) :: info
          !!! return 0 on success
          CLASS(ppm_t_discr_kind),OPTIONAL,TARGET, INTENT(IN):: Discr_src
          !!! Discretization where the data comes from (source). By default, this
          !!! is the same as Discr_to
          CLASS(ppm_t_discr_data),OPTIONAL,        INTENT(IN):: nn_sq
          !!! Particle property where is stored the nearest-neighbour distance
          !!! (squared) of the source data point (only makes sense for interpolating
          !!! operators)

          !Local variables
          CLASS(ppm_t_discr_kind), POINTER :: Discr2

          start_subroutine("operator_discretize_on")

          IF (PRESENT(Discr_src)) THEN
             Discr2 => Discr_src
          ELSE
             Discr2 => Discr_to
          ENDIF

          SELECT TYPE(opts)
          CLASS IS (ppm_t_options_op)
             SELECT TYPE(Discr_to)
             CLASS IS (ppm_t_equi_mesh_)
                fail("mesh methods not yet implemented")
                SELECT CASE (opts%method)
                CASE (ppm_param_op_pse)
                   fail("Cannot use PSE on a mesh.")
                CASE (ppm_param_op_fd)
                   fail("Finite differences not yet implemented")
                CASE DEFAULT
                   fail("Method not recognized")
                END SELECT

             CLASS IS (ppm_t_particles_s_)
                SELECT CASE (opts%method)
                CASE (ppm_param_op_dcpse)
                   ALLOCATE(ppm_t_dcop_s::op_discr,stat=info)
                   fail("single-precision DCops not yet implemented (no big deal: just waiting for templating features...")

                CASE (ppm_param_op_pse)
                   !allocate(ppm_t_pseop::op_discr,stat=info)
                   fail("PSE method not yet implemented")

                CASE DEFAULT
                   fail("Method not recognized")

                END SELECT

             CLASS IS (ppm_t_particles_d_)
                SELECT CASE (opts%method)
                CASE (ppm_param_op_dcpse)
                   ALLOCATE(ppm_t_dcop_d::op_discr,stat=info)
                   SELECT TYPE(op_discr)
                   TYPE IS (ppm_t_dcop_d)
                      CALL op_discr%create(this,Discr_to,Discr2,info, &
                      &    opts%with_ghosts,opts%vector,opts%interp,  &
                      &    order=opts%order,order_v=opts%order_v,prop=nn_sq)
                      or_fail("op_discr%create failed")

                      CALL op_discr%comp_weights(info,c=opts%c)
                      or_fail("Failed to compute the weights of the DC-PSE operator")
                   END SELECT

                CASE (ppm_param_op_pse)
                   !allocate(ppm_t_pseop::op_discr,stat=info)
                   fail("PSE method not yet implemented")

                CASE (ppm_param_op_fd)
                   fail("cannot do finite differences on particles")

                CASE DEFAULT
                   fail("Method not recognized")

                END SELECT

             CLASS DEFAULT
                fail("Class not recognized for discretization argument")

             END SELECT !TYPE(Discr_to)

          CLASS DEFAULT
             fail("Wrong ppm_t_options type for operator option argument")

          END SELECT !TYPE(opts)

          end_subroutine()
      END SUBROUTINE

      END MODULE ppm_module_operator_typedef

