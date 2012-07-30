      MODULE ppm_module_options
      !!! Declares options data types
      !!!
         !----------------------------------------------------------------------
         !  Modules
         !----------------------------------------------------------------------
         USE ppm_module_alloc
         USE ppm_module_data
         USE ppm_module_error
         USE ppm_module_write
         USE ppm_module_substart
         USE ppm_module_substop

         IMPLICIT NONE

TYPE,ABSTRACT  ::  ppm_t_options
END TYPE

TYPE,EXTENDS(ppm_t_options)        :: ppm_t_options_op
    INTEGER                        :: method = 0
    !!! one of
    !!!  - ppm_param_op_fd
    !!!  - ppm_param_op_pse
    !!!  - ppm_param_op_dcpse
    LOGICAL                        :: with_ghosts = .FALSE.
    LOGICAL                        :: vector = .FALSE.
    LOGICAL                        :: interp = .FALSE.
    INTEGER,DIMENSION(:),POINTER   :: order_v => NULL()
    INTEGER                        :: order = 2
    REAL(ppm_kind_double)          :: c = 0.5D0
    !!! - for dcpse operators

    CONTAINS
        PROCEDURE :: create   => options_op_create
        PROCEDURE :: destroy  => options_op_destroy
END TYPE


         !----------------------------------------------------------------------
         ! Type-bound procedures
         !----------------------------------------------------------------------


         CONTAINS
         
SUBROUTINE options_op_create(this,method,info,with_ghosts,vector,&
        interp,order,order_v,c)
    CLASS(ppm_t_options_op)        :: this
    !!! option data type for discretized operator
    INTEGER                               :: method
    INTEGER                               :: info
    LOGICAL, OPTIONAL                     :: with_ghosts
    LOGICAL, OPTIONAL                     :: vector
    LOGICAL, OPTIONAL                     :: interp
    INTEGER,DIMENSION(:), OPTIONAL        :: order_v
    !!! order of approximation (passed as a vector)
    !!! one number per term in the operator
    INTEGER,              OPTIONAL        :: order
    !!! order of approximation (passed as a scalar)
    !!! same number for all terms in the operator
    REAL(ppm_kind_double),OPTIONAL        :: c

    start_subroutine("options_op_create")

    this%method = method

    IF (PRESENT(with_ghosts)) this%with_ghosts=with_ghosts
    IF (PRESENT(vector)) this%vector=vector
    IF (PRESENT(interp)) this%interp=interp
    check_false(<#PRESENT(order_v).AND.PRESENT(order)#>,&
        "provide values for either order_v or order, not both")
    IF (PRESENT(order_v)) THEN
        allocate(this%order_v(size(order_v)),STAT=info)
            or_fail_alloc("order")
        this%order_v=order_v
    ENDIF
    IF (PRESENT(order)) THEN
        this%order=order
        dealloc_pointer("this%order_v")
    ENDIF
    IF (PRESENT(c)) this%c=c


    end_subroutine()
END SUBROUTINE

SUBROUTINE options_op_destroy(this,info)
    CLASS(ppm_t_options_op)              :: this
    INTEGER                              :: info
    start_subroutine("options_op_destroy")

    dealloc_pointer("this%order_v")
    this%order = 2
    this%method = 0
    this%c = 0.5D0
    this%vector =.FALSE.
    this%interp =.FALSE.
    this%with_ghosts =.FALSE.

    end_subroutine()
END SUBROUTINE



 END MODULE ppm_module_options

