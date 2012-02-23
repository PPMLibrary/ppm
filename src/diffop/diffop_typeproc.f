!BEGIN
FUNCTION DTYPE(begin_op)(this) RESULT (iterator)
    CLASS(DTYPE(ppm_c_operators)),  INTENT(INOUT)   :: this
    TYPE(DTYPE(ppm_t_operator)),   POINTER      :: iterator

    this%iter_id = this%min_id
    IF (this%nb.GT.0) THEN
        iterator => this%vec(this%min_id)%t
    ELSE
        iterator => NULL()
    ENDIF

END FUNCTION DTYPE(begin_op)
!NEXT
FUNCTION DTYPE(next_op)(this) RESULT (iterator)
    CLASS(DTYPE(ppm_c_operators)),  INTENT(INOUT)   :: this
    TYPE(DTYPE(ppm_t_operator)),   POINTER      :: iterator

    DO WHILE(this%iter_id.LT.this%max_id)
        this%iter_id = this%iter_id + 1
        IF (ASSOCIATED(this%vec(this%iter_id)%t)) THEN 
            iterator => this%vec(this%iter_id)%t
            RETURN
        ENDIF
    ENDDO
    iterator => NULL()
    RETURN
       
END FUNCTION DTYPE(next_op)
!PREVIOUS
FUNCTION DTYPE(prev_op)(this) RESULT (iterator)
    CLASS(DTYPE(ppm_c_operators)),  INTENT(INOUT)   :: this
    TYPE(DTYPE(ppm_t_operator)),   POINTER      :: iterator

    DO WHILE(this%iter_id.GT.this%min_id)
        this%iter_id = this%iter_id - 1
        IF (ASSOCIATED(this%vec(this%iter_id)%t)) THEN 
            iterator => this%vec(this%iter_id)%t
            RETURN
        ENDIF
    ENDDO
    iterator => NULL()
    RETURN
       
END FUNCTION DTYPE(prev_op)
!LAST
FUNCTION DTYPE(last_op)(this) RESULT (iterator)
    CLASS(DTYPE(ppm_c_operators)),  INTENT(INOUT)   :: this
    TYPE(DTYPE(ppm_t_operator)),   POINTER      :: iterator

    this%iter_id = this%max_id
    IF (this%nb.GT.0) THEN
        iterator => this%vec(this%max_id)%t
    ELSE
        iterator => NULL()
    ENDIF

END FUNCTION DTYPE(last_op)
!DESTROY CONTAINER
SUBROUTINE DTYPE(op_container_destroy)(this,info)
    CLASS(DTYPE(ppm_c_operators)), INTENT(INOUT)   :: this
    TYPE(DTYPE(ppm_t_operator)),POINTER :: p => NULL()
    INTEGER, INTENT(OUT) :: info

    info = 0
    p => this%begin()
    DO WHILE(ASSOCIATED(p))
        CALL p%destroy(info)
        p => this%next()
    ENDDO
    IF (ASSOCIATED(this%vec))  DEALLOCATE(this%vec)
    NULLIFY(this%vec)
    this%min_id = HUGE(1)
    this%max_id = 0
    this%nb = 0
    this%vec_size=0

END SUBROUTINE DTYPE(op_container_destroy)
!
FUNCTION DTYPE(op_exists)(cont,id,caller) RESULT(exists)
    !!!------------------------------------------------------------------------!
    !!! Check whether a neighbor list exists and can be accessed at this id
    !!!------------------------------------------------------------------------!
    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    CLASS(DTYPE(ppm_c_operators))                       :: cont
    !!! Data structure containing the particles
    INTEGER,                            INTENT(IN   )   :: id
    !!! id where the data is stored
    LOGICAL                                             :: exists
    !!! Return status, on success 0.
    CHARACTER(LEN = *),OPTIONAL                         :: caller
    !!! Calling routine

#include "../cont/container_exists.inc"

END FUNCTION DTYPE(op_exists)
!CREATE ENTRY
SUBROUTINE DTYPE(op_create)(op,nterms,coeffs,degree,order,&
        name,with_ghosts,vector,interp,pid,nlid,info)
    !!! Create a differential operator
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_operator))          :: op
    INTEGER,                INTENT(IN   ) :: nterms
    !!! Number of terms in the linear combination
    REAL(MK),DIMENSION(:),  INTENT(IN   ) :: coeffs
    !!! Multiplicative coefficients of each term in the linear combination of
    !!! differential operators
    INTEGER,DIMENSION(:),   INTENT(IN   ) :: degree
    !!! Degree of differentiation of each term
    INTEGER,DIMENSION(:),   INTENT(IN   ) :: order
    !!! Order of approxmiation for each term
    LOGICAL,                INTENT(IN   ) :: with_ghosts
    !!! True if the operator should be computed for ghost particles too. 
    !!! Note that the resulting values will be wrong for the ghost particles
    !!! that have some neighbours outside the ghost layers. Default is false.
    LOGICAL,                INTENT(IN   ) :: vector
    !!! True if the operator is a vector field. Default is false.
    LOGICAL,                INTENT(IN   ) :: interp
    !!! True if the operator interpolates data from one set of particles to
    !!! another. Default is false.
    INTEGER,                INTENT(IN   ) :: pid
    !!! Id of the set of particles that this operator takes data from.
    !!! The default, 0, stands for "self" (the operator is computed
    !!! on the same set of particles than the one which contains the data).
    INTEGER,                INTENT(IN   ) :: nlid
    !!! Id of the neighbour list that should be used
    !!! The default, 1, refers to "self": the list of neighbours within
    !!! the same set of particles. 
    CHARACTER(LEN=*)                      :: name
    !!! name for this operator
    INTEGER,                INTENT(OUT)   :: info
    !!! Returns status, 0 upon success.

    CHARACTER(LEN=ppm_char)               :: caller = 'op_create'
    REAL(KIND(1.D0))                      :: t0
    
    CALL substart(caller,t0,info)

    op%flags = .FALSE.
    op%flags(ppm_ops_inc_ghosts) = with_ghosts
    op%flags(ppm_ops_interp) = interp
    op%flags(ppm_ops_vector) = vector
    op%flags(ppm_ops_isdefined) = .TRUE.
    op%P_id = pid
    op%neigh_id = nlid

    IF (ASSOCIATED(op%desc)) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_sub_failed,caller,   &
            &       'operator struct not clean. Use destroy first ',&
            &       __LINE__,info)
        GOTO 9999
    ENDIF

    ALLOCATE(op%desc,STAT=info)
    IF (info.NE.0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &       'allocation of ker or desc failed',__LINE__,info)
        GOTO 9999
    ENDIF

    CALL op%desc%create(nterms,coeffs,degree,order,name,info)

    CALL substop(caller,t0,info)
    
    9999 CONTINUE

END SUBROUTINE DTYPE(op_create)
!DESTROY ENTRY
SUBROUTINE DTYPE(op_destroy)(op,info)
    !!! Destroy the description for a differential operator
    CLASS(DTYPE(ppm_t_operator))              :: op
    INTEGER                                   :: i
    INTEGER,                   INTENT(  OUT)  :: info
    !!! Returns status, 0 upon success.
    REAL(KIND(1.D0))                   :: t0
    CHARACTER(LEN=ppm_char)            :: caller = 'op_destroy'
    
    CALL substart(caller,t0,info)

    CALL ppm_alloc(op%ker,ldc,ppm_param_dealloc,info)
    IF (info.NE.0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_dealloc,caller,   &
            &       'ker deallocate failed ',__LINE__,info)
        GOTO 9999
    ENDIF
    CALL op%desc%destroy(info)
    IF (info.NE.0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_sub_failed,caller,   &
            &       'desc destroy failed ',__LINE__,info)
        GOTO 9999
    ENDIF

    op%flags = .FALSE.
    op%P_id = -1
    op%neigh_id = 1

    CALL substop(caller,t0,info)

    9999 CONTINUE

END SUBROUTINE DTYPE(op_destroy)

#undef DEFINE_MK
#undef DTYPE


