#define CONTAINER DTYPE(ppm_c_operators)
#define __CONTAINER(a) DTYPE(ppm_c_operators)_/**/a
#define VEC_TYPE DTYPE(ppm_t_ptr_ops)
#define ITERATOR_TYPE DTYPE(ppm_t_operator)
#include "cont/extended_container_typeproc.f"

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

SUBROUTINE DTYPE(desc_create)(desc,nterms,coeffs,degree,order,name,info)
    !!! Create a description for a differential operator
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_opdesc))            :: desc
    INTEGER,                INTENT(IN   ) :: nterms
    !!! Number of terms in the linear combination
    REAL(MK),DIMENSION(:),  INTENT(IN   ) :: coeffs
    !!! Multiplicative coefficients of each term in the linear combination of
    !!! differential operators
    INTEGER,DIMENSION(:),   INTENT(IN   ) :: degree
    !!! Degree of differentiation of each term
    INTEGER,DIMENSION(:),   INTENT(IN   ) :: order
    !!! Order of approxmiation for each term
    CHARACTER(LEN=*)                      :: name
    !!! name for this operator
    INTEGER,                INTENT(OUT)   :: info
    !!! Returns status, 0 upon success.

    CHARACTER(LEN=ppm_char)               :: caller = 'desc_create'
    REAL(KIND(1.D0))                      :: t0
    
    CALL substart(caller,t0,info)

    !Check arguments
    IF (MINVAL(degree).LT.0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &       'invalid degree: must be positive',__LINE__,info)
        GOTO 9999
    ENDIF
    IF (MINVAL(order).LT.0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &     'invalid approx order: must be positive',__LINE__,info)
        GOTO 9999
    ENDIF
    IF (SIZE(degree).NE.ppm_dim*nterms) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &     'wrong number of terms in degree argument',__LINE__,info)
        GOTO 9999
    ENDIF
    IF (SIZE(order).NE.nterms) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &      'wrong number of terms in order argument',__LINE__,info)
        GOTO 9999
    ENDIF
    IF (SIZE(coeffs).NE.nterms) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &      'wrong number of terms in coeffs argument',__LINE__,info)
        GOTO 9999
    ENDIF


    !allocate operators descriptors
    ldc(1) = ppm_dim * nterms
    CALL ppm_alloc(desc%degree,ldc,ppm_param_alloc_fit,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &            'failed to allocate ops%desc',__LINE__,info)
        GOTO 9999
    ENDIF
    ldc(1) = nterms
    CALL ppm_alloc(desc%order,ldc,ppm_param_alloc_fit,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &            'failed to allocate ops%desc',__LINE__,info)
        GOTO 9999
    ENDIF
    ldc(1) = nterms
    CALL ppm_alloc(desc%coeffs,ldc,ppm_param_alloc_fit,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &            'failed to allocate ops%desc',__LINE__,info)
        GOTO 9999
    ENDIF
    desc%order = order 
    desc%coeffs = coeffs 
    desc%degree = degree 
    desc%nterms = nterms 
    desc%name = name


    CALL substop(caller,t0,info)

    9999 CONTINUE

END SUBROUTINE DTYPE(desc_create)

SUBROUTINE DTYPE(desc_destroy)(desc,info)
    CLASS(DTYPE(ppm_t_opdesc))              :: desc
    INTEGER,                 INTENT(  OUT)  :: info
    !!! Returns status, 0 upon success.
    
    IF (ASSOCIATED(desc%degree)) DEALLOCATE(desc%degree,STAT=info)
    IF (ASSOCIATED(desc%order))  DEALLOCATE(desc%order,STAT=info)
    IF (ASSOCIATED(desc%coeffs)) DEALLOCATE(desc%coeffs,STAT=info)

END SUBROUTINE DTYPE(desc_destroy)

SUBROUTINE DTYPE(part_op_create)(Op,Pc,id,nterms,coeffs,degree,order,info,&
        name,with_ghosts,vector,interp,P_id_from,P_from,neigh_id)
    !!! Adds a differential operator to a particle set
    !!!------------------------------------------------------------------------!
    !!! Define a DC operator as a linear combination (with scalar coefficients)
    !!! of nterms partial derivatives of arbitrary degrees. 
    !!! These are given by a matrix
    !!! of integers where each row represents one term of the linear combination
    !!! and each of the ppm_dim columns is the order of differentiation in that
    !!! dimension.
    !!! The definition of the operator is stored in the ppm_t_operator derived 
    !!! type under the index eta_id. 
    !!! The operator itself is computed elsewhere and will be stored in 
    !!! the same data structure.
    !!!
    !!! Usage example:
    !!!
    !!!   The differential operator:
    !!!   3.0 df/dx -7.0 d^4f/dxdydz^2 + 8.0 d^3f/dx^2dz
    !!!   would be defined by calling particles_dcop_define with
    !!!   coeffs = (/3.0, -7.0, 8.0/)
    !!!   degree = (/1,0,0,  1,1,2,  2,0,1 /)
    !!!   order =  (/2,      1,      3     /)
    !!!   nterms = 3
    !!!------------------------------------------------------------------------!
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_operator))                        :: Op
    CLASS(DTYPE(ppm_t_particles))                       :: Pc
    INTEGER,                            INTENT(  OUT)   :: id
    !!! id for this operator 
    INTEGER,                            INTENT(IN   )   :: nterms
    !!! Number of terms in the linear combination
    REAL(MK),DIMENSION(:),              INTENT(IN   )   :: coeffs
    !!! Multiplicative coefficients of each term in the linear combination of
    !!! differential operators
    INTEGER,DIMENSION(:),               INTENT(IN   )   :: degree
    !!! Degree of differentiation of each term
    INTEGER,DIMENSION(:),               INTENT(IN   )   :: order
    !!! Order of approxmiation for each term
    LOGICAL,OPTIONAL,                   INTENT(IN   )   :: with_ghosts
    !!! True if the operator should be computed for ghost particles too. 
    !!! Note that the resulting values will be wrong for the ghost particles
    !!! that have some neighbours outside the ghost layers. Default is false.
    LOGICAL,OPTIONAL,                   INTENT(IN   )   :: vector
    !!! True if the operator is a vector field. Default is false.
    LOGICAL,OPTIONAL,                   INTENT(IN   )   :: interp
    !!! True if the operator interpolates data from one set of particles to
    !!! another. Default is false.
    CLASS(DTYPE(ppm_t_particles)),OPTIONAL,INTENT(IN)   :: P_from
    !!! The set of particles that this operator takes data from.
    !!! The default is "self" (P_from = Pc and the operator is computed
    !!! on the same set of particles than the one which contains the data).
    INTEGER,OPTIONAL,                   INTENT(IN   )   :: neigh_id
    !!! Id of the neighbour list that should be used
    !!! The default, 1, refers to "self": the list of neighbours within
    !!! the same set of particles. 
    CHARACTER(LEN=*) , OPTIONAL           :: name
    !!! name for this operator
    INTEGER,               INTENT(OUT)    :: info

    INTEGER                               :: i,vec_size,npart,lpid,lnlid
    CHARACTER(LEN=ppm_char)               :: caller = 'particle_op_create'
    CHARACTER(LEN=ppm_char)               :: lname
    LOGICAL                               :: lwith_ghosts,lvector,linterp
    REAL(KIND(1.D0))                      :: t0
    TYPE(DTYPE(ppm_t_ptr_ops)),DIMENSION(:),POINTER  :: vec_tmp => NULL()
    TYPE(DTYPE(ppm_t_operator)),          POINTER  :: op => NULL()

    CALL substart(caller,t0,info)

    !Generate a new id (we should use templating here...)
    ASSOCIATE (cont => Pc%ops )
        id = 0
        IF (cont%nb.LT.cont%vec_size) THEN
            !there is at least one empty slot in the array
            ! of property pointers
            id = id + 1
            DO WHILE (ASSOCIATED(cont%vec(id)%t))
                id = id + 1
            ENDDO
        ELSE
            IF (.NOT. ASSOCIATED(cont%vec)) THEN
                !need to allocate the array of property pointers 
                vec_size=20
                ALLOCATE(cont%vec(vec_size),STAT=info)
                IF (info .NE. 0) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_alloc,caller,&
                        'allocating property array failed',__LINE__,info)
                    GOTO 9999
                ENDIF
                id = 1
            ELSE
                !need to resize the array of property pointers 
                vec_size=MAX(2*cont%vec_size,20)
                ALLOCATE(vec_tmp(vec_size),STAT=info)
                IF (info .NE. 0) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_alloc,caller,&
                        'allocating property array failed',__LINE__,info)
                    GOTO 9999
                ENDIF
                DO i=1,cont%vec_size
                    vec_tmp(i)%t => cont%vec(i)%t
                ENDDO
                DEALLOCATE(cont%vec)
                cont%vec => vec_tmp
            ENDIF
            cont%vec_size = vec_size
            id = cont%nb + 1
        ENDIF
        cont%nb = cont%nb + 1
            

        IF (id .GT. cont%max_id) cont%max_id = id
        IF (id .LT. cont%min_id) cont%min_id = id

    END ASSOCIATE

    !Allocate operator struct
    IF (.NOT. ASSOCIATED(Pc%ops%vec(id)%t)) THEN
        ALLOCATE(Pc%ops%vec(id)%t,STAT=info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,&
                'allocating operator pointer failed',__LINE__,info)
            GOTO 9999
        ENDIF
    ENDIF

    op => Pc%ops%vec(id)%t

    IF (PRESENT(name)) THEN
        lname = name
    ELSE
        lname = particles_dflt_opname(id)
    ENDIF
    IF (PRESENT(with_ghosts)) THEN
        lwith_ghosts = with_ghosts
    ELSE
        lwith_ghosts = .FALSE.
    ENDIF
    IF (PRESENT(interp)) THEN
        linterp = interp
    ELSE
        linterp = .FALSE.
    ENDIF
    IF (PRESENT(vector)) THEN
        lvector = vector
    ELSE
        lvector = .FALSE.
    ENDIF
    IF (PRESENT(P_id)) THEN
        lpid = P_id
    ELSE
        lpid = 0
    ENDIF
    IF (PRESENT(neigh_id)) THEN
        lnlid = neigh_id
    ELSE
        lnlid = ppm_param_default_nlID
    ENDIF

    IF (.NOT. Pc%neighs%exists(lnlid)) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,&
            'Invalid neighbour list. Use comp_neigh() first.',&
            __LINE__,info)
        GOTO 9999
    ENDIF
    IF (Pc%neighs%vec(lnlid)%t%P_id .NE. lpid) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,&
            'incompatible P_id and neigh_id',__LINE__,info)
        GOTO 9999
    ENDIF

    ! Create/Initialize operator
    CALL op%create(nterms,coeffs,degree,order,&
        lname,lwith_ghosts,lvector,linterp,lpid,lnlid,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_sub_failed,caller,&
            'creating operator object failed',__LINE__,info)
        GOTO 9999
    ENDIF

    CALL substop(caller,t0,info)
    9999  CONTINUE

END SUBROUTINE DTYPE(part_op_create)


SUBROUTINE DTYPE(part_op_destroy)(Pc,id,info)
    !!! Destroy a property from an existing particle set
    CLASS(DTYPE(ppm_t_particles))         :: Pc
    INTEGER,                INTENT(INOUT) :: id
    INTEGER,               INTENT(OUT)    :: info

    CHARACTER(LEN=ppm_char)               :: caller = 'particle_op_destroy'
    REAL(KIND(1.D0))                      :: t0

    CALL substart(caller,t0,info)

    ASSOCIATE (cont => Pc%ops)
        IF (id .LE. 0 .OR. id .GT. cont%vec_size) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,&
                &    'property id larger than size of properties array',&
                __LINE__,info)
            GOTO 9999
        ENDIF

        CALL cont%vec(id)%t%destroy(info)
        NULLIFY(cont%vec(id)%t)

        cont%nb = cont%nb - 1
        IF (id .EQ. cont%max_id) THEN
            cont%max_id = cont%max_id - 1
            IF (cont%max_id .GT. 0) THEN
                DO WHILE(.NOT.ASSOCIATED(cont%vec(cont%max_id)%t))
                    cont%max_id = cont%max_id - 1
                    IF (cont%max_id .EQ. 0) EXIT
                ENDDO
            ENDIF
        ENDIF
        IF (cont%nb.EQ.0) THEN
            cont%min_id = HUGE(1)
        ELSE IF (id .EQ. cont%min_id) THEN
            cont%min_id = cont%min_id + 1
            IF (cont%min_id .LE. cont%vec_size) THEN
                DO WHILE(.NOT.ASSOCIATED(cont%vec(cont%min_id)%t))
                    cont%min_id = cont%min_id + 1
                    IF (cont%min_id .GT. cont%vec_size) THEN
                        info = ppm_error_error
                        CALL ppm_error(ppm_err_alloc,caller,&
                            &    'coding error in the data structure',&
                            __LINE__,info)
                        GOTO 9999
                    ENDIF
                ENDDO
            ENDIF
        ENDIF
    END ASSOCIATE

    CALL substop(caller,t0,info)

    9999  CONTINUE

END SUBROUTINE DTYPE(part_op_destroy)

SUBROUTINE DTYPE(part_op_compute)(Pc,op_id,info,c,min_sv)

    USE ppm_module_write
    IMPLICIT NONE
#ifdef __MPI
    INCLUDE 'mpif.h'
#endif

    DEFINE_MK()
    !---------------------------------------------------------
    ! arguments
    !---------------------------------------------------------
    CLASS(DTYPE(ppm_t_particles))                        :: Pc
    !!! particles
    INTEGER,                             INTENT(IN   )   :: op_id
    !!! id of the operator 
    INTEGER,                             INTENT(  OUT)   :: info
    !!! non-zero on output if some error occurred
    !---------------------------------------------------------
    ! Optional arguments
    !---------------------------------------------------------
    REAL(MK),OPTIONAL                       :: c
    !!! ratio h/epsilon (default is 1.0)
    REAL(MK),OPTIONAL   ,  INTENT(  OUT)    :: min_sv
    !!! smallest singular value
    !---------------------------------------------------------
    ! local variables
    !---------------------------------------------------------
    CHARACTER(LEN = ppm_char)               :: caller = 'part_dcop_compute'
    CHARACTER(LEN = ppm_char)               :: cbuf
    REAL(KIND(1.D0))                        :: t0,t1,t2
    TYPE(DTYPE(ppm_t_operator)), POINTER    :: op => NULL()

    !-------------------------------------------------------------------------
    ! Initialize
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)
#ifdef __MPI
    t1 = MPI_WTIME(info)
#endif

    !-------------------------------------------------------------------------
    ! Check arguments
    !-------------------------------------------------------------------------
    IF (.NOT. ASSOCIATED(Pc%xp)) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,'Particles not defined',&
            __LINE__,info)
        GOTO 9999
    ENDIF
    IF (.NOT. Pc%ops%exists(op_id)) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            & 'No operator data structure found, use create_op() first',&
            __LINE__,info)
        GOTO 9999
    ENDIF
    op => Pc%ops%vec(op_id)%t
    IF (.NOT. op%flags(ppm_ops_isdefined)) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            & 'Operator not found, use create_op() first',&
            __LINE__,info)
        GOTO 9999
    ENDIF
    IF (op%flags(ppm_ops_iscomputed)) THEN
        WRITE(cbuf,*) 'WARNING: The operator with id ',op_id,&
            & ' and name *',TRIM(ADJUSTL(op%desc%name)),&
            &'* seems to have already been computed. Unnecessary call to',&
            &' particles_dcop_compute()'
        CALL ppm_write(ppm_rank,caller,cbuf,info)
    ENDIF

    !-------------------------------------------------------------------------
    ! Compute the DC operator
    !-------------------------------------------------------------------------

    Pc%stats%nb_dc_comp = Pc%stats%nb_dc_comp + 1

    IF (ppm_dim .EQ. 2) THEN
        CALL Pc%DTYPE(ppm_dcop_compute2d)(op_id,info,c,min_sv)
    ELSE
        CALL Pc%DTYPE(ppm_dcop_compute3d)(op_id,info,c,min_sv)
    ENDIF
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            & 'ppm_dcop_compute failed',&
            __LINE__,info)
        GOTO 9999
    ENDIF

    !-------------------------------------------------------------------------
    ! Update states
    !-------------------------------------------------------------------------
    op%flags(ppm_ops_iscomputed) = .TRUE.
#ifdef __MPI
    t2 = MPI_WTIME(info)
    Pc%stats%t_dc_comp = Pc%stats%t_dc_comp + (t2-t1)
#endif

    !-------------------------------------------------------------------------
    ! Finalize
    !-------------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999 CONTINUE ! jump here upon error

END SUBROUTINE DTYPE(part_op_compute)

SUBROUTINE DTYPE(part_op_apply)(Pc,from_id,to_id,op_id,info)
    !!!------------------------------------------------------------------------!
    !!! Apply DC kernel stored in op_id to the scalar property stored
    !!! prop_from_id and store the results in prop_to_id
    !!!------------------------------------------------------------------------!

#ifdef __MPI
    INCLUDE "mpif.h"
#endif

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_particles))                       :: Pc
    !!! Data structure containing the particles
    INTEGER,                            INTENT(IN   )   :: from_id
    !!! id where the data is stored
    INTEGER,                            INTENT(INOUT)   :: to_id
    !!! id where the result should be stored (0 if it needs to be allocated)
    INTEGER,                            INTENT(IN   )   :: op_id
    !!! id where the DC kernel has been stored
    INTEGER,                            INTENT(  OUT)   :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    ! local variables
    !-------------------------------------------------------------------------
    CHARACTER(LEN = ppm_char)                  :: filename
    CHARACTER(LEN = ppm_char)                  :: caller = 'part_dcop_apply'
    INTEGER                                    :: ip,iq,ineigh,lda,np_target
    REAL(KIND(1.D0))                           :: t0,t1,t2
    REAL(MK),DIMENSION(:,:),POINTER            :: eta => NULL()
    REAL(MK),DIMENSION(:),  POINTER            :: wps1 => NULL(),wps2=>NULL()
    REAL(MK),DIMENSION(:,:),POINTER            :: wpv1 => NULL(),wpv2=>NULL()
    REAL(MK),DIMENSION(:),  POINTER            :: dwps => NULL()
    REAL(MK),DIMENSION(:,:),POINTER            :: dwpv => NULL()
    INTEGER, DIMENSION(:),  POINTER            :: nvlist => NULL()
    INTEGER, DIMENSION(:,:),POINTER            :: vlist => NULL()
    REAL(MK)                                   :: sig
    LOGICAL                                    :: vector_output
    LOGICAL                                    :: vector_input
    LOGICAL                                    :: with_ghosts,isinterp

    TYPE(DTYPE(ppm_t_sop)),POINTER             :: Pc2 => NULL()
    TYPE(DTYPE(ppm_t_neighlist)),POINTER       :: Nlist => NULL()
    TYPE(DTYPE(ppm_t_operator)), POINTER       :: op => NULL()
    TYPE(DTYPE(ppm_t_part_prop)), POINTER      :: prop_from => NULL()
    !-------------------------------------------------------------------------
    ! Initialize
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)
#ifdef __MPI
    t1 = MPI_WTIME(info)
#endif

    !-------------------------------------------------------------------------
    ! Check arguments
    !-------------------------------------------------------------------------
    IF (.NOT. ASSOCIATED(Pc%xp)) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,'Particles not defined',&
            __LINE__,info)
        GOTO 9999
    ENDIF
    IF (.NOT. Pc%ops%exists(op_id)) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            & 'No operator data structure found, use create_op() first',&
            __LINE__,info)
        GOTO 9999
    ENDIF
    op => Pc%ops%vec(op_id)%t
    IF (.NOT. op%flags(ppm_ops_isdefined)) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            & 'Operator not found, use create_op() first',&
            __LINE__,info)
        GOTO 9999
    ENDIF
    IF (.NOT.op%flags(ppm_ops_iscomputed)) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            & 'Operator not computed, use comp_op() first',&
            __LINE__,info)
        GOTO 9999
    ENDIF

    isinterp = op%flags(ppm_ops_interp)
    vector_output =  op%flags(ppm_ops_vector)
    !if true, then each term of the differential opearator is stored as one
    !component in eta. This is used when computing e.g. the gradient opearator.
    !if false, the same input parameters would yield an operator approximating
    ! the divergence operator.
    with_ghosts = op%flags(ppm_ops_inc_ghosts)
    !if true, then the operator should be computed for ghost particles too. 
    !Note that the resulting values will be wrong for the ghost particles
    !that have some neighbours outside the ghost layers. Some of these particles
    !may also not have enough neighbours for the Vandermonde matrix to be
    !invertible. These particles will be skipped without raising a warning.

    IF (with_ghosts) THEN
        np_target = Pc%Mpart
    ELSE
        np_target = Pc%Npart
    ENDIF

    lda = op%desc%nterms

    IF (.NOT. Pc%neighs%exists(op%neigh_id)) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            & 'Neighbour lists have not been created',&
            __LINE__,info)
        GOTO 9999
    ENDIF

    Nlist => Pc%neighs%vec(op%neigh_id)%t
    IF (.NOT. Nlist%uptodate) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,&
            'Neighbour lists are not up to date',__LINE__,info)
        GOTO 9999
    ENDIF
    nvlist => Nlist%nvlist
    vlist => Nlist%vlist


    SELECT TYPE(Pc)
    TYPE IS (DTYPE(ppm_t_sop))
        IF (isinterp) THEN
            Pc2 => DTYPE(ppm_Particles)%vec(Pc%set_aPc%vec(op%P_id))
            IF (.NOT. Pc2%props%exists(from_id)) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,caller,   &
                    & 'The operator input is not allocated.',&
                    __LINE__,info)
                GOTO 9999
            ELSE
                prop_from => Pc2%props%vec(from_id)%t
            ENDIF
        ELSE
            IF (.NOT. Pc%props%exists(from_id)) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,caller,   &
                    & 'The operator input is not allocated.',&
                    __LINE__,info)
                GOTO 9999
            ELSE
                prop_from => Pc%props%vec(from_id)%t
            ENDIF
        ENDIF
    CLASS DEFAULT
        IF (.NOT. Pc%props%exists(from_id)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,   &
                & 'The operator input is not allocated.',&
                __LINE__,info)
            GOTO 9999
        ELSE
            prop_from => Pc%props%vec(from_id)%t
        ENDIF
    END SELECT

    IF (.NOT.prop_from%flags(ppm_ppt_ghosts)) THEN
        WRITE(cbuf,*) 'Ghost values of ',TRIM(ADJUSTL(&
            prop_from%name)),' are needed.'
        CALL ppm_write(ppm_rank,caller,cbuf,info)
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,&
            'Please call particles_mapping_ghosts first',&
            __LINE__,info)
        GOTO 9999
    ENDIF

    IF (vector_output .AND. prop_from%lda .NE. lda) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            & 'Incompatible dimensions between operator and input data',&
            __LINE__,info)
        GOTO 9999
    ENDIF
    vector_input = (prop_from%lda .GE.2)

    !allocate output field if needed
    !otherwise simply check that the output array had been allocated
    !to the right size
    IF (to_id.EQ.0) THEN
        IF (vector_output) THEN
            CALL Pc%create_prop(to_id,ppm_type_real,info,lda=lda,&  
                name="dflt_dcop_apply",with_ghosts=with_ghosts)
        ELSE
            CALL Pc%create_prop(to_id,ppm_type_real,info,&
                name="dflt_dcop_apply",with_ghosts=with_ghosts)
        ENDIF
    ELSE
        ASSOCIATE (prop_to => Pc%props%vec(to_id)%t)
        !Destroy and reallocate the target property data structure
        ! if its type/dimension do not match that of the operator
        IF (      vector_output.AND.prop_to%lda.LT.2 .OR. &
             .NOT.vector_output.AND.prop_to%lda.NE.1 .OR. &
             prop_to%data_type.NE.ppm_type_real) THEN 
                CALL Pc%realloc_prop(to_id,info,with_ghosts=with_ghosts,&
                    datatype=ppm_type_real,lda=lda)
        ENDIF
        !Resize the target property array if its size does not match
        !that of the operators output.
        IF (.NOT.Pc%props%vec(to_id)%t%flags(ppm_ppt_partial).OR. &
            &  with_ghosts .AND. &
            &  .NOT.Pc%props%vec(to_id)%t%flags(ppm_ppt_ghosts)) THEN
            CALL Pc%realloc_prop(to_id,info,with_ghosts=with_ghosts)
        ENDIF
        END ASSOCIATE
    ENDIF
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,&
            'ppm_prop_(re)allocate failed',__LINE__,info)
        GOTO 9999
    ENDIF

    !zero the output array
    IF (vector_output) THEN
        CALL Pc%get(dwpv,to_id,with_ghosts=with_ghosts)
        DO ip = 1,np_target
            dwpv(1:lda,ip) = 0._MK
        ENDDO
    ELSE
        CALL Pc%get(dwps,to_id,with_ghosts=with_ghosts)
        DO ip = 1,np_target
            dwps(ip) = 0._MK
        ENDDO
    ENDIF
    eta => Pc%get_dcop(op_id,with_ghosts=with_ghosts)


    IF (isinterp) THEN
        IF (vector_output) THEN
            IF(vector_input) THEN
                CALL Pc2%get(wpv2,from_id,with_ghosts=.TRUE.)
                DO ip = 1,np_target
                    DO ineigh = 1,nvlist(ip)
                        iq = vlist(ineigh,ip)
                        dwpv(1:lda,ip) = dwpv(1:lda,ip) + &
                            wpv2(1:lda,iq) * eta(1+(ineigh-1)*lda:ineigh*lda,ip)
                    ENDDO
                ENDDO
                CALL Pc2%set(wpv2,from_id,read_only=.TRUE.)
            ELSE
                CALL Pc2%get(wps2,from_id,with_ghosts=.TRUE.)
                DO ip = 1,np_target
                    DO ineigh = 1,nvlist(ip)
                        iq = vlist(ineigh,ip)
                        dwpv(1:lda,ip) = dwpv(1:lda,ip) + &
                            wps2(iq) * eta(1+(ineigh-1)*lda:ineigh*lda,ip)
                    ENDDO
                ENDDO
                CALL Pc2%set(wps2,from_id,read_only=.TRUE.)
            ENDIF
        ELSE
            CALL Pc2%get(wps2,from_id,with_ghosts=.TRUE.)
            DO ip = 1,np_target
                DO ineigh = 1,nvlist(ip)
                    iq = vlist(ineigh,ip)
                    dwps(ip) = dwps(ip) + wps2(iq) * eta(ineigh,ip)
                ENDDO
            ENDDO
            CALL Pc2%set(wps2,from_id,read_only=.TRUE.)
        ENDIF
    ELSE
        sig = -1._mk 
        IF (vector_output) THEN
            IF(vector_input) THEN
                CALL Pc%get(wpv1,from_id,with_ghosts=.TRUE.)
                DO ip = 1,np_target
                    DO ineigh = 1,nvlist(ip)
                        iq = vlist(ineigh,ip)
                        dwpv(1:lda,ip) = dwpv(1:lda,ip) + &
                            (wpv1(1:lda,iq) + sig*(wpv1(1:lda,ip)))* &
                            eta(1+(ineigh-1)*lda:ineigh*lda,ip)
                    ENDDO
                ENDDO
                CALL Pc%set(wpv1,from_id,read_only=.TRUE.)
            ELSE
                CALL Pc%get(wps1,from_id,with_ghosts=.TRUE.)
                DO ip = 1,np_target
                    DO ineigh = 1,nvlist(ip)
                        iq = vlist(ineigh,ip)
                        dwpv(1:lda,ip) = dwpv(1:lda,ip) + &
                            (wps1(iq) + sig*(wps1(ip)))* &
                            eta(1+(ineigh-1)*lda:ineigh*lda,ip)
                    ENDDO
                ENDDO
                CALL Pc%set(wps1,from_id,read_only=.TRUE.)
            ENDIF
        ELSE
            CALL Pc%get(wps1,from_id,with_ghosts=.TRUE.)
            DO ip = 1,np_target
                DO ineigh = 1,nvlist(ip)
                    iq = vlist(ineigh,ip)
                    dwps(ip) = dwps(ip) + &
                        (wps1(iq)+sig*(wps1(ip))) * eta(ineigh,ip)
                ENDDO
            ENDDO
            CALL Pc%set(wps1,from_id,read_only=.TRUE.)
        ENDIF
    ENDIF

    eta => Pc%set_dcop(op_id)
    IF (vector_output) THEN
        IF (with_ghosts) THEN
            !we assume that the ghosts are up-to-date even though
            !they clearly are not. we assume you know what you are
            !doing when using this option.
            CALL Pc%set(dwpv,to_id,ghosts_ok=.TRUE.)
        ELSE
            CALL Pc%set(dwpv,to_id)
        ENDIF
    ELSE
        IF (with_ghosts) THEN
            CALL Pc%set(dwps,to_id,ghosts_ok=.TRUE.)
        ELSE
            CALL Pc%set(dwps,to_id)
        ENDIF
    ENDIF
    nvlist => NULL()
    vlist => NULL()

    Pc%stats%nb_dc_apply = Pc%stats%nb_dc_apply + 1
#ifdef __MPI
    t2 = MPI_WTIME(info)
    Pc%stats%t_dc_apply = Pc%stats%t_dc_apply+(t2-t1)
#endif

    CALL substop(caller,t0,info)

    9999  CONTINUE ! jump here upon error


END SUBROUTINE DTYPE(part_op_apply)


#undef DEFINE_MK
#undef DTYPE


