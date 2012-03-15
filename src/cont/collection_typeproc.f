#define WRAP(a) a
!!----------------------------------------------------------------
!! Template for Procedures for container class
!!----------------------------------------------------------------
!BEGIN
FUNCTION __CONTAINER(begin)(this) RESULT (iterator)
    CLASS(CONTAINER), INTENT(INOUT)   :: this
    CLASS(VEC_TYPE),POINTER   :: iterator

    this%iter_id = this%min_id
    IF (this%min_id.GT.0) THEN
        SELECT TYPE (t => this%vec(this%min_id))
        TYPE IS (VEC_TYPE)
            iterator => t
            RETURN
        END SELECT
    ENDIF
    iterator => NULL()
END FUNCTION __CONTAINER(begin)
!NEXT
FUNCTION __CONTAINER(next)(this) RESULT (iterator)
    CLASS(CONTAINER), INTENT(INOUT)   :: this
    CLASS(VEC_TYPE),POINTER   :: iterator

    IF (this%iter_id.LT.this%max_id) THEN
        this%iter_id = this%iter_id + 1
        SELECT TYPE (t => this%vec(this%iter_id))
        TYPE IS (VEC_TYPE)
            iterator => t
            RETURN
        END SELECT
    ENDIF
    iterator => NULL()
END FUNCTION __CONTAINER(next)
!PREVIOUS
FUNCTION __CONTAINER(prev)(this) RESULT (iterator)
    CLASS(CONTAINER), INTENT(INOUT)   :: this
    CLASS(VEC_TYPE),POINTER   :: iterator

    IF (this%iter_id.GT.this%min_id) THEN
        this%iter_id = this%iter_id - 1
        SELECT TYPE (t => this%vec(this%iter_id))
        TYPE IS (VEC_TYPE)
            iterator => t
            RETURN
        END SELECT
    ENDIF
    iterator => NULL()
END FUNCTION __CONTAINER(prev)
!LAST
FUNCTION __CONTAINER(last)(this) RESULT (iterator)
    CLASS(CONTAINER), INTENT(INOUT)   :: this
    CLASS(VEC_TYPE),POINTER   :: iterator

    this%iter_id = this%max_id
    IF (this%max_id.GT.0) THEN
        SELECT TYPE (t => this%vec(this%max_id))
        TYPE IS (VEC_TYPE)
            iterator => t
            RETURN
        END SELECT
    ENDIF
    iterator => NULL()

END FUNCTION __CONTAINER(last)

!DESTROY CONTAINER
SUBROUTINE __CONTAINER(destroy)(this,info)
    CLASS(CONTAINER), INTENT(INOUT)   :: this
    CLASS(VEC_TYPE),POINTER :: p => NULL()
    INTEGER, INTENT(OUT) :: info

    info = 0
    p => this%begin()
    DO WHILE(ASSOCIATED(p))
        CALL p%destroy(info)
        p => this%next()
    ENDDO
    IF (ASSOCIATED(this%vec))  DEALLOCATE(this%vec)
    this%min_id = HUGE(1)
    this%max_id = 0
    this%nb = 0
    this%vec_size=0

END SUBROUTINE __CONTAINER(destroy)

!EXISTS
FUNCTION __CONTAINER(exists)(this,id,caller) RESULT(exists)
    !!!------------------------------------------------------------------------!
    !!! Check whether a neighbor list exists and can be accessed at this id
    !!!------------------------------------------------------------------------!
    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    CLASS(CONTAINER)                               :: this
    !!! Data structure containing the particles
    INTEGER,                       INTENT(IN   )   :: id
    !!! id where the data is stored
    LOGICAL                                        :: exists
    !!! Return status, on success 0.
    CHARACTER(LEN = *),OPTIONAL                    :: caller
    !!! Calling routine

#include "../cont/container_exists.inc"

END FUNCTION __CONTAINER(exists)

!PUSH
SUBROUTINE __CONTAINER(push)(this,element,info)
    !!! add an element into the collection
    CLASS(CONTAINER)                   :: this
    CLASS(WRAP(VEC_TYPE)_),POINTER            :: element
    INTEGER,               INTENT(OUT) :: info
    TYPE(VEC_TYPE),DIMENSION(:),POINTER :: vec_temp

    INTEGER                            :: i

    info = 0
    !add the element at the end of the array
    this%max_id = this%max_id + 1
    this%nb = this%nb + 1

    !if the array is full, double its size
    IF (this%max_id.GT.this%vec_size) THEN
        SELECT TYPE (t => this%vec)
        TYPE IS (VEC_TYPE)
            ALLOCATE(vec_temp(this%vec_size))
            FORALL (i=1:this%vec_size)
                vec_temp(i) = t(i)
            END FORALL
        END SELECT

        DEALLOCATE(this%vec)
        ALLOCATE(VEC_TYPE::this%vec(2*this%vec_size))

        SELECT TYPE (t => this%vec)
        TYPE IS (VEC_TYPE)
            FORALL (i=1:this%vec_size)
                t(i) = vec_temp(i)
            END FORALL
            DEALLOCATE(vec_temp)
            this%vec_size = 2 * this%vec_size
        END SELECT
    ENDIF

    SELECT TYPE (t => this%vec(this%max_id))
    TYPE IS (VEC_TYPE)
        SELECT TYPE (e => element)
        TYPE IS (VEC_TYPE)
            t = e
        END SELECT
    END SELECT

END SUBROUTINE __CONTAINER(push)
!REMOVE
SUBROUTINE __CONTAINER(remove)(this,id,info)
    !!! remove an element (identified by its id) from the collection
    CLASS(CONTAINER)                   :: this
    INTEGER,          INTENT(IN   )    :: id
    INTEGER,          INTENT(  OUT)    :: info

    !deallocate the element
    CALL this%vec(id)%destroy(info)
    !swap with the last non-empty element of the collection
    IF (this%max_id.GT.this%min_id) THEN
        SELECT TYPE(v => this%vec)
        TYPE IS (VEC_TYPE)
            v(id) = v(this%max_id) 
        CLASS DEFAULT
        END SELECT
    ENDIF
    this%nb = this%nb - 1
    this%max_id = this%max_id - 1

END SUBROUTINE __CONTAINER(remove)
#undef CONTAINER
#undef __CONTAINER
#undef VEC_TYPE
#undef WRAP
