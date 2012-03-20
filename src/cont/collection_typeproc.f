#define WRAP(a) a
!!----------------------------------------------------------------
!! Template for Procedures for container class
!!----------------------------------------------------------------
!BEGIN
FUNCTION __CONTAINER(begin)(this) RESULT (iterator)
    CLASS(CONTAINER), INTENT(INOUT)   :: this
    CLASS(VEC_TYPE),POINTER   :: iterator

    this%iter_id = this%min_id
    IF (this%iter_id.GT.0 .AND. this%iter_id.LE.this%vec_size) THEN
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

    IF (this%iter_id.GE.this%min_id .AND. this%iter_id.LT.this%max_id) THEN
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

    IF (this%iter_id.GT.this%min_id .AND. this%iter_id.LE.this%max_id) THEN
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
    IF (this%iter_id.GT.0 .AND. this%iter_id.LE.this%vec_size) THEN
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

    CHARACTER(LEN=ppm_char)            :: caller = "__CONTAINER(destroy)"
    REAL(KIND(1.D0))                   :: t0

    CALL substart(caller,t0,info)

    p => this%begin()
    DO WHILE(ASSOCIATED(p))
        CALL p%destroy(info)
        or_fail_dealloc("Could not deallocate collection element")
        p => this%next()
    ENDDO
    IF (ASSOCIATED(this%vec))  DEALLOCATE(this%vec,STAT=info)
    or_fail_dealloc("Could not deallocate collection array")
    this%iter_id = 0
    this%min_id = 0
    this%max_id = 0
    this%nb = 0
    this%vec_size=0

    CALL substart(caller,t0,info)
    9999 CONTINUE
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
SUBROUTINE __CONTAINER(push)(this,element,info,id)
    !!! add an element into the collection
    CLASS(CONTAINER)                   :: this
    CLASS(WRAP(VEC_TYPE)_),POINTER            :: element
    INTEGER,               INTENT(OUT) :: info
    INTEGER,OPTIONAL,      INTENT(OUT) :: id
    !!! index of the element in the collection

    TYPE(VEC_TYPE),DIMENSION(:),POINTER :: vec_temp
    INTEGER                            :: i
    CHARACTER(LEN=ppm_char)            :: caller = "__CONTAINER(push)"
    REAL(KIND(1.D0))                   :: t0

    CALL substart(caller,t0,info)

    !add the element at the end of the array
    this%min_id = 1
    this%max_id = this%max_id + 1
    this%nb = this%nb + 1
    IF (PRESENT(id)) id = this%max_id

    IF (this%max_id.GT.this%vec_size) THEN
        !if the array is empty, allocate with a reasonable size
        IF (this%vec_size.EQ.0 .OR. .NOT.ASSOCIATED(this%vec)) THEN
            this%vec_size = 10
            ALLOCATE(VEC_TYPE::this%vec(this%vec_size),STAT=info)
        ELSE
        !if the array is full, double its size
            SELECT TYPE (t => this%vec)
            TYPE IS (VEC_TYPE)
                ALLOCATE(vec_temp(this%vec_size))
                FORALL (i=1:this%vec_size)
                    vec_temp(i) = t(i)
                END FORALL
            END SELECT

            DEALLOCATE(this%vec)
            ALLOCATE(VEC_TYPE::this%vec(2*this%vec_size),STAT=info)
            or_fail_alloc("could not allocate collection array")

            SELECT TYPE (t => this%vec)
            TYPE IS (VEC_TYPE)
                FORALL (i=1:this%vec_size)
                    t(i) = vec_temp(i)
                END FORALL
                DEALLOCATE(vec_temp)
                this%vec_size = 2 * this%vec_size
            END SELECT
        ENDIF
    ENDIF

    SELECT TYPE (t => this%vec(this%max_id))
    TYPE IS (VEC_TYPE)
        SELECT TYPE (e => element)
        TYPE IS (VEC_TYPE)
            t = e
        END SELECT
    END SELECT

    CALL substart(caller,t0,info)
    9999 CONTINUE

END SUBROUTINE __CONTAINER(push)
!REMOVE
SUBROUTINE __CONTAINER(remove)(this,id,info)
    !!! remove an element (identified by its id) from the collection
    CLASS(CONTAINER)                   :: this
    INTEGER,          INTENT(IN   )    :: id
    INTEGER,          INTENT(  OUT)    :: info

    !deallocate the element
    !CALL this%vec(id)%destroy(info)
    SELECT TYPE (t => this%vec(id))
    TYPE IS (VEC_TYPE)
        CALL t%destroy(info)
    END SELECT
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
    IF (this%nb.EQ.0 .OR. this%max_id.EQ.0) this%min_id = 0

END SUBROUTINE __CONTAINER(remove)
#undef CONTAINER
#undef __CONTAINER
#undef VEC_TYPE
#undef WRAP
