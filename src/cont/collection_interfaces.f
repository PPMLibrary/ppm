!! INTERFACES for container class
!BEGIN
FUNCTION __CONTAINER(begin)(this) RESULT (iterator)
    IMPORT CONTAINER,VEC_TYPE
    CLASS(CONTAINER), INTENT(INOUT)   :: this
    CLASS(VEC_TYPE),POINTER   :: iterator
END FUNCTION
!NEXT
FUNCTION __CONTAINER(next)(this) RESULT (iterator)
    IMPORT CONTAINER,VEC_TYPE
    CLASS(CONTAINER), INTENT(INOUT)   :: this
    CLASS(VEC_TYPE),POINTER   :: iterator
END FUNCTION
!PREVIOUS
FUNCTION __CONTAINER(prev)(this) RESULT (iterator)
    IMPORT CONTAINER,VEC_TYPE
    CLASS(CONTAINER), INTENT(INOUT)   :: this
    CLASS(VEC_TYPE),POINTER   :: iterator
END FUNCTION
!LAST
FUNCTION __CONTAINER(last)(this) RESULT (iterator)
    IMPORT CONTAINER,VEC_TYPE
    CLASS(CONTAINER), INTENT(INOUT)   :: this
    CLASS(VEC_TYPE),POINTER   :: iterator
END FUNCTION
!DESTROY CONTAINER
SUBROUTINE __CONTAINER(destroy)(this,info)
    IMPORT CONTAINER,VEC_TYPE
    CLASS(CONTAINER), INTENT(INOUT)   :: this
    CLASS(VEC_TYPE),POINTER :: p => NULL()
    INTEGER, INTENT(OUT) :: info
END SUBROUTINE
!EXISTS
FUNCTION __CONTAINER(exists)(this,id,caller) RESULT(exists)
    IMPORT CONTAINER
    CLASS(CONTAINER)                               :: this
    INTEGER,                       INTENT(IN   )   :: id
    LOGICAL                                        :: exists
    CHARACTER(LEN = *),OPTIONAL                    :: caller
END FUNCTION
!PUSH
SUBROUTINE __CONTAINER(push)(this,element,info,id)
    IMPORT CONTAINER,VEC_TYPE
    CLASS(CONTAINER)                   :: this
    CLASS(VEC_TYPE),POINTER            :: element
    INTEGER,               INTENT(OUT) :: info
    INTEGER,OPTIONAL,      INTENT(OUT) :: id
END SUBROUTINE
!REMOVE
SUBROUTINE __CONTAINER(remove)(this,id,info)
    IMPORT CONTAINER
    CLASS(CONTAINER)                   :: this
    INTEGER,          INTENT(IN   )    :: id
    INTEGER,          INTENT(  OUT)    :: info
END SUBROUTINE

#undef CONTAINER
#undef __CONTAINER
#undef VEC_TYPE
