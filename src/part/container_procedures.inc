!BEGIN
FUNCTION DTYPE(begin_op)(this) RESULT (iterator)
    CLASS(DTYPE(ppm_c_operators)),  INTENT(INOUT)   :: this
    CLASS(DTYPE(ppm_t_operator)),   POINTER      :: iterator

    this%iter_id = this%min_id
    IF (this%min_id.NE.0) THEN
        iterator => this%vec(this%min_id)%t
    ELSE
        iterator => NULL()
    ENDIF

END FUNCTION DTYPE(begin_op)

FUNCTION DTYPE(begin_neigh)(this) RESULT (iterator)
    CLASS(DTYPE(ppm_c_neighlists)),  INTENT(INOUT)   :: this
    CLASS(DTYPE(ppm_t_neighlist)),   POINTER      :: iterator

    this%iter_id = this%min_id
    IF (this%min_id.NE.0) THEN
        iterator => this%vec(this%min_id)%t
    ELSE
        iterator => NULL()
    ENDIF

END FUNCTION DTYPE(begin_neigh)

FUNCTION DTYPE(begin_prop)(this) RESULT (iterator)
    CLASS(DTYPE(ppm_c_props)), INTENT(INOUT)   :: this
    CLASS(DTYPE(ppm_t_part_prop)),POINTER   :: iterator

    this%iter_id = this%min_id
    IF (this%min_id.NE.0) THEN
        iterator => this%vec(this%min_id)%t
    ELSE
        iterator => NULL()
    ENDIF

END FUNCTION DTYPE(begin_prop)

!LAST
FUNCTION DTYPE(last_op)(this) RESULT (iterator)
    CLASS(DTYPE(ppm_c_operators)),  INTENT(INOUT)   :: this
    CLASS(DTYPE(ppm_t_operator)),   POINTER      :: iterator

    this%iter_id = this%max_id
    IF (this%max_id.NE.0) THEN
        iterator => this%vec(this%max_id)%t
    ELSE
        iterator => NULL()
    ENDIF

END FUNCTION DTYPE(last_op)

FUNCTION DTYPE(last_neigh)(this) RESULT (iterator)
    CLASS(DTYPE(ppm_c_neighlists)),  INTENT(INOUT)   :: this
    CLASS(DTYPE(ppm_t_neighlist)),   POINTER      :: iterator

    this%iter_id = this%max_id
    IF (this%max_id.NE.0) THEN
        iterator => this%vec(this%max_id)%t
    ELSE
        iterator => NULL()
    ENDIF

END FUNCTION DTYPE(last_neigh)

FUNCTION DTYPE(last_prop)(this) RESULT (iterator)
    CLASS(DTYPE(ppm_c_props)), INTENT(INOUT)   :: this
    CLASS(DTYPE(ppm_t_part_prop)),POINTER   :: iterator

    this%iter_id = this%max_id
    IF (this%max_id.NE.0) THEN
        iterator => this%vec(this%max_id)%t
    ELSE
        iterator => NULL()
    ENDIF

END FUNCTION DTYPE(last_prop)


!NEXT
FUNCTION DTYPE(next_op)(this) RESULT (iterator)
    CLASS(DTYPE(ppm_c_operators)),  INTENT(INOUT)   :: this
    CLASS(DTYPE(ppm_t_operator)),   POINTER      :: iterator

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

FUNCTION DTYPE(next_neigh)(this) RESULT (iterator)
    CLASS(DTYPE(ppm_c_neighlists)),  INTENT(INOUT):: this
    CLASS(DTYPE(ppm_t_neighlist)),   POINTER      :: iterator

    DO WHILE(this%iter_id.LT.this%max_id)
        this%iter_id = this%iter_id + 1
        IF (ASSOCIATED(this%vec(this%iter_id)%t)) THEN 
            iterator => this%vec(this%iter_id)%t
            RETURN
        ENDIF
    ENDDO
    iterator => NULL()
    RETURN

END FUNCTION DTYPE(next_neigh)

FUNCTION DTYPE(next_prop)(this) RESULT (iterator)
    CLASS(DTYPE(ppm_c_props)), INTENT(INOUT)   :: this
    CLASS(DTYPE(ppm_t_part_prop)),POINTER   :: iterator

    DO WHILE(this%iter_id.LT.this%max_id)
        this%iter_id = this%iter_id + 1
        IF (ASSOCIATED(this%vec(this%iter_id)%t)) THEN 
            iterator => this%vec(this%iter_id)%t
            RETURN
        ENDIF
    ENDDO
    iterator => NULL()
    RETURN

END FUNCTION DTYPE(next_prop)


!PREVIOUS
FUNCTION DTYPE(prev_op)(this) RESULT (iterator)
    CLASS(DTYPE(ppm_c_operators)),  INTENT(INOUT)   :: this
    CLASS(DTYPE(ppm_t_operator)),   POINTER      :: iterator

    DO WHILE(this%iter_id.GE.this%min_id)
        this%iter_id = this%iter_id - 1
        IF (ASSOCIATED(this%vec(this%iter_id)%t)) THEN 
            iterator => this%vec(this%iter_id)%t
            RETURN
        ENDIF
    ENDDO
    iterator => NULL()
    RETURN
       
END FUNCTION DTYPE(prev_op)

FUNCTION DTYPE(prev_neigh)(this) RESULT (iterator)
    CLASS(DTYPE(ppm_c_neighlists)),  INTENT(INOUT)   :: this
    CLASS(DTYPE(ppm_t_neighlist)),   POINTER      :: iterator

    DO WHILE(this%iter_id.GE.this%min_id)
        this%iter_id = this%iter_id - 1
        IF (ASSOCIATED(this%vec(this%iter_id)%t)) THEN 
            iterator => this%vec(this%iter_id)%t
            RETURN
        ENDIF
    ENDDO
    iterator => NULL()
    RETURN

END FUNCTION DTYPE(prev_neigh)

FUNCTION DTYPE(prev_prop)(this) RESULT (iterator)
    CLASS(DTYPE(ppm_c_props)), INTENT(INOUT)   :: this
    CLASS(DTYPE(ppm_t_part_prop)),POINTER   :: iterator

    DO WHILE(this%iter_id.GE.this%min_id)
        this%iter_id = this%iter_id - 1
        IF (ASSOCIATED(this%vec(this%iter_id)%t)) THEN 
            iterator => this%vec(this%iter_id)%t
            RETURN
        ENDIF
    ENDDO
    iterator => NULL()
    RETURN

END FUNCTION DTYPE(prev_prop)
