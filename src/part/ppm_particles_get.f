#define WRAP(a) a


#define __FUNCNAME DTYPE(WRAP(DATANAME)_check)
SUBROUTINE __FUNCNAME(this,wp,info) 
    !!!------------------------------------------------------------------------!
    !!! Check whether a Data Structure exists and can be accessed at this id
    !!!------------------------------------------------------------------------!

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    CLASS(DTYPE(ppm_t_part_prop))                       :: this
    !!! Data structure containing the particles
#if   __DIM == 1
    __TYPE,DIMENSION(:),POINTER,        INTENT(IN   )   :: wp
#elif __DIM == 2
    __TYPE,DIMENSION(:,:),POINTER,      INTENT(IN   )   :: wp
#endif
    INTEGER,                            INTENT(   OUT)  :: info
    !!! Return status, on success 0.
    INTEGER, DIMENSION(:),POINTER :: nullv=>NULL()
    
    start_subroutine(__FUNCNAME)
    !-------------------------------------------------------------------------
    ! Check arguments
    !-------------------------------------------------------------------------

#if   __DIM == 1
    IF (this%lda.NE.1) THEN
#elif __DIM == 2
    IF (this%lda.LT.2) THEN
#endif
        fail ("Argument has wrong dimension for this property")
    ENDIF

    IF (this%data_type.NE. &
#if   __MYTYPE == __INTEGER
        ppm_type_int&
#elif __MYTYPE == __LONGINT
        ppm_type_longint& 
#elif __MYTYPE == __REAL
        ppm_type_real& 
#elif __MYTYPE == __COMPLEX
        ppm_type_comp& 
#elif __MYTYPE == __LOGICAL
        ppm_type_logical& 
#else
        ppm_type_none&
#endif
       &  ) THEN

       fail ("Argument has wrong dimension for this property")
    ENDIF

    end_subroutine()

END SUBROUTINE __FUNCNAME
#undef __FUNCNAME

#define __FUNCNAME DTYPE(WRAP(DATANAME)_get_prop)
SUBROUTINE __FUNCNAME(this,discr_data,wp,info,with_ghosts,read_only)
    CLASS(DTYPE(ppm_t_particles))   :: this
    !CLASS(DTYPE(ppm_t_part_prop)_),POINTER  :: discr_data
    CLASS(ppm_t_discr_data)          :: discr_data
#if   __DIM == 1
    __TYPE,DIMENSION(:),POINTER     :: wp
#elif __DIM == 2
    __TYPE,DIMENSION(:,:),POINTER   :: wp
#endif
    INTEGER                         :: info
    !!! Return status, on success 0.
    LOGICAL,OPTIONAL                :: with_ghosts
    !!! returns array between 1:Mpart (default is 1:Npart)
    LOGICAL,OPTIONAL                :: read_only
    !!! swear on your favourite book that you will not modify the data
    !!! that you access. Default is false (which implies that the state
    !!! variables will be modify with the assumption that the data
    !!! accessed has been changed, e.g. so that a subsequent ghost update
    !!! will also update this property)

    INTEGER                         :: np

    start_subroutine(__FUNCNAME)

    wp => NULL()

    SELECT TYPE(discr_data)
    CLASS IS (DTYPE(ppm_t_part_prop)_)

    IF (ppm_debug.GE.1) THEN
        CALL discr_data%checktype(wp,info)
        or_fail("Argument has the wrong type for this property")
    ENDIF

    np = this%Npart
    IF (PRESENT(with_ghosts)) THEN
        IF (with_ghosts) THEN
            IF (discr_data%flags(ppm_ppt_ghosts)) THEN
                np = this%Mpart
            ELSE
                WRITE(cbuf,*)"ERROR: tried to get DATANAME (name = ",&
                    & TRIM(ADJUSTL(discr_data%name)),&
                    & ") with ghosts when ghosts are not up-to-date. ",&
                    & "Returning NULL pointer"
                CALL ppm_write(ppm_rank,caller,cbuf,info)
                    fail("ghosts are not up-to-date")
            ENDIF
        ENDIF
    ENDIF

    IF (discr_data%flags(ppm_ppt_partial)) THEN
        wp => &
#if   __DIM == 1
            discr_data%WRAP(DATANAME)(1:np)
#elif __DIM == 2
            discr_data%WRAP(DATANAME)(:,1:np)
#endif
    ELSE
        WRITE(cbuf,*) 'ERROR: tried to get DATANAME (name = ',&
            & TRIM(ADJUSTL(discr_data%name)),&
            & ') when mapping is not up-to-date. ',&
            & 'Returning NULL pointer', &
            'Run with traceback option to debug'
        CALL ppm_write(ppm_rank,cbuf,caller,info)
            fail("unmapped particles")
    ENDIF

    !Assume that the ghost values are now incorrect, unless explicitely
    !told otherwise
    IF (PRESENT(read_only)) THEN
        IF (.NOT.read_only) THEN
            discr_data%flags(ppm_ppt_ghosts) = .FALSE.
        ENDIF
    ELSE
        discr_data%flags(ppm_ppt_ghosts) = .FALSE.
    ENDIF

    END SELECT

    check_associated(wp,"Get_Prop returned a NULL pointer")

    end_subroutine()
END SUBROUTINE __FUNCNAME
#undef __FUNCNAME

#define __FUNCNAME DTYPE(WRAP(DATANAME)_set_prop)
SUBROUTINE __FUNCNAME(this,discr_data,wp,info,read_only,ghosts_ok)
    CLASS(DTYPE(ppm_t_particles))    :: this
    CLASS(DTYPE(ppm_t_part_prop)_)   :: discr_data
    INTEGER                          :: info
    !!! Return status, on success 0.
    LOGICAL,OPTIONAL                 :: read_only
    LOGICAL,OPTIONAL                 :: ghosts_ok
#if   __DIM == 1
    __TYPE,DIMENSION(:),POINTER      :: wp
#elif __DIM == 2
    __TYPE,DIMENSION(:,:),POINTER    :: wp
#endif

    start_subroutine(__FUNCNAME)


    !If read_only was not explicitely set to true, then assume
    !that ghosts are no longer up to date, unless ghosts_ok was
    ! explicitely set to true
    IF (PRESENT(ghosts_ok)) THEN
        IF (.NOT.ghosts_ok) THEN
            !Assume that the ghost values are now incorrect
            discr_data%flags(ppm_ppt_ghosts) = .FALSE.
        ENDIF
    ENDIF

    IF (PRESENT(read_only)) THEN
        IF (.NOT.read_only) THEN
            !Assume that the ghost values are now incorrect
            discr_data%flags(ppm_ppt_ghosts) = .FALSE.
        ENDIF
    ENDIF

    wp => NULL()

    end_subroutine()
END SUBROUTINE __FUNCNAME
#undef __FUNCNAME

#define __FUNCNAME DTYPE(WRAP(DATANAME)_get_field)
SUBROUTINE __FUNCNAME(this,Field,wp,info,with_ghosts,read_only)
    !!! Returns a pointer to the data array where that contains
    !!! the discretized elements of Field on this particle set.
    CLASS(DTYPE(ppm_t_particles))   :: this
    CLASS(ppm_t_field_)             :: Field
#if   __DIM == 1
    __TYPE,DIMENSION(:),POINTER     :: wp
    !!! data array 
#elif __DIM == 2
    __TYPE,DIMENSION(:,:),POINTER   :: wp
    !!! data array
#endif
    INTEGER                         :: info
    !!! Return status, on success 0.
    LOGICAL,OPTIONAL                :: with_ghosts
    !!! returns array between 1:Mpart (default is 1:Npart)
    LOGICAL,OPTIONAL                :: read_only
    !!! swear on your favourite book that you will not modify the data
    !!! that you access. Default is false (which implies that the state
    !!! variables will be modify with the assumption that the data
    !!! accessed has been changed, e.g. so that a subsequent ghost update
    !!! will also update this property)

    INTEGER                         :: np
    CLASS(ppm_t_discr_data),      POINTER :: discr_data => NULL()

    start_subroutine(__FUNCNAME)

    wp => NULL()

    CALL Field%get_discr(this,discr_data,info)
        or_fail("could not get discr data for this field on that particle set")

    SELECT TYPE(prop => discr_data)
    CLASS IS (DTYPE(ppm_t_part_prop))
        IF (ppm_debug.GE.1) THEN
            CALL prop%checktype(wp,info)
            or_fail("Argument has the wrong type for this property")
        ENDIF

        np = this%Npart
        IF (PRESENT(with_ghosts)) THEN
            IF (with_ghosts) THEN
                IF (prop%flags(ppm_ppt_ghosts)) THEN
                    np = this%Mpart
                ELSE
                    WRITE(cbuf,*)"ERROR: tried to get DATANAME (name = ",&
                        & TRIM(ADJUSTL(prop%name)),&
                        & ") with ghosts when ghosts are not up-to-date. ",&
                        & "Returning NULL pointer"
                    CALL ppm_write(ppm_rank,caller,cbuf,info)
                    fail("ghosts are not up-to-date")
                ENDIF
            ENDIF
        ENDIF

        IF (prop%flags(ppm_ppt_partial)) THEN
            wp => &
#if   __DIM == 1
                prop%WRAP(DATANAME)(1:np)
#elif __DIM == 2
                prop%WRAP(DATANAME)(:,1:np)
#endif
        ELSE
            WRITE(cbuf,*) 'ERROR: tried to get DATANAME (name = ',&
                & TRIM(ADJUSTL(prop%name)),&
                & ') when mapping is not up-to-date. ',&
                & 'Returning NULL pointer', &
                'Run with traceback option to debug'
            CALL ppm_write(ppm_rank,cbuf,caller,info)
            fail("unmapped particles")
        ENDIF
        !Assume that the ghost values are now incorrect, unless explicitely
        !told otherwise
        IF (PRESENT(read_only)) THEN
            IF (.NOT.read_only) THEN
                prop%flags(ppm_ppt_ghosts) = .FALSE.
            ENDIF
        ELSE
            prop%flags(ppm_ppt_ghosts) = .FALSE.
        ENDIF
    CLASS DEFAULT
        fail("wrong type. Discretized data should be a ppm_t_part_prop")
    END SELECT

    check_associated(wp,"Get_Field returned a NULL pointer")

    end_subroutine()
END SUBROUTINE __FUNCNAME
#undef __FUNCNAME

#define __FUNCNAME DTYPE(WRAP(DATANAME)_set_field)
SUBROUTINE __FUNCNAME(this,Field,wp,info,read_only,ghosts_ok)
    CLASS(DTYPE(ppm_t_particles))    :: this
    CLASS(ppm_t_field_)              :: Field
    INTEGER                          :: info
    !!! Return status, on success 0.
    LOGICAL,OPTIONAL                 :: read_only
    LOGICAL,OPTIONAL                 :: ghosts_ok
#if   __DIM == 1
    __TYPE,DIMENSION(:),POINTER      :: wp
#elif __DIM == 2
    __TYPE,DIMENSION(:,:),POINTER    :: wp
#endif

    CLASS(ppm_t_discr_data),      POINTER :: discr_data => NULL()

    start_subroutine(__FUNCNAME)

    CALL Field%get_discr(this,discr_data,info)
        or_fail("could not get discr data for this field on that particle set")

    SELECT TYPE(prop => discr_data)
    CLASS IS (DTYPE(ppm_t_part_prop))

    !If read_only was not explicitely set to true, then assume
    !that ghosts are no longer up to date, unless ghosts_ok was
    ! explicitely set to true
    IF (PRESENT(ghosts_ok)) THEN
        IF (.NOT.ghosts_ok) THEN
            !Assume that the ghost values are now incorrect
            prop%flags(ppm_ppt_ghosts) = .FALSE.
        ENDIF
    ENDIF

    IF (PRESENT(read_only)) THEN
        IF (.NOT.read_only) THEN
            !Assume that the ghost values are now incorrect
            prop%flags(ppm_ppt_ghosts) = .FALSE.
        ENDIF
    ENDIF
    END SELECT

    wp => NULL()

    end_subroutine()
END SUBROUTINE __FUNCNAME
#undef __FUNCNAME

#define __FUNCNAME DTYPE(WRAP(DATANAME)_get)
#define __CHECKTYPE DTYPE(WRAP(DATANAME)_check)
SUBROUTINE __FUNCNAME(Pc,wp,ppt_id,with_ghosts,read_only)
    CLASS(DTYPE(ppm_t_particles))   :: Pc
    INTEGER                         :: ppt_id
#if   __DIM == 1
    __TYPE,DIMENSION(:),POINTER     :: wp
#elif __DIM == 2
    __TYPE,DIMENSION(:,:),POINTER   :: wp
#endif
    INTEGER                         :: info
    !!! Return status, on success 0.
    LOGICAL,OPTIONAL                :: with_ghosts
    !!! returns array between 1:Mpart (default is 1:Npart)
    LOGICAL,OPTIONAL                :: read_only
    !!! swear on your favourite book that you will not modify the data
    !!! that you access. Default is false (which implies that the state
    !!! variables will be modify with the assumption that the data
    !!! accessed has been changed, e.g. so that a subsequent ghost update
    !!! will also update this property)
    LOGICAL   :: lghosts


    start_subroutine(__FUNCNAME)

    wp => NULL()
    lghosts = .FALSE.
    IF (PRESENT(with_ghosts)) THEN
        IF (with_ghosts) lghosts = .TRUE.
    ENDIF

    IF (ppt_id .LE. 0) THEN
        stdout("ERROR: failed to get DATANAME for property with ppt_id = ",ppt_id)
        fail("Cannot get property. Returning Null pointer")
    ENDIF


    IF (ppt_id .LE. Pc%props%max_id) THEN
        ASSOCIATE (prop => Pc%props%vec(ppt_id)%t)
        IF (prop%flags(ppm_ppt_partial)) THEN
            IF (lghosts) THEN
                IF (prop%flags(ppm_ppt_ghosts)) THEN
                    wp => &
#if   __DIM == 1
                        prop%WRAP(DATANAME)(1:Pc%Mpart)
#elif __DIM == 2
                        prop%WRAP(DATANAME)(:,1:Pc%Mpart)
#endif
                ELSE
                    stdout("ERROR: tried to get DATANAME (name = ",&
                        & 'TRIM(ADJUSTL(prop%name))',&
                        & ") with ghosts when ghosts are not up-to-date. ",&
                        & "Returning NULL pointer")
                    fail("Ghosts not up-to-date. Call map_ghosts()?")
                ENDIF
            ELSE
                wp => &
#if   __DIM == 1
                    prop%WRAP(DATANAME)(1:Pc%Npart)
#elif __DIM == 2
                    prop%WRAP(DATANAME)(:,1:Pc%Npart)
#endif
            ENDIF
        ELSE
            stdout("ERROR: tried to get DATANAME (name = ",&
                & 'TRIM(ADJUSTL(Pc%props%vec(ppt_id)%t%name))',&
                & ") when mapping is not up-to-date. ",&
                & "Returning NULL pointer")
            fail("unmapped particles")
        ENDIF

        IF (PRESENT(read_only)) THEN
            IF (.NOT.read_only) prop%flags(ppm_ppt_ghosts) = .FALSE.
        ELSE
            prop%flags(ppm_ppt_ghosts) = .FALSE.
        ENDIF

        END ASSOCIATE
    ELSE
        fail("Invalid id for particle property, returning NULL pointer")
    ENDIF

    end_subroutine()
END SUBROUTINE __FUNCNAME
#undef __FUNCNAME


#define __FUNCNAME DTYPE(WRAP(DATANAME)_set)
SUBROUTINE __FUNCNAME(Pc,wp,ppt_id,read_only,ghosts_ok)
    CLASS(DTYPE(ppm_t_particles))    :: Pc
    INTEGER                          :: ppt_id
    LOGICAL,OPTIONAL                 :: read_only
    LOGICAL,OPTIONAL                 :: ghosts_ok
#if   __DIM == 1
    __TYPE,DIMENSION(:),POINTER      :: wp
#elif __DIM == 2
    __TYPE,DIMENSION(:,:),POINTER    :: wp
#endif

    !If read_only was not explicitely set to true, then assume
    !that ghosts are no longer up to date, unless ghosts_ok was
    ! explicitely set to true
    IF (PRESENT(ghosts_ok)) THEN
        IF (ghosts_ok) THEN
            wp => NULL()
            RETURN
        ENDIF
    ENDIF

    IF (PRESENT(read_only)) THEN
        IF (read_only) THEN
            wp => NULL()
            RETURN
        ENDIF
    ENDIF

    !Assume that the ghost values are now incorrect
    Pc%props%vec(ppt_id)%t%flags(ppm_ppt_ghosts) = .FALSE.
    wp => NULL()

END SUBROUTINE __FUNCNAME
#undef __FUNCNAME

#undef DATANAME
#undef __TYPE
#undef __MYTYPE
