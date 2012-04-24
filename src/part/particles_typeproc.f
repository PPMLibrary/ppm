minclude define_collection_procedures(DTYPE(ppm_t_part_prop))
minclude define_collection_procedures(DTYPE(ppm_t_neighlist))
minclude define_collection_procedures(DTYPE(ppm_t_particles))
!minclude define_collection_procedures(DTYPE(ppm_t_sop))

SUBROUTINE DTYPE(prop_create)(prop,datatype,npart,lda,name,flags,info,field,zero)
    !!! Constructor for particle property data structure
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_part_prop))      :: prop
    INTEGER,                INTENT(IN) :: datatype
    INTEGER,                INTENT(IN) :: npart
    INTEGER,                INTENT(IN) :: lda
    CHARACTER(LEN=*),       INTENT(IN) :: name
    !!! name to this property
    LOGICAL, DIMENSION(ppm_param_length_pptflags),INTENT(IN) :: flags
    INTEGER,               INTENT(OUT) :: info
    CLASS(ppm_t_field_),OPTIONAL,TARGET, INTENT(IN) :: field
    LOGICAL, OPTIONAL,     INTENT( IN) :: zero
    !!! if true, then initialize the data to zero

    INTEGER                            :: iopt
    LOGICAL                            :: is2d
    LOGICAL                            :: zero_data

    start_subroutine("prop_create")

    prop%lda       = lda
    prop%data_type = datatype
    IF (PRESENT(field)) THEN
        prop%field_ptr => field
    ENDIF
    prop%name      = name
    prop%flags     = flags

    IF (PRESENT(zero)) THEN
        zero_data = zero
    ELSE
        zero_data = .FALSE.
    ENDIF


    iopt   = ppm_param_alloc_grow

    IF (lda.GE.2) THEN
        ldc(1) = lda
        ldc(2) = npart
        is2d = .TRUE.
    ELSE
        ldc(1) = npart
        is2d = .FALSE.
    ENDIF


    IF (is2d) THEN
        SELECT CASE (datatype)
        CASE (ppm_type_int)
            CALL ppm_alloc(prop%data_2d_i,ldc,iopt,info)
            IF (zero_data) prop%data_2d_i(1:lda,1:npart) = 0
        CASE (ppm_type_longint)
            CALL ppm_alloc(prop%data_2d_li,ldc,iopt,info)
            IF (zero_data) prop%data_2d_li(1:lda,1:npart) = 0
        CASE (ppm_type_real)
            CALL ppm_alloc(prop%data_2d_r,ldc,iopt,info)
            IF (zero_data) prop%data_2d_r(1:lda,1:npart) = 0._MK
        CASE (ppm_type_comp)
            CALL ppm_alloc(prop%data_2d_c,ldc,iopt,info)
            IF (zero_data) prop%data_2d_c(1:lda,1:npart) = 0._MK
        CASE (ppm_type_logical )
            CALL ppm_alloc(prop%data_2d_l,ldc,iopt,info)
            IF (zero_data) prop%data_2d_l(1:lda,1:npart) = .FALSE.
        CASE DEFAULT
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_argument,caller,   &
                &        'invalid type for particle property',__LINE__,info)
        END SELECT
    ELSE
        SELECT CASE (datatype)
        CASE (ppm_type_int)
            CALL ppm_alloc(prop%data_1d_i,ldc,iopt,info)
            IF (zero_data) prop%data_1d_i(1:npart) = 0
        CASE (ppm_type_longint)
            CALL ppm_alloc(prop%data_1d_li,ldc,iopt,info)
            IF (zero_data) prop%data_1d_li(1:npart) = 0
        CASE (ppm_type_real)
            CALL ppm_alloc(prop%data_1d_r,ldc,iopt,info)
            IF (zero_data) prop%data_1d_r(1:npart) = 0._MK
        CASE (ppm_type_comp)
            CALL ppm_alloc(prop%data_1d_c,ldc,iopt,info)
            IF (zero_data) prop%data_1d_c(1:npart) = 0._MK
        CASE (ppm_type_logical )
            CALL ppm_alloc(prop%data_1d_l,ldc,iopt,info)
            IF (zero_data) prop%data_1d_l(1:npart) = .FALSE.
        CASE DEFAULT
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_argument,caller,   &
                &        'invalid type for particle property',__LINE__,info)
        END SELECT
    ENDIF

    or_fail_dealloc("allocating property failed")

    end_subroutine()
END SUBROUTINE DTYPE(prop_create)
!DESTROY ENTRY
SUBROUTINE DTYPE(prop_destroy)(prop,info)
    CLASS(DTYPE(ppm_t_part_prop))      :: prop
    INTEGER,                                INTENT(  OUT)  :: info
    !!! Returns status, 0 upon success.

    start_subroutine("prop_destroy")

    dealloc_pointer(prop%data_1d_i)
    dealloc_pointer(prop%data_2d_i)
    dealloc_pointer(prop%data_1d_li)
    dealloc_pointer(prop%data_2d_li)
    dealloc_pointer(prop%data_1d_r)
    dealloc_pointer(prop%data_2d_r)
    dealloc_pointer(prop%data_1d_c)
    dealloc_pointer(prop%data_2d_c)
    dealloc_pointer(prop%data_1d_l)
    dealloc_pointer(prop%data_2d_l)

    prop%data_type = ppm_type_none
    prop%lda = 0
    prop%flags = .FALSE.
    prop%field_ptr => NULL()
    prop%name = ""

    end_subroutine()
END SUBROUTINE DTYPE(prop_destroy)

SUBROUTINE DTYPE(prop_print_info)(prop,info,level,fileunit,propid)
    !-----------------------------------------------------------------------
    ! Print out summary information about this property
    !-----------------------------------------------------------------------

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_part_prop))                          :: prop
    !!! Data structure containing the particles
    INTEGER,                            INTENT(  OUT)      :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    !  Optional arguments
    !-------------------------------------------------------------------------
    INTEGER, OPTIONAL,                  INTENT(IN   )      :: level
    !!! indentation level at which to printout the info. Default = 0
    INTEGER, OPTIONAL,                  INTENT(IN   )      :: fileunit
    !!! Already open file unit for printout. Default = stdout
    INTEGER, OPTIONAL,                  INTENT(IN   )      :: propid
    !!! id of this property in the parent struct
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    INTEGER                              :: lev,fileu,id
    CHARACTER(LEN = ppm_char)            :: myformat

    start_subroutine("prop_print_info")


    IF (PRESENT(fileunit)) THEN
        fileu = fileunit
    ELSE
        fileu = 6
    ENDIF
    IF (PRESENT(level)) THEN
        lev = MAX(level,1)
    ELSE
        lev = 1
    ENDIF
    IF (PRESENT(propid)) THEN
        id = propid
    ELSE
        id = 1
    ENDIF

    WRITE(myformat,'(A,I0,A)') '(',4*lev,'X,A,I0,A,A,A,I0)'

    WRITE(fileu,myformat) 'Property ',id,': ',&
        TRIM(color_print(prop%name,33)),&
        ' Type: ',prop%data_type

    lev = lev + 1

    WRITE(myformat,'(A,I0,A)') '(',4*lev,'X,A,I0)'
    WRITE(fileu,myformat) 'lda: ',prop%lda

    WRITE(myformat,'(A,I0,A,I0,A)') '(',4*lev,'X,A,',&
        ppm_param_length_pptflags,'L)'
    WRITE(fileu,myformat) 'flags: ',prop%flags


    end_subroutine()
END SUBROUTINE DTYPE(prop_print_info)

!!----------------------------------------------------------------
!! Procedures for Particle Sets DS
!!----------------------------------------------------------------
SUBROUTINE DTYPE(get_xp)(Pc,xp,with_ghosts)
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_particles))         :: Pc
    LOGICAL,OPTIONAL                      :: with_ghosts
    REAL(MK),DIMENSION(:,:),     POINTER  :: xp
    INTEGER                               :: info

    IF (PRESENT(with_ghosts)) THEN
        IF (with_ghosts) THEN
            IF (Pc%flags(ppm_part_ghosts)) THEN
                xp => Pc%xp(1:ppm_dim,1:Pc%Mpart)
            ELSE
                write(cbuf,*) 'WARNING: tried to get xp with ghosts ',&
                    'when ghosts are not up-to-date'
                CALL ppm_write(ppm_rank,'get_xp',cbuf,info)
                xp => NULL()
            ENDIF
            RETURN
        ENDIF
    ENDIF

    xp => Pc%xp(1:ppm_dim,1:Pc%Npart)

END SUBROUTINE DTYPE(get_xp)

SUBROUTINE DTYPE(set_xp)(Pc,xp,read_only,ghosts_ok)
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_particles))    :: Pc
    LOGICAL,OPTIONAL                 :: read_only
    LOGICAL,OPTIONAL                 :: ghosts_ok
    REAL(MK),DIMENSION(:,:),POINTER  :: xp
    INTEGER                          :: i

    CLASS(DTYPE(ppm_t_part_prop)_), POINTER :: prop => NULL()
    CLASS(DTYPE(ppm_t_neighlist)_), POINTER :: nl => NULL()
    CLASS(ppm_t_operator_discr_),   POINTER :: op => NULL()

    IF (PRESENT(ghosts_ok)) THEN
        IF (ghosts_ok) THEN
            xp => NULL()
            RETURN
        ENDIF
    ENDIF

    IF (PRESENT(read_only)) THEN
        IF (read_only) THEN
            xp => NULL()
            RETURN
        ENDIF
    ENDIF

    Pc%flags(ppm_part_areinside) = .FALSE.
    Pc%flags(ppm_part_partial) = .FALSE.
    Pc%flags(ppm_part_ghosts) = .FALSE.
    Pc%flags(ppm_part_cartesian) = .FALSE.

    prop => Pc%props%begin()
    DO WHILE (ASSOCIATED(prop))
        prop%flags(ppm_ppt_ghosts) = .FALSE.
        prop%flags(ppm_ppt_partial) = .FALSE.
        prop => Pc%props%next()
    ENDDO

    nl => Pc%neighs%begin()
    DO WHILE (ASSOCIATED(nl))
        nl%uptodate = .FALSE.
        nl => Pc%neighs%next()
    ENDDO

    IF (ASSOCIATED(Pc%ops)) THEN
        op => Pc%ops%begin()
        DO WHILE (ASSOCIATED(op))
            op%flags(ppm_ops_iscomputed) = .FALSE.
            op => Pc%ops%next()
        ENDDO
    ENDIF

    xp => NULL()

END SUBROUTINE DTYPE(set_xp)

SUBROUTINE DTYPE(part_prop_create)(Pc,id,datatype,info,&
        field,name,lda,zero,with_ghosts)
    !!! Adds a property to an existing particle set
    CLASS(DTYPE(ppm_t_particles))         :: Pc
    INTEGER,                INTENT(  OUT) :: id
    INTEGER,                INTENT(IN   ) :: datatype
    INTEGER, OPTIONAL,      INTENT(IN   ) :: lda
    CLASS(ppm_t_field_),OPTIONAL,INTENT(IN   ) :: field
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN  ) :: name
    !!! name to this property
    LOGICAL, OPTIONAL                     :: zero
    !!! if true, then initialise the data to zero
    LOGICAL, OPTIONAL                     :: with_ghosts
    !!! if true, then allocate with Mpart instead of the default size of Npart
    INTEGER,               INTENT(OUT)    :: info
    CHARACTER(LEN=ppm_char)               :: name2
    INTEGER                               :: lda2,vec_size,npart,i
    CLASS(DTYPE(ppm_t_part_prop)_),POINTER:: prop => NULL()
    LOGICAL, DIMENSION(ppm_param_length_pptflags):: flags

    start_subroutine("particle_prop_create")

    lda2 = 1
    IF (PRESENT(lda)) THEN
        IF (lda.GE.2) THEN
            lda2 = lda
        ENDIF
    ENDIF
    IF (PRESENT(name)) THEN
        name2 = name
    ELSE
        IF (PRESENT(field)) THEN
            name2 = field%name
        ELSE
            name2 = "default_ppt_name"
        ENDIF
    ENDIF

    npart = Pc%Npart
    flags = .FALSE.
    IF (PRESENT(with_ghosts)) THEN
        IF (with_ghosts) THEN
            IF (Pc%flags(ppm_part_ghosts)) THEN
                npart = Pc%Mpart
                flags(ppm_ppt_ghosts) = .TRUE.
            ELSE
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,caller,&
            'trying to init property for ghosts when ghosts are not computed',&
                    __LINE__,info)
                GOTO 9999
            ENDIF
        ENDIF
    ENDIF
    flags(ppm_ppt_partial) = .TRUE.
    flags(ppm_ppt_map_ghosts) = .TRUE.
    flags(ppm_ppt_map_parts) = .TRUE.


    ALLOCATE(DTYPE(ppm_t_part_prop)::prop,STAT=info)
        or_fail_alloc("prop")
    ! Create the property
    CALL prop%create(datatype,npart,lda2,name2,flags,info,field,zero)
        or_fail("creating property array failed")

    CALL Pc%props%push(prop,info,id)
        or_fail("pushing new property into collection failed")

    end_subroutine()
END SUBROUTINE DTYPE(part_prop_create)

SUBROUTINE DTYPE(part_prop_destroy)(Pc,id,info)
    !!! Destroy a property from an existing particle set
    CLASS(DTYPE(ppm_t_particles))         :: Pc
    INTEGER,                INTENT(INOUT) :: id
    INTEGER,               INTENT(OUT)    :: info

    start_subroutine("part_prop_destroy")

    ASSOCIATE (cont => Pc%props)
        IF (id .LE. 0 .OR. id .GT. cont%vec_size) THEN
            fail("property id larger than size of properties array")
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
                        fail("coding error in the data structure...")
                    ENDIF
                ENDDO
            ENDIF
        ENDIF
    END ASSOCIATE

    end_subroutine()
END SUBROUTINE DTYPE(part_prop_destroy)

SUBROUTINE DTYPE(part_prop_realloc)(Pc,id,info,with_ghosts,datatype,lda)
    !!! Reallocate the property array to the correct size
    !!! (e.g. if the number of particles has changed or if the type
    !!! of the data changes)
    CLASS(DTYPE(ppm_t_particles))         :: Pc
    INTEGER,                INTENT(IN   ) :: id
    INTEGER,               INTENT(OUT)    :: info
    LOGICAL, OPTIONAL                     :: with_ghosts
    !!! if true, then allocate with Mpart instead of the default size of Npart
    INTEGER, OPTIONAL                     :: datatype
    !!! deallocate the old data array and allocate a new one,
    !!! possibly of a different data type.
    INTEGER, OPTIONAL                     :: lda
    !!! deallocate the old data array and allocate a new one,
    !!! possibly of a different dimension

    CHARACTER(LEN=ppm_char)               :: name2
    INTEGER                               :: lda2,vec_size,npart,i,dtype
    CLASS(DTYPE(ppm_t_part_prop)_),POINTER:: prop => NULL()
    LOGICAL, DIMENSION(ppm_param_length_pptflags):: flags
    CLASS(ppm_t_field_),POINTER           :: field => NULL()

    start_subroutine("realloc_prop")

    IF (.NOT. ASSOCIATED(Pc%props%vec(id)%t)) THEN
        ALLOCATE(DTYPE(ppm_t_part_prop)::Pc%props%vec(id)%t,STAT=info)
            or_fail_alloc("allocating property pointer failed")
    ENDIF

    prop => Pc%props%vec(id)%t
    flags = prop%flags
    name2 = prop%name
    SELECT TYPE(f => prop%field_ptr)
    CLASS IS (ppm_t_field_)
        field => f
    END SELECT

    npart = Pc%Npart
    IF (PRESENT(with_ghosts)) THEN
        IF (with_ghosts) THEN
            IF (Pc%flags(ppm_part_ghosts)) THEN
                npart = Pc%Mpart
                flags(ppm_ppt_ghosts) = .TRUE.
            ELSE
    fail("trying to init property for ghosts when ghosts are not computed")
            ENDIF
        ENDIF
    ENDIF
    flags(ppm_ppt_partial) = .TRUE.

    IF (PRESENT(lda)) THEN
        lda2 = lda
    ELSE
        lda2 = prop%lda
    ENDIF
    IF (PRESENT(datatype)) THEN
        dtype = datatype
    ELSE
        dtype = prop%data_type
    ENDIF
    IF (lda2.NE.prop%lda .OR. dtype.NE.prop%data_type) THEN
        CALL prop%destroy(info)
    ENDIF

    ! Create the property
    CALL prop%create(dtype,npart,lda2,name2,flags,info,field)
        or_fail("reallocating property array failed")

    end_subroutine()
END SUBROUTINE DTYPE(part_prop_realloc)


SUBROUTINE DTYPE(part_neigh_create)(Pc,id,info,&
        P_id,name,skin,symmetry,cutoff)
    !!! Create a data structure to store a neighbour list
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_particles))         :: Pc
    INTEGER,               INTENT(  OUT)  :: id
    INTEGER,               INTENT(OUT)    :: info
    INTEGER, OPTIONAL                     :: P_id    
    !!! Id of the set of particles that this neighbor list refers to
    !!! The default, 0, stands for "self"
    CHARACTER(LEN=*) , OPTIONAL           :: name
    !!! name of this neighbour list
    REAL(MK), OPTIONAL                    :: skin
    REAL(MK), OPTIONAL                    :: cutoff
    LOGICAL, OPTIONAL                     :: symmetry    

    INTEGER                               :: vec_size,i
    CLASS(DTYPE(ppm_t_neighlist)_),POINTER :: Nlist=>NULL()
    REAL(MK), DIMENSION(:), POINTER       :: rcp => NULL()

    start_subroutine("particle_neigh_create")

    ! Create the neighbour list
    ALLOCATE(DTYPE(ppm_t_neighlist)::Nlist,STAT=info)
        or_fail_alloc("Nlist")

    IF (PRESENT(name)) THEN
        Nlist%name = name
    ELSE
        Nlist%name = particles_dflt_nlname(id)
    ENDIF

    IF (PRESENT(P_id)) THEN
        Nlist%P_id = P_id
    ELSE
        Nlist%P_id = 0
    ENDIF

    SELECT TYPE(Pc)
    TYPE IS (DTYPE(ppm_t_sop))
        ASSOCIATE (ghosts => Pc%flags(ppm_part_ghosts))
            IF (Pc%rcp_id.LE.0) THEN
                CALL Pc%create_prop(Pc%rcp_id,ppm_type_real,info,&
                    name='rcp',with_ghosts=ghosts) 
            ENDIF
            CALL Pc%get(rcp,Pc%rcp_id,with_ghosts=ghosts)
            IF (PRESENT(cutoff)) THEN
                rcp = cutoff
            ELSE
                rcp = Pc%ghostlayer
            ENDIF
            CALL Pc%set(rcp,Pc%rcp_id,ghosts_ok=ghosts)
        END ASSOCIATE
        Nlist%cutoff = -1._MK 
        !this field should not be used with adaptive particles
    CLASS DEFAULT
        IF (PRESENT(cutoff)) THEN
            Nlist%cutoff = cutoff
        ELSE
            Nlist%cutoff = Pc%ghostlayer
        ENDIF
    END SELECT

    IF (PRESENT(skin)) THEN
        Nlist%skin = skin
    ELSE
        Nlist%skin = 0._mk
    ENDIF

    IF (PRESENT(symmetry)) THEN
        IF (symmetry) THEN
            Nlist%isymm = 1
        ELSE
            Nlist%isymm = 0
        ENDIF
    ELSE
        Nlist%isymm = 0
    ENDIF

    Nlist%uptodate = .FALSE.
    Nlist%nneighmin = 0
    Nlist%nneighmax = 0

    CALL Pc%neighs%push(Nlist,info)
        or_fail("pushing new neighbour list into collection failed")

    end_subroutine()
END SUBROUTINE DTYPE(part_neigh_create)

SUBROUTINE DTYPE(part_neigh_destroy)(Pc,id,info)
    !!! Destroy a property from an existing particle set
    CLASS(DTYPE(ppm_t_particles))         :: Pc
    INTEGER,                INTENT(INOUT) :: id
    INTEGER,               INTENT(OUT)    :: info

    start_subroutine("part_neigh_destroy")

    ASSOCIATE (cont => Pc%neighs)
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

    end_subroutine()
END SUBROUTINE DTYPE(part_neigh_destroy)


SUBROUTINE DTYPE(part_create)(Pc,Npart,info,name)

    !!! create a ppm_t_particles data type

    !-------------------------------------------------------------------------
    !  Includes
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    !  Modules
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_particles))                          :: Pc
    !!! Data structure containing the particles
    INTEGER,                                INTENT(IN   )  :: Npart
    !!! Number of particles 
    INTEGER,                                INTENT(  OUT)  :: info
    !!! Returns status, 0 upon success.
    !-------------------------------------------------------------------------
    !  Optional arguments
    !-------------------------------------------------------------------------
    CHARACTER(LEN=*) , OPTIONAL                            :: name
    !!! give a name to this Particle set
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------

    LOGICAL                                         :: lalloc,ldealloc
    INTEGER                                         :: i

    start_subroutine("part_create")

    !-----------------------------------------------------------------
    !  Destroy the DS if it already exists
    !-----------------------------------------------------------------
    IF (ASSOCIATED(Pc%xp)) THEN
        CALL Pc%destroy(info)
            or_fail_dealloc("Pc%destroy")
    ENDIF

    !dumb way of creating a global ID for this mesh
    !TODO find something better? (needed if one creates and destroy
    ! many meshes)
    ppm_nb_part_sets = ppm_nb_part_sets + 1
    Pc%ID = ppm_nb_part_sets 

    !-----------------------------------------------------------------
    !  Allocate memory for the positions
    !-----------------------------------------------------------------
    ldc(1) = ppm_dim
    ldc(2) = Npart
    CALL ppm_alloc(Pc%xp,ldc(1:2),ppm_param_alloc_fit,info)
        or_fail_alloc("Pc%xp")
    Pc%Npart = Npart
    Pc%Mpart = Npart
    Pc%flags(ppm_part_ghosts) = .FALSE.
    Pc%flags(ppm_part_areinside) = .FALSE.
    Pc%flags(ppm_part_partial) = .FALSE.
    Pc%flags(ppm_part_reqput) = .FALSE.
    Pc%flags(ppm_part_cartesian) = .FALSE.
    Pc%flags(ppm_part_neighlists) = .FALSE.
    Pc%flags(ppm_part_global_index) = .FALSE.
    Pc%active_topoid = -1
    ! No active topology yet

    ! Give a default name to this Particle set
    IF (PRESENT(name)) THEN
        Pc%name = ADJUSTL(TRIM(name))
    ELSE
        Pc%name = particles_dflt_partname()
    ENDIF

    ! Particles have not been initialised yet
    Pc%h_avg = -1._MK
    Pc%h_min = -1._MK

    Pc%time = 0._MK
    Pc%itime = 0

    Pc%gi_id = 0

    SELECT TYPE(Pc)
    CLASS IS (DTYPE(ppm_t_sop))
        !-----------------------------------------------------------------
        !  Initialize fields of the extended SOP type
        !-----------------------------------------------------------------
        ! Particles are by default not adaptive
        Pc%adaptive = .FALSE.
        Pc%adapt_wpid = 0
        Pc%rcp_id = 0
        Pc%D_id = 0
        Pc%Dtilde_id = 0
        Pc%nn_sq_id = 0
        ! Particles do not represent a level-set function
        Pc%level_set = .FALSE.
        Pc%level_id = 0
        !        Pc%level_old_id = 0
        Pc%level_grad_id = 0
        !        Pc%level_grad_old_id = 0
        ! Particles are by default isotropic
        Pc%anisotropic = .FALSE.
    CLASS DEFAULT
    END SELECT

    IF (.NOT.ASSOCIATED(Pc%neighs)) THEN
        ALLOCATE(DTYPE(ppm_c_neighlist)::Pc%neighs,STAT=info)
        or_fail_alloc("could not allocate Pc%neighs")
    ELSE
        fail("neighbour list collection is already allocated. Call destroy() first?")
    ENDIF

    IF (ASSOCIATED(Pc%ops)) THEN
        fail("collection of operator pointers is already associated. Call destroy() first?")
    ENDIF

    IF (.NOT.ASSOCIATED(Pc%props)) THEN
        ALLOCATE(DTYPE(ppm_c_part_prop)::Pc%props,STAT=info)
        or_fail_alloc("could not allocate Pc%props")
    ELSE
        fail("property collection is already allocated. Call destroy() first?")
    ENDIF

    IF (.NOT.ALLOCATED(Pc%stats)) THEN
        ALLOCATE(DTYPE(particles_stats)::Pc%stats,STAT=info)
        or_fail_alloc("could not allocate Pc%stats")
    ELSE
        fail("stats structure is already allocated. Call destroy() first?")
    ENDIF



    end_subroutine()
END SUBROUTINE DTYPE(part_create)


SUBROUTINE DTYPE(part_destroy)(Pc,info)
    !!! Deallocate a ppm_t_particles data type

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_particles))                          :: Pc
    !!! Data structure containing the Pc
    INTEGER,                                INTENT(  OUT)  :: info
    !!! Returns status, 0 upon success.

    LOGICAL                                         :: lalloc,ldealloc
    INTEGER                                         :: i

    start_subroutine("part_destroy")

    Pc%ID = 0
    ! first deallocate all content of Pc
    dealloc_pointer(Pc%xp)

    dealloc_pointer(Pc%pcost)

    dealloc_allocatable(Pc%stats)

    !Deallocate neighbour lists
    destroy_collection_ptr(Pc%neighs)

    !Deallocate properties
    destroy_collection_ptr(Pc%props)

    !Deallocate operators
    destroy_collection_ptr(Pc%ops)

    !-------------------------------------------------------------------------
    !  Finalize
    !-------------------------------------------------------------------------
    end_subroutine()
END SUBROUTINE DTYPE(part_destroy)

#define __DIM 2
#include "part/ppm_particles_initialize.f"
#define __DIM 3
#include "part/ppm_particles_initialize.f"

!!temporary hack to deal with both 2d and 3d
SUBROUTINE DTYPE(part_initialize)(Pc,Npart_global,info,&
        distrib,topoid,minphys,maxphys,cutoff,name)
    !-----------------------------------------------------------------------
    ! Set initial particle positions
    !-----------------------------------------------------------------------

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_particles))                          :: Pc
    !!! Data structure containing the particles
    INTEGER,                            INTENT(INOUT)      :: Npart_global
    !!! total number of particles that will be initialized
    INTEGER,                            INTENT(  OUT)      :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    !  Optional arguments
    !-------------------------------------------------------------------------
    INTEGER,OPTIONAL,                   INTENT(IN   )      :: distrib
    !!! type of initial distribution. One of
    !!! ppm_param_part_init_cartesian (default)
    !!! ppm_param_part_init_random
    INTEGER,OPTIONAL,                   INTENT(IN   )      :: topoid
    !!! topology id (used only to get the extent of the physical domain)
    REAL(MK),DIMENSION(ppm_dim),OPTIONAL,INTENT(IN   )     :: minphys
    !!! extent of the physical domain. Only if topoid is not present.
    REAL(MK),DIMENSION(ppm_dim),OPTIONAL,INTENT(IN   )     :: maxphys
    !!! extent of the physical domain. Only if topoid is not present.
    REAL(MK),                   OPTIONAL,INTENT(IN   )     :: cutoff
    !!! cutoff of the particles
    CHARACTER(LEN=*),           OPTIONAL,INTENT(IN   )     :: name
    !!! name for this set of particles

    start_subroutine("part_initialize")

    IF (ppm_dim .eq. 2) THEN
        CALL  DTYPE(particles_initialize2d)(Pc,Npart_global,info,&
            distrib,topoid,minphys,maxphys,cutoff,name=name)
    ELSE
        CALL  DTYPE(particles_initialize3d)(Pc,Npart_global,info,&
            distrib,topoid,minphys,maxphys,cutoff,name=name)
    ENDIF

    end_subroutine()
END SUBROUTINE DTYPE(part_initialize)

SUBROUTINE DTYPE(part_print_info)(Pc,info,level,fileunit)
    !-----------------------------------------------------------------------
    ! Print out summary information about this Particle set
    ! (list of properties, operators, etc...)
    !-----------------------------------------------------------------------

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_particles))                          :: Pc
    !!! Data structure containing the particles
    INTEGER,                            INTENT(  OUT)      :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    !  Optional arguments
    !-------------------------------------------------------------------------
    INTEGER, OPTIONAL,                  INTENT(IN   )      :: level
    !!! indentation level at which to printout the info. Default = 0
    INTEGER, OPTIONAL,                  INTENT(IN   )      :: fileunit
    !!! Already open file unit for printout. Default = stdout
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    INTEGER                              :: lev,fileu
    CHARACTER(LEN = ppm_char)            :: myformat
    CLASS(DTYPE(ppm_t_part_prop)_),POINTER :: prop => NULL()

    start_subroutine("part_print_info")


    IF (PRESENT(fileunit)) THEN
        fileu = fileunit
    ELSE
        fileu = 6
    ENDIF
    IF (PRESENT(level)) THEN
        lev = MAX(level,1)
    ELSE
        lev = 1
    ENDIF

    WRITE(myformat,'(A,I0,A)') '(',4*lev,'X,A,A,2X,2(A,I0),A)'

    WRITE(fileu,myformat) 'Particle set: ',TRIM(color_print(Pc%name,31)),&
       '(N = ',Pc%Npart,' M = ',Pc%Mpart,')'

    lev = lev + 1

    WRITE(myformat,'(A,I0,A,I0,A)') '(',4*lev,'X,A,',&
        ppm_param_length_partflags,'L)'
    WRITE(fileu,myformat) 'flags: ',Pc%flags

    prop => Pc%props%begin()
    DO WHILE (ASSOCIATED(prop))
        CALL prop%print_info(info,lev,fileunit,Pc%props%iter_id)
        prop => Pc%props%next()
    ENDDO

    end_subroutine()
END SUBROUTINE DTYPE(part_print_info)

SUBROUTINE DTYPE(part_del_parts)(Pc,list_del_parts,nb_del,info)
    !!! remove some particles from a Particle set
    !!! WARNING: this implementation is NOT efficient
    !!! if the number of particles to delete is large.

    !-------------------------------------------------------------------------
    !  Includes
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    !  Modules
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_particles))                          :: Pc
    !!! Data structure containing the particles
    INTEGER,DIMENSION(:),POINTER,           INTENT(IN   )  :: list_del_parts
    !!! list of particles to be deleted
    INTEGER,                                INTENT(IN   )  :: nb_del
    !!! number of particles to be deleted
    INTEGER,                                INTENT(  OUT)  :: info
    !!! Returns status, 0 upon success.
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------

    INTEGER                              :: i,ip,Npart,del_part,lda
    CLASS(DTYPE(ppm_t_part_prop)_),POINTER :: prop => NULL()

    start_subroutine("part_del_parts")

    !-----------------------------------------------------------------
    !  check arguments
    !-----------------------------------------------------------------
    IF (.NOT.ASSOCIATED(Pc%xp)) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'Pc structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    prop => Pc%props%begin()
    DO WHILE (ASSOCIATED(prop))
        IF (.NOT.prop%flags(ppm_ppt_partial)) THEN
            IF (prop%flags(ppm_ppt_map_parts)) THEN
                info = ppm_error_error
                CALL ppm_error(999,caller,   &
                    & 'property not mapped, data will be lost',&
                    &  __LINE__,info)
                GOTO 9999
            ENDIF
        ENDIF
        prop => Pc%props%next()
    ENDDO

    !-----------------------------------------------------------------
    !  Delete particles
    !-----------------------------------------------------------------
    Npart = Pc%Npart

    del_part = 1

    DO i=1,nb_del
        ip = list_del_parts(i)

        ! copying particles from the end of xp to the index that has
        ! to be removed
        Pc%xp(1:ppm_dim,ip) = Pc%xp(1:ppm_dim,Npart-i+1)

        prop => Pc%props%begin()
        DO WHILE (ASSOCIATED(prop))
            IF (prop%flags(ppm_ppt_partial)) THEN

                    lda = prop%lda

                    IF (lda.GE.2) THEN
                        SELECT CASE (prop%data_type)
                        CASE (ppm_type_int)
                            prop%data_2d_i(1:lda,ip) = &
                                prop%data_2d_i(1:lda,Npart-i+1)
                        CASE (ppm_type_longint)
                            prop%data_2d_li(1:lda,ip) = &
                                prop%data_2d_li(1:lda,Npart-i+1)
                        CASE (ppm_type_real)
                            prop%data_2d_r(1:lda,ip) = &
                                prop%data_2d_r(1:lda,Npart-i+1)
                        CASE (ppm_type_comp)
                            prop%data_2d_c(1:lda,ip) = &
                                prop%data_2d_c(1:lda,Npart-i+1)
                        CASE (ppm_type_logical )
                            prop%data_2d_l(1:lda,ip) = &
                                prop%data_2d_l(1:lda,Npart-i+1)
                        END SELECT

                    ELSE
                        SELECT CASE (prop%data_type)
                        CASE (ppm_type_int)
                            prop%data_1d_i(ip) = &
                                prop%data_1d_i(Npart-i+1)
                        CASE (ppm_type_longint)
                            prop%data_1d_li(ip) = &
                                prop%data_1d_li(Npart-i+1)
                        CASE (ppm_type_real)
                            prop%data_1d_r(ip) = &
                                prop%data_1d_r(Npart-i+1)
                        CASE (ppm_type_comp)
                            prop%data_1d_c(ip) = &
                                prop%data_1d_c(Npart-i+1)
                        CASE (ppm_type_logical )
                            prop%data_1d_l(ip) = &
                                prop%data_1d_l(Npart-i+1)
                        END SELECT
                    ENDIF
                prop => Pc%props%next()
            ENDIF
        ENDDO
    ENDDO
    !New number of particles, after deleting some
    Pc%Npart = Npart - del_part

    end_subroutine()
END SUBROUTINE DTYPE(part_del_parts)

SUBROUTINE DTYPE(part_set_rel)(this,field,info)
    !!! Create bookkeeping data structure to log the relationship between
    !!! the particle set and a field
    CLASS(DTYPE(ppm_t_particles))      :: this
    CLASS(ppm_t_field_)                :: field
    INTEGER,               INTENT(OUT) :: info
    TYPE(ppm_t_field_info),POINTER     :: fdinfo => NULL()

    start_subroutine("part_set_rel")


    ALLOCATE(fdinfo,STAT=info)
        or_fail_alloc("could not allocate new field_info pointer")


    CALL fdinfo%create(field,info)
        or_fail("could not create new field_info object")

    IF (.NOT.ASSOCIATED(this%field_ptr)) THEN
        ALLOCATE(ppm_c_field_info::this%field_ptr,STAT=info)
            or_fail_alloc("could not allocate field_info collection")
    ENDIF

    IF (this%field_ptr%vec_size.LT.field%ID) THEN
        CALL this%field_ptr%grow_size(info)
            or_fail("could not grow_size this%field_ptr")
    ENDIF
    !if the collection is still too small, it means we should consider
    ! using a hash table for particle set IDs...
    IF (this%field_ptr%vec_size.LT.field%ID) THEN
        fail("The id of the partile set that we are trying to store in a Field collection seems to large. Why not implement a hash function?")
    ENDIF

    IF (this%field_ptr%exists(field%ID)) THEN
        fail("It seems like this particle set already has a pointer to that field. Are you sure you want to do that a second time?")
    ELSE
        !add field_info to the collection, at index fieldID
        !Update the counters nb,max_id and min_id accordingly
        !(this is not very clean without hash table)
        this%field_ptr%vec(field%ID)%t => fdinfo
        this%field_ptr%nb = this%field_ptr%nb + 1
        IF (this%field_ptr%nb.EQ.1) THEN
            this%field_ptr%max_id = field%ID
            this%field_ptr%min_id = field%ID
        ELSE
            this%field_ptr%max_id = MAX(this%field_ptr%max_id,field%ID)
            this%field_ptr%min_id = MIN(this%field_ptr%min_id,field%ID)
        ENDIF
    ENDIF

    end_subroutine()
END SUBROUTINE DTYPE(part_set_rel)

SUBROUTINE DTYPE(part_prop_push)(Pc,prop_id,info)

    !!! wrapper for ppm_map_part_push

    !-------------------------------------------------------------------------
    !  Includes
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    !  Modules
    !-------------------------------------------------------------------------
    USE ppm_module_map

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_particles))                          :: Pc
    !!! Data structure containing the particles
    INTEGER,                                INTENT(IN   )  :: prop_id
    !!! id of the property to be pushed
    INTEGER,                                INTENT(  OUT)  :: info
    !!! Returns status, 0 upon success.

    INTEGER                              :: lda
    CLASS(DTYPE(ppm_t_part_prop)_),POINTER :: prop => NULL()
    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    start_subroutine("map_part_push")

    !-----------------------------------------------------------------
    !  Call ppm_map_part_push
    !-----------------------------------------------------------------
    prop => Pc%props%vec(prop_id)%t
    lda = prop%lda

    IF (lda.GE.2) THEN
        SELECT CASE (prop%data_type)
        CASE (ppm_type_int)
            CALL ppm_map_part_push(&
            prop%data_2d_i,lda,Pc%Npart,info)
        CASE (ppm_type_longint)
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,   &
                &  'Type not supported for mappings.',&
                &  __LINE__,info)
            GOTO 9999
        CASE (ppm_type_real)
            CALL ppm_map_part_push(&
            prop%data_2d_r,lda,Pc%Npart,info)
        CASE (ppm_type_comp)
            CALL ppm_map_part_push(&
            prop%data_2d_c,lda,Pc%Npart,info)
        CASE (ppm_type_logical )
            CALL ppm_map_part_push(&
            prop%data_2d_l,lda,Pc%Npart,info)
        END SELECT

    ELSE

        SELECT CASE (prop%data_type)
        CASE (ppm_type_int)
            CALL ppm_map_part_push(&
            prop%data_1d_i,Pc%Npart,info)
        CASE (ppm_type_longint)
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,   &
                &  'Type not supported for mappings.',&
                &  __LINE__,info)
            GOTO 9999
        CASE (ppm_type_real)
            CALL ppm_map_part_push(&
            prop%data_1d_r,Pc%Npart,info)
        CASE (ppm_type_comp)
            CALL ppm_map_part_push(&
            prop%data_1d_c,Pc%Npart,info)
        CASE (ppm_type_logical )
            CALL ppm_map_part_push(&
            prop%data_1d_l,Pc%Npart,info)
        END SELECT
    ENDIF

    end_subroutine()
END SUBROUTINE DTYPE(part_prop_push)

SUBROUTINE DTYPE(part_prop_pop)(Pc,prop_id,Npart_new,info)

    !!! wrapper for ppm_map_part_pop

    !-------------------------------------------------------------------------
    !  Includes
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    !  Modules
    !-------------------------------------------------------------------------
    USE ppm_module_map

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_particles))                          :: Pc
    !!! Data structure containing the particles
    INTEGER,                                INTENT(IN   )  :: prop_id
    !!! id of the property to be pushed
    INTEGER,                                INTENT(IN   )  :: Npart_new
    !!! number of particles to pop from the buffer
    INTEGER,                                INTENT(  OUT)  :: info
    !!! Returns status, 0 upon success.

    INTEGER                                :: lda
    CLASS(DTYPE(ppm_t_part_prop)_),POINTER :: prop => NULL()

    start_subroutine("map_part_pop")

    !-----------------------------------------------------------------
    !  Call ppm_map_part_pop
    !-----------------------------------------------------------------
    prop => Pc%props%vec(prop_id)%t
    lda = prop%lda

    IF (lda.GE.2) THEN
        SELECT CASE (prop%data_type)
        CASE (ppm_type_int)
            CALL ppm_map_part_pop(&
            prop%data_2d_i,lda,Pc%Npart,Npart_new,info)
        CASE (ppm_type_longint)
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,   &
                &  'Type not supported for mappings.',&
                &  __LINE__,info)
            GOTO 9999
        CASE (ppm_type_real)
            CALL ppm_map_part_pop(&
            prop%data_2d_r,lda,Pc%Npart,Npart_new,info)
        CASE (ppm_type_comp)
            CALL ppm_map_part_pop(&
            prop%data_2d_c,lda,Pc%Npart,Npart_new,info)
        CASE (ppm_type_logical )
            CALL ppm_map_part_pop(&
            prop%data_2d_l,lda,Pc%Npart,Npart_new,info)
        END SELECT

    ELSE

        SELECT CASE (prop%data_type)
        CASE (ppm_type_int)
            CALL ppm_map_part_pop(&
            prop%data_1d_i,Pc%Npart,Npart_new,info)
        CASE (ppm_type_longint)
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,   &
                &  'Type not supported for mappings.',&
                &  __LINE__,info)
            GOTO 9999
        CASE (ppm_type_real)
            CALL ppm_map_part_pop(&
            prop%data_1d_r,Pc%Npart,Npart_new,info)
        CASE (ppm_type_comp)
            CALL ppm_map_part_pop(&
            prop%data_1d_c,Pc%Npart,Npart_new,info)
        CASE (ppm_type_logical )
            CALL ppm_map_part_pop(&
            prop%data_1d_l,Pc%Npart,Npart_new,info)
        END SELECT
    ENDIF

    end_subroutine()
END SUBROUTINE DTYPE(part_prop_pop)


SUBROUTINE DTYPE(part_mapping)(Pc,info,debug,global,topoid)

    !!!  Partial/Global mapping for particles
    !!!  Assumptions:
    !!! * All the particles have to be inside the domain
    !!!   (otherwise -> "unassigned particle error")

    USE ppm_module_map
#ifdef __MPI
    INCLUDE "mpif.h"
#endif

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    CLASS(DTYPE(ppm_t_particles))                          :: Pc
    DEFINE_MK()
    !!! Data structure containing the particles
    INTEGER,                            INTENT(  OUT)      :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    !  Optional Arguments
    !-------------------------------------------------------------------------
    LOGICAL, OPTIONAL                                   :: debug
    !!! IF true, printout more
    LOGICAL, OPTIONAL                                   :: global
    !!! does a global mapping. Default is false (i.e. partial mapping)
    INTEGER, OPTIONAL                                   :: topoid
    !!! topology id
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    INTEGER                   :: Npart_new
    !!! new number of particles on this processor
    INTEGER                   :: ltopoid
    !!! index variable
    REAL(KIND(1.D0))          :: t1,t2
    LOGICAL                   :: dbg,partial
    CLASS(DTYPE(ppm_t_part_prop)_), POINTER :: prop => NULL()
    CLASS(DTYPE(ppm_t_neighlist)_), POINTER :: nl => NULL()
    CLASS(ppm_t_operator_discr_),   POINTER :: op => NULL()

    start_subroutine("particles_mapping")

#ifdef __MPI
    t1 = MPI_WTIME(info)
#endif
    dbg = .FALSE.
    IF (PRESENT(debug)) dbg = debug
    !-----------------------------------------------------------------
    !  Checks
    !-----------------------------------------------------------------
    IF (.NOT.ASSOCIATED(Pc%xp)) THEN
        fail("Pc structure had not been defined. Call allocate first")
    ENDIF
    IF (.NOT.Pc%flags(ppm_part_areinside)) THEN
        fail("some Pc may be outside the domain. Apply BC first")
    ENDIF

    IF (PRESENT(global)) THEN
        IF(.NOT.PRESENT(topoid)) THEN
            fail("need the topoid parameter for global mapping")
        ENDIF
        IF (global) partial = .FALSE.
    ELSE
        partial = .TRUE.
    ENDIF

    IF (partial) THEN
        ltopoid = Pc%active_topoid
    ELSE
        ltopoid = topoid
    ENDIF

    !-----------------------------------------------------------------------
    !  Map the particles onto the topology
    !-----------------------------------------------------------------------
    IF (partial .AND. Pc%flags(ppm_part_partial)) THEN
        !Particles have already been mapped onto this topology
        !nothing to do
    ELSE
#ifdef __MPI
        t1 = MPI_WTIME(info)
#endif
        IF (partial) THEN
            CALL ppm_map_part_partial(ltopoid,Pc%xp,Pc%Npart,info) 
            Pc%stats%nb_part_map = Pc%stats%nb_part_map + 1
            IF (info .NE. 0) THEN
                info = ppm_error_error
                CALL ppm_error(999,caller,&
                    'ppm_map_part_partial failed',__LINE__,info)
                GOTO 9999
            ENDIF
        ELSE
            CALL ppm_map_part_global(ltopoid,Pc%xp,Pc%Npart,info) 
            Pc%stats%nb_global_map = Pc%stats%nb_global_map + 1
            IF (info .NE. 0) THEN
                info = ppm_error_error
                CALL ppm_error(999,caller,&
                    'ppm_map_part_global failed',__LINE__,info)
                GOTO 9999
            ENDIF
        ENDIF

        prop => Pc%props%begin()
        DO WHILE (ASSOCIATED(prop))
            IF (prop%flags(ppm_ppt_map_parts)) THEN
                IF(dbg) &
                    write(*,*) 'pushing property ',Pc%props%iter_id
                CALL Pc%map_part_push_legacy(Pc%props%iter_id,info)
                IF (info .NE. 0) THEN
                    info = ppm_error_error
                    CALL ppm_error(999,caller,&
                        'ppm_map_part_prop_push failed',__LINE__,info)
                    GOTO 9999
                ENDIF
            ENDIF
            prop => Pc%props%next()
        ENDDO

        CALL ppm_map_part_send(Pc%Npart,Npart_new,info)
        IF (info .NE. 0) THEN
            CALL ppm_error(0,caller,&
                'ppm_map_part_send failed',__LINE__,info)
            GOTO 9999
        ENDIF

        prop => Pc%props%last()
        DO WHILE (ASSOCIATED(prop))
            IF (prop%flags(ppm_ppt_map_parts)) THEN
                IF(dbg) &
                    write(*,*) 'poping property ',Pc%props%iter_id
                CALL Pc%map_part_pop_legacy(Pc%props%iter_id,Npart_new,info)
                IF (info .NE. 0) THEN
                    info = ppm_error_error
                    CALL ppm_error(999,caller,&
                        'ppm_map_part_prop_pop failed',__LINE__,info)
                    GOTO 9999
                ENDIF
                prop%flags(ppm_ppt_partial) = .TRUE.
            ENDIF
            prop => Pc%props%prev()
        ENDDO

        IF(dbg) &
            write(*,*) 'popping xp'
        CALL ppm_map_part_pop(Pc%xp,ppm_dim,Pc%Npart,Npart_new,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(999,caller,&
                'ppm_map_part_pop failed',__LINE__,info)
            GOTO 9999
        ENDIF

        ! Update states
        ! Number of particles on this processor
        Pc%Npart = Npart_new
        Pc%Mpart = Pc%Npart

        ! This is the active topology for these particles
        IF (.NOT.partial) Pc%active_topoid = topoid

        ! Particles are now mapped on the active topology
        Pc%flags(ppm_part_partial) = .TRUE.
        ! Pc have been re-indexed and ghosts have not been computed
        Pc%flags(ppm_part_ghosts) = .FALSE.

        !   values for poperty arrays have been mapped and ghosts
        !   are no longer up-to-date
        prop => Pc%props%begin()
        DO WHILE (ASSOCIATED(prop))
            prop%flags(ppm_ppt_ghosts) = .FALSE.
            prop => Pc%props%next()
        ENDDO

        ! particles have been re-indexed and neighbour lists not updated
        nl => Pc%neighs%begin()
        DO WHILE (ASSOCIATED(nl))
            nl%uptodate = .FALSE.
            nl => Pc%neighs%next()
        ENDDO
        Pc%flags(ppm_part_neighlists) = .FALSE.

        ! particles have been re-indexed and operators need be recomputed
        IF (ASSOCIATED(Pc%ops)) THEN
            op => Pc%ops%begin()
            DO WHILE (ASSOCIATED(op))
                op%flags(ppm_ops_iscomputed) = .FALSE.
                op => Pc%ops%next()
            ENDDO
        ENDIF

    ENDIF

#ifdef __MPI
    t2 = MPI_WTIME(info)
    IF (partial) THEN
        Pc%stats%t_part_map = Pc%stats%t_part_map + (t2-t1)
    ELSE
        Pc%stats%t_global_map = Pc%stats%t_global_map + (t2-t1)
    ENDIF
#endif

    end_subroutine()
END SUBROUTINE DTYPE(part_mapping)

SUBROUTINE DTYPE(part_mapping_ghosts)(Pc,info,ghostsize,debug)

    !!!  Ghost mapping for particles
    !!!  Assumptions:
    !!! * Particles positions need to have been mapped onto the topology
    !!!
    USE ppm_module_map
#ifdef __MPI
    INCLUDE "mpif.h"
#endif

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    CLASS(DTYPE(ppm_t_particles))                       :: Pc
    DEFINE_MK()
    !!! Data structure containing the particles
    INTEGER,                            INTENT(  OUT)   :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    !  Optional Arguments
    !-------------------------------------------------------------------------
    REAL(MK), OPTIONAL                                  :: ghostsize
    !!! size of the ghost layers. Default is to use the particles cutoff
    LOGICAL, OPTIONAL                                   :: debug
    !!! IF true, printout more
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    INTEGER                                   :: topoid
    !!! index variable
    REAL(MK)                                  :: cutoff
    !!! cutoff radius
    TYPE(ppm_t_topo),POINTER                  :: topo => NULL()
    REAL(KIND(1.D0))                          :: t1,t2
    LOGICAL                                   :: dbg
    LOGICAL                                   :: skip_ghost_get
    LOGICAL                                   :: skip_send
    CLASS(DTYPE(ppm_t_part_prop)_), POINTER :: prop => NULL()
    CLASS(DTYPE(ppm_t_neighlist)_), POINTER :: nl => NULL()
    CLASS(ppm_t_operator_discr_),   POINTER :: op => NULL()

    start_subroutine("part_mapping_ghosts")

    dbg = .FALSE.
    IF (PRESENT(debug)) dbg=debug
    skip_ghost_get = .FALSE.
    skip_send = .TRUE.
    !we must not call ppm_map_part_send unless ppm_map_part_push (or ghost_get)
    ! has been called (in which case, skip_send is set to FALSE)

    !-----------------------------------------------------------------
    !  Checks
    !-----------------------------------------------------------------
    IF (.NOT.ASSOCIATED(Pc%xp)) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (.NOT.Pc%flags(ppm_part_partial)) THEN
        !Particles have not been mapped onto this topology
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            &  'Do a partial/global mapping before doing a ghost mapping',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
    IF (.NOT.Pc%flags(ppm_part_areinside)) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            &  'some particles may be outside the domain. Apply BC first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF


    topoid = Pc%active_topoid
    topo=>ppm_topo(topoid)%t

    cutoff = Pc%ghostlayer
    IF (PRESENT(ghostsize)) THEN
        IF (ghostsize .LT. cutoff) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,   &
                &  'using ghostsize < cutoff+skin. Increase ghostsize.',&
                &  __LINE__,info)
            GOTO 9999
        ELSE
            cutoff = ghostsize
        ENDIF
    ENDIF

#if   __KIND == __SINGLE_PRECISION
    IF (cutoff .GT. topo%ghostsizes) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            &  'ghostsize of topology may be smaller than that of particles',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
#elif   __KIND == __DOUBLE_PRECISION
    IF (cutoff .GT. topo%ghostsized) THEN

        write(*,*) 'cutoff = ',cutoff
        write(*,*) 'cutoff used to create topology = ',topo%ghostsized

        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            &  'ghostsize of topology may be smaller than that of particles',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
#endif
    IF (cutoff .GT. 0._MK) THEN
        IF (Pc%flags(ppm_part_ghosts)) THEN
            IF (dbg) THEN
                write(*,*) 'ghosts have already been updated'
            ENDIF

            IF (ppm_map_type_isactive(ppm_param_map_ghost_get)) THEN
                IF (dbg) THEN
                    write(*,*) 'we skip the ghost_get and go straight to'
                    write(*,*) 'push/send/pop'
                ENDIF
                skip_ghost_get = .TRUE.
            ENDIF
        ENDIF

        IF (.NOT.skip_ghost_get) THEN
            Pc%stats%nb_ghost_get = Pc%stats%nb_ghost_get + 1
#ifdef __MPI
            t1 = MPI_WTIME(info)
#endif
            IF(dbg) &
                write(*,*) 'ghost-get '
            CALL ppm_map_part_ghost_get(topoid,Pc%xp,ppm_dim,&
                Pc%Npart,Pc%isymm,cutoff,info)
            IF (info .NE. 0) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_sub_failed,caller,&
                    'ppm_map_part_ghost_get failed',__LINE__,info)
                GOTO 9999
            ENDIF
#ifdef __MPI
            t2 = MPI_WTIME(info)
            Pc%stats%t_ghost_get = Pc%stats%t_ghost_get + (t2-t1)
#endif
            skip_send = .FALSE.
        ELSE
            IF(dbg) &
                write(*,*) 'skipping ghost-get '
        ENDIF

        !Update the ghost for the properties if
        ! 1) they have been mapped to this topology,
        ! 2) the ghosts have not yet been updated, and
        ! 3) the user wants them to be updated
        prop => Pc%props%begin()
        DO WHILE (ASSOCIATED(prop))

            IF (prop%flags(ppm_ppt_map_ghosts)) THEN
                IF (.NOT.prop%flags(ppm_ppt_ghosts)) THEN
                    IF (prop%flags(ppm_ppt_partial)) THEN

                        IF(dbg) &
                            write(*,*) 'pushing property ',Pc%props%iter_id,&
                            TRIM(prop%name)
                            Pc%stats%nb_ghost_push = &
                                Pc%stats%nb_ghost_push + 1
#ifdef __MPI
                        t1 = MPI_WTIME(info)
#endif
                        CALL Pc%map_part_push_legacy(Pc%props%iter_id,info)
                        IF (info .NE. 0) THEN
                            info = ppm_error_error
                            CALL ppm_error(ppm_err_sub_failed,caller,&
                                'ppm_map_part_push failed',__LINE__,info)
                            GOTO 9999
                        ENDIF
#ifdef __MPI
                        t2 = MPI_WTIME(info)
                        Pc%stats%t_ghost_push = &
                            Pc%stats%t_ghost_push + (t2-t1)
#endif
                        skip_send = .FALSE.
                    ELSE
                        write(*,*) 'pushing property ',&
                            Pc%props%iter_id,TRIM(prop%name)
                        info = ppm_error_error
                        CALL ppm_error(ppm_err_argument,caller,&
                            'getting ghosts for a property thats not mapped',&
                            __LINE__,info)
                        GOTO 9999
                    ENDIF
                ENDIF
            ENDIF
            prop => Pc%props%next()
        ENDDO

        IF (.NOT. skip_send) THEN
            CALL ppm_map_part_send(Pc%Npart,Pc%Mpart,info)
            IF (info .NE. 0) THEN
                write(*,*) 'ppm_map_part_send failed with info = ',info
                info = ppm_error_error
                CALL ppm_error(ppm_err_sub_failed,caller,&
                    'ppm_map_part_send failed',__LINE__,info)
                GOTO 9999
            ENDIF

            prop => Pc%props%last()
            DO WHILE (ASSOCIATED(prop))

                IF (prop%flags(ppm_ppt_map_ghosts)) THEN
                    IF (.NOT.prop%flags(ppm_ppt_ghosts)) THEN
                        IF (prop%flags(ppm_ppt_partial)) THEN

                            IF(dbg) &
                                write(*,*) 'popping property ',Pc%props%iter_id,&
                                TRIM(prop%name)
                            CALL Pc%map_part_pop_legacy(Pc%props%iter_id,Pc%Mpart,info)
                            IF (info .NE. 0) THEN
                                write(*,*) 'popping property ',Pc%props%iter_id,&
                                    TRIM(prop%name)
                                info = ppm_error_error
                                CALL ppm_error(ppm_err_sub_failed,caller,&
                                    'ppm_map_part_pop failed',__LINE__,info)
                                GOTO 9999
                            ENDIF
                            prop%flags(ppm_ppt_ghosts) = .TRUE.
                        ENDIF
                    ENDIF
                ENDIF
                prop => Pc%props%prev()
            ENDDO

            IF (.NOT.skip_ghost_get) THEN
                IF(dbg) &
                    write(*,*) 'popping-xp '
                CALL ppm_map_part_pop(Pc%xp,ppm_dim,Pc%Npart,&
                    Pc%Mpart,info)
                IF (info .NE. 0) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_sub_failed,caller,&
                        'ppm_map_part_pop failed',__LINE__,info)
                    GOTO 9999
                ENDIF
            ENDIF
        ENDIF !.NOT.skip_send

    ELSE ! if cutoff .le. 0

        IF(dbg) THEN
            write(*,*) 'cutoff = 0, nothing to do'
            write(*,*) 'setting all %has_ghost properties to true'
        ENDIF

        prop => Pc%props%begin()
        DO WHILE (ASSOCIATED(prop))
            IF (prop%flags(ppm_ppt_map_ghosts)) THEN
                prop%flags(ppm_ppt_ghosts) = .TRUE.
            ENDIF
            prop => Pc%props%next()
        ENDDO
    ENDIF


    ! Update states
    !   ghosts have been computed
    Pc%flags(ppm_part_ghosts) = .TRUE.
    ! the states for the properties have already been updated above

    end_subroutine()
END SUBROUTINE DTYPE(part_mapping_ghosts)

SUBROUTINE DTYPE(part_apply_bc)(Pc,info)

    !!!  Apply boundary conditions for particles positions
    !!!  Assumptions:
    !!! * Particles positions need to have been mapped onto the topology

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    CLASS(DTYPE(ppm_t_particles))                         :: Pc
    DEFINE_MK()
    !!! Data structure containing the particles
    INTEGER,                            INTENT(  OUT)     :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    REAL(MK), DIMENSION(:,:),POINTER                      :: xp => NULL()
    !!! pointer to positions
    TYPE(ppm_t_topo),POINTER                              :: topo => NULL()
    !!! pointer to topology
    REAL(MK), DIMENSION(ppm_dim)                          :: min_phys,max_phys
    !!! computational domain corners
    REAL(MK), DIMENSION(ppm_dim)                          :: len_phys
    !!! length of the computational domain
    INTEGER                                               :: di,ip
    INTEGER                                               :: topoid
    INTEGER                                               :: Npart,del_part
    INTEGER,DIMENSION(:),POINTER                          :: list_del_parts
    REAL(MK)                                              :: almostone
    CLASS(DTYPE(ppm_t_part_prop)_),POINTER                :: prop => NULL()

    start_subroutine("particles_apply_bc")

    !-----------------------------------------------------------------
    !  checks
    !-----------------------------------------------------------------
    IF (.NOT.ASSOCIATED(Pc%xp)) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    topoid = Pc%active_topoid
    topo=>ppm_topo(topoid)%t
    xp=>Pc%xp
    Npart = Pc%Npart
    almostone = 1._MK - EPSILON(1._MK)

    !-----------------------------------------------------------------
    !  Move particles if needed
    !-----------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
    min_phys = topo%min_physs
    max_phys = topo%max_physs
#elif __KIND == __DOUBLE_PRECISION
    min_phys = topo%min_physd
    max_phys = topo%max_physd
#endif
    len_phys=max_phys-min_phys

    del_part = 0
    DO di=1,ppm_dim
        IF (topo%bcdef(di) .EQ. ppm_param_bcdef_periodic) THEN
            DO ip=1,Npart
                IF (xp(di,ip) .EQ. max_phys(di)) &
                    xp(di,ip) = xp(di,ip) - len_phys(di)*almostone
                IF (xp(di,ip) .GT. max_phys(di)) &
                    xp(di,ip) = xp(di,ip) - len_phys(di)
                IF (xp(di,ip) .LT. min_phys(di)) &
                    xp(di,ip) = xp(di,ip) + len_phys(di)
            ENDDO
        ELSE IF (topo%bcdef(di) .EQ. ppm_param_bcdef_freespace) THEN
            !delete particles that have crossed the boundary
            DO ip=Npart,1,-1
                IF (xp(di,ip).GE.max_phys(di).OR.xp(di,ip).LT.min_phys(di)) THEN
                    del_part = del_part+1
                ENDIF
            ENDDO
        ELSE
            info = ppm_error_error
            CALL ppm_error(999,caller,   &
                & 'this type of BC is not implemented/tested in this version',&
                &  __LINE__,info)
            GOTO 9999
        ENDIF
    ENDDO

    IF (del_part .GT. 0 ) THEN
        ldc(1) = del_part
        CALL ppm_alloc(list_del_parts,ldc,ppm_param_alloc_fit,info)
        del_part = 0
        DO di=1,ppm_dim
            IF (topo%bcdef(di) .EQ. ppm_param_bcdef_freespace) THEN
                DO ip=Npart,1,-1
                    IF (xp(di,ip).GE.max_phys(di).OR.xp(di,ip).LT.min_phys(di)) THEN
                        del_part = del_part+1
                        list_del_parts(del_part)=ip
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
        CALL Pc%del_parts(list_del_parts,del_part,info)
        IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,caller,   &
                & 'could not delete particles',&
                &  __LINE__,info)
            GOTO 9999
        ENDIF
        DEALLOCATE(list_del_parts)
    ENDIF

    ! Update states
    Pc%Npart = Npart
    ! Particles are no longer on the right processors
    Pc%flags(ppm_part_partial) = .FALSE.
    ! But they are now all inside the computational domain
    Pc%flags(ppm_part_areinside) = .TRUE.
    ! Dangereous to use the ghosts
    Pc%flags(ppm_part_ghosts) = .FALSE.
    ! ghosts values for properties are also dangerous to use
    prop => Pc%props%begin()
    DO WHILE (ASSOCIATED(prop))
        prop%flags(ppm_ppt_ghosts) = .FALSE.
        prop => Pc%props%next()
    ENDDO


    !-----------------------------------------------------------------------
    ! Finalize
    !-----------------------------------------------------------------------
    xp => NULL()
    topo => NULL()

    end_subroutine()
END SUBROUTINE DTYPE(part_apply_bc)

SUBROUTINE DTYPE(part_move)(Pc,disp,info)

    !!!  Move all particles according to some displacement field
    !!!  The size of disp must match the size of xp

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_particles))                         :: Pc
    !!! Data structure containing the particles
    REAL(MK), DIMENSION(:,:), POINTER,  INTENT(IN   )     :: disp
    !!! Data structure containing the particles
    INTEGER,                            INTENT(  OUT)     :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    INTEGER                                   :: ip
    REAL(MK),DIMENSION(:,:),       POINTER    :: xp=>NULL()
    CLASS(DTYPE(ppm_t_part_prop)_),POINTER    :: prop => NULL()
    CLASS(ppm_t_operator_discr_),  POINTER    :: op => NULL()

    start_subroutine("particle_move")

    !-----------------------------------------------------------------
    !  checks
    !-----------------------------------------------------------------
    IF (.NOT.ASSOCIATED(Pc%xp)) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF


    CALL Pc%get_xp(xp)
    FORALL (ip=1:Pc%Npart) &
            xp(1:ppm_dim,ip) = xp(1:ppm_dim,ip) + disp(1:ppm_dim,ip)
    CALL Pc%set_xp(xp)

    !-----------------------------------------------------------------
    !  update states
    !-----------------------------------------------------------------
    Pc%flags(ppm_part_ghosts) = .FALSE.

    prop => Pc%props%begin()
    DO WHILE (ASSOCIATED(prop))
        prop%flags(ppm_ppt_ghosts) = .FALSE.
        prop => Pc%props%next()
    ENDDO

    IF (ASSOCIATED(Pc%ops)) THEN
        op => Pc%ops%begin()
        DO WHILE (ASSOCIATED(op))
            op%flags(ppm_ops_iscomputed) = .FALSE.
            op => Pc%ops%next()
        ENDDO
    ENDIF

    Pc%flags(ppm_part_partial) = .FALSE.
    Pc%flags(ppm_part_cartesian) = .FALSE.

    end_subroutine()
END SUBROUTINE DTYPE(part_move)

SUBROUTINE DTYPE(part_neighlist)(Pc,info,nlid,lstore,incl_ghosts,knn)
    !-----------------------------------------------------------------
    !  Neighbor lists for particles
    !-----------------------------------------------------------------
    !  Assumptions:
    ! * Particles positions need to have been mapped onto the topology
    ! * Ghost positions have been computed
    !
    USE ppm_module_neighlist
#ifdef __WITH_CNL
    USE ppm_module_cnl
#endif
    USE ppm_module_inl_vlist
#ifdef __WITH_KDTREE
    USE ppm_module_inl_k_vlist
    USE ppm_module_kdtree
#endif
    
#ifdef __MPI
    INCLUDE "mpif.h"
#endif

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    CLASS(DTYPE(ppm_t_particles))                          :: Pc
    DEFINE_MK()
    !!! Data structure containing the particles
    INTEGER,                            INTENT(  OUT)      :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    !  Optional arguments
    !-------------------------------------------------------------------------
    INTEGER, OPTIONAL,                  INTENT(INOUT)      :: nlid
    !!! which neighbour list are we computing. Default is 1
    LOGICAL, OPTIONAL,                  INTENT(IN   )      :: lstore
    !!! store verlet lists
    LOGICAL, OPTIONAL,                  INTENT(IN   )      :: incl_ghosts
    !!! if true, then verlet lists are computed for all particles, incl. ghosts.
    !!! Default is false.
    INTEGER, OPTIONAL,                  INTENT(IN   )      :: knn
    !!! if present, neighbour lists are constructed such that each particle
    !!! has at least knn neighbours.
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    INTEGER                                   :: op_id,np_target,i
    INTEGER                                   :: ip,ineigh
    !!! index variable
    LOGICAL                                   :: symmetry
    !!! backward compatibility
    LOGICAL                                   :: ensure_knn
    !!! uses a neighbour-finding algorithm that finds enough neighbours
    REAL(MK),DIMENSION(2*ppm_dim):: ghostlayer
    !!!
    REAL(KIND(1.D0))                          :: t1,t2
    TYPE(ppm_t_topo), POINTER                 :: topo => NULL()
#ifdef __WITH_KDTREE
    TYPE(DTYPE(kdtree2)),POINTER              :: tree => NULL()
    TYPE(DTYPE(kdtree2_result)),ALLOCATABLE   :: results(:)
#endif
    INTEGER                                   :: neigh_id,topoid
    INTEGER                                   :: nneighmin,nneighmax
    CLASS(DTYPE(ppm_t_neighlist)_), POINTER   :: Nlist => NULL()
    REAL(MK)                                  :: skin

    REAL(MK),DIMENSION(:),POINTER             :: rcp  => NULL()
    CLASS(ppm_t_operator_discr_),POINTER      :: op => NULL()

    start_subroutine("part_comp_neighlist")

    !-----------------------------------------------------------------
    !  Checks
    !-----------------------------------------------------------------
    IF (.NOT.ASSOCIATED(Pc%xp)) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (.NOT.Pc%flags(ppm_part_partial)) THEN
        !Particles have not been mapped onto this topology
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            &  'Do a partial/global mapping before',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (.NOT.Pc%flags(ppm_part_ghosts)) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            &  'Ghosts have not been updated. They are needed for neighlists',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (Pc%neighs%vec_size.LE.0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            &  'neighlist DS not allocated. Call create_neighlist() first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    !checks that the data structure for the neighbour list has been
    !defined already
    IF (PRESENT(nlid)) THEN
        neigh_id = nlid
        IF (.NOT. Pc%neighs%exists(neigh_id)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,   &
                &  'invalid id for neighbour list.',&
                &  __LINE__,info)
            GOTO 9999
        ENDIF
    ELSE
        neigh_id = 1
        IF (Pc%neighs%vec(neigh_id)%t%P_id.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,   &
        &  'We assumed default neighlist, but the P_id is not default',&
                &  __LINE__,info)
            GOTO 9999
        ENDIF
    ENDIF

    !define an alias for the neighbour list DS
    Nlist => Pc%neighs%vec(neigh_id)%t
    IF (.NOT.ASSOCIATED(Nlist)) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            &  'neighbour list not allocated. ',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    !check that we have a cutoff radius
    SELECT TYPE (Pc)
    CLASS IS (DTYPE(ppm_t_sop))
        IF (Pc%rcp_id.LE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,   &
                &  'cutoff radii for adaptive particles have not been defined',&
                &  __LINE__,info)
            GOTO 9999
        ENDIF
        CALL Pc%get(rcp,Pc%rcp_id,with_ghosts=.TRUE.)
        IF (.NOT.ASSOCIATED(rcp)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,   &
                &  'DS for cutoff radii is not associated',&
                &  __LINE__,info)
            GOTO 9999
        ENDIF
    CLASS DEFAULT
        IF (Nlist%cutoff.LE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,   &
            &  'cutoff is negative or zero - do we really want neighbour lists?',&
                &  __LINE__,info)
            GOTO 9999
        ENDIF
    END SELECT

    IF (Nlist%isymm.EQ.1) THEN
        symmetry =.TRUE.
    ELSE
        symmetry = .FALSE.
    ENDIF
    IF (PRESENT(knn)) THEN
        ensure_knn = .TRUE.
    ELSE
        ensure_knn = .FALSE.
    ENDIF
    skin = Nlist%skin
    topoid = Pc%active_topoid

    do_something: IF (Nlist%uptodate .OR. Pc%Npart.EQ.0) THEN
        !neighbor lists are already up-to-date, or no particles on this proc
        !nothing to do
        IF (Nlist%uptodate) THEN
            info = ppm_error_notice
            CALL ppm_error(999,caller,   &
                &  'neighlists are already up-to-date, NOTHING to do',&
                &  __LINE__,info)
            info = 0
        ELSE
            Nlist%nneighmin = 0
            Nlist%nneighmax = 0
        ENDIF
    ELSE
        !hack to build (potentially incomplete) neighbour lists even 
        !for ghost particles
        np_target = Pc%Npart
        IF (PRESENT(incl_ghosts)) THEN
            IF (incl_ghosts) THEN
                np_target = Pc%Mpart
                topo => ppm_topo(topoid)%t
                IF (MK.EQ.ppm_kind_single) THEN
                    topo%min_subs(:,:) = topo%min_subs(:,:) - topo%ghostsizes
                    topo%max_subs(:,:) = topo%max_subs(:,:) + topo%ghostsizes
                ELSE IF (MK.EQ.ppm_kind_double) THEN
                    topo%min_subd(:,:) = topo%min_subd(:,:) - topo%ghostsized
                    topo%max_subd(:,:) = topo%max_subd(:,:) + topo%ghostsized
                ENDIF
            ENDIF
        ENDIF

        IF (ensure_knn) THEN
#ifdef __WITH_KDTREE

                Pc%stats%nb_kdtree = Pc%stats%nb_kdtree+1
#ifdef __MPI
                t1 = MPI_WTIME(info)
#endif

                tree => kdtree2_create(Pc%xp(1:ppm_dim,1:Pc%Mpart),&
                    sort=.true.,rearrange=.true.)
                allocate(results(knn+1))
                ldc(1) = knn
                ldc(2) = Pc%Npart
                CALL ppm_alloc(Nlist%vlist,ldc,ppm_param_alloc_grow,info)
                IF (info .NE. 0) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_alloc,caller,   &
                        &            'failed to allocate vlist',__LINE__,info)
                    GOTO 9999
                ENDIF
                ldc(1) = Pc%Npart
                CALL ppm_alloc(Nlist%nvlist,ldc,ppm_param_alloc_grow,info)
                IF (info .NE. 0) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_alloc,caller,   &
                        &            'failed to allocate vlist',__LINE__,info)
                    GOTO 9999
                ENDIF
                DO ip=1,Pc%Npart
                    call kdtree2_n_nearest(tp=tree,qv=Pc%xp(1:ppm_dim,ip),&
                        nn=knn+1,results=results)
                    !remove ip from the list
                    ineigh=0
                    DO i=1,knn+1
                        IF(results(i)%idx.ne.ip) THEN
                            ineigh=ineigh+1
                            Nlist%vlist(ineigh,ip)=results(i)%idx
                        ENDIF
                    ENDDO
                    Nlist%nvlist(ip)=knn
                ENDDO
                call kdtree2_destroy(tree)
                deallocate(results)
#ifdef __MPI
                t2 = MPI_WTIME(info)
                Pc%stats%t_kdtree = Pc%stats%t_kdtree+(t2-t1)
#endif
#else
                info = ppm_error_error
                CALL ppm_error(ppm_err_alloc,caller,   &
                    &   'option required the kdtree module.',__LINE__,info)
                GOTO 9999
#endif
!__WITH_KDTREE
        ELSE  

            SELECT TYPE (Pc)
            CLASS IS (DTYPE(ppm_t_sop))

                !FIXME: when adaptive ghost layers are available
                ghostlayer(1:2*ppm_dim)=Pc%ghostlayer

#ifdef __WITH_CNL
                conventionalinl: IF (Pc%conventionalinl) THEN
                    Pc%stats%nb_cinl = Pc%stats%nb_cinl+1

#ifdef __MPI
                    t1 = MPI_WTIME(info)
#endif
                    !HUGLY HACK to make CNL routines work on a topology with
                    !several subdomains
#if   __KIND == __SINGLE_PRECISION
                    CALL cnl_vlist(Pc%xp,&
                        rcp,Pc%Npart,Pc%Mpart,&
                        ppm_topo(topoid)%t%min_subs(:,1)-Pc%ghostlayer,&
                        ppm_topo(topoid)%t%max_subs(:,1)+Pc%ghostlayer,&
                        Pc%nvlist,Pc%vlist,ppm_dim,info)
#elif __KIND == __DOUBLE_PRECISION
                    CALL cnl_vlist(Pc%xp,&
                        rcp,Pc%Npart,Pc%Mpart,&
                        ppm_topo(topoid)%t%min_subd(:,1)-Pc%ghostlayer,&
                        ppm_topo(topoid)%t%max_subd(:,1)+Pc%ghostlayer,&
                        Pc%nvlist,Pc%vlist,ppm_dim,info)
#endif
                    IF (info .NE. 0) THEN
                        info = ppm_error_error
                        CALL ppm_error(ppm_err_sub_failed,caller,&
                            'ppm_cinl_vlist failed',__LINE__,info)
                        GOTO 9999
                    ENDIF
                    !end HUGLY HACK
#ifdef __MPI
                    t2 = MPI_WTIME(info)
                    Pc%stats%t_cinl = Pc%stats%t_cinl + (t2 - t1)
#endif
                ELSE
#endif
!__WITH_CNL
                    Pc%stats%nb_inl = Pc%stats%nb_inl+1

#ifdef __MPI
                    t1 = MPI_WTIME(info)
#endif
                    CALL ppm_inl_vlist(topoid,Pc%xp,np_target,&
                        Pc%Mpart,rcp,skin,symmetry,ghostlayer,info,&
                        Nlist%vlist,Nlist%nvlist)
                    IF (info .NE. 0) THEN
                        info = ppm_error_error
                        CALL ppm_error(ppm_err_sub_failed,caller,&
                            'ppm_inl_vlist failed',__LINE__,info)
                        GOTO 9999
                    ENDIF
#ifdef __MPI
                    t2 = MPI_WTIME(info)
                    Pc%stats%t_inl = Pc%stats%t_inl + (t2 - t1)
#endif
#ifdef __WITH_CNL
                ENDIF conventionalinl
#endif

            CLASS DEFAULT

                Pc%stats%nb_nl = Pc%stats%nb_nl+1

#ifdef __MPI
                t1 = MPI_WTIME(info)
#endif
                CALL ppm_neighlist_vlist(topoid,Pc%xp,Pc%Mpart,&
                    Nlist%cutoff,skin,symmetry,Nlist%vlist,&
                    Nlist%nvlist,info,lstore=lstore)
                IF (info .NE. 0) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_sub_failed,caller,&
                        'ppm_neighlist_vlist failed',__LINE__,info)
                    GOTO 9999
                ENDIF
#ifdef __MPI
                t2 = MPI_WTIME(info)
                Pc%stats%t_nl = Pc%stats%t_nl + (t2 - t1)
#endif
            END SELECT
        ENDIF

        !restore subdomain sizes (revert hack)
        IF (PRESENT(incl_ghosts)) THEN
            IF (incl_ghosts) THEN
                IF (MK.EQ.ppm_kind_single) THEN
                    topo%min_subs(:,:) = topo%min_subs(:,:) + topo%ghostsizes
                    topo%max_subs(:,:) = topo%max_subs(:,:) - topo%ghostsizes
                ELSE IF (MK.EQ.ppm_kind_double) THEN
                    topo%min_subd(:,:) = topo%min_subd(:,:) + topo%ghostsized
                    topo%max_subd(:,:) = topo%max_subd(:,:) - topo%ghostsized
                ENDIF
                topo => NULL()
            ENDIF
        ENDIF

        !-----------------------------------------------------------------------
        !Update state
        !-----------------------------------------------------------------------
        Nlist%uptodate = .TRUE.

        Nlist%nneighmin = MINVAL(Nlist%nvlist(1:Pc%Npart))
        Nlist%nneighmax = MAXVAL(Nlist%nvlist(1:np_target))

        ! DC operators that do not use a xset neighbour list, if they exist, 
        ! are no longer valid (they depend on the neighbour lists)
        IF (ASSOCIATED(Pc%ops)) THEN
            op => Pc%ops%begin()
            DO WHILE (ASSOCIATED(op))
                IF (.NOT.op%flags(ppm_ops_interp)) THEN
                    op%flags(ppm_ops_iscomputed) = .FALSE.
                ENDIF
                op => Pc%ops%next()
            ENDDO
        ENDIF

        
        Nlist => NULL()

        !FIXME
        ! We want to distinguish between "self" neighbour lists
        ! and cross-set ones.
        IF (neigh_id .EQ. 1) THEN
            Pc%flags(ppm_part_neighlists) = .TRUE.
        ENDIF

    ENDIF do_something

    end_subroutine()
END SUBROUTINE DTYPE(part_neighlist)


SUBROUTINE DTYPE(part_set_cutoff)(Pc,cutoff,info,Nlist)
    !!! Set a cutoff radius for a Particle set and update the
    !!! the ghostlayer sizes.
    !!! The cutoff radius concerns the default neighbor list, unless
    !!! specified otherwise.
    !-------------------------------------------------------------------------
    ! Arguments
    !-------------------------------------------------------------------------
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_particles))            :: Pc
    REAL(MK),                 INTENT(IN   )  :: cutoff
    !!! cutoff radius
    INTEGER,                  INTENT(   OUT) :: info
    !!! return status. On success, 0
    CLASS(DTYPE(ppm_t_neighlist)_),OPTIONAL,INTENT(INOUT) :: NList
    !!! Neighbor list for which this cutoff radius
    !!! applies. Default is ppm_param_default_nlID
    
    !-------------------------------------------------------------------------
    ! local variables
    !-------------------------------------------------------------------------
    CLASS(DTYPE(ppm_t_neighlist)_),   POINTER :: nl => NULL()
    INTEGER                                   :: neigh_id 

    start_subroutine("part_set_cutoff")

    !-------------------------------------------------------------------------
    !  Set new cutoff
    !-------------------------------------------------------------------------
    IF (PRESENT(NList)) THEN
        check_true("Pc%neighs%has(NList)",&
            "Neighbour list does not concern this particle set")
        NList%cutoff = cutoff
    ELSE
        nl => Pc%get_neighlist()
        check_associated(nl,"Compute neighbour lists first")
        nl%cutoff = cutoff
    ENDIF

    
    ! Compute ghostlayer sizes
    IF (cutoff.GT.Pc%ghostlayer) THEN
        !If the new cutoff is larger than the current ghostsize
        ! then the new ghostsize is the new cutoff
        Pc%ghostlayer = cutoff
        ! update states
        Pc%flags(ppm_part_ghosts) = .FALSE.
    ELSE IF (cutoff .LT. Pc%ghostlayer) THEN
        !Else, we find the new maximum cutoff radius amongst
        !all existing neighbor lists on this Particle set
        Pc%ghostlayer = 0._mk
        nl => Pc%neighs%begin()
        DO WHILE (ASSOCIATED(nl))
            IF (nl%cutoff .GT. Pc%ghostlayer) THEN
                Pc%ghostlayer = nl%cutoff
            ENDIF
            nl => Pc%neighs%next()
        ENDDO
        !no need to update states: ghosts are still ok.
    ENDIF

    end_subroutine()
END SUBROUTINE DTYPE(part_set_cutoff)

SUBROUTINE DTYPE(part_comp_global_index)(Pc,info)
    !!! Compute a global index for particles
    !!! (Uses MPI communications)
    !-------------------------------------------------------------------------
    ! Arguments
    !-------------------------------------------------------------------------
#ifdef __MPI
    INCLUDE "mpif.h"
#endif

    DEFINE_MK()
    CLASS(DTYPE(ppm_t_particles))            :: Pc
    INTEGER,                  INTENT(   OUT) :: info
    !!! return status. On success, 0

    !-------------------------------------------------------------------------
    ! local variables
    !-------------------------------------------------------------------------
    INTEGER                        :: offset
    INTEGER                        :: i
    INTEGER, DIMENSION(:), POINTER :: wp
    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    start_subroutine("part_comp_global_index")

    IF (.NOT. Pc%flags(ppm_part_global_index)) THEN
        CALL Pc%create_prop(Pc%gi_id,ppm_type_int,info,name="GlobalIndex")
        Pc%flags(ppm_part_global_index) = .TRUE.
    END IF
#ifdef __MPI
    CALL MPI_Scan(Pc%Npart,offset,1,MPI_INTEGER,MPI_SUM,ppm_comm,info)
    offset = offset - Pc%Npart
#else
    offset = 0
#endif
    CALL Pc%get(wp,Pc%gi_id)
    FORALL (i=1:Pc%Npart) wp(i) = offset + i !- 1 !uncomment if index from 0
    CALL Pc%set(wp,Pc%gi_id)

    end_subroutine()
END SUBROUTINE DTYPE(part_comp_global_index)




SUBROUTINE DTYPE(part_map_create)(Pc,id,source_topoid,target_topoid,info)
    !!! Adds a property to an existing particle set
    CLASS(DTYPE(ppm_t_particles))         :: Pc
    INTEGER,                INTENT(  OUT) :: id
    INTEGER,                INTENT(IN   ) :: source_topoid
    INTEGER,                INTENT(IN   ) :: target_topoid
    INTEGER,               INTENT(OUT)    :: info

    INTEGER                               :: vec_size,npart,i
    TYPE(DTYPE(ppm_t_ptr_part_mapping)),DIMENSION(:),POINTER:: vec_tmp => NULL()
    CLASS(DTYPE(ppm_t_part_mapping)_),               POINTER:: map => NULL()
    LOGICAL, DIMENSION(ppm_param_length_pptflags):: flags

    start_subroutine("part_map_create")

    !Generate a new id (we should use templating here...)
    ASSOCIATE (cont => Pc%maps )
        id = 0
        IF (cont%nb.LT.cont%vec_size) THEN
            !there is at least one empty slot in the array
            ! of mapping pointers
            id = id + 1
            DO WHILE (ASSOCIATED(cont%vec(id)%t))
                id = id + 1
            ENDDO
        ELSE
            IF (.NOT. ASSOCIATED(cont%vec)) THEN
                !need to allocate the array of mapping pointers 
                vec_size=20
                ALLOCATE(cont%vec(vec_size),STAT=info)
                IF (info .NE. 0) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_alloc,caller,&
                        'allocating mapping array failed',__LINE__,info)
                    GOTO 9999
                ENDIF
                id = 1
            ELSE
                !need to resize the array of mapping pointers 
                vec_size=MAX(2*cont%vec_size,20)
                ALLOCATE(vec_tmp(vec_size),STAT=info)
                IF (info .NE. 0) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_alloc,caller,&
                        'allocating mapping array failed',__LINE__,info)
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

    IF (.NOT. ASSOCIATED(Pc%maps%vec(id)%t)) THEN
        ALLOCATE(DTYPE(ppm_t_part_mapping)::Pc%maps%vec(id)%t,STAT=info)
            or_fail_alloc("Pc%maps%vec(id)%t")
    ENDIF

    map => Pc%maps%vec(id)%t

    ! Create the mapping
    CALL map%create(source_topoid,target_topoid,info)
        or_fail("map%create")

    end_subroutine()
END SUBROUTINE DTYPE(part_map_create)

SUBROUTINE DTYPE(part_map_destroy)(Pc,id,info)
    !!! Destroy a mapping from an existing particle set
    CLASS(DTYPE(ppm_t_particles))         :: Pc
    INTEGER,                INTENT(INOUT) :: id
    INTEGER,               INTENT(OUT)    :: info

    start_subroutine("part_map_destroy")

    ASSOCIATE (cont => Pc%maps)
        IF (id .LE. 0 .OR. id .GT. cont%vec_size) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,&
                &    'mapping id larger than size of mappings array',&
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


    end_subroutine()
END SUBROUTINE DTYPE(part_map_destroy)


!DESTROY ENTRY
SUBROUTINE DTYPE(neigh_destroy)(neigh,info)
    CLASS(DTYPE(ppm_t_neighlist))      :: neigh
    INTEGER,                                INTENT(  OUT)  :: info
    !!! Returns status, 0 upon success.

    start_subroutine("neigh_destroy")

    IF(ASSOCIATED(neigh%nvlist)) DEALLOCATE(neigh%nvlist,STAT=info)
    IF(ASSOCIATED(neigh%vlist))  DEALLOCATE(neigh%vlist,STAT=info)

    end_subroutine()
END SUBROUTINE DTYPE(neigh_destroy)

FUNCTION DTYPE(has_ghosts)(this,Field) RESULT(res)
    !!! Returns true if the discretization of the field on this
    !!! particle set has its ghosts up-to-date.
    !!! If the Field argument is not present, then the function
    !!! returns whether the particles themselves have their
    !!! ghosts up-to-date.
    CLASS(DTYPE(ppm_t_particles))                  :: this
    CLASS(ppm_t_field_),OPTIONAL                   :: Field
    LOGICAL                                        :: res

    !Local variables
    CLASS(DTYPE(ppm_t_part_prop)_),POINTER         :: prop => NULL()

    start_function("has_ghosts")

    IF (PRESENT(Field)) THEN
        CALL this%get_prop(Field,prop,info)
            or_fail("this field is not discretized on this particle set")
        res = prop%flags(ppm_ppt_ghosts)
    ELSE
        res = this%flags(ppm_part_ghosts)
    ENDIF

    end_function()

END FUNCTION

SUBROUTINE DTYPE(get_prop)(this,Field,prop,info)
    CLASS(DTYPE(ppm_t_particles))                  :: this
    CLASS(ppm_t_field_),TARGET,             INTENT(IN   )  :: Field
    CLASS(DTYPE(ppm_t_part_prop)_),POINTER, INTENT(  OUT)  :: prop
    INTEGER,                                INTENT(  OUT)  :: info

    start_subroutine("get_prop")
    
    prop => this%props%begin()
    loop: DO WHILE(ASSOCIATED(prop))
        IF (ASSOCIATED(prop%field_ptr,Field)) THEN
            EXIT loop
        ENDIF
        prop => this%props%next()
    ENDDO loop

    check_associated(prop,&
        "Failed to find a discretization of this field for this particle set")

    end_subroutine()
END SUBROUTINE

#undef DEFINE_MK


