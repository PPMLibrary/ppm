minclude ppm_create_collection_procedures(DTYPE(part_prop),DTYPE(part_prop)_)
minclude ppm_create_collection_procedures(DTYPE(neighlist),DTYPE(neighlist)_)
minclude ppm_create_collection_procedures(DTYPE(particles),DTYPE(particles)_)
!minclude ppm_create_collection_procedures(DTYPE(ppm_t_sop))

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
SUBROUTINE DTYPE(get_xp)(this,xp,info,with_ghosts)
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_particles))                :: this
    REAL(MK),DIMENSION(:,:),POINTER,INTENT(OUT)  :: xp
    INTEGER,                        INTENT(OUT)  :: info
    LOGICAL,OPTIONAL                             :: with_ghosts

    start_subroutine("get_xp")

    check_associated("this%xp")

    IF (PRESENT(with_ghosts)) THEN
        IF (with_ghosts) THEN
            IF (this%flags(ppm_part_ghosts)) THEN
                xp => this%xp(1:ppm_dim,1:this%Mpart)
            ELSE
                write(cbuf,*) 'WARNING: tried to get xp with ghosts ',&
                    'when ghosts are not up-to-date'
                CALL ppm_write(ppm_rank,'get_xp',cbuf,info)
                xp => NULL()
            ENDIF
            RETURN
        ENDIF
    ENDIF

    xp => this%xp(1:ppm_dim,1:this%Npart)

    end_subroutine()
END SUBROUTINE DTYPE(get_xp)

SUBROUTINE DTYPE(set_xp)(this,xp,info,read_only,ghosts_ok)
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_particles))    :: this
    LOGICAL,OPTIONAL                 :: read_only
    LOGICAL,OPTIONAL                 :: ghosts_ok
    REAL(MK),DIMENSION(:,:),POINTER  :: xp
    INTEGER,            INTENT(OUT)  :: info

    INTEGER                          :: i
    CLASS(DTYPE(ppm_t_part_prop)_), POINTER :: prop => NULL()
    CLASS(DTYPE(ppm_t_neighlist)_), POINTER :: nl => NULL()
    CLASS(ppm_t_operator_discr_),   POINTER :: op => NULL()

    start_subroutine("set_xp") 

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

    this%flags(ppm_part_areinside) = .FALSE.
    this%flags(ppm_part_partial) = .FALSE.
    this%flags(ppm_part_ghosts) = .FALSE.
    this%flags(ppm_part_cartesian) = .FALSE.

    prop => this%props%begin()
    DO WHILE (ASSOCIATED(prop))
        prop%flags(ppm_ppt_ghosts) = .FALSE.
        prop%flags(ppm_ppt_partial) = .FALSE.
        prop => this%props%next()
    ENDDO

    nl => this%neighs%begin()
    DO WHILE (ASSOCIATED(nl))
        nl%uptodate = .FALSE.
        nl => this%neighs%next()
    ENDDO

    IF (ASSOCIATED(this%ops)) THEN
        op => this%ops%begin()
        DO WHILE (ASSOCIATED(op))
            op%flags(ppm_ops_iscomputed) = .FALSE.
            op => this%ops%next()
        ENDDO
    ENDIF

    xp => NULL()

    end_subroutine()
END SUBROUTINE DTYPE(set_xp)

SUBROUTINE DTYPE(part_prop_create)(this,info,field,part_prop,discr_data,&
        dtype,name,lda,zero,with_ghosts)
    !!! Adds a property to an existing particle set
    CLASS(DTYPE(ppm_t_particles))         :: this
    INTEGER,               INTENT(OUT)    :: info
    CLASS(ppm_t_field_),OPTIONAL,INTENT(IN   ) :: field
    CLASS(DTYPE(ppm_t_part_prop)_),OPTIONAL,POINTER,INTENT(OUT):: part_prop
    !!! Pointer to the ppm_t_part_prop object for that property
    CLASS(ppm_t_discr_data),OPTIONAL,POINTER,       INTENT(OUT):: discr_data
    !!! Pointer to the ppm_t_discr_data object for that property
    INTEGER, OPTIONAL,      INTENT(IN   ) :: dtype
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN ) :: name
    INTEGER, OPTIONAL,      INTENT(IN   ) :: lda
    !!! name to this property
    LOGICAL, OPTIONAL                     :: zero
    !!! if true, then initialise the data to zero
    LOGICAL, OPTIONAL                     :: with_ghosts
    !!! if true, then allocate with Mpart instead of the default size of Npart

    CLASS(DTYPE(ppm_t_part_prop)_),POINTER :: prop => NULL()
    CHARACTER(LEN=ppm_char)               :: name2
    INTEGER                               :: lda2,vec_size,npart,i,datatype
    LOGICAL, DIMENSION(ppm_param_length_pptflags):: flags

    start_subroutine("particle_prop_create")

    IF (PRESENT(field)) THEN
        lda2 = field%lda
        name2 = field%name
        datatype = field%data_type
    ELSE
        lda2 = 1
        IF (PRESENT(lda)) THEN
            IF (lda.GE.2) THEN
                lda2 = lda
            ENDIF
        ENDIF
        IF (PRESENT(name)) THEN
            name2 = name
        ELSE
            name2 = "default_ppt_name"
        ENDIF
        IF (PRESENt(dtype)) THEN
            datatype = dtype
        ELSE
            datatype = ppm_type_real
        ENDIF
    ENDIF

    npart = this%Npart
    flags = .FALSE.
    IF (PRESENT(with_ghosts)) THEN
        IF (with_ghosts) THEN
            IF (this%flags(ppm_part_ghosts)) THEN
                npart = this%Mpart
                flags(ppm_ppt_ghosts) = .TRUE.
            ELSE
        fail("trying to init property for ghosts when ghosts are not computed")
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

    IF (PRESENT(part_prop)) THEN
        part_prop => prop
    ENDIF
    IF (PRESENT(discr_data)) THEN
        discr_data => prop
    ENDIF

    CALL this%props%push(prop,info)
        or_fail("pushing new property into collection failed")

    end_subroutine()
END SUBROUTINE DTYPE(part_prop_create)

SUBROUTINE DTYPE(part_prop_destroy)(this,prop,info)
    !!! Destroy a property from an existing particle set
    CLASS(DTYPE(ppm_t_particles))         :: this
    CLASS(ppm_t_discr_data),INTENT(INOUT) :: prop
    INTEGER,               INTENT(OUT)    :: info

    start_subroutine("part_prop_destroy")

    check_associated("this%props",&
        "Particle set does not contain any discretized data")

    SELECT TYPE(prop)
    CLASS IS (DTYPE(ppm_t_part_prop))
        CALL this%props%remove(info,prop)
            or_fail("could not remove property from its container")
    CLASS DEFAULT
        fail("discretization data has to be of class ppm_t_part_prop")
    END SELECT

    

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


SUBROUTINE DTYPE(part_neigh_create)(this,Part_src,info,&
        name,skin,symmetry,cutoff,Nlist)
    !!! Create a data structure to store a neighbour list
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_particles))                      :: this
    CLASS(DTYPE(ppm_t_particles)_),TARGET, INTENT(IN)  :: Part_src
    !!! Particle set to which the neighbours belong (can be the same as this)
    INTEGER,                               INTENT(OUT) :: info
    CHARACTER(LEN=*) , OPTIONAL                                 :: name
    !!! name of this neighbour list
    REAL(MK), OPTIONAL                                          :: skin
    REAL(MK), OPTIONAL                                          :: cutoff
    LOGICAL, OPTIONAL                                           :: symmetry    
    CLASS(DTYPE(ppm_t_neighlist)_),POINTER,OPTIONAL,INTENT(OUT) :: Nlist
    !!! returns a pointer to the newly created verlet list

    INTEGER                                :: vec_size,i
    CLASS(DTYPE(ppm_t_neighlist)_),POINTER :: Nl=>NULL()
    REAL(MK), DIMENSION(:),        POINTER :: rcp => NULL()

    start_subroutine("particle_neigh_create")

    ! Create the neighbour list
    ALLOCATE(DTYPE(ppm_t_neighlist)::Nl,STAT=info)
        or_fail_alloc("Nl")

    IF (PRESENT(name)) THEN
        Nl%name = name
    ELSE
        WRITE(Nl%name,*) 'Nl',TRIM(ADJUSTL(this%name)),'_',&
            TRIM(ADJUSTL(Part_src%name))
    ENDIF

    check_associated("Part_src%xp","Invalid particle set Part_src")

    Nl%Part => Part_src

    SELECT TYPE(this)
    TYPE IS (DTYPE(ppm_t_sop))
        ASSOCIATE (ghosts => this%flags(ppm_part_ghosts))
            IF (.NOT.ASSOCIATED(this%rcp)) THEN

                CALL this%create_prop(info,part_prop=this%rcp,&
                    dtype=ppm_type_real,name='rcp',with_ghosts=ghosts) 
                    or_fail("Creating property for rcp failed")
            ENDIF

            CALL this%get(this%rcp,rcp,info,with_ghosts=ghosts)
                or_fail("Cannot access this%rcp")

            IF (PRESENT(cutoff)) THEN
                rcp = cutoff
            ELSE
                rcp = this%ghostlayer
            ENDIF
            CALL this%set(this%rcp,rcp,info,ghosts_ok=ghosts)
        END ASSOCIATE
        Nl%cutoff = -1._MK 
        !this field should not be used with adaptive particles
    CLASS DEFAULT
        IF (PRESENT(cutoff)) THEN
            Nl%cutoff = cutoff
        ELSE
            Nl%cutoff = this%ghostlayer
        ENDIF
    END SELECT

    IF (PRESENT(skin)) THEN
        Nl%skin = skin
    ELSE
        Nl%skin = 0._mk
    ENDIF

    IF (PRESENT(symmetry)) THEN
        IF (symmetry) THEN
            Nl%isymm = 1
        ELSE
            Nl%isymm = 0
        ENDIF
    ELSE
        Nl%isymm = 0
    ENDIF

    Nl%uptodate = .FALSE.
    Nl%nneighmin = 0
    Nl%nneighmax = 0

    !returning a pointer to the neighbour list, before it is pushed
    ! into the collection.
    IF (PRESENT(Nlist)) Nlist => Nl

    CALL this%neighs%push(Nl,info)
        or_fail("pushing new neighbour list into collection failed")

    end_subroutine()
END SUBROUTINE DTYPE(part_neigh_create)

SUBROUTINE DTYPE(part_neigh_destroy)(this,Nlist,info)
    !!! Destroy a property from an existing particle set
    CLASS(DTYPE(ppm_t_particles))                        :: this
    CLASS(DTYPE(ppm_t_neighlist)_),POINTER,INTENT(INOUT) :: Nlist
    INTEGER,                               INTENT(OUT)   :: info

    start_subroutine("part_neigh_destroy")

    CALL this%neighs%remove(info,Nlist)
        or_fail("could not remove Nlist from its collection")

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

    Pc%gi => NULL()

    IF (.NOT. ASSOCIATED(Pc%field_ptr)) THEN
        ALLOCATE(Pc%field_ptr,STAT=info)
            or_fail_alloc("Pc%field_ptr")
    ELSE
        fail("Pc%field_ptr was already associated. Use destroy() first")
    ENDIF

    SELECT TYPE(Pc)
    CLASS IS (DTYPE(ppm_t_sop))
        !-----------------------------------------------------------------
        !  Initialize fields of the extended SOP type
        !-----------------------------------------------------------------
        ! Particles are by default not adaptive
        Pc%adaptive = .FALSE.
        Pc%adapt_wp => NULL()
        Pc%rcp => NULL()
        Pc%D => NULL()
        Pc%Dtilde => NULL()
        ! Particles do not represent a level-set function
        Pc%level_set = .FALSE.
        Pc%level => NULL()
        !        Pc%level_old_id = 0
        Pc%level_grad => NULL()
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


SUBROUTINE DTYPE(part_map)(Pc,info,debug,global,topoid)

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
END SUBROUTINE DTYPE(part_map)

SUBROUTINE DTYPE(part_map_ghosts)(Pc,info,ghostsize,debug)

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
    CLASS(ppm_t_discr_data), POINTER          :: prop => NULL()

    start_subroutine("part_map_ghosts")

    dbg = .FALSE.
    IF (PRESENT(debug)) dbg=debug
    skip_ghost_get = .FALSE.
    skip_send = .TRUE.
    !we must not call ppm_map_part_send unless ppm_map_part_push (or ghost_get)
    ! has been called (in which case, skip_send is set to FALSE)

    !-----------------------------------------------------------------
    !  Checks
    !-----------------------------------------------------------------
    !check that particles are allocated
    check_associated("Pc%xp",&
        "Particles structure had not been defined. Call allocate first")
    !check that particles are mapped onto this topology
    check_true("Pc%flags(ppm_part_partial)",&
        "Do a partial/global mapping before doing a ghost mapping")
    !check that particles are inside the domain
    check_true("Pc%flags(ppm_part_areinside)",&
        "some particles may be outside the domain. Apply BC first")

    topoid = Pc%active_topoid
    topo=>ppm_topo(topoid)%t

    cutoff = Pc%ghostlayer
    IF (PRESENT(ghostsize)) THEN
        IF (ghostsize .LT. cutoff) THEN
            fail("using ghostsize < cutoff+skin. Increase ghostsize.")
        ELSE
            cutoff = ghostsize
        ENDIF
    ENDIF

#if   __KIND == __SINGLE_PRECISION
    IF (cutoff .GT. topo%ghostsizes) THEN
        fail("ghostsize of topology may be smaller than that of particles")
    ENDIF
#elif   __KIND == __DOUBLE_PRECISION
    IF (cutoff .GT. topo%ghostsized) THEN
        stdout("cutoff for ghost mapping       = ",cutoff)
        stdout("cutoff used to create topology = ",'topo%ghostsized')
        fail("ghostsize of topology may be smaller than that of particles")
    ENDIF
#endif
    IF (cutoff .GT. 0._MK) THEN
        IF (Pc%flags(ppm_part_ghosts)) THEN
            IF (dbg) THEN
                stdout("ghosts have already been updated")
            ENDIF

            IF (ppm_map_type_isactive(ppm_param_map_ghost_get)) THEN
                IF (dbg) THEN
                    stdout("we skip the ghost_get and go straight to",&
                        " push/send/pop")
                ENDIF
                skip_ghost_get = .TRUE.
            ENDIF
        ENDIF

        IF (.NOT.skip_ghost_get) THEN
            Pc%stats%nb_ghost_get = Pc%stats%nb_ghost_get + 1
#ifdef __MPI
            t1 = MPI_WTIME(info)
#endif
            CALL ppm_map_part_ghost_get(topoid,Pc%xp,ppm_dim,&
                Pc%Npart,Pc%isymm,cutoff,info)
                or_fail("ppm_map_part_ghost_get failed")
#ifdef __MPI
            t2 = MPI_WTIME(info)
            Pc%stats%t_ghost_get = Pc%stats%t_ghost_get + (t2-t1)
#endif
            skip_send = .FALSE.
        ELSE
            IF(dbg) THEN
                stdout("skipping ghost-get")
            ENDIF
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

                        IF(dbg) THEN
                            stdout("pushing property ",'Pc%props%iter_id',&
                                'TRIM(prop%name)')
                        ENDIF
                        Pc%stats%nb_ghost_push = Pc%stats%nb_ghost_push + 1
#ifdef __MPI
                        t1 = MPI_WTIME(info)
#endif
                        CALL Pc%map_part_push_legacy(Pc%props%iter_id,info)
                            or_fail("map_part_push")
#ifdef __MPI
                        t2 = MPI_WTIME(info)
                        Pc%stats%t_ghost_push = &
                            Pc%stats%t_ghost_push + (t2-t1)
#endif
                        skip_send = .FALSE.
                    ELSE
                        stdout("pushing property ",&
                            'Pc%props%iter_id','TRIM(prop%name)')
                        fail("getting ghosts for a property thats not mapped")
                    ENDIF
                ENDIF
            ENDIF
            prop => Pc%props%next()
        ENDDO

        IF (.NOT. skip_send) THEN
            CALL ppm_map_part_send(Pc%Npart,Pc%Mpart,info)
                or_fail("ppm_map_part_send")

            prop => Pc%props%last()
            DO WHILE (ASSOCIATED(prop))

                IF (prop%flags(ppm_ppt_map_ghosts)) THEN
                    IF (.NOT.prop%flags(ppm_ppt_ghosts)) THEN
                        IF (prop%flags(ppm_ppt_partial)) THEN

                            IF(dbg) THEN
                                stdout("popping property ",'Pc%props%iter_id',&
                                    'TRIM(prop%name)')
                            ENDIF
                            CALL Pc%map_part_pop_legacy(Pc%props%iter_id,&
                                Pc%Mpart,info)
                                or_fail("map_part_pop")
                            prop%flags(ppm_ppt_ghosts) = .TRUE.
                        ENDIF
                    ENDIF
                ENDIF
                prop => Pc%props%prev()
            ENDDO

            IF (.NOT.skip_ghost_get) THEN
                IF(dbg) THEN
                    stdout("popping xp")
                ENDIF
                CALL ppm_map_part_pop(Pc%xp,ppm_dim,Pc%Npart,&
                    Pc%Mpart,info)
                    or_fail("map_part_pop")
            ENDIF
        ENDIF !.NOT.skip_send

    ELSE ! if cutoff .le. 0

        IF(dbg) THEN
            stdout("cutoff = 0, nothing to do")
            stdout("setting all %has_ghost properties to true")
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
END SUBROUTINE DTYPE(part_map_ghosts)

SUBROUTINE DTYPE(part_map_ghost_get)(Pc,info,ghostsize,debug)
    !!! Push ghost particles properties into the send buffer.
    !!! If ghost positions are not up-to-date, they are re-computed
    !!! with a call to ghost_get. Otherwise ghost_get is skipped.
    !!! If the Field argument is present, only the discretization of that
    !!! Field on this particle set is pushed. Otherwise, all the properties
    !!! (that are not already up-to-date) are pushed, except those that
    !!! have the ppm_ppt_map_ghosts flag set to false.
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

    start_subroutine("part_map_ghost_get")

    dbg = .FALSE.
    IF (PRESENT(debug)) dbg=debug

    !-----------------------------------------------------------------
    !  Checks
    !-----------------------------------------------------------------
    !check that particles are allocated
    check_associated("Pc%xp",&
        "Particles structure had not been defined. Call allocate first")
    !check that particles are mapped onto this topology
    check_true("Pc%flags(ppm_part_partial)",&
        "Partial/global mapping required before doing a ghost mapping")
    !check that particles are inside the domain
    check_true("Pc%flags(ppm_part_areinside)",&
        "Some particles may be outside the domain. Apply BC first")

    topoid = Pc%active_topoid
    topo=>ppm_topo(topoid)%t

    cutoff = Pc%ghostlayer
    IF (PRESENT(ghostsize)) THEN
        IF (ghostsize .LT. cutoff) THEN
            fail("using ghostsize < cutoff+skin. Increase ghostsize.")
        ELSE
            cutoff = ghostsize
        ENDIF
    ENDIF

#if   __KIND == __SINGLE_PRECISION
    IF (cutoff .GT. topo%ghostsizes) THEN
        stdout("cutoff for ghost mapping       = ",cutoff)
        stdout("cutoff used to create topology = ",'topo%ghostsizes')
        fail("ghostsize of topology may be smaller than that of particles")
    ENDIF
#elif   __KIND == __DOUBLE_PRECISION
    IF (cutoff .GT. topo%ghostsized) THEN
        stdout("cutoff for ghost mapping       = ",cutoff)
        stdout("cutoff used to create topology = ",'topo%ghostsized')
        fail("ghostsize of topology may be smaller than that of particles")
    ENDIF
#endif
    IF (cutoff .GT. 0._MK) THEN
        IF (Pc%flags(ppm_part_ghosts)) THEN
            IF (dbg) THEN
                stdout("ghosts have already been updated")
            ENDIF
        ENDIF

        IF (Pc%flags(ppm_part_ghosts)) THEN
            IF (dbg) THEN
                stdout("ghosts have already been updated")
            ENDIF

            IF (ppm_map_type_isactive(ppm_param_map_ghost_get)) THEN
                IF (dbg) THEN
                    stdout("we skip the ghost_get and go straight to",&
                        " push/send/pop")
                ENDIF

                Pc%stats%nb_ghost_get = Pc%stats%nb_ghost_get + 1
#ifdef __MPI
                t1 = MPI_WTIME(info)
#endif
                CALL ppm_map_part_ghost_get(topoid,Pc%xp,ppm_dim,&
                    Pc%Npart,Pc%isymm,cutoff,info)
                    or_fail("ppm_map_part_ghost_get failed")
#ifdef __MPI
                t2 = MPI_WTIME(info)
                Pc%stats%t_ghost_get = Pc%stats%t_ghost_get + (t2-t1)
#endif
            ELSE
                IF(dbg) THEN
                    stdout("skipping ghost-get")
                ENDIF
            ENDIF
        ELSE
            IF(dbg) THEN
        stdout("flags(ppm_part_ghosts) set to false. Skipping ghost-get")
            ENDIF
        ENDIF
    ELSE ! if cutoff .le. 0
        IF(dbg) THEN
            stdout("cutoff = 0, nothing to do")
        ENDIF
    ENDIF

    end_subroutine()
END SUBROUTINE DTYPE(part_map_ghost_get)

SUBROUTINE DTYPE(part_map_ghost_push)(Pc,info,Field,ghostsize,debug)
    !!! Push ghost particles properties into the send buffer.
    !!! If ghost positions are not up-to-date, they are re-computed
    !!! with a call to ghost_get. Otherwise ghost_get is skipped.
    !!! If the Field argument is present, only the discretization of that
    !!! Field on this particle set is pushed. Otherwise, all the properties
    !!! (that are not already up-to-date) are pushed, except those that
    !!! have the ppm_ppt_map_ghosts flag set to false.
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
    CLASS(ppm_t_field_), OPTIONAL                       :: Field
    !!! Push only this field. Default is to push all of them.
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
    CLASS(ppm_t_discr_data), POINTER          :: prop => NULL()

    start_subroutine("part_map_ghost_push")

    dbg = .FALSE.
    IF (PRESENT(debug)) dbg=debug

    !-----------------------------------------------------------------
    !  Checks
    !-----------------------------------------------------------------
    !check that particles are allocated
    check_associated("Pc%xp",&
        "Particles structure had not been defined. Call allocate first")
    !check that particles are mapped onto this topology
    check_true("Pc%flags(ppm_part_partial)",&
        "Do a partial/global mapping before doing a ghost mapping")
    !check that particles are inside the domain
    check_true("Pc%flags(ppm_part_areinside)",&
        "some particles may be outside the domain. Apply BC first")

    topoid = Pc%active_topoid
    topo=>ppm_topo(topoid)%t

    cutoff = Pc%ghostlayer
    IF (PRESENT(ghostsize)) THEN
        IF (ghostsize .LT. cutoff) THEN
            fail("using ghostsize < cutoff+skin. Increase ghostsize.")
        ELSE
            cutoff = ghostsize
        ENDIF
    ENDIF

#if   __KIND == __SINGLE_PRECISION
    IF (cutoff .GT. topo%ghostsizes) THEN
        stdout("cutoff for ghost mapping       = ",cutoff)
        stdout("cutoff used to create topology = ",'topo%ghostsizes')
        fail("ghostsize of topology may be smaller than that of particles")
    ENDIF
#elif   __KIND == __DOUBLE_PRECISION
    IF (cutoff .GT. topo%ghostsized) THEN
        stdout("cutoff for ghost mapping       = ",cutoff)
        stdout("cutoff used to create topology = ",'topo%ghostsized')
        fail("ghostsize of topology may be smaller than that of particles")
    ENDIF
#endif
    IF (cutoff .GT. 0._MK) THEN
        IF (Pc%flags(ppm_part_ghosts)) THEN
            IF (dbg) THEN
                stdout("ghosts have already been updated")
            ENDIF
        ENDIF
        check_true("ppm_map_type_isactive(ppm_param_map_ghost_get)",&
                "Need to call ghost_get before ghost_push")

        IF (PRESENT(Field)) THEN
            !Get discretization
            !Points the iterator (props%iter_id) to the discretization
            !of Field on the particle set Pc.
            CALL Pc%get_discr(Field,prop,info)
                or_fail("could not get discretization for this Field")
            Pc%stats%nb_ghost_push = Pc%stats%nb_ghost_push + 1
#ifdef __MPI
            t1 = MPI_WTIME(info)
#endif
            CALL Pc%map_part_push_legacy(Pc%props%iter_id,info)
                or_fail("map_part_push")
#ifdef __MPI
            t2 = MPI_WTIME(info)
            Pc%stats%t_ghost_push = Pc%stats%t_ghost_push + (t2-t1)
#endif
        ELSE
            !Update the ghost for the properties if
            ! 1) they have been mapped to this topology,
            ! 2) the ghosts have not yet been updated, and
            ! 3) the user wants them to be updated
            prop => Pc%props%begin()
            DO WHILE (ASSOCIATED(prop))

            IF (prop%flags(ppm_ppt_map_ghosts)) THEN
            IF (.NOT.prop%flags(ppm_ppt_ghosts)) THEN
                IF (prop%flags(ppm_ppt_partial)) THEN

                    IF(dbg) THEN
                        stdout("pushing property ",'Pc%props%iter_id',&
                            'TRIM(prop%name)')
                    ENDIF
                    Pc%stats%nb_ghost_push = Pc%stats%nb_ghost_push + 1
#ifdef __MPI
                    t1 = MPI_WTIME(info)
#endif
                    CALL Pc%map_part_push_legacy(Pc%props%iter_id,info)
                        or_fail("map_part_push")
#ifdef __MPI
                    t2 = MPI_WTIME(info)
                    Pc%stats%t_ghost_push = Pc%stats%t_ghost_push + (t2-t1)
#endif
                ELSE
                    stdout("pushing property ",&
                            'Pc%props%iter_id','TRIM(prop%name)')
                    fail("getting ghosts for a property thats not mapped")
                ENDIF
            ENDIF
            ENDIF
            prop => Pc%props%next()
        ENDDO
        ENDIF

    ELSE ! if cutoff .le. 0
        IF(dbg) THEN
            stdout("cutoff = 0, nothing to do")
        ENDIF
    ENDIF

    end_subroutine()
END SUBROUTINE DTYPE(part_map_ghost_push)

SUBROUTINE DTYPE(part_map_ghost_send)(Pc,info)

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
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    REAL(KIND(1.D0))                          :: t1,t2

    start_subroutine("part_map_ghost_send")

    !-----------------------------------------------------------------
    !  Checks
    !-----------------------------------------------------------------
    !check that particles are allocated
    check_associated("Pc%xp",&
        "Particles structure had not been defined. Call allocate first")
    !check that particles are mapped onto this topology
    check_true("Pc%flags(ppm_part_partial)",&
        "Do a partial/global mapping before doing a ghost mapping")

    !-----------------------------------------------------------------
    !  Send the buffer
    !-----------------------------------------------------------------
    CALL ppm_map_part_send(Pc%Npart,Pc%Mpart,info)
        or_fail("ppm_map_part_send")

    end_subroutine()
END SUBROUTINE DTYPE(part_map_ghost_send)

SUBROUTINE DTYPE(part_map_ghost_pop)(Pc,info,Field,ghostsize,debug)
    !!! Push ghost particles properties into the send buffer.
    !!! If ghost positions are not up-to-date, they are re-computed
    !!! with a call to ghost_get. Otherwise ghost_get is skipped.
    !!! If the Field argument is present, only the discretization of that
    !!! Field on this particle set is pushed. Otherwise, all the properties
    !!! (that are not already up-to-date) are pushed, except those that
    !!! have the ppm_ppt_map_ghosts flag set to false.
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
    CLASS(ppm_t_field_), OPTIONAL                       :: Field
    !!! Pop only this field. Default is to pop all of them
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
    CLASS(ppm_t_discr_data), POINTER          :: prop => NULL()

    start_subroutine("part_map_ghost_pop")

    dbg = .FALSE.
    IF (PRESENT(debug)) dbg=debug

    !-----------------------------------------------------------------
    !  Checks
    !-----------------------------------------------------------------
    !check that particles are allocated
    check_associated("Pc%xp",&
        "Particles structure had not been defined. Call allocate first")
    !check that particles are mapped onto this topology
    check_true("Pc%flags(ppm_part_partial)",&
        "Do a partial/global mapping before doing a ghost mapping")

    topoid = Pc%active_topoid
    topo=>ppm_topo(topoid)%t

    cutoff = Pc%ghostlayer
    IF (PRESENT(ghostsize)) THEN
        IF (ghostsize .LT. cutoff) THEN
            fail("using ghostsize < cutoff+skin. Increase ghostsize.")
        ELSE
            cutoff = ghostsize
        ENDIF
    ENDIF

#if   __KIND == __SINGLE_PRECISION
    IF (cutoff .GT. topo%ghostsizes) THEN
        stdout("cutoff for ghost mapping       = ",cutoff)
        stdout("cutoff used to create topology = ",'topo%ghostsizes')
        fail("ghostsize of topology may be smaller than that of particles")
    ENDIF
#elif   __KIND == __DOUBLE_PRECISION
    IF (cutoff .GT. topo%ghostsized) THEN
        stdout("cutoff for ghost mapping       = ",cutoff)
        stdout("cutoff used to create topology = ",'topo%ghostsized')
        fail("ghostsize of topology may be smaller than that of particles")
    ENDIF
#endif
    IF (cutoff .GT. 0._MK) THEN
        IF (.NOT.Pc%flags(ppm_part_ghosts).OR. &
            .NOT.ppm_map_type_isactive(ppm_param_map_ghost_get)) THEN
            fail("Ghost buffer invalid. Correct sequence is ghost_get, ghost_push, ghost_send and ghost_pop.")
        ENDIF

        IF (PRESENT(Field)) THEN
            !Get discretization
            !Points the iterator (props%iter_id) to the discretization
            !of Field on the particle set Pc.
            CALL Pc%get_discr(Field,prop,info)
                or_fail("could not get discretization for this Field")
            CALL Pc%map_part_pop_legacy(Pc%props%iter_id,&
                Pc%Mpart,info)
                or_fail("map_part_pop")
            prop%flags(ppm_ppt_ghosts) = .TRUE.
        ELSE
            prop => Pc%props%last()
            DO WHILE (ASSOCIATED(prop))

                IF (prop%flags(ppm_ppt_map_ghosts)) THEN
                IF (.NOT.prop%flags(ppm_ppt_ghosts)) THEN
                IF (prop%flags(ppm_ppt_partial)) THEN

                    IF(dbg) THEN
                        stdout("popping property ",'Pc%props%iter_id',&
                            'TRIM(prop%name)')
                    ENDIF
                    CALL Pc%map_part_pop_legacy(Pc%props%iter_id,&
                        Pc%Mpart,info)
                        or_fail("map_part_pop")
                    prop%flags(ppm_ppt_ghosts) = .TRUE.
                ENDIF
                ENDIF
                ENDIF
                prop => Pc%props%prev()
            ENDDO
        ENDIF

    ELSE ! if cutoff .le. 0

        IF(dbg) THEN
            stdout("cutoff = 0, nothing to do")
            stdout("setting all %has_ghost properties to true")
        ENDIF

        ! Update states
        prop => Pc%props%begin()
        DO WHILE (ASSOCIATED(prop))
            IF (prop%flags(ppm_ppt_map_ghosts)) THEN
                prop%flags(ppm_ppt_ghosts) = .TRUE.
            ENDIF
            prop => Pc%props%next()
        ENDDO
    ENDIF

    end_subroutine()
END SUBROUTINE DTYPE(part_map_ghost_pop)

SUBROUTINE DTYPE(part_map_ghost_pop_pos)(Pc,info,ghostsize,debug)
    !!! Pop ghost particles positions from the send buffer.
    !!!  Assumptions:
    !!! * Needs to be called at the end of the 
    !!!  ghost_get, ghost_push, send and ghost_pop call sequence.
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

    start_subroutine("part_map_ghost_pop_pos")

    dbg = .FALSE.
    IF (PRESENT(debug)) dbg=debug

    !-----------------------------------------------------------------
    !  Checks
    !-----------------------------------------------------------------
    !check that particles are allocated
    check_associated("Pc%xp",&
        "Particles structure had not been defined. Call allocate first")
    !check that particles are mapped onto this topology
    check_true("Pc%flags(ppm_part_partial)",&
        "Do a partial/global mapping before doing a ghost mapping")

    topoid = Pc%active_topoid
    topo=>ppm_topo(topoid)%t

    cutoff = Pc%ghostlayer
    IF (PRESENT(ghostsize)) THEN
        IF (ghostsize .LT. cutoff) THEN
            fail("using ghostsize < cutoff+skin. Increase ghostsize.")
        ELSE
            cutoff = ghostsize
        ENDIF
    ENDIF

#if   __KIND == __SINGLE_PRECISION
    IF (cutoff .GT. topo%ghostsizes) THEN
        stdout("cutoff for ghost mapping       = ",cutoff)
        stdout("cutoff used to create topology = ",'topo%ghostsizes')
        fail("ghostsize of topology may be smaller than that of particles")
    ENDIF
#elif   __KIND == __DOUBLE_PRECISION
    IF (cutoff .GT. topo%ghostsized) THEN
        stdout("cutoff for ghost mapping       = ",cutoff)
        stdout("cutoff used to create topology = ",'topo%ghostsized')
        fail("ghostsize of topology may be smaller than that of particles")
    ENDIF
#endif
    IF (cutoff .GT. 0._MK) THEN
        check_true("Pc%flags(ppm_part_ghosts)",&
            "flags(ppm_part_ghosts) need to be set to .true.")
        check_true("ppm_map_type_isactive(ppm_param_map_ghost_get)",&
            "Ghost buffer invalid. Correct sequence is ghost_get, ghost_push, ghost_send and ghost_pop.")

        CALL ppm_map_part_pop(Pc%xp,ppm_dim,Pc%Npart,&
            Pc%Mpart,info)
            or_fail("map_part_pop")

    ELSE ! if cutoff .le. 0

        IF(dbg) THEN
            stdout("cutoff = 0, nothing to do")
            stdout("setting all %has_ghost properties to true")
        ENDIF

    ENDIF

    ! Update states
    !   ghosts have been computed
    Pc%flags(ppm_part_ghosts) = .TRUE.

    end_subroutine()
END SUBROUTINE DTYPE(part_map_ghost_pop_pos)

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
            fail("this type of BC is not implemented/tested in this version")
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
            or_fail("Pc%del_parts: could not delete particles")
        DEALLOCATE(list_del_parts,STAT=info)
            or_fail_dealloc("list_del_parts")
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
    CALL Pc%get_xp(xp,info)
        or_fail("Particle positions cannot be accessed")
    FORALL (ip=1:Pc%Npart) &
            xp(1:ppm_dim,ip) = xp(1:ppm_dim,ip) + disp(1:ppm_dim,ip)
    CALL Pc%set_xp(xp,info)
        or_fail("set_xp")

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

SUBROUTINE DTYPE(part_neighlist)(this,info,P_xset,skin,symmetry,cutoff,name,&
        lstore,incl_ghosts,knn)
    !!!  Neighbor lists for particles
    !!!  Compute the Verlet lists for the target particles, using neighbours
    !!!  from the particle set P_set (the default is that P_xset is the
    !!!  same set as the target particles)
    !!!-----------------------------------------------------------------
    !!!  Assumptions:
    !!! * Particles positions need to have been mapped onto the topology
    !!! * Ghost positions have been computed
    USE ppm_module_neighlist
#ifdef __WITH_CNL
    USE ppm_module_cnl
#endif
    USE ppm_module_inl_vlist
    USE ppm_module_inl_xset_vlist
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
    CLASS(DTYPE(ppm_t_particles)),TARGET                   :: this
    DEFINE_MK()
    !!! Data structure containing the particles
    INTEGER,                            INTENT(  OUT)      :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    !  Optional arguments
    !-------------------------------------------------------------------------
    CLASS(DTYPE(ppm_t_particles)_),OPTIONAL,TARGET         :: P_xset
    !!! Particle set from which the neighbours are sought
    REAL(MK), OPTIONAL                                     :: skin
    !!! skin
    LOGICAL, OPTIONAL                                      :: symmetry    
    !!! if using symmetry
    REAL(MK), OPTIONAL                                     :: cutoff
    !!! cutoff radius
    CHARACTER(LEN=*) , OPTIONAL                            :: name
    !!! name of this neighbour list
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
    LOGICAL                                   :: ensure_knn,lsymm
    !!! uses a neighbour-finding algorithm that finds enough neighbours
    REAL(MK),DIMENSION(2*ppm_dim):: ghostlayer
    !!!
    REAL(KIND(1.D0))                          :: t1,t2
    TYPE(ppm_t_topo), POINTER                 :: topo => NULL()
#ifdef __WITH_KDTREE
    TYPE(DTYPE(kdtree2)),POINTER              :: tree => NULL()
    TYPE(DTYPE(kdtree2_result)),ALLOCATABLE   :: results(:)
#endif
    INTEGER                                   :: topoid
    INTEGER                                   :: nneighmin,nneighmax
    CLASS(DTYPE(ppm_t_neighlist)_), POINTER   :: Nlist => NULL()
    REAL(MK)                                  :: lskin

    REAL(MK),DIMENSION(:),POINTER             :: rcp  => NULL()
    CLASS(ppm_t_operator_discr_),POINTER      :: op => NULL()
    CLASS(DTYPE(ppm_t_particles)_), POINTER   :: Part_src => NULL()
    LOGICAL                                   :: xset_neighlists

    start_subroutine("part_comp_neighlist")

    !-----------------------------------------------------------------
    !  Checks
    !-----------------------------------------------------------------
    check_associated("this%xp",&
        "Particles structure had not been defined. Call allocate first")

    check_true("this%flags(ppm_part_partial)",&
        "Particles not mapped. Do a partial/global mapping")

    xset_neighlists = .FALSE.
    IF (PRESENT(P_xset)) THEN
        Part_src => P_xset
        check_associated("Part_src%xp",&
            "Cross-Set particles have not been defined. Call allocate first")
        check_true("Part_src%flags(ppm_part_partial)",&
            "Particles not mapped. Do a partial/global mapping")
        IF (.NOT.ASSOCIATED(Part_src,this)) THEN
            xset_neighlists = .TRUE.
        ENDIF
    ELSE    
        Part_src => this
    ENDIF


    check_true("Part_src%flags(ppm_part_ghosts)",&
            "Ghosts have not been updated. They are needed for neighlists")


    check_associated("this%neighs")

    !check whether the neighbour list already exists
    IF (this%has_neighlist(Part_src)) THEN
        Nlist => this%get_neighlist(Part_src)
        IF (PRESENT(skin).OR.PRESENT(symmetry).OR.PRESENT(cutoff)) THEN
            stdout("the optional arguments skin,",&
            "symmetry or cutoff will not be used",&
            " because the neighbour list already exists. We",&
            " should perhaps change the API?  ")
            fail("Need to destroy/re-create this neighbour list first")
        ENDIF
    ELSE
        CALL this%create_neighlist(Part_src,info,name=name,skin=skin,&
            symmetry=symmetry,cutoff=cutoff,Nlist=Nlist)
            or_fail("failed to create neighbour list")
    ENDIF

    check_associated(Nlist)

    !check that we have a cutoff radius
    SELECT TYPE (this)
    CLASS IS (DTYPE(ppm_t_sop))
        check_associated("this%rcp",&
            "cutoff radii for adaptive particles have not been defined")
        CALL this%get(this%rcp,rcp,info,with_ghosts=.TRUE.,read_only=.true.)
            or_fail("could not access cutoff radii")
    CLASS DEFAULT
        check_true("Nlist%cutoff.GT.0",&
            "cutoff is negative or zero - do we really want neighbour lists?")
    END SELECT

    IF (Nlist%isymm.EQ.1) THEN
        lsymm =.TRUE.
    ELSE
        lsymm = .FALSE.
    ENDIF

    IF (PRESENT(knn)) THEN
        ensure_knn = .TRUE.
    ELSE
        ensure_knn = .FALSE.
    ENDIF

    lskin = Nlist%skin
    topoid = this%active_topoid

    do_something: IF (Nlist%uptodate .OR. this%Npart.EQ.0) THEN
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
        np_target = this%Npart
        IF (PRESENT(incl_ghosts)) THEN
            IF (incl_ghosts) THEN
                np_target = this%Mpart
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

                this%stats%nb_kdtree = this%stats%nb_kdtree+1
#ifdef __MPI
                t1 = MPI_WTIME(info)
#endif

                tree => kdtree2_create(Part_src%xp(1:ppm_dim,1:Part_src%Mpart),&
                    sort=.true.,rearrange=.true.)
                allocate(results(knn+1),STAT=info)
                    or_fail_alloc("results")
                ldc(1) = knn
                ldc(2) = this%Npart
                CALL ppm_alloc(Nlist%vlist,ldc,ppm_param_alloc_grow,info)
                    or_fail_alloc("Nlist%vlist")
                ldc(1) = this%Npart
                CALL ppm_alloc(Nlist%nvlist,ldc,ppm_param_alloc_grow,info)
                    or_fail_alloc("Nlist%nvlist")
                DO ip=1,this%Npart
                    call kdtree2_n_nearest(tp=tree,qv=this%xp(1:ppm_dim,ip),&
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
                deallocate(results,STAT=info)
                    or_fail_dealloc("results")
#ifdef __MPI
                t2 = MPI_WTIME(info)
                this%stats%t_kdtree = this%stats%t_kdtree+(t2-t1)
#endif
#else
                fail("option required the kdtree module.")
#endif
!__WITH_KDTREE
        ELSE  

            SELECT TYPE (this)
            CLASS IS (DTYPE(ppm_t_sop))

                !FIXME: when adaptive ghost layers are available
                ghostlayer(1:2*ppm_dim)=Part_src%ghostlayer

#ifdef __WITH_CNL
                conventionalinl: IF (this%conventionalinl) THEN
                    this%stats%nb_cinl = this%stats%nb_cinl+1

#ifdef __MPI
                    t1 = MPI_WTIME(info)
#endif
                    !HUGLY HACK to make CNL routines work on a topology with
                    !several subdomains
#if   __KIND == __SINGLE_PRECISION
                    CALL cnl_vlist(this%xp,&
                        rcp,this%Npart,this%Mpart,&
                        ppm_topo(topoid)%t%min_subs(:,1)-this%ghostlayer,&
                        ppm_topo(topoid)%t%max_subs(:,1)+this%ghostlayer,&
                        this%nvlist,this%vlist,ppm_dim,info)
#elif __KIND == __DOUBLE_PRECISION
                    CALL cnl_vlist(this%xp,&
                        rcp,this%Npart,this%Mpart,&
                        ppm_topo(topoid)%t%min_subd(:,1)-this%ghostlayer,&
                        ppm_topo(topoid)%t%max_subd(:,1)+this%ghostlayer,&
                        this%nvlist,this%vlist,ppm_dim,info)
#endif
                    or_fail("ppm_cinl_vlist failed")
                    !end HUGLY HACK
#ifdef __MPI
                    t2 = MPI_WTIME(info)
                    this%stats%t_cinl = this%stats%t_cinl + (t2 - t1)
#endif
                ELSE
#endif
!__WITH_CNL
                IF (xset_neighlists) THEN

                    this%stats%nb_xset_nl = this%stats%nb_xset_nl + 1
#ifdef __MPI
                    t1 = MPI_WTIME(info)
#endif
                    CALL ppm_inl_xset_vlist(topoid,this%xp,&
                        this%Npart,this%Mpart,Part_src%xp,Part_src%Npart,&
                        Part_src%Mpart,rcp,&
                        lskin,ghostlayer,info,Nlist%vlist,&
                        Nlist%nvlist,lstore)
                        or_fail("ppm_inl_xset_vlist failed")
#ifdef __MPI
                    t2 = MPI_WTIME(info)
                    this%stats%t_xset_nl = this%stats%t_xset_nl + (t2 - t1)
#endif
                ELSE
                    this%stats%nb_inl = this%stats%nb_inl+1
#ifdef __MPI
                    t1 = MPI_WTIME(info)
#endif
                    CALL ppm_inl_vlist(topoid,this%xp,np_target,&
                        this%Mpart,rcp,lskin,lsymm,ghostlayer,info,&
                        Nlist%vlist,Nlist%nvlist)
                        or_fail("ppm_inl_vlist failed")
#ifdef __MPI
                    t2 = MPI_WTIME(info)
                    this%stats%t_inl = this%stats%t_inl + (t2 - t1)
#endif
                ENDIF ! XSET
#ifdef __WITH_CNL
                ENDIF conventionalinl
#endif

            CLASS DEFAULT
                IF (xset_neighlists) THEN

                    this%stats%nb_xset_nl = this%stats%nb_xset_nl + 1
#ifdef __MPI
                    t1 = MPI_WTIME(info)
#endif
                    CALL ppm_inl_xset_vlist(topoid,this%xp,&
                        this%Npart,this%Mpart,Part_src%xp,Part_src%Npart,&
                        Part_src%Mpart,Nlist%cutoff,&
                        lskin,ghostlayer,info,Nlist%vlist,&
                        Nlist%nvlist,lstore)
                        or_fail("ppm_inl_xset_vlist failed")
#ifdef __MPI
                    t2 = MPI_WTIME(info)
                    this%stats%t_xset_nl = this%stats%t_xset_nl + (t2 - t1)
#endif
                ELSE
                    this%stats%nb_nl = this%stats%nb_nl+1

#ifdef __MPI
                    t1 = MPI_WTIME(info)
#endif
                    CALL ppm_neighlist_vlist(topoid,this%xp,this%Mpart,&
                        Nlist%cutoff,lskin,lsymm,Nlist%vlist,&
                        Nlist%nvlist,info,lstore=lstore)
                    or_fail("ppm_neighlist_vlist failed")
#ifdef __MPI
                    t2 = MPI_WTIME(info)
                    this%stats%t_nl = this%stats%t_nl + (t2 - t1)
#endif
                ENDIF ! XSET
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

        Nlist%nneighmin = MINVAL(Nlist%nvlist(1:this%Npart))
        Nlist%nneighmax = MAXVAL(Nlist%nvlist(1:np_target))

        ! DC operators that do not use a xset neighbour list, if they exist, 
        ! are no longer valid (they depend on the neighbour lists)
        IF (ASSOCIATED(this%ops)) THEN
            op => this%ops%begin()
            DO WHILE (ASSOCIATED(op))
                IF (.NOT.op%flags(ppm_ops_interp)) THEN
                    op%flags(ppm_ops_iscomputed) = .FALSE.
                ENDIF
                op => this%ops%next()
            ENDDO
        ENDIF

        ! We want to distinguish between "self" neighbour lists
        ! and cross-set ones.
        IF (ASSOCIATED(Nlist%Part,this)) THEN
            this%flags(ppm_part_neighlists) = .TRUE.
        ENDIF
        
        Nlist => NULL()


    ENDIF do_something

    end_subroutine()
END SUBROUTINE DTYPE(part_neighlist)


SUBROUTINE DTYPE(part_set_cutoff)(Pc,cutoff,info,Nlist)
    !!! Set a cutoff radius for a Particle set and update the
    !!! the ghostlayer sizes.
    !!! The cutoff radius concerns the default neighbor list, unless
    !!! specified otherwise.
    !!! If the cutoff is increased from its previous value, the neighbour
    !!! list is flagged as "not up-to-date" and will have to be recomputed
    !!! before it is used again.
    !!! If the ghostlayer sizes are increased from their previous values,
    !!! the ghosts are flagged as "not up-to-date" and will have to be
    !!! recomputed before they are used again.
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
    !!! applies. By default, this is the "standard" Verlet list, with neighbours
    !!! sought within the particle set itself.
    
    !-------------------------------------------------------------------------
    ! local variables
    !-------------------------------------------------------------------------
    CLASS(DTYPE(ppm_t_neighlist)_),   POINTER :: nl => NULL()

    start_subroutine("part_set_cutoff")

    !-------------------------------------------------------------------------
    !  Set new cutoff
    !-------------------------------------------------------------------------
    IF (PRESENT(NList)) THEN
        check_true("Pc%neighs%has(NList)",&
            "Neighbour list does not concern this particle set")
        IF (cutoff .LT. NList%cutoff) NList%uptodate = .FALSE.
        NList%cutoff = cutoff
    ELSE
        nl => Pc%get_neighlist()
        check_associated(nl,"Compute neighbour lists first")
        IF (cutoff .LT. nl%cutoff) nl%uptodate = .FALSE.
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
        CALL Pc%create_prop(info,part_prop=Pc%gi,dtype=ppm_type_int,&
            name="GlobalIndex")
        Pc%flags(ppm_part_global_index) = .TRUE.
    END IF
#ifdef __MPI
    CALL MPI_Scan(Pc%Npart,offset,1,MPI_INTEGER,MPI_SUM,ppm_comm,info)
    offset = offset - Pc%Npart
#else
    offset = 0
#endif
    CALL Pc%get(Pc%gi,wp,info)
    FORALL (i=1:Pc%Npart) wp(i) = offset + i !- 1 !uncomment if index from 0
    CALL Pc%set(Pc%gi,wp,info)

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
                    or_fail_alloc("cont%vec")
                id = 1
            ELSE
                !need to resize the array of mapping pointers 
                vec_size=MAX(2*cont%vec_size,20)
                ALLOCATE(vec_tmp(vec_size),STAT=info)
                    or_fail_alloc("vec_tmp")
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
            fail("mapping id larger than size of mappings array")
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
                        fail("fatal error in the data structure")
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
    CLASS(ppm_t_discr_data), POINTER               :: prop => NULL()

    start_function("has_ghosts")

    IF (PRESENT(Field)) THEN
        CALL this%get_discr(Field,prop,info)
            or_fail("this field is not discretized on this particle set")
        res = prop%flags(ppm_ppt_ghosts)
    ELSE
        res = this%flags(ppm_part_ghosts)
    ENDIF

    end_function()

END FUNCTION

SUBROUTINE DTYPE(part_get_discr)(this,Field,prop,info)
    !!! Returns a pointer to the ppm_t_discr_data object that is
    !!! the discretization of that Field on this Particle set.
    !!! Fails with an error if the Field is not discretized here.
    CLASS(DTYPE(ppm_t_particles))                          :: this
    CLASS(ppm_t_field_),TARGET,             INTENT(IN   )  :: Field
    CLASS(ppm_t_discr_data),POINTER,        INTENT(  OUT)  :: prop
    INTEGER,                                INTENT(  OUT)  :: info

    start_subroutine("part_get_discr")
    
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


SUBROUTINE DTYPE(part_prop_zero)(this,Field,info)
    !!! Reset values of a property to zero
    !!! (a bit unrolled)
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_particles))                          :: this
    CLASS(ppm_t_field_),TARGET,             INTENT(IN   )  :: Field
    INTEGER,                                INTENT(  OUT)  :: info

    INTEGER                        :: lda

    start_subroutine("part_prop_zero")
    
    lda = Field%lda
    SELECT CASE(lda)
    CASE (1)
        foreach p in particles(this) with sca_fields(w=Field) prec(DTYPE(prec))
            w_p = 0._mk
        end foreach
    CASE (2)
        foreach p in particles(this) with vec_fields(w=Field) prec(DTYPE(prec))
            w_p(1) = 0._mk
            w_p(2) = 0._mk
        end foreach
    CASE (3)
        foreach p in particles(this) with vec_fields(w=Field) prec(DTYPE(prec))
            w_p(1) = 0._mk
            w_p(2) = 0._mk
            w_p(3) = 0._mk
        end foreach
    CASE (4)
        foreach p in particles(this) with vec_fields(w=Field) prec(DTYPE(prec))
            w_p(1) = 0._mk
            w_p(2) = 0._mk
            w_p(3) = 0._mk
            w_p(4) = 0._mk
        end foreach
    CASE (5)
        foreach p in particles(this) with vec_fields(w=Field) prec(DTYPE(prec))
            w_p(1) = 0._mk
            w_p(2) = 0._mk
            w_p(3) = 0._mk
            w_p(4) = 0._mk
            w_p(5) = 0._mk
        end foreach
    CASE DEFAULT
        foreach p in particles(this) with vec_fields(w=Field) prec(DTYPE(prec))
            w_p(1:lda) = 0._mk
        end foreach
    END SELECT


    end_subroutine()
END SUBROUTINE

