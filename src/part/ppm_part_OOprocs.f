SUBROUTINE DTYPE(get_xp)(Particles,xp,with_ghosts)
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_particles))         :: Particles
    LOGICAL,OPTIONAL                      :: with_ghosts
    REAL(MK),DIMENSION(:,:),     POINTER  :: xp
    INTEGER                               :: info

    IF (PRESENT(with_ghosts)) THEN
        IF (with_ghosts) THEN
            IF (Particles%flags(ppm_part_ghosts)) THEN
                xp => Particles%xp(1:ppm_dim,1:Particles%Mpart)
            ELSE
                write(cbuf,*) 'WARNING: tried to get xp with ghosts ',&
                    'when ghosts are not up-to-date'
                CALL ppm_write(ppm_rank,'get_xp',cbuf,info)
                xp => NULL()
            ENDIF
            RETURN
        ENDIF
    ENDIF

    xp => Particles%xp(1:ppm_dim,1:Particles%Npart)

END SUBROUTINE DTYPE(get_xp)

SUBROUTINE DTYPE(set_xp)(Particles,xp,read_only,ghosts_ok)
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_particles))    :: Particles
    LOGICAL,OPTIONAL                 :: read_only
    LOGICAL,OPTIONAL                 :: ghosts_ok
    REAL(MK),DIMENSION(:,:),POINTER  :: xp
    INTEGER                          :: i

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

    Particles%flags(ppm_part_areinside) = .FALSE.
    Particles%flags(ppm_part_partial) = .FALSE.
    Particles%flags(ppm_part_ghosts) = .FALSE.
    Particles%flags(ppm_part_cartesian) = .FALSE.
    DO i = 1,Particles%max_wpid
        Particles%props(i)%t%flags(ppm_ppt_ghosts) = .FALSE.
        Particles%props(i)%t%flags(ppm_ppt_partial) = .FALSE.
    ENDDO
    DO i = 1,Particles%max_nlid
        Particles%neighs(i)%t%uptodate = .FALSE.
    ENDDO

    xp => NULL()

END SUBROUTINE DTYPE(set_xp)

SUBROUTINE DTYPE(part_prop_create)(Particles,propid,datatype,info,&
        lda,name,zero,with_ghosts)
    !!! Adds a property to an existing particle set
    CLASS(DTYPE(ppm_t_particles))         :: Particles
    INTEGER,                INTENT(INOUT) :: propid
    INTEGER,                INTENT(IN)    :: datatype
    INTEGER, OPTIONAL,      INTENT(IN)    :: lda
    CHARACTER(LEN=*) , OPTIONAL           :: name
    !!! name to this property
    LOGICAL, OPTIONAL                     :: zero
    !!! if true, then initialise the data to zero
    LOGICAL, OPTIONAL                     :: with_ghosts
    !!! if true, then allocate with Mpart instead of the default size of Npart
    INTEGER,               INTENT(OUT)    :: info

    INTEGER                               :: lda2,nprops,npart,i
    CHARACTER(LEN=ppm_char)               :: caller = 'particle_prop_create'
    CHARACTER(LEN=ppm_char)               :: name2
    REAL(KIND(1.D0))                      :: t0
    TYPE(DTYPE(ppm_ptr_part_prop)),DIMENSION(:),POINTER  :: prop_tmp => NULL()
    LOGICAL, DIMENSION(ppm_param_length_pptflags):: flags

    CALL substart(caller,t0,info)

    !Generate a new propid
    IF (propid .EQ. 0) THEN
        IF (Particles%nwp.LT.Particles%swp) THEN
            !there is at least one empty slot in the array
            ! of property pointers
            propid = propid + 1
            DO WHILE (Particles%props(propid)%t%data_type .EQ. ppm_type_none)
                propid = propid + 1

!--------not a necessary check-------
                IF (propid .GT. Particles%swp) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_alloc,caller,&
         &    'number of properties greater than allocated size',&
                        __LINE__,info)
                    GOTO 9999
                ENDIF
!------------------------------------

            ENDDO
        ELSE
            IF (.NOT. ASSOCIATED(Particles%props)) THEN
            !need to allocate the array of property pointers 
                nprops=20
                ALLOCATE(Particles%props(nprops),STAT=info)
                IF (info .NE. 0) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_alloc,caller,&
                        'allocating property array failed',__LINE__,info)
                    GOTO 9999
                ENDIF
                propid = 1
            ELSE
            !need to resize the array of property pointers 
                nprops=MAX(2*Particles%swp,20)
                ALLOCATE(prop_tmp(nprops),STAT=info)
                IF (info .NE. 0) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_alloc,caller,&
                        'allocating property array failed',__LINE__,info)
                    GOTO 9999
                ENDIF
                DO i=1,Particles%swp
                    prop_tmp(i)%t => Particles%props(i)%t
                ENDDO
                DEALLOCATE(Particles%props)
                Particles%props => prop_tmp
            ENDIF
            Particles%swp = nprops
            propid = Particles%nwp+1
        ENDIF
        Particles%nwp = Particles%nwp + 1
    ELSE
        !overwriting an existing property
    ENDIF
        

    IF (propid .GT. Particles%max_wpid) Particles%max_wpid = propid

    !
    lda2 = 1
    IF (PRESENT(lda)) THEN
        IF (lda.GE.2) THEN
            lda2 = lda
        ENDIF
    ENDIF

    IF (PRESENT(name)) THEN
        name2 = name
    ELSE
        name2 = particles_dflt_pptname(propid,1)
    ENDIF

    npart = Particles%Npart
    flags = .FALSE.
    IF (PRESENT(with_ghosts)) THEN
        IF (with_ghosts) THEN
            npart = Particles%Mpart
            flags(ppm_ppt_ghosts) = .TRUE.
        ENDIF
    ENDIF
    flags(ppm_ppt_partial) = .TRUE.

    IF (.NOT. ASSOCIATED(Particles%props(propid)%t)) THEN
        ALLOCATE(Particles%props(propid)%t,STAT=info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,&
                'allocating property pointer failed',__LINE__,info)
            GOTO 9999
        ENDIF
    ENDIF
    CALL Particles%props(propid)%t%create(datatype,npart,lda2,name2,flags,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_sub_failed,caller,&
            'creating property array failed',__LINE__,info)
        GOTO 9999
    ENDIF

    CALL substop(caller,t0,info)
    9999  CONTINUE

END SUBROUTINE DTYPE(part_prop_create)

SUBROUTINE DTYPE(prop_create)(prop,datatype,npart,lda,name,flags,info)
    !!! Constructor for particle property data structure
    CLASS(DTYPE(part_prop))            :: prop
    INTEGER,                INTENT(IN) :: datatype
    INTEGER,                INTENT(IN) :: npart
    INTEGER,                INTENT(IN) :: lda
    CHARACTER(LEN=*)                   :: name
    !!! name to this property
    LOGICAL, DIMENSION(ppm_param_length_pptflags),OPTIONAL,INTENT(IN) :: flags
    INTEGER,               INTENT(OUT) :: info

    REAL(KIND(1.D0))                   :: t0
    INTEGER                            :: iopt
    CHARACTER(LEN=ppm_char)            :: caller = 'prop_create'
    LOGICAL                            :: is2d



    CALL substart(caller,t0,info)

    prop%lda       = lda
    prop%data_type = datatype
    prop%name      = name
    prop%flags     = flags

    ldc(1) = npart
    iopt   = ppm_param_alloc_fit

    IF (lda.GE.2) THEN
        ldc(2) = lda
        is2d = .TRUE.
    ELSE
        is2d = .FALSE.
    ENDIF

    IF (is2d) THEN
        SELECT CASE (datatype)
        CASE (ppm_type_int)
            CALL ppm_alloc(prop%data_2d_i,ldc,iopt,info)
        CASE (ppm_type_longint)
            CALL ppm_alloc(prop%data_2d_li,ldc,iopt,info)
#if   __KIND == __SINGLE_PRECISION
        CASE (ppm_type_real_single )
            CALL ppm_alloc(prop%data_2d_r,ldc,iopt,info)
        CASE (ppm_type_comp_single)
            CALL ppm_alloc(prop%data_2d_c,ldc,iopt,info)
#elif __KIND ==__DOUBLE_PRECISION
        CASE (ppm_type_real_double)
            CALL ppm_alloc(prop%data_2d_r,ldc,iopt,info)
        CASE (ppm_type_comp_double)
            CALL ppm_alloc(prop%data_2d_c,ldc,iopt,info)
#endif
        CASE (ppm_type_logical )
            CALL ppm_alloc(prop%data_2d_l,ldc,iopt,info)
        CASE DEFAULT
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_argument,caller,   &
                &        'invalid type for particle property',__LINE__,info)
        END SELECT
    ELSE
        SELECT CASE (datatype)
        CASE (ppm_type_int)
            CALL ppm_alloc(prop%data_1d_i,ldc,iopt,info)
        CASE (ppm_type_longint)
            CALL ppm_alloc(prop%data_1d_li,ldc,iopt,info)
#if   __KIND == __SINGLE_PRECISION
        CASE (ppm_type_real_single )
            CALL ppm_alloc(prop%data_1d_r,ldc,iopt,info)
        CASE (ppm_type_comp_single)
            CALL ppm_alloc(prop%data_1d_c,ldc,iopt,info)
#elif __KIND ==__DOUBLE_PRECISION
        CASE (ppm_type_real_double)
            CALL ppm_alloc(prop%data_1d_r,ldc,iopt,info)
        CASE (ppm_type_comp_double)
            CALL ppm_alloc(prop%data_1d_c,ldc,iopt,info)
#endif
        CASE (ppm_type_logical )
            CALL ppm_alloc(prop%data_1d_l,ldc,iopt,info)
        CASE DEFAULT
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_argument,caller,   &
                &        'invalid type for particle property',__LINE__,info)
        END SELECT
    ENDIF

    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,&
            'allocating property failed',__LINE__,info)
        GOTO 9999
    ENDIF

    CALL substop(caller,t0,info)

    9999  CONTINUE

END SUBROUTINE DTYPE(prop_create)

SUBROUTINE DTYPE(prop_destroy)(prop)
    CLASS(DTYPE(part_prop))      :: prop

    IF(ASSOCIATED(prop%data_1d_i)) DEALLOCATE(prop%data_1d_i)
    IF(ASSOCIATED(prop%data_2d_i)) DEALLOCATE(prop%data_2d_i)
    IF(ASSOCIATED(prop%data_1d_li)) DEALLOCATE(prop%data_1d_li)
    IF(ASSOCIATED(prop%data_2d_li)) DEALLOCATE(prop%data_2d_li)
    IF(ASSOCIATED(prop%data_1d_r)) DEALLOCATE(prop%data_1d_r)
    IF(ASSOCIATED(prop%data_2d_r)) DEALLOCATE(prop%data_2d_r)
    IF(ASSOCIATED(prop%data_1d_c)) DEALLOCATE(prop%data_1d_c)
    IF(ASSOCIATED(prop%data_2d_c)) DEALLOCATE(prop%data_2d_c)
    IF(ASSOCIATED(prop%data_1d_l)) DEALLOCATE(prop%data_1d_l)
    IF(ASSOCIATED(prop%data_2d_l)) DEALLOCATE(prop%data_2d_l)

    prop%data_type = ppm_type_none
    prop%lda = 0

END SUBROUTINE DTYPE(prop_destroy)

SUBROUTINE DTYPE(neigh_destroy)(neigh)
    CLASS(DTYPE(ppm_t_neighlist))      :: neigh

    IF(ASSOCIATED(neigh%nvlist)) DEALLOCATE(neigh%nvlist)
    IF(ASSOCIATED(neigh%vlist))  DEALLOCATE(neigh%vlist)

END SUBROUTINE DTYPE(neigh_destroy)

SUBROUTINE DTYPE(op_destroy)(op)
    CLASS(DTYPE(ppm_t_operator))              :: op
    INTEGER                                   :: i
    
    DO i=1,op%max_opsid
        CALL op%ker(i)%t%destroy()
        CALL op%desc(i)%t%destroy()
    ENDDO
    op%max_opsid = 0
    op%nb_ops = 0

END SUBROUTINE DTYPE(op_destroy)

SUBROUTINE DTYPE(desc_destroy)(desc)
    CLASS(DTYPE(ppm_t_opdesc))              :: desc
    
    IF (ASSOCIATED(desc%degree)) DEALLOCATE(desc%degree)
    IF (ASSOCIATED(desc%order))  DEALLOCATE(desc%order)
    IF (ASSOCIATED(desc%coeffs)) DEALLOCATE(desc%coeffs)

END SUBROUTINE DTYPE(desc_destroy)

SUBROUTINE DTYPE(part_create)(Particles,Npart,info,name)

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
    CLASS(DTYPE(ppm_t_particles))                          :: Particles
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
    REAL(KIND(1.D0))                                :: t0
    CHARACTER(LEN = ppm_char)            :: caller = 'ppm_create_particles'
    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)

    NULLIFY(Particles%xp)
    NULLIFY(Particles%neighs)
    NULLIFY(Particles%props)
    NULLIFY(Particles%ops)
    !-----------------------------------------------------------------
    !  Allocate memory for the positions
    !-----------------------------------------------------------------
    ldc(1) = ppm_dim
    ldc(2) = Npart
    CALL ppm_alloc(Particles%xp,ldc(1:2),ppm_param_alloc_fit,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &        'Could not allocate Particles elements',__LINE__,info)
        GOTO 9999
    ENDIF
    Particles%Npart = Npart
    Particles%Mpart = Npart
    Particles%flags(ppm_part_ghosts) = .FALSE.
    Particles%flags(ppm_part_areinside) = .FALSE.
    Particles%flags(ppm_part_partial) = .FALSE.
    Particles%flags(ppm_part_reqput) = .FALSE.
    Particles%active_topoid = -1
    ! No active topology yet

    ! Give a default name to this Particle set
    IF (PRESENT(name)) THEN
        Particles%name = ADJUSTL(TRIM(name))
    ELSE
        Particles%name = particles_dflt_partname()
    ENDIF

    ! No properties defined
    Particles%nwp = 0
    Particles%swp = 0
    Particles%max_wpid = 0
    Particles%props => NULL()

    ! No neighbor lists defined
    Particles%nnl = 0
    Particles%snl = 0
    Particles%max_nlid = 0
    Particles%neighs => NULL()

    ! Particles are by default not adaptive
    Particles%adaptive = .FALSE.
    Particles%adapt_wpid = 0
    Particles%gi_id = 0
    Particles%rcp_id = 0
    Particles%D_id = 0
    Particles%Dtilde_id = 0
    Particles%nn_sq_id = 0
    ! Particles do not represent a level-set function
    Particles%level_set = .FALSE.
    Particles%level_id = 0
    !        Particles%level_old_id = 0
    Particles%level_grad_id = 0
    !        Particles%level_grad_old_id = 0
    ! Particles are by default isotropic
    Particles%anisotropic = .FALSE.
    ! Particles are by default not a Cartesian grid
    Particles%cartesian = .FALSE.
    ! Particles have not been initialised yet
    Particles%h_avg = -1._MK
    Particles%h_min = -1._MK

    Particles%time = 0._MK
    Particles%itime = 0


    !-------------------------------------------------------------------------
    !  Finalize
    !-------------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999 CONTINUE

END SUBROUTINE DTYPE(part_create)

SUBROUTINE DTYPE(part_destroy)(Particles,info)

    !!! Deallocate a ppm_t_particles data type

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    DEFINE_MK()
    CLASS(DTYPE(ppm_t_particles))                          :: Particles
    !!! Data structure containing the particles
    INTEGER,                                INTENT(  OUT)  :: info
    !!! Returns status, 0 upon success.

    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    LOGICAL                                         :: lalloc,ldealloc
    INTEGER                                         :: i
    REAL(KIND(1.D0))                                :: t0
    CHARACTER(LEN = ppm_char)            :: caller = 'ppm_destroy_particles'
    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)

    !-------------------------------------------------------------------------
    !  Check arguments
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    !  If reallocating, deallocate old data first
    !-------------------------------------------------------------------------
    !----------------------------------------------------------------------
    !  deallocate
    !----------------------------------------------------------------------
    ! first deallocate all content of Particles
    IF (ASSOCIATED(Particles%xp)) DEALLOCATE(Particles%xp,STAT=info)
    IF (ASSOCIATED(Particles%pcost)) DEALLOCATE(Particles%pcost,STAT=info)

    !Deallocate neighbour lists
    IF (ASSOCIATED(Particles%neighs)) THEN
        DO i=1,Particles%max_nlid
            CALL Particles%neighs(i)%t%destroy()
        ENDDO
    ENDIF

    !Deallocate properties
    IF (ASSOCIATED(Particles%props)) THEN
        DO i=1,Particles%max_wpid
            CALL Particles%props(i)%t%destroy()
        ENDDO
    ENDIF

    IF (ASSOCIATED(Particles%ops)) THEN
        CALL Particles%ops%destroy()
!        CALL particles_dcop_deallocate(Particles,info)
!        IF (info .NE. 0) THEN
!            info = ppm_error_error
!            CALL ppm_error(ppm_err_dealloc,caller,   &
!                &          'Deallocating ops',__LINE__,info)
!            GOTO 9999
!        ENDIF
    ENDIF

    !-------------------------------------------------------------------------
    !  Finalize
    !-------------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999 CONTINUE

END SUBROUTINE DTYPE(part_destroy)

#undef DEFINE_MK


