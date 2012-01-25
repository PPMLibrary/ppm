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
    INTEGER                          :: info

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

    info = 0
    CALL particles_updated_positions(Particles,info)
    xp => NULL()

END SUBROUTINE DTYPE(set_xp)

SUBROUTINE DTYPE(ppm_alloc_particles)(Particles,Npart,iopt,info,name)

    !!! (Re)allocates the memory of ppm_t_particles data type

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
    TYPE(DTYPE(ppm_t_particles)),POINTER,   INTENT(INOUT)  :: Particles
    !!! Data structure containing the particles
    INTEGER,                                INTENT(IN   )  :: Npart
    !!! Number of particles 
    INTEGER                   ,             INTENT(IN   )  :: iopt
    !!! Allocation mode. One of:
    !!!
    !!! * ppm_param_alloc_fit
    !!! * ppm_param_dealloc
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
    CHARACTER(LEN = ppm_char)            :: caller = 'ppm_alloc_particles'
    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)

    !-------------------------------------------------------------------------
    !  Check arguments
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    !  Check the allocation type
    !-------------------------------------------------------------------------
    lalloc   = .FALSE.
    ldealloc = .FALSE.
    IF (iopt.EQ.ppm_param_alloc_fit) THEN
        !----------------------------------------------------------------------
        !  fit memory but skip the present contents
        !----------------------------------------------------------------------
        IF (ASSOCIATED(Particles)) ldealloc = .TRUE.
        lalloc   = .TRUE.
    ELSEIF (iopt.EQ.ppm_param_dealloc) THEN
        ldealloc = .TRUE.
    ELSE
        !----------------------------------------------------------------------
        !  Unknown iopt
        !----------------------------------------------------------------------
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,                 &
            &                  'unknown iopt',__LINE__,info)
        GOTO 9999
    ENDIF

    !-------------------------------------------------------------------------
    !  If reallocating, deallocate old data first
    !-------------------------------------------------------------------------
    IF (ldealloc) THEN
        !----------------------------------------------------------------------
        !  deallocate
        !----------------------------------------------------------------------
        IF (ASSOCIATED(Particles)) THEN
            ! first deallocate all content of Particles
            IF (ASSOCIATED(Particles%xp)) DEALLOCATE(Particles%xp,STAT=info)
            IF (ASSOCIATED(Particles%pcost)) DEALLOCATE(Particles%pcost,STAT=info)
            IF (ASSOCIATED(Particles%nvlist)) DEALLOCATE(Particles%nvlist,STAT=info)

            IF (ASSOCIATED(Particles%props)) THEN
                DO i=1,Particles%max_wpid
                    IF (ASSOCIATED(Particles%props(i))) &
                        CALL Particles%props(i)%destroy()
                    NULLIFY(Particles%props(i))
                ENDDO
                DEALLOCATE(Particles%props,STAT=info)
            ENDIF
            IF (ASSOCIATED(Particles%ops)) THEN
                CALL particles_dcop_deallocate(Particles,info)
                IF (info .NE. 0) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_dealloc,caller,   &
                        &          'Deallocating ops',__LINE__,info)
                    GOTO 9999
                ENDIF
            ENDIF
            DEALLOCATE(Particles,stat=info)
            IF (info .NE. 0) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_dealloc,caller,   &
                    &          'Deallocating Particles',__LINE__,info)
                GOTO 9999
            ENDIF
            NULLIFY(Particles)
        ENDIF
    ENDIF

    !-------------------------------------------------------------------------
    !  Allocate new memory
    !-------------------------------------------------------------------------
    IF (lalloc) THEN

        ALLOCATE(Particles,STAT=info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,   &
                &           'Allocating Particles',__LINE__,info)
            GOTO 9999
        ENDIF
        NULLIFY(Particles%xp)
        NULLIFY(Particles%nvlist)
        NULLIFY(Particles%vlist)
        NULLIFY(Particles%wpi)
        NULLIFY(Particles%wps)
        NULLIFY(Particles%wpv)
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
        Particles%has_ghosts = .FALSE.

        ! Give a default name to this Particle set
        IF (PRESENT(name)) THEN
            Particles%name = ADJUSTL(TRIM(name))
        ELSE
            Particles%name = particles_dflt_partname()
        ENDIF
        ! No properties defined yet
        Particles%nwpi = 0
        Particles%nwps = 0
        Particles%nwpv = 0
        Particles%max_wpiid = 0
        Particles%max_wpsid = 0
        Particles%max_wpvid = 0
        ! No active topology yet
        Particles%active_topoid = -1
        ! Particles are not yet mapped on any topology
        Particles%ontopology = .FALSE.
        ! Neighbour lists are not yet computed
        Particles%neighlists = .FALSE.
        Particles%conventionalinl = .FALSE.
        Particles%nneighmin = 0
        Particles%neighlists_cross = .FALSE.
        Particles%nneighmin_cross = 0
        Particles%nneighmax_cross = 0
        ! Default cutoff is zero
        Particles%cutoff = 0._MK
        ! Default skin is zero
        Particles%skin = 0._MK
        ! Particles interactions are by default not assumed symmetric
        Particles%isymm = 0
        ! Particles are by default not adaptive
        Particles%adaptive = .FALSE.
        Particles%adapt_wpid = 0
!        Particles%adapt_wpgradid = 0
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
        Particles%areinside = .FALSE.

        Particles%time = 0._MK
        Particles%itime = 0

        IF (verbose) &
            write(*,*) 'Allocated particles with Np=',Npart


    ENDIF

    !-------------------------------------------------------------------------
    !  Finalize
    !-------------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999 CONTINUE

END SUBROUTINE DTYPE(ppm_alloc_particles)

#undef DEFINE_MK


