MODULE ppm_module_particles

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

#define __KIND __DOUBLE_PRECISION

USE ppm_module_typedef
USE ppm_module_alloc
USE ppm_module_error
USE ppm_module_data, ONLY: ppm_dim

IMPLICIT NONE

#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER,PRIVATE :: prec = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER,PRIVATE :: prec = ppm_kind_double
#endif

    INTEGER, PARAMETER :: ppm_param_part_init_cartesian = 1
    INTEGER, PARAMETER :: ppm_param_part_init_random = 2

!for debugging
LOGICAL           :: verbose = .FALSE.

TYPE pnt_array_1d
    REAL(ppm_kind_double), DIMENSION(:), POINTER :: vec => NULL()
END TYPE pnt_array_1d

TYPE pnt_array_2d
    REAL(ppm_kind_double), DIMENSION(:,:), POINTER :: vec =>NULL()
END TYPE pnt_array_2d

TYPE ppm_t_particles

    REAL(prec), DIMENSION(:,:), POINTER             :: xp => NULL()
    !!! positions of the particles
    INTEGER                                         :: xp_g 
    !!! pseudo-boolean for the positions
    !!! takes value 1 if the ghost values have been computed
    !!! 0 if they have not, and -1 if they do not need to be updated
    INTEGER                                         :: xp_m 
    !!! pseudo-boolean for the positions
    !!! takes value 1 if the real particles are mapped one-to-one
    !!! onto the xp array, 0 if they are not, and -1 if they do 
    !!! not need to be updated
    INTEGER                                         :: Npart
    !!! Number of real particles on this processor
    INTEGER                                         :: Mpart
    !!! Number of particles (including ghosts) on this processor
    LOGICAL                                         :: areinside
    !!! true if all the particles are inside the comp. domain


    ! Particles properties
    !   scalar
    TYPE(pnt_array_1d),    DIMENSION(:),   POINTER  :: wps  => NULL()
    !!! array of scalar properties
    INTEGER              , DIMENSION(:  ), POINTER  :: wps_g=> NULL()
    !!! array of pseudo-booleans for the scalar properties
    !!! takes values 1 if the ghost values have been computed
    !!! 0 if they have not, and -1 if they do not need to be updated
    INTEGER              , DIMENSION(:  ), POINTER  :: wps_m=> NULL()
    !!! array of pseudo-booleans for the scalar properties
    !!! takes values 1 if the values are used on the active topology
    !!! 0 if they are not, and -1 if they do not need to be updated
    !!! Note that the value of 0 indicates that the array is allocated 
    !!! to the correct whereas a value of -1 does not.
    INTEGER                                         :: nwps
    !!! number of scalar properties
    INTEGER                                         :: max_wpsid
    !!! maximum index for scalar properties


    !   vectors
    TYPE(pnt_array_2d),    DIMENSION(:),   POINTER  :: wpv  => NULL()
    !!! array of vector properties
    INTEGER              , DIMENSION(:  ), POINTER  :: wpv_s=> NULL()
    !!! leading dimension of each vector property
    INTEGER              , DIMENSION(:  ), POINTER  :: wpv_g=> NULL()
    !!! array of pseudo-booleans for the vector properties
    !!! takes values 1 if the ghost values have been computed
    !!! 0 if they have not, and -1 if they do not need to be updated
    INTEGER              , DIMENSION(:  ), POINTER  :: wpv_m=> NULL()
    !!! array of pseudo-booleans for the vector properties
    !!! takes values 1 if the values are used on the active topology
    !!! 0 if they are not, and -1 if they do not need to be updated
    !!! Note that the value of 0 indicates that the array is allocated 
    !!! to the correct whereas a value of -1 does not.
    INTEGER                                         :: nwpv
    !!! number of vector properties
    INTEGER                                         :: max_wpvid
    !!! maximum index for vector properties


    ! Load balancing
    REAL(prec), DIMENSION(:  ), POINTER             :: pcost=> NULL()
    !!! cost associated to each particle 


    ! Neighbor lists
    REAL(prec)                                      :: cutoff
    !!! cutoff radius
    REAL(prec)                                      :: skin
    !!! skin layer around the particles
    INTEGER                                         :: isymm
    !!! using symmetry
    INTEGER              , DIMENSION(:  ), POINTER  :: nvlist=> NULL()
    !!! Number of neighbors of each particles
    INTEGER              , DIMENSION(:,:), POINTER  :: vlist=> NULL()
    !!! Neighbor lists
    LOGICAL                                         :: neighlists
    !!! true if the neighbor lists have been computed
    INTEGER                                         :: nneighmin
    !!! smallest number of neighbors
    INTEGER                                         :: nneighmax
    !!! highest number of neighbors

    
    ! xset Neighbor lists
    INTEGER              , DIMENSION(:  ), POINTER  :: nvlist_cross=> NULL()
    !!! Number of neighbors of each particles
    INTEGER              , DIMENSION(:,:), POINTER  :: vlist_cross=> NULL()
    !!! Neighbor lists
    LOGICAL                                         :: neighlists_cross
    !!! true if the neighbor lists have been computed
    INTEGER                                         :: nneighmin_cross
    !!! smallest number of neighbors
    INTEGER                                         :: nneighmax_cross
    !!! highest number of neighbors


    ! Links between Particles and topologies
    INTEGER                                         :: active_topoid
    !!! active topology for these particles
    LOGICAL                                         :: ontopology
    !!! true if the particles have been mapped on the active topology


    ! DC-PSE
    INTEGER                                         :: eta_id
    !!! index of the wpv array where the current DC operators are stored

    ! Adaptive particles
    LOGICAL                                         :: adaptive
    !!! true if the particles have their own cutoff radii
    !!! in this case, the cutoff will be stored in wps(rcp_id)%vec
    INTEGER                                         :: rcp_id
    !!! index of the wps array where the cutoff radius is stored
    INTEGER                                         :: D_id
    !!! index of the wps array where D is stored
    INTEGER                                         :: Dtilde_id
    !!! index of the wps array where D_tilde is stored
    INTEGER                                         :: adapt_wpid
    !!! index of the wps array where is stored the property on 
    !!! which adaptation is based 
    !!! default is first 1d property that is not rcp_id (if any)
    !!! otherwise, it is rcp_id
    INTEGER                                         :: adapt_wpgradid
    !!! index of the wpv array where is stored the gradient of the property 
    !!! on which adaptation is based (if needed)
    LOGICAL                                         :: level_set
    !!! true if particles carry a level-set function
    INTEGER                                         :: level_id
    !!! index of the wps array where the level-set is stored
!    INTEGER                                         :: level_old_id
    !!! index of the wps array where the level-set is backed up before adapt

    INTEGER                                         :: level_grad_id
    !!! index of the wps array where the gradient of the level-set is stored
!    INTEGER                                         :: level_grad_old_id
    !!! index of the wps array where the gradient of the level-set 
    !!! is backed up before adapt


    ! Anisotropic particles
    LOGICAL                                         :: anisotropic
    !!! true if the particles have their own cutoff radii
    !!! in this case, the G tensor will be stored in wpv(G_id)%vec
    INTEGER                                         :: G_id
    !!! index where G is stored


    ! spacings between particles
    LOGICAL                                         :: cartesian
    !!! true if particles are located exactly on the nodes of a uniform grid
    REAL(prec)                                      :: h_avg
    !!! (global) average distance between particles
    REAL(prec)                                      :: h_min
    !!! (global) minimum distance between particles

END TYPE ppm_t_particles

!----------------------------------------------------------------------
! Wrapper type to be able to have a pointer array to hold Particles
!----------------------------------------------------------------------
TYPE ppm_ptr_t_Particles
    TYPE(ppm_t_particles), POINTER  :: t => NULL()
END TYPE ppm_ptr_t_Particles

CONTAINS

FUNCTION get_xp(Particles,with_ghosts)
    TYPE(ppm_t_particles)            :: Particles
    LOGICAL,OPTIONAL                 :: with_ghosts
    REAL(prec),DIMENSION(:,:),POINTER:: get_xp

    IF (PRESENT(with_ghosts)) THEN
        IF (with_ghosts) THEN
            IF (Particles%xp_g.EQ.1) THEN
                get_xp => Particles%xp
            ELSE
                write(*,*) 'WARNING: tried to get xp with ghosts &
                    & when ghosts are not up-to-date'
                get_xp => NULL()
            ENDIF
            RETURN
        ENDIF
    ENDIF

    get_xp => Particles%xp

END FUNCTION get_xp

FUNCTION Set_xp(Particles,read_only,ghosts_ok)
    TYPE(ppm_t_particles)            :: Particles
    LOGICAL,OPTIONAL                 :: read_only
    LOGICAL,OPTIONAL                 :: ghosts_ok
    REAL(prec),DIMENSION(:,:),POINTER:: Set_xp

    IF (PRESENT(ghosts_ok)) THEN
        IF (ghosts_ok) THEN
            Set_xp => NULL()
            RETURN
        ENDIF
    ENDIF

    IF (PRESENT(read_only)) THEN
        IF (.NOT.read_only) THEN
            IF (Particles%xp_g.EQ.1) THEN
                Particles%xp_g = 0
            ENDIF
        ENDIF
    ELSE
        IF (Particles%xp_g.EQ.1) THEN
            Particles%xp_g = 0
        ENDIF
    ENDIF
    Set_xp => NULL()

END FUNCTION Set_xp

FUNCTION get_wps(Particles,wps_id,with_ghosts)
    TYPE(ppm_t_particles)            :: Particles
    INTEGER                          :: wps_id
    LOGICAL,OPTIONAL                 :: with_ghosts
    REAL(prec),DIMENSION(:),POINTER  :: get_wps

    IF (wps_id .LE. 0) THEN
        write(*,*) 'ERROR: failed to get wps for property &
            & wps_id = ',wps_id
        get_wps => NULL()
        RETURN
    ENDIF

    IF (wps_id .LE. Particles%max_wpsid) THEN
        IF (Particles%wps_m(wps_id).EQ.1) THEN
            IF (PRESENT(with_ghosts)) THEN
                IF (with_ghosts) THEN
                    IF (Particles%wps_g(wps_id).EQ.1) THEN
                        get_wps => Particles%wps(wps_id)%vec
                    ELSE
                        write(*,*) 'ERROR: tried to get wps with ghosts &
                            & when ghosts are not up-to-date. Returning NULL pointer'
                        write(*,*) 'Run with traceback option to debug'
                        get_wps => NULL()
                    ENDIF
                    RETURN
                ENDIF
            ENDIF
            get_wps => Particles%wps(wps_id)%vec
            RETURN
        ENDIF
    ENDIF
    write(*,*) 'ERROR: tried to get wps when mapping&
        & is not up-to-date. Returning NULL pointer'
    write(*,*) 'Run with traceback option to debug'
    get_wps => NULL()

END FUNCTION get_wps

FUNCTION set_wps(Particles,wps_id,read_only,ghosts_ok)
    TYPE(ppm_t_particles)            :: Particles
    INTEGER                          :: wps_id
    LOGICAL,OPTIONAL                 :: read_only
    LOGICAL,OPTIONAL                 :: ghosts_ok
    REAL(prec),DIMENSION(:),POINTER  :: set_wps

    !If read_only was not explicitely set to true, then assume
    !that ghosts are no longer up to date, unless ghosts_ok was
    ! explicitely set to true
    IF (PRESENT(ghosts_ok)) THEN
        IF (ghosts_ok) THEN
            set_wps => NULL()
            RETURN
        ENDIF
    ENDIF

    IF (PRESENT(read_only)) THEN
        IF (read_only) THEN
            set_wps => NULL()
            RETURN
        ENDIF
    ENDIF

    !Assume that the ghost values are now incorrect
    IF (Particles%wps_g(wps_id).EQ.1) THEN
        Particles%wps_g(wps_id) = 0
    ENDIF
    set_wps => NULL()

END FUNCTION set_wps

FUNCTION get_wpv(Particles,wpv_id,with_ghosts)
    TYPE(ppm_t_particles)            :: Particles
    INTEGER                          :: wpv_id
    LOGICAL,OPTIONAL                 :: with_ghosts
    REAL(prec),DIMENSION(:,:),POINTER:: get_wpv

    IF (wpv_id .LE. 0) THEN
        write(*,*) 'ERROR: failed to get wpv for property &
            & wpv_id = ',wpv_id
        get_wpv => NULL()
        RETURN
    ENDIF

    IF (wpv_id .LE. Particles%max_wpvid) THEN
        IF (Particles%wpv_m(wpv_id).EQ.1) THEN
            IF (PRESENT(with_ghosts)) THEN
                IF (with_ghosts) THEN
                    IF (Particles%wpv_g(wpv_id).EQ.1) THEN
                        get_wpv => Particles%wpv(wpv_id)%vec
                    ELSE
                        write(*,*) 'ERROR: tried to get wpv with ghosts &
                            & when ghosts are not up-to-date. Returning NULL pointer'
                        write(*,*) 'Run with traceback option to debug'
                        get_wpv => NULL()
                    ENDIF
                    RETURN
                ENDIF
            ENDIF
            get_wpv => Particles%wpv(wpv_id)%vec
            RETURN
        ENDIF
    ENDIF
    write(*,*) 'ERROR: tried to get wpv when mapping&
        & is not up-to-date. Returning NULL pointer'
    write(*,*) 'Run with traceback option to debug'
    get_wpv => NULL()

END FUNCTION get_wpv

FUNCTION set_wpv(Particles,wpv_id,read_only,ghosts_ok)
    TYPE(ppm_t_particles)            :: Particles
    INTEGER                          :: wpv_id
    LOGICAL,OPTIONAL                 :: read_only
    LOGICAL,OPTIONAL                 :: ghosts_ok
    REAL(prec),DIMENSION(:,:),POINTER:: set_wpv

    !If read_only was not explicitely set to true, then assume
    !that ghosts are no longer up to date, unless ghosts_ok
    ! was explicitely set to true

    IF (PRESENT(ghosts_ok)) THEN
        IF (ghosts_ok) THEN
            set_wpv => NULL()
            RETURN
        ENDIF
    ENDIF

    IF (PRESENT(read_only)) THEN
        IF (read_only) THEN
            set_wpv => NULL()
            RETURN
        ENDIF
    ENDIF

    !Assume that the ghost values are now incorrect
    IF (Particles%wpv_g(wpv_id).EQ.1) THEN
        Particles%wpv_g(wpv_id) = 0
    ENDIF
    set_wpv => NULL()

END FUNCTION set_wpv

SUBROUTINE ppm_alloc_particles(Particles,Npart,iopt,info)
    !!! (Re)allocates the memory of ppm_t_particles data type
    !-------------------------------------------------------------------------
    !  Includes
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    !  Modules
    !-------------------------------------------------------------------------
    USE ppm_module_substart
    USE ppm_module_substop
    USE ppm_module_error
#if   __KIND == __SINGLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    TYPE(ppm_t_particles), POINTER,  INTENT(INOUT)  :: Particles
    !!! Data structure containing the particles
    INTEGER,                         INTENT(IN   )  :: Npart
    !!! Number of particles 
    INTEGER                   ,      INTENT(IN   )  :: iopt
    !!! Allocation mode. One of:
    !!!
    !!! * ppm_param_alloc_fit
    !!! * ppm_param_dealloc
    INTEGER,                         INTENT(  OUT)  :: info
    !!! Returns status, 0 upon success.
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------

    INTEGER, DIMENSION(2)                           :: ldc
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
    ! maybe add some sanity checks later


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
            IF (ASSOCIATED(Particles%vlist)) DEALLOCATE(Particles%vlist,STAT=info)
            IF (ASSOCIATED(Particles%nvlist)) DEALLOCATE(Particles%nvlist,STAT=info)
            IF (ASSOCIATED(Particles%wps_g)) DEALLOCATE(Particles%wps_g,STAT=info)
            IF (ASSOCIATED(Particles%wps_m)) DEALLOCATE(Particles%wps_m,STAT=info)
            IF (ASSOCIATED(Particles%wpv_g)) DEALLOCATE(Particles%wpv_g,STAT=info)
            IF (ASSOCIATED(Particles%wpv_m)) DEALLOCATE(Particles%wpv_m,STAT=info)
            IF (ASSOCIATED(Particles%wpv_s)) DEALLOCATE(Particles%wpv_s,STAT=info)

            IF (ASSOCIATED(Particles%wps)) THEN
                DO i=1,Particles%max_wpsid
                    IF (ASSOCIATED(Particles%wps(i)%vec)) &
                        DEALLOCATE(Particles%wps(i)%vec,STAT=info)
                        NULLIFY(Particles%wps(i)%vec)
                ENDDO
                DEALLOCATE(Particles%wps,STAT=info)
            ENDIF
            IF (ASSOCIATED(Particles%wpv)) THEN
                DO i=1,Particles%max_wpvid
                    IF (ASSOCIATED(Particles%wpv(i)%vec)) &
                        DEALLOCATE(Particles%wpv(i)%vec,STAT=info)
                        NULLIFY(Particles%wpv(i)%vec)
                ENDDO
                DEALLOCATE(Particles%wpv,STAT=info)
            ENDIF
            DEALLOCATE(Particles,stat=info)
            NULLIFY(Particles)
            IF (info .NE. 0) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_dealloc,caller,   &
                    &          'Deallocating Particles',__LINE__,info)
            ENDIF
        ENDIF
    ENDIF

    !-------------------------------------------------------------------------
    !  Allocate new memory
    !-------------------------------------------------------------------------
    IF (lalloc) THEN

        ALLOCATE(Particles,STAT=info)
        IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,caller,   &
                &           'Allocating Particles',__LINE__,info)
            GOTO 9999
        ENDIF
        NULLIFY(Particles%xp)
        NULLIFY(Particles%nvlist)
        NULLIFY(Particles%vlist)
        NULLIFY(Particles%wps)
        NULLIFY(Particles%wps_g)
        NULLIFY(Particles%wps_m)
        NULLIFY(Particles%wpv)
        NULLIFY(Particles%wpv_g)
        NULLIFY(Particles%wpv_m)
        NULLIFY(Particles%wpv_s)
        !-----------------------------------------------------------------
        !  Allocate memory for the positions
        !-----------------------------------------------------------------
        ldc(1) = ppm_dim
        ldc(2) = Npart
        CALL ppm_alloc(Particles%xp,ldc,ppm_param_alloc_fit,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,   &
                &        'Could not allocate Particles elements',__LINE__,info)
            GOTO 9999
        ENDIF
        Particles%Npart = Npart
        Particles%Mpart = Npart
        Particles%xp_g = 0

        ! No properties defined yet
        Particles%nwps = 0
        Particles%nwpv = 0
        Particles%max_wpsid = 0
        Particles%max_wpvid = 0
        ! No active topology yet
        Particles%active_topoid = -1
        ! Particles are not yet mapped on any topology
        Particles%ontopology = .FALSE.
        Particles%xp_m = 0
        ! Neighbour lists are not yet computed
        Particles%neighlists = .FALSE.
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
        Particles%adapt_wpgradid = 0
        Particles%rcp_id = 0
        Particles%eta_id = 0
        Particles%D_id = 0
        Particles%Dtilde_id = 0
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


        IF (verbose) &
            write(*,*) 'Allocated particles with Np=',Npart


    ENDIF
    !-------------------------------------------------------------------------
    !  Return
    !-------------------------------------------------------------------------
    9999 CONTINUE
    CALL substop(caller,t0,info)
    RETURN
END SUBROUTINE ppm_alloc_particles

SUBROUTINE particles_allocate_wps(Particles,wp_id,info,with_ghosts,zero,iopt)
    USE ppm_module_substart
    USE ppm_module_substop
#if   __KIND == __SINGLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    TYPE(ppm_t_particles), POINTER,     INTENT(INOUT)      :: Particles
    !!! Data structure containing the particles
    INTEGER,                            INTENT(INOUT)      :: wp_id
    !!! id of this property
    INTEGER,                            INTENT(  OUT)      :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    !  Optional arguments
    !-------------------------------------------------------------------------
    INTEGER , OPTIONAL                                     :: iopt
    !!! Allocation mode. One of:
    !!!
    !!! * ppm_param_alloc_fit (default)
    !!! * ppm_param_alloc_grow
    !!! * ppm_param_alloc_grow_preserve
    !!! * ppm_param_dealloc
    LOGICAL , OPTIONAL                                     :: with_ghosts
    !!! if true, then allocate with Mpart instead of the default size of Npart
    LOGICAL , OPTIONAL                                     :: zero
    !!! if true, then initialise the array to zero
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    INTEGER , DIMENSION(3)    :: ldc
    !!! Number of elements in all dimensions for allocation
    LOGICAL                   :: lalloc,ldealloc,incr_nb_wp
    INTEGER                   :: iopt_l
    CHARACTER(LEN = ppm_char) :: caller = 'particles_allocate_wps'
    REAL(KIND(1.D0))          :: t0
    INTEGER                   :: i
    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)
    info = 0
    !-----------------------------------------------------------------
    !  Some checks
    !-----------------------------------------------------------------
    IF (.NOT.ASSOCIATED(Particles).OR..NOT.ASSOCIATED(Particles%xp)) THEN
        CALL ppm_error(ppm_err_alloc,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
    IF(.NOT.ASSOCIATED(Particles%wps)) THEN
        !allocate memory
        ALLOCATE(Particles%wps(20),STAT=info)
        ALLOCATE(Particles%wps_g(20),STAT=info)
        ALLOCATE(Particles%wps_m(20),STAT=info)
        !set defaults
        Particles%wps_g = -1
        Particles%wps_m = -1
    ENDIF
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &            'allocation failed for wps',__LINE__,info)
        GOTO 9999
    ENDIF

    !-------------------------------------------------------------------------
    !  Check the allocation type
    !-------------------------------------------------------------------------
    IF (PRESENT(iopt)) THEN
        iopt_l = iopt
    ELSE
        iopt_l = ppm_param_alloc_fit
    ENDIF
    lalloc   = .FALSE.
    ldealloc = .FALSE.

    IF (iopt_l.EQ.ppm_param_alloc_fit) THEN
        !----------------------------------------------------------------------
        !  fit memory but skip the present contents
        !----------------------------------------------------------------------
        IF (wp_id.NE.0) THEN
            IF (ASSOCIATED(Particles%wps(wp_id)%vec)) ldealloc = .TRUE.
        ENDIF
        lalloc   = .TRUE.
    ELSEIF (iopt_l.EQ.ppm_param_alloc_grow_preserve) THEN
        lalloc = .TRUE.
    ELSEIF (iopt_l.EQ.ppm_param_alloc_grow) THEN
        lalloc = .TRUE.
    ELSEIF (iopt_l.EQ.ppm_param_dealloc) THEN
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

    IF (ldealloc) THEN
        ldc = 0
        CALL ppm_alloc(Particles%wps(wp_id)%vec,ldc,ppm_param_dealloc,info)
        IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_dealloc,caller,&
                'ppm_alloc (deallocate) failed',__LINE__,info)
            GOTO 9999
        ENDIF

        !----------------------------------------------------------------------
        ! Update state
        !----------------------------------------------------------------------
        Particles%wps_m(wp_id) = -1
        Particles%wps_g(wp_id) = -1
        Particles%wps(wp_id)%vec => NULL()
        !Decrement number of properties
        Particles%nwps = Particles%nwps-1
        IF (wp_id .GE. Particles%max_wpsid) THEN
            Particles%max_wpsid = Particles%max_wpsid - 1
            DO WHILE (.NOT.ASSOCIATED(Particles%wps(Particles%max_wpsid)%vec) )
                Particles%max_wpsid = Particles%max_wpsid - 1
                IF (Particles%max_wpsid.LE.0) EXIT
            ENDDO
        ENDIF

        IF (verbose) &
            write(*,*) 'De-allocated wps with id=', wp_id
        IF (.NOT. lalloc) THEN
            wp_id = 0
        ENDIF
    ENDIF

    IF (lalloc) THEN
        !-----------------------------------------------------------------
        !  If not provided by the user, generate new id for the property
        !-----------------------------------------------------------------
        IF (wp_id.EQ.0) THEN
            !create a new property id (first free index)
            wp_id = 1
            DO WHILE (ASSOCIATED(Particles%wps(wp_id)%vec))
                wp_id = wp_id + 1
                IF (wp_id.GT.SIZE(Particles%wps,1)) THEN
                    write(*,*) 'need to increase the size of the property array'
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_alloc,caller,   &
                        &            'allocation failed for wps',__LINE__,info)
                    GOTO 9999
                ENDIF
            ENDDO
            incr_nb_wp = .TRUE.
        ELSE
            IF (ldealloc) THEN
                !number of properties has decreased in dealloc and
                ! needs to be incremented back again
                incr_nb_wp = .TRUE.
            ELSE
                !number of properties will not change
                incr_nb_wp = .FALSE.
            ENDIF
        ENDIF

        !-----------------------------------------------------------------
        !  Allocate memory for the properties
        !-----------------------------------------------------------------
        ldc(1) = Particles%Npart
        IF (PRESENT(with_ghosts)) THEN
            IF (with_ghosts) ldc(1) = Particles%Mpart
        ENDIF
        CALL ppm_alloc(Particles%wps(wp_id)%vec,ldc,iopt_l,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,   &
                &            'Could not allocate ppm_topo elements',__LINE__,info)
            GOTO 9999
        ENDIF

        IF (PRESENT(zero)) THEN
            IF (zero) THEN
                DO i=1,ldc(1)
                    Particles%wps(wp_id)%vec(i)=0._MK
                ENDDO
            ENDIF
        ENDIF

        !-----------------------------------------------------------------
        ! Update state
        !-----------------------------------------------------------------
        !Set its state to "mapped" (every index corresponds to exactly one
        !real particle)
        Particles%wps_m(wp_id) = 1

        !if Mpart was up-to-date, then set the for this property ghosts to "updated"
        ! (because to every index between Npart+1 and Mpart corresponds exactly
        ! one ghost value)
        Particles%wps_g(wp_id) = 0
        IF (PRESENT(with_ghosts)) THEN
            IF (with_ghosts .AND. Particles%xp_g .EQ. 1) THEN
                Particles%wps_g(wp_id) = 1
            ENDIF
        ENDIF

        IF (incr_nb_wp) THEN
            !Increment number of properties
            Particles%nwps = Particles%nwps+1
        ENDIF

        IF (wp_id .GT. Particles%max_wpsid) THEN
            Particles%max_wpsid = wp_id
        ENDIF
        IF (Particles%max_wpsid .GT. SIZE(Particles%wps)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,&
                'Pointer array wps is growing and should be reallocated',&
                __LINE__,info)
            GOTO 9999
        ENDIF

        !We assume that if the user allocates an array for the cutoff radius,
        ! then it means the particles are meant to be adaptive.
        IF (Particles%rcp_id.NE.0) Particles%adaptive=.TRUE.

        IF (verbose) &
            write(*,*) 'Allocated wps with size=',ldc(1),' id=',wp_id
    ENDIF

    9999 CONTINUE

END SUBROUTINE particles_allocate_wps


SUBROUTINE particles_allocate_wpv(Particles,wp_id,lda,info, & 
        with_ghosts,zero,iopt)
    USE ppm_module_substart
    USE ppm_module_substop
#if   __KIND == __SINGLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    TYPE(ppm_t_particles), POINTER,     INTENT(INOUT)      :: Particles
    !!! Data structure containing the particles
    INTEGER,                            INTENT(INOUT)      :: wp_id
    !!! id of this property
    INTEGER,                            INTENT(IN   )      :: lda
    !!! leading dimension of the property array
    INTEGER,                            INTENT(  OUT)      :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    !  Optional arguments
    !-------------------------------------------------------------------------
    INTEGER , OPTIONAL                                     :: iopt
    !!! Allocation mode. One of:
    !!!
    !!! * ppm_param_alloc_fit (default)
    !!! * ppm_param_alloc_grow
    !!! * ppm_param_alloc_grow_preserve
    !!! * ppm_param_dealloc
    LOGICAL , OPTIONAL                                     :: with_ghosts
    !!! if true, then allocate with Mpart instead of the default size of Npart
    LOGICAL , OPTIONAL                                     :: zero
    !!! if true, then initialise the array to zero
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    INTEGER , DIMENSION(3)    :: ldc
    !!! Number of elements in all dimensions for allocation
    LOGICAL                   :: lalloc,ldealloc,incr_nb_wp
    INTEGER                   :: iopt_l
    !!! allocation mode, see one of ppm_alloc_* subroutines.
    CHARACTER(LEN = ppm_char) :: caller = 'particles_allocate_wpv'
    REAL(KIND(1.D0))          :: t0
    INTEGER                   :: i
    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)
    info = 0
    !-----------------------------------------------------------------
    !  Checks
    !-----------------------------------------------------------------
    IF (.NOT.ASSOCIATED(Particles).OR..NOT.ASSOCIATED(Particles%xp)) THEN
        CALL ppm_error(ppm_err_alloc,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
    IF(.NOT.ASSOCIATED(Particles%wpv)) THEN
        !allocate memory
        ALLOCATE(Particles%wpv(20),STAT=info)
        ALLOCATE(Particles%wpv_g(20),STAT=info)
        ALLOCATE(Particles%wpv_m(20),STAT=info)
        ALLOCATE(Particles%wpv_s(20),STAT=info)
        !set defaults
        Particles%wpv_g = -1
        Particles%wpv_m = -1
        Particles%wpv_s = 0
    ENDIF
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &            'allocation failed for wpv',__LINE__,info)
        GOTO 9999
    ENDIF


    !-------------------------------------------------------------------------
    !  Check the allocation type
    !-------------------------------------------------------------------------
    IF (PRESENT(iopt)) THEN
        iopt_l = iopt
    ELSE
        iopt_l = ppm_param_alloc_fit
    ENDIF
    lalloc   = .FALSE.
    ldealloc = .FALSE.
    IF (iopt_l.EQ.ppm_param_alloc_fit) THEN
        !----------------------------------------------------------------------
        !  fit memory but skip the present contents
        !----------------------------------------------------------------------
        IF (wp_id.NE.0) THEN
            IF (ASSOCIATED(Particles%wpv(wp_id)%vec)) ldealloc = .TRUE.
        ENDIF
        lalloc   = .TRUE.
    ELSEIF (iopt_l.EQ.ppm_param_alloc_grow_preserve) THEN
        lalloc = .TRUE.
    ELSEIF (iopt_l.EQ.ppm_param_alloc_grow) THEN
        lalloc = .TRUE.
    ELSEIF (iopt_l.EQ.ppm_param_dealloc) THEN
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

    IF (ldealloc) THEN
        ldc(1) = 0
        ldc(2) = 0
        CALL ppm_alloc(Particles%wpv(wp_id)%vec,ldc,ppm_param_dealloc,info)
        IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_dealloc,caller,&
                'ppm_alloc (deallocate) failed',__LINE__,info)
            GOTO 9999
        ENDIF

        ! Update state
        Particles%wpv_m(wp_id) = -1
        Particles%wpv_g(wp_id) = -1
        Particles%wpv_s(wp_id) = 0
        Particles%wpv(wp_id)%vec => NULL()
        !Decrement number of properties
        Particles%nwpv = Particles%nwpv-1
        IF (wp_id .GE. Particles%max_wpvid) THEN
            Particles%max_wpvid = Particles%max_wpvid - 1
            DO WHILE (.NOT.ASSOCIATED(Particles%wpv(Particles%max_wpvid)%vec) )
                Particles%max_wpvid = Particles%max_wpvid - 1
                IF (Particles%max_wpvid.LE.0) EXIT
            ENDDO
        ENDIF

        IF (verbose) &
            write(*,*) 'De-allocated wpv with id=', wp_id
        IF (.NOT. lalloc) THEN
            wp_id = 0
        ENDIF
    ENDIF

    IF (lalloc) THEN
        !-----------------------------------------------------------------
        !  If not provided by the user, generate new id for the property
        !-----------------------------------------------------------------
        IF (wp_id.EQ.0) THEN
            !create a new property id (first free index)
            wp_id = 1
            DO WHILE (ASSOCIATED(Particles%wpv(wp_id)%vec))
                wp_id = wp_id + 1
                IF (wp_id.GT.SIZE(Particles%wpv,1)) THEN
                    write(*,*) 'need to increase the size of the property array'
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_alloc,caller,   &
                        &            'allocation failed for wpv',__LINE__,info)
                    GOTO 9999
                ENDIF
            ENDDO
            incr_nb_wp = .TRUE.
        ELSE
            !number of properties will not change
            incr_nb_wp = .FALSE.
        ENDIF

        !-----------------------------------------------------------------
        !  Allocate memory for the properties
        !-----------------------------------------------------------------
        iopt_l = ppm_param_alloc_fit
        ldc(1) = lda
        ldc(2) = Particles%Npart
        IF (PRESENT(with_ghosts)) THEN
            IF (with_ghosts) ldc(2) = Particles%Mpart
        ENDIF
        CALL ppm_alloc(Particles%wpv(wp_id)%vec,ldc(1:2),iopt_l,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,   &
                &            'Could not allocate ppm_topo elements',__LINE__,info)
            GOTO 9999
        ENDIF
        IF (PRESENT(zero)) THEN
            IF (zero) THEN
                DO i=1,ldc(2)
                    Particles%wpv(wp_id)%vec(:,i)=0._MK
                ENDDO
            ENDIF
        ENDIF

        ! Update state
        !Set its state to "mapped" (every index corresponds to exactly one
        !real particle)
        Particles%wpv_m(wp_id) = 1
        !if Mpart was up-to-date, then set the state for this property ghosts to 
        ! "updated" (because to every index between Npart+1 and Mpart 
        ! corresponds exactly one ghost value)
        Particles%wpv_g(wp_id) = 0
        IF (PRESENT(with_ghosts)) THEN
            IF (with_ghosts .AND. Particles%xp_g .EQ. 1) THEN
                Particles%wpv_g(wp_id) = 1
            ENDIF
        ENDIF
        !Set its size to lda
        Particles%wpv_s(wp_id) = lda

        IF (incr_nb_wp) THEN
            !Increment number of properties
            Particles%nwpv = Particles%nwpv+1
        ENDIF

        IF (wp_id .GT. Particles%max_wpvid) THEN
            Particles%max_wpvid = wp_id
        ENDIF

        IF (Particles%max_wpvid .GT. SIZE(Particles%wpv)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,&
                'Pointer array wpv is growing and should be reallocated',&
                __LINE__,info)
            GOTO 9999
        ENDIF

        IF (verbose) &
            write(*,*) 'Allocated wpv with size=',ldc,' id=',wp_id

    ENDIF

    !-----------------------------------------------------------------------
    ! Finalize
    !-----------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999 CONTINUE ! jump here upon error

END SUBROUTINE particles_allocate_wpv

SUBROUTINE particles_mapping_global(Particles,topoid,info)
    !-----------------------------------------------------------------
    !  Global mapping for particles
    !-----------------------------------------------------------------
    !  Assumptions:
    ! * All the particles have to be inside the domain
    !   (otherwise -> "unassigned particle error")
    !
    USE ppm_module_map
    USE ppm_module_substart
    USE ppm_module_substop
#if   __KIND == __SINGLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    TYPE(ppm_t_particles), POINTER,  INTENT(INOUT)      :: Particles
    !!! Data structure containing the particles
    INTEGER,                            INTENT(IN   )      :: topoid
    !!! topology id
    INTEGER,                            INTENT(  OUT)      :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    INTEGER , DIMENSION(3)    :: ldc
    !!! Number of elements in all dimensions for allocation
    INTEGER                   :: iopt
    !!! allocation mode, see one of ppm_alloc_* subroutines.
    INTEGER                   :: Npart_new
    !!! new number of particles on this processor
    INTEGER                   :: prop_id
    !!! index variable
    CHARACTER(LEN = ppm_char) :: caller = 'particles_mapping_global'
    REAL(KIND(1.D0))          :: t0
    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)
    !-----------------------------------------------------------------
    !  Checks
    !-----------------------------------------------------------------
    IF (.NOT.ASSOCIATED(Particles).OR..NOT.ASSOCIATED(Particles%xp)) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
    IF (.NOT.Particles%areinside) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'some particles may be outside the domain. Apply BC first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    !-----------------------------------------------------------------------
    !  Map the particles onto the topology
    !-----------------------------------------------------------------------
    IF (Particles%ontopology.AND.Particles%active_topoid.EQ.topoid) THEN
        !Particles have already been mapped onto this topology
        !nothing to do
    ELSE
        CALL ppm_map_part_global(topoid,Particles%xp,Particles%Npart,info) 
        IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(999,caller,&
                'ppm_map_part_global failed',__LINE__,info)
            GOTO 9999
        ENDIF

        DO prop_id = 1,Particles%max_wpsid
            IF(Particles%wps_m(prop_id) .EQ. 0) THEN
                CALL ppm_map_part_push(Particles%wps(prop_id)%vec,&
                    Particles%Npart,info)
                IF (info .NE. 0) THEN
                    info = ppm_error_fatal
                    CALL ppm_error(999,caller,&
                        'ppm_map_part_push failed',__LINE__,info)
                    GOTO 9999
                ENDIF
            ENDIF
        ENDDO

        DO prop_id = 1,Particles%max_wpvid
            IF(Particles%wpv_m(prop_id) .EQ. 0) THEN
                CALL ppm_map_part_push(Particles%wpv(prop_id)%vec,&
                    Particles%wpv_s(prop_id), Particles%Npart,info)
                IF (info .NE. 0) THEN
                    info = ppm_error_fatal
                    CALL ppm_error(999,caller,&
                        'ppm_map_part_push failed',__LINE__,info)
                    GOTO 9999
                ENDIF
            ENDIF
        ENDDO

        CALL ppm_map_part_send(Particles%Npart,Npart_new,info)
        IF (info .NE. 0) THEN
            CALL ppm_error(0,caller,&
                'ppm_map_part_send failed',__LINE__,info)
            GOTO 9999
        ENDIF

        DO prop_id = Particles%max_wpvid,1,-1
            IF(Particles%wpv_m(prop_id) .EQ. 0) THEN
                CALL ppm_map_part_pop(Particles%wpv(prop_id)%vec,&
                    Particles%wpv_s(prop_id), Particles%Npart,Npart_new,info)
                IF (info .NE. 0) THEN
                    info = ppm_error_fatal
                    CALL ppm_error(999,caller,&
                        'ppm_map_part_pop failed',__LINE__,info)
                    GOTO 9999
                ENDIF
            ENDIF
        ENDDO
        DO prop_id = Particles%max_wpsid,1,-1
            IF(Particles%wps_m(prop_id) .EQ. 0) THEN
                CALL ppm_map_part_pop(Particles%wps(prop_id)%vec,&
                    Particles%Npart,Npart_new,info)
                IF (info .NE. 0) THEN
                    info = ppm_error_fatal
                    CALL ppm_error(999,caller,&
                        'ppm_map_part_pop failed',__LINE__,info)
                    GOTO 9999
                ENDIF
            ENDIF
        ENDDO

        CALL ppm_map_part_pop(Particles%xp,ppm_dim,Particles%Npart,&
            Npart_new,info)
        IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(999,caller,&
                'ppm_map_part_pop failed',__LINE__,info)
            GOTO 9999
        ENDIF

        ! Update states
        ! Number of particles on this processor
        Particles%Npart = Npart_new
        ! This is the active topology for these particles
        Particles%active_topoid = topoid
        ! Particles are now mapped on the active topology
        Particles%ontopology = .TRUE.
        Particles%xp_m = 1
        !   values for some scalar arrays have been mapped and ghosts
        !   are no longer up-to-date
        DO prop_id = 1,Particles%max_wpsid
            IF (Particles%wps_m(prop_id) .EQ. 0) Particles%wps_m(prop_id) = 1
            IF (Particles%wps_g(prop_id) .EQ. 1) Particles%wps_g(prop_id) = 0
        ENDDO
        !   values for some vector arrays have been mapped and ghosts 
        !   are no longer up-to-date
        DO prop_id = 1,Particles%max_wpvid
            IF (Particles%wpv_m(prop_id) .EQ. 0) Particles%wpv_m(prop_id) = 1
            IF (Particles%wpv_g(prop_id) .EQ. 1) Particles%wpv_g(prop_id) = 0
        ENDDO
        ! particles have been re-indexed and ghosts have not been computed
        Particles%xp_g = 0
        Particles%Mpart = Particles%Npart
        ! particles have been re-indexed and neighbour lists not updated
        Particles%neighlists = .FALSE.

    IF (verbose) &
        write(*,*) 'Finished Global mapping'

    ENDIF

    !-----------------------------------------------------------------------
    ! Finalize
    !-----------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999 CONTINUE ! jump here upon error

END SUBROUTINE particles_mapping_global

SUBROUTINE particles_mapping_partial(Particles,topoid,info,debug)
    !-----------------------------------------------------------------
    !  Partial mapping for particles
    !-----------------------------------------------------------------
    !  Assumptions:
    ! * All the particles have to be inside the domain
    !   (otherwise -> "unassigned particle error")
    !
    USE ppm_module_map
    USE ppm_module_substart
    USE ppm_module_substop
#if   __KIND == __SINGLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    TYPE(ppm_t_particles), POINTER,  INTENT(INOUT)      :: Particles
    !!! Data structure containing the particles
    INTEGER,                            INTENT(IN   )      :: topoid
    !!! topology id
    INTEGER,                            INTENT(  OUT)      :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    !  Optional Arguments
    !-------------------------------------------------------------------------
    LOGICAL, OPTIONAL                                   :: debug
    !!! IF true, printout more
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    INTEGER , DIMENSION(3)    :: ldc
    !!! Number of elements in all dimensions for allocation
    INTEGER                   :: iopt
    !!! allocation mode, see one of ppm_alloc_* subroutines.
    INTEGER                   :: Npart_new
    !!! new number of particles on this processor
    INTEGER                   :: prop_id
    !!! index variable
    CHARACTER(LEN = ppm_char) :: caller = 'particles_mapping_partial'
    REAL(KIND(1.D0))          :: t0
    LOGICAL                   :: dbg
    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)
    dbg = .FALSE.
    IF (PRESENT(debug)) dbg = debug
    !-----------------------------------------------------------------
    !  Checks
    !-----------------------------------------------------------------
    IF (.NOT.ASSOCIATED(Particles).OR..NOT.ASSOCIATED(Particles%xp)) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
    IF (.NOT.Particles%areinside) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'some particles may be outside the domain. Apply BC first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    !-----------------------------------------------------------------------
    !  Map the particles onto the topology
    !-----------------------------------------------------------------------
    IF (Particles%ontopology.AND.Particles%active_topoid.EQ.topoid) THEN
        !Particles have already been mapped onto this topology
        !nothing to do
        info = ppm_error_notice
        CALL ppm_error(999,caller,   &
            'Particles have not moved since last mapping - not doing anything',&
            &  __LINE__,info)
    ELSE
        CALL ppm_map_part_partial(topoid,Particles%xp,Particles%Npart,info) 
        IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(999,caller,&
                'ppm_map_part_global failed',__LINE__,info)
            GOTO 9999
        ENDIF

        DO prop_id = 1,Particles%max_wpsid
            IF(Particles%wps_m(prop_id) .EQ. 0) THEN
                IF(dbg) &
                    write(*,*) 'pushing-1 ',prop_id
                CALL ppm_map_part_push(Particles%wps(prop_id)%vec,&
                    Particles%Npart,info)
                IF (info .NE. 0) THEN
                    info = ppm_error_fatal
                    CALL ppm_error(999,caller,&
                        'ppm_map_part_push failed',__LINE__,info)
                    GOTO 9999
                ENDIF
            ENDIF
        ENDDO

        DO prop_id = 1,Particles%max_wpvid
            IF(Particles%wpv_m(prop_id) .EQ. 0) THEN
                IF(dbg) &
                    write(*,*) 'pushing-2 ',prop_id
                CALL ppm_map_part_push(Particles%wpv(prop_id)%vec,&
                    Particles%wpv_s(prop_id), Particles%Npart,info)
                IF (info .NE. 0) THEN
                    info = ppm_error_fatal
                    CALL ppm_error(999,caller,&
                        'ppm_map_part_push failed',__LINE__,info)
                    GOTO 9999
                ENDIF
            ENDIF
        ENDDO

        CALL ppm_map_part_send(Particles%Npart,Npart_new,info)
        IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(999,caller,&
                'ppm_map_part_send failed',__LINE__,info)
            GOTO 9999
        ENDIF

        DO prop_id = Particles%max_wpvid,1,-1
            IF(Particles%wpv_m(prop_id) .EQ. 0) THEN
                IF(dbg) &
                    write(*,*) 'popping-2 ',prop_id
                CALL ppm_map_part_pop(Particles%wpv(prop_id)%vec,&
                    Particles%wpv_s(prop_id), Particles%Npart,Npart_new,info)
                IF (info .NE. 0) THEN
                    info = ppm_error_fatal
                    CALL ppm_error(999,caller,&
                        'ppm_map_part_pop failed',__LINE__,info)
                    GOTO 9999
                ENDIF
            ENDIF
        ENDDO
        DO prop_id = Particles%max_wpsid,1,-1
            IF(Particles%wps_m(prop_id) .EQ. 0) THEN
                IF(dbg) &
                    write(*,*) 'popping-1 ',prop_id
                CALL ppm_map_part_pop(Particles%wps(prop_id)%vec,&
                    Particles%Npart,Npart_new,info)
                IF (info .NE. 0) THEN
                    info = ppm_error_fatal
                    CALL ppm_error(999,caller,&
                        'ppm_map_part_pop failed',__LINE__,info)
                    GOTO 9999
                ENDIF
            ENDIF
        ENDDO

        CALL ppm_map_part_pop(Particles%xp,ppm_dim,Particles%Npart,&
            Npart_new,info)
        IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(999,caller,&
                'ppm_map_part_pop failed',__LINE__,info)
            GOTO 9999
        ENDIF

        ! Update states
        ! Number of particles on this processor
        Particles%Npart = Npart_new
        ! This is the active topology for these particles
        Particles%active_topoid = topoid
        ! Particles are now mapped on the active topology
        Particles%ontopology = .TRUE.
        Particles%xp_m = 1
        !   values for some scalar arrays have been mapped and ghosts
        !   are no longer up-to-date
        DO prop_id = 1,Particles%max_wpsid
            IF (Particles%wps_m(prop_id) .EQ. 0) Particles%wps_m(prop_id) = 1
            IF (Particles%wps_g(prop_id) .EQ. 1) Particles%wps_g(prop_id) = 0
        ENDDO
        !   values for some vector arrays have been mapped and ghosts 
        !   are no longer up-to-date
        DO prop_id = 1,Particles%max_wpvid
            IF (Particles%wpv_m(prop_id) .EQ. 0) Particles%wpv_m(prop_id) = 1
            IF (Particles%wpv_g(prop_id) .EQ. 1) Particles%wpv_g(prop_id) = 0
        ENDDO
        ! particles have been re-indexed and ghosts have not been computed
        Particles%xp_g = 0
        Particles%Mpart = Particles%Npart
        ! particles have been re-indexed and neighbour lists not updated
        Particles%neighlists = .FALSE.

    ENDIF

    IF (verbose) &
        write(*,*) 'Finished partial mapping'

    !-----------------------------------------------------------------------
    ! Finalize
    !-----------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999 CONTINUE ! jump here upon error

END SUBROUTINE particles_mapping_partial

SUBROUTINE particles_mapping_ghosts(Particles,topoid,info,debug)
    !-----------------------------------------------------------------
    !  Ghost mapping for particles
    !-----------------------------------------------------------------
    !  Assumptions:
    ! * Particles positions need to have been mapped onto the topology
    !
    USE ppm_module_data, ONLY: ppm_topo
    USE ppm_module_map
    USE ppm_module_substart
    USE ppm_module_substop
#if   __KIND == __SINGLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    TYPE(ppm_t_particles), POINTER,  INTENT(INOUT)      :: Particles
    !!! Data structure containing the particles
    INTEGER,                            INTENT(IN   )   :: topoid
    !!! topology id
    INTEGER,                            INTENT(  OUT)   :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    !  Optional Arguments
    !-------------------------------------------------------------------------
    LOGICAL, OPTIONAL                                   :: debug
    !!! IF true, printout more
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    INTEGER , DIMENSION(3)                    :: ldc
    !!! Number of elements in all dimensions for allocation
    INTEGER                                   :: iopt
    !!! allocation mode, see one of ppm_alloc_* subroutines.
    INTEGER                                   :: prop_id
    !!! index variable
    REAL(ppm_kind_double)                     :: cutoff
    !!! cutoff radius
    TYPE(ppm_t_topo),POINTER                  :: topo => NULL()
    CHARACTER(LEN = ppm_char) :: caller = 'particles_mapping_ghosts'
    REAL(KIND(1.D0))                          :: t0
    LOGICAL                                   :: dbg
    LOGICAL                                   :: skip_ghost_get
    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)
    dbg = .FALSE.
    IF (PRESENT(debug)) dbg=debug
    skip_ghost_get = .FALSE.
    !-----------------------------------------------------------------
    !  Checks
    !-----------------------------------------------------------------
    IF (.NOT.ASSOCIATED(Particles).OR..NOT.ASSOCIATED(Particles%xp)) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (.NOT.Particles%ontopology .OR. Particles%active_topoid.NE.topoid) THEN
        !Particles have not been mapped onto this topology
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'Do a partial/global mapping before doing a ghost mapping',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
    IF (Particles%xp_m .NE. 1) THEN
        !Particles have not been mapped onto this topology
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'Real particles are not mapped one-to-one with xp',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
    IF (.NOT.Particles%areinside) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'some particles may be outside the domain. Apply BC first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF


    topo=>ppm_topo(topoid)%t

    cutoff = Particles%cutoff + Particles%skin

#if   __KIND == __SINGLE_PRECISION
    IF (cutoff .GT. topo%ghostsizes) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'ghostsize of topology may be smaller than that of particles',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
#elif   __KIND == __DOUBLE_PRECISION
    IF (cutoff .GT. topo%ghostsized) THEN

        write(*,*) cutoff
        write(*,*) topo%ghostsized

        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'ghostsize of topology may be smaller than that of particles',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
#endif
    IF (cutoff .GT. 0.0D0) THEN
        IF (Particles%xp_g .EQ. 1) THEN
            IF (verbose) THEN
                write(*,*) 'ghosts have already been updated'
            ENDIF

            IF (ppm_map_type_isactive(ppm_param_map_ghost_get)) THEN
                IF (verbose) THEN
                    write(*,*) 'we skip the ghost_get and go straight to'
                    write(*,*) 'push/send/pop'
                ENDIF
                skip_ghost_get = .TRUE.
            ENDIF
        ENDIF

        IF (.NOT.skip_ghost_get) THEN
            CALL ppm_map_part_ghost_get(topoid,Particles%xp,ppm_dim,&
                Particles%Npart,Particles%isymm,cutoff,info)
            IF (info .NE. 0) THEN
                info = ppm_error_fatal
                CALL ppm_error(999,caller,&
                    'ppm_map_part_ghost_get failed',__LINE__,info)
                GOTO 9999
            ENDIF
        ENDIF

        !Update the ghost for the properties if
        ! 1) they have been mapped to this topology,
        ! 2) the ghosts have not yet been updated, and
        ! 3) the user wants them to be updated
        DO prop_id = 1,Particles%max_wpsid
            IF(Particles%wps_m(prop_id) .EQ. 1) THEN
                IF(Particles%wps_g(prop_id) .EQ. 0) THEN
                    IF(dbg) &
                        write(*,*) 'pushing-1 ',prop_id
                    CALL ppm_map_part_push(Particles%wps(prop_id)%vec,&
                        Particles%Npart,info)
                    IF (info .NE. 0) THEN
                        info = ppm_error_fatal
                        CALL ppm_error(999,caller,&
                            'ppm_map_part_push failed',__LINE__,info)
                        GOTO 9999
                    ENDIF
                ENDIF
            ENDIF
        ENDDO
        DO prop_id = 1,Particles%max_wpvid
            IF(Particles%wpv_m(prop_id) .EQ. 1) THEN
                IF(Particles%wpv_g(prop_id) .EQ. 0) THEN
                    IF (dbg) &
                        write(*,*) 'pushing-2 ',prop_id
                    CALL ppm_map_part_push(Particles%wpv(prop_id)%vec,&
                        Particles%wpv_s(prop_id), Particles%Npart,info)
                    IF (info .NE. 0) THEN
                        info = ppm_error_fatal
                        CALL ppm_error(999,caller,&
                            'ppm_map_part_push failed',__LINE__,info)
                        GOTO 9999
                    ENDIF
                ENDIF
            ENDIF
        ENDDO

        CALL ppm_map_part_send(Particles%Npart,Particles%Mpart,info)
        IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(999,caller,&
                'ppm_map_part_send failed',__LINE__,info)
            GOTO 9999
        ENDIF

        DO prop_id = Particles%max_wpvid,1,-1
            IF(Particles%wpv_m(prop_id) .EQ. 1) THEN
                IF(Particles%wpv_g(prop_id) .EQ. 0) THEN
                    IF(dbg) &
                        write(*,*) 'popping-2 ',prop_id
                    CALL ppm_map_part_pop(Particles%wpv(prop_id)%vec,&
                        Particles%wpv_s(prop_id),Particles%Npart, &
                        Particles%Mpart,info)
                    IF (info .NE. 0) THEN
                        info = ppm_error_fatal
                        CALL ppm_error(999,caller,&
                            'ppm_map_part_pop failed',__LINE__,info)
                        GOTO 9999
                    ENDIF
                ENDIF
            ENDIF
        ENDDO
        DO prop_id = Particles%max_wpsid,1,-1
            IF(Particles%wps_m(prop_id) .EQ. 1) THEN
                IF(Particles%wps_g(prop_id) .EQ. 0) THEN
                    IF(dbg) &
                        write(*,*) 'popping-1 ',prop_id
                    CALL ppm_map_part_pop(Particles%wps(prop_id)%vec,&
                        Particles%Npart,Particles%Mpart,info)
                    IF (info .NE. 0) THEN
                        info = ppm_error_fatal
                        CALL ppm_error(999,caller,&
                            'ppm_map_part_pop failed',__LINE__,info)
                        GOTO 9999
                    ENDIF
                ENDIF
            ENDIF
        ENDDO

        IF (.NOT.skip_ghost_get) THEN
            CALL ppm_map_part_pop(Particles%xp,ppm_dim,Particles%Npart,&
                Particles%Mpart,info)
            IF (info .NE. 0) THEN
                info = ppm_error_fatal
                CALL ppm_error(999,caller,&
                    'ppm_map_part_pop failed',__LINE__,info)
                GOTO 9999
            ENDIF
        ENDIF
    ENDIF


    ! Update states
    !   ghosts have been computed
    Particles%xp_g = 1
    !   ghosts values for some scalar arrays have been computed
    DO prop_id = 1,Particles%max_wpsid
        IF (Particles%wps_m(prop_id) .EQ. 1) THEN
            IF (Particles%wps_g(prop_id) .EQ. 0) Particles%wps_g(prop_id) = 1
        ENDIF
    ENDDO
    !   ghosts values for some vector arrays have been computed
    DO prop_id = 1,Particles%max_wpvid
        IF (Particles%wpv_m(prop_id) .EQ. 1) THEN
            IF (Particles%wpv_g(prop_id) .EQ. 0) Particles%wpv_g(prop_id) = 1
        ENDIF
    ENDDO

    IF (verbose) &
        write(*,*) 'Finished ghost mapping'
    !-----------------------------------------------------------------------
    ! Finalize
    !-----------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999 CONTINUE ! jump here upon error

END SUBROUTINE particles_mapping_ghosts

SUBROUTINE particles_apply_bc(Particles,topoid,info)
    !-----------------------------------------------------------------
    !  Apply boundary conditions for particles positions
    !-----------------------------------------------------------------
    !  Assumptions:
    ! * Particles positions need to have been mapped onto the topology
    !
    USE ppm_module_data, ONLY: ppm_topo,ppm_rank
    USE ppm_module_substart
    USE ppm_module_substop
#if   __KIND == __SINGLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    TYPE(ppm_t_particles), POINTER,  INTENT(INOUT)        :: Particles
    !!! Data structure containing the particles
    INTEGER,                            INTENT(IN   )     :: topoid
    !!! topology id
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
    INTEGER                                               :: di,ip,prop_id
    INTEGER                                               :: Npart,del_part
    REAL(ppm_kind_double)                                 :: t0
    CHARACTER(LEN = ppm_char)                   :: caller = 'particles_apply_bc'
    REAL(MK)                                              :: almostone
    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)
    !-----------------------------------------------------------------
    !  checks
    !-----------------------------------------------------------------
    IF (.NOT.ASSOCIATED(Particles).OR..NOT.ASSOCIATED(Particles%xp)) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (Particles%active_topoid.NE.topoid) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'WARNING: this topoid is not the one for which &
            & these particles were mapped',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (Particles%xp_m .NE. 1) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'Real particles are not mapped one-to-one with xp',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    topo=>ppm_topo(topoid)%t
    xp=>Particles%xp
    Npart = Particles%Npart
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

    DO di=1,ppm_dim
        IF (topo%bcdef(di) .EQ. ppm_param_bcdef_periodic) THEN
            DO ip=1,Npart
                IF (xp(di,ip) .GE. max_phys(di)) &
                    xp(di,ip) = xp(di,ip) - len_phys(di)*almostone
                IF (xp(di,ip) .LT. min_phys(di)) &
                    xp(di,ip) = xp(di,ip) + len_phys(di)*almostone
            ENDDO
        ELSE IF (topo%bcdef(di) .EQ. ppm_param_bcdef_freespace) THEN
            del_part = 0
            DO ip=Npart,1,-1
                IF (xp(di,ip).GE.max_phys(di).OR.xp(di,ip).LT.min_phys(di)) THEN
                    !delete particles that have crossed the boundary
                    IF (ip .EQ. Npart - del_part) THEN
                        ! no need to copy anything
                        ! this particle will be deleted when Npart is decreased
                    ELSE    
                        ! copying particles from the end of xp to the index that has
                        ! to be removed
                        xp(1:ppm_dim,ip) = xp(1:ppm_dim,Npart-del_part)
                        DO prop_id = 1,Particles%max_wpsid
                            IF (Particles%wps_m(prop_id) .EQ. 1) THEN
                                Particles%wps(prop_id)%vec(ip) = &
                                    Particles%wps(prop_id)%vec(Npart-del_part)
                            ENDIF
                        ENDDO
                        DO prop_id = 1,Particles%max_wpvid
                            IF (Particles%wpv_m(prop_id) .EQ. 1) THEN
                                Particles%wpv(prop_id)%vec(:,ip) = &
                                    Particles%wpv(prop_id)%vec(:,Npart-del_part)
                            ENDIF
                        ENDDO
                    ENDIF
                    del_part = del_part+1
                ENDIF
            ENDDO
            !New number of particles, after deleting some
            Npart = Npart - del_part
        ELSE
            CALL pwrite(ppm_rank,caller,&
                & 'this type of BC is not implemented/tested in this version',info)
            info=-1
            GOTO 9999
        ENDIF
    ENDDO

    ! Update states
    Particles%Npart = Npart
    ! Particles are no longer on the right processors
    Particles%ontopology = .FALSE.
    ! But they are now all inside the computational domain
    Particles%areinside = .TRUE.
    ! Dangereous to use the ghosts
    Particles%xp_g = 0
    ! ghosts values for properties are also dangerous to use
    DO prop_id = 1,Particles%max_wpsid
        IF (Particles%wps_g(prop_id) .EQ. 1) Particles%wps_g(prop_id) = 0
    ENDDO
    DO prop_id = 1,Particles%max_wpvid
        IF (Particles%wpv_g(prop_id) .EQ. 1) Particles%wpv_g(prop_id) = 0
    ENDDO

    IF (verbose) &
        write(*,*) 'applied periodic boundary conditions'
    !-----------------------------------------------------------------------
    ! Finalize
    !-----------------------------------------------------------------------
    xp => NULL()
    topo => NULL()
    CALL substop(caller,t0,info)

    9999 CONTINUE ! jump here upon error

END SUBROUTINE particles_apply_bc

SUBROUTINE particles_have_moved(Particles,info)
    !-----------------------------------------------------------------
    !  Update states when particles have moved 
    !   (assuming they have moved across processor boundaries)
    !-----------------------------------------------------------------
    USE ppm_module_data, ONLY: ppm_topo,ppm_rank
    USE ppm_module_substart
    USE ppm_module_substop
#if   __KIND == __SINGLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    TYPE(ppm_t_particles), POINTER,  INTENT(INOUT)         :: Particles
    !!! Data structure containing the particles
    INTEGER,                            INTENT(  OUT)      :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    INTEGER                                               :: di,prop_id
    REAL(ppm_kind_double)                                 :: t0
    CHARACTER(LEN = ppm_char)                 :: caller ='particles_have_moved'
    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)
    !-----------------------------------------------------------------
    !  checks
    !-----------------------------------------------------------------
    IF (.NOT.ASSOCIATED(Particles)) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    Particles%xp_g = 0

    DO prop_id = 1,Particles%max_wpsid
        IF(Particles%wps_g(prop_id) .EQ. 1) THEN
            Particles%wps_g(prop_id) = 0
        ENDIF
    ENDDO
    DO prop_id = 1,Particles%max_wpvid
        IF(Particles%wpv_g(prop_id) .EQ. 1) THEN
            Particles%wpv_g(prop_id) = 0
        ENDIF
    ENDDO

    Particles%ontopology = .FALSE.

    CALL substop(caller,t0,info)

    9999 CONTINUE ! jump here upon error

END SUBROUTINE particles_have_moved

SUBROUTINE particles_neighlists(Particles,topoid,info,lstore)
    !-----------------------------------------------------------------
    !  Neighbor lists for particles
    !-----------------------------------------------------------------
    !  Assumptions:
    ! * Particles positions need to have been mapped onto the topology
    ! * Ghost positions have been computed
    !
    USE ppm_module_substart
    USE ppm_module_substop
    USE ppm_module_neighlist
    USE ppm_module_inl_vlist
#if   __KIND == __SINGLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    TYPE(ppm_t_particles), POINTER,  INTENT(INOUT)         :: Particles
    !!! Data structure containing the particles
    INTEGER,                            INTENT(IN   )      :: topoid
    !!! topology id
    INTEGER,                            INTENT(  OUT)      :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    !  Optional arguments
    !-------------------------------------------------------------------------
    LOGICAL, OPTIONAL,                  INTENT(IN   )      :: lstore
    !!! store verlet lists
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    INTEGER , DIMENSION(3)                    :: ldc
    !!! Number of elements in all dimensions for allocation
    INTEGER                                   :: iopt
    !!! allocation mode, see one of ppm_alloc_* subroutines.
    INTEGER                                   :: prop_id
    !!! index variable
    LOGICAL                                   :: symmetry
    !!! backward compatibility
    REAL(MK),DIMENSION(2*ppm_dim):: ghostlayer
    !!!
    CHARACTER(LEN = ppm_char)                 :: caller = 'particles_neighlists'
    REAL(KIND(1.D0))                          :: t0
    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)
    !-----------------------------------------------------------------
    !  Checks
    !-----------------------------------------------------------------
    IF (.NOT.ASSOCIATED(Particles).OR..NOT.ASSOCIATED(Particles%xp)) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (.NOT.Particles%ontopology .OR. Particles%active_topoid.NE.topoid) THEN
        !Particles have not been mapped onto this topology
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'Do a partial/global mapping before',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (Particles%xp_g.NE.1) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'Ghosts have not been updated. They are needed for neighlists',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (Particles%isymm.EQ.1) THEN
        symmetry =.TRUE.
    ELSE
        symmetry = .FALSE.
    ENDIF

    IF ( Particles%neighlists ) THEN
        !neighbor lists are already up-to-date, nothing to do
        write(*,*) 'neighlists are supposedly already up-to-date, NOTHING to do'
    ELSE
        IF (Particles%adaptive) THEN
            !FIXME: when adaptive ghost layers are available
            ghostlayer(1:2*ppm_dim)=Particles%cutoff
            CALL ppm_inl_vlist(topoid,Particles%xp,Particles%Npart,&
                Particles%Mpart,Particles%wps(Particles%rcp_id)%vec,Particles%skin,&
                symmetry,ghostlayer,info,Particles%vlist,&
                Particles%nvlist)
            IF (info .NE. 0) THEN
                info = ppm_error_fatal
                CALL ppm_error(999,caller,&
                    'ppm_inl_vlist failed',__LINE__,info)
                GOTO 9999
            ENDIF
        ELSE
            CALL ppm_neighlist_vlist(topoid,Particles%xp,Particles%Mpart,&
                Particles%cutoff,Particles%skin,symmetry,Particles%vlist,&
                Particles%nvlist,info,lstore=lstore)
            IF (info .NE. 0) THEN
                info = ppm_error_fatal
                CALL ppm_error(999,caller,&
                    'ppm_neighlist_vlist failed',__LINE__,info)
                GOTO 9999
            ENDIF
        ENDIF
    ENDIF

    !Update state
    Particles%neighlists = .TRUE.
    !

    ! TODO:
    !WARNING: does not work with several processors!
    ! would require an MPI_GATHER
    ! Since this is mainly for debugging/warnings, I am not sure
    ! if it is really needed
    Particles%nneighmin = MINVAL(Particles%nvlist(1:Particles%Npart))
    Particles%nneighmax = MAXVAL(Particles%nvlist(1:Particles%Npart))

    IF (verbose) &
        write(*,*) 'computed neighlists'
    !-----------------------------------------------------------------------
    ! Finalize
    !-----------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999 CONTINUE ! jump here upon error

END SUBROUTINE particles_neighlists

SUBROUTINE particles_neighlists_xset(Particles_1,Particles_2,topoid,info,lstore)
    !------------------------------------------------------------------------
    !  Cross-set neighbor lists for particles
    !  Find the neighbours of each (real) Particles_1 within all Particles_2
    !------------------------------------------------------------------------
    !  Assumptions:
    ! * Particles (1 and 2) need to have been mapped onto the topology
    ! * Ghost positions for Particles_2 have been computed
    !
    USE ppm_module_substart
    USE ppm_module_substop
    USE ppm_module_inl_xset_vlist
#if   __KIND == __SINGLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    TYPE(ppm_t_particles), POINTER,  INTENT(INOUT)         :: Particles_1
    !!! Data structure containing the particles
    TYPE(ppm_t_particles), POINTER,  INTENT(IN   )         :: Particles_2
    !!! Data structure containing the particles
    INTEGER,                            INTENT(IN   )      :: topoid
    !!! topology id
    INTEGER,                            INTENT(  OUT)      :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    !  Optional arguments
    !-------------------------------------------------------------------------
    LOGICAL, OPTIONAL,                  INTENT(IN   )      :: lstore
    !!! store verlet lists
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    INTEGER , DIMENSION(3)                    :: ldc
    !!! Number of elements in all dimensions for allocation
    INTEGER                                   :: iopt
    !!! allocation mode, see one of ppm_alloc_* subroutines.
    INTEGER                                   :: prop_id
    !!! index variable
    LOGICAL                                   :: symmetry
    !!! backward compatibility
    REAL(MK),DIMENSION(2*ppm_dim):: ghostlayer
    !!!
    CHARACTER(LEN = ppm_char)          :: caller = 'particles_neighlists_xset'
    REAL(KIND(1.D0))                          :: t0
    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)
    !-----------------------------------------------------------------
    !  Checks
    !-----------------------------------------------------------------
    IF (.NOT.ASSOCIATED(Particles_1).OR..NOT.ASSOCIATED(Particles_2)) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
    IF (.NOT.ASSOCIATED(Particles_1%xp).OR..NOT.ASSOCIATED(Particles_2%xp)) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'Particles positions had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (.NOT.Particles_1%ontopology .OR. Particles_1%active_topoid.NE.topoid) THEN
        !Particles have not been mapped onto this topology
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'Do a partial/global mapping for Particles_1 before',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
    IF (.NOT.Particles_2%ontopology .OR. Particles_2%active_topoid.NE.topoid) THEN
        !Particles have not been mapped onto this topology
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'Do a partial/global mapping for Particles_2 before',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (Particles_2%xp_g.NE.1) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'Ghosts for Particles_2have not been updated.&
            & They are needed for neighlists',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF ( Particles_2%neighlists_cross ) THEN
        !neighbor lists are already up-to-date, nothing to do
        write(*,*) 'xset neighlists are supposedly already up-to-date, NOTHING to do'
    ELSE
        IF (Particles_2%adaptive) THEN
            !FIXME: when adaptive ghost layers are available
            ghostlayer(1:2*ppm_dim)=Particles_2%cutoff
            CALL ppm_inl_xset_vlist(topoid,Particles_1%xp,Particles_1%Npart,&
                Particles_1%Mpart,Particles_2%xp,Particles_2%Npart,&
                Particles_2%Mpart,Particles_2%wps(Particles_2%rcp_id)%vec,&
                Particles_2%skin,ghostlayer,info,Particles_1%vlist_cross,&
                Particles_1%nvlist_cross,lstore)
            IF (info .NE. 0) THEN
                info = ppm_error_fatal
                CALL ppm_error(999,caller,&
                    'ppm_inl_xset_vlist failed',__LINE__,info)
                GOTO 9999
            ENDIF
        ELSE
            ghostlayer(1:2*ppm_dim)=Particles_2%cutoff
            CALL ppm_inl_xset_vlist(topoid,Particles_1%xp,Particles_1%Npart,&
                Particles_1%Mpart,Particles_2%xp,Particles_2%Npart,&
                Particles_2%Mpart,Particles_2%cutoff,&
                Particles_2%skin,ghostlayer,info,Particles_1%vlist_cross,&
                Particles_1%nvlist_cross,lstore)
            IF (info .NE. 0) THEN
                info = ppm_error_fatal
                CALL ppm_error(999,caller,&
                    'ppm_neighlist_vlist failed',__LINE__,info)
                GOTO 9999
            ENDIF
        ENDIF
    ENDIF

    !Update state
    Particles_1%neighlists_cross = .TRUE.
    !

    ! TODO:
    !WARNING: does not work with several processors!
    ! would require an MPI_GATHER
    ! Since this is mainly for debugging/warnings, I am not sure
    ! if it is really needed
    Particles_1%nneighmin_cross = &
        MINVAL(Particles_1%nvlist_cross(1:Particles_1%Npart))
    Particles_1%nneighmax_cross = &
        MAXVAL(Particles_1%nvlist_cross(1:Particles_1%Npart))

    IF (verbose) &
        write(*,*) 'computed xset neighlists'
    !-----------------------------------------------------------------------
    ! Finalize
    !-----------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999 CONTINUE ! jump here upon error

END SUBROUTINE particles_neighlists_xset

SUBROUTINE particles_updated_positions(Particles,info)
    !-----------------------------------------------------------------
    !  routine to call when positions have changed
    !-----------------------------------------------------------------
    !  Effects:
    ! * discard ghosts
    ! * require
    ! * discard neighbour lists
    !  Assumptions:
    ! * Particles positions need to have been mapped onto the topology
    ! * cutoff radii are up-to-date
    !
    ! Note:
    ! * Does an MPI_Allreduce for to get the maximum cutoff
    USE ppm_module_data, ONLY: ppm_comm,ppm_mpi_kind
    USE ppm_module_substart
    USE ppm_module_substop
    INCLUDE "mpif.h"
#if   __KIND == __SINGLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    TYPE(ppm_t_particles), POINTER,  INTENT(INOUT)         :: Particles
    !!! Data structure containing the particles
    INTEGER,                            INTENT(  OUT)      :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    !!!
    CHARACTER(LEN = ppm_char)           :: caller='particles_updated_positions'
    REAL(MK)                              :: cutoff
    REAL(KIND(1.D0))                      :: t0
    INTEGER                               :: prop_id
    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)
    !-----------------------------------------------------------------
    !  Checks
    !-----------------------------------------------------------------
    IF (.NOT.ASSOCIATED(Particles).OR..NOT.ASSOCIATED(Particles%xp)) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    ! Update states
    ! particles may no longer be on the right processor
    Particles%ontopology=.FALSE.
    ! note that there is still a one-to-one relationship between real particles
    !  and indices. They can still be used for computing stuff or for
    ! being written to output (so we keep Particles%xp_m as it is)
    ! ghosts are now wrong (some particles may have entered the ghost layers)
     Particles%xp_g = 0
    ! and we have to throw away the neighbour lists (we dont destroy
    ! the pointers, but we flag the lists as being outdated)
    Particles%neighlists = .FALSE.
    ! ghosts values for properties are also wrong
    DO prop_id = 1,Particles%max_wpsid
        IF (Particles%wps_g(prop_id) .EQ. 1) Particles%wps_g(prop_id) = 0
    ENDDO
    DO prop_id = 1,Particles%max_wpvid
        IF (Particles%wpv_g(prop_id) .EQ. 1) Particles%wpv_g(prop_id) = 0
    ENDDO

    !-----------------------------------------------------------------------
    ! Finalize
    !-----------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999 CONTINUE ! jump here upon error

END SUBROUTINE particles_updated_positions

SUBROUTINE particles_updated_nb_part(Particles,info,&
        preserve_wps,preserve_wpv,preserve_neighlist)
    !-----------------------------------------------------------------
    !  routine to call when the number of particles has been
    !  manually changed
    !-----------------------------------------------------------------
    !  Effects:
    ! * discard ghosts
    ! * set properties to unmapped, unless they are specified
    !   in the ok_list (this means that the arrays for these
    !   properties have been reallocated to their new sizes
    !   and that a value has been assigned to each new particle)
    !  Assumptions:
    ! * The xp array has been allocated to the correct size
    !    and the positions of all particles are set properly
    !-----------------------------------------------------------------
    USE ppm_module_data, ONLY: ppm_comm,ppm_mpi_kind
    USE ppm_module_substart
    USE ppm_module_substop
    INCLUDE "mpif.h"
#if   __KIND == __SINGLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    TYPE(ppm_t_particles), POINTER,  INTENT(INOUT)         :: Particles
    !!! Data structure containing the particles
    INTEGER,                            INTENT(  OUT)      :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    !  Optional arguments
    !-------------------------------------------------------------------------
    INTEGER, DIMENSION(:), OPTIONAL,    INTENT(IN   )   :: preserve_wps
    !!! list of 1d properties that should not be discarded (default is empty)
    INTEGER, DIMENSION(:), OPTIONAL,    INTENT(IN   )   :: preserve_wpv
    !!! list of 2d properties that should not be discarded (default is empty)
    LOGICAL,               OPTIONAL,    INTENT(IN   )   :: preserve_neighlist
    !!! if true, do not discard neighbour lists (default is FALSE) 
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    !!!
    CHARACTER(LEN = ppm_char)           :: caller='particles_updated_nb_part'
    REAL(MK)                              :: cutoff
    REAL(KIND(1.D0))                      :: t0
    INTEGER                               :: prop_id,id
    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)
    !-----------------------------------------------------------------
    !  Checks
    !-----------------------------------------------------------------
    IF (.NOT.ASSOCIATED(Particles).OR..NOT.ASSOCIATED(Particles%xp)) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF


    ! Update states
    ! particles may no longer be on the right processor
    Particles%ontopology=.FALSE.
    ! Ghosts are now wrong (some particles may have been inserted in 
    !  the ghost layers)
     Particles%xp_g = 0
    ! Throw away the neighbour lists  (unless explicitely specified)
    IF (PRESENT(preserve_neighlist)) THEN
        IF (.NOT.(preserve_neighlist .AND. Particles%neighlists)) THEN
            Particles%neighlists = .FALSE.
        ENDIF
    ELSE
        Particles%neighlists = .FALSE.
    ENDIF
    Particles%neighlists = .FALSE.
    ! Mappings for some properties are now wrong
    ! 1D properties 
    !    increment wps_m for user-selected properties
    DO prop_id = 1,SIZE(preserve_wps)
        id = preserve_wps(prop_id)
        IF (id.GT.Particles%max_wpsid) THEN
            info = ppm_error_fatal
            CALL ppm_error(999,caller,   &
                &  'id in preserve_wps is larger than max id for Particles',&
                &  __LINE__,info)
            GOTO 9999
        ENDIF
        Particles%wps_m(id) = Particles%wps_m(id) + 1
    ENDDO
    !    decrement all wps_m (that had values 1 or 2)
    DO prop_id = 1,Particles%max_wpsid
        IF (Particles%wps_m(prop_id) .GE. 1) &
            Particles%wps_m(prop_id) =  Particles%wps_m(prop_id) - 1
    ENDDO

    ! 2D properties 
    !    increment wpv_m for user-selected properties
    DO prop_id = 1,SIZE(preserve_wpv)
        id = preserve_wpv(prop_id)
        IF (id.GT.Particles%max_wpvid) THEN
            info = ppm_error_fatal
            CALL ppm_error(999,caller,   &
                &  'id in preserve_wpv is larger than max id for Particles',&
                &  __LINE__,info)
            GOTO 9999
        ENDIF
        Particles%wpv_m(id) = Particles%wpv_m(id) + 1
    ENDDO
    !    decrement all wpv_m (that had values 1 or 2)
    DO prop_id = 1,Particles%max_wpvid
        IF (Particles%wpv_m(prop_id) .GE. 1) &
            Particles%wpv_m(prop_id) =  Particles%wpv_m(prop_id) - 1
    ENDDO

    ! Grow arrays with wp_m = 0 to size Npart, if needed
    DO prop_id = 1,Particles%max_wpsid
        IF (Particles%wps_m(prop_id) .EQ. 0) THEN
            id = prop_id
            CALL particles_allocate_wps(Particles,&
                id,info,&
                iopt=ppm_param_alloc_grow)
        ENDIF
    ENDDO
    DO prop_id = 1,Particles%max_wpvid
        IF (Particles%wpv_m(prop_id) .EQ. 0) THEN
            id = prop_id
            CALL particles_allocate_wpv(Particles,&
                id,Particles%wpv_s(id),info,&
                iopt=ppm_param_alloc_grow)
        ENDIF
    ENDDO

    ! ghosts values for properties are also wrong
    DO prop_id = 1,Particles%max_wpsid
        IF (Particles%wps_g(prop_id) .EQ. 1) Particles%wps_g(prop_id) = 0
    ENDDO
    DO prop_id = 1,Particles%max_wpvid
        IF (Particles%wpv_g(prop_id) .EQ. 1) Particles%wpv_g(prop_id) = 0
    ENDDO

    !-----------------------------------------------------------------------
    ! Finalize
    !-----------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999 CONTINUE ! jump here upon error

END SUBROUTINE particles_updated_nb_part

SUBROUTINE particles_updated_cutoff(Particles,info)
    !-----------------------------------------------------------------
    !  routine to call when cutoffs have changed
    !-----------------------------------------------------------------
    !  Effects:
    ! * update cutoff to maxval(rcp)
    ! * discard neighbour lists
    !  Assumptions:
    ! * Particles positions need to have been mapped onto the topology
    ! * cutoff radii are up-to-date
    !
    ! Note:
    ! * Does an MPI_Allreduce for to get the maximum cutoff
    USE ppm_module_data, ONLY: ppm_comm,ppm_mpi_kind
    USE ppm_module_substart
    USE ppm_module_substop
    INCLUDE "mpif.h"
#if   __KIND == __SINGLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    TYPE(ppm_t_particles), POINTER,  INTENT(INOUT)         :: Particles
    !!! Data structure containing the particles
    INTEGER,                            INTENT(  OUT)      :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    !!!
    CHARACTER(LEN = ppm_char)             :: caller = 'particles_updated_cutoff'
    REAL(MK)                              :: cutoff
    REAL(KIND(1.D0))                      :: t0
    INTEGER                               :: prop_id
    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)
    !-----------------------------------------------------------------
    !  Checks
    !-----------------------------------------------------------------
    IF (.NOT.ASSOCIATED(Particles).OR..NOT.ASSOCIATED(Particles%xp)) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (.NOT.Particles%adaptive .OR. Particles%rcp_id.EQ.0) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'Particles are not adaptive, cannot update cutoff',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    cutoff = MAXVAL(Particles%wps(Particles%rcp_id)%vec(1:Particles%Npart))

    CALL MPI_Allreduce(cutoff,cutoff,1,ppm_mpi_kind,MPI_MAX,ppm_comm,info)
    IF (info .NE. 0) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'MPI_Allreduce failed',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF


    ! Update states
    ! If cutoff has increased, then ghosts may no longer be 
    ! up to date
    IF (Particles%cutoff .LT. cutoff) THEN
        Particles%xp_g = 0
        !   ghosts values for some scalar arrays have been computed
        DO prop_id = 1,Particles%max_wpsid
            IF (Particles%wps_g(prop_id) .EQ. 1) Particles%wps_g(prop_id) = 0
        ENDDO
        !   ghosts values for some vector arrays have been computed
        DO prop_id = 1,Particles%max_wpvid
            IF (Particles%wpv_g(prop_id) .EQ. 1) Particles%wpv_g(prop_id) = 0
        ENDDO
    ENDIF
    ! neighbour lists have to be thrown away
    Particles%neighlists=.FALSE.
    !store new cutoff value
    Particles%cutoff = cutoff

    !-----------------------------------------------------------------------
    ! Finalize
    !-----------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999 CONTINUE ! jump here upon error

END SUBROUTINE particles_updated_cutoff

SUBROUTINE particles_compute_hmin(Particles,info)
    !-----------------------------------------------------------------
    !  compute minimum distance between particles
    !-----------------------------------------------------------------
    !  Effects:
    !
    ! Note:
    ! * Does an MPI_Allreduce to synchronize between processors
    !  (this could/should be modified such that the synchronization
    !   for such variables is done only when the variable is actually needed)
    !  (maybe set a flag Particles%synchronized=.FALSE./.TRUE. ? )
    USE ppm_module_data, ONLY: ppm_comm,ppm_mpi_kind
    USE ppm_module_substart
    USE ppm_module_substop
    INCLUDE "mpif.h"
#if   __KIND == __SINGLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    TYPE(ppm_t_particles), POINTER,  INTENT(INOUT)         :: Particles
    !!! Data structure containing the particles
    INTEGER,                            INTENT(  OUT)      :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    !!!
    CHARACTER(LEN = ppm_char)             :: caller = 'particles_compute_hmin'
    REAL(MK)                              :: hmin2
    REAL(KIND(1.D0))                      :: t0
    INTEGER                               :: ip,ineigh,iq
    INTEGER, DIMENSION(:),    POINTER     :: nvlist
    INTEGER, DIMENSION(:,:),  POINTER     :: vlist
    REAL(MK), DIMENSION(:,:), POINTER     :: xp
    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)
    !-----------------------------------------------------------------
    !  Checks
    !-----------------------------------------------------------------
    IF (.NOT.ASSOCIATED(Particles).OR..NOT.ASSOCIATED(Particles%xp)) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (.NOT.Particles%neighlists) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'Please construct neighbour lists first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
    IF (Particles%xp_g.NE.1) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'Please update ghost positions first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    nvlist => Particles%nvlist
    vlist => Particles%vlist
    xp => Particles%xp
    hmin2 = HUGE(1._MK)
    DO ip=1,Particles%Npart
        DO ineigh=1,nvlist(ip)
            iq=vlist(ineigh,ip)
            hmin2 = MIN(hmin2,SUM((xp(1:ppm_dim,ip)-xp(1:ppm_dim,iq))**2))
        ENDDO
    ENDDO
    xp=>NULL()
    nvlist=>NULL()
    vlist=>NULL()

    CALL MPI_Allreduce(hmin2,hmin2,1,ppm_mpi_kind,MPI_MIN,ppm_comm,info)
    IF (info .NE. 0) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'MPI_Allreduce failed',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF


    ! Update states
    Particles%h_min = SQRT(hmin2)

    !-----------------------------------------------------------------------
    ! Finalize
    !-----------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999 CONTINUE ! jump here upon error

END SUBROUTINE particles_compute_hmin

SUBROUTINE particles_initialise(Particles,Npart_global,info,&
        distrib,topoid,minphys,maxphys)
    !-----------------------------------------------------------------------
    ! Set initial particle positions
    !-----------------------------------------------------------------------
    USE ppm_module_substart
    USE ppm_module_substop
    USE ppm_module_data, ONLY: ppm_rank,ppm_nproc,ppm_topo
#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    TYPE(ppm_t_particles), POINTER,     INTENT(INOUT)      :: Particles
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
    !!! ppm_param_part_init_cartesian
    !!! ppm_param_part_init_random
    INTEGER,OPTIONAL,                   INTENT(IN   )      :: topoid
    !!! topology id (used only to get the extent of the physical domain)
    REAL(MK),DIMENSION(ppm_dim),OPTIONAL,INTENT(IN   )     :: minphys
    !!! extent of the physical domain. Only if topoid is not present.
    REAL(MK),DIMENSION(ppm_dim),OPTIONAL,INTENT(IN   )     :: maxphys
    !!! extent of the physical domain. Only if topoid is not present.
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------

    INTEGER                               :: ip,i,j,k,Npart
    INTEGER                               :: nijk(ppm_dim),nijk_global(ppm_dim)
    CHARACTER(LEN = ppm_char)              :: filename,cbuf
    CHARACTER(LEN = ppm_char)              :: caller = 'particles_initialise'
    REAL(MK)                              :: y,z,h
    REAL(KIND(1.D0))                      :: t0
    INTEGER                               :: remaining_rows

    REAL(MK)                              :: shift
    INTEGER                               :: distribution
    TYPE(ppm_t_topo),POINTER              :: topo => NULL()
    REAL(MK), DIMENSION(ppm_dim)          :: min_phys,max_phys,len_phys

    REAL(MK), DIMENSION(:,:), POINTER     :: xp
    REAL(MK), DIMENSION(:  ), ALLOCATABLE :: randnb
    INTEGER,  DIMENSION(:  ), ALLOCATABLE :: seed
    INTEGER                               :: seedsize


    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)

    IF(PRESENT(distrib)) THEN
        distribution=distrib
    ELSE
        distribution=ppm_param_part_init_cartesian
    ENDIF

    !Get boundaries of computational domain
    IF (PRESENT(topoid) .AND. (PRESENT(minphys).OR.PRESENT(maxphys))) THEN
            info = ppm_error_fatal
            CALL ppm_error(999,caller,&
               'probable conflict of optional arguments. Use topoid OR minphys'&
               ,__LINE__,info)
            GOTO 9999
    ENDIF
    IF (PRESENT(topoid)) THEN
        topo => ppm_topo(topoid)%t
        min_phys = topo%min_physd
        max_phys = topo%max_physd
    ELSE IF (PRESENT(minphys).AND.PRESENT(maxphys)) THEN
        min_phys = minphys
        max_phys = maxphys
    ELSE
        info = ppm_error_fatal
        CALL ppm_error(999,caller,&
            'optional arguments needed to define the domain boundaries'&
            ,__LINE__,info)
        GOTO 9999
    ENDIF
    len_phys=max_phys-min_phys


    h = (PRODUCT(len_phys)/REAL(Npart_global))**(1./REAL(ppm_dim))
    nijk_global = FLOOR(len_phys/h)
    Npart_global = PRODUCT(nijk_global)
    remaining_rows = MOD(nijk_global(ppm_dim),ppm_nproc)

    !number of particles along x 
    nijk(1:ppm_dim-1) = nijk_global(1:ppm_dim-1)
    !number of particles along y 
    nijk(ppm_dim) = nijk_global(ppm_dim)/ppm_nproc

    !number of particles on this processor
    Npart = PRODUCT(nijk)
    !proc 0 takes care of the additional rows (remainder)
    IF (ppm_rank.EQ.0) THEN
        Npart = Npart + remaining_rows * PRODUCT(nijk(1:ppm_dim-1))
    ENDIF

    !Deallocate Particles if already allocated
    IF (ASSOCIATED(Particles)) THEN
        CALL ppm_alloc_particles(Particles,Npart,ppm_param_dealloc,info)
        IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_dealloc,caller,&
                'ppm_alloc_particles (deallocate) failed',__LINE__,info)
            GOTO 9999
        ENDIF
    ENDIF

    CALL ppm_alloc_particles(Particles,Npart,ppm_param_alloc_fit,info)
    IF (info .NE. 0) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_alloc,caller,&
            'ppm_alloc_particles (allocate) failed',__LINE__,info)
        GOTO 9999
    ENDIF
    !use a shortcut, for convenience
    xp => Particles%xp

#ifdef same_random_sequence_nproc
    ALLOCATE(randnb(ppm_dim*Npart_global))
#else if
    ALLOCATE(randnb(ppm_dim*Npart))
#endif
    CALL RANDOM_SEED(SIZE=seedsize)
    ALLOCATE(seed(seedsize))
    DO i=1,seedsize
        seed(i)=i*i*i*i
    ENDDO
    CALL RANDOM_SEED(PUT=seed)
    CALL RANDOM_NUMBER(randnb)

    !-----------------------------------------------------------------------
    ! set particles
    !-----------------------------------------------------------------------
    ip = 0
    shift = 0._MK !shifts positions of the particles. Set to 0.5 to place
    !particles in the middle of each cell, set to 0 to place them in the lower
    !left corner

    if_cartesian: IF (distribution .EQ. ppm_param_part_init_cartesian) THEN
#ifdef __3D
        DO k = 1,nijk(3)
            h = len_phys(3)/REAL(nijk(3),MK)
            z = min_phys(3) + h*(k-1) + shift*h
#endif
        DO j = 1,nijk(2)
            h = len_phys(2)/REAL(nijk_global(2),MK)
            y = min_phys(2) + h*(j-1 + ppm_rank*nijk(2)) + shift*h
            DO i = 1,nijk(1)
                h = len_phys(1)/REAL(nijk(1),MK)

                ip = ip + 1
                xp(1,ip) = min_phys(1) + h*(i-1) + shift*h
                xp(2,ip) = y                  
#ifdef __3D
                xp(3,ip) = z 
#endif

                ! impose periodic boundaries:
                IF (xp(1,ip) .GE. max_phys(1)) xp(1,ip) = xp(1,ip) - len_phys(1)
                IF (xp(2,ip) .GE. max_phys(2)) xp(2,ip) = xp(2,ip) - len_phys(2)
                IF (xp(1,ip) .LT. min_phys(1)) xp(1,ip) = xp(1,ip) + len_phys(1)
                IF (xp(2,ip) .LT. min_phys(2)) xp(2,ip) = xp(2,ip) + len_phys(2)
#ifdef __3D
                IF (xp(3,ip) .GE. max_phys(3)) xp(3,ip) = xp(3,ip) - len_phys(3)
                IF (xp(3,ip) .LT. min_phys(3)) xp(3,ip) = xp(3,ip) + len_phys(3)
#endif

            ENDDO
        ENDDO
#ifdef __3D
        ENDDO
#endif

        IF(ppm_rank.EQ.0) THEN
#ifdef __3D
        DO k = 1,nijk(3)
            h = len_phys(3)/REAL(nijk(3),MK)
            z = min_phys(3) + h*(k-1) + shift*h
#endif
            DO j = 1,remaining_rows
                h = len_phys(2)/REAL(nijk_global(2),MK)
                y = min_phys(2) + h*(j-1 + ppm_nproc*nijk(2)) + shift*h
                DO i = 1,nijk(1)
                    h = len_phys(1)/REAL(nijk(1),MK)

                    ip = ip + 1
                    xp(1,ip) = min_phys(1) + h*(i-1) + shift*h
                    xp(2,ip) = y                  
#ifdef __3D
                    xp(3,ip) = z 
#endif

                    ! impose periodic boundaries:
                    IF (xp(1,ip) .GE. max_phys(1)) xp(1,ip) = xp(1,ip) - len_phys(1)
                    IF (xp(2,ip) .GE. max_phys(2)) xp(2,ip) = xp(2,ip) - len_phys(2)
                    IF (xp(1,ip) .LT. min_phys(1)) xp(1,ip) = xp(1,ip) + len_phys(1)
                    IF (xp(2,ip) .LT. min_phys(2)) xp(2,ip) = xp(2,ip) + len_phys(2)
#ifdef __3D
                    IF (xp(3,ip) .GE. max_phys(3)) xp(3,ip) = xp(3,ip) - len_phys(3)
                    IF (xp(3,ip) .LT. min_phys(3)) xp(3,ip) = xp(3,ip) + len_phys(3)
#endif

                ENDDO
            ENDDO
#ifdef __3D
        ENDDO
#endif
        ENDIF
        Particles%cartesian = .TRUE.
    ELSE !random distribution
#ifdef __3D
        DO k = 1,nijk(3)
            h = len_phys(3)/REAL(nijk(3),MK)
            z = min_phys(3) + h*(k-1) + shift*h
#endif
        DO j = 1,nijk(2)
            h = len_phys(2)/REAL(nijk_global(2),MK)
            y = min_phys(2) + h*(j-1 + ppm_rank*nijk(2))  + shift* h

            DO i = 1,nijk(1)
                h = len_phys(1)/REAL(nijk(1),MK)

                ip = ip + 1
                ! uniformly random in cells
#ifdef same_random_sequence_nproc
                xp(1,ip) = min_phys(1) + h*(i-1) + shift + &
                    randnb(ppm_dim*ppm_rank*PRODUCT(nijk)+ppm_dim*ip - 1)*h
                xp(2,ip) = y                   + &
                    randnb(ppm_dim*ppm_rank*PRODUCT(nijk)+ppm_dim*ip    )*h
#ifdef __3D
                xp(3,ip) = z                   + &
                    randnb(ppm_dim*ppm_rank*PRODUCT(nijk)+ppm_dim*ip - 2)*h
#endif
#else if
                xp(1,ip) = min_phys(1) + h*(i-1) + shift + &
                    randnb(ppm_dim*ip - 1)*h
                xp(2,ip) = y                   + &
                    randnb(ppm_dim*ip    )*h
#ifdef __3D
                xp(3,ip) = z                   + &
                    randnb(ppm_dim*ip - 2)*h
#endif
#endif
                ! impose periodic boundaries:
                IF (xp(1,ip) .GE. max_phys(1)) xp(1,ip) = xp(1,ip) - len_phys(1)
                IF (xp(2,ip) .GE. max_phys(2)) xp(2,ip) = xp(2,ip) - len_phys(2)
                IF (xp(1,ip) .LT. min_phys(1)) xp(1,ip) = xp(1,ip) + len_phys(1)
                IF (xp(2,ip) .LT. min_phys(2)) xp(2,ip) = xp(2,ip) + len_phys(2)
#ifdef __3D
                IF (xp(3,ip) .GE. max_phys(3)) xp(3,ip) = xp(3,ip) - len_phys(3)
                IF (xp(3,ip) .LT. min_phys(3)) xp(3,ip) = xp(3,ip) + len_phys(3)
#endif
            ENDDO
        ENDDO
#ifdef __3D
        ENDDO
#endif
        IF(ppm_rank.EQ.0) THEN
#ifdef __3D
        DO k = 1,nijk(3)
            h = len_phys(3)/REAL(nijk(3),MK)
            z = min_phys(3) + h*(k-1) + shift*h
#endif
            DO j = 1,remaining_rows
                h = len_phys(2)/REAL(nijk_global(2),MK)
                y = min_phys(2) + h*(j-1 + ppm_nproc*nijk(2)) + shift*h
                DO i = 1,nijk(1)
                    h = len_phys(1)/REAL(nijk(1),MK)

                    ip = ip + 1
                    ! uniformly random in cells
#ifdef same_random_sequence_nproc
                    xp(1,ip) = min_phys(1) + h*(i-1) + shift + &
                        randnb(ppm_dim*(ppm_nproc-1)*&
                        PRODUCT(nijk)+ppm_dim*ip - 1)*h
                    xp(2,ip) = y                   + &
                        randnb(ppm_dim*(ppm_nproc-1)*&
                        PRODUCT(nijk)+ppm_dim*ip    )*h
#ifdef __3D
                    xp(3,ip) = z                   + &
                        randnb(ppm_dim*(ppm_nproc-1)*&
                        PRODUCT(nijk)+ppm_dim*ip - 2)*h
#endif
#else if
                    xp(1,ip) = min_phys(1) + h*(i-1) + shift + &
                        randnb(ppm_dim*ip - 1)*h
                    xp(2,ip) = y                   + &
                        randnb(ppm_dim*ip    )*h
#ifdef __3D
                    xp(3,ip) = z                   + &
                        randnb(ppm_dim*ip - 2)*h
#endif
#endif
                    ! impose periodic boundaries:
                    IF (xp(1,ip) .GE. max_phys(1)) xp(1,ip) = xp(1,ip) - len_phys(1)
                    IF (xp(2,ip) .GE. max_phys(2)) xp(2,ip) = xp(2,ip) - len_phys(2)
                    IF (xp(1,ip) .LT. min_phys(1)) xp(1,ip) = xp(1,ip) + len_phys(1)
                    IF (xp(2,ip) .LT. min_phys(2)) xp(2,ip) = xp(2,ip) + len_phys(2)
#ifdef __3D
                    IF (xp(3,ip) .GE. max_phys(3)) xp(3,ip) = xp(3,ip) - len_phys(3)
                    IF (xp(3,ip) .LT. min_phys(3)) xp(3,ip) = xp(3,ip) + len_phys(3)
#endif
                ENDDO
            ENDDO
#ifdef __3D
        ENDDO
#endif
        ENDIF
        Particles%cartesian = .FALSE.
    ENDIF if_cartesian

    xp=>NULL()
    DEALLOCATE(randnb,seed)

    ! (global) average interparticle spacing
    Particles%h_avg = (PRODUCT(len_phys)/REAL(Npart_global))**(1./REAL(ppm_dim)) 
    ! min interparticle spacing (not needed now)
    Particles%h_min = -1._MK

    Particles%areinside = .TRUE.
    ! neighbour lists not updated
    Particles%neighlists = .FALSE.

    !-----------------------------------------------------------------------
    ! Finalize
    !-----------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999 CONTINUE ! jump here upon error

END SUBROUTINE particles_initialise


SUBROUTINE particles_io_xyz(Particles,itnum,writedir,info,&
        with_ghosts,wps_list,wpv_list)
    !!!------------------------------------------------------------------------!
    !!! Dump data files to disk
    !!! All data is gathered to rank 0 and then written as one single file.
    !!!
    !!! NOT to be used with large number of processors!...
    !!!------------------------------------------------------------------------!
    USE ppm_module_substart
    USE ppm_module_substop
    USE ppm_module_data, ONLY: ppm_rank,ppm_nproc,ppm_comm

    INCLUDE "mpif.h"

#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    TYPE(ppm_t_particles), POINTER,     INTENT(IN   )   :: Particles
    !!! Data structure containing the particles
    INTEGER,                            INTENT(IN   )   :: itnum
    !!! iteration number
    CHARACTER(LEN=ppm_char),            INTENT(IN   )   :: writedir
    !!! directory name
    INTEGER,                            INTENT(  OUT)   :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    !  Optional arguments
    !-------------------------------------------------------------------------
    LOGICAL, OPTIONAL,                  INTENT(IN   )   :: with_ghosts
    !!! also print ghost particles (and the values of their properties)
    INTEGER, DIMENSION(:), OPTIONAL,    INTENT(IN   )   :: wps_list
    !!! list of 1d properties that will be written to file
    INTEGER, DIMENSION(:), OPTIONAL,    INTENT(IN   )   :: wpv_list
    !!! list of 2d properties that will be written to file
    !-------------------------------------------------------------------------
    ! local variables
    !-------------------------------------------------------------------------
    CHARACTER(LEN = ppm_char)                  :: cbuf,filename
    CHARACTER(LEN = ppm_char)                  :: caller = 'particles_io_xyz'
    CHARACTER(LEN=128)                         :: myformat
    INTEGER,DIMENSION(:),ALLOCATABLE           :: Npart_vec
    REAL(MK),DIMENSION(:,:),ALLOCATABLE        :: propp_iproc
    REAL(MK),DIMENSION(:,:),ALLOCATABLE        :: xp_iproc
    INTEGER                                    :: count,iproc,i,j,ip
    INTEGER,DIMENSION(MPI_STATUS_SIZE)         :: status
    REAL(KIND(1.D0))                           :: t0
    INTEGER                                    :: ndim2,wpv_s,Npart_local
    LOGICAL                                    :: ghosts
    INTEGER                                    :: nb_wps,nb_wpv
    INTEGER,DIMENSION(:),ALLOCATABLE           :: wps_l,wpv_l


    !-------------------------------------------------------------------------
    ! Initialize
    !-------------------------------------------------------------------------
    info = 0 ! change if error occurs
    CALL substart(caller,t0,info)

    IF(PRESENT(with_ghosts)) THEN
        ghosts=with_ghosts
    ELSE
        ghosts=.FALSE.
    ENDIF

    !-------------------------------------------------------------------------
    ! Checks
    !-------------------------------------------------------------------------
    IF (.NOT.ASSOCIATED(Particles).OR..NOT.ASSOCIATED(Particles%xp)) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    !-------------------------------------------------------------------------
    ! List of properties (1d and 2d) to writeout (default is
    ! all those that are mapped with the particles, i.e. for which
    ! wps_m = 1 (resp. wpv_m = 1)
    !-------------------------------------------------------------------------
    IF(PRESENT(wps_list)) THEN
        nb_wps=SIZE(wps_list)
        ALLOCATE(wps_l(nb_wps),STAT=info)
        wps_l=wps_list
        DO i=1,nb_wps
            IF (wps_l(i).GT.Particles%max_wpsid) THEN
                info = ppm_error_fatal
                CALL ppm_error(999,caller,   &
                    &  'property index exceeds size of property array',&
                    &  __LINE__,info)
                GOTO 9999
            ENDIF
            IF (Particles%wps_m(wps_l(i)).NE.1) THEN
                info = ppm_error_fatal
                CALL ppm_error(999,caller,   &
                    &  'trying to printout a property that is not mapped &
                    & to the particles',&
                    &  __LINE__,info)
                GOTO 9999
            ENDIF
        ENDDO
    ELSE
        !printout all properties i for which wps_m(i) = 1
        nb_wps = 0
        DO i=1,Particles%max_wpsid
            IF (Particles%wps_m(i).EQ.1) &
                nb_wps = nb_wps + 1
        ENDDO
        ALLOCATE(wps_l(nb_wps),STAT=info)
        nb_wps = 0
        DO i=1,Particles%max_wpsid
            IF (Particles%wps_m(i).EQ.1) &
                nb_wps = nb_wps + 1
                wps_l(nb_wps) = i
        ENDDO
    ENDIF
    IF(PRESENT(wpv_list)) THEN
        nb_wpv=SIZE(wpv_list)
        ALLOCATE(wpv_l(nb_wpv),STAT=info)
        wpv_l=wpv_list
        DO i=1,nb_wpv
            IF (wpv_l(i).GT.Particles%max_wpvid) THEN
                info = ppm_error_fatal
                CALL ppm_error(999,caller,   &
                    &  'property index exceeds size of property array',&
                    &  __LINE__,info)
                GOTO 9999
            ENDIF
            IF (Particles%wpv_m(wpv_l(i)).NE.1) THEN
                info = ppm_error_fatal
                CALL ppm_error(999,caller,   &
                    &  'trying to printout a property that is not mapped &
                    & to the particles',&
                    &  __LINE__,info)
                GOTO 9999
            ENDIF
        ENDDO
    ELSE
        !printout all properties i for which wpv_m(i) = 1
        nb_wpv = 0
        DO i=1,Particles%max_wpvid
            IF (Particles%wpv_m(i).EQ.1) &
                nb_wpv = nb_wpv + 1
        ENDDO
        ALLOCATE(wpv_l(nb_wpv),STAT=info)
        nb_wpv = 0
        DO i=1,Particles%max_wpvid
            IF (Particles%wpv_m(i).EQ.1) &
                nb_wpv = nb_wpv + 1
                wpv_l(nb_wpv) = i
        ENDDO
    ENDIF

    !count number of properties to be printed
    ndim2 = nb_wps
    DO i=1,nb_wpv
        ndim2 = ndim2 + Particles%wpv_s(wpv_l(i))
    ENDDO
    ! add one column that gives the rank for ghost particles
    ! and -1 for real particles
    IF (ghosts) ndim2 = ndim2 + 1

#ifdef highprec_writeout
    WRITE(myformat,'(A,I2,A)') '( ',ppm_dim+ndim2,'E24.13E3)'
#else
    WRITE(myformat,'(A,I2,A)') '( ',ppm_dim+ndim2,'E17.7E3)'
#endif

    !====================================================================!
    ! open file
    WRITE(filename,'(A,A,I6.6,A)') TRIM(writedir),'Particles_',itnum,'.xyz'
    IF (ppm_rank .EQ. 0) THEN
        OPEN(23,FILE=filename,FORM='FORMATTED',STATUS='REPLACE',IOSTAT=info)
        IF (info .NE. 0) THEN
            CALL pwrite(ppm_rank,caller,'opening file failed.',info)
            info = -1
            GOTO 9999
        ENDIF
    ENDIF

    !====================================================================!
    ! rank 0 collects data and writes out sequentially
    !FIXME these lines should not need to be commented
    !had to do this otherwise ifort complains at run time
    !when all the debug flags are turned on...
    ALLOCATE(Npart_vec(ppm_nproc),STAT=info)
    IF (info .NE. 0) THEN
        CALL pwrite(ppm_rank,caller,'allocation failed.',info)
        info = -1
        GOTO 9999
    ENDIF

    IF(ghosts) THEN
        Npart_local = Particles%Mpart
    ELSE
        Npart_local = Particles%Npart
    ENDIF

    CALL MPI_Gather(Npart_local,1,MPI_INTEGER,Npart_vec,1,&
        MPI_INTEGER,0,ppm_comm,info)
    IF (info .NE. 0) THEN
        CALL pwrite(ppm_rank,caller,'MPI_Gather failed.',info)
        info = -1
        GOTO 9999
    ENDIF

    IF (ppm_rank .EQ. 0) THEN

        ip = 0
        DO iproc = 1,ppm_nproc

            IF (iproc .GT. 1) THEN

                ALLOCATE(xp_iproc(ppm_dim,Npart_vec(iproc)),STAT=info)
                IF (info .NE. 0) THEN
                    CALL pwrite(ppm_rank,caller,'allocation failed.',info)
                    info = -1
                    GOTO 9999
                ENDIF
                ALLOCATE(propp_iproc(ndim2,Npart_vec(iproc)),STAT=info)
                IF (info .NE. 0) THEN
                    CALL pwrite(ppm_rank,caller,'allocation failed.',info)
                    info = -1
                    GOTO 9999
                ENDIF
                IF (MK .EQ. 4) THEN
                    count = ppm_dim*Npart_vec(iproc)
                    CALL MPI_Recv(xp_iproc,count,MPI_REAL,&
                        iproc-1,iproc-1,ppm_comm,status,info)
                    count = ndim2*Npart_vec(iproc)
                    CALL MPI_Recv(propp_iproc,count,MPI_REAL,&
                        iproc-1,iproc-1,ppm_comm,status,info)
                ELSE
                    count = ppm_dim*Npart_vec(iproc)
                    CALL MPI_Recv(xp_iproc,count,MPI_DOUBLE_PRECISION,  &
                        iproc-1,iproc-1,ppm_comm,status,info)
                    count = ndim2*Npart_vec(iproc)
                    CALL MPI_Recv(propp_iproc,count,MPI_DOUBLE_PRECISION,  &
                        iproc-1,iproc-1,ppm_comm,status,info)
                ENDIF
                DO i=1,Npart_vec(iproc)
                    ip = ip + 1
                    WRITE (23, myformat) xp_iproc(1:ppm_dim,i),&
                        (propp_iproc(j,i), j=1,ndim2) 
                END DO

                DEALLOCATE(propp_iproc)
                DEALLOCATE(xp_iproc)

            ELSE

                ALLOCATE(propp_iproc(ndim2,Npart_local),STAT=info)
                IF (info .NE. 0) THEN
                    CALL pwrite(ppm_rank,caller,'allocation failed.',info)
                    info = -1
                    GOTO 9999
                ENDIF

                ndim2= 0
                DO i=1,nb_wps
                    DO j=1,Npart_local
                        propp_iproc(i,j)=Particles%wps(wps_l(i))%vec(j)
                    ENDDO
                    ndim2 = ndim2 + 1
                ENDDO
                DO i=1,nb_wpv
                    wpv_s=Particles%wpv_s(wpv_l(i))
                    DO j=1,Npart_local
                        propp_iproc(ndim2+1:ndim2+wpv_s,j)=&
                            Particles%wpv(wpv_l(i))%vec(1:wpv_s,j)
                    ENDDO
                    ndim2=ndim2+wpv_s
                ENDDO
                IF (ghosts) THEN
                    DO j=1,Particles%Npart
                        propp_iproc(ndim2+1,j)= -1._MK
                    ENDDO
                    DO j=Particles%Npart+1,Particles%Mpart
                        propp_iproc(ndim2+1,j)= REAL(ppm_rank,MK)
                    ENDDO
                    ndim2=ndim2+1
                ENDIF

                DO i=1,Npart_vec(iproc)
                    ip = ip + 1
                    WRITE (23, myformat) Particles%xp(1:ppm_dim,i),&
                        (propp_iproc(j,i), j=1,ndim2) 
                END DO

                DEALLOCATE(propp_iproc)
            ENDIF

        ENDDO

        !==================================================================!
        ! close file

        CLOSE(23,IOSTAT=info)
        IF (info .NE. 0) THEN
            CALL pwrite(ppm_rank,caller,'closing file failed.',info)
            info = -1
            GOTO 9999
        ENDIF

        DEALLOCATE(Npart_vec)

    ELSE

        count = ppm_dim*Npart_local
        IF (MK .EQ. 4) THEN
            CALL MPI_Send(Particles%xp,count,MPI_REAL,0,ppm_rank,ppm_comm,info)
            IF (info .NE. 0) THEN
                CALL pwrite(ppm_rank,caller,'MPI_Send failed.',info)
                info = -1
                GOTO 9999
            ENDIF
        ELSE
            CALL MPI_Send(Particles%xp,count,MPI_DOUBLE_PRECISION,0,&
                ppm_rank,ppm_comm,info)
            IF (info .NE. 0) THEN
                CALL pwrite(ppm_rank,caller,'MPI_Send failed.',info)
                info = -1
                GOTO 9999
            ENDIF
        ENDIF

        count = ndim2*Npart_local
        ALLOCATE(propp_iproc(ndim2,Npart_local),STAT=info)
        IF (info .NE. 0) THEN
            CALL pwrite(ppm_rank,caller,'allocation failed.',info)
            info = -1
            GOTO 9999
        ENDIF

        ndim2 = 0
        DO i=1,nb_wps
            DO j=1,Npart_local
                propp_iproc(i,j)=Particles%wps(wps_l(i))%vec(j)
            ENDDO
            ndim2 = ndim2 + 1
        ENDDO
        DO i=1,nb_wpv
            wpv_s=Particles%wpv_s(wpv_l(i))
            DO j=1,Npart_local
                propp_iproc(ndim2+1:ndim2+wpv_s,j) = &
                    Particles%wpv(wpv_l(i))%vec(1:wpv_s,j)
            ENDDO
            ndim2=ndim2+wpv_s
        ENDDO
        IF (ghosts) THEN
            DO j=1,Particles%Npart
                propp_iproc(ndim2+1,j)= -1._MK
            ENDDO
            DO j=Particles%Npart+1,Particles%Mpart
                propp_iproc(ndim2+1,j)= REAL(ppm_rank,MK)
            ENDDO
        ENDIF
        IF (MK .EQ. 4) THEN
            CALL MPI_Send(propp_iproc,count,MPI_REAL,0,ppm_rank,ppm_comm,info)
            IF (info .NE. 0) THEN
                CALL pwrite(ppm_rank,caller,'MPI_Send failed.',info)
                info = -1
                GOTO 9999
            ENDIF
        ELSE
            CALL MPI_Send(propp_iproc,count,MPI_DOUBLE_PRECISION,0,&
                ppm_rank,ppm_comm,info)
            IF (info .NE. 0) THEN
                CALL pwrite(ppm_rank,caller,'MPI_Send failed.',info)
                info = -1
                GOTO 9999
            ENDIF
        ENDIF

        DEALLOCATE(propp_iproc)
    ENDIF

    call MPI_barrier(ppm_comm,info)
    If (ppm_rank.EQ.0) write(*,*) 'Wrote ',ip,' particles to file'

    !-------------------------------------------------------------------------
    !  Finalize
    !-------------------------------------------------------------------------
    DEALLOCATE(wps_l,wpv_l)
    CALL substop(caller,t0,info)

    9999  CONTINUE ! jump here upon error

END SUBROUTINE particles_io_xyz

SUBROUTINE particles_apply_dcops(Particles,from_id,to_id,eta_id,sig,&
        info,Particles_old,from_old_id)
    !!!------------------------------------------------------------------------!
    !!! Apply DC kernel stored in eta_id to the scalar property stored
    !!! prop_from_id and store the results in prop_to_id
    !!!------------------------------------------------------------------------!
    USE ppm_module_substart
    USE ppm_module_substop
    USE ppm_module_data, ONLY: ppm_rank

    INCLUDE "mpif.h"

#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    TYPE(ppm_t_particles), POINTER,     INTENT(INOUT)   :: Particles
    !!! Data structure containing the particles
    INTEGER,                            INTENT(IN   )   :: from_id
    !!! id where the data is stored
    INTEGER,                            INTENT(INOUT)   :: to_id
    !!! id where the result should be stored (0 if it needs to be allocated)
    INTEGER,                            INTENT(IN   )   :: eta_id
    !!! id where the DC kernel has been stored
    INTEGER,                            INTENT(IN   )   :: sig
    !!! value of the plus/minus sign in the formula for DC ops
    !!! valid choices are 1,-1 or 0
    INTEGER,                            INTENT(  OUT)   :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    ! Optional argument
    !-------------------------------------------------------------------------
    TYPE(ppm_t_particles), POINTER,OPTIONAL,INTENT(IN   )   :: Particles_old
    !!! (for interpolation kernels) Data structure containing another set
    !!! of particles particles
    INTEGER,OPTIONAL,                   INTENT(IN   )   :: from_old_id
    !!! id where the data is stored
    !-------------------------------------------------------------------------
    ! local variables
    !-------------------------------------------------------------------------
    CHARACTER(LEN = ppm_char)                  :: cbuf,filename
    CHARACTER(LEN = ppm_char)               :: caller = 'particles_apply_dcops'
    INTEGER                                    :: ip,iq,ineigh
    REAL(KIND(1.D0))                           :: t0
    REAL(MK),DIMENSION(:,:),POINTER            :: xp1 => NULL()
    REAL(MK),DIMENSION(:,:),POINTER            :: xp2 => NULL()
    REAL(MK),DIMENSION(:,:),POINTER            :: eta => NULL()
    REAL(MK),DIMENSION(:),  POINTER            :: wp1 => NULL()
    REAL(MK),DIMENSION(:),  POINTER            :: wp2 => NULL()
    REAL(MK),DIMENSION(:),  POINTER            :: dwp => NULL()
    INTEGER, DIMENSION(:),  POINTER            :: nvlist => NULL()
    INTEGER, DIMENSION(:,:),POINTER            :: vlist => NULL()


    !-------------------------------------------------------------------------
    ! Initialize
    !-------------------------------------------------------------------------
    info = 0 ! change if error occurs
    CALL substart(caller,t0,info)

    IF (to_id.EQ.0) THEN
        CALL particles_allocate_wps(Particles,to_id,info)
        IF (info .NE. 0) THEN
            CALL pwrite(ppm_rank,caller,'particles_allocate_wps failed.',info)
            info = -1
            GOTO 9999
        ENDIF
    ENDIF


    wp1 => Get_wps(Particles,from_id)
    xp1 => Get_xp(Particles)
    IF (PRESENT(Particles_old)) THEN
        wp2 => Get_wps(Particles_old,from_old_id,with_ghosts=.TRUE.)
        xp2 => Get_xp(Particles,with_ghosts=.TRUE.)
        nvlist => Particles%nvlist_cross
        vlist => Particles%vlist_cross
    ELSE
        wp2 => wp1
        xp2 => xp1
        nvlist => Particles%nvlist
        vlist => Particles%vlist
    ENDIF

    dwp => Get_wps(Particles,to_id)
    eta => Get_wpv(Particles,eta_id)

    DO ip = 1,Particles%Npart 
        dwp(ip) = 0._MK
    ENDDO
    DO ip = 1,Particles%Npart
        DO ineigh = 1,nvlist(ip)
            iq = vlist(ineigh,ip)
            dwp(ip) = dwp(ip) + (wp2(iq)+sig*(wp1(ip))) * eta(ineigh,ip)
        ENDDO
    ENDDO

    eta => Set_wpv(Particles,eta_id,read_only=.TRUE.)
    wp1 => Set_wps(Particles,from_id,read_only=.TRUE.)
    xp1 => Set_xp(Particles,read_only=.TRUE.)
    IF (PRESENT(Particles_old)) THEN
        wp2 => Set_wps(Particles_old,from_old_id,read_only=.TRUE.)
        xp2 => Set_xp(Particles,read_only=.TRUE.)
    ELSE
        wp2 => NULL()
        xp2 => NULL()
    ENDIF
    nvlist => NULL()
    vlist => NULL()

    CALL substop(caller,t0,info)

    9999  CONTINUE ! jump here upon error


END SUBROUTINE particles_apply_dcops


END MODULE ppm_module_particles
