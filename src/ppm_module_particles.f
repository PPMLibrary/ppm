MODULE ppm_module_particles

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

#define __KIND __DOUBLE_PRECISION

USE ppm_module_particles_typedef
USE ppm_module_typedef
USE ppm_module_alloc
USE ppm_module_substart
USE ppm_module_substop
USE ppm_module_error
USE ppm_module_write
USE ppm_module_data, ONLY: ppm_dim

IMPLICIT NONE
#include "ppm_define.h"

#if   __KIND == __SINGLE_PRECISION
INTEGER, PARAMETER,PRIVATE :: prec = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
INTEGER, PARAMETER,PRIVATE :: prec = ppm_kind_double
#endif

    INTEGER, PARAMETER :: ppm_param_part_init_cartesian = 1
    INTEGER, PARAMETER :: ppm_param_part_init_random = 2


    !----------------------------------------------------------------------
    ! Private variables for the module
    !----------------------------------------------------------------------
    !for debugging
    LOGICAL,PRIVATE           :: verbose = .FALSE.

    INTEGER , PRIVATE, DIMENSION(3)    :: ldc
    !!! Number of elements in all dimensions for allocation

CONTAINS

FUNCTION get_xp(Particles,with_ghosts)
    TYPE(ppm_t_particles)            :: Particles
    LOGICAL,OPTIONAL                 :: with_ghosts
    REAL(prec),DIMENSION(:,:),POINTER:: get_xp

    IF (PRESENT(with_ghosts)) THEN
        IF (with_ghosts) THEN
            IF (Particles%xp_g.EQ.1) THEN
                get_xp => Particles%xp(1:ppm_dim,1:Particles%Mpart)
            ELSE
                write(*,*) 'WARNING: tried to get xp with ghosts &
                    & when ghosts are not up-to-date'
                get_xp => NULL()
            ENDIF
            RETURN
        ENDIF
    ENDIF

    get_xp => Particles%xp(1:ppm_dim,1:Particles%Npart)

END FUNCTION get_xp

FUNCTION set_xp(Particles,read_only,ghosts_ok)
    TYPE(ppm_t_particles)            :: Particles
    LOGICAL,OPTIONAL                 :: read_only
    LOGICAL,OPTIONAL                 :: ghosts_ok
    REAL(prec),DIMENSION(:,:),POINTER:: set_xp

    IF (PRESENT(ghosts_ok)) THEN
        IF (ghosts_ok) THEN
            set_xp => NULL()
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
    set_xp => NULL()

END FUNCTION set_xp

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
                        get_wps => Particles%wps(wps_id)%vec(1:Particles%Mpart)
                    ELSE
                        write(*,*) 'ERROR: tried to get wps with ghosts &
                            & when ghosts are not up-to-date. Returning NULL pointer'
                        write(*,*) 'Run with traceback option to debug'
                        get_wps => NULL()
                    ENDIF
                    RETURN
                ENDIF
            ENDIF
            get_wps => Particles%wps(wps_id)%vec(1:Particles%Npart)
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
    INTEGER                          :: wpv_id,lda
    LOGICAL,OPTIONAL                 :: with_ghosts
    REAL(prec),DIMENSION(:,:),POINTER:: get_wpv

    IF (wpv_id .LE. 0) THEN
        write(*,*) 'ERROR: failed to get wpv for property &
            & wpv_id = ',wpv_id
        get_wpv => NULL()
        stop
        RETURN
    ENDIF

    lda = Particles%wpv_s(wpv_id)

    IF (wpv_id .LE. Particles%max_wpvid) THEN
        IF (Particles%wpv_m(wpv_id).EQ.1) THEN
            IF (PRESENT(with_ghosts)) THEN
                IF (with_ghosts) THEN
                    IF (Particles%wpv_g(wpv_id).EQ.1) THEN
                        get_wpv => &
                            Particles%wpv(wpv_id)%vec(1:lda,1:Particles%Mpart)
                    ELSE
                        write(*,*) 'ERROR: tried to get wpv with ghosts &
                            & when ghosts are not up-to-date. Returning NULL pointer'
                        write(*,*) 'Run with traceback option to debug'
                        get_wpv => NULL()
                    ENDIF
                    RETURN
                ENDIF
            ENDIF
            get_wpv => Particles%wpv(wpv_id)%vec(1:lda,1:Particles%Npart)
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

FUNCTION get_dcop(Particles,eta_id)
    TYPE(ppm_t_particles)              :: Particles
    INTEGER                            :: eta_id
    REAL(prec),DIMENSION(:,:),POINTER  :: get_dcop

    IF (eta_id .LE. 0 .OR. eta_id .GT. Particles%ops%max_opsid) THEN
        write(*,*) 'ERROR: failed to get operator for id ',eta_id
        get_dcop => NULL()
        RETURN
    ENDIF

    get_dcop => Particles%ops%ker(eta_id)%vec(:,1:Particles%Npart)

END FUNCTION get_dcop

FUNCTION set_dcop(Particles,eta_id)
    TYPE(ppm_t_particles)             :: Particles
    INTEGER                           :: eta_id
    REAL(prec),DIMENSION(:,:),POINTER :: set_dcop

    set_dcop => NULL()

END FUNCTION set_dcop

SUBROUTINE ppm_alloc_particles(Particles,Npart,iopt,info,name)

    !!! (Re)allocates the memory of ppm_t_particles data type

    !-------------------------------------------------------------------------
    !  Includes
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    !  Modules
    !-------------------------------------------------------------------------
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
    !  Optional arguments
    !-------------------------------------------------------------------------
    CHARACTER(LEN=*) , OPTIONAL                     :: name
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
            NULLIFY(Particles)
            IF (info .NE. 0) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_dealloc,caller,   &
                    &          'Deallocating Particles',__LINE__,info)
                GOTO 9999
            ENDIF
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
        NULLIFY(Particles%wps)
        NULLIFY(Particles%wps_g)
        NULLIFY(Particles%wps_m)
        NULLIFY(Particles%wpv)
        NULLIFY(Particles%wpv_g)
        NULLIFY(Particles%wpv_m)
        NULLIFY(Particles%wpv_s)
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
        Particles%xp_g = 0

        ! Give a default name to this Particle set
        IF (PRESENT(name)) THEN
            Particles%name = name
        ELSE
            Particles%name = particles_dflt_partname()
        ENDIF
        ! No properties defined yet
        Particles%nwps = 0
        Particles%nwpv = 0
        Particles%max_wpsid = 0
        Particles%max_wpvid = 0
        ! No active topology yet
        Particles%active_topoid = -1
        ! Particles are not yet mapped on any topology
        Particles%ontopology = .FALSE.
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


        IF (verbose) &
            write(*,*) 'Allocated particles with Np=',Npart


    ENDIF

    !-------------------------------------------------------------------------
    !  Finalize
    !-------------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999 CONTINUE

END SUBROUTINE ppm_alloc_particles

SUBROUTINE particles_dcop_deallocate(Particles,info)

    !!! Deallocates the memory of ppm_t_operator data type

    !-------------------------------------------------------------------------
    !  Includes
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    !  Modules
    !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    TYPE(ppm_t_particles), POINTER,  INTENT(INOUT)  :: Particles
    !!! Data structure containing the operators
    INTEGER,                         INTENT(  OUT)  :: info
    !!! Returns status, 0 upon success.
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------

    INTEGER                                         :: i
    REAL(KIND(1.D0))                                :: t0
    CHARACTER(LEN = ppm_char)            :: caller = 'ppm_dcop_deallocate'
    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)

    IF (ASSOCIATED(Particles%ops)) THEN
        IF (ASSOCIATED(Particles%ops%desc)) THEN
            DO i=1,Particles%ops%max_opsid
                IF (ASSOCIATED(Particles%ops%desc(i)%degree)) THEN
                    CALL particles_dcop_free(Particles,i,info)
                    IF (info .NE. 0) THEN
                        info = ppm_error_error
                        CALL ppm_error(ppm_err_dealloc,caller,   &
                            &    'freeing Particles%ops%desc(i)',__LINE__,info)
                        GOTO 9999
                    ENDIF
                ENDIF
            ENDDO
            DEALLOCATE(Particles%ops%desc,STAT=info)
            IF (info .NE. 0) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_dealloc,caller,   &
                    &          'Deallocating Particles%ops%desc',__LINE__,info)
                GOTO 9999
            ENDIF
            DEALLOCATE(Particles%ops%ker,STAT=info)
            IF (info .NE. 0) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_dealloc,caller,   &
                    &          'Deallocating Particles%ops%ker',__LINE__,info)
                GOTO 9999
            ENDIF
            NULLIFY(Particles%ops%ker,Particles%ops%desc)
        ENDIF
        DEALLOCATE(Particles%ops,STAT=info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_dealloc,caller,   &
                &          'Deallocating Particles%ops',__LINE__,info)
            GOTO 9999
        ENDIF
        NULLIFY(Particles%ops)
    ENDIF
    !-------------------------------------------------------------------------
    !  Finalize
    !-------------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999 CONTINUE

END SUBROUTINE particles_dcop_deallocate

FUNCTION particles_dflt_partname(i)
    !!! Default name for a set of particles
    CHARACTER(LEN=ppm_char)   :: particles_dflt_partname
    INTEGER,OPTIONAL          :: i
    CHARACTER(LEN=ppm_char)   :: buf

    IF (PRESENT(i)) THEN
        WRITE(buf,*) i
        WRITE(particles_dflt_partname,*) 'Particles_',ADJUSTL(TRIM(buf))
    ELSE
        WRITE(particles_dflt_partname,*) 'Particles'
    ENDIF
    RETURN
END FUNCTION

FUNCTION particles_dflt_pptname(i,ndim)
    !!! Default name for a scalar or vector property
    CHARACTER(LEN=ppm_char)   :: particles_dflt_pptname
    INTEGER                   :: i,ndim
    CHARACTER(LEN=ppm_char)   :: buf

    WRITE(buf,*) i
    IF (ndim .EQ. 1) THEN
        WRITE(particles_dflt_pptname,*) 'property_s',ADJUSTL(TRIM(buf))
    ELSE
        WRITE(particles_dflt_pptname,*) 'property_v',ADJUSTL(TRIM(buf))
    ENDIF
    RETURN
END FUNCTION

FUNCTION particles_dflt_opname(i)
    !!! Default name for an operator
    CHARACTER(LEN=ppm_char)   :: particles_dflt_opname
    INTEGER                   :: i
    CHARACTER(LEN=ppm_char)   :: buf

    WRITE(buf,*) i
    WRITE(particles_dflt_opname,*) 'operator_',ADJUSTL(TRIM(buf))
    RETURN
END FUNCTION


SUBROUTINE particles_allocate_wps(Particles,wp_id,info,with_ghosts,&
        zero,iopt,name)

    !!! Allocate a scalar property in this Particles data structure
    !!! If _wp_id_ is zero, a new id is generated
    !!! otherwise, the property wp_id is overwritten
    !!! The array can be of size Npart or Mpart (if _with_ghosts_ is
    !!! true) and initialized to zero (if _zero_ is true) or not.

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
    CHARACTER(LEN=*) , OPTIONAL                     :: name
    !!! give a name to this scalar-valued property
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
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
        DO i=1,20
            Particles%wps(i)%name = particles_dflt_pptname(i,1)
        ENDDO
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
            info = ppm_error_error
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
        !reset the label of this property to its default value
        !Commented out for now - not sure whether we want this or no
        !Particles%wps(wp_id)%name = particles_dflt_pptname(wp_id,1)
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
        CALL ppm_alloc(Particles%wps(wp_id)%vec,ldc(1:1),iopt_l,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,   &
                &            'Could not allocate wps elements',__LINE__,info)
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
        !Set the name of the property
        IF (PRESENT(name)) Particles%wps(wp_id)%name = TRIM(name)
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
            write(*,*) 'Allocated wps with size=',ldc(1),&
            ' id=',wp_id, 'name = ', TRIM(ADJUSTL(Particles%wps(wp_id)%name))
    ENDIF

    !-----------------------------------------------------------------------
    ! Finalize
    !-----------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999 CONTINUE

END SUBROUTINE particles_allocate_wps


SUBROUTINE particles_allocate_wpv(Particles,wp_id,lda,info, & 
        with_ghosts,zero,iopt,name)

    !!! Allocate a vector valued property in this Particles data structure
    !!! If _wp_id_ is zero, a new id is generated
    !!! otherwise, the property wp_id is overwritten
    !!! The array can be of size Npart or Mpart (if _with_ghosts_ is
    !!! true) and initialized to zero (if _zero_ is true) or not.

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
    CHARACTER(LEN=*) , OPTIONAL                     :: name
    !!! give a name to this vector-valued property
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
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
        DO i=1,20
            Particles%wpv(i)%name = particles_dflt_pptname(i,2)
        ENDDO
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
        CALL ppm_alloc(Particles%wpv(wp_id)%vec,ldc(1:2),ppm_param_dealloc,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_dealloc,caller,&
                'ppm_alloc (deallocate) failed',__LINE__,info)
            GOTO 9999
        ENDIF

        ! Update state
        Particles%wpv_m(wp_id) = -1
        Particles%wpv_g(wp_id) = -1
        Particles%wpv_s(wp_id) = 0
        Particles%wpv(wp_id)%vec => NULL()
        !reset the label of this property to its default value
        !Commented out for now - not sure whether we want this or no
        !Particles%wpv(wp_id)%name = particles_dflt_pptname(wp_id,2)
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
                &            'Could not allocate wpv elements',__LINE__,info)
            GOTO 9999
        ENDIF
        IF (PRESENT(zero)) THEN
            IF (zero) THEN
                DO i=1,ldc(2)
                    Particles%wpv(wp_id)%vec(:,i)=0._MK
                ENDDO
            ENDIF
        ENDIF

        !-----------------------------------------------------------------
        ! Update state
        !-----------------------------------------------------------------
        !Set the name of the property
        IF (PRESENT(name)) Particles%wpv(wp_id)%name = TRIM(name)
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
            write(*,*) 'Allocated wpv with size=',ldc,&
            ' id=',wp_id, 'name = ', TRIM(ADJUSTL(Particles%wpv(wp_id)%name))

    ENDIF

    !-----------------------------------------------------------------------
    ! Finalize
    !-----------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999 CONTINUE ! jump here upon error

END SUBROUTINE particles_allocate_wpv

SUBROUTINE particles_mapping_global(Particles,topoid,info)

    !!!  Global mapping for particles
    !!!  Assumptions:
    !!! * All the particles have to be inside the domain
    !!!   (otherwise -> "unassigned particle error")

    USE ppm_module_map
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
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
    IF (.NOT.Particles%areinside) THEN
        info = ppm_error_error
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
            info = ppm_error_error
            CALL ppm_error(999,caller,&
                'ppm_map_part_global failed',__LINE__,info)
            GOTO 9999
        ENDIF

        DO prop_id = 1,Particles%max_wpsid
            IF(Particles%wps_m(prop_id) .EQ. 0) THEN
                CALL ppm_map_part_push(Particles%wps(prop_id)%vec,&
                    Particles%Npart,info)
                IF (info .NE. 0) THEN
                    info = ppm_error_error
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
                    info = ppm_error_error
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
                    info = ppm_error_error
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
                    info = ppm_error_error
                    CALL ppm_error(999,caller,&
                        'ppm_map_part_pop failed',__LINE__,info)
                    GOTO 9999
                ENDIF
            ENDIF
        ENDDO

        CALL ppm_map_part_pop(Particles%xp,ppm_dim,Particles%Npart,&
            Npart_new,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
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

    !!!  Partial mapping for particles
    !!!  Assumptions:
    !!! * All the particles have to be inside the domain
    !!!   (otherwise -> "unassigned particle error")

    USE ppm_module_map
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
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
    IF (.NOT.Particles%areinside) THEN
        info = ppm_error_error
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
            info = ppm_error_error
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
                    info = ppm_error_error
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
                    info = ppm_error_error
                    CALL ppm_error(999,caller,&
                        'ppm_map_part_push failed',__LINE__,info)
                    GOTO 9999
                ENDIF
            ENDIF
        ENDDO

        CALL ppm_map_part_send(Particles%Npart,Npart_new,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
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
                    info = ppm_error_error
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
                    info = ppm_error_error
                    CALL ppm_error(999,caller,&
                        'ppm_map_part_pop failed',__LINE__,info)
                    GOTO 9999
                ENDIF
            ENDIF
        ENDDO

        CALL ppm_map_part_pop(Particles%xp,ppm_dim,Particles%Npart,&
            Npart_new,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
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

    !!!  Ghost mapping for particles
    !!!  Assumptions:
    !!! * Particles positions need to have been mapped onto the topology
    !!!
    USE ppm_module_data, ONLY: ppm_topo
    USE ppm_module_map
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
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (.NOT.Particles%ontopology .OR. Particles%active_topoid.NE.topoid) THEN
        !Particles have not been mapped onto this topology
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'Do a partial/global mapping before doing a ghost mapping',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
    IF (.NOT.Particles%areinside) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'some particles may be outside the domain. Apply BC first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF


    topo=>ppm_topo(topoid)%t

    cutoff = Particles%cutoff + Particles%skin 

#if   __KIND == __SINGLE_PRECISION
    IF (cutoff .GT. topo%ghostsizes) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'ghostsize of topology may be smaller than that of particles',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
#elif   __KIND == __DOUBLE_PRECISION
    IF (cutoff .GT. topo%ghostsized) THEN

        write(*,*) cutoff
        write(*,*) topo%ghostsized

        info = ppm_error_error
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
                info = ppm_error_error
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
                        info = ppm_error_error
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
                        info = ppm_error_error
                        CALL ppm_error(999,caller,&
                            'ppm_map_part_push failed',__LINE__,info)
                        GOTO 9999
                    ENDIF
                ENDIF
            ENDIF
        ENDDO

        CALL ppm_map_part_send(Particles%Npart,Particles%Mpart,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
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
                        info = ppm_error_error
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
                        info = ppm_error_error
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
                info = ppm_error_error
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

    !!!  Apply boundary conditions for particles positions
    !!!  Assumptions:
    !!! * Particles positions need to have been mapped onto the topology

    USE ppm_module_data, ONLY: ppm_topo,ppm_rank
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
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (Particles%active_topoid.NE.topoid) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'WARNING: this topoid is not the one for which &
            & these particles were mapped',&
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
                IF (xp(di,ip) .EQ. max_phys(di)) &
                    xp(di,ip) = xp(di,ip) - len_phys(di)*almostone
                IF (xp(di,ip) .GT. max_phys(di)) &
                    xp(di,ip) = xp(di,ip) - len_phys(di)
                IF (xp(di,ip) .LT. min_phys(di)) &
                    xp(di,ip) = xp(di,ip) + len_phys(di)
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
            info = ppm_error_error
            CALL ppm_error(999,caller,   &
                & 'this type of BC is not implemented/tested in this version',&
                &  __LINE__,info)
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

SUBROUTINE particles_move(Particles,disp,info)

    !!!  Move all particles according to some displacement field
    !!!  The size of disp must match the size of xp

    USE ppm_module_data, ONLY: ppm_topo,ppm_rank
#if   __KIND == __SINGLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    TYPE(ppm_t_particles),    POINTER,  INTENT(INOUT)     :: Particles
    !!! Data structure containing the particles
    REAL(MK), DIMENSION(:,:), POINTER,  INTENT(IN   )     :: disp
    !!! Data structure containing the particles
    INTEGER,                            INTENT(  OUT)     :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    REAL(ppm_kind_double)                                 :: t0
    INTEGER                                               :: ip
    CHARACTER(LEN = ppm_char)                 :: caller ='particles_move'
    REAL(MK),DIMENSION(:,:),POINTER                       :: xp=>NULL()
    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)
    !-----------------------------------------------------------------
    !  checks
    !-----------------------------------------------------------------
    IF (.NOT.ASSOCIATED(Particles)) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF


    xp => get_xp(Particles)
    FORALL (ip=1:Particles%Npart) &
            xp(1:ppm_dim,ip) = xp(1:ppm_dim,ip) + disp(1:ppm_dim,ip)
    xp => set_xp(Particles)

    !-----------------------------------------------------------------
    !  update states
    !-----------------------------------------------------------------
    CALL particles_have_moved(Particles,info)

    !-----------------------------------------------------------------
    !  Finalize
    !-----------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999 CONTINUE ! jump here upon error

END SUBROUTINE particles_move

SUBROUTINE particles_have_moved(Particles,info)
    !-----------------------------------------------------------------
    !  Update states when particles have moved 
    !   (assuming they have moved across processor boundaries)
    !-----------------------------------------------------------------
    USE ppm_module_data, ONLY: ppm_topo,ppm_rank
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
        info = ppm_error_error
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
    Particles%cartesian = .FALSE.

    !-----------------------------------------------------------------
    !  Finalize
    !-----------------------------------------------------------------
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
    USE ppm_module_neighlist
    USE ppm_module_inl_vlist
#ifdef __MPI
    INCLUDE "mpif.h"
#endif
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
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (.NOT.Particles%ontopology .OR. Particles%active_topoid.NE.topoid) THEN
        !Particles have not been mapped onto this topology
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'Do a partial/global mapping before',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (Particles%xp_g.NE.1) THEN
        info = ppm_error_error
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
                info = ppm_error_error
                CALL ppm_error(999,caller,&
                    'ppm_inl_vlist failed',__LINE__,info)
                GOTO 9999
            ENDIF
        ELSE
            CALL ppm_neighlist_vlist(topoid,Particles%xp,Particles%Mpart,&
                Particles%cutoff,Particles%skin,symmetry,Particles%vlist,&
                Particles%nvlist,info,lstore=lstore)
            IF (info .NE. 0) THEN
                info = ppm_error_error
                CALL ppm_error(999,caller,&
                    'ppm_neighlist_vlist failed',__LINE__,info)
                GOTO 9999
            ENDIF
        ENDIF
    ENDIF

    !Update state
    Particles%neighlists = .TRUE.
    !

    Particles%nneighmin = MINVAL(Particles%nvlist(1:Particles%Npart))
    Particles%nneighmax = MAXVAL(Particles%nvlist(1:Particles%Npart))

    ! TODO: maybe the MPI_Allreduce is not really needed for production runs
    ! This is mainly used for debugging/warnings
#ifdef __MPI
    CALL MPI_Allreduce(Particles%nneighmin,Particles%nneighmin,1,&
        MPI_INTEGER,MPI_MIN,ppm_comm,info)
    CALL MPI_Allreduce(Particles%nneighmax,Particles%nneighmax,1,&
        MPI_INTEGER,MPI_MAX,ppm_comm,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,'MPI_Allreduce failed',__LINE__,info)
        GOTO 9999
    ENDIF
#endif

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
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
    IF (.NOT.ASSOCIATED(Particles_1%xp).OR..NOT.ASSOCIATED(Particles_2%xp)) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'Particles positions had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (.NOT.Particles_1%ontopology .OR. Particles_1%active_topoid.NE.topoid) THEN
        !Particles have not been mapped onto this topology
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'Do a partial/global mapping for Particles_1 before',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
    IF (.NOT.Particles_2%ontopology .OR. Particles_2%active_topoid.NE.topoid) THEN
        !Particles have not been mapped onto this topology
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'Do a partial/global mapping for Particles_2 before',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (Particles_2%xp_g.NE.1) THEN
        info = ppm_error_error
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
                info = ppm_error_error
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
                info = ppm_error_error
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
    Particles_1%Particles_cross => Particles_2

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
#ifdef __MPI
    INCLUDE "mpif.h"
#endif
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
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    ! Update states
    ! particles may be outside of the computational domain
    Particles%areinside=.FALSE.
    ! particles may no longer be on the right processor
    Particles%ontopology=.FALSE.
    Particles%cartesian = .FALSE.
    ! note that there is still a one-to-one relationship between real particles
    !  and indices. They can still be used for computing stuff or for
    ! being written to output.
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
#ifdef __MPI
    INCLUDE "mpif.h"
#endif
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
        info = ppm_error_error
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
            info = ppm_error_error
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
            info = ppm_error_error
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

SUBROUTINE particles_updated_cutoff(Particles,info,cutoff)
    !-----------------------------------------------------------------
    !  routine to call when cutoffs have changed
    !-----------------------------------------------------------------
    !  Effects:
    ! * update cutoff to maxval(rcp) or to cutoff (if present)
    ! * discard neighbour lists
    !  Assumptions:
    ! * Particles positions need to have been mapped onto the topology
    ! * cutoff radii are up-to-date
    !
    ! Note:
    ! * Does an MPI_Allreduce for to get the maximum cutoff
    USE ppm_module_data, ONLY: ppm_comm,ppm_mpi_kind
#ifdef __MPI
    INCLUDE "mpif.h"
#endif
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
    REAL(MK),OPTIONAL                                      :: cutoff
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    !!!
    CHARACTER(LEN = ppm_char)             :: caller = 'particles_updated_cutoff'
    REAL(MK)                              :: cutoff_new
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
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (PRESENT(cutoff)) THEN
        cutoff_new = cutoff
        IF (Particles%adaptive) THEN
            Particles%wps(Particles%rcp_id)%vec = cutoff
        ENDIF
    ELSE
        IF (Particles%adaptive) THEN

            cutoff_new = MAXVAL(Particles%wps(Particles%rcp_id)%vec(1:Particles%Npart))

#ifdef __MPI
            CALL MPI_Allreduce(cutoff_new,cutoff_new,1,ppm_mpi_kind,&
                MPI_MAX,ppm_comm,info)
            IF (info .NE. 0) THEN
                info = ppm_error_error
                CALL ppm_error(999,caller,'MPI_Allreduce failed',__LINE__,info)
                GOTO 9999
            ENDIF
#endif
        ELSE
            info = ppm_error_error
            CALL ppm_error(999,caller,'need to provide new cutoff value',&
                &  __LINE__,info)
            GOTO 9999
        ENDIF
    ENDIF

    ! Update states
    ! If cutoff has increased, then ghosts may no longer be 
    ! up to date
    IF (Particles%cutoff .LT. cutoff_new) THEN
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
    Particles%cutoff = cutoff_new

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
#ifdef __MPI
    INCLUDE "mpif.h"
#endif
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
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (.NOT.Particles%neighlists) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'Please construct neighbour lists first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
    IF (Particles%xp_g.NE.1) THEN
        info = ppm_error_error
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

#ifdef __MPI
    CALL MPI_Allreduce(hmin2,hmin2,1,ppm_mpi_kind,MPI_MIN,ppm_comm,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'MPI_Allreduce failed',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
#endif


    ! Update states
    Particles%h_min = SQRT(hmin2)

    !-----------------------------------------------------------------------
    ! Finalize
    !-----------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999 CONTINUE ! jump here upon error

END SUBROUTINE particles_compute_hmin

#define __DIM 2
#include "part/ppm_particles_initialize.f"
#define __DIM 3
#include "part/ppm_particles_initialize.f"


!!temporary hack to deal with both 2d and 3d
SUBROUTINE particles_initialize(Particles,Npart_global,info,&
        distrib,topoid,minphys,maxphys,cutoff)
    !-----------------------------------------------------------------------
    ! Set initial particle positions
    !-----------------------------------------------------------------------
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

    IF (ppm_dim .eq. 2) THEN
        CALL  particles_initialize2d(Particles,Npart_global,info,&
            distrib,topoid,minphys,maxphys,cutoff)
    ELSE
        CALL  particles_initialize3d(Particles,Npart_global,info,&
            distrib,topoid,minphys,maxphys,cutoff)
    ENDIF
END SUBROUTINE particles_initialize


SUBROUTINE particles_io_xyz(Particles,itnum,writedir,info,&
        with_ghosts,wps_list,wpv_list)
    !!!------------------------------------------------------------------------!
    !!! Dump data files to disk
    !!! All data is gathered to rank 0 and then written as one single file.
    !!!
    !!! NOT to be used with large number of processors!...
    !!!  use the vtk module with ppm_io instead
    !!!------------------------------------------------------------------------!
    USE ppm_module_data, ONLY: ppm_rank,ppm_nproc,ppm_comm

#ifdef __MPI
    INCLUDE "mpif.h"
#endif

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
    CHARACTER(LEN=ppm_char)                    :: cbuf
    CHARACTER(LEN=2*ppm_char)                  :: filename
    CHARACTER(LEN=ppm_char)                    :: caller = 'particles_io_xyz'
    CHARACTER(LEN=128)                         :: myformat
    INTEGER,DIMENSION(:),ALLOCATABLE           :: Npart_vec
    REAL(MK),DIMENSION(:,:),ALLOCATABLE        :: propp_iproc
    REAL(MK),DIMENSION(:,:),ALLOCATABLE        :: xp_iproc
    INTEGER                                    :: count,iproc,i,j,ip,di
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
        info = ppm_error_error
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
                info = ppm_error_error
                CALL ppm_error(999,caller,   &
                    &  'property index exceeds size of property array',&
                    &  __LINE__,info)
                GOTO 9999
            ENDIF
            IF (Particles%wps_m(wps_l(i)).NE.1) THEN
                info = ppm_error_error
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
            IF (Particles%wps_m(i).EQ.1) THEN
                nb_wps = nb_wps + 1
                wps_l(nb_wps) = i
            ENDIF
        ENDDO
    ENDIF
    IF(PRESENT(wpv_list)) THEN
        nb_wpv=SIZE(wpv_list)
        ALLOCATE(wpv_l(nb_wpv),STAT=info)
        wpv_l=wpv_list
        DO i=1,nb_wpv
            IF (wpv_l(i).GT.Particles%max_wpvid) THEN
                info = ppm_error_error
                CALL ppm_error(999,caller,   &
                    &  'property index exceeds size of property array',&
                    &  __LINE__,info)
                GOTO 9999
            ENDIF
            IF (Particles%wpv_m(wpv_l(i)).NE.1) THEN
                info = ppm_error_error
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
            IF (Particles%wpv_m(i).EQ.1) THEN
                nb_wpv = nb_wpv + 1
                wpv_l(nb_wpv) = i
            ENDIF
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
    WRITE(myformat,'(A,I4,A)') '( ',ppm_dim+ndim2,'E24.13E3)'
#else
    WRITE(myformat,'(A,I4,A)') '( ',ppm_dim+ndim2,'E17.7E3)'
#endif

    !====================================================================!
    ! open file
    IF (ppm_rank .EQ. 0) THEN
        WRITE(filename,'(A,A,A,I6.6,A)') TRIM(ADJUSTL(writedir)),&
            TRIM(ADJUSTL(Particles%name)),'_',itnum,'.xyz'
        OPEN(23,FILE=TRIM(ADJUSTL(filename)),FORM='FORMATTED',STATUS='REPLACE',IOSTAT=info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(999,caller,'failed to open file 23',__LINE__,info)
            GOTO 9999
        ENDIF
        WRITE(filename,'(A,A,A,I6.6,A)') TRIM(ADJUSTL(writedir)),&
            TRIM(ADJUSTL(Particles%name)),'_',itnum,'.txt'
        OPEN(24,FILE=TRIM(filename),FORM='FORMATTED',STATUS='REPLACE',IOSTAT=info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(999,caller,'failed to open file 24',__LINE__,info)
            GOTO 9999
        ENDIF
    ENDIF

    !====================================================================!
    ! rank 0 collects data and writes out sequentially
    ! FIXME
    ALLOCATE(Npart_vec(ppm_nproc),STAT=info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,'allocation failed',__LINE__,info)
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
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,'MPI_Gather failed',__LINE__,info)
        GOTO 9999
    ENDIF

    IF (ppm_rank .EQ. 0) THEN

        ip = 0
        DO iproc = 1,ppm_nproc

            IF (iproc .GT. 1) THEN

                ALLOCATE(xp_iproc(ppm_dim,Npart_vec(iproc)),STAT=info)
                IF (info .NE. 0) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_alloc,caller,'allocation failed',__LINE__,info)
                    GOTO 9999
                ENDIF
                ALLOCATE(propp_iproc(ndim2,Npart_vec(iproc)),STAT=info)
                IF (info .NE. 0) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_alloc,caller,'allocation failed',__LINE__,info)
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
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_alloc,caller,'allocation failed',__LINE__,info)
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

            !Write a small description file
            WRITE (24, '(A)') 'Column| Variable name'
            WRITE (24, '(A)') '---------------------'
            DO di=1,ppm_dim
                WRITE (24, '(I2,A,I1)') di,': xp_',di
            ENDDO
            ndim2= ppm_dim+1
            DO i=1,nb_wps
                WRITE (24, '(I2,A,A)') ndim2,': ',&
                    TRIM(ADJUSTL(Particles%wps(wps_l(i))%name))
                ndim2 = ndim2 + 1
            ENDDO
            DO i=1,nb_wpv
                wpv_s=Particles%wpv_s(wpv_l(i))
                DO di=1,wpv_s
                    WRITE (24, '(I2,A,A,A,I2.2)') ndim2,': ',&
                        TRIM(ADJUSTL(Particles%wpv(wpv_l(i))%name)),'_',di
                    ndim2=ndim2+1
                ENDDO
            ENDDO
            IF (ghosts) THEN
                WRITE (24, '(I2,A,I1)') ndim2,': is_a_ghost'
            ENDIF
        ENDDO

        !==================================================================!
        ! close file

        CLOSE(23,IOSTAT=info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,'failed to close file',__LINE__,info)
            GOTO 9999
        ENDIF
        CLOSE(24,IOSTAT=info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,'failed to close file',__LINE__,info)
            GOTO 9999
        ENDIF

    ELSE

        count = ppm_dim*Npart_local
        IF (MK .EQ. 4) THEN
            CALL MPI_Send(Particles%xp,count,MPI_REAL,0,ppm_rank,ppm_comm,info)
            IF (info .NE. 0) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_alloc,caller,'MPI_Send failed',__LINE__,info)
                GOTO 9999
            ENDIF
        ELSE
            CALL MPI_Send(Particles%xp,count,MPI_DOUBLE_PRECISION,0,&
                ppm_rank,ppm_comm,info)
            IF (info .NE. 0) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_alloc,caller,'MPI_Send failed',__LINE__,info)
                GOTO 9999
            ENDIF
        ENDIF

        count = ndim2*Npart_local
        ALLOCATE(propp_iproc(ndim2,Npart_local),STAT=info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,'allocation failed',__LINE__,info)
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
                info = ppm_error_error
                CALL ppm_error(ppm_err_alloc,caller,'MPI_Send failed',__LINE__,info)
                GOTO 9999
            ENDIF
        ELSE
            CALL MPI_Send(propp_iproc,count,MPI_DOUBLE_PRECISION,0,&
                ppm_rank,ppm_comm,info)
            IF (info .NE. 0) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_alloc,caller,'MPI_Send failed',__LINE__,info)
                GOTO 9999
            ENDIF
        ENDIF

        DEALLOCATE(propp_iproc)
    ENDIF

    DEALLOCATE(Npart_vec)


    IF (ppm_rank.EQ.0) THEN
        WRITE(cbuf,*) 'Wrote ',ip,' particles to file'
        CALL ppm_write(ppm_rank,caller,cbuf,info)
    ENDIF

    !-------------------------------------------------------------------------
    !  Finalize
    !-------------------------------------------------------------------------
    DEALLOCATE(wps_l,wpv_l)
    CALL substop(caller,t0,info)

    9999  CONTINUE ! jump here upon error

END SUBROUTINE particles_io_xyz

SUBROUTINE particles_dcop_define(Particles,eta_id,coeffs,degree,order,nterms,&
        info,name,interp,vector)
    !!!------------------------------------------------------------------------!
    !!! Define a DC operator as a linear combination (with scalar coefficients)
    !!! of nterms partial derivatives of arbitrary degrees. These are given by a matrix
    !!! of integers where each row represents one term of the linear combination
    !!! and each of the ppm_dim columns is the order of differentiation in that
    !!! dimension.
    !!! The definition of the operator is stored in the ppm_t_operator derived 
    !!! type under the index eta_id. 
    !!! The operator itself is computed elsewhere and will be stored in 
    !!! the same data structure.
    !!!
    !!! Usage example:
    !!!
    !!!   The differential operator:
    !!!   3.0 df/dx -7.0 d^4f/dxdydz^2 + 8.0 d^3f/dx^2dz
    !!!   would be defined by calling particles_dcop_define with
    !!!   coeffs = (/3.0, -7.0, 8.0/)
    !!!   degree = (/1,0,0,  1,1,2,  2,0,1 /)
    !!!   order =  (/2,      1,      3     /)
    !!!   nterms = 3
    !!!------------------------------------------------------------------------!
#ifdef __MPI
    INCLUDE "mpif.h"
#endif

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
    INTEGER,                            INTENT(INOUT)   :: eta_id
    !!! id where the data is stored
    REAL(MK),DIMENSION(:),              INTENT(IN   )   :: coeffs
    !!! Multiplicative coefficients of each term in the linear combination of
    !!! differential operators
    INTEGER,DIMENSION(:),               INTENT(IN   )   :: degree
    !!! Degree of differentiation of each term
    INTEGER,DIMENSION(:),               INTENT(IN   )   :: order
    !!! Order of approxmiation for each term
    INTEGER,                            INTENT(IN   )   :: nterms
    !!! Number of terms in the linear combination
    INTEGER,                            INTENT(  OUT)   :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    ! Optional argument
    !-------------------------------------------------------------------------
    CHARACTER(LEN=*),OPTIONAL,          INTENT(IN   )   :: name
    !!! Descriptive name for this operator
    LOGICAL,OPTIONAL,                   INTENT(IN   )   :: interp
    !!! True if the operator interpolates data from one set of particles to
    !!! another. Default is false.
    LOGICAL,OPTIONAL,                   INTENT(IN   )   :: vector
    !!! True if the operator is a vector field. Default is false.
    !-------------------------------------------------------------------------
    ! local variables
    !-------------------------------------------------------------------------
    CHARACTER(LEN = ppm_char)               :: caller = 'particles_dcop_define'
    INTEGER,DIMENSION(3)                       :: ldc
    REAL(KIND(1.D0))                           :: t0
    INTEGER                                    :: i


    !-------------------------------------------------------------------------
    ! Initialize
    !-------------------------------------------------------------------------
    info = 0 ! change if error occurs
    CALL substart(caller,t0,info)

    IF (.NOT. ASSOCIATED(Particles%ops)) THEN
    !allocate operators data structure
        ALLOCATE(Particles%ops,STAT=info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,   &
                &            'failed to allocate Particles%ops',__LINE__,info)
            GOTO 9999
        ENDIF
        Particles%ops%nb_ops = 0
        Particles%ops%max_opsid = 0
        ALLOCATE(Particles%ops%ker(20),STAT=info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,   &
                &            'failed to allocate Particles%ops%ker',__LINE__,info)
            GOTO 9999
        ENDIF
        ALLOCATE(Particles%ops%desc(20),STAT=info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,   &
                &            'failed to allocate Particles%ops%desc',__LINE__,info)
            GOTO 9999
        ENDIF
        DO i=1,20
            Particles%ops%desc(i)%name = particles_dflt_opname(i)
            Particles%ops%desc(i)%interp = .FALSE.
            Particles%ops%desc(i)%nterms = 0
            Particles%ops%desc(i)%vector = .FALSE.
        ENDDO
    ENDIF
    IF (MINVAL(degree).LT.0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &            'invalid degree: must be positive',__LINE__,info)
        GOTO 9999
    ENDIF
    IF (MINVAL(order).LT.0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &            'invalid approx order: must be positive',__LINE__,info)
        GOTO 9999
    ENDIF
    IF (SIZE(degree).NE.ppm_dim*nterms) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &            'wrong number of terms in degree argument',__LINE__,info)
        GOTO 9999
    ENDIF
    IF (SIZE(order).NE.nterms) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &            'wrong number of terms in order argument',__LINE__,info)
        GOTO 9999
    ENDIF
    IF (SIZE(coeffs).NE.nterms) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &            'wrong number of terms in coeffs argument',__LINE__,info)
        GOTO 9999
    ENDIF

    IF (eta_id.EQ.0) THEN
        eta_id = 1
        DO WHILE (ASSOCIATED(Particles%ops%desc(eta_id)%degree))
            eta_id = eta_id + 1
        ENDDO
    ELSE
        IF (ASSOCIATED(Particles%ops%desc(eta_id)%degree)) THEN
            CALL particles_dcop_free(Particles,eta_id,info)
            IF (info .NE. 0) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_alloc,caller,   &
                    &            'failed to free Particles%ops',__LINE__,info)
                GOTO 9999
            ENDIF
        ENDIF
    ENDIF


    !allocate operators descriptors
    ldc(1) = ppm_dim * nterms
    CALL ppm_alloc(Particles%ops%desc(eta_id)%degree,ldc(1:1),&
        ppm_param_alloc_fit,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &            'failed to allocate ops%desc',__LINE__,info)
        GOTO 9999
    ENDIF
    ldc(1) = nterms
    CALL ppm_alloc(Particles%ops%desc(eta_id)%order,ldc(1:1),&
        ppm_param_alloc_fit,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &            'failed to allocate ops%desc',__LINE__,info)
        GOTO 9999
    ENDIF
    ldc(1) = nterms
    CALL ppm_alloc(Particles%ops%desc(eta_id)%coeffs,ldc(1:1),&
        ppm_param_alloc_fit,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &            'failed to allocate ops%desc',__LINE__,info)
        GOTO 9999
    ENDIF
    Particles%ops%desc(eta_id)%degree = degree 
    Particles%ops%desc(eta_id)%order = order 
    Particles%ops%desc(eta_id)%coeffs = coeffs 
    Particles%ops%desc(eta_id)%nterms = nterms 
    IF (PRESENT(name)) THEN
        Particles%ops%desc(eta_id)%name = name
    ENDIF
    IF (PRESENT(interp)) THEN
        Particles%ops%desc(eta_id)%interp = interp
    ELSE
        Particles%ops%desc(eta_id)%interp = .FALSE.
    ENDIF
    IF (PRESENT(vector)) THEN
        Particles%ops%desc(eta_id)%vector = vector
    ELSE
        Particles%ops%desc(eta_id)%vector = .FALSE.
    ENDIF

    !update states
    Particles%ops%nb_ops = Particles%ops%nb_ops+1
    IF (eta_id .GT. Particles%ops%max_opsid) &
        Particles%ops%max_opsid = eta_id

    !-------------------------------------------------------------------------
    ! Finalize
    !-------------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999  CONTINUE ! jump here upon error


END SUBROUTINE particles_dcop_define

SUBROUTINE particles_dcop_free(Particles,eta_id,info)
    !!!------------------------------------------------------------------------!
    !!! Remove an operator from the Particles data structure
    !!!------------------------------------------------------------------------!
#ifdef __MPI
    INCLUDE "mpif.h"
#endif

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
    INTEGER,                            INTENT(IN   )   :: eta_id
    !!! id where the data is stored
    INTEGER,                            INTENT(  OUT)   :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    ! local variables
    !-------------------------------------------------------------------------
    CHARACTER(LEN = ppm_char)               :: caller = 'particles_dcop_free'
    INTEGER,DIMENSION(3)                       :: ldc
    REAL(KIND(1.D0))                           :: t0


    !-------------------------------------------------------------------------
    ! Initialize
    !-------------------------------------------------------------------------
    info = 0 ! change if error occurs
    CALL substart(caller,t0,info)

    IF (.NOT. ASSOCIATED(Particles%ops)) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            & 'No operator data structure found - cannot free element',&
            __LINE__,info)
        GOTO 9999
    ENDIF

    IF (eta_id.LE.0 .OR. eta_id .GT. Particles%ops%max_opsid) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            & 'invalid value for ops id - cannot free element',&
            __LINE__,info)
        GOTO 9999
    ENDIF

    !free memory
    ldc = 0
    CALL ppm_alloc(Particles%ops%desc(eta_id)%degree,ldc,&
        ppm_param_dealloc,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &            'failed to deallocate ops%desc',__LINE__,info)
        GOTO 9999
    ENDIF
    CALL ppm_alloc(Particles%ops%desc(eta_id)%order,ldc,&
        ppm_param_dealloc,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &            'failed to deallocate ops%desc',__LINE__,info)
        GOTO 9999
    ENDIF
    CALL ppm_alloc(Particles%ops%desc(eta_id)%coeffs,ldc,&
        ppm_param_dealloc,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &            'failed to deallocate ops%desc',__LINE__,info)
        GOTO 9999
    ENDIF
    CALL ppm_alloc(Particles%ops%ker(eta_id)%vec,ldc,&
        ppm_param_dealloc,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &            'failed to deallocate ops%ker',__LINE__,info)
        GOTO 9999
    ENDIF
    Particles%ops%desc(eta_id)%degree => NULL()
    Particles%ops%desc(eta_id)%order => NULL()
    Particles%ops%desc(eta_id)%coeffs => NULL()
    Particles%ops%desc(eta_id)%name = particles_dflt_opname(eta_id)
    Particles%ops%desc(eta_id)%interp = .FALSE.
    Particles%ops%desc(eta_id)%nterms = 0
    Particles%ops%desc(eta_id)%vector = .FALSE.
    Particles%ops%ker(eta_id)%vec=> NULL()

    !update indices
    Particles%ops%nb_ops = Particles%ops%nb_ops - 1
    DO WHILE (.NOT.ASSOCIATED(Particles%ops%desc(Particles%ops%max_opsid)%degree))
        Particles%ops%max_opsid = Particles%ops%max_opsid - 1
        IF (Particles%ops%max_opsid .LE. 0) EXIT
    ENDDO

    !-------------------------------------------------------------------------
    ! Finalize
    !-------------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999  CONTINUE ! jump here upon error


END SUBROUTINE particles_dcop_free

SUBROUTINE particles_dcop_apply(Particles,from_id,to_id,eta_id,info)
    !!!------------------------------------------------------------------------!
    !!! NEW version
    !!! Apply DC kernel stored in eta_id to the scalar property stored
    !!! prop_from_id and store the results in prop_to_id
    !!!------------------------------------------------------------------------!
    USE ppm_module_data, ONLY: ppm_rank

#ifdef __MPI
    INCLUDE "mpif.h"
#endif

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
    INTEGER,                            INTENT(  OUT)   :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    ! local variables
    !-------------------------------------------------------------------------
    CHARACTER(LEN = ppm_char)                  :: cbuf,filename
    CHARACTER(LEN = ppm_char)               :: caller = 'particles_dcop_apply'
    INTEGER                                    :: ip,iq,ineigh,lda
    REAL(KIND(1.D0))                           :: t0
    REAL(MK),DIMENSION(:,:),POINTER            :: eta => NULL()
    REAL(MK),DIMENSION(:),  POINTER            :: wps1 => NULL(),wps2=>NULL()
    REAL(MK),DIMENSION(:,:),POINTER            :: wpv1 => NULL(),wpv2=>NULL()
    REAL(MK),DIMENSION(:),  POINTER            :: dwps => NULL()
    REAL(MK),DIMENSION(:,:),POINTER            :: dwpv => NULL()
    INTEGER, DIMENSION(:),  POINTER            :: nvlist => NULL()
    INTEGER, DIMENSION(:,:),POINTER            :: vlist => NULL()
    REAL(MK)                                   :: sig
    LOGICAL                                    :: vector

    !-------------------------------------------------------------------------
    ! Initialize
    !-------------------------------------------------------------------------
    info = 0 ! change if error occurs
    CALL substart(caller,t0,info)

    IF (.NOT. ASSOCIATED(Particles%ops)) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            & 'No operator data structure found',&
            __LINE__,info)
        GOTO 9999
    ENDIF
    IF (eta_id.LE.0 .OR. eta_id .GT. Particles%ops%max_opsid) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            & 'invalid value for ops id',&
            __LINE__,info)
        GOTO 9999
    ENDIF
    vector =  Particles%ops%desc(eta_id)%vector
    IF (to_id.EQ.0) THEN
        IF (vector) THEN
            CALL particles_allocate_wpv(Particles,to_id,&
                Particles%ops%desc(eta_id)%nterms,info,name="dflt_dcop_apply")
            IF (info .NE. 0) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_alloc,caller,&
                    'particles_allocate_wps failed',__LINE__,info)
                GOTO 9999
            ENDIF
        ELSE
            CALL particles_allocate_wps(Particles,to_id,&
                info,name="dflt_dcop_apply")
            IF (info .NE. 0) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_alloc,caller,&
                    'particles_allocate_wps failed',__LINE__,info)
                GOTO 9999
            ENDIF
        ENDIF
    ENDIF


    IF (vector) THEN
        dwpv => Get_wpv(Particles,to_id)
        lda = Particles%ops%desc(eta_id)%nterms
        DO ip = 1,Particles%Npart 
            dwpv(1:lda,ip) = 0._MK
        ENDDO
    ELSE
        dwps => Get_wps(Particles,to_id)
        DO ip = 1,Particles%Npart 
            dwps(ip) = 0._MK
        ENDDO
    ENDIF
    eta => Get_dcop(Particles,eta_id)


    IF (Particles%ops%desc(eta_id)%interp) THEN
        IF (.NOT. ASSOCIATED(Particles%Particles_cross)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,&
                'Need to specify which set of particles &
                &   (particles_cross) should be used for interpolation',&
                & __LINE__,info)
            GOTO 9999
        ENDIF
        IF (.NOT. (Particles%neighlists_cross)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,&
                'Please compute xset neighbor lists first',&
                __LINE__,info)
            GOTO 9999
        ENDIF
        nvlist => Particles%nvlist_cross
        vlist => Particles%vlist_cross
        IF (vector) THEN
            wpv2 => Get_wpv(Particles%Particles_cross,from_id,with_ghosts=.TRUE.)
            DO ip = 1,Particles%Npart
                DO ineigh = 1,nvlist(ip)
                    iq = vlist(ineigh,ip)
                    dwpv(1:lda,ip) = dwpv(1:lda,ip) + &
                        wpv2(1:lda,iq) * eta(1+(ineigh-1)*lda:ineigh*lda,ip)
                ENDDO
            ENDDO
            wpv2 => Set_wpv(Particles%Particles_cross,from_id,read_only=.TRUE.)
        ELSE
            wps2 => Get_wps(Particles%Particles_cross,from_id,with_ghosts=.TRUE.)
            DO ip = 1,Particles%Npart
                DO ineigh = 1,nvlist(ip)
                    iq = vlist(ineigh,ip)
                    dwps(ip) = dwps(ip) + wps2(iq) * eta(ineigh,ip)
                ENDDO
            ENDDO
            wps2 => Set_wps(Particles%Particles_cross,from_id,read_only=.TRUE.)
        ENDIF
    ELSE
        nvlist => Particles%nvlist
        vlist => Particles%vlist
        sig = -1._mk !TODO FIXME
        IF (vector) THEN
            wpv1 => Get_wpv(Particles%Particles_cross,from_id,with_ghosts=.TRUE.)
            DO ip = 1,Particles%Npart
                DO ineigh = 1,nvlist(ip)
                    iq = vlist(ineigh,ip)
                    dwpv(1:lda,ip) = dwpv(1:lda,ip) + &
                        (wpv1(1:lda,iq) + sig*(wpv1(1:lda,ip)))* &
                        eta(1+(ineigh-1)*lda:ineigh*lda,ip)
                ENDDO
            ENDDO
            wpv1 => Set_wpv(Particles%Particles_cross,from_id,read_only=.TRUE.)
        ELSE
            wps1 => Get_wps(Particles,from_id,with_ghosts=.TRUE.)
            DO ip = 1,Particles%Npart
                DO ineigh = 1,nvlist(ip)
                    iq = vlist(ineigh,ip)
                    dwps(ip) = dwps(ip) + (wps1(iq)+sig*(wps1(ip))) * eta(ineigh,ip)
                ENDDO
            ENDDO
            wps1 => Set_wps(Particles,from_id,read_only=.TRUE.)
        ENDIF
    ENDIF

    eta => Set_dcop(Particles,eta_id)
    IF (vector) THEN
        dwpv => Set_wpv(Particles,to_id)
    ELSE
        dwps => Set_wps(Particles,to_id)
    ENDIF
    nvlist => NULL()
    vlist => NULL()

    CALL substop(caller,t0,info)

    9999  CONTINUE ! jump here upon error


END SUBROUTINE particles_dcop_apply

SUBROUTINE particles_apply_dcops(Particles,from_id,to_id,eta_id,sig,&
        info,Particles_old,from_old_id)
    !!!------------------------------------------------------------------------!
    !!! OLD version
    !!! Apply DC kernel stored in eta_id to the scalar property stored
    !!! prop_from_id and store the results in prop_to_id
    !!!------------------------------------------------------------------------!
    USE ppm_module_data, ONLY: ppm_rank

#ifdef __MPI
    INCLUDE "mpif.h"
#endif

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
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,&
                'particles_allocate_wps failed',__LINE__,info)
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

#undef __KIND
