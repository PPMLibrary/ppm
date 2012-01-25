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
    DO i = 1,Pc%max_wpid
        Pc%props(i)%t%flags(ppm_ppt_ghosts) = .FALSE.
        Pc%props(i)%t%flags(ppm_ppt_partial) = .FALSE.
    ENDDO
    DO i = 1,Pc%max_nlid
        Pc%neighs(i)%t%uptodate = .FALSE.
    ENDDO

    xp => NULL()

END SUBROUTINE DTYPE(set_xp)

SUBROUTINE DTYPE(part_prop_create)(Pc,propid,datatype,info,&
        lda,name,zero,with_ghosts)
    !!! Adds a property to an existing particle set
    CLASS(DTYPE(ppm_t_particles))         :: Pc
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
        IF (Pc%nwp.LT.Pc%swp) THEN
            !there is at least one empty slot in the array
            ! of property pointers
            propid = propid + 1
            DO WHILE (ASSOCIATED(Pc%props(propid)%t))
                propid = propid + 1

!--------not a necessary check-------
                IF (propid .GT. Pc%swp) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_alloc,caller,&
         &    'number of properties greater than allocated size',&
                        __LINE__,info)
                    GOTO 9999
                ENDIF
!------------------------------------

            ENDDO
        ELSE
            IF (.NOT. ASSOCIATED(Pc%props)) THEN
            !need to allocate the array of property pointers 
                nprops=20
                ALLOCATE(Pc%props(nprops),STAT=info)
                IF (info .NE. 0) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_alloc,caller,&
                        'allocating property array failed',__LINE__,info)
                    GOTO 9999
                ENDIF
                propid = 1
            ELSE
            !need to resize the array of property pointers 
                nprops=MAX(2*Pc%swp,20)
                ALLOCATE(prop_tmp(nprops),STAT=info)
                IF (info .NE. 0) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_alloc,caller,&
                        'allocating property array failed',__LINE__,info)
                    GOTO 9999
                ENDIF
                DO i=1,Pc%swp
                    prop_tmp(i)%t => Pc%props(i)%t
                ENDDO
                DEALLOCATE(Pc%props)
                Pc%props => prop_tmp
            ENDIF
            Pc%swp = nprops
            propid = Pc%nwp+1
        ENDIF
    ELSE
        !using a given propid
        !check that this id is not already used
        IF (ASSOCIATED(Pc%props(propid)%t)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,&
                'property id already in use. Use prop_destroy() first',&
                __LINE__,info)
            GOTO 9999
        ENDIF
    ENDIF
    Pc%nwp = Pc%nwp + 1
        

    IF (propid .GT. Pc%max_wpid) Pc%max_wpid = propid

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

    npart = Pc%Npart
    flags = .FALSE.
    IF (PRESENT(with_ghosts)) THEN
        IF (with_ghosts) THEN
            npart = Pc%Mpart
            flags(ppm_ppt_ghosts) = .TRUE.
        ENDIF
    ENDIF
    flags(ppm_ppt_partial) = .TRUE.

    IF (.NOT. ASSOCIATED(Pc%props(propid)%t)) THEN
        ALLOCATE(Pc%props(propid)%t,STAT=info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,&
                'allocating property pointer failed',__LINE__,info)
            GOTO 9999
        ENDIF
    ENDIF

    ! Create the property
    CALL Pc%props(propid)%t%create(&
       &                datatype,npart,lda2,name2,flags,info,zero)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_sub_failed,caller,&
            'creating property array failed',__LINE__,info)
        GOTO 9999
    ENDIF

    CALL substop(caller,t0,info)
    9999  CONTINUE

END SUBROUTINE DTYPE(part_prop_create)

SUBROUTINE DTYPE(part_prop_destroy)(Pc,propid,info)
    !!! Destroy a property from an existing particle set
    CLASS(DTYPE(ppm_t_particles))         :: Pc
    INTEGER,                INTENT(INOUT) :: propid
    INTEGER,               INTENT(OUT)    :: info

    CHARACTER(LEN=ppm_char)               :: caller = 'particle_prop_destroy'
    REAL(KIND(1.D0))                      :: t0

    CALL substart(caller,t0,info)

    IF (propid .LE. 0 .OR. propid .GT. Pc%swp) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,&
            &    'property id larger than size of properties array',&
            __LINE__,info)
        GOTO 9999
    ENDIF

    CALL Pc%props(propid)%t%destroy(info)
    NULLIFY(Pc%props(propid)%t)

    Pc%nwp = Pc%nwp - 1
    IF (propid .EQ. Pc%max_wpid) THEN
        Pc%max_wpid = Pc%max_wpid - 1
        IF (Pc%max_wpid .GT. 0) THEN
            DO WHILE(.NOT.ASSOCIATED(Pc%props(Pc%max_wpid)%t))
                Pc%max_wpid = Pc%max_wpid - 1
                IF (Pc%max_wpid .EQ. 0) EXIT
            ENDDO
        ENDIF
    ENDIF

    CALL substop(caller,t0,info)

    9999  CONTINUE

END SUBROUTINE DTYPE(part_prop_destroy)

SUBROUTINE DTYPE(prop_create)(prop,datatype,npart,lda,name,flags,info,zero)
    !!! Constructor for particle property data structure
    DEFINE_MK()
    CLASS(DTYPE(part_prop))            :: prop
    INTEGER,                INTENT(IN) :: datatype
    INTEGER,                INTENT(IN) :: npart
    INTEGER,                INTENT(IN) :: lda
    CHARACTER(LEN=*)                   :: name
    !!! name to this property
    LOGICAL, DIMENSION(ppm_param_length_pptflags),INTENT(IN) :: flags
    INTEGER,               INTENT(OUT) :: info
    LOGICAL, OPTIONAL,     INTENT( IN) :: zero
    !!! if true, then initialize the data to zero

    REAL(KIND(1.D0))                   :: t0
    INTEGER                            :: iopt
    CHARACTER(LEN=ppm_char)            :: caller = 'prop_create'
    LOGICAL                            :: is2d
    LOGICAL                            :: zero_data



    CALL substart(caller,t0,info)

    prop%lda       = lda
    prop%data_type = datatype
    prop%name      = name
    prop%flags     = flags

    IF (PRESENT(zero)) THEN
        zero_data = zero
    ELSE
        zero_data = .FALSE.
    ENDIF


    iopt   = ppm_param_alloc_fit

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
#if   __KIND == __SINGLE_PRECISION
        CASE (ppm_type_real_single )
            CALL ppm_alloc(prop%data_2d_r,ldc,iopt,info)
            IF (zero_data) prop%data_2d_r(1:lda,1:npart) = 0._MK
        CASE (ppm_type_comp_single)
            CALL ppm_alloc(prop%data_2d_c,ldc,iopt,info)
            IF (zero_data) prop%data_2d_c(1:lda,1:npart) = 0._MK
#elif __KIND ==__DOUBLE_PRECISION
        CASE (ppm_type_real_double)
            CALL ppm_alloc(prop%data_2d_r,ldc,iopt,info)
            IF (zero_data) prop%data_2d_r(1:lda,1:npart) = 0._MK
        CASE (ppm_type_comp_double)
            CALL ppm_alloc(prop%data_2d_c,ldc,iopt,info)
            IF (zero_data) prop%data_2d_c(1:lda,1:npart) = 0._MK
#endif
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
#if   __KIND == __SINGLE_PRECISION
        CASE (ppm_type_real_single )
            CALL ppm_alloc(prop%data_1d_r,ldc,iopt,info)
            IF (zero_data) prop%data_1d_r(1:npart) = 0._MK
        CASE (ppm_type_comp_single)
            CALL ppm_alloc(prop%data_1d_c,ldc,iopt,info)
            IF (zero_data) prop%data_1d_c(1:npart) = 0._MK
#elif __KIND ==__DOUBLE_PRECISION
        CASE (ppm_type_real_double)
            CALL ppm_alloc(prop%data_1d_r,ldc,iopt,info)
            IF (zero_data) prop%data_1d_r(1:npart) = 0._MK
        CASE (ppm_type_comp_double)
            CALL ppm_alloc(prop%data_1d_c,ldc,iopt,info)
            IF (zero_data) prop%data_1d_c(1:npart) = 0._MK
#endif
        CASE (ppm_type_logical )
            CALL ppm_alloc(prop%data_1d_l,ldc,iopt,info)
            IF (zero_data) prop%data_1d_l(1:npart) = .FALSE.
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

SUBROUTINE DTYPE(prop_destroy)(prop,info)
    CLASS(DTYPE(part_prop))      :: prop
    INTEGER,                                INTENT(  OUT)  :: info
    !!! Returns status, 0 upon success.

    IF(ASSOCIATED(prop%data_1d_i)) DEALLOCATE(prop%data_1d_i,STAT=info)
    IF(ASSOCIATED(prop%data_2d_i)) DEALLOCATE(prop%data_2d_i,STAT=info)
    IF(ASSOCIATED(prop%data_1d_li)) DEALLOCATE(prop%data_1d_li,STAT=info)
    IF(ASSOCIATED(prop%data_2d_li)) DEALLOCATE(prop%data_2d_li,STAT=info)
    IF(ASSOCIATED(prop%data_1d_r)) DEALLOCATE(prop%data_1d_r,STAT=info)
    IF(ASSOCIATED(prop%data_2d_r)) DEALLOCATE(prop%data_2d_r,STAT=info)
    IF(ASSOCIATED(prop%data_1d_c)) DEALLOCATE(prop%data_1d_c,STAT=info)
    IF(ASSOCIATED(prop%data_2d_c)) DEALLOCATE(prop%data_2d_c,STAT=info)
    IF(ASSOCIATED(prop%data_1d_l)) DEALLOCATE(prop%data_1d_l,STAT=info)
    IF(ASSOCIATED(prop%data_2d_l)) DEALLOCATE(prop%data_2d_l,STAT=info)

    prop%data_type = ppm_type_none
    prop%lda = 0
    prop%flags = .FALSE.

END SUBROUTINE DTYPE(prop_destroy)

SUBROUTINE DTYPE(neigh_destroy)(neigh,info)
    CLASS(DTYPE(ppm_t_neighlist))      :: neigh
    INTEGER,                                INTENT(  OUT)  :: info
    !!! Returns status, 0 upon success.

    IF(ASSOCIATED(neigh%nvlist)) DEALLOCATE(neigh%nvlist,STAT=info)
    IF(ASSOCIATED(neigh%vlist))  DEALLOCATE(neigh%vlist,STAT=info)

END SUBROUTINE DTYPE(neigh_destroy)

SUBROUTINE DTYPE(op_destroy)(op,info)
    CLASS(DTYPE(ppm_t_operator))              :: op
    INTEGER                                   :: i
    INTEGER,                                INTENT(  OUT)  :: info
    !!! Returns status, 0 upon success.
    
    DO i=1,op%max_opsid
        CALL op%ker(i)%t%destroy(info)
        CALL op%desc(i)%t%destroy(info)
    ENDDO
    op%max_opsid = 0
    op%nb_ops = 0

END SUBROUTINE DTYPE(op_destroy)

SUBROUTINE DTYPE(desc_destroy)(desc,info)
    CLASS(DTYPE(ppm_t_opdesc))              :: desc
    INTEGER,                                INTENT(  OUT)  :: info
    !!! Returns status, 0 upon success.
    
    IF (ASSOCIATED(desc%degree)) DEALLOCATE(desc%degree,STAT=info)
    IF (ASSOCIATED(desc%order))  DEALLOCATE(desc%order,STAT=info)
    IF (ASSOCIATED(desc%coeffs)) DEALLOCATE(desc%coeffs,STAT=info)

END SUBROUTINE DTYPE(desc_destroy)

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
    REAL(KIND(1.D0))                                :: t0
    CHARACTER(LEN = ppm_char)            :: caller = 'ppm_create_particles'
    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)

    !-----------------------------------------------------------------
    !  Destroy the DS if it already exists
    !-----------------------------------------------------------------
    IF (ASSOCIATED(Pc%xp)) THEN
        CALL Pc%destroy(info)
    ENDIF
    !-----------------------------------------------------------------
    !  Allocate memory for the positions
    !-----------------------------------------------------------------
    ldc(1) = ppm_dim
    ldc(2) = Npart
    CALL ppm_alloc(Pc%xp,ldc(1:2),ppm_param_alloc_fit,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &        'Could not allocate Particles elements',__LINE__,info)
        GOTO 9999
    ENDIF
    Pc%Npart = Npart
    Pc%Mpart = Npart
    Pc%flags(ppm_part_ghosts) = .FALSE.
    Pc%flags(ppm_part_areinside) = .FALSE.
    Pc%flags(ppm_part_partial) = .FALSE.
    Pc%flags(ppm_part_reqput) = .FALSE.
    Pc%flags(ppm_part_cartesian) = .FALSE.
    Pc%active_topoid = -1
    ! No active topology yet

    ! Give a default name to this Particle set
    IF (PRESENT(name)) THEN
        Pc%name = ADJUSTL(TRIM(name))
    ELSE
        Pc%name = particles_dflt_partname()
    ENDIF

    ! No properties defined
    Pc%nwp = 0
    Pc%swp = 0
    Pc%max_wpid = 0
    Pc%props => NULL()

    ! No neighbor lists defined
    Pc%nnl = 0
    Pc%snl = 0
    Pc%max_nlid = 0
    Pc%neighs => NULL()

    ! Particles are by default not adaptive
    Pc%adaptive = .FALSE.
    Pc%adapt_wpid = 0
    Pc%gi_id = 0
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
    ! Particles are by default not a Cartesian grid
    Pc%cartesian = .FALSE.
    ! Particles have not been initialised yet
    Pc%h_avg = -1._MK
    Pc%h_min = -1._MK

    Pc%time = 0._MK
    Pc%itime = 0


    !-------------------------------------------------------------------------
    !  Finalize
    !-------------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999 CONTINUE

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
    ! first deallocate all content of Pc
    IF (ASSOCIATED(Pc%xp)) THEN 
        DEALLOCATE(Pc%xp,STAT=info)
        NULLIFY(Pc%xp)
    ENDIF
    IF (ASSOCIATED(Pc%pcost)) THEN
        DEALLOCATE(Pc%pcost,STAT=info)
        NULLIFY(Pc%pcost)
    ENDIF

    !Deallocate neighbour lists
    IF (ASSOCIATED(Pc%neighs)) THEN
        DO i=1,Pc%max_nlid
            CALL Pc%neighs(i)%t%destroy(info)
        ENDDO
        DEALLOCATE(Pc%neighs)
        NULLIFY(Pc%neighs)
    ENDIF

    !Deallocate properties
    IF (ASSOCIATED(Pc%props)) THEN
        DO i=1,Pc%max_wpid
            CALL Pc%props(i)%t%destroy(info)
        ENDDO
        DEALLOCATE(Pc%props)
        NULLIFY(Pc%props)
    ENDIF

    IF (ASSOCIATED(Pc%ops)) THEN
        CALL Pc%ops%destroy(info)
        DEALLOCATE(Pc%ops)
        NULLIFY(Pc%ops)
    ENDIF

    !-------------------------------------------------------------------------
    !  Finalize
    !-------------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999 CONTINUE

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
    USE ppm_module_data, ONLY: ppm_rank,ppm_nproc,ppm_topo

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

    IF (ppm_dim .eq. 2) THEN
        CALL  DTYPE(particles_initialize2d)(Pc,Npart_global,info,&
            distrib,topoid,minphys,maxphys,cutoff,name=name)
    ELSE
        CALL  DTYPE(particles_initialize3d)(Pc,Npart_global,info,&
            distrib,topoid,minphys,maxphys,cutoff,name=name)
    ENDIF
END SUBROUTINE DTYPE(part_initialize)

#undef DEFINE_MK


