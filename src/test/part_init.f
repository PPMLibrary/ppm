#if   __KIND == __SINGLE_PRECISION
SUBROUTINE part_inits(xp,Npart_global,min_phys,max_phys,info,&
           &distrib,grid_shift)
#elif __KIND == __DOUBLE_PRECISION
SUBROUTINE part_initd(xp,Npart_global,min_phys,max_phys,info,&
           &distrib,grid_shift)
#endif
    !-----------------------------------------------------------------------
    ! Set initial particle positions
    !-----------------------------------------------------------------------
    USE ppm_module_substart
    USE ppm_module_substop
    USE ppm_module_typedef
    USE ppm_module_data, ONLY: ppm_dim,ppm_rank,ppm_nproc,ppm_topo,ppm_comm
    USE ppm_module_write
    USE ppm_module_alloc
    USE ppm_module_error
#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    REAL(MK),DIMENSION(:,:), POINTER                       :: xp
    !!! Data structure containing the particles
    REAL(MK),DIMENSION(ppm_dim),         INTENT(IN   )     :: min_phys
    !!! extent of the physical domain.
    REAL(MK),DIMENSION(ppm_dim),         INTENT(IN   )     :: max_phys
    !!! extent of the physical domain.
    INTEGER,                            INTENT(INOUT)      :: Npart_global
    !!! total number of particles that will be initialized
    INTEGER,                            INTENT(  OUT)      :: info
    !!! Return status, on success 0.
    !-------------------------------------------------------------------------
    !  Optional arguments
    !-------------------------------------------------------------------------
    INTEGER,OPTIONAL,                   INTENT(IN   )      :: distrib
    !!! type of initial distribution. One of
    !!! * ppm_param_part_init_cartesian (default)
    !!! * ppm_param_part_init_random
    REAL(MK),OPTIONAL,                   INTENT(IN   )     :: grid_shift
    !!! shifts positions of the particles. Set to 0.5 to place
    !!! particles in the middle of each cell, set to 0 to place 
    !!! them in the lower left corner


    !-------------------------------------------------------------------------
    !  Local variables
    !------------------------------------------------------------------------- 
    INTEGER                               :: ip,i,j,k,Npart,iopt
    INTEGER, DIMENSION(3)                 :: nijk
    INTEGER, DIMENSION(3)                 :: nijk_global
    CHARACTER(LEN = ppm_char)             :: caller = 'part_init'
    REAL(MK)                              :: h
    REAL(MK)                              :: shift
    REAL(MK),DIMENSION(3)                 :: pos
    REAL(KIND(1.D0))                      :: t0
    INTEGER                               :: remaining_rows

    INTEGER                               :: distribution
    INTEGER                               :: part_seedsize
    INTEGER,  DIMENSION(3)                :: ldc
    INTEGER,  DIMENSION(:), POINTER       :: part_seed => NULL()

    REAL(MK), DIMENSION(:  ), POINTER     :: randnb => NULL()
    REAL(MK),DIMENSION(3)                 :: minphys
    REAL(MK),DIMENSION(3)                 :: maxphys
    REAL(MK),DIMENSION(3)                 :: lenphys


    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    CALL substart(caller,t0,info)

    IF(PRESENT(distrib)) THEN
        distribution=distrib
    ELSE
        distribution=ppm_param_part_init_cartesian
    ENDIF
    IF (PRESENT(grid_shift)) THEN
        shift = grid_shift
    ELSE
        shift = 0.0_mk
    ENDIF
    minphys = 0.0_mk
    maxphys = 0.0_mk
    minphys(1:ppm_dim) = min_phys(1:ppm_dim)
    maxphys(1:ppm_dim) = max_phys(1:ppm_dim)
    lenphys = maxphys - minphys


    h = (PRODUCT(lenphys(1:ppm_dim))/REAL(Npart_global))**(1./REAL(ppm_dim))
    nijk_global = 1
    nijk_global(1:ppm_dim) = FLOOR(lenphys(1:ppm_dim)/h)
    Npart_global = PRODUCT(nijk_global(1:ppm_dim))
    remaining_rows = MOD(nijk_global(ppm_dim),ppm_nproc)
    !number of particles along x and z
    nijk = 1
    nijk(1:ppm_dim) = nijk_global(1:ppm_dim)
    !number of particles along y 
    nijk(2) = nijk_global(ppm_dim)/ppm_nproc

    !number of particles on this processor
    Npart = PRODUCT(nijk)

    !proc 0 takes care of the additional rows (remainder)
    IF (ppm_rank.EQ.0) THEN
        Npart = Npart + remaining_rows * nijk(1)*nijk(3)
    ENDIF
    iopt = ppm_param_alloc_fit
    ldc(1) = ppm_dim
    ldc(2) = Npart
    CALL ppm_alloc(xp,ldc,iopt,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,&
            'ppm_alloc_particles (allocate) failed',__LINE__,info)
        GOTO 9999
    ENDIF

    !-----------------------------------------------------------------------
    ! set particles
    !-----------------------------------------------------------------------
    ip = 0
    if_cartesian: IF (distribution .EQ. ppm_param_part_init_cartesian) THEN
        DO k = 1,nijk(3)
            h = lenphys(3)/REAL(nijk(3),MK)
            pos(3) = minphys(3) + h*(k-1) + shift*h
        DO j = 1,nijk(2)
            h = lenphys(2)/REAL(nijk_global(2),MK)
            pos(2) = minphys(2) + h*(j-1 + ppm_rank*nijk(2)) + shift*h
            DO i = 1,nijk(1)
                h = lenphys(1)/REAL(nijk(1),MK)
                pos(1) = minphys(1) + h*(i-1) + shift*h

                ip = ip + 1
                xp(1:ppm_dim,ip) = pos(1:ppm_dim)

                ! impose periodic boundaries:
                IF (xp(1,ip) .GE. maxphys(1)) xp(1,ip) = xp(1,ip) - lenphys(1)
                IF (xp(2,ip) .GE. maxphys(2)) xp(2,ip) = xp(2,ip) - lenphys(2)
                IF (xp(1,ip) .LT. minphys(1)) xp(1,ip) = xp(1,ip) + lenphys(1)
                IF (xp(2,ip) .LT. minphys(2)) xp(2,ip) = xp(2,ip) + lenphys(2)
            ENDDO
        ENDDO
        ENDDO
        IF (ppm_dim.EQ.3) THEN
            DO ip=1,Npart
                IF (xp(3,ip) .GE. maxphys(3)) xp(3,ip) = xp(3,ip) - lenphys(3)
                IF (xp(3,ip) .LT. minphys(3)) xp(3,ip) = xp(3,ip) + lenphys(3)
            ENDDO
        ENDIF


        IF(ppm_rank.EQ.0) THEN
        DO k = 1,nijk(3)
            h = lenphys(3)/REAL(nijk(3),MK)
            pos(3) = minphys(3) + h*(k-1) + shift*h
            DO j = 1,remaining_rows
                h = lenphys(2)/REAL(nijk_global(2),MK)
                pos(2) = minphys(2) + h*(j-1 + ppm_nproc*nijk(2)) + shift*h
                DO i = 1,nijk(1)
                    h = lenphys(1)/REAL(nijk(1),MK)
                    pos(1) = minphys(1) + h*(i-1) + shift*h

                    ip = ip + 1
                    xp(1:ppm_dim,ip) = pos(1:ppm_dim)

                    ! impose periodic boundaries:
                    IF (xp(1,ip) .GE. maxphys(1)) xp(1,ip) = xp(1,ip) - lenphys(1)
                    IF (xp(2,ip) .GE. maxphys(2)) xp(2,ip) = xp(2,ip) - lenphys(2)
                    IF (xp(1,ip) .LT. minphys(1)) xp(1,ip) = xp(1,ip) + lenphys(1)
                    IF (xp(2,ip) .LT. minphys(2)) xp(2,ip) = xp(2,ip) + lenphys(2)
                ENDDO
            ENDDO
        ENDDO
        IF (ppm_dim.EQ.3) THEN
            DO ip=1,Npart
                IF (xp(3,ip) .GE. maxphys(3)) xp(3,ip) = xp(3,ip) - lenphys(3)
                IF (xp(3,ip) .LT. minphys(3)) xp(3,ip) = xp(3,ip) + lenphys(3)
            ENDDO
        ENDIF
        ENDIF
    ELSE !random distribution
        iopt = ppm_param_alloc_fit
        ldc(1) = ppm_dim*Npart
        CALL ppm_alloc(randnb,ldc,iopt,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,   &
                &            'allocation failed',__LINE__,info)
            GOTO 9999
        ENDIF
        IF (.NOT.ASSOCIATED(part_seed)) THEN
            CALL RANDOM_SEED(SIZE=part_seedsize)
            ldc(1) = part_seedsize
            CALL ppm_alloc(part_seed,ldc,iopt,info)
            IF (info .NE. 0) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_alloc,caller,   &
                    &            'allocation failed',__LINE__,info)
                GOTO 9999
            ENDIF
            DO i=1,part_seedsize
                part_seed(i)=i*i*i*i
            ENDDO
            CALL RANDOM_SEED(PUT=part_seed)
        ENDIF
        CALL RANDOM_NUMBER(randnb)
        ip = 1
        DO k = 1,nijk(3)
            h = lenphys(3)/REAL(nijk(3),MK)
            pos(3) = minphys(3) + &
            &        ((k-1) + shift+randnb(ppm_dim*ip - 2))*h
        DO j = 1,nijk(2)
            h = lenphys(2)/REAL(nijk_global(2),MK)
            pos(2) = minphys(2) + &
            &        ((j-1 + ppm_rank*nijk(2)) + shift + randnb(ppm_dim*ip))*h
            DO i = 1,nijk(1)
                h = lenphys(1)/REAL(nijk(1),MK)
                pos(1) = minphys(1) + &
                &        ((i-1) + shift + randnb(ppm_dim*ip - 1))*h

                ip = ip + 1
                ! uniformly random in cells
                xp(1:ppm_dim,ip) = pos(1:ppm_dim)
                ! impose periodic boundaries:
                IF (xp(1,ip) .GE. maxphys(1)) xp(1,ip) = xp(1,ip) - lenphys(1)
                IF (xp(2,ip) .GE. maxphys(2)) xp(2,ip) = xp(2,ip) - lenphys(2)
                IF (xp(1,ip) .LT. minphys(1)) xp(1,ip) = xp(1,ip) + lenphys(1)
                IF (xp(2,ip) .LT. minphys(2)) xp(2,ip) = xp(2,ip) + lenphys(2)
            ENDDO
        ENDDO
        ENDDO
        IF (ppm_dim.EQ.3) THEN
            DO ip=1,Npart
                IF (xp(3,ip) .GE. maxphys(3)) xp(3,ip) = xp(3,ip) - lenphys(3)
                IF (xp(3,ip) .LT. minphys(3)) xp(3,ip) = xp(3,ip) + lenphys(3)
            ENDDO
        ENDIF
        ip = 1
        IF(ppm_rank.EQ.0) THEN
        DO k = 1,nijk(3)
            h = lenphys(3)/REAL(nijk(3),MK)
            pos(3) = minphys(3) + &
            &        ((k-1) + shift + randnb(ppm_dim*ip - 2))*h
            DO j = 1,remaining_rows
                h = lenphys(2)/REAL(nijk_global(2),MK)
                pos(2) = minphys(2) + &
                &        h*((j-1 + ppm_nproc*nijk(2)) + shift + randnb(ppm_dim*ip))*h
                DO i = 1,nijk(1)
                    h = lenphys(1)/REAL(nijk(1),MK)
                    pos(1) = minphys(1) + &
                    &        ((i-1) + shift + randnb(ppm_dim*ip - 1))*h

                    ip = ip + 1
                    ! uniformly random in cells
                    xp(1:ppm_dim,ip) = pos(1:ppm_dim)
                       
                        
                    ! impose periodic boundaries:
                    IF (xp(1,ip) .GE. maxphys(1)) xp(1,ip) = xp(1,ip) - lenphys(1)
                    IF (xp(2,ip) .GE. maxphys(2)) xp(2,ip) = xp(2,ip) - lenphys(2)
                    IF (xp(1,ip) .LT. minphys(1)) xp(1,ip) = xp(1,ip) + lenphys(1)
                    IF (xp(2,ip) .LT. minphys(2)) xp(2,ip) = xp(2,ip) + lenphys(2)
#if __DIM == __3D
                    IF (xp(3,ip) .GE. maxphys(3)) xp(3,ip) = xp(3,ip) - lenphys(3)
                    IF (xp(3,ip) .LT. minphys(3)) xp(3,ip) = xp(3,ip) + lenphys(3)
#endif
                ENDDO
            ENDDO
        ENDDO
        IF (ppm_dim.EQ.3) THEN
            DO ip=1,Npart
                IF (xp(3,ip) .GE. maxphys(3)) xp(3,ip) = xp(3,ip) - lenphys(3)
                IF (xp(3,ip) .LT. minphys(3)) xp(3,ip) = xp(3,ip) + lenphys(3)
            ENDDO
        ENDIF
        ENDIF

        DEALLOCATE(randnb)

    ENDIF if_cartesian

    Npart_global = Npart

    !-----------------------------------------------------------------------
    ! Finalize
    !-----------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999 CONTINUE ! jump here upon error

#if   __KIND == __SINGLE_PRECISION
END SUBROUTINE part_inits
#elif __KIND == __DOUBLE_PRECISION
END SUBROUTINE part_initd
#endif

