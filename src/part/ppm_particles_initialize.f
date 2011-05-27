#if __DIM == 2
SUBROUTINE particles_initialize2d(Particles,Npart_global,info,&
        distrib,topoid,minphys,maxphys,cutoff)
#elif __DIM == 3
SUBROUTINE particles_initialize3d(Particles,Npart_global,info,&
        distrib,topoid,minphys,maxphys,cutoff)
#endif
    !-----------------------------------------------------------------------
    ! Set initial particle positions
    !-----------------------------------------------------------------------
    USE ppm_module_substart
    USE ppm_module_substop
    USE ppm_module_data, ONLY: ppm_rank,ppm_nproc,ppm_topo,ppm_comm
    USE ppm_module_write
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
    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------

    INTEGER                               :: ip,i,j,k,Npart,iopt
    INTEGER                               :: nijk(ppm_dim),nijk_global(ppm_dim)
    CHARACTER(LEN = ppm_char)              :: filename,cbuf
    CHARACTER(LEN = ppm_char)              :: caller = 'particles_initialize'
    REAL(MK)                              :: y,z,h
    REAL(KIND(1.D0))                      :: t0
    INTEGER                               :: remaining_rows

    REAL(MK)                              :: shift
    INTEGER                               :: distribution
    TYPE(ppm_t_topo),POINTER              :: topo => NULL()
    REAL(MK), DIMENSION(ppm_dim)          :: min_phys,max_phys,len_phys

    REAL(MK), DIMENSION(:,:), POINTER     :: xp
    REAL(MK), DIMENSION(:  ), POINTER     :: randnb
    INTEGER,  DIMENSION(:  ), POINTER     :: seed
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
            info = ppm_error_error
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
        info = ppm_error_error
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
    nijk(1:ppm_dim) = nijk_global(1:ppm_dim-1)
    !number of particles along y 
    nijk(2) = nijk_global(ppm_dim)/ppm_nproc

    !number of particles on this processor
    Npart = PRODUCT(nijk)

    !proc 0 takes care of the additional rows (remainder)
    IF (ppm_rank.EQ.0) THEN
#if   __DIM == 2
        Npart = Npart + remaining_rows * nijk(1)
#elif __DIM == 3
        Npart = Npart + remaining_rows * nijk(1)*nijk(3)
#endif
    ENDIF

    !Deallocate Particles if already allocated
    IF (ASSOCIATED(Particles)) THEN
        CALL ppm_alloc_particles(Particles,Npart,ppm_param_dealloc,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_dealloc,caller,&
                'ppm_alloc_particles (deallocate) failed',__LINE__,info)
            GOTO 9999
        ENDIF
    ENDIF

    CALL ppm_alloc_particles(Particles,Npart,ppm_param_alloc_fit,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,&
            'ppm_alloc_particles (allocate) failed',__LINE__,info)
        GOTO 9999
    ENDIF
    !use a shortcut, for convenience
    xp => Particles%xp

    !-----------------------------------------------------------------------
    ! set particles
    !-----------------------------------------------------------------------
    ip = 0
    shift = 0._MK !shifts positions of the particles. Set to 0.5 to place
    !particles in the middle of each cell, set to 0 to place them in the lower
    !left corner

    if_cartesian: IF (distribution .EQ. ppm_param_part_init_cartesian) THEN
#if __DIM == 3
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
#if __DIM == 3
                xp(3,ip) = z 
#endif

                ! impose periodic boundaries:
                IF (xp(1,ip) .GE. max_phys(1)) xp(1,ip) = xp(1,ip) - len_phys(1)
                IF (xp(2,ip) .GE. max_phys(2)) xp(2,ip) = xp(2,ip) - len_phys(2)
                IF (xp(1,ip) .LT. min_phys(1)) xp(1,ip) = xp(1,ip) + len_phys(1)
                IF (xp(2,ip) .LT. min_phys(2)) xp(2,ip) = xp(2,ip) + len_phys(2)
#if __DIM == 3
                IF (xp(3,ip) .GE. max_phys(3)) xp(3,ip) = xp(3,ip) - len_phys(3)
                IF (xp(3,ip) .LT. min_phys(3)) xp(3,ip) = xp(3,ip) + len_phys(3)
#endif

            ENDDO
        ENDDO
#if __DIM == 3
        ENDDO
#endif

        IF(ppm_rank.EQ.0) THEN
#if __DIM == 3
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
#if __DIM == 3
                    xp(3,ip) = z 
#endif

                    ! impose periodic boundaries:
                    IF (xp(1,ip) .GE. max_phys(1)) xp(1,ip) = xp(1,ip) - len_phys(1)
                    IF (xp(2,ip) .GE. max_phys(2)) xp(2,ip) = xp(2,ip) - len_phys(2)
                    IF (xp(1,ip) .LT. min_phys(1)) xp(1,ip) = xp(1,ip) + len_phys(1)
                    IF (xp(2,ip) .LT. min_phys(2)) xp(2,ip) = xp(2,ip) + len_phys(2)
#if __DIM == 3
                    IF (xp(3,ip) .GE. max_phys(3)) xp(3,ip) = xp(3,ip) - len_phys(3)
                    IF (xp(3,ip) .LT. min_phys(3)) xp(3,ip) = xp(3,ip) + len_phys(3)
#endif

                ENDDO
            ENDDO
#if __DIM == 3
        ENDDO
#endif
        ENDIF
        Particles%cartesian = .TRUE.
    ELSE !random distribution
        iopt = ppm_param_alloc_fit
#ifdef same_random_sequence_nproc
        ldc(1) = ppm_dim*Npart_global
#else
        ldc(1) = ppm_dim*Npart
#endif
        CALL ppm_alloc(randnb,ldc(1:1),iopt,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,   &
                &            'allocation failed',__LINE__,info)
            GOTO 9999
        ENDIF
        CALL RANDOM_SEED(SIZE=seedsize)
        ldc(1) = seedsize
        CALL ppm_alloc(seed,ldc(1:1),iopt,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,   &
                &            'allocation failed',__LINE__,info)
            GOTO 9999
        ENDIF
        DO i=1,seedsize
            seed(i)=i*i*i*i
        ENDDO
        CALL RANDOM_SEED(PUT=seed)
        CALL RANDOM_NUMBER(randnb)

#if __DIM == 3
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
#if __DIM == 3
                xp(3,ip) = z                   + &
                    randnb(ppm_dim*ppm_rank*PRODUCT(nijk)+ppm_dim*ip - 2)*h
#endif
#else
                xp(1,ip) = min_phys(1) + h*(i-1) + shift + &
                    randnb(ppm_dim*ip - 1)*h
                xp(2,ip) = y                   + &
                    randnb(ppm_dim*ip    )*h
#if __DIM == 3
                xp(3,ip) = z                   + &
                    randnb(ppm_dim*ip - 2)*h
#endif
#endif
                ! impose periodic boundaries:
                IF (xp(1,ip) .GE. max_phys(1)) xp(1,ip) = xp(1,ip) - len_phys(1)
                IF (xp(2,ip) .GE. max_phys(2)) xp(2,ip) = xp(2,ip) - len_phys(2)
                IF (xp(1,ip) .LT. min_phys(1)) xp(1,ip) = xp(1,ip) + len_phys(1)
                IF (xp(2,ip) .LT. min_phys(2)) xp(2,ip) = xp(2,ip) + len_phys(2)
#if __DIM == 3
                IF (xp(3,ip) .GE. max_phys(3)) xp(3,ip) = xp(3,ip) - len_phys(3)
                IF (xp(3,ip) .LT. min_phys(3)) xp(3,ip) = xp(3,ip) + len_phys(3)
#endif
            ENDDO
        ENDDO
#if __DIM == 3
        ENDDO
#endif
        IF(ppm_rank.EQ.0) THEN
#if __DIM == 3
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
#if __DIM == 3
                    xp(3,ip) = z                   + &
                        randnb(ppm_dim*(ppm_nproc-1)*&
                        PRODUCT(nijk)+ppm_dim*ip - 2)*h
#endif
#else
                    xp(1,ip) = min_phys(1) + h*(i-1) + shift + &
                        randnb(ppm_dim*ip - 1)*h
                    xp(2,ip) = y                   + &
                        randnb(ppm_dim*ip    )*h
#if __DIM == 3
                    xp(3,ip) = z                   + &
                        randnb(ppm_dim*ip - 2)*h
#endif
#endif
                    ! impose periodic boundaries:
                    IF (xp(1,ip) .GE. max_phys(1)) xp(1,ip) = xp(1,ip) - len_phys(1)
                    IF (xp(2,ip) .GE. max_phys(2)) xp(2,ip) = xp(2,ip) - len_phys(2)
                    IF (xp(1,ip) .LT. min_phys(1)) xp(1,ip) = xp(1,ip) + len_phys(1)
                    IF (xp(2,ip) .LT. min_phys(2)) xp(2,ip) = xp(2,ip) + len_phys(2)
#if __DIM == 3
                    IF (xp(3,ip) .GE. max_phys(3)) xp(3,ip) = xp(3,ip) - len_phys(3)
                    IF (xp(3,ip) .LT. min_phys(3)) xp(3,ip) = xp(3,ip) + len_phys(3)
#endif
                ENDDO
            ENDDO
#if __DIM == 3
        ENDDO
#endif
        ENDIF
        Particles%cartesian = .FALSE.

        DEALLOCATE(randnb,seed)

    ENDIF if_cartesian

    xp=>NULL()

    ! (global) average interparticle spacing
    Particles%h_avg = (PRODUCT(len_phys)/REAL(Npart_global))**(1./REAL(ppm_dim)) 
    ! min interparticle spacing (not needed now)
    Particles%h_min = -1._MK

    Particles%areinside = .TRUE.
    ! neighbour lists not updated
    Particles%neighlists = .FALSE.
    ! set cutoff to a default value
    IF (PRESENT(cutoff)) THEN
        Particles%cutoff = cutoff
    ELSE
        Particles%cutoff = 2.1_MK * Particles%h_avg
    ENDIF

    !-----------------------------------------------------------------------
    ! Finalize
    !-----------------------------------------------------------------------
    CALL substop(caller,t0,info)

    9999 CONTINUE ! jump here upon error

#if __DIM == 2
END SUBROUTINE particles_initialize2d
#elif __DIM == 3
END SUBROUTINE particles_initialize3d
#endif

#undef __DIM
