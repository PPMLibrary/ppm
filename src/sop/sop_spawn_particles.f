SUBROUTINE sop_spawn_particles(Particles,opts,info,nb_part_added,&
        nneigh_threshold,wp_fun,printp)
    !!!----------------------------------------------------------------------------!
    !!!
    !!! Insert particles around those with too few neighbours
    !!!
    !!! Uses ppm_alloc to grow/shrink arrays 
    !!!
    !!! Warning: on output, some particles may be outside the computational domain
    !!! Call impose_boundary_conditions to fix this.
    !!!----------------------------------------------------------------------------!

    USE ppm_module_alloc, ONLY: ppm_alloc

    IMPLICIT NONE
#ifdef __MPI
    INCLUDE 'mpif.h'
#endif
#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    ! arguments
    TYPE(ppm_t_particles), POINTER,       INTENT(INOUT)   :: Particles
    TYPE(sop_t_opts), POINTER,            INTENT(IN   )   :: opts
    INTEGER,                              INTENT(  OUT)   :: info

    INTEGER, OPTIONAL,                    INTENT(  OUT)   :: nb_part_added
    INTEGER, OPTIONAL                                     :: nneigh_threshold
    OPTIONAL                                              :: wp_fun
    INTEGER, OPTIONAL                                     :: printp
    !!! printout particles that are created into file fort.(6000+printp)
    ! argument-functions need an interface
    INTERFACE
        FUNCTION wp_fun(pos)
            USE ppm_module_data, ONLY: ppm_dim
            USE ppm_module_typedef
#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
            REAL(MK),DIMENSION(ppm_dim),INTENT(IN)           :: pos
            REAL(MK)                                      :: wp_fun
        END FUNCTION wp_fun

    END INTERFACE

    ! local variables
    REAL(MK),     DIMENSION(:,:),POINTER   :: xp => NULL()
    REAL(MK),     DIMENSION(:,:),  POINTER :: inv => NULL()
    REAL(MK),     DIMENSION(:,:),  POINTER :: D => NULL()
    REAL(MK),     DIMENSION(:),  POINTER   :: wp => NULL()
    REAL(MK),     DIMENSION(:),  POINTER   :: level => NULL()
    INTEGER                                :: Npart, Mpart
    INTEGER,      DIMENSION(:),  POINTER   :: nvlist => NULL()
    INTEGER                                :: ip,iq,ineigh,i,di,j
    CHARACTER(LEN=256)                     :: cbuf
    CHARACTER(LEN=256)                     :: caller='sop_spawn_particles'
    REAL(KIND(1.D0))                       :: t0
    REAL(MK)                               :: lev
    REAL(MK)                               :: theta1,theta2
    INTEGER                                :: nvlist_theoretical
    INTEGER                                :: add_part, num_try
    INTEGER,        DIMENSION(2)           :: lda
    INTEGER,        DIMENSION(1)           :: lda1
    INTEGER, PARAMETER                     :: nb_new_part = 5
    REAL(MK)                               :: angle, dist, dist1, dist2, leng, min_dist, minmin_d
    REAL(MK), DIMENSION(ppm_dim)           :: displace
#ifdef __USE_RANDOMNUMBERS
    LOGICAL                                :: alloc_rand
#endif
    REAL(MK),     DIMENSION(:),   POINTER      :: Matrix_A => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: Matrix_B => NULL()
    LOGICAL                                :: too_close
    !!! number of new particles that are generated locally

    !!-------------------------------------------------------------------------!
    !! Initialize
    !!-------------------------------------------------------------------------!
    info = 0
#if debug_verbosity > 0
    CALL substart(caller,t0,info)
#endif

    IF (PRESENT(nneigh_threshold)) THEN
        nvlist_theoretical = nneigh_threshold
    ELSE
        nvlist_theoretical = opts%nneigh_theo
    ENDIF

    add_part = 0

    nvlist => Particles%nvlist
    Npart = Particles%Npart
    Mpart = Particles%Mpart
   
    !!-------------------------------------------------------------------------!
    !! Count number of new particles to insert
    !!-------------------------------------------------------------------------!
    ! counting how many particles have to be added
    ! if not enough neighbours, add nb_new_part particle

    ! HAECKIC: add opts% which sampling method is used: random, closestconsider and quadrant

    DO ip=1,Npart
       !write(*,*) nvlist(ip), nvlist_theoretical
       IF (nvlist(ip) .LT. nvlist_theoretical) &
             add_part = add_part + nvlist_theoretical - nvlist(ip)
    ENDDO

    write(*,*) 'new particles ', add_part

    nvlist => NULL()

#if debug_verbosity > 1
        IF (PRESENT(printp)) THEN
            OPEN(6000+printp)
            OPEN(7000+printp)
        ENDIF
#endif


    ! HAECKIC: arrays are extended compared to Mpart, because maybe needed for comparison

    !!-------------------------------------------------------------------------!
    !! Re-allocate (grow) arrays if necessary
    !!-------------------------------------------------------------------------!
    IF (add_part .GT. 0) THEN
        lda = (/ppm_dim,Mpart+add_part/)
        CALL ppm_alloc(Particles%xp,lda,ppm_param_alloc_grow_preserve,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,   &
                &    'allocation of xp failed',__LINE__,info)
            GOTO 9999
        ENDIF
        lda = (/Particles%tensor_length,Mpart+add_part/)
        CALL ppm_alloc(Particles%wpv(Particles%D_id)%vec,lda,&
            ppm_param_alloc_grow_preserve,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,   &
                &    'allocation of D failed',__LINE__,info)
            GOTO 9999
        ENDIF
        CALL ppm_alloc(Particles%wpv(Particles%G_id)%vec,lda,&
            ppm_param_alloc_grow_preserve,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,   &
                &    'allocation of rcp failed',__LINE__,info)
            GOTO 9999
        ENDIF

#ifdef __USE_RANDOMNUMBERS

    !!-------------------------------------------------------------------------!
    !! Draw random numbers to add noise on positions of new particles
    !!-------------------------------------------------------------------------!
        !generate new random numbers if we are running out of them
        alloc_rand = .FALSE.
        IF (.NOT. ALLOCATED(randnb)) THEN
            alloc_rand = .TRUE.
        ELSE IF((ppm_dim-1)*randnb_i+nb_new_part*ppm_dim*(add_part+1) .GE. SIZE(randnb)) THEN
            alloc_rand = .TRUE.
            DEALLOCATE(randnb)
        ENDIF
        IF(alloc_rand) THEN
            ALLOCATE(randnb(add_part*ppm_dim*(2*(Npart))),STAT=info)
            IF (info .NE. 0) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_alloc,caller,   &
                    &    'allocation of randnb failed',__LINE__,info)
                GOTO 9999
            ENDIF
            IF (.NOT.ASSOCIATED(ppm_particles_seed)) THEN
                CALL RANDOM_SEED(SIZE=ppm_particles_seedsize)
                ldc(1) = ppm_particles_seedsize
                CALL ppm_alloc(ppm_particles_seed,ldc(1:1),ppm_param_alloc_fit,info)
                IF (info .NE. 0) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_alloc,caller,   &
                        &            'allocation failed',__LINE__,info)
                    GOTO 9999
                ENDIF
                DO i=1,ppm_particles_seedsize
                    ppm_particles_seed(i)=ppm_rank*i*i*i*i
                ENDDO
                CALL RANDOM_SEED(PUT=ppm_particles_seed)
            ENDIF
            CALL RANDOM_NUMBER(randnb)
            randnb_i = 0
        ENDIF
#endif

    !!-------------------------------------------------------------------------!
    !! Add particles
    !! (either randomly, or using fake random numbers based on particles'
    !! positions. This is a quick hack to ensure exact consistency between
    !! simulations run an different number of processors)
    !!-------------------------------------------------------------------------!
    add_part = 0

    xp => Particles%xp !Cannot use get_xp, because we need access to elements
    ! of the array xp that are beyond Particles%Npart
    D => Particles%wpv(Particles%D_id)%vec      !same reason
    inv => Particles%wpv(Particles%G_id)%vec  !same reason
    nvlist => Particles%nvlist

    min_dist = 1.0_mk

   ! HAECKIC: TODO add sampling methods

    IF (ppm_dim .EQ. 2) THEN
        add_particles2d: DO ip=1,Npart

            IF (nvlist(ip) .LT. nvlist_theoretical) THEN
            !IF (Particles%nvlist(ip).GT.0) THEN

               ! HAECKIC: double sampling problem!
                DO i=1,(nvlist_theoretical-nvlist(ip))
                    add_part = add_part + 1

                    ! get the isotropic -> anisotropic transofmration matrix
                    CALL particles_inverse_matrix(D(:,ip),Matrix_A,info)

                    angle = 0._MK
#ifdef __USE_RANDOMNUMBERS
                    randnb_i = randnb_i + 1
                    angle = 2*PI*(randnb(randnb_i)+0.5_mk) 
#endif
                    ! Pseudo random
                    angle = angle + PI/2._mk * &
                        REAL(Particles%nvlist(ip),MK)


                    ! haeckic: do here the correct sampling
                    !          simply use the same inv and D?
                    ! sample in unit circle then transform it onto ellipse
                    

                    ! haeckic: add the minimum distance to other particles (=attractive radius)
                    ! we only compare with neighbors of ip
                    ! maybe: add a checking of neighbors of neighbors of ip
                    ! If after 100 tries no success, then do decrease min_dist
!                     too_close = .TRUE.
!                     num_try = 1
                    !write(*,*) 'start',xp(1:ppm_dim,ip), inv(2,ip), inv(3,ip)
!                     DO WHILE (too_close .AND. (num_try .LT. Npart))
! 
!                         angle = 0._MK
! 
! #ifdef __USE_RANDOMNUMBERS
!                         angle = 2*PI*randnb(num_try) 
! #endif
!                         ! Pseudo random
!                         !angle = angle + PI/2._mk * REAL(Particles%nvlist(ip),MK)
! 
!                         !where in ellipse: between attractive radius and 1
!                         ! attractive + rand*(1-attractive)
!                         leng = opts%attractive_radius0 + & 
!                                 & randnb(Npart+num_try)*(1.0_mk-opts%attractive_radius0)
!                         displace = leng*(/COS(angle),SIN(angle)/)
!                         displace = (/Matrix_A(1)*displace(1) + Matrix_A(2)*displace(2),&
!                            &           Matrix_A(3)*displace(1) + Matrix_A(4)*displace(2)/)
!                         xp(1:ppm_dim,Mpart + add_part) = xp(1:ppm_dim,ip) + displace
! 
!                      !write(*,*) xp(1:ppm_dim,Npart + add_part)!, angle, leng, (/COS(angle),SIN(angle)/)
!                         DO j=1,Particles%tensor_length
!                               D(j,Mpart + add_part)   = D(j,ip)
!                               inv(j,Mpart + add_part) = inv(j,ip)
!                         ENDDO
!                         
!                         ! Go through all -> optimize by going through neihbor lists
!                         too_close = .FALSE.
!                         minmin_d = 1.0_mk
!                         DO j=1,Particles%Mpart
!                            !IF (.NOT.(j .EQ. Npart + add_part)) THEN
!                               CALL particles_anisotropic_distance(Particles,Mpart + add_part,j,dist1,info)
!                               CALL particles_anisotropic_distance(Particles,j,Mpart + add_part,dist2,info)
!                               dist = MIN(dist1,dist2)
!                               IF (dist .LT. minmin_d) THEN
!                                     minmin_d = dist
!                               ENDIF
!                               IF(dist .LT. opts%attractive_radius0) THEN
!                                  too_close = .TRUE.
!                               ENDIF
!                            !ENDIF
!                         ENDDO
!       
!                         !against the newly added we need to consider periodicity
!                         DO j=Particles%Mpart+1,Particles%Mpart+add_part-1
!                            !IF (.NOT.(j .EQ. Npart + add_part)) THEN
! 
!                               !temporarily set particle's position, then undo
! 
!                               IF (ppm_dim .EQ. 2) THEN
!                                  
!                                  CALL particles_anisotropic_distance(Particles,Mpart + add_part,j,dist1,info)
!                                  CALL particles_anisotropic_distance(Particles,j,Mpart + add_part,dist2,info)
!                                  dist = MIN(dist1,dist2)
!                                  IF (dist .LT. minmin_d) THEN
!                                        minmin_d = dist
!                                  ENDIF
!                                  IF(dist .LT. opts%attractive_radius0) THEN
!                                     too_close = .TRUE.
!                                  ENDIF
! 
!                                  displace = (/1.0_mk ,0.0_mk /)
!                                  xp(1:ppm_dim,j) = xp(1:ppm_dim,j) + displace
! 
!                                  CALL particles_anisotropic_distance(Particles,Mpart + add_part,j,dist1,info)
!                                  CALL particles_anisotropic_distance(Particles,j,Mpart + add_part,dist2,info)
!                                  dist = MIN(dist1,dist2)
!                                  IF (dist .LT. minmin_d) THEN
!                                        minmin_d = dist
!                                  ENDIF
!                                  IF(dist .LT. opts%attractive_radius0) THEN
!                                     too_close = .TRUE.
!                                  ENDIF
! 
!                                  xp(1:ppm_dim,j) = xp(1:ppm_dim,j) - displace
! 
!                                  displace = (/0.0_mk ,1.0_mk /)
!                                  xp(1:ppm_dim,j) = xp(1:ppm_dim,j) + displace
! 
!                                  CALL particles_anisotropic_distance(Particles,Mpart + add_part,j,dist1,info)
!                                  CALL particles_anisotropic_distance(Particles,j,Mpart + add_part,dist2,info)
!                                  dist = MIN(dist1,dist2)
!                                  IF (dist .LT. minmin_d) THEN
!                                        minmin_d = dist
!                                  ENDIF
!                                  IF(dist .LT. opts%attractive_radius0) THEN
!                                     too_close = .TRUE.
!                                  ENDIF
! 
!                                  xp(1:ppm_dim,j) = xp(1:ppm_dim,j) - displace
! 
!                                  displace = (/-1.0_mk ,0.0_mk /)
!                                  xp(1:ppm_dim,j) = xp(1:ppm_dim,j) + displace
! 
!                                  CALL particles_anisotropic_distance(Particles,Mpart + add_part,j,dist1,info)
!                                  CALL particles_anisotropic_distance(Particles,j,Mpart + add_part,dist2,info)
!                                  dist = MIN(dist1,dist2)
!                                  IF (dist .LT. minmin_d) THEN
!                                        minmin_d = dist
!                                  ENDIF
!                                  IF(dist .LT. opts%attractive_radius0) THEN
!                                     too_close = .TRUE.
!                                  ENDIF
! 
!                                  xp(1:ppm_dim,j) = xp(1:ppm_dim,j) - displace
! 
!                                  displace = (/0.0_mk ,-1.0_mk /)
!                                  xp(1:ppm_dim,j) = xp(1:ppm_dim,j) + displace
! 
!                                  CALL particles_anisotropic_distance(Particles,Mpart + add_part,j,dist1,info)
!                                  CALL particles_anisotropic_distance(Particles,j,Mpart + add_part,dist2,info)
!                                  dist = MIN(dist1,dist2)
!                                  IF (dist .LT. minmin_d) THEN
!                                        minmin_d = dist
!                                  ENDIF
!                                  IF(dist .LT. opts%attractive_radius0) THEN
!                                     too_close = .TRUE.
!                                  ENDIF
! 
!                                  xp(1:ppm_dim,j) = xp(1:ppm_dim,j) - displace
! 
!                               ELSE
! 
!                               ENDIF
!    
!                               
!                            !ENDIF
!                         ENDDO
! 
!                         !Go through neighbors
! !                         too_close = .FALSE.
! !                         DO ineigh=1,nvlist(ip)
! !                            iq=Particles%vlist(ineigh,ip)
! ! 
! !                            CALL particles_anisotropic_distance(Particles,Npart + add_part,iq,dist1,info)
! !                            CALL particles_anisotropic_distance(Particles,iq,Npart + add_part,dist2,info)
! !                            dist = MIN(dist1,dist2)
! 
! !                            IF(dist .LT. opts%attractive_radius0) THEN
! !                               too_close = .TRUE.
! !                            ENDIF
! ! 
! !                         ENDDO
! 
!                         num_try = num_try + 1
!                     ENDDO
!                     
!                     !write(*,*) 'end'
!                     IF (minmin_d .LT. min_dist) THEN
!                         min_dist = minmin_d
!                     ENDIF
! 
!                     IF (num_try .GE. Npart) THEN
!                      write(*,*) 'Npart samples none found' 
!                     ENDIF

! (opts%attractive_radius0 + randnb(randnb_i)* &
!                     & (1.0_mk-2*opts%attractive_radius0))

                     ! the only random one
                    displace = 0.749_mk*(/COS(angle),SIN(angle)/)
                    displace = (/Matrix_A(1)*displace(1) + Matrix_A(2)*displace(2),&
                     &           Matrix_A(3)*displace(1) + Matrix_A(4)*displace(2)/)
                    xp(1:ppm_dim,Mpart + add_part) = xp(1:ppm_dim,ip) + displace

!                     write(*,*) ' From sample: ',  xp(1:ppm_dim,ip)
!                     write(*,*) ' angle: ', angle 
!                     write(*,*) ' axes: (', Matrix_A(1), ',',Matrix_A(3),') (',Matrix_A(2), ',',Matrix_A(4),')'
!                      write(*,*) 'A new sample: ', xp(1:ppm_dim,Npart + add_part)
!                      write(*,*) ' '

!                         0.723_MK*D(ip)&       !radius
!                         * (/COS(angle),SIN(angle)/) !direction 
                    DO j=1,Particles%tensor_length
                        ! why D? we update it anyway...
                        D(j,Mpart + add_part)   = D(j,ip)
                        inv(j,Mpart + add_part) = inv(j,ip)
                    ENDDO
                    
#if debug_verbosity > 1
                    IF (PRESENT(printp)) THEN
                        write(6000+printp,*) xp(1:ppm_dim,Npart+add_part)
                        write(7000+printp,'(4(E12.4,2X))') xp(1:ppm_dim,ip),&
                            xp(1:ppm_dim,ip)- xp(1:ppm_dim,Npart+add_part)
                    ENDIF
#endif
                ENDDO
            ENDIF
        ENDDO add_particles2d
    ELSE ! if ppm_dim .eq. 3
        add_particles3d: DO ip=1,Npart

            IF (nvlist(ip) .LT. nvlist_theoretical) THEN

                DO i=1,nb_new_part
                    add_part = add_part + 1

#ifdef __USE_RANDOMNUMBERS
                    randnb_i = randnb_i + 1
#endif

#ifdef __USE_RANDOMNUMBERS
                    theta1 = ACOS(1._MK - 2._MK*randnb(2*randnb_i))
                    theta2 = 2._MK * PI * randnb(2*randnb_i-1)
#else
                    theta1 = ACOS(SIN(1000._MK*xp(1,ip)/D(ip)))
                    theta2 = PI * (1._MK+COS(1000._MK*xp(2,ip)/D(ip)))
#endif

                    ! HAECKIC: TODO 3D sampling for anisotropic case

!                     xp(1:ppm_dim,Npart + add_part) = xp(1:ppm_dim,ip) + &
!                         !random 3D points on a sphere
!                     0.7_MK*D(ip)&       !radius
!                         * (/COS(theta1), &
!                         &   SIN(theta1) * COS(theta2), &
!                         &   SIN(theta1) * SIN(theta2)  &
!                         &   /) 
!                     D(Npart + add_part)   = D(ip)
!                     inv(Npart + add_part) = inv(ip)

#if debug_verbosity > 1
                    IF (PRESENT(printp)) THEN
                        write(6000+printp,*) xp(1:ppm_dim,Npart+add_part)
                    ENDIF
#endif
                ENDDO
            ENDIF
        ENDDO add_particles3d
    ENDIF

!write(*,*) 'closest: ', min_dist 

    !!-------------------------------------------------------------------------!
    !!new size Npart
    !!-------------------------------------------------------------------------!

    ! Set new Npart and append added particles
    DO ip=1,add_part
         xp(1:ppm_dim,Npart+ip) = xp(1:ppm_dim,Mpart+ip)
         DO j=1,Particles%tensor_length
               D(j,Npart+ip)   = D(j,Mpart+ip)
               inv(j,Npart+ip) = inv(j,Mpart+ip)
         ENDDO
    ENDDO

    xp => Set_xp(Particles)
    D => Set_wpv(Particles,Particles%D_id)
    inv => Set_wpv(Particles,Particles%G_id)
    nvlist => NULL()
    Particles%Npart = Npart + add_part

    CALL particles_updated_nb_part(Particles,info,&
        preserve_wpv=(/Particles%D_id,Particles%G_id/),&
        preserve_wps= (/ (i, i=1,0) /)) !F90-friendly way to init an empty array
        !preserve_wpv=(/ INTEGER :: /)) !valid only in F2003
    IF (info .NE.0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_sub_failed,caller,   &
            &    'particles_updated_nb_part failed',__LINE__,info)
        GOTO 9999
    ENDIF

ENDIF !add_part .NE.0

#if debug_verbosity > 1
IF (PRESENT(printp)) THEN
    CLOSE(6000+printp)
    CLOSE(7000+printp)
ENDIF
#endif

!!-------------------------------------------------------------------------!
!! Finalize
!!-------------------------------------------------------------------------!
#if debug_verbosity > 1
#ifdef __MPI
CALL MPI_Allreduce(add_part,add_part,1,MPI_INTEGER,MPI_SUM,ppm_comm,info)
#endif
IF (ppm_rank .EQ.0) THEN
    WRITE(cbuf,'(A,I8,A)') 'Adding ', add_part,' particles'
    CALL ppm_write(ppm_rank,caller,cbuf,info)
ENDIF
#endif
IF (PRESENT(nb_part_added)) THEN
    nb_part_added = add_part
ENDIF

#if debug_verbosity > 0
CALL substop(caller,t0,info)
#endif
9999 CONTINUE ! jump here upon error

END SUBROUTINE sop_spawn_particles

!haeckic: do the quadrant?
SUBROUTINE check_quadrants(Particles,info)
    IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    ! arguments
    TYPE(ppm_t_particles), POINTER,       INTENT(INOUT)   :: Particles
    INTEGER,                              INTENT(  OUT)   :: info
    !local variables

    REAL(MK),DIMENSION(:,:), POINTER                      :: xp=>NULL()
    INTEGER,DIMENSION(:),    POINTER                      :: nvlist=>NULL()
    INTEGER,DIMENSION(:,:),  POINTER                      :: vlist=>NULL()
    INTEGER                                               :: ip,iq,ineigh
    LOGICAL,DIMENSION(4)                                  :: needs_neigh_l

    info = 0
    nvlist => Particles%nvlist
    vlist => Particles%vlist

    xp => Get_xp(Particles,with_ghosts=.TRUE.)


    DO ip=1,Particles%Npart
        ineigh = 1
        needs_neigh_l = .TRUE.
        IF (nvlist(ip).GT.0) THEN
            DO WHILE(ANY(needs_neigh_l) .AND. ineigh.LE.nvlist(ip))
                iq = vlist(ineigh,ip)
                IF (xp(1,iq) .GT. xp(1,ip)) THEN
                    IF (xp(2,iq) .GT. xp(2,ip)) THEN
                        needs_neigh_l(1) = .FALSE.
                    ELSE
                        needs_neigh_l(4) = .FALSE.
                    ENDIF
                ELSE
                    IF (xp(2,iq) .GT. xp(2,ip)) THEN
                        needs_neigh_l(2) = .FALSE.
                    ELSE
                        needs_neigh_l(3) = .FALSE.
                    ENDIF
                ENDIF
                ineigh = ineigh + 1
            ENDDO
        ENDIF

        IF (ANY(needs_neigh_l)) THEN
            loop_quadrants: DO iq=1,4
                IF(needs_neigh_l(iq)) THEN
                    nvlist(ip) = iq
                    EXIT loop_quadrants
                ENDIF
            ENDDO loop_quadrants
        ELSE
            nvlist(ip) = 0
        ENDIF

    ENDDO
    xp => Set_xp(Particles,read_only=.TRUE.)
    vlist => NULL()
    nvlist => NULL()

END SUBROUTINE check_quadrants

! haeckic: do we need that?
SUBROUTINE check_duplicates(Particles)
    IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    ! arguments
    TYPE(ppm_t_particles), POINTER,       INTENT(INOUT)   :: Particles
    !local variables

    REAL(MK),DIMENSION(:,:), POINTER                      :: xp=>NULL()
    INTEGER,DIMENSION(:,:),  POINTER                      :: vlist=>NULL()
    INTEGER                                               :: ip,iq,ineigh
    INTEGER                                               :: info
    
    info = 0

    CALL particles_apply_bc(Particles,Particles%active_topoid,info)
    CALL particles_mapping_partial(Particles,Particles%active_topoid,info)
    CALL particles_mapping_ghosts(Particles,Particles%active_topoid,info)
    CALL particles_neighlists(Particles,Particles%active_topoid,info)

    vlist => Particles%vlist
    xp => Get_xp(Particles,with_ghosts=.TRUE.)

    DO ip=1,Particles%Npart
        DO ineigh = 1,Particles%nvlist(ip)
            iq = vlist(ineigh,ip)
            IF (SUM((xp(1:ppm_dim,iq) - xp(1:ppm_dim,ip))**2).LT.1E-10) THEN
                write(*,*) 'duplicate particles'
                write(*,*) 'ip = ',ip,' iq = ',iq
            stop
            ENDIF
        ENDDO
    ENDDO
    xp => Set_xp(Particles,read_only=.TRUE.)
    vlist => NULL()

END SUBROUTINE check_duplicates



#undef __KIND
