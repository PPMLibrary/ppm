test_suite ppm_module_inl_vlist

use ppm_module_typedef
use ppm_module_mktopo
use ppm_module_map
use ppm_module_topo_check
use ppm_module_neighlist
use ppm_module_particles

#ifdef __MPI
    INCLUDE "mpif.h"
#endif

integer, parameter              :: debug = 0
integer, parameter              :: mk = kind(1.0d0) !kind(1.0e0)
real(mk),parameter              :: pi = 3.1415926535897931_mk
real(mk),parameter              :: skin = 0._mk
integer,parameter               :: ndim=2
integer,parameter               :: pdim=2
integer                         :: decomp,assig,tolexp
real(mk)                        :: tol,min_rcp,max_rcp
integer                         :: info,comm,rank,nproc, num_neigh
integer                         :: topoid
integer                         :: np_global = 100
integer                         :: mp
integer                         :: newnp
real(mk),dimension(:,:),pointer :: xp => NULL()
real(mk),dimension(:,:),pointer :: disp => NULL()
real(mk),dimension(:  ),pointer :: rcp => NULL()
real(mk),dimension(:,:),pointer :: wp => NULL()
real(mk),dimension(:  ),pointer :: min_phys => NULL()
real(mk),dimension(:  ),pointer :: max_phys => NULL()
real(mk),dimension(:  ),pointer :: h => NULL()
real(mk),dimension(:  ),pointer :: p_h => NULL()
real(mk),dimension(:  ),pointer :: len_phys => NULL()
real(mk),dimension(:  ),pointer :: ghostlayer => NULL()
integer, dimension(:  ),pointer :: ghostsize => NULL()
integer                         :: i,j,k,sum1,sum2,n_ind
integer                         :: p_i, n_nei
integer, dimension(6)           :: bcdef
real(mk),dimension(:  ),pointer :: cost => NULL()
integer, dimension(:  ),pointer :: nm => NULL()
integer                         :: seedsize
integer,  dimension(:),allocatable :: seed
real(mk), dimension(:),allocatable :: randnb
integer                          :: isymm = 0
logical                          :: lsymm = .false.,ok
real(mk)                         :: t0,t1,t2,t3
real(mk)                         :: eps, dist,dist1,dist2
real(mk)                         :: cutoff = 0.15_mk
type(ppm_t_particles),pointer    :: Particles=>NULL()
type(ppm_t_particles),pointer    :: Particles2=>NULL()
REAL(MK),DIMENSION(:,:),POINTER  :: inv => NULL()
REAL(MK),DIMENSION(:),Pointer    :: Matrix_A => NULL()
REAL(MK),DIMENSION(:),Pointer    :: Matrix_B => NULL()

    init

        use ppm_module_typedef
        use ppm_module_init
        
        allocate(min_phys(ndim),max_phys(ndim),len_phys(ndim),&
            &         ghostsize(ndim),ghostlayer(2*ndim),&
            &         nm(ndim),h(ndim),p_h(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        ghostsize(1:ndim) = 2
        ghostlayer(1:2*ndim) = max_rcp
        bcdef(1:6) = ppm_param_bcdef_periodic
        
        eps = epsilon(1.0_mk)
        tolexp = int(log10(epsilon(1.0_mk)))
        
        nullify(xp,rcp,wp)

#ifdef __MPI
        comm = mpi_comm_world
        call mpi_comm_rank(comm,rank,info)
        call mpi_comm_size(comm,nproc,info)
#else
        rank = 0
        nproc = 1
#endif
        call ppm_init(ndim,mk,tolexp,0,debug,info,99)

    end init


    finalize
        use ppm_module_finalize

        call ppm_finalize(info)

        deallocate(min_phys,max_phys,len_phys,ghostsize,nm)

    end finalize


    setup

        call random_seed(size=seedsize)
        allocate(seed(seedsize))
        do i=1,seedsize
            seed(i)=10+i*i*(rank+1)
        enddo
        call random_seed(put=seed)

    end setup
        

    teardown
        
        deallocate(xp,rcp,wp,stat=info)
        deallocate(seed)

    end teardown

! test inl_neighlist_isotropic_2d
! 
!         use ppm_module_typedef
!         use ppm_module_mktopo
!         use ppm_module_map
!         use ppm_module_topo_check
!         use ppm_module_neighlist
!         use ppm_module_particles
! 
!         decomp = ppm_param_decomp_cuboid
!         assig  = ppm_param_assign_internal
! 
!         topoid = 0
! 
!         call ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,bcdef,&
!             cutoff,cost,info)
! 
!         call particles_initialize(Particles,np_global,info,&
!                 ppm_param_part_init_cartesian,topoid,cutoff=cutoff)
! 
!         call particles_mapping_global(Particles,topoid,info)
! 
!         call particles_allocate_wps(Particles,Particles%rcp_id,info,name='rcp')
!         call particles_update_cutoff(Particles,cutoff,info)
! 
!         allocate(disp(ndim,Particles%Npart),stat=info)
!         call random_number(disp)
!         disp = 0.03_mk * (disp-0.5_mk)
!         
!         call particles_move(Particles,disp,info)
!         deallocate(disp)
!         ! consider displacement
!         call particles_apply_bc(Particles,topoid,info)
! 
!         ! give particles to other proc which changed due to displacement
!         call particles_mapping_partial(Particles,topoid,info)
! 
!         ! get ghosts
!         call particles_mapping_ghosts(Particles,topoid,info)
! 
!         ! build neighlists
!         call particles_neighlists(Particles,topoid,info)
! 
!         ! check that we have all neighbors
!         ! it is assumed that the ghost mapping works
!         ! get cutoff radii of particles
! 
!         DO i=1,Particles%Npart
!                n_nei = 0
!                DO j=1,Particles%Mpart
! 
!                   dist = sqrt((Particles%xp(1,i)-Particles%xp(1,j))**2 + (Particles%xp(2,i)-Particles%xp(2,j))**2)
! 
!                   IF (dist .LE. Particles%wps(Particles%rcp_id)%vec(i) .AND. &
!                   &  dist .LE. Particles%wps(Particles%rcp_id)%vec(j) .AND. .NOT.(i.EQ.j)) THEN
! 
!                      !Check if we have neighbor in list
!                      ok = .FALSE.
!                      DO k=1,Particles%nvlist(i)
!                         n_ind = Particles%vlist(k,i)
!                         IF (abs(Particles%xp(1,j)-Particles%xp(1,n_ind)) .LT. 1e-12&
!                               & .AND. abs(Particles%xp(2,j)-Particles%xp(2,n_ind)) .LT. 1e-12) THEN
!                            ok = .TRUE.
!                         ENDIF
!                      ENDDO
!                      assert_true(ok)
!                      n_nei = n_nei + 1
! 
!                      ok = .FALSE.
!                      !check that neighbor has myself as neighbor
!                      IF (j.LE.Particles%Npart) THEN
!                         DO k=1,Particles%nvlist(j)
!                            n_ind = Particles%vlist(k,j)
!                            IF (abs(Particles%xp(1,i)-Particles%xp(1,n_ind)) .LT. 1e-12&
!                                  & .AND. abs(Particles%xp(2,i)-Particles%xp(2,n_ind)) .LT. 1e-12) THEN
!                               ok = .TRUE.
!                            ENDIF
!                         ENDDO
!                      ENDIF
!                      assert_true(ok)
! 
!                   ENDIF
! 
!                ENDDO
!               
!                !Check if n_nei is equal neighbor count
!                ok = n_nei .EQ. Particles%nvlist(i)
!                assert_true(ok)
! 
!         ENDDO
! 
!         !make random isotropic radii
!         allocate(disp(ndim,Particles%Npart),stat=info)
!         call random_number(disp)
!         DO i=1,Particles%Npart
!             Particles%wps(Particles%rcp_id)%vec(i) = 2*cutoff*disp(1,i)
!         ENDDO
!         deallocate(disp)
! 
!         ! updated ghosts...
!         call particles_updated_cutoff(Particles,info)
! 
!         ! make a new topology because cutoff max is bigger now
!         call ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,bcdef,&
!             Particles%cutoff,cost,info)
! 
!         ! get ghosts
!         call particles_mapping_ghosts(Particles,topoid,info)
! 
!         ! build neighlists
!         call particles_neighlists(Particles,topoid,info)
! 
!         ! check that we have all neighbors
!         ! it is assumed that the ghost mapping works
!         DO i=1,Particles%Npart
!                n_nei = 0
!                DO j=1,Particles%Mpart
! 
!                   dist = sqrt((Particles%xp(1,i)-Particles%xp(1,j))**2 + (Particles%xp(2,i)-Particles%xp(2,j))**2)
! 
!                   IF (dist .LE. Particles%wps(Particles%rcp_id)%vec(i) .AND. &
!                   &  dist .LE. Particles%wps(Particles%rcp_id)%vec(j) .AND. .NOT.(i.EQ.j)) THEN
! 
!                      !Check if we have neighbor in list
!                      ok = .FALSE.
!                      DO k=1,Particles%nvlist(i)
!                         n_ind = Particles%vlist(k,i)
!                         IF (abs(Particles%xp(1,j)-Particles%xp(1,n_ind)) .LT. 1e-12&
!                               & .AND. abs(Particles%xp(2,j)-Particles%xp(2,n_ind)) .LT. 1e-12) THEN
!                            ok = .TRUE.
!                         ENDIF
!                      ENDDO
!                      assert_true(ok)
!                      n_nei = n_nei + 1
!                   ENDIF
! 
!                ENDDO
!               
!                !Check if n_nei is equal neighbor count
!                ok = n_nei .EQ. Particles%nvlist(i)
!                assert_true(ok)
! 
!         ENDDO
! 
!         ! create another particles set and check cross sets anisotropic
!         call random_seed(size=seedsize)
!         do i=1,seedsize
!             seed(i)=10*2.21312+i*i*(rank+1)*0.5 + 21
!         enddo
!         call random_seed(put=seed)
!         
!         call particles_initialize(Particles2,np_global,info,&
!                 ppm_param_part_init_cartesian,topoid,cutoff=cutoff)
!         
!         call particles_mapping_global(Particles2,topoid,info)
!         
!         allocate(disp(ndim,Particles2%Npart),stat=info)
!         call random_number(disp)
!         disp = 0.03_mk * (disp-0.5_mk)
!         call particles_move(Particles2,disp,info)
!         deallocate(disp)
! 
!         ! consider displacement
!         call particles_apply_bc(Particles2,topoid,info)
! 
!         ! give particles to other proc which changed due to displacement
!         call particles_mapping_partial(Particles2,topoid,info)
! 
!         ! make the cross neighbor list
!         CALL particles_neighlists_xset(Particles2,Particles,topoid,info)
!       
!         DO i=1,Particles2%Npart
!                n_nei = 0
!                DO j=1,Particles%Mpart
! 
!                   ! If j_th particle should be a neighbor of i_th, check if it is in vlist_cross
!                   dist = sqrt((Particles2%xp(1,i)-Particles%xp(1,j))**2 + (Particles2%xp(2,i)-Particles%xp(2,j))**2)
! 
!                   IF (dist .LE. Particles%wps(Particles%rcp_id)%vec(j)) THEN
! 
!                      !Check if we have neighbor in list
!                      ok = .FALSE.
!                      DO k=1,Particles2%nvlist_cross(i)
!                         n_ind = Particles2%vlist_cross(k,i)
!                         IF (abs(Particles%xp(1,j)-Particles%xp(1,n_ind)) .LT. 1e-12&
!                               & .AND. abs(Particles%xp(2,j)-Particles%xp(2,n_ind)) .LT. 1e-12) THEN
!                            ok = .TRUE.
!                         ENDIF
!                      ENDDO
!                      assert_true(ok)
!                      n_nei = n_nei + 1
!                   ENDIF
! 
!                ENDDO
!               
!                !Check if n_nei is equal neighbor count
!                ok = n_nei .EQ. Particles2%nvlist_cross(i)
!                assert_true(ok)
!         ENDDO
! 
! 
!    end test

    test inl_neighlist_anisotropic_2d

        use ppm_module_typedef
        use ppm_module_mktopo
        use ppm_module_map
        use ppm_module_topo_check
	use ppm_module_neighlist
        use ppm_module_particles

        decomp = ppm_param_decomp_cuboid
        assig  = ppm_param_assign_internal

        topoid = 0

        call ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,bcdef,&
            cutoff,cost,info)

        call particles_initialize_anisotropic(Particles,np_global,cutoff,info,&
                ppm_param_part_init_cartesian,topoid)
        
        call particles_mapping_global(Particles,topoid,info)
        
        allocate(disp(ndim,Particles%Npart),stat=info)
        call random_number(disp)
        disp = 0.03_mk * (disp-0.5_mk)
        call particles_move(Particles,disp,info)
        deallocate(disp)

        ! consider displacement
        call particles_apply_bc(Particles,topoid,info)

        ! give particles to other proc which changed due to displacement
        call particles_mapping_partial(Particles,topoid,info)

        ! get ghosts
        call particles_mapping_ghosts(Particles,topoid,info)

        ! build neighlists
        call particles_neighlists(Particles,topoid,info)

        ! check that we have all neighbors
        ! it is assumed that the ghost mapping works
        DO i=1,Particles%Npart
               n_nei = 0
               DO j=1,Particles%Mpart

                  ! If j_th particle should be a neighbor, check if it is in vlist
                  CALL particles_anisotropic_distance(Particles,i,j,dist1,info)
                  CALL particles_anisotropic_distance(Particles,j,i,dist2,info)

                  dist = MAX(dist1,dist2)

                  IF (dist .LE. 1.0_mk .AND. .NOT.(i.EQ.j)) THEN

                     !Check if we have neighbor in list
                     ok = .FALSE.
                     DO k=1,Particles%nvlist(i)
                        n_ind = Particles%vlist(k,i)
                        IF (abs(Particles%xp(1,j)-Particles%xp(1,n_ind)) .LT. 1e-12&
                              & .AND. abs(Particles%xp(2,j)-Particles%xp(2,n_ind)) .LT. 1e-12) THEN
                           ok = .TRUE.
                        ENDIF
                     ENDDO
                     assert_true(ok)
                     n_nei = n_nei + 1
                     
                     ok = .FALSE.
                     !check that neighbor has myself as neighbor
                     IF (j.LE.Particles%Npart) THEN
                        DO k=1,Particles%nvlist(j)
                           n_ind = Particles%vlist(k,j)
                           IF (abs(Particles%xp(1,i)-Particles%xp(1,n_ind)) .LT. 1e-12&
                                 & .AND. abs(Particles%xp(2,i)-Particles%xp(2,n_ind)) .LT. 1e-12) THEN
                              ok = .TRUE.
                           ENDIF
                        ENDDO
                     ENDIF
                     !assert_true(ok)
                     
                  ENDIF

               ENDDO
              
               !Check if n_nei is equal neighbor count
               ok = n_nei .EQ. Particles%nvlist(i)
               assert_true(ok)

        ENDDO

        ! set the tensors randomly and remap and test everything again
        inv => get_wpv(Particles,Particles%G_id)
        CALL ppm_alloc(Matrix_A,(/ 4 /),ppm_param_alloc_fit,info)

        allocate(disp(ndim,Particles%Npart),stat=info)
        call random_number(disp)
        disp = 2*disp
        DO i=1,Particles%Npart
         
            ! let particles not get too small
            IF (disp(1,i) .LT. 0.5_mk) THEN
               disp(1,i) = 0.5_mk
            ENDIF

            IF (disp(2,i) .LT. 0.5_mk) THEN
               disp(2,i) = 0.5_mk
            ENDIF

            Matrix_A(1) = cutoff*disp(1,i)
            Matrix_A(3) = cutoff*disp(2,i)

            Matrix_A(2) = -cutoff*disp(2,i)/3
            Matrix_A(4) = cutoff*disp(1,i)/3

            CALL particles_inverse_matrix(Matrix_A, Matrix_B,info)

            inv(1,i) = Matrix_B(1)
            inv(2,i) = Matrix_B(2)
            inv(3,i) = Matrix_B(3)
            inv(4,i) = Matrix_B(4)

        ENDDO
        deallocate(disp)

        ! updated ghosts...
        call particles_updated_cutoff(Particles,info)

        ! make new topology
        call ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,bcdef,&
            Particles%cutoff,cost,info)

        ! get ghosts
        call particles_mapping_ghosts(Particles,topoid,info)

        ! build neighlists
        call particles_neighlists(Particles,topoid,info)

        ! check that we have all neighbors
        ! it is assumed that the ghost mapping works
        DO i=1,Particles%Npart
               n_nei = 0
               DO j=1,Particles%Mpart

                  ! If j_th particle should be a neighbor, check if it is in vlist
                  CALL particles_anisotropic_distance(Particles,i,j,dist1,info)
                  CALL particles_anisotropic_distance(Particles,j,i,dist2,info)

                  dist = MAX(dist1,dist2)

                  IF (dist .LE. 1.0_mk .AND. .NOT.(i.EQ.j)) THEN

                     !Check if we have neighbor in list
                     ok = .FALSE.
                     DO k=1,Particles%nvlist(i)
                        n_ind = Particles%vlist(k,i)
                        IF (abs(Particles%xp(1,j)-Particles%xp(1,n_ind)) .LT. 1e-12&
                              & .AND. abs(Particles%xp(2,j)-Particles%xp(2,n_ind)) .LT. 1e-12) THEN
                           ok = .TRUE.
                        ENDIF
                     ENDDO
                     assert_true(ok)
                     n_nei = n_nei + 1
                  ENDIF

               ENDDO
              
               !Check if n_nei is equal neighbor count
               ok = n_nei .EQ. Particles%nvlist(i)
               assert_true(ok)

        ENDDO

        ! create another particles set and check cross sets anisotropic
        
        call random_seed(size=seedsize)
        do i=1,seedsize
            seed(i)=10*2.21312+i*i*(rank+1)*0.5 + 21
        enddo
        call random_seed(put=seed)
        
        call particles_initialize_anisotropic(Particles2,np_global,cutoff,info,&
                ppm_param_part_init_cartesian,topoid)
        
        call particles_mapping_global(Particles2,topoid,info)
        
        allocate(disp(ndim,Particles2%Npart),stat=info)
        call random_number(disp)
        disp = 0.03_mk * (disp-0.5_mk)
        call particles_move(Particles2,disp,info)
        deallocate(disp)

        ! consider displacement
        call particles_apply_bc(Particles2,topoid,info)

        ! give particles to other proc which changed due to displacement
        call particles_mapping_partial(Particles2,topoid,info)

        ! make the cross neighbor list
        CALL particles_neighlists_xset(Particles2,Particles,topoid,info)
      
        DO i=1,Particles2%Npart
               n_nei = 0
               DO j=1,Particles%Mpart

                  ! If j_th particle should be a neighbor of i_th, check if it is in vlist_cross
                  CALL particles_sep_anisotropic_distance(Particles,Particles2,j,i,dist2,info)

                  IF (dist2 .LE. 1.0_mk) THEN

                     !Check if we have neighbor in list
                     ok = .FALSE.
                     DO k=1,Particles2%nvlist_cross(i)
                        n_ind = Particles2%vlist_cross(k,i)
                        IF (abs(Particles%xp(1,j)-Particles%xp(1,n_ind)) .LT. 1e-12&
                              & .AND. abs(Particles%xp(2,j)-Particles%xp(2,n_ind)) .LT. 1e-12) THEN
                           ok = .TRUE.
                        ENDIF
                     ENDDO
                     assert_true(ok)
                     n_nei = n_nei + 1
                  ENDIF

               ENDDO
              
               !Check if n_nei is equal neighbor count
               ok = n_nei .EQ. Particles2%nvlist_cross(i)
               assert_true(ok)
        ENDDO


    end test


end test_suite
