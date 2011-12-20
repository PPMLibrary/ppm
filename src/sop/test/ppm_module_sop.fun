test_suite ppm_module_sop
use ppm_module_particles
use ppm_module_sop_typedef
use ppm_module_dcops
use ppm_module_io_vtk
use ppm_module_data !, ONLY: ppm_mpi_kind

#ifdef __MPI
    INCLUDE "mpif.h"
#endif

integer, parameter              :: debug = 0
integer, parameter              :: mk = kind(1.0d0) !kind(1.0e0)
real(mk),parameter              :: tol=epsilon(1._mk)*100
real(mk),parameter              :: pi = 3.1415926535897931_mk
real(mk),parameter              :: skin = 0._mk
integer,parameter               :: ndim=2
integer                         :: decomp,assig,tolexp
integer                         :: info,comm,rank,nproc
integer                         :: topoid,nneigh_theo
integer                         :: np_global = 300
integer                         :: npart_g
real(mk),parameter              :: cutoff = 0.45_mk
real(mk),dimension(:,:),pointer :: xp=>NULL(),disp=>NULL()
real(mk),dimension(:,:),pointer :: grad_wp=>NULL()
real(mk),dimension(:),pointer   :: err_x=>NULL(),err_y=>NULL(),err_z=>NULL()
real(mk),dimension(:  ),pointer :: min_phys,max_phys
real(mk),dimension(:  ),pointer :: len_phys
real(mk),dimension(:  ),pointer :: rcp,wp
integer                         :: i,j,k,isum1,isum2,ip,ineigh,iq
integer                         :: wp_id,wp1_id,wpv_id,eta_id
integer                         :: grad_wp_id,grad_exact_id
integer                         :: err_x_id,err_y_id,err_z_id
real(mk)                        :: rsum1,rsum2,dt
integer                         :: nstep
real(mk),dimension(:),pointer   :: delta
integer,dimension(3)            :: ldc
integer, dimension(6)           :: bcdef
real(mk),dimension(:  ),pointer :: cost=>NULL()
character(len=ppm_char)         :: dirname,mesg
logical                         :: ok
real(mk)                        :: t0,t1,t2,t3
type(ppm_t_particles_d),pointer :: Particles=>NULL()
type(sop_t_opts_d),pointer      :: opts=>NULL()
type(sop_t_stats_d),pointer     :: sop_stats=>NULL()
integer                         :: seedsize
integer,  dimension(:),allocatable :: seed
integer, dimension(:),pointer   :: nvlist=>NULL()
integer, dimension(:,:),pointer :: vlist=>NULL()
integer,dimension(:),allocatable:: degree,order
real(mk),dimension(:),allocatable:: coeffs
integer                         :: nterms
real(mk)                        :: linf,min_sv
real(mk),dimension(:),allocatable :: exact_vec,err_vec


    init

        use ppm_module_typedef
        use ppm_module_init
        use ppm_module_mktopo
        
        allocate(min_phys(ndim),max_phys(ndim),len_phys(ndim),&
            &         delta(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic
        
#ifdef __MPI
        comm = mpi_comm_world
        call mpi_comm_rank(comm,rank,info)
        call mpi_comm_size(comm,nproc,info)
#else
        rank = 0
        nproc = 1
#endif
        tolexp = int(log10(epsilon(1._mk)))+10
        call ppm_init(ndim,mk,tolexp,0,debug,info,99)

        call random_seed(size=seedsize)
        allocate(seed(seedsize))
        do i=1,seedsize
            seed(i)=10+i*i*(rank+1)
        enddo
        call random_seed(put=seed)

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        !decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0

        call ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,bcdef,&
            cutoff,cost,info)
    end init


    finalize
        use ppm_module_finalize

        call ppm_finalize(info)

        deallocate(min_phys,max_phys,len_phys,delta)

    end finalize


    setup


    end setup
        

    teardown
        
        call ppm_alloc_particles(Particles,np_global,ppm_param_dealloc,info)

    end teardown

    test adapt_particles
        ! test particle adaptation

        use ppm_module_typedef
        use ppm_module_topo_check

        use ppm_module_util_commopt


        !start with slightly perturbed cartesian particles
        call particles_initialize(Particles,np_global,info,&
                ppm_param_part_init_random,topoid)
        Assert_Equal(info,0)
        call particles_mapping_global(Particles,topoid,info)
        Assert_Equal(info,0)
        allocate(disp(ndim,Particles%Npart),stat=info)
        call random_number(disp)
        disp = Particles%h_avg * 1.5_mk * disp
        call particles_move(Particles,disp,info)
        Assert_Equal(info,0)
        deallocate(disp)
        call particles_apply_bc(Particles,topoid,info)
        Assert_Equal(info,0)
        call particles_mapping_partial(Particles,topoid,info)
        Assert_Equal(info,0)
        !init cutoff radii
        call particles_allocate_wps(Particles,Particles%rcp_id,info,&
            with_ghosts=.false.,name='rcp')
        Assert_Equal(info,0)
        wp1_id=0

        !init one property (not used, it will just be carried around) 
        call particles_allocate_wps(Particles,wp1_id,info,&
            with_ghosts=.false.,name='wp_test1')
        Assert_Equal(info,0)

        ! initialisation
        rcp => get_wps(Particles,Particles%rcp_id)
        wp => get_wps(Particles,wp1_id)
        xp => get_xp(Particles)
        FORALL(ip=1:Particles%Npart) 
            rcp(ip) = MIN(1.9_mk*Particles%h_avg,cutoff) !uniform cutoffs
            wp(ip) = 7._mk * f0_fun(xp(1:ndim,ip))
        END FORALL
        rcp => set_wps(Particles,Particles%rcp_id)
        wp => set_wps(Particles,wp1_id)
        xp => set_xp(Particles,read_only=.true.)

        ! call this routine when the cutoffs have been modified manually
        ! (it will check whether the ghosts have to be re-computed)
        call particles_updated_cutoff(Particles,info)
        Assert_Equal(info,0)

        !compute ghosts and neighbour lists
        call particles_mapping_ghosts(Particles,topoid,info)
        Assert_True(info.eq.0)
        call particles_neighlists(Particles,topoid,info)
        Assert_True(info.eq.0)

        ! printout 
        call ppm_vtk_particle_cloud('before_adapt0',Particles,info)
        Assert_Equal(info,0)

        ! Now, we will adapt the particles positions such that their
        ! resolution matches the monitor function (D_fun). D_fun
        ! is a function of the fields gradient. When these are
        ! known analytically, they are passed as the optional
        ! argument wp_grad_fun (and the field itself is passed as
        ! the optional argument wp_grad).

        ! We first adapt the particles in the case where the field and
        ! its gradient are known analytically. Then when they are not.

        !setup options for sop

        ! init data structure for options
        call sop_init_opts(opts,info)
        call sop_init_stats(sop_stats,info)
        Assert_Equal(info,0)

        ! set parameters
        opts%D_needs_gradients = .true.
        opts%add_parts = .true.
        opts%add_parts = .true.
        opts%param_morse = 2.5_mk
        opts%rcp_over_D = 2.0_mk
        opts%attractive_radius0 = 0.5_mk !0.4_mk
        opts%adaptivity_criterion = 9.5_mk
        opts%fuse_radius = 0.05_mk
        opts%c = 1.4_mk
        opts%nneigh_critical = 15
        opts%spawn_radius = 0.5_mk

        if (ndim .eq. 2) then
            opts%scale_D = 0.05_mk
            opts%minimum_D = 0.00001_mk
            opts%maximum_D = 0.10_mk
        else
            opts%scale_D = 0.05_mk
            opts%minimum_D = 0.01_mk
            opts%maximum_D = 0.05_mk
        endif

        Particles%itime = 0

        !one can print-to-file the interaction potential that will be used
        !call sop_plot_potential(opts,'potential.dat',info)

        !--------------------------------------------------
        ! adapt particles using analytical initial condition
        !--------------------------------------------------
        call sop_adapt_particles(topoid,Particles,D_fun,opts,info,&
            wp_fun=f0_fun,wp_grad_fun=f0_grad_fun)
        Assert_Equal(info,0)

        !printout
        call ppm_vtk_particle_cloud('after_adapt0',Particles,info,&
            with_nvlist=.true.)

        ! Now we adapt the particles in the case where the field and
        ! its gradient are NOT known analytically.

        !Define a property, which will be used for adaptation
        ! SOP expects this property to be the field on which it has
        ! to adapt the particles positions. It will evaluate its gradients
        ! numerically using DC operators.
        call particles_allocate_wps(Particles,Particles%adapt_wpid,&
            info,name='wp_adapt')
        wp_id = Particles%adapt_wpid
        Assert_Equal(info,0)

        xp => get_xp(Particles)
        wp => get_wps(Particles,wp_id)
        FORALL(ip=1:Particles%Npart) wp(ip) = f0_fun(xp(1:ndim,ip)) 
        xp => set_xp(Particles,read_only=.true.)
        wp => set_wps(Particles,wp_id)

        call particles_mapping_ghosts(Particles,topoid,info)
        !call particles_neighlists(Particles,topoid,info,knn=opts%nneigh_critical+15)
        call particles_neighlists(Particles,topoid,info)

        !check that Derivatives can be computed accurately
        nterms=ndim
        allocate(degree(nterms*ndim),coeffs(nterms),order(nterms))
        if (ndim .eq. 2) then
               degree =  (/1,0,   0,1/)
        else 
               degree =  (/1,0,0, 0,1,0, 0,0,1/)
        endif
        coeffs = 1.0_mk; order =  4; eta_id = 0
        call particles_dcop_define(Particles,eta_id,coeffs,degree,order,nterms,&
                info,name="gradient",vector=.true.)
        deallocate(degree,coeffs,order)

        min_sv = HUGE(1._MK)
        call particles_dcop_compute(Particles,eta_id,info,&
                        c=opts%c,min_sv=min_sv)
        Assert_Equal(info,0)
        WRITE(*,*) 'Smallest singular value for TEST gradient was ',min_sv

        grad_wp_id=0
        call particles_allocate_wpv(Particles,grad_wp_id,ppm_dim,info,name="grad_wp")
        call particles_dcop_apply(Particles,wp_id,grad_wp_id,eta_id,info)
        Assert_Equal(info,0)
        grad_exact_id=0
        call particles_allocate_wpv(Particles,grad_exact_id,ppm_dim,info,&
            name="grad_exact")

        grad_wp => Get_wpv(Particles,grad_exact_id)
        xp => Get_xp(Particles)
        FORALL(ip=1:Particles%Npart) grad_wp(1:ndim,ip) = f0_grad_fun(xp(1:ndim,ip)) 
        grad_wp => Set_wpv(Particles,grad_exact_id)
        xp => Set_xp(Particles,read_only=.true.)


        wp => Get_wps(Particles,wp_id)
        grad_wp => Get_wpv(Particles,grad_wp_id)
        xp => Get_xp(Particles)
        linf = 0._mk
        allocate(exact_vec(nterms),err_vec(nterms))
        err_vec = 0._mk

        err_x_id=0
        call particles_allocate_wps(Particles,err_x_id,info,name="error_grad_x")
        err_y_id=0
        call particles_allocate_wps(Particles,err_y_id,info,name="error_grad_y")
        !err_z_id=0
        !call particles_allocate_wps(Particles,err_z_id,info,name="error_grad_z")
        err_x => Get_wps(Particles,err_x_id)
        err_y => Get_wps(Particles,err_y_id)
        !err_z => Get_wps(Particles,err_z_id)
        DO ip=1,Particles%Npart
            exact_vec = f0_grad_fun(xp(1:ndim,ip))
            DO i=1,nterms 
                err_vec(i) = MAX(err_vec(i),abs(grad_wp(i,ip) - exact_vec(i)))
            ENDDO
            err_x(ip) = abs(grad_wp(1,ip) - exact_vec(1))
            err_y(ip) = abs(grad_wp(2,ip) - exact_vec(2))
            !err_z(ip) = abs(grad_wp(3,ip) - exact_vec(3))
            linf = MAX(linf,MAXVAL(abs(exact_vec)))
        ENDDO

        !removme
        !visualize stencil that gave rise to the biggest error
        !wp(1:Particles%Npart) = 0._mk
        !ip = MAXLOC(err_x,1)
        !wp(ip) = 2._mk
        !DO ineigh=1,Particles%nvlist(ip)
        !    iq=Particles%vlist(ineigh,ip)
        !    wp(iq) = 1._mk
        !ENDDO
        !removme

        err_x => Set_wps(Particles,err_x_id)
        err_y => Set_wps(Particles,err_y_id)
        !err_z => Set_wps(Particles,err_z_id)
#ifdef __MPI
        call MPI_Allreduce(linf,linf,1,ppm_mpi_kind,MPI_MAX,comm,info)
#endif
        wp => Set_wps(Particles,wp_id,read_only=.TRUE.)
        grad_wp => Set_wpv(Particles,grad_wp_id,read_only=.TRUE.)
        xp => Set_xp(Particles,read_only=.TRUE.)
        call ppm_vtk_particle_cloud('grad',Particles,info)

        Assert_Equal(info,0)
        write(*,*) 'L-infinity error is ', MAXVAL(err_vec)/linf
        Assert_True(MAXVAL(err_vec)/linf.LT.0.1)

        !printout errors
        call ppm_vtk_particle_cloud('errors_in_gradients',Particles,info,&
            with_nvlist=.true.)

        deallocate(exact_vec,err_vec)
        call particles_dcop_free(Particles,eta_id,info)
        call particles_allocate_wps(Particles,err_x_id,info,&
            iopt=ppm_param_dealloc)
        call particles_allocate_wps(Particles,err_y_id,info,&
            iopt=ppm_param_dealloc)
        call particles_allocate_wpv(Particles,grad_wp_id,ndim,info,&
            iopt=ppm_param_dealloc)
        call particles_allocate_wpv(Particles,grad_exact_id,ndim,info,&
            iopt=ppm_param_dealloc)


        Particles%itime = 1
        opts%check_dcops=.true.
        !--------------------------------------------------
        ! adapt particles using discretized fields
        !--------------------------------------------------
        call sop_adapt_particles(topoid,Particles,D_fun,opts,info,stats=sop_stats)
        Assert_Equal(info,0)

        !printout
        call ppm_vtk_particle_cloud('after_adapt1',Particles,info,&
            with_nvlist=.true.)

        !do nothing and adapt again
        call particles_mapping_ghosts(Particles,topoid,info)
        call particles_neighlists(Particles,topoid,info)!,knn=opts%nneigh_critical+5)
        write(*,*) 'NOW, neighmin = ',Particles%nneighmin

        Particles%itime = 2
        !--------------------------------------------------
        ! adapt particles using discretized fields
        !--------------------------------------------------
        call sop_adapt_particles(topoid,Particles,D_fun,opts,info)
        Assert_Equal(info,0)

        !printout
        call ppm_vtk_particle_cloud('after_adapt2',Particles,info,&
            with_nvlist=.true.)

        !move particles and adapt one more time
        allocate(disp(ndim,Particles%Npart),stat=info)
        xp => get_xp(Particles)
        dt = Particles%h_min
        FORALL (ip=1:Particles%Npart)
            disp(1:2,ip) = dt*2._mk*&
                (/-sin(pi*xp(1,ip))**2*sin(pi*xp(2,ip))*cos(pi*xp(2,ip)),&
                   sin(pi*xp(2,ip))**2*sin(pi*xp(1,ip))*cos(pi*xp(1,ip)) /)
        END FORALL
        if (ndim .eq. 3) then
            FORALL (ip=1:Particles%Npart)
                disp(3,ip) = dt*2._mk*&
                    sin(pi*xp(2,ip))**2*sin(pi*xp(1,ip))*cos(pi*xp(1,ip))
            END FORALL
        endif
!write(*,*) 'dt = ',dt
!write(*,*) 'max_disp = ',MAXVAL(abs(disp(:,1:Particles%Npart)))
!write(*,*) 'min_xp = ',MINVAL(xp(:,1:Particles%Npart))
!write(*,*) 'max_xp = ',MAXVAL(xp(:,1:Particles%Npart))
        xp => set_xp(Particles,read_only=.true.)
        call particles_move(Particles,disp,info)
        Assert_Equal(info,0)
        deallocate(disp)
        call particles_apply_bc(Particles,topoid,info)
        Assert_Equal(info,0)
        call particles_mapping_partial(Particles,topoid,info)
        Assert_Equal(info,0)
        call particles_mapping_ghosts(Particles,topoid,info)
        Assert_Equal(info,0)
        call particles_neighlists(Particles,topoid,info)!,knn=opts%nneigh_critical+5)
        Assert_Equal(info,0)

        Particles%itime = 3

        call sop_adapt_particles(topoid,Particles,D_fun,opts,info)
        Assert_Equal(info,0)

        call ppm_vtk_particle_cloud('after_adapt3',Particles,info)

    end test

pure function f0_fun(pos)

    use ppm_module_data, ONLY: ppm_dim
    real(mk)                                 :: f0_fun
    real(mk), dimension(ppm_dim), intent(in) :: pos
    real(mk), dimension(ppm_dim)             :: centre
    real(mk)                                 :: radius,eps

    centre = 0.5_mk
    centre(2) = 0.5_mk
    radius=0.15_mk
    eps = 0.003_mk

    f0_fun = tanh((sqrt(sum((pos(1:ppm_dim)-centre)**2)) - radius)/eps)

    !f0_fun = 1._mk

end function f0_fun

pure function f0_grad_fun(pos)

    use ppm_module_data, ONLY: ppm_dim
    real(mk), dimension(ppm_dim)             :: f0_grad_fun
    real(mk), dimension(ppm_dim), intent(in) :: pos
    real(mk), dimension(ppm_dim)             :: centre
    real(mk)                                 :: radius,eps,f0,d

    centre = 0.5_mk
    centre(2) = 0.5_mk
    radius=0.15_mk
    eps = 0.003_mk

    d = sqrt(sum((pos(1:ppm_dim)-centre)**2))
    f0 = tanh((d - radius)/eps)
    f0_grad_fun = (1._mk - f0**2) * (pos(1:ppm_dim)-centre)/(eps*d)

    !f0_grad_fun = 0._mk

end function f0_grad_fun

pure function level0_fun(pos)

    use ppm_module_data, ONLY: ppm_dim
    real(mk)                              :: level0_fun
    real(mk), dimension(ppm_dim), intent(in) :: pos
    real(mk), dimension(ppm_dim)             :: centre
    real(mk)                              :: radius

    centre = 0.5_mk
    centre(2) = 0.5_mk
    radius=0.15_mk

    level0_fun = sqrt(sum((pos(1:ppm_dim)-centre)**2)) - radius

end function level0_fun

pure function D_fun(wp,wp_grad,opts,level)
    use ppm_module_data, ONLY: ppm_dim
    real(mk)                               :: D_fun
    real(mk),                   intent(in) :: wp
    real(mk),dimension(ppm_dim),intent(in) :: wp_grad
    type(sop_t_opts_d),pointer,   intent(in) :: opts
    real(mk), optional,         intent(in) :: level
    real(mk)                               :: lengthscale

    D_fun =  opts%scale_d / sqrt(1_mk + sqrt(SUM(wp_grad**2)))

end function D_fun


end test_suite
