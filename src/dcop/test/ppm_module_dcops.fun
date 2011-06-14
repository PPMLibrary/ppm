test_suite ppm_module_dcops
use ppm_module_particles
use ppm_module_particles_typedef
use ppm_module_io_vtk
#include "../../ppm_define.h"

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
integer                         :: np_global = 1000
integer                         :: npart_g
real(mk),parameter              :: cutoff = 0.15_mk
real(mk),dimension(:,:),pointer :: xp=>NULL(),disp=>NULL()
real(mk),dimension(:  ),pointer :: min_phys=>NULL(),max_phys=>NULL()
real(mk),dimension(:  ),pointer :: len_phys=>NULL()
real(mk),dimension(:  ),pointer :: rcp=>NULL(),wp=>NULL(),dwp=>NULL()
real(mk),dimension(:,:),pointer :: grad_wp=>NULL()
integer                         :: i,j,k,isum1,isum2,ip,iq,ineigh,vect_wp_id
integer                         :: wp_id,dwp_id,grad_wp_id,eta_id,eta2_id
real(mk)                        :: coeff,err,exact,linf
real(mk),dimension(:),allocatable:: exact_vec,err_vec
integer                         :: nterms
real(mk),dimension(:),pointer   :: delta=>NULL()
integer,dimension(3)            :: ldc
integer,dimension(ndim)         :: dg
integer,dimension(:),allocatable:: degree,degree2,order,order2
real(mk),dimension(:),allocatable:: coeffs,coeffs2
integer, dimension(6)           :: bcdef
real(mk),dimension(:  ),pointer :: cost=>NULL()
character(len=ppm_char)         :: dirname
integer                         :: isymm = 0
logical                         :: lsymm = .false.,ok
real(mk)                        :: t0,t1,t2,t3
type(ppm_t_particles),pointer   :: Particles=>NULL()
type(ppm_t_particles),pointer   :: Particles2=>NULL()
integer                         :: seedsize
integer,  dimension(:),allocatable :: seed
integer, dimension(:),pointer   :: nvlist=>NULL()
integer, dimension(:,:),pointer :: vlist=>NULL()

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

        call ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,bcdef,cutoff,cost,info)
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

    test allocate_operator
        ! test data structure (mostly define and free)

        call particles_initialize(Particles,np_global,info,ppm_param_part_init_cartesian,topoid)
        call particles_mapping_global(Particles,topoid,info)

        allocate(degree(3*ndim),coeffs(3),order(3),degree2(7*ndim),coeffs2(7),order2(7))
        if (ndim .eq. 2) then
               degree =  (/1,0,    1,1,     0,1  /)
               degree2=  (/1,8,   0,1,   3,3,   1,1,   0,7,   0,2,   3,3/)
        else 
               degree =  (/1,0,0,  1,1,1,  0,0,1 /)
               degree2=  (/1,2,9, 1,1,8, 3,3,1, 1,1,1, 3,0,7, 2,0,2, 3,3,3/)
        endif
        coeffs = (/2.0_mk, 1.0_mk, -3.2_mk/)
        order =  (/2,       1,       1  /)
        coeffs2 = (/0.1_mk, -1.0_mk, -3.8_mk, -3.3_mk, 0.001_mk, 10._mk, 4._mk/)
        order2=  (/2,       2,       2,        2,      1,        2,      1/)

        call particles_dcop_deallocate(Particles,info) !this should do nothing
        Assert_Equal(info,0)

        eta_id = 0
        call particles_dcop_define(Particles,eta_id,coeffs,degree,order,3,info,name="test")
        Assert_Equal(info,0)
        Assert_False(eta_id.le.0)
        Assert_Equal(sum(abs(Particles%ops%desc(eta_id)%degree-degree)),0)
        Assert_Equal_Within(sum(abs(Particles%ops%desc(eta_id)%coeffs-coeffs)),0,1e-5)
        Assert_Equal(Particles%ops%nb_ops,1)

        eta2_id = 0
        call particles_dcop_define(Particles,eta2_id,coeffs2,degree2,order,7,info,name="test2")
        Assert_False(info.eq.0)
        call particles_dcop_define(Particles,eta2_id,coeffs2,degree2,order2,7,info,name="test2")
        Assert_Equal(info,0)
        Assert_False(eta2_id.le.0)
        Assert_Equal(sum(abs(Particles%ops%desc(eta2_id)%degree-degree2)),0)
        Assert_Equal_Within(sum(abs(Particles%ops%desc(eta2_id)%coeffs-coeffs2)),0,1e-5)
        Assert_Equal(Particles%ops%nb_ops,2)

!free the first operator and reallocate a new one
        call particles_dcop_free(Particles,eta_id,info)
        Assert_Equal(info,0)
        Assert_Equal(Particles%ops%nb_ops,1)
        Assert_Equal(Particles%ops%max_opsid,2)
        
        eta_id = 0
        call particles_dcop_define(Particles,eta_id,coeffs,degree,order,3,info,name="test1")
        Assert_Equal(info,0)
        Assert_False(eta_id.le.0)
        Assert_Equal(sum(abs(Particles%ops%desc(eta_id)%degree-degree)),0)
        Assert_Equal_Within(sum(abs(Particles%ops%desc(eta_id)%coeffs-coeffs)),0,1e-5)
        Assert_Equal(Particles%ops%nb_ops,2)
        Assert_Equal(Particles%ops%max_opsid,2)

        call particles_dcop_free(Particles,eta_id,info)
        Assert_Equal(info,0)
        call particles_dcop_free(Particles,eta2_id,info)
        Assert_Equal(info,0)
        
!deallocate all operators when there are no operators defined
        call particles_dcop_deallocate(Particles,info)
        Assert_Equal(info,0)
!defining an operator with an id that was not used before
        eta_id = 3
        call particles_dcop_define(Particles,eta_id,coeffs,degree,order,3,info,name="test1")
        Assert_Equal(info,0)
        Assert_Equal(eta_id,3)
!redefining an operator with an id that was already used before (overwritting)
        call particles_dcop_define(Particles,eta_id,coeffs,degree,order,3,info,name="test1")
        Assert_Equal(info,0)
        Assert_Equal(eta_id,3)

        deallocate(degree,coeffs,order,degree2,coeffs2,order2)
    end test

    test compute_operator
        ! test if we can compute the dc operators, then evaluate them on some test functions

        call particles_initialize(Particles,np_global,info,ppm_param_part_init_cartesian,topoid)
        allocate(disp(ndim,Particles%Npart))
        call random_number(disp)
        disp=0.15_mk*Particles%h_avg*disp
        call particles_move(Particles,disp,info)
        call particles_apply_bc(Particles,topoid,info)
        Assert_Equal(info,0)
        Particles%cutoff = Particles%h_avg * 3.3_mk
        call particles_mapping_global(Particles,topoid,info)
        wp_id=0
        call particles_allocate_wps(Particles,wp_id,info,name='wp')
        wp => Get_wps(Particles,wp_id)
        xp => Get_xp(Particles)
        FORALL(ip=1:Particles%Npart) wp(ip) = f0_fun(xp(1:ndim,ip),ndim)
        wp => Set_wps(Particles,wp_id)
        xp => Set_xp(Particles,read_only=.TRUE.)
        call particles_mapping_ghosts(Particles,topoid,info)
        call particles_neighlists(Particles,topoid,info)
        dwp_id=0
        call particles_allocate_wps(Particles,dwp_id,info,name='dwp')

        call ppm_vtk_particle_cloud('testvtk1',Particles,info)
        Assert_Equal(info,0)



!check d2dx2
        nterms=1
        allocate(degree(nterms*ndim),coeffs(nterms),order(nterms))
        if (ndim .eq. 2) then
               degree =  (/2,0/)
        else 
               degree =  (/2,0,0/)
        endif
        coeffs = 1.0_mk
        order =  2

        eta_id = 0
        call particles_dcop_define(Particles,eta_id,coeffs,degree,order,nterms,info,name="d2dx2")
        Assert_Equal(info,0)
        call particles_dcop_compute(Particles,eta_id,info)
        Assert_Equal(info,0)

        call particles_dcop_apply(Particles,wp_id,dwp_id,eta_id,info)
        Assert_Equal(info,0)

        wp => Get_wps(Particles,wp_id)
        dwp => Get_wps(Particles,dwp_id)
        xp => Get_xp(Particles)
        err = 0._mk
        linf = 0._mk
        DO ip=1,Particles%Npart
            exact = 0._mk
            DO i=1,nterms
                coeff = Particles%ops%desc(eta_id)%coeffs(i)
                dg = Particles%ops%desc(eta_id)%degree(1+(i-1)*ndim:i*ndim)
                exact = exact + coeff*df0_fun(xp(1:ndim,ip),dg,ndim)
            ENDDO
            write(101,'(5(E20.10,2X))') xp(1:ndim,ip),wp(ip),exact,dwp(ip)
            err = MAX(err,abs(dwp(ip) - exact))
            linf = MAX(linf,abs(exact))
        ENDDO
        wp => Set_wps(Particles,wp_id,read_only=.TRUE.)
        dwp => Set_wps(Particles,dwp_id,read_only=.TRUE.)
        xp => Set_xp(Particles,read_only=.TRUE.)
write(*,*) 'error is ', err/linf
        Assert_True(err/linf.LT.0.1)
        deallocate(degree,coeffs,order)
        call particles_dcop_free(Particles,eta_id,info)
        Assert_Equal(info,0)

!check Laplacian
        nterms=ndim
        allocate(degree(nterms*ndim),coeffs(nterms),order(nterms))
        if (ndim .eq. 2) then
               degree =  (/2,0,   0,2/)
        else 
               degree =  (/2,0,0, 0,2,0, 0,0,2/)
        endif
        coeffs = 1.0_mk
        order =  2

        eta_id = 0
        call particles_dcop_define(Particles,eta_id,coeffs,degree,order,nterms,info,name="laplacian")
        Assert_Equal(info,0)
        call particles_dcop_compute(Particles,eta_id,info)
        Assert_Equal(info,0)

        call particles_dcop_apply(Particles,wp_id,dwp_id,eta_id,info)
        Assert_Equal(info,0)

        wp => Get_wps(Particles,wp_id)
        dwp => Get_wps(Particles,dwp_id)
        xp => Get_xp(Particles)
        err = 0._mk
        linf = 0._mk
        DO ip=1,Particles%Npart
            exact = 0._mk
            DO i=1,nterms
                coeff = Particles%ops%desc(eta_id)%coeffs(i)
                dg = Particles%ops%desc(eta_id)%degree(1+(i-1)*ndim:i*ndim)
                exact = exact + coeff*df0_fun(xp(1:ndim,ip),dg,ndim)
            ENDDO
            write(102,'(5(E20.10,2X))') xp(1:ndim,ip),wp(ip),exact,dwp(ip)
            err = MAX(err,abs(dwp(ip) - exact))
            linf = MAX(linf,abs(exact))
        ENDDO
        wp => Set_wps(Particles,wp_id,read_only=.TRUE.)
        dwp => Set_wps(Particles,dwp_id,read_only=.TRUE.)
        xp => Set_xp(Particles,read_only=.TRUE.)
write(*,*) 'error is ', err/linf
        Assert_True(err/linf.LT.0.1)
        deallocate(degree,coeffs,order)
        call particles_dcop_free(Particles,eta_id,info)
        Assert_Equal(info,0)


!check something more complicated
        nterms= ndim
        allocate(degree(nterms*ndim),coeffs(nterms),order(nterms))
        if (ndim .eq. 2) then
               degree =  (/1,1,   4,1/)
        else 
               degree =  (/2,0,0, 1,2,0, 0,1,2/)
        endif
        coeffs = 1.0_mk
        coeffs(2) = 4.0_mk
        order =  2

        eta_id = 0
        call particles_dcop_define(Particles,eta_id,coeffs,degree,order,nterms,info,name="test")
        Assert_Equal(info,0)
        call particles_dcop_compute(Particles,eta_id,info)
        Assert_Equal(info,0)

        call particles_dcop_apply(Particles,wp_id,dwp_id,eta_id,info)
        Assert_Equal(info,0)

        wp => Get_wps(Particles,wp_id)
        dwp => Get_wps(Particles,dwp_id)
        xp => Get_xp(Particles)
        err = 0._mk
        linf = 0._mk
        DO ip=1,Particles%Npart
            exact = 0._mk
            DO i=1,nterms
                coeff = Particles%ops%desc(eta_id)%coeffs(i)
                dg = Particles%ops%desc(eta_id)%degree(1+(i-1)*ndim:i*ndim)
                exact = exact + coeff*df0_fun(xp(1:ndim,ip),dg,ndim)
            ENDDO
            write(103,'(5(E20.10,2X))') xp(1:ndim,ip),wp(ip),exact,dwp(ip)
            err = MAX(err,abs(dwp(ip) - exact))
            linf = MAX(linf,abs(exact))
        ENDDO
        wp => Set_wps(Particles,wp_id,read_only=.TRUE.)
        dwp => Set_wps(Particles,dwp_id,read_only=.TRUE.)
        xp => Set_xp(Particles,read_only=.TRUE.)
write(*,*) 'error is ', err/linf
        Assert_True(err/linf.LT.0.1)
        deallocate(degree,coeffs,order)
        call particles_dcop_free(Particles,eta_id,info)
        Assert_Equal(info,0)

!check vector-valued operators (gradient)
        nterms= ndim
        allocate(degree(nterms*ndim),coeffs(nterms),order(nterms))
        if (ndim .eq. 2) then
               degree =  (/1,0,   0,1/)
        else 
               degree =  (/1,0,0, 0,1,0, 0,0,1/)
        endif
        coeffs = 1.0_mk
        order =  2

        eta_id = 0
        call particles_dcop_define(Particles,eta_id,coeffs,degree,order,nterms,&
                info,name="gradient",vector=.true.)
        Assert_Equal(info,0)
        call particles_dcop_compute(Particles,eta_id,info)
        Assert_Equal(info,0)

        grad_wp_id=0
        call particles_allocate_wpv(Particles,grad_wp_id,ppm_dim,info,name="grad_wp")
        Assert_Equal(info,0)
        call particles_dcop_apply(Particles,wp_id,grad_wp_id,eta_id,info)
        Assert_Equal(info,0)

        wp => Get_wps(Particles,wp_id)
        grad_wp => Get_wpv(Particles,grad_wp_id)
        xp => Get_xp(Particles)
        err = 0._mk
        linf = 0._mk
        allocate(exact_vec(nterms),err_vec(nterms))
        DO ip=1,Particles%Npart
            DO i=1,nterms
                coeff = Particles%ops%desc(eta_id)%coeffs(i)
                dg = Particles%ops%desc(eta_id)%degree(1+(i-1)*ndim:i*ndim)
                exact_vec(i) = coeff*df0_fun(xp(1:ndim,ip),dg,ndim)
            ENDDO
            DO i=1,nterms 
                err_vec(i) = MAX(err_vec(i),abs(grad_wp(i,ip) - exact_vec(i)))
            ENDDO
            linf = MAX(linf,MAXVAL(abs(exact_vec)))
        ENDDO
        wp => Set_wps(Particles,wp_id,read_only=.TRUE.)
        grad_wp => Set_wpv(Particles,grad_wp_id,read_only=.TRUE.)
        xp => Set_xp(Particles,read_only=.TRUE.)
write(*,*) 'error is ', MAXVAL(err_vec)/linf
        Assert_True(MAXVAL(err_vec)/linf.LT.0.1)
        deallocate(degree,coeffs,order,exact_vec,err_vec)
        call particles_dcop_free(Particles,eta_id,info)
        Assert_Equal(info,0)

    end test

    test compute_operator_interp
        ! test if we can compute the dc operators with interpolating properties
        ! then evaluate them on some test functions

        call particles_initialize(Particles,np_global,info,ppm_param_part_init_cartesian,topoid)
        call particles_initialize(Particles2,np_global,info,ppm_param_part_init_cartesian,topoid)
        allocate(disp(ndim,Particles%Npart))
        call random_number(disp)
        disp=0.15_mk*Particles%h_avg*disp
        call particles_move(Particles,disp,info)
        call random_number(disp)
        disp=0.15_mk*Particles%h_avg*disp
        call particles_move(Particles2,disp,info)
        call particles_apply_bc(Particles,topoid,info)
        call particles_apply_bc(Particles2,topoid,info)
        Particles%cutoff = Particles%h_avg * 2.1_mk
        Particles2%cutoff = Particles2%h_avg * 2.1_mk
        call particles_mapping_global(Particles,topoid,info)
        call particles_mapping_global(Particles2,topoid,info)

        wp_id=0
        call particles_allocate_wps(Particles,wp_id,info,name='wp')
        wp => Get_wps(Particles,wp_id)
        xp => Get_xp(Particles)
        FORALL(ip=1:Particles%Npart) wp(ip) = f0_fun(xp(1:ndim,ip),ndim)
        wp => Set_wps(Particles,wp_id)
        xp => Set_xp(Particles,read_only=.TRUE.)
        call particles_mapping_ghosts(Particles,topoid,info)
        call particles_neighlists(Particles,topoid,info)

!compute nearest-neigbour distances
        call particles_allocate_wps(Particles,Particles%nn_sq_id,info,name='nn_sq')
        wp => Get_wps(Particles,Particles%nn_sq_id)
        xp => Get_xp(Particles,with_ghosts=.true.)
         forall(ip=1:Particles%Npart)  wp(ip) = huge(1._mk)
        do ip=1,Particles%Npart
            do ineigh=1,Particles%nvlist(ip)
                iq=Particles%vlist(ineigh,ip)
                wp(ip) = min(wp(ip),sum((xp(1:ndim,ip)-xp(1:ndim,iq))**2))
            enddo
        enddo
        wp => Set_wps(Particles,Particles%nn_sq_id)
        xp => Set_xp(Particles,read_only=.TRUE.)
        call particles_mapping_ghosts(Particles,topoid,info)
        call particles_mapping_ghosts(Particles2,topoid,info)
        call particles_neighlists(Particles2,topoid,info)
        call particles_neighlists_xset(Particles2,Particles,topoid,info)

        dwp_id=0
        call particles_allocate_wps(Particles2,dwp_id,info,name='dwp')

!check data interpolation
        nterms=1
        allocate(degree(nterms*ndim),coeffs(nterms),order(nterms))
        if (ndim .eq. 2) then
               degree =  (/0,0/)
        else 
               degree =  (/0,0,0/)
        endif
        coeffs = 1.0_mk
        order =  2

        eta_id = 0
        call particles_dcop_define(Particles2,eta_id,coeffs,degree,order,nterms,&
                info,name="interpolation",interp=.true.)
        Assert_Equal(info,0)
        call particles_dcop_compute(Particles2,eta_id,info)
        Assert_Equal(info,0)

        call particles_dcop_apply(Particles2,wp_id,dwp_id,eta_id,info)
        Assert_Equal(info,0)

        wp => Get_wps(Particles,wp_id)
        dwp => Get_wps(Particles2,dwp_id)
        xp => Get_xp(Particles2)
        err = 0._mk
        linf = 0._mk
        DO ip=1,Particles2%Npart
            exact = 0._mk
            DO i=1,nterms
                coeff = Particles2%ops%desc(eta_id)%coeffs(i)
                dg = Particles2%ops%desc(eta_id)%degree(1+(i-1)*ndim:i*ndim)
                exact = exact + coeff*df0_fun(xp(1:ndim,ip),dg,ndim)
            ENDDO
            err = MAX(err,abs(dwp(ip) - exact))
            linf = MAX(linf,abs(exact))
        ENDDO
        wp => Set_wps(Particles,wp_id,read_only=.TRUE.)
        dwp => Set_wps(Particles2,dwp_id,read_only=.TRUE.)
        xp => Set_xp(Particles2,read_only=.TRUE.)
write(*,*) 'error is ', err/linf
        Assert_True(err/linf.LT.0.1)
        deallocate(degree,coeffs,order)
        call particles_dcop_free(Particles2,eta_id,info)
        Assert_Equal(info,0)
        

!check data interpolation with derivatives
        call particles_updated_cutoff(Particles,info,Particles%h_avg*3.1_mk)
        call particles_updated_cutoff(Particles2,info,Particles2%h_avg*3.1_mk)
        call particles_mapping_ghosts(Particles,topoid,info)
        call particles_mapping_ghosts(Particles2,topoid,info)
        call particles_neighlists_xset(Particles2,Particles,topoid,info)
        nterms=2
        allocate(degree(nterms*ndim),coeffs(nterms),order(nterms))
        if (ndim .eq. 2) then
               degree =  (/0,0, 1,0/)
        else 
               degree =  (/0,0,0, 1,0,0/)
        endif
        coeffs = 1.0_mk
        order =  2
        eta_id = 0
        call particles_dcop_define(Particles2,eta_id,coeffs,degree,order,nterms,info,name="interp",interp=.true.)
        Assert_Equal(info,0)
        call particles_dcop_compute(Particles2,eta_id,info)
        Assert_Equal(info,0)

        call particles_dcop_apply(Particles2,wp_id,dwp_id,eta_id,info)
        Assert_Equal(info,0)

        wp => Get_wps(Particles,wp_id)
        dwp => Get_wps(Particles2,dwp_id)
        xp => Get_xp(Particles2)
        err = 0._mk
        linf = 0._mk
        DO ip=1,Particles2%Npart
            exact = 0._mk
            DO i=1,nterms
                coeff = Particles2%ops%desc(eta_id)%coeffs(i)
                dg = Particles2%ops%desc(eta_id)%degree(1+(i-1)*ndim:i*ndim)
                exact = exact + coeff*df0_fun(xp(1:ndim,ip),dg,ndim)
            ENDDO
            err = MAX(err,abs(dwp(ip) - exact))
            linf = MAX(linf,abs(exact))
        ENDDO
        wp => Set_wps(Particles,wp_id,read_only=.TRUE.)
        dwp => Set_wps(Particles2,dwp_id,read_only=.TRUE.)
        xp => Set_xp(Particles2,read_only=.TRUE.)
write(*,*) 'error is ', err/linf
        Assert_True(err/linf.LT.0.1)
        deallocate(degree,coeffs,order)
        call particles_dcop_free(Particles2,eta_id,info)
        Assert_Equal(info,0)

!check data interpolation with derivatives of several degrees
        nterms=5
        allocate(degree(nterms*ndim),coeffs(nterms),order(nterms))
        if (ndim .eq. 2) then
               degree =  (/0,0, 1,0, 2,0, 1,2, 0,4/)
        else 
               degree =  (/0,0,0, 1,0,0, 2,0,0, 2,0,2, 1,1,2/)
        endif
        coeffs = (/1._mk, -3._mk, 1.3_mk, 2.1_mk, 0.1_mk/)
        order =  (/2, 1, 3, 1, 2/)
        eta_id = 0
        call particles_dcop_define(Particles2,eta_id,coeffs,degree,order,nterms,info,name="everything",interp=.true.)
        Assert_Equal(info,0)
        call particles_dcop_compute(Particles2,eta_id,info)
        Assert_Equal(info,0)

        call particles_dcop_apply(Particles2,wp_id,dwp_id,eta_id,info)
        Assert_Equal(info,0)

        wp => Get_wps(Particles,wp_id)
        dwp => Get_wps(Particles2,dwp_id)
        xp => Get_xp(Particles2)
        err = 0._mk
        linf = 0._mk
        DO ip=1,Particles2%Npart
            exact = 0._mk
            DO i=1,nterms
                coeff = Particles2%ops%desc(eta_id)%coeffs(i)
                dg = Particles2%ops%desc(eta_id)%degree(1+(i-1)*ndim:i*ndim)
                exact = exact + coeff*df0_fun(xp(1:ndim,ip),dg,ndim)
            ENDDO
            err = MAX(err,abs(dwp(ip) - exact))
            linf = MAX(linf,abs(exact))
        ENDDO
        wp => Set_wps(Particles,wp_id,read_only=.TRUE.)
        dwp => Set_wps(Particles2,dwp_id,read_only=.TRUE.)
        xp => Set_xp(Particles2,read_only=.TRUE.)
write(*,*) 'error is ', err/linf
        Assert_True(err/linf.LT.0.1)
        deallocate(degree,coeffs,order)
        call particles_dcop_free(Particles2,eta_id,info)
        Assert_Equal(info,0)

!check data interpolation with vector-valued operators
        nterms=4
        allocate(degree(nterms*ndim),coeffs(nterms),order(nterms))
        if (ndim .eq. 2) then
               degree =  (/1,0, 2,0, 1,2, 0,4/)
        else 
               degree =  (/1,0,0, 2,0,0, 2,0,2, 1,1,2/)
        endif
        coeffs = (/-3._mk, 1.3_mk, 2.1_mk, 0.1_mk/)
        order =  (/1, 3, 1, 2/)
        eta_id = 0
        call particles_dcop_define(Particles2,eta_id,coeffs,degree,order,&
                nterms,info,name="vectvalued",interp=.true.,vector=.true.)
        Assert_Equal(info,0)
        call particles_dcop_compute(Particles2,eta_id,info)
        Assert_Equal(info,0)

        vect_wp_id = 0 !let dcop_apply allocate the array for the result
        call particles_dcop_apply(Particles2,wp_id,vect_wp_id,eta_id,info)
        Assert_Equal(info,0)


        call ppm_vtk_particle_cloud('testvtk',Particles,info)
        Assert_Equal(info,0)

        wp => Get_wps(Particles,wp_id)
        grad_wp => Get_wpv(Particles2,vect_wp_id)
        xp => Get_xp(Particles2)
        err = 0._mk
        linf = 0._mk
        allocate(exact_vec(nterms),err_vec(nterms))
        DO ip=1,Particles2%Npart
            DO i=1,nterms
                coeff = Particles2%ops%desc(eta_id)%coeffs(i)
                dg = Particles2%ops%desc(eta_id)%degree(1+(i-1)*ndim:i*ndim)
                exact_vec(i) = coeff*df0_fun(xp(1:ndim,ip),dg,ndim)
            ENDDO
            DO i=1,nterms 
                err_vec(i) = MAX(err_vec(i),abs(grad_wp(i,ip) - exact_vec(i)))
            ENDDO
            linf = MAX(linf,MAXVAL(abs(exact_vec)))
        ENDDO
        wp => Set_wps(Particles,wp_id,read_only=.TRUE.)
        grad_wp => Set_wpv(Particles2,vect_wp_id,read_only=.TRUE.)
        xp => Set_xp(Particles2,read_only=.TRUE.)
write(*,*) 'error is ', MAXVAL(err_vec)/linf
        Assert_True(MAXVAL(err_vec)/linf.LT.0.1)
        deallocate(degree,coeffs,order)
        call particles_dcop_free(Particles2,eta_id,info)
        Assert_Equal(info,0)

    end test

!-------------------------------------------------------------
! test function
!-------------------------------------------------------------
pure function f0_fun(pos,ndim)

    real(mk)                              :: f0_fun
    integer                 ,  intent(in) :: ndim
    real(mk), dimension(ndim), intent(in) :: pos

    f0_fun =  sin(2._mk*pi*pos(1)) * cos(2._mk*pi*pos(2)) 
end function f0_fun

!-------------------------------------------------------------
! derivatives of the test function
!-------------------------------------------------------------
function df0_fun(pos,order_deriv,ndim)

    real(mk)                                  :: df0_fun
    integer                     , intent(in)  :: ndim
    real(mk), dimension(ppm_dim), intent(in)  :: pos
    integer,  dimension(ppm_dim), intent(in)  :: order_deriv

    select case (order_deriv(1))
    case (0)
        df0_fun =    sin (2._mk*pi * pos(1)) 
    case (1)
        df0_fun =    2._mk*pi*cos(2._mk*pi * pos(1))
    case (2)
        df0_fun =  -(2._mk*pi)**2*sin(2._mk*pi * pos(1))
    case (3)
        df0_fun =  -(2._mk*pi)**3*cos(2._mk*pi * pos(1))
    case (4)
        df0_fun =   (2._mk*pi)**4*sin(2._mk*pi * pos(1))
    case default
        df0_fun =  0._mk
    endselect

    select case (order_deriv(2))
    case (0)
        df0_fun =   df0_fun * cos (2._mk*pi * pos(2)) 
    case (1)
        df0_fun =   df0_fun * (-2._mk*pi)*sin(2._mk*pi * pos(2))
    case (2)
        df0_fun =  df0_fun * (-(2._mk*pi)**2)*cos(2._mk*pi * pos(2))
    case (3)
        df0_fun =  df0_fun * ( (2._mk*pi)**3)*sin(2._mk*pi * pos(2))
    case (4)
        df0_fun =  df0_fun * ( (2._mk*pi)**4)*cos(2._mk*pi * pos(2))
    case default
        df0_fun =  0._mk
    endselect

    if (ndim .eq. 3 ) then
        select case (order_deriv(3))
        case (0)
            df0_fun =   df0_fun * sin (2._mk*pi * pos(3)) 
        case (1)
            df0_fun =   df0_fun * (2._mk*pi)*cos(2._mk*pi * pos(3))
        case (2)
            df0_fun =  df0_fun * (-(2._mk*pi)**2)*sin(2._mk*pi * pos(3))
        case (3)
            df0_fun =  df0_fun * (-(2._mk*pi)**3)*cos(2._mk*pi * pos(3))
        case (4)
            df0_fun =  df0_fun * ( (2._mk*pi)**4)*sin(2._mk*pi * pos(3))
        case default
            df0_fun =  0._mk
        endselect
    endif

end function df0_fun

end test_suite
