test_suite ppm_module_field_typedef

use ppm_module_mesh_typedef
use ppm_module_topo_typedef
use ppm_module_data
use ppm_module_mktopo
use ppm_module_finalize
use ppm_module_interfaces

#ifdef __MPI
    INCLUDE "mpif.h"
#endif

    integer, parameter              :: debug = 0
    integer, parameter              :: mk = kind(1.0d0) !kind(1.0e0)
#ifdef __MPI
    integer, parameter              :: comm = mpi_comm_world
#endif
    integer                         :: ndim
    integer                         :: nspec
    integer                         :: rank
    integer                         :: nproc
    integer                         :: decomp
    integer                         :: assig
    integer                         :: tolexp
    real(mk)                        :: tol
    integer                         :: info
    integer                         :: topoid=-1
    real(mk),dimension(:,:),pointer :: xp => NULL()
    real(mk),dimension(:,:),pointer :: wp => NULL()
    real(mk),dimension(:  ),pointer :: min_phys => NULL()
    real(mk),dimension(:  ),pointer :: max_phys => NULL()
    real(mk),dimension(:  ),pointer :: h => NULL()
    integer, dimension(:),  pointer :: ighostsize => NULL()
    real(mk), dimension(:), pointer :: ghostsize => NULL()
    integer                         :: i,j,k,p_i,ai,aj,it,isub
    integer, dimension(:  ),pointer :: bcdef => NULL()
    real(mk),dimension(:  ),pointer :: cost => NULL()
    integer, dimension(:,:),pointer :: istart => NULL()
    integer, dimension(:,:),pointer :: ndata => NULL()
    integer, dimension(:  ),pointer :: nm => NULL()
    integer                         :: np,mp
    integer                         :: kernel
    integer                         :: seedsize
    integer,  dimension(:),allocatable :: seed
        real(mk),dimension(:,:),pointer  :: field2d_1,field2d_2
        real(mk),dimension(:,:,:),pointer:: field3d_1,field3d_2
        real(mk),dimension(:,:,:),pointer:: field4d_1,field4d_2
        integer, dimension(2)            :: maxndata
        integer, dimension(:  ), pointer :: isublist => NULL()
        integer                          :: nsublist
        integer                          :: ipatch,mypatchid
        type(ppm_t_topo),      POINTER   :: topo => NULL()
        real(mk)                         :: sca_ghostsize


!-------------------------- init testsuit -------------------------------------
    init
        use ppm_module_data
        use ppm_module_init

        tol = 100.0_mk*epsilon(1.0_mk)
        tolexp = int(log10(epsilon(1.0_mk)))
        ndim = 2
        nspec = 1
        allocate(bcdef(2*ndim))
        allocate(min_phys(ndim),max_phys(ndim),ghostsize(ndim),&
        &        ighostsize(ndim),nm(ndim),h(ndim),stat=info)


#ifdef __MPI
        call mpi_comm_rank(comm,rank,info)
        call mpi_comm_size(comm,nproc,info)
#else
        nproc = 1
        rank = 0
#endif

        call ppm_init(ndim,mk,tolexp,0,debug,info,99)

    end init
!------------------------------------------------------------------------------


!------------------------- finalize testsuit ----------------------------------
    finalize
        use ppm_module_finalize

        call ppm_finalize(info)
        
        deallocate(min_phys,max_phys,ghostsize,nm)

    end finalize
!------------------------------------------------------------------------------


!------------------------------ test setup ------------------------------------
    setup
        
        nullify(xp)
        nullify(wp)
        nullify(cost)

        np = 400*nproc
        mp = 0
        
        call random_seed(size=seedsize)
        allocate(seed(seedsize))
        do i=1,seedsize
            seed(i)=10+i*i*(rank+1)
        enddo
        call random_seed(put=seed)

        bcdef(1:2*ndim) = ppm_param_bcdef_freespace
        kernel = ppm_param_rmsh_kernel_mp4
        do i=1,ndim
            min_phys(i) = 0.0_mk
            max_phys(i) = 1.0_mk
            ighostsize(i) = 2
            ghostsize(i) = 0.05_mk
        enddo
        
    end setup
!------------------------------------------------------------------------------
        

!--------------------------- test teardown ------------------------------------
    teardown
        
        if (associated(xp)) deallocate(xp)
        if (associated(wp)) deallocate(wp)
        if (associated(cost)) deallocate(cost)
        if (allocated(seed)) deallocate(seed)

    end teardown
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    test field_uniform_basics
        real(mk),dimension(ndim)         :: offset

        type(ppm_t_field)                :: Vort,Veloc
        type(ppm_t_equi_mesh)            :: Mesh1,Mesh2
        class(ppm_t_subpatch_),POINTER   :: p => NULL()

        integer,dimension(:),pointer :: ilist => NULL()
        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        assig  = ppm_param_assign_internal
        topoid = 0

        allocate(nm(ndim),stat=info)
        nm(1:ndim) = 16*nproc

        sca_ghostsize = 0.05_mk 

        call ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,    &
        &               bcdef,sca_ghostsize,cost,info)
        Assert_Equal(info,0)
        
        !--------------------------
        !Define Fields
        !--------------------------
        call Vort%create(1,info,name="Vorticity") !scalar field
        Assert_Equal(info,0)
        call Veloc%create(ndim,info,name="Velocity") !vector field
        Assert_Equal(info,0)

        !--------------------------
        !Create Mesh
        !--------------------------
        offset = 0._mk
        call Mesh1%create(topoid,offset,info,Nm=Nm)
        Assert_Equal(info,0)

        call Mesh1%def_uniform(info) 
        Assert_Equal(info,0)
        Assert_True(associated(Mesh1%subpatch))

        p => Mesh1%subpatch%begin()
        do while(associated(p))
           Assert_True(associated(p%subpatch_data))
           p => Mesh1%subpatch%next()
        enddo

        !--------------------------
        !Create data arrays on the mesh for the vorticity and velocity fields
        !--------------------------
        call Vort%discretize_on(Mesh1,info)
        Assert_Equal(info,0)

        call Veloc%discretize_on(Mesh1,info)
        Assert_Equal(info,0)

        !--------------------------
        ! Iterate through patches and initialize the data arrays
        !--------------------------
        p => Mesh1%subpatch%begin()

        do while (ASSOCIATED(p))
            call p%get_field(Vort,field2d_1,info)
            call p%get_field(Veloc,field3d_1,info)

            do i = 1,p%nnodes(1)
                do j = 1,p%nnodes(2)
                    field2d_1(i,j) = cos(i*h(1)+j)
                    field3d_1(1,i,j) = sin(field2d_1(i,j))
                    field3d_1(2,i,j) = cos(field2d_1(i,j))
                enddo
            enddo
            p => Mesh1%subpatch%next()
        enddo

        !Second version
        do ipatch = 1,Mesh1%subpatch%nb
            p => Mesh1%subpatch%vec(ipatch)%t
            call p%get_field(Vort,field2d_1,info)
            call p%get_field(Veloc,field3d_1,info)

            do i = 1,p%nnodes(1)
                do j = 1,p%nnodes(2)
                    field2d_1(i,j) = cos(i*h(1)+j)
                    field3d_1(1,i,j) = sin(field2d_1(i,j))
                    field3d_1(2,i,j) = cos(field2d_1(i,j))
                enddo
            enddo
        enddo

        field2d_1 => NULL()
        field3d_1 => NULL()


        ilist => Mesh1%list_of_fields(info)
        if(associated(ilist)) then
            write(*,*) 'mesh1%list: ', ilist
        endif
        ilist => Mesh2%list_of_fields(info)
        write(*,*) 'mesh2%list: ', associated(ilist)

        !--------------------
        ! Remove a patch
        !--------------------
        !Mesh1%remove_patch(patch_ID = 2)
        
        call Vort%destroy(info)
        Assert_equal(info,0)
        call Veloc%destroy(info)
        Assert_equal(info,0)
        call Mesh1%destroy(info)
        Assert_equal(info,0)
        call Mesh2%destroy(info)
        Assert_equal(info,0)

        p=>NULL()

    end test
!------------------------------------------------------------------------------
    test field_basics
        implicit none
        real(mk),dimension(2*ndim)       :: my_patch
        real(mk),dimension(ndim)         :: offset

        type(ppm_t_field)                :: Vort,Veloc
        type(ppm_t_equi_mesh)            :: Mesh1,Mesh2
        class(ppm_t_subpatch_),POINTER   :: p => NULL()

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        assig  = ppm_param_assign_internal
        topoid = 0

        allocate(nm(ndim),stat=info)
        nm(1:ndim) = 16*nproc

        sca_ghostsize = 0.05_mk 

        call ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,    &
        &               bcdef,sca_ghostsize,cost,info)
        Assert_Equal(info,0)
        
        !--------------------------
        !Define Fields
        !--------------------------
        call Vort%create(1,info,name="Vorticity")
        Assert_Equal(info,0)
        call Veloc%create(ndim,info,name="Velocity")
        Assert_Equal(info,0)

        !--------------------------
        !Create Mesh
        !--------------------------
        offset = 0._mk
        call Mesh1%create(topoid,offset,info,Nm=Nm)
        Assert_Equal(info,0)

        mypatchid = 1
        my_patch = (/0.5,0.1,5.1,10.0/)

        !topo => ppm_topo(Mesh1%topoid)%t
        !do isub = 1,topo%nsublist
        !    write(*,*) 'On sub ',isub
        !    write(*,*) '   ',topo%min_subd(1:ndim,isub)
        !    write(*,*) '   ',topo%max_subd(1:ndim,isub)
        !enddo   

        call Mesh1%def_patch(my_patch,info,mypatchid) 
        Assert_Equal(info,0)
        Assert_True(associated(Mesh1%subpatch))

        p => Mesh1%subpatch%begin()
        do while(associated(p))
           Assert_True(associated(p%subpatch_data))
           p => Mesh1%subpatch%next()
        enddo


        !--------------------------
        !Create data arrays on the mesh for the vorticity and velocity fields
        !--------------------------
        call Vort%discretize_on(Mesh1,info)
        Assert_Equal(info,0)

        call Veloc%discretize_on(Mesh1,info)
        Assert_Equal(info,0)

        !--------------------------
        ! Iterate through patches and initialize the data arrays
        !--------------------------
        p => Mesh1%subpatch%begin()

        do while (ASSOCIATED(p))
            call p%get_field(Vort,field2d_1,info)
            call p%get_field(Veloc,field3d_1,info)

            do i = 1,p%nnodes(1)
                do j = 1,p%nnodes(2)
                    field2d_1(i,j) = cos(i*h(1)+j)
                    field3d_1(1:ndim,i,j) = 17.4
                enddo
            enddo
            p => Mesh1%subpatch%next()
        enddo

        !Second version
        do ipatch = 1,Mesh1%subpatch%nb
            p => Mesh1%subpatch%vec(ipatch)%t
            call p%get_field(Vort,field2d_1,info)
            call p%get_field(Veloc,field3d_1,info)

            do i = 1,p%nnodes(1)
                do j = 1,p%nnodes(2)
                    field2d_1(i,j) = cos(i*h(1)+j)
                    field3d_1(1:ndim,i,j) = 17.4
                enddo
            enddo
        enddo

        !--------------------
        ! Remove a patch
        !--------------------
        !Mesh1%remove_patch(patch_ID = 2)

    end test

!------------------------------------------------------------------------------

end test_suite
