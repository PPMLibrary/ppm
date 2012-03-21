test_suite ppm_module_field_typedef

use ppm_module_mesh_typedef
use ppm_module_topo_typedef

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
        
        use ppm_module_finalize

        call ppm_finalize(info)
        
        if (associated(xp)) deallocate(xp)
        if (associated(wp)) deallocate(wp)
        if (associated(cost)) deallocate(cost)
        if (allocated(seed)) deallocate(seed)

    end teardown
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    test create_field
        use ppm_module_data
        use ppm_module_mktopo
        use ppm_module_finalize
        use ppm_module_interfaces

        implicit none
        integer, dimension(2)            :: maxndata
        integer, dimension(:  ), pointer :: isublist => NULL()
        integer                          :: nsublist
        integer                          :: ipatch,isub

        integer                          :: mypatchid
        real(mk),dimension(2*ndim)       :: my_patch
        real(mk),dimension(ndim)         :: offset
        real(mk),dimension(ndim)         :: u_maxsub
        real(mk)                         :: sca_ghostsize

        type(ppm_t_field)                :: Vort,Veloc
        type(ppm_t_equi_mesh)            :: Mesh1,Mesh2
    !type(ppm_t_subpatch),POINTER     :: p => NULL()
        class(ppm_t_subpatch_),POINTER   :: p => NULL()
        type(ppm_t_topo),      POINTER   :: topo => NULL()

        real(mk),dimension(:,:),pointer  :: field2d_1,field2d_2
        real(mk),dimension(:,:,:),pointer:: field3d_1,field3d_2

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        assig  = ppm_param_assign_internal
        topoid = 0

        allocate(nm(ndim),stat=info)
        do i=1,ndim
            nm(i) = 16*nproc
        enddo

        u_maxsub = 0.2_mk
        sca_ghostsize = 0.05_mk 

        call ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,    &
        &               bcdef,sca_ghostsize,cost,info)
        Assert_Equal(info,0)
        
        !--------------------------
        !Create Mesh
        !--------------------------
        offset = 0._mk
        call Mesh1%create(topoid,offset,info,Nm=Nm)
        Assert_Equal(info,0)

        mypatchid = 1
        my_patch = (/0.5,0.1,5.1,10.0/)


        topo => ppm_topo(Mesh1%topoid)%t
        do isub = 1,topo%nsublist
            write(*,*) 'On sub ',isub
            write(*,*) '   ',topo%min_subd(1:ndim,isub)
            write(*,*) '   ',topo%max_subd(1:ndim,isub)
        enddo   

        call Mesh1%add_patch(my_patch,info,mypatchid) 
        Assert_Equal(info,0)


        Assert_True(allocated(Mesh1%subpatch))

        i=1
        p => Mesh1%subpatch%begin()
        do while(associated(p))
           write(*,*) 'subpatch no ',i,' status = ',allocated(p%subpatch_data)
           i=i+1
           p => Mesh1%subpatch%next()
        enddo


        !--------------------------
        !Create data arrays on the mesh for the vorticity and velocity fields
        !--------------------------
        write(*,*) 'in discretize'
        call Vort%discretize_on(Mesh1,info)
        Assert_Equal(info,0)
        write(*,*) 'ok till here'

        call Veloc%discretize_on(Mesh1,info)
        Assert_Equal(info,0)

        !--------------------------
        ! Iterate through patches and initialize the data arrays
        !--------------------------
        p => Mesh1%subpatch%begin()

        DO WHILE (ASSOCIATED(p))

            call p%get_field(field2d_1,Vort,info)
            call p%get_field(field2d_2,Veloc,info)

            DO i = p%istart(1),p%iend(1)
                DO j = p%istart(2),p%iend(2)
                    field2d_1(i,j) = cos(i*h(1)+j)
                    field2d_2(i,j) = sqrt(field2d_1(i,j))
                ENDDO
            ENDDO
            p => Mesh1%subpatch%next()

        ENDDO

        !Second version
        DO ipatch = 1,Mesh1%subpatch%nb
            p => Mesh1%subpatch%vec(ipatch) 
            call p%get_field(field2d_1,Vort,info)
            call p%get_field(field2d_2,Veloc,info)

            DO i = p%istart(1),p%iend(1)
                DO j = p%istart(2),p%iend(2)
                    field2d_1(i,j) = cos(i*h(1)+j)
                    field2d_2(i,j) = sqrt(field2d_1(i,j))
                ENDDO
            ENDDO
        ENDDO

        !--------------------
        ! Remove a patch
        !--------------------
        !Mesh1%remove_patch(patch_ID = 2)


        
 
    end test
!------------------------------------------------------------------------------

end test_suite
